"""Validate saved Qollab reports and export matched optimizer data and figures.

Reads HTML as text, never executes it. CSV collection uses standard Python.
Optional figures use installed base R, with no new Python plotting dependencies.
"""

import argparse
import csv
from html.parser import HTMLParser
import json
from math import isclose, isfinite
from pathlib import Path
import shutil
import subprocess
import sys

from qphylo_benchmark import benchmark_metrics
from qphylo_qaoa import QAOAConfig, summarize_counts


def read_report(path):
    path = Path(path)
    source = path.read_text()
    if path.suffix.lower() == ".json":
        return json.loads(source)
    if path.suffix.lower() == ".html":
        class Messages(HTMLParser):
            def __init__(self): super().__init__(); self.parts = []
            def handle_data(self, data): self.parts.append(data)
        parser = Messages()
        parser.feed(source)
        source = "\n".join(parser.parts)
    matches = [line.split("REPORT_JSON:", 1)[1].strip()
               for line in source.splitlines() if "REPORT_JSON:" in line]
    if not matches:
        raise ValueError(f"No REPORT_JSON in {path.name}")
    return json.loads(matches[-1])


def run_inventory(paths):
    """Quickly show recorded settings; filenames are not experimental metadata."""
    rows = []
    for path in paths:
        report = read_report(path)
        config = report["config"]
        rows.append({"file": str(path), "optimizer": config.get("optimizer", "seeded_coordinate_search"),
                     "stored_seed": config["seed"], "budget": config["max_evaluations"],
                     "simulator_seed": report.get("simulator_seed"), "status": report["status"]})
    return rows


def summaries_agree(expected, stored):
    """Allow floating-point roundoff only; discrete sample fields stay exact."""
    if isinstance(expected, float):
        return (type(stored) in (int, float) and isfinite(stored) and
                isclose(expected, stored, rel_tol=1e-12, abs_tol=1e-12))
    if type(expected) is not type(stored):
        return False
    if isinstance(expected, dict):
        return expected.keys() == stored.keys() and all(
            summaries_agree(value, stored[key]) for key, value in expected.items())
    if isinstance(expected, list):
        return len(expected) == len(stored) and all(
            summaries_agree(a, b) for a, b in zip(expected, stored))
    return expected == stored


def validated_report(path):
    report = read_report(path)
    if report.get("mode") not in ("toy_tree", "protein_tree"):
        raise ValueError(f"{Path(path).name}: expected a toy_tree or protein_tree run.")
    if report.get("execution_mode", "multi_job") != "multi_job":
        raise ValueError("Comparison collector currently accepts completed/failed local multi_job runs.")
    if report["status"] not in ("completed", "failed"):
        raise ValueError("Run is unfinished; do not treat a checkpoint as a replicate.")
    config = QAOAConfig(**report["config"])
    for row in report.get("evaluations", []):
        stats = summarize_counts(row["weights"], row["target_size"], row["penalty"],
                                 row["counts"], row["shots"])
        if not summaries_agree(json.loads(json.dumps(stats)), row["summary"]):
            raise ValueError(f"Stored scores disagree with counts in {Path(path).name}.")
        if row["backend"].get("status") != "DONE":
            raise ValueError("A scored evaluation has no completed job.")
        if row.get("optimizer", "seeded_coordinate_search") != config.optimizer:
            raise ValueError("Evaluation optimizer disagrees with run configuration.")
    report["config"] = vars(config)
    benchmark_metrics(report)
    report["source_file"] = str(Path(path).resolve())
    return report


def comparison_tables(paths):
    reports = [validated_report(path) for path in paths]
    if not reports:
        raise ValueError("No saved result files were supplied.")
    signatures, backend_signatures, versions, seen = set(), set(), set(), set()
    for report in reports:
        config = dict(report["config"])
        optimizer, seed = config.pop("optimizer"), config.pop("seed")
        pair_key = (optimizer, seed)
        if pair_key in seen:
            conflicting = [r["source_file"] for r in reports
                           if (r["config"]["optimizer"], r["config"]["seed"]) == pair_key]
            raise ValueError(f"Duplicate optimizer/seed {pair_key} in these reports:\n  " +
                             "\n  ".join(conflicting) +
                             "\nSeeds come from REPORT_JSON, not filenames. Select one replicate "
                             "explicitly, or rerun missing seeds after changing the top-level seed setting.")
        seen.add(pair_key)
        jobs = report.get("jobs", [])
        actual = sorted(set(j.get("backend_type", "unknown") for j in jobs))
        noise = sorted(set(str(j.get("noise_model")) for j in jobs))
        qiskit = sorted(set(j.get("qiskit_version", "unknown") for j in jobs))
        signatures.add(json.dumps({"graph": report["graph_signature"], "config": config}, sort_keys=True))
        if jobs:
            if actual != ["BasicSimulator"] or noise != ["None"]:
                raise ValueError("This first comparison requires the built-in ideal BasicSimulator; keep noise runs separate.")
            if report.get("simulator_seed") is not None and any(
                    j.get("seed_simulator") != report["simulator_seed"] for j in jobs):
                raise ValueError("Recorded simulator seed disagrees with actual backend options.")
            backend_signatures.add(json.dumps({"backends": actual, "noise": noise,
                                               "qiskit": qiskit}, sort_keys=True))
        if config["max_evaluations"] < 2:
            raise ValueError("Invalid budget.")
        versions.update(e["scipy_version"] for e in report.get("evaluations", [])
                        if e.get("scipy_version"))
    if len(signatures) != 1 or len(backend_signatures) > 1:
        raise ValueError("Unmatched graphs/settings/actual backends: keep smoke runs and comparison runs separate.")
    if len(versions) > 1:
        raise ValueError("Mixed SciPy versions; do not pool these COBYLA runs.")
    by_seed = {}
    for r in reports:
        by_seed.setdefault(r["config"]["seed"], []).append(r)
    paired_seeds = []
    for seed, pair in by_seed.items():
        if len(pair) == 2:
            if pair[0].get("simulator_seed") != pair[1].get("simulator_seed"):
                raise ValueError(f"Paired seed {seed} has different simulator seeds.")
            initial = []
            for report in pair:
                initial.append({(e["target_size"], e["attempt"]): e.get("initial_parameters")
                                for e in report.get("evaluations", [])
                                if len(e["weights"]) == len(report["labels"]) and
                                e["phase"] == "training" and e["evaluation"] == 0})
            for key in initial[0].keys() & initial[1].keys():
                if initial[0][key] != initial[1][key]:
                    raise ValueError(f"Paired seed {seed} has different root initial angles.")
            paired_seeds.append(seed)
    runs, evaluations, finals = [], [], []
    for report in reports:
        cfg = report["config"]
        root_n = len(report["labels"])
        base = {"optimizer": cfg["optimizer"], "seed": cfg["seed"],
                "dataset_kind": "protein_derived" if report["mode"] == "protein_tree" else "synthetic",
                "source_file": report["source_file"]}
        resources = report["resource_metrics"]
        topology = report["topology_metrics"] or {}
        runs.append({**base, "status": report["status"], "n": root_n,
            "budget_per_attempt": cfg["max_evaluations"], "shots_per_job": cfg["shots"],
            "p": cfg["p"],
            "scipy_version": ",".join(sorted(set(e["scipy_version"] for e in report.get("evaluations", []) if e.get("scipy_version")))),
            "simulator_seed": report.get("simulator_seed"),
            "actual_backend": ",".join(sorted(set(j.get("backend_type", "unknown") for j in report.get("jobs", [])))),
            "rf_to_exact": topology.get("unrooted_rf_to_exact"),
            "matches_exact": topology.get("matches_exact_splits"),
            "split_nodes": topology.get("split_nodes"), **resources,
            "reinforcements": sum(a["action"] == "reinforce" for a in report.get("attempts", [])),
            "tree": report.get("tree"), "exact_tree": report["exact_tree"],
            "error": report.get("error")})
        best_values = {}
        for row in report.get("evaluations", []):
            value = row["summary"]["expectation"] / row["scale"]
            key = (row["call"], row["attempt"])
            if row["phase"] == "training":
                best_values[key] = min(best_values.get(key, float("inf")), value)
            evaluations.append({**base, "call": row["call"], "attempt": row["attempt"],
                "evaluation": row["evaluation"], "phase": row["phase"],
                "n": len(row["weights"]), "root_n": root_n,
                "target_size": row["target_size"], "penalty": row["penalty"],
                "scaled_expectation": value, "best_training_so_far": best_values.get(key),
                "feasible_rate": row["summary"]["feasible_rate"],
                "job_id": row["backend"].get("job_id"),
                "scipy_version": row.get("scipy_version"),
                "initial_parameters": json.dumps(row.get("initial_parameters")),
                "parameters": json.dumps(row["parameters"])})
        for row in report["final_cut_metrics"]:
            finals.append({**base, "root_n": root_n, **row,
                           "taxa": ",".join(row["taxa"]) if row["taxa"] else ""})
    return reports, {"runs": runs, "evaluations": evaluations, "final_cuts": finals}, sorted(paired_seeds)


def export_tables(paths, output):
    reports, tables, paired_seeds = comparison_tables(paths)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for name, rows in tables.items():
        if not rows:
            # Remove only our derived table, avoiding stale measurements if this
            # export replaces a previous completed study with failed runs.
            (output / f"{name}.csv").unlink(missing_ok=True)
            continue
        with (output / f"{name}.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    manifest = {"source_files": [r["source_file"] for r in reports],
        "graph_signature": reports[0]["graph_signature"], "paired_seeds": paired_seeds,
        "replicates": len(reports), "failed_runs": sum(r["status"] == "failed" for r in reports),
        "measured_evaluations": len(tables["evaluations"]),
        "final_measurements": len(tables["final_cuts"]),
        "dataset_provenance": reports[0].get("dataset_provenance"),
        "scope": ("Protein-derived affinity" if reports[0]["mode"] == "protein_tree" else "Synthetic") +
            " ideal-simulator optimizer pilot; no biological-accuracy or quantum-advantage claim.",
        "curve": "Root c=n/2, attempt 0, scaled penalized training energy versus measured evaluations.",
        "timing": "Sum of recorded job elapsed times, not full program wall time.",
        "uniform_baseline": "Analytic uniform bitstring sampling, not a measured random-solver run."}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", help="Saved HTML/JSON reports or a folder of such reports")
    parser.add_argument("--output", default="results/eight_taxa_comparison")
    parser.add_argument("--plot", action="store_true", help="Export PNG/PDF with installed Rscript")
    parser.add_argument("--inspect", action="store_true", help="List recorded seeds/settings without exporting or plotting")
    args = parser.parse_args()
    paths = []
    for item in args.inputs:
        path = Path(item)
        paths.extend(sorted(p for p in path.iterdir() if p.suffix.lower() in (".html", ".json"))
                     if path.is_dir() else [path])
    try:
        if args.inspect:
            rows = run_inventory(paths)
            if rows:
                writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            raise SystemExit(0)
        manifest = export_tables(paths, args.output)
    except ValueError as exc:
        parser.exit(2, f"Collection stopped: {exc}\n")
    print(json.dumps(manifest, indent=2))
    if args.plot:
        if not manifest["measured_evaluations"] or not manifest["final_measurements"]:
            raise SystemExit("Tables saved; there are no final measurements to plot yet.")
        rscript = shutil.which("Rscript")
        if not rscript:
            raise SystemExit("CSV files are saved. Rscript is required for --plot.")
        subprocess.run([rscript, str(Path(__file__).with_name("plot_qaoa_results.R")),
                        str(Path(args.output).resolve())], check=True)

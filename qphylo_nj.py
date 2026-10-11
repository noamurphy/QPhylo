"""Run the local WAG+Gamma neighbor-joining baseline and tree comparison.

Requires R with ape/phangorn in the qphylo Conda environment.
From the project root: python qphylo_nj.py
No quantum jobs are submitted. Reuses the supplied alignment, not new MUSCLE output.
"""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from collect_qaoa_results import validated_report
from qphylo_proteins import read_protein_fasta

ROOT = Path(__file__).resolve().parent
DEFAULT_FASTA = ROOT / "PhyloTree_ReConst/alignments/8/assembled_sequences_1822.fasta"
DEFAULT_REPORT = ROOT / "qollab_smoke/protein1822_cobyla_run1.html"
DEFAULT_REFERENCE = ROOT / "PhyloTree_ReConst/trees/8/tree_best_1822.newick"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_alignment(path):
    """Validate FASTA while retaining aligned columns and gap positions."""
    proteins = read_protein_fasta(path)
    chunks, sequences = [], []
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        if not line.strip():
            continue
        if line.strip().startswith(">"):
            if chunks:
                sequences.append("".join(chunks).upper().replace(".", "-"))
            chunks = []
        else:
            chunks.append("".join(line.split()))
    if chunks:
        sequences.append("".join(chunks).upper().replace(".", "-"))
    if len(sequences) != len(proteins):
        raise ValueError("Invalid alignment records.")
    if len({len(s) for s in sequences}) != 1:
        raise ValueError("NJ needs an alignment with equal column counts; align sequences first.")
    if len(proteins) < 3:
        raise ValueError("NJ requires at least three taxa.")
    for protein, sequence in zip(proteins, sequences):
        if set(sequence) - set("ACDEFGHIKLMNPQRSTVWY-"):
            raise ValueError("This pilot accepts canonical amino acids and alignment gaps only.")
        if sequence.replace("-", "") != protein.sequence:
            raise ValueError("Alignment and ungapped sequence disagree.")
    return [p.label for p in proteins], sequences


def find_rscript(explicit=None):
    if explicit:
        path = shutil.which(explicit) or explicit
        if Path(path).is_file():
            return str(Path(path).resolve())
        raise FileNotFoundError(f"Rscript not found: {explicit}")
    # Prefer project-environment R over the unrelated system installation.
    candidates = [Path(sys.executable).resolve().parent / "Rscript"]
    if os.environ.get("CONDA_PREFIX"):
        candidates.insert(0, Path(os.environ["CONDA_PREFIX"]) / "bin/Rscript")
    conda = os.environ.get("CONDA_EXE") or shutil.which("conda")
    if conda:
        candidates.append(Path(conda).resolve().parents[1] / "envs/qphylo/bin/Rscript")
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    path = shutil.which("Rscript")
    if path:
        return path
    raise FileNotFoundError("Install R, ape and phangorn in qphylo, or supply --rscript.")


def run(fasta=DEFAULT_FASTA, report_path=DEFAULT_REPORT, reference=DEFAULT_REFERENCE,
        output=ROOT / "results/protein1822_comparison", bootstraps=100, seed=42, rscript=None):
    if not isinstance(bootstraps, int) or isinstance(bootstraps, bool) or bootstraps < 1:
        raise ValueError("At least one bootstrap replicate is required.")
    if type(seed) is not int or not 0 <= seed <= 2147483647:
        raise ValueError("Use an integer bootstrap seed from 0 to 2147483647.")
    labels, sequences = read_alignment(fasta)
    report = validated_report(report_path)
    if report["mode"] != "protein_tree" or report["status"] != "completed":
        raise ValueError("Use a completed protein-derived QAOA report.")
    provenance = report.get("dataset_provenance", {}).get("preprocessing", {})
    if labels != report["labels"] or sha256(fasta) != provenance.get("source_sha256"):
        raise ValueError("Alignment taxa/order or source hash differs from the QAOA input.")
    reference_text = Path(reference).read_text().strip()
    executable = find_rscript(rscript)
    folder = Path(output).resolve()
    folder.mkdir(parents=True, exist_ok=True)
    owned = ("alignment.fasta", "labels.txt", "qaoa_settings.csv", "qaoa_tree.newick", "exact_tree.newick",
             "reference_tree.newick", "nj_tree.newick", "nj_tree_with_support.newick",
             "nj_bootstrap_trees.newick", "distances.csv", "rf_comparison.csv", "rf_matrix.csv",
             "splits.csv", "tree_comparison.png", "versions.csv", "nj_diagnostics.csv",
             "sessionInfo.txt", "analysis.stdout.txt", "analysis.stderr.txt", "summary.md", "manifest.json")
    inputs = {Path(p).resolve() for p in (fasta, report_path, reference)}
    if any(folder / name in inputs for name in owned):
        raise ValueError("Output would overwrite a source input.")
    for name in owned:
        (folder / name).unlink(missing_ok=True)
    (folder / "alignment.fasta").write_text("".join(f">{label}\n{seq}\n" for label, seq in zip(labels, sequences)))
    (folder / "labels.txt").write_text("\n".join(labels) + "\n")
    with (folder / "qaoa_settings.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["p", "seed"])
        writer.writeheader()
        writer.writerow({key: report["config"][key] for key in ("p", "seed")})
    for name, text in (("qaoa_tree", report["tree"]), ("exact_tree", report["exact_tree"]),
                       ("reference_tree", reference_text)):
        (folder / f"{name}.newick").write_text(text + "\n")
    manifest = {"schema": "qphylo-model-nj-v1", "status": "prepared",
        "alignment": {"source": str(Path(fasta).resolve()), "source_sha256": sha256(fasta),
            "prepared_sha256": sha256(folder / "alignment.fasta"), "labels": labels,
            "n_taxa": len(labels), "columns": len(sequences[0]),
            "gap_characters": sum(s.count("-") for s in sequences),
            "method": "Reuse supplied FO alignment; no new MUSCLE alignment."},
        "model": {"name": "WAG", "gamma_categories": 5, "gamma_shape": 0.5,
            "frequencies": "WAG model frequencies, not empirical +F", "invariant_sites": "No separate +I component",
            "gap_handling": "phangorn dist.ml exclude=pairwise"},
        "bootstrap": {"replicates": bootstraps, "seed": seed, "method": "Site resampling and NJ refitting",
            "jumble": True, "multicore": False, "support_units": "Percent; phangorn FBP proportions multiplied by 100"},
        "implementation_sha256": {name: sha256(ROOT / name) for name in ("qphylo_nj.py", "build_model_nj.R")},
        "qaoa": {"source": str(Path(report_path).resolve()), "source_sha256": sha256(report_path),
            "graph_signature": report["graph_signature"], "config": report["config"],
            "resources": report["resource_metrics"], "dataset_provenance": report["dataset_provenance"]},
        "reference": {"source": str(Path(reference).resolve()), "sha256": sha256(reference),
            "meaning": "Provided FO tree_best_1822; generating-tree provenance not established by this pilot."},
        "comparison": "Unrooted RF with identical taxa; branch lengths ignored; display anchored on last taxon only",
        "sources": ["Onodera paper section 2.8: WAG+Gamma5, shape0.5, NJ and100 bootstraps",
            "https://klausvigo.github.io/phangorn/reference/dist.hamming.html",
            "https://klausvigo.github.io/phangorn/reference/bootstrap.pml.html"],
        "limits": "Provided alignment/newer libraries differ from MUSCLE3.8.31/phangorn2.8.1. One graph and QAOA seed; no quantum advantage claim."}
    manifest_path = folder / "manifest.json"
    def save():
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    save()
    started = time.monotonic()
    try:
        command = [executable, str(ROOT / "build_model_nj.R"), str(folder), str(bootstraps), str(seed)]
        manifest["command"] = command
        result = subprocess.run(command, capture_output=True, text=True)
        (folder / "analysis.stdout.txt").write_text(result.stdout)
        (folder / "analysis.stderr.txt").write_text(result.stderr)
        if result.returncode:
            raise RuntimeError(result.stderr.strip() or f"R analysis exited {result.returncode}")
        for name in ("nj_tree.newick", "nj_tree_with_support.newick", "nj_bootstrap_trees.newick",
                     "distances.csv", "rf_comparison.csv", "rf_matrix.csv", "splits.csv", "tree_comparison.png",
                     "versions.csv", "nj_diagnostics.csv", "sessionInfo.txt"):
            if not (folder / name).is_file():
                raise RuntimeError(f"R analysis did not produce {name}.")
        with (folder / "versions.csv").open() as stream:
            manifest["versions"] = {row["name"]: row["value"] for row in csv.DictReader(stream)}
        with (folder / "rf_comparison.csv").open() as stream:
            comparisons = list(csv.DictReader(stream))
        with (folder / "nj_diagnostics.csv").open() as stream:
            manifest["nj_diagnostics"] = next(csv.DictReader(stream))
        with (folder / "splits.csv").open() as stream:
            nj_splits = [row for row in csv.DictReader(stream) if row["method"] == "Model_NJ"]
        manifest["pairwise_comparisons"] = comparisons
        manifest["status"] = "completed"
        manifest["output_sha256"] = {name: sha256(folder / name) for name in owned
                                     if (folder / name).is_file() and name not in ("manifest.json", "summary.md")}
        summary = ["# Protein1822: QAOA, exact NMcut and model-based NJ", "",
            f"NJ used WAG+Gamma (5 categories, shape0.5), {bootstraps} site bootstraps, seed{seed}.",
            f"Input: the same {len(labels)} proteins / {len(sequences[0])} aligned columns used for the QAOA graph.", "",
            "| Method 1 | Method 2 | Unrooted RF | Normalized RF |", "| --- | --- | --- | --- |"]
        summary += [f"| {row['method1']} | {row['method2']} | {row['rf']} | {row['normalized_rf']} |" for row in comparisons]
        pair_rf = {(row["method1"], row["method2"]): float(row["rf"]) for row in comparisons}
        if pair_rf[("QAOA", "Exact_NMcut")] == 0 and pair_rf[("QAOA", "Model_NJ")] > 0:
            summary += ["", "QAOA reproduces the exact affinity-tree topology on this run. Its disagreement with NJ is also present in exact NMcut, so it is not a QAOA topology approximation error on this input."]
        summary += ["", "RF0 means identical unrooted splits. Reference agreement is conditional on the supplied reference; its generating-tree provenance is unverified.",
            "", "## NJ site-bootstrap support", "", "| Split sides | Support (%) |", "| --- | --- |"]
        summary += [f"| {row['side1']} / {row['side2']} | {row['bootstrap_percent']} |" for row in nj_splits]
        summary += ["", f"NJ produced {manifest['nj_diagnostics']['negative_nj_edges']} negative branch lengths. These are retained in the Newick files; topology comparisons and the figure ignore branch lengths.",
            "", "NJ bootstrap percentages describe site-resampling support, not probabilities that clades are correct.",
            "", "The figure compares topology only: branch lengths are hidden and the last taxon is a display anchor, not a biological outgroup.",
            "", "This reuses the provided alignment and records current R/library versions; it is not a complete reproduction of Onodera's alignment/tool versions.",
            "", "![Tree comparison](tree_comparison.png)"]
        (folder / "summary.md").write_text("\n".join(summary) + "\n")
        print(result.stdout.strip())
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        manifest["elapsed_seconds"] = time.monotonic() - started
        save()
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--qaoa-report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--output", type=Path, default=ROOT / "results/protein1822_comparison")
    parser.add_argument("--bootstraps", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--rscript", help="Optional Rscript executable; prefers R in qphylo")
    args = parser.parse_args()
    run(args.fasta, args.qaoa_report, args.reference, args.output, args.bootstraps, args.seed, args.rscript)
    print("Saved:", args.output / "summary.md")

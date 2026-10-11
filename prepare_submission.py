"""Freeze audited protein data and comparison into a self-contained Qollab main.

Run: python prepare_submission.py
This builds code only. No quantum jobs are submitted; local reference files stay intact.
"""

import argparse
import ast
import csv
import hashlib
import json
from pathlib import Path
import shutil

from collect_qaoa_results import validated_report
from prepare_qollab import protein_graph, without_local_imports_and_cli
from qphylo_benchmark import graph_signature
from qphylo_qaoa import hamiltonian_scale

ROOT = Path(__file__).resolve().parent


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def submission_data(graph_path, report_path, comparison_path):
    weights, labels, provenance = protein_graph(graph_path)
    report = validated_report(report_path)
    signature = graph_signature(weights, labels)
    if (report["status"] != "completed" or report["mode"] != "protein_tree" or
            report["graph_signature"] != signature or report["dataset_provenance"] != provenance):
        raise ValueError("Training report does not match the audited protein graph/provenance.")
    jobs = report["jobs"]
    if not jobs or any(job.get("backend_type") != "BasicSimulator" or
                       job.get("status") != "DONE" or job.get("noise_model") is not None for job in jobs):
        raise ValueError("Freeze angles only from a completed built-in ideal-simulator run.")
    target_size = len(labels)//2
    finals = [row for row in report["evaluations"] if row["call"] == target_size-1 and
              row["phase"] == "final" and row["weights"] == weights and row["target_size"] == target_size]
    accepted = [a for a in report["attempts"] if a["call"] == target_size-1 and a["action"] == "accept"]
    if len(accepted) != 1:
        raise ValueError("No unique accepted balanced-root cut in training report.")
    attempt = accepted[0]
    final = next((row for row in finals if row["attempt"] == attempt["attempt"]), None)
    if (final is None or not final["summary"]["best_energy_sample"]["feasible"] or
            final["parameters"] != attempt["parameters"] or final["penalty"] != attempt["penalty"] or
            not abs(final["scale"]-hamiltonian_scale(weights,target_size,final["penalty"])) < 1e-9):
        raise ValueError("Accepted cut angles/scale are not supported by the saved final measurement.")
    training = [row for row in report["evaluations"] if row["call"] == final["call"] and
                row["attempt"] == final["attempt"] and row["phase"] == "training"]
    best_training = min(training, key=lambda row: row["summary"]["expectation"]/row["scale"])
    if final["parameters"] != best_training["parameters"]:
        raise ValueError("Final angles do not match the best measured training point.")
    folder = Path(comparison_path)
    manifest = json.loads((folder/"manifest.json").read_text())
    if (manifest["status"] != "completed" or manifest["qaoa"]["source_sha256"] != file_hash(report_path) or
            manifest["qaoa"]["graph_signature"] != signature or manifest["alignment"]["labels"] != labels):
        raise ValueError("NJ comparison does not match the training report and protein taxa.")
    for name, digest in manifest["output_sha256"].items():
        if Path(name).name != name or file_hash(folder/name) != digest:
            raise ValueError(f"Comparison output hash mismatch: {name}")
    if ((folder/"qaoa_tree.newick").read_text().strip() != report["tree"] or
            (folder/"exact_tree.newick").read_text().strip() != report["exact_tree"]):
        raise ValueError("Saved comparison trees disagree with the verified run.")
    with (folder/"splits.csv").open() as stream:
        splits = list(csv.DictReader(stream))
    def canonical_splits(method):
        return sorted([list(min((tuple(sorted(row["side1"].split(','))),
                                     tuple(sorted(row["side2"].split(',')))), key=lambda x:(len(x),x)))
                       for row in splits if row["method"] == method])
    return {"schema": "qphylo-submission-data-v1", "weights": weights, "labels": labels,
        "graph_signature": signature, "dataset_provenance": provenance,
        "config": report["config"], "training_simulator_seed": report["simulator_seed"],
        "frozen_cut": {"target_size": target_size, "penalty": final["penalty"],
            "parameters": final["parameters"], "scale": final["scale"], "graph_signature": signature,
            "source_file": Path(report_path).name, "source_sha256": file_hash(report_path),
            "training_config": report["config"], "training_simulator_seed": report["simulator_seed"],
            "call": final["call"], "attempt": final["attempt"],
            "meaning": "Angles trained on ideal BasicSimulator; one new cut measurement, not a full tree or noisy reoptimization"},
        "saved_comparison": {"origin": "Previously verified local protein1822 run and model NJ",
            "n_taxa": len(labels), "columns": manifest["alignment"]["columns"],
            "submitted_jobs": len(jobs), "reinforcements": sum(a["action"] == "reinforce" for a in report["attempts"]),
            "qaoa_source_sha256": file_hash(report_path), "manifest_sha256": file_hash(folder/"manifest.json"),
            "pairwise_comparisons": manifest["pairwise_comparisons"],
            "qaoa_tree": report["tree"], "exact_tree": report["exact_tree"],
            "nj_tree": (folder/"nj_tree_with_support.newick").read_text().strip(),
            "reference_tree": (folder/"reference_tree.newick").read_text().strip(),
            "nj_splits": canonical_splits("Model_NJ"), "reference_splits": canonical_splits("Provided_reference"),
            "nj_support": [{"side": row["side1"].split(','), "percent": float(row["bootstrap_percent"])}
                           for row in splits if row["method"] == "Model_NJ"],
            "model": manifest["model"], "bootstrap": manifest["bootstrap"], "versions": manifest["versions"],
            "reference_meaning": manifest["reference"]["meaning"]}}


def build_submission(output=None, graph=None, report=None, comparison=None):
    graph = Path(graph or ROOT/"results/protein1822/graph.json")
    report = Path(report or ROOT/"qollab_smoke/protein1822_cobyla_run1.html")
    comparison = Path(comparison or ROOT/"results/protein1822_comparison")
    data = submission_data(graph, report, comparison)
    parts = ["# GENERATED by prepare_submission.py; upload as Qollab main.py.\n"
             "# Copyright (c) 2026 QPhylo contributors\n"
             "# SPDX-License-Identifier: MIT\n"
             "# See LICENSE for terms and CREDITS.md for attribution/upstream notices.\n"
             "# Code/data are self-contained; no main() wrapper is needed.\n"
             "# Auto: Built-in simulator -> fresh full tree; IonQ/others -> one frozen-angle cut.\n"
             "# Saved comparison is labeled historical; current jobs are always separately reported.\n"
             "# This effort is supported by Qollab & IonQ\n"
             "# Method: https://doi.org/10.1016/j.ympev.2022.107636\n"
             "# Prior QAOA work/dataset: https://github.com/Nicofero/PhyloTree_ReConst\n\n"
             "RUN_MODE = 'auto'  # auto, full_tree (built-in only), or single_cut.\n\n",
             f"SUBMISSION_DATA = {data!r}\n"]
    for source in ("qphylo_core.py", "qphylo_qaoa.py", "qphylo_benchmark.py", "qollab_smoke/submission_entry.py"):
        parts.append(f"\n# ---- Source: {source} ----\n")
        parts.append(without_local_imports_and_cli(ROOT/source))
    code = "".join(parts)
    ast.parse(code)
    output = Path(output or ROOT/"qollab_smoke/submission_main.py")
    # Preserve the input snapshots and maintainable sources if an output override is supplied.
    if (output.resolve() in (graph.resolve(), report.resolve()) or
            output.resolve().parent == comparison.resolve() or
            output.resolve() in [(ROOT/name).resolve() for name in
                ("qphylo_core.py", "qphylo_qaoa.py", "qphylo_benchmark.py", "qollab_smoke/submission_entry.py", "prepare_submission.py")]):
        raise ValueError("Output would overwrite a source or comparison result.")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(code)
    return output


def build_upload(main_path, folder=None):
    """Prepare code/notices plus the maintained qollab_upload README/assets.

    The default README/assets are source artifacts, copied for custom targets.
    Regenerate plots separately with prepare_project_header.R when data changes.
    """
    folder = Path(folder or ROOT/"qollab_upload")
    folder.mkdir(parents=True, exist_ok=True)
    for source, name in ((Path(main_path), "main.py"), (ROOT/"LICENSE", "LICENSE"),
                         (ROOT/"CREDITS.md", "CREDITS.md")):
        target = folder/name
        if source.resolve() != target.resolve():
            shutil.copy2(source, target)
    card = ROOT/"qollab_upload/README.md"
    if card.resolve() != (folder/"README.md").resolve():
        shutil.copy2(card, folder/"README.md")
    assets = ROOT/"qollab_upload/assets"
    if assets.resolve() != (folder/"assets").resolve():
        shutil.copytree(assets, folder/"assets", dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns(".ipynb_checkpoints", ".DS_Store"))
    return folder


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--upload-dir", type=Path, help="Default: qollab_upload with main.py, notices, README and assets")
    args = parser.parse_args()
    main_path = build_submission(output=args.output)
    print(main_path)
    print("Project upload set:", build_upload(main_path, args.upload_dir))

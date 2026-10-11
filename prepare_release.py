"""Copy a curated, reviewable local release; never publish or submit jobs.

Run: python prepare_release.py
The target must not exist, so previously edited release copies are preserved.
"""

import argparse
import ast
import hashlib
import json
from pathlib import Path
import shutil

from prepare_submission import build_submission, build_upload, ROOT


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def release_files():
    files = [ROOT/name for name in (
        "README.md", "LICENSE", "CREDITS.md", "DESIGN_DECISIONS.md", "PROJECT_HANDOFF.md", "AGENTS.md",
        "environment.yml", ".gitignore",
        "qphylo_core.py", "qphylo_qaoa.py", "qphylo_run_steps.py", "qphylo_benchmark.py",
        "qphylo_proteins.py", "qphylo_nj.py", "build_model_nj.R", "collect_qaoa_results.py",
        "plot_qaoa_results.R", "prepare_project_header.R", "prepare_qollab.py", "prepare_submission.py", "prepare_release.py",
        "qollab_smoke/README.md", "qollab_smoke/EIGHT_TAXA_GUIDE.md",
        "qollab_smoke/submission_entry.py", "qollab_smoke/submission_main.py",
        "qollab_smoke/cut_pilot.py", "qollab_smoke/qaoa_cut_standalone.py",
        "qollab_smoke/protein1822_cobyla.py", "qollab_smoke/eight_smoke.py",
        "qollab_smoke/eight_coordinate.py", "qollab_smoke/eight_cobyla.py", "qollab_smoke/next_checks.py",
        "qollab_smoke/protein1822_cobyla_run1.html",
        "qollab_smoke/qaoa_toyTestResults_noamurphy-qaoa-phylogenetic-tree-reconst-output-completed-2026-10-10T19-58-08.397Z.html",
        "qollab_smoke/qaoa_cut_standaloneResults_noamurphy-qaoa-phylogenetic-tree-reconst-output-completed-2026-10-10T18-13-10.946Z.html",
        "qollab_smoke/noamurphy-qaoa-phylogenetic-tree-reconst-output-completed-2026-10-10T19-32-54.722Z.html",
        "PhyloTree_ReConst/LICENSE",
        "PhyloTree_ReConst/alignments/8/assembled_sequences_1822.fasta",
        "PhyloTree_ReConst/trees/8/tree_best_1822.newick")]
    for folder in ("docs", "tests", "qollab_upload", "results/protein1822", "results/protein1822_cobyla_run1",
                   "results/protein1822_comparison", "results/eight_taxa_comparison", "results/toy4_validation",
                   "results/submission_validation",
                   "qollab_smoke/eight_taxa_results/comparison", "qollab_smoke/Next Check Results"):
        files.extend(p for p in (ROOT/folder).rglob("*") if p.is_file() and
                     "__pycache__" not in p.parts and ".ipynb_checkpoints" not in p.parts and
                     p.name not in {".DS_Store", ".Rhistory"})
    # Include final live output only if it has actually been supplied.
    files.extend((ROOT/"qollab_smoke").glob("submission_*run*.html"))
    return sorted(set(files))


def prepare_release(output=None):
    destination = Path(output or ROOT/"submission_release").resolve()
    if destination.exists():
        raise FileExistsError("Release target exists. Choose a new --output folder to preserve edits.")
    # Regenerate and revalidate the submitted entry before making a snapshot.
    build_upload(build_submission())
    files = release_files()
    if any(not path.is_file() for path in files):
        raise ValueError("A required source/evidence file is missing.")
    for path in files:
        if path.suffix == ".py":
            ast.parse(path.read_text(), filename=str(path))
    inventory = {str(path.relative_to(ROOT)): digest(path) for path in files}
    destination.mkdir(parents=True)
    for source in files:
        target = destination/source.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        if digest(target) != inventory[str(source.relative_to(ROOT))]:
            raise RuntimeError(f"Copy verification failed: {target}")
    (destination/"PhyloTree_ReConst/README_SOURCE.md").write_text(
        "# Selected upstream material\n\n"
        "Only the chosen alignment, supplied reference tree and original MIT license\n"
        "are copied here, unchanged. This is not the full upstream code repository.\n"
        "Source: https://github.com/Nicofero/PhyloTree_ReConst\n"
        "Checkout: 2fba9f1904a13fa41246af9cef25574f367f0239\n"
        "See ../CREDITS.md for attribution, RAxMLGrove ODbL/contents notices and\n"
        "the unresolved selected-data/generating-tree provenance.\n")
    notice = destination/"PhyloTree_ReConst/README_SOURCE.md"
    inventory[str(notice.relative_to(destination))] = digest(notice)
    final_checks = "Pending; inspect any subsequently supplied final-entry reports"
    audit_path = ROOT/"results/submission_validation/audit.json"
    if audit_path.exists():
        audit = json.loads(audit_path.read_text())
        if (audit.get("status") == "verified" and
                audit.get("local_prepared_main_sha256") == digest(ROOT/"qollab_upload/main.py") and
                all((ROOT/"qollab_smoke"/name).is_file() and
                    digest(ROOT/"qollab_smoke"/name) == expected
                    for name, expected in audit["sources"].items())):
            final_checks = "Verified: BasicSimulator full tree and IonQ Aria-1 noisy cut; see results/submission_validation/audit.json"
    (destination/"RELEASE_MANIFEST.json").write_text(json.dumps({
        "schema": "qphylo-release-v1", "status": "local_preparation_only",
        "entry_point": "qollab_smoke/submission_main.py",
        "final_live_checks": final_checks,
        "excluded": ["full upstream fork/Git internals", "research PDFs", "scratch outputs", "environments"],
        "source_file_sha256": inventory}, indent=2, sort_keys=True)+"\n")
    return destination, len(inventory)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    folder, count = prepare_release(args.output)
    print(f"Prepared {count} verified files in {folder}; nothing published.")

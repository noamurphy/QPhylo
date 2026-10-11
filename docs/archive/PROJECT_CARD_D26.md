> Historical D26 field draft, superseded by D29. Use ../../qollab_upload/README.md.

# Qollab Project Card: ready-to-copy content

Keep the current project's identity if preferred; the title below is a proposed
display title. Add the actual public repository and demo URLs after publication.
Final-entry live checks passed on the built-in ideal simulator and IonQ Aria-1
noise simulator; their saved counts and execution records were independently audited.

## Title

QPhylo: protein phylogeny with QAOA normalized minimum cut

## Short description (under 240 characters)

Reconstruct an eight-protein tree with Qiskit QAOA, recursive normalized cuts and constraint reinforcement. Compare sampled cuts with exact solving and model-based neighbor joining, with auditable results.

## What it does and why

Phylogenetic trees summarize evolutionary relationships. Recursive graph cuts
offer a different approach from bottom-up methods such as neighbor joining, but
both the affinity construction and cut solver can influence the result.

QPhylo asks whether Fernández-Otero and colleagues' QAOA approach can follow
Onodera-style constraint reinforcement and normalization more closely. This
eight-taxon pilot separates sampled optimization error from differences between
graph-cut and model-based phylogenetic methods.

## How it works

Proteins become a BLASTP affinity graph using self-score normalization and
explicit directed-score averaging. At each node, QAOA solves separate fixed-size
cuts. A wrong-size lowest-energy sample triggers one permitted penalty increase;
accepted feasible candidates are compared by normalized cut on the current
induced graph. Recursion returns a topology.

The same alignment supplies a WAG+Gamma model-distance NJ baseline with 100 site
bootstraps. Exact fixed-size recursion diagnoses the cut solver; it is never a
fallback that replaces a failed measured partition.

## Press Run

This is a Python/Qiskit project. Leave `RUN_MODE = 'auto'` in `main.py`.

- **Built-in simulator (24q):** fresh full-tree QAOA/COBYLA optimization, p=1,
  1,000 shots per circuit. The output shows grouping, Newick, exact/NJ comparisons,
  reinforcement and current job counts.
- **IonQ simulator:** exactly one eight-qubit root-cut measurement using angles
  trained in the saved ideal run. This demonstrates selected-backend circuit
  execution; it does not retrain angles or infer a new full tree.

No supporting code/data uploads or optional algorithm packages are required. The selected
backend is used directly. **SAVED COMPARISON** is clearly labeled historical;
**CURRENT RUN** and `REPORT_JSON` identify new measurements and actual job metadata.
The uploaded `LICENSE` and `CREDITS.md` contain the full terms and attribution;
`main.py` uses a short copyright/SPDX header and runs independently of them.

Both execution paths are verified: full-tree mode completed 403 jobs and matched
the exact NMcut topology. Aria-1 noise single-cut mode completed one job, with
218 feasible and 27 optimal-cut shots/1,000 and an exact best observed cut.
This noisy check samples previously trained angles; it is not a new full tree.

## Findings so far

The synthetic eight-taxon optimizer pilot used five paired initialization seeds.
All ten trees matched exact graph-cut topology. COBYLA used fewer mean circuit
jobs in this fixture; broader optimizer superiority is not established.

The verified protein run completed 403 ideal-simulator jobs and five successful
penalty reinforcements. Its topology matched exact NMcut (unrooted RF=0).
Model NJ matched the provided reference (RF=0), while both NMcut trees differed
from NJ by one internal split (RF=2, normalized 0.2). NJ's taxon1/taxon8 split had
94% bootstrap support. Because the discrepancy also exists with exact NMcut,
it is not attributable to QAOA topology approximation on this run.

## Contribution, limits and next work

Our contribution is an auditable Qiskit adaptation and comparison workflow:
local recursive normalization, bounded reinforcement, protein affinity provenance,
matched-input NJ and backend-aware execution. Preprocessing, retry limits and
optimizer budgets are explicit adaptations rather than a full paper reproduction.

This is one protein input/initialization on ideal simulation. The supplied
reference's generating history is unverified. No biological superiority,
quantum advantage, hardware speedup or large-scale performance is claimed.
Next work would use more seeds and controlled inputs with verified generating
trees, then test circuit depth and noise under matched budgets.

## Reproduce or extend

The accompanying repository contains reusable modules, tests, selected input,
raw reports, tables, figures, hashes and local Conda instructions. Start with
README and `docs/REPRODUCING.md`. `prepare_submission.py` regenerates the one-file
Qollab entry from audited snapshots. Preserve fresh outputs separately and use
the collector to rescore them before making new claims.

## Credit and disclosure

NMcut methodology: Onodera, Hara, Aoki, Asahi and Sawamura,
[2023 paper](https://doi.org/10.1016/j.ympev.2022.107636).
Prior quantum phylogeny work: Fernández-Otero, Pena and Pichel,
[2026 paper](https://doi.org/10.1007/s11227-026-08465-x) and
[PhyloTree_ReConst](https://github.com/Nicofero/PhyloTree_ReConst).
QAOA: Farhi, Goldstone and Gutmann. Data source credited by FO: RAxMLGrove.
Qiskit, SciPy, BLAST+, R, ape, phangorn and Qollab tutorials are credited in full
in `CREDITS.md`, with preserved licenses and data-source notices.

The project began from earlier planning, an initial reference core and the
upstream fork. Hackathon additions and AI-assisted Codex development under Noah's
direction are disclosed in `CREDITS.md`. Onodera's Digital Annealer was
quantum-inspired; our simulator outputs are not quantum-hardware results.
Original project software/documentation is MIT; third-party terms are preserved.

**This effort is supported by Qollab & IonQ**

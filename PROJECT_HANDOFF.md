# Qollab phylogenetic reconstruction project — handoff

Updated: 11 October 2026, submission preparation.

## Start here

This document transfers the relevant context from the research/planning chat to
a new coding chat in the project folder. It is a project snapshot, not a full
transcript. Read it before changing the implementation, then inspect the actual
repository and files: work done after this snapshot takes precedence.

Noah wants to continue development from the project interface, with the source,
decisions, and results together. The immediate objective is a small working
experiment within the limited hackathon time. Implement and test concrete steps;
keep scientific claims tied to evidence. Distinguish completed work, proposals,
and unresolved choices. Explain changes in plain language.

**Current milestone:** gate/ordering checks and the fixed-angle IonQ cut probe
passed; the sampled optimized cut and four-taxon tree have completed on Qollab's
built-in BasicSimulator. Eight-taxon smoke and live COBYLA now also pass.
The user corrected the coordinate seed outputs and completed the five-pair
comparison: all ten trees match the exact unrooted topology, with no failed runs.
See D18 and `results/eight_taxa_comparison/`. The protein FASTA/BLASTP bridge
in D19 is now implemented and verified on the existing eight-protein input.
The first protein Qollab run is verified (D22): 403 DONE jobs, RF0 against the
exact affinity tree and five successful alpha100->600 reinforcements.
Model-based NJ on the same proteins and the three-method comparison are complete
(D23): NJ matches the supplied reference (RF0); QAOA/exact NMcut differ from NJ
by one internal split (RF2, normalized0.2). NJ's1/8 split has94% site-bootstrap
support. See `results/protein1822_comparison/summary.md` and `tree_comparison.png`.
The controlled biological benchmark and reference provenance remain pending.
Protein-derived reconstruction works; biological accuracy and quantum advantage
remain unestablished.

**Submission priority (D24):** User wants almost everything complete tonight
because tomorrow is travel/family time. Deadline: October11,19:00 Eastern /
20:00 Halifax. Official rules/FAQ were read: runnable MIT Qiskit project on
Qollab, disclosed starting work, exact sponsor credit, public2–3 minute demo.
Judges use free local simulation; IonQ execution must fit one job per run.
The final entry is now prepared (D25): `qollab_smoke/submission_main.py` runs
fresh full-tree COBYLA on recognized ideal BasicSimulator and one frozen-angle
protein root cut on IonQ/other backends. Both final live checks are verified
(D28): 403 DONE ideal jobs / RF0 to exact, and one DONE Aria-1 noisy cut job with
218 feasible/27 optimal shots out of 1,000, best observed gap 0.
README, credits, Project Card, reproduction guide, demo outline and curated
`submission_release/` are prepared (D26–D30):131 inventoried files, valid links,
and55 tests from inside the copy (54 pass/one skip). The user explicitly deferred more experiments and
wants video near last. Curated source is now published on GitHub (D34);
hackathon submission and video remain pending.
D31: the user confirms Qollab uploaded files, title, about, README and header
are ready. They accept the remaining platform header rendering; no more header
iterations are needed. Next is publishing the curated release contents as a
public GitHub repository, adding verified links, recording/uploading the demo,
and captain submission. Actual repository/project/video URLs and public access
have not yet been verified.
D32: user authorizes a new public GitHub repository named `QPhylo`. Connector
profile is `noamurphy`, but no repository-creation tool is exposed; browser
access to GitHub was denied. No remote repository/commit was created. Ask the
user to create public QPhylo with an initial README and provide its URL, then
resume publication using the prepared release and existing-repository tools.
D33: user created `https://github.com/noamurphy/QPhylo`. Public flag and initial
main commit `0b3816ec9708137fbb7bea53af54c45fecd70dc8` verified; only a stub
README exists. Connector blob upload returned403 `Resource not accessible by
integration`; repository and installation listings are empty despite GitHub
plugin being installed. No files/commit published. Await user finishing or
reconnecting GitHub access to QPhylo, then upload the131-file inventory plus
manifest, preserve initial commit, and verify all Git blob hashes. Do not use
other credentials/browser access as a workaround. Local release remains valid.
After user-reported reconnection, the write retry still returned403 with empty
repository/installation listings; no files were uploaded. Next ask the user to
inspect the installed GitHub App's repository access and include QPhylo, or
complete app installation if absent. Avoid another reconnect-only loop.
D34 supersedes the above access blocker: after another user access change,
license upload succeeds and QPhylo appears in the authorized repository list.
Source is published at `https://github.com/noamurphy/QPhylo`, commit
`1b1324e58ad85ce9dc74dfd35a3b132867856bdf`. All132 source/evidence/manifest
files match the prepared bytes; initial commit is retained as parent. Follow-up
tracking status is synchronized in the working files, release and repository.
`qollab_upload/README.md` now contains the source/reproduction links; reupload
only that README to Qollab and add the repository URL to any dedicated field.
Next is public Qollab access check, demo recording/upload, then captain submission.

D35: user has applied the source-linked README, reviewed GitHub and connected
it to Qollab. Qollab is intentionally private; the signed-out page returns
"Project not found", consistent with that setting. Public review is deferred
at the user's request. The user reports folders are not viable in Qollab and
inline figures fail with relative paths. The upload README now uses direct
`raw.githubusercontent.com` URLs for the two verified public GitHub figures.
Only reupload that README and check rendering; main/notices/images are unchanged.
This supersedes the asset-folder upload requirement for inline figures below.

**Upload layout (D27/D29):** The user prefers separate notices instead of full license
text embedded in code. `python prepare_submission.py` generates the bundle and
`qollab_upload/{main.py,LICENSE,CREDITS.md}`, plus the maintained
`qollab_upload/README.md` and `assets/`. The README is the finished Project Card
(user reports bidirectional synchronization). Upload the full set at the Qollab project
root, replacing its main.py. The program's short copyright/SPDX header points
to the separate notices; its executable code and embedded graph remain unchanged.
D29 adds original tree header `assets/project_header.png`, full tree and optimizer
figures inside the README, and `prepare_project_header.R` for regeneration from
saved data. The historical field draft is archived; Qollab presentation is now
accepted by the user (D31), while public access is not yet verified.
D30 supersedes the tree-header recommendation after the
user's screenshot showed an approximately 5.3:1 crop losing its labels/context.
Use `assets/project_banner_v2.png`, conceptual AI protein/quantum artwork with
subjects in the center band. Prompt/disclosure are in `docs/HEADER_ART.md` and
credits. The earlier tree header and actual results plots are retained unchanged.

Read [DESIGN_DECISIONS.md](DESIGN_DECISIONS.md) for the running decision record
and reasons. Update it when making or revising choices; retain superseded choices.
[AGENTS.md](AGENTS.md) carries that maintenance instruction for future chats.
[qollab_smoke/README.md](qollab_smoke/README.md) is the current run guide.

### Verified execution history

- The early one-X/two-X tests produced the expected 1,000 measurements but then
  reported `main` undefined. A main-function wrapper did not work for the user.
  Later top-level Bell execution was clean: 00=499, 01=3, 10=6, 11=492. Its exact
  noise selection was not recorded; the earlier runner error remains unexplained.
- Qubit ordering checks on IonQ Aria 1 (Noise Model) returned q0:001=999/000=1
  and q2:100=1,000. RX/RZ/ZZ checks returned their expected strings in all shots.
  These support the adapter's reversed count-string decoding and basic gate plan.
- The four-qubit fixed-angle probe completed cleanly on IonQ Aria 1 noise model,
  Qiskit 2.5.2: 196 feasible and 51 optimal-cut shots out of 1,000. Chosen AB|CD,
  cut=4, gap=0, Ncut=4/11. This was not angle optimization. Saved output:
  `qollab_smoke/qaoa_cut_standaloneResults_noamurphy-qaoa-phylogenetic-tree-reconst-output-completed-2026-10-10T18-13-10.946Z.html`.
- An IonQ optimized run hit "only one job is allowed per code run" on its second
  submission. The checkpoint implementation preserves sampled optimization across
  invocations; it is locally tested, with live IonQ continuation still pending.
- The supplied exact tutorial uses NumPyMinimumEigensolver (classical solving,
  zero IonQ jobs). Its QAOA tutorial uses StatevectorSampler with COBYLA locally
  and also submits zero IonQ jobs. That tutorial demonstrates working SciPy and
  optional algorithm imports in Qollab. Its printed probabilities sum to 0.4617;
  their square roots sum to 1, suggesting a conversion issue. Do not treat its
  percentages as validated; our adapter independently scores raw integer counts.
- The original optimized-cut code completed eight training jobs and one fresh
  final job on BasicSimulator, Qiskit 2.5.2. All nine count dictionaries were
  independently rescored: AB|CD, cut=4, gap=0, Ncut=4/11. Final feasible rate=70.2%,
  optimum shots=305/1,000, final penalized mean=86.166 versus initial 110.16.
  Alpha=100 accepted on first attempt; no reinforcement needed. Output:
  `qollab_smoke/noamurphy-qaoa-phylogenetic-tree-reconst-output-completed-2026-10-10T19-32-54.722Z.html`.
  Its stale user-entered label says IonQ Aria 1. Every job's actual backend metadata
  says BasicSimulator and no noise model. Preserve the original and classify it
  as an ideal local run; fix the manual label on future runs.
- The updated multi_job toy-tree run completed 18 DONE jobs/18,000 shots on
  BasicSimulator with the correct ideal-backend label. Both c=1/c=2 attempts
  accepted alpha=100 on their first attempt. Final feasible rates were 47.1%
  and 46.2%; the candidate Ncuts were 4/3 and 4/11. The selected tree and exact
  tree were both `(('A','B'),('C','D'));`. Every count dictionary, optimizer
  request/decision, attempt and trace was replayed/verified without new jobs.
  Output: `qollab_smoke/qaoa_toyTestResults_noamurphy-qaoa-phylogenetic-tree-reconst-output-completed-2026-10-10T19-58-08.397Z.html`.
  Only the root is split; the pair children terminate directly. Deeper recursive
  optimization was subsequently checked on eight taxa; live penalty
  reinforcement was later verified on proteins (D22).
- Eight-taxon smoke output under `eight_taxa_results/smoke/` completed 72 DONE
  jobs on BasicSimulator, with root plus two child split nodes, exact topology
  and RF=0. All comparison outputs were independently rescored and traces
  validated. COBYLA seeds 42–46 completed 169/220/222/194/214 jobs, SciPy 1.18.0,
  each with RF=0 and no retries. Coordinate files labeled seed42–seed46 all
  record seed=42, simulator_seed=20261052, identical measured parameters/counts
  and 248 jobs. Preserve them; only one is an initialization replicate.
  `results/eight_taxa_seed42_preview/` contains a valid single-pair CSV/PNG/PDF.

### Current implementation and execution choice

- `qphylo_core.py`: exact fixed-size reference, local induced-graph normalized cut,
  recursive topology and QUBO coefficients. Core remains unchanged.
- `qphylo_qaoa.py`: Ising conversion, explicit basic-gate QAOA, backend sampling,
  count validation/scoring, seeded coordinate search, selectable SciPy COBYLA
  and bounded reinforcement. Initial alpha=100, increment=500, one retry.
  Default/smoke budget is nine jobs per attempt; comparison budget is 31,
  including fresh final sampling. Scaling preserves minima but changes gamma units.
  COBYLA uses the algorithm underneath FO's wrapper, with documented differences
  in budget and settings; live execution is verified on Qollab for five seeds.
- `qphylo_run_steps.py`: replay saved measurements and submit only the next missing
  one. No automatic resubmission on errors/timeouts, bit projection or exact fallback.
- `qollab_smoke/cut_pilot.py`: maintainable source. `execution_mode="multi_job"`
  completes the optimizer/recursive calls on local simulators in one invocation.
  `execution_mode="one_job"` checkpoints for restricted IonQ runs. Both use the
  injected selected backend, identical objectives, optimizer and scoring. The
  previous assumption that the limit applied to all Qollab runs is superseded.
- `qollab_smoke/qaoa_cut_standalone.py`: generated paste-ready bundle, default
  fixed_probe and one_job. `prepare_qollab.py` rebuilds it; `--resume <output>`
  restores a checkpoint and explicitly selects one_job. Rebuilding without
  --resume returns to the initial source settings.
- `qphylo_benchmark.py`: synthetic four/eight-taxon fixtures, exact diagnostics,
  unrooted split metrics and resource bounds. Exact scoring never replaces a
  failed quantum split.
- `qphylo_proteins.py`: standard-library FASTA reader, local external BLASTP,
  mean-self-score normalization and explicit bidirectional averaging. Exports
  raw/audit matrices, zero-diagonal core graph, provenance and connectivity.
  For connected <=12-taxon inputs, also exports an exact affinity-tree reference.
  Default input is the existing eight-protein FASTA. `--validate-only` requires
  no BLAST installation. The first protein-derived QAOA run is verified (D22);
  model-based NJ is now verified (D23).
- `qphylo_nj.py` / `build_model_nj.R`: reuse the supplied alignment, validate its
  hash/taxa against QAOA, calculate WAG+Gamma5/shape0.5 distances, NJ and100 site
  bootstraps, then export four trees, splits, unrooted RF, PNG, versions and hashes.
  Run `python qphylo_nj.py` inside qphylo. R/ape/phangorn are installed there.
  Quote decoding and FBP proportion-to-percent conversion are checked; display
  rooting treats support as edge labels. Negative NJ branches are retained.
- `prepare_qollab.py --preset eight_smoke|eight_coordinate|eight_cobyla`:
  separate paste-ready eight-taxon scripts. See `qollab_smoke/EIGHT_TAXA_GUIDE.md`.
- `collect_qaoa_results.py` and `plot_qaoa_results.R`: validate saved raw counts,
  require matching ideal-backend settings, export CSV/manifest and four-panel
  PNG/PDF. `results/toy4_validation/` is a real single-method preview from the
  saved four-taxon result, not an eight-taxon comparison.
- No StatevectorSampler version was created; it was discussed as an alternative.
  The new successful built-in run makes replacing the sampler unnecessary.
  D25 now implements local-trained angles followed by a single selected-backend
  cut measurement in the final entry; Aria-1 noise validation passed (D28). This remains
  separate from sampled-IonQ optimization/checkpoint continuation.

### Next live session

**Current task:** Review the completed first protein QAOA/exact/NJ/reference
comparison and finish submission. Extra experiments are deferred until all
submission materials are ready. The final `submission_main.py` has now passed
built-in and IonQ Aria-1 noise execution and independent output audits (D28).
See `results/submission_validation/summary.md`; original outputs are preserved.
Next apply the Project Card and finish publication/video/submission. Follow
`docs/QOLLAB_FINAL_CHECK.md` for repeats and
`docs/SUBMISSION_CHECKLIST.md`. The final runner embeds historical NJ/results
without claiming they are new measurements. Publish/apply prepared materials
only when authorized; record the2–3 minute video near last. The user's
`qollab_smoke/protein1822_cobyla_run1.html` is verified: 403 DONE ideal jobs,
403,000 shots, four recursive split nodes, nine accepted candidates and five
successful alpha100->600 retries. All candidate cut gaps are zero within
roundoff and the tree matches the exact affinity reference (RF0). Read
`results/protein1822_cobyla_run1/summary.md`; collector exports are adjacent.
Floating-point score audits now allow1e-12 relative/absolute roundoff while
discrete sample fields stay exact. Original output is untouched. Reused shot
seeds make identical-angle final jobs deterministic rather than independent.
One protein seed does not establish reliability or biological accuracy.

The model-NJ comparison is complete in `results/protein1822_comparison/`.
WAG+Gamma5, shape0.5, 100 site bootstraps, seed42; supplied389-column alignment
reused rather than realigned with MUSCLE. All100 trees, distances and supports
are saved. Independent split parsing confirms all RF values and bootstrap counts.
NJ matches the supplied reference and differs from both NMcut trees by one split:
{taxon1,taxon8} instead of {taxon1,taxon6,taxon7}. Support for1/8 is94%; shared
pairs2/3,4/5,6/7 each100%, other internal split89%. Two negative NJ terminal
branch lengths remain in exported trees; figures/RF ignore lengths. This
NMcut/NJ discrepancy exists with exact solving and cannot be attributed to QAOA
topology approximation on this run. The supplied reference is not established
as a generating tree. D23 records adaptations and limits.

For repeat runs, use `qollab_smoke/protein1822_cobyla.py` on Qollab's built-in
ideal simulator. It embeds the verified protein graph and provenance, so only
the standalone code needs to be pasted into `main.py`. It uses multi_job,
COBYLA, p=1, 1,000 shots, seed42, budget31, alpha100/+500 and one retry.
The arbitrary-recursion/retry ceiling remains 930 jobs; actual cost depends on
optimization termination, selected splits and reinforcement. Save the full
protein output separately from synthetic runs. The collector accepts
`protein_tree` reports and preserves dataset provenance; figures label these
as protein-derived. The first protein QAOA result is recorded in D22.
The corrected five-seed comparison is complete (D18), and the protein bridge
has passed a real BLASTP run (D19). BLAST+2.17.0 is installed in qphylo.
Repeat with `python qphylo_proteins.py` inside that environment.
Default output is `results/protein1822/`. `graph.json` contains labels and weights
accepted by the unchanged core. The exact tree is an affinity-graph reference,
not biological truth. Inspect connectivity and scoring provenance before the
next Qollab experiment. Symmetry averaging and zero diagonal are adaptations,
not established copies of every Onodera preprocessing convention (D19).

The original eight-taxon protocol below remains available for repeat experiments;
the smoke test itself and actual COBYLA execution have passed.

The user accepted the first-draft workflow and chose eight taxa. Paste
`qollab_smoke/eight_smoke.py`, select Built-in simulator (24q), and save the
complete output in `qollab_smoke/eight_taxa_results/smoke/`. Expect three split
nodes and unrooted RF=0 against the exact reference. Root cut=16/Ncut=8/19;
each four-taxon child cut=12/Ncut=0.8. Exact-topology/no-retry cost is 72 jobs;
strict worst case is 270. A failed run is evidence, not permission for fallback.

Then run `eight_coordinate.py` and `eight_cobyla.py` for paired seeds 42–46,
changing only the seed assignment. Both allow 30 training plus one fresh-final
evaluation per attempt, with separate matching simulator seeds. Exact-topology
no-retry cost is up to 248 jobs/run; worst case 930. COBYLA requires multi_job;
only coordinate search supports current IonQ checkpoint replay. FO used 300
COBYLA evaluations, so this is a matched-budget optimizer pilot inside our
adapted NMcut, not a full FO reproduction. Preserve failures and all original
reports in `qollab_smoke/eight_taxa_results/comparison/`.

Collect and plot with:

```bash
python collect_qaoa_results.py qollab_smoke/eight_taxa_results/comparison --output results/eight_taxa_comparison --plot
```

D14/D15 and the eight-taxon guide specify settings, expected results, pairing,
data exports and figure interpretation. Protein preprocessing and model-based
NJ for one documented dataset follow this synthetic validation.

Local validation: 55 tests, 54 pass and one real-Qiskit Statevector test skipped
because Qiskit is absent locally. Tests cover math, scoring, job protocol,
bounded optimization, multi-job cut/tree packaging and checkpoint continuation.
Core self-test passes. The eight-taxon packaging and COBYLA budgeting checks
use explicit backend/minimizer test doubles, not measured cloud results. The
collector/figure pipeline was exercised on real four-taxon output and the
eight-taxon paired study. Protein parser/scoring checks use explicit subprocess
doubles; distinguish these from real BLAST execution. Live reinforcement is
verified by D22; IonQ checkpoint continuation remains pending. The NJ test suite
also includes a real R protein-baseline regression when ape/phangorn are present.

Earlier quantum preparation installed no packages. The protein stage installed
BLAST+2.17.0 and its dependencies only in qphylo, with filesystem approval.
The user then authorized R/ape/phangorn in the same environment (D23): R4.5.3,
ape5.8.1, phangorn2.12.1 are installed and exercised. The reviewed installation
added packages without replacing existing ones.
Python3.11.17, JupyterLab4.6.4 and ipykernel7.4.0 remain unchanged. Local Qiskit,
NumPy, SciPy, Aer and IonQ provider remain absent. Use qphylo's Rscript for NJ;
the unrelated system R lacks ape/phangorn. MUSCLE,
MAFFT, IQ-TREE and INDELible were not found on PATH. Browser packages are separate.

## 1. Goal and intended experiment

Noah is building on Fernández-Otero and colleagues' quantum graph-cut
phylogenetic reconstruction work (called FO below), using RAxMLGrove as a source
of benchmark material. The revised primary focus is to adapt FO's Qiskit QAOA
approach to more closely follow Onodera's constraint reinforcement, recursive
normalization, and model-based neighbor-joining comparisons. Increasing protein
sequence divergence remains a possible evaluation axis rather than the defining
goal. Reported performance differences motivate this comparison; they do not
yet establish which implementation differences caused those outcomes.

Working research question:

> Can FO's Qiskit QAOA approach more closely reproduce Onodera's NMcut procedure
> when constraint reinforcement, recursive normalization, and the model-based
> neighbor-joining baseline are aligned with the original method?

`qphylo_core.py` provides fixed-size cut validation, induced-subgraph normalized
cut scoring, QUBO coefficients, and the recursive structure for this work.
Adaptive penalty reinforcement/retries are implemented in the QAOA solver adapter;
they are outside the core. The biological NJ baseline is now implemented (D23).

Immediate agreed sequence: review Noah's Qollab tutorial image, establish an
isolated Conda environment and compatible package versions, execute one small
circuit on Qollab, then connect the fixed-size QAOA cut adapter to the core.

High divergence was already studied by Onodera. Do not claim that testing
divergence itself is new. The useful contribution would be a controlled QAOA
implementation and comparison that separates biological preprocessing,
recursive algorithm behavior, and optimization error.

Long-branch conditions were also discussed as a possible direction. A global
branch-length multiplier tests increasing divergence; a long-branch-attraction
test requires deliberately heterogeneous branches/topology. These are distinct
experiments. Divergence is a possible benchmark for the revised fidelity goal;
an explicit long-branch experiment is an optional extension, not an implemented
feature.

Suggested initial scope, still a proposal:

- Eight taxa, so every fixed-size cut can be checked exactly.
- One known generating topology, four branch-length/divergence settings, and
  three independent sequence replicates per setting as a small first batch.
- QAOA depth p=1 initially, with explicit seeds, shots, and optimizer budget.
- Compare model-based neighbor joining (NJ), exact fixed-size graph-cut
  recursion, and QAOA using the same graph-cut recursion and affinity matrices.
- Establish one successful Qollab/IonQ cut execution before expanding to full
  trees or a batch. Backend latency may determine the feasible experiment size.

Final numerical settings and sequence simulator have not been chosen.

## 2. Platform constraints and workflow

The planning conversation identified Qollab and Qiskit as hackathon
requirements, with IonQ simulator execution as the intended backend. These are
project constraints carried forward from that discussion; inspect the actual
hackathon rules and supplied tutorial before making submission-compliance claims.

Keep Qiskit responsible for constructing/executing QAOA circuits. Run the final
experiment notebook through the required Qollab environment/backend. The
classical reference and biological preprocessing are supporting Python modules.
Codex is the development assistant: it can edit modules/notebooks, verify math,
debug the backend adapter, automate comparisons, and prepare figures and a demo.
Its involvement does not replace the required quantum framework/platform.

FO has separate quantum annealing (QA) and QAOA paths. Its QA ecosystem includes
D-Wave Ocean components; that path is not the proposed hackathon execution route.
Onodera used Fujitsu's Digital Annealer, a **quantum-inspired** optimizer. Do not
describe those experiments as gate-based quantum hardware executions.

The prior notes recorded this attribution phrase:

> This effort is supported by Qollab & IonQ.

Confirm its required wording and placement against the current submission
instructions. No credentials, private API keys, or authenticated backend setup
have been supplied to the coding session.

## 3. Sources and files to bring into the folder

### Upstream projects

- FO code: https://github.com/Nicofero/PhyloTree_ReConst
- RAxMLGrove: https://github.com/angtft/RAxMLGrove

The user fork is cloned at `PhyloTree_ReConst`, remote
`git@github.com:noamurphy/PhyloTree_ReConst.git`, branch `main`, checked-out commit
`2fba9f1904a13fa41246af9cef25574f367f0239`. The existing local modification to
`qaoa_phylo_tree.ipynb` changes notebook metadata only. No fork files were edited
during this coding-session preparation, and no commits or pushes were made.

Keep the upstream code as a reference/reproducibility source. Preserve its MIT
license notices wherever its code is reused. RAxMLGrove's database/data licensing
is separate; the previous review identified ODbL/DBCL terms. Check those terms
before redistributing database-derived material.

### Original paper and companion paper supplied by Noah

1. **Onodera et al., "Phylogenetic tree reconstruction via graph cut presented
   using a quantum-inspired computer."** Molecular Phylogenetics and Evolution
   178 (2023), 107636. DOI: https://doi.org/10.1016/j.ympev.2022.107636
   Available online in October 2022; journal volume is dated 2023.
2. **Onodera et al., "Data on the solution and processing time reached when
   constructing a phylogenetic tree using a quantum-inspired computer."**
   Data in Brief 47 (2023), 108970. This is the companion data article, not the
   original methodological paper. It was initially supplied before the original.

Available PDF filenames, which should be copied/attached if the next chat needs
to inspect the papers rather than rely on this summary:

- `Phylogenetic tree reconstruction via graph cut presented using a quantum-inspired computer.pdf`
- `Data on the solution and processing time reached when constructing a phylogenetic tree using a quantum-inspired computer.pdf`

The FO paper was discussed and its repository inspected previously. Its full
bibliographic information is not reproduced here; recover and verify it from
the repository/paper before writing formal references.

### A small upstream smoke-test dataset

- https://github.com/Nicofero/PhyloTree_ReConst/blob/main/alignments/8/assembled_sequences_1822.fasta
- https://github.com/Nicofero/PhyloTree_ReConst/blob/main/trees/8/tree_best_1822.newick

These were identified as an eight-taxon protein/tree pair suitable for checking
file parsing and pipeline execution. They have now been processed by BLASTP,
QAOA/exact NMcut and model-based NJ (D19/D22/D23). This pair alone is not the
controlled divergence study.

### Current implementation

Copy `qphylo_core.py` alongside this handoff. The earlier chat supplied it as a
separate downloadable file. The handoff does not embed the Python module or PDFs.

## 4. What the original algorithm actually does

The original method is NMcutDA (Normalized-Minimum cut by Digital Annealer);
FO's quantum adaptation is referred to as QNMcut in our discussion.

For the graph currently being split, with vertex set S and nonnegative
similarities w_ij:

1. For each target size c=1,...,floor(|S|/2), solve a **minimum cut with exactly
   c selected vertices**.
2. Evaluate normalized cut for the resulting candidate for each c.
3. Choose the candidate with the smallest normalized-cut score.
4. Recurse on the two induced subgraphs, stopping at groups of one or two taxa.

For S=A union B:

```text
cut(A,B) = sum of w_ij for i in A, j in B
vol_S(A) = sum over i in A of sum over j in S of w_ij
Ncut_S(A,B) = cut(A,B) * (1/vol_S(A) + 1/vol_S(B))
```

The fixed-cardinality QUBO is:

```text
E(x) = sum over i<j of w_ij * (x_i-x_j)^2
       + alpha * (sum_i x_i-c)^2
```

For binary x, the expanded coefficients are:

```text
constant = alpha*c^2
linear[i] = sum_j w_ij + alpha*(1-2*c)
quadratic[i,j] = 2*(alpha-w_ij), for i<j
```

This architecture is **not** direct global minimization of normalized cut.
The exact reference solves the same fixed-size cut subproblems exactly, then
compares their Ncut values. It must not silently substitute a different algorithm
that enumerates all partitions directly to minimize Ncut.

Minimum-cut ties can produce multiple valid candidates with different Ncut.
The starter uses a deterministic first subset in combinations order. Record
that tie convention; alternative tie treatment is a separate implementation
choice, not a demonstrated reproduction of Onodera's optimizer behavior.

## 5. Findings from the prior FO/Onodera audit

These findings summarize the earlier paper/code comparison. Recheck the exact
upstream commit before modifying its code or publishing claims about it. The
repository changed after publication, so current-code findings cannot by
themselves establish how the published experiments were run.

| Area | Prior finding | Action for this project |
| --- | --- | --- |
| Recursive cut architecture | FO retains the fixed-size min-cut candidates followed by Ncut selection and recursion. | Keep that architecture. |
| Affinity preprocessing | Onodera uses BLASTP local-alignment bit scores, BLOSUM62, E-value cutoff 10, normalized by mean self-scores times 100; no qualifying hit gives zero. FO uses global Biopython raw alignment scores, BLOSUM62, gap-open -10 and extension -0.5, with similar self-score normalization. Global raw scores can be negative. | Prefer an original-style BLAST affinity pipeline for the main fidelity experiment. Treat an FO-style pipeline as a labeled alternative. Do not quietly clip negative weights. |
| Recursive normalization | In the inspected FO path, cuts are computed on the current subgraph but Ncut denominators use rows from the original full matrix. External taxa can change the chosen split. | Compute numerator and both volumes from the same current induced graph. |
| Feasibility | A nonempty partition does not establish that its size equals c. Fixed penalties can leave invalid solutions. | Validate exact cardinality. Report failure/retry explicitly. |
| Scoring quantum outputs | A penalized QUBO expectation, or an Ising eigenvalue without its offset, is not the actual unpenalized cut of a returned partition. An inspected fallback risks mixing these. | Decode a measured bitstring, validate it, then rescore its actual cut and local Ncut classically. Preserve offsets when checking energy equivalence. |
| NJ baseline | Onodera uses model-based evolutionary distances from an MSA. The inspected FO baseline uses 100 minus similarity. | Use model-based NJ as the main biological baseline; label an FO-style NJ separately if retained. |
| Imports | FO's QAOA module imports shared QA utilities, bringing D-Wave-related dependencies into the QAOA path. | Separate graph/biology utilities from backend-specific code. |
| Benchmark provenance | RAxMLGrove's inferred empirical trees can seed sequence simulations; the repository does not provide the original empirical alignments as the assumed benchmark source. | Track empirical template, branch modifications, sequence simulator, and known generating tree separately. |

One inspected current QAOA implementation, `qaoa_phylo_tree_qiskit`, used a
penalty of 500*n, a caller-supplied depth, and COBYLA with 300 iterations. Earlier
paper notes referred to a different fixed penalty. Onodera described an adaptive
penalty/retry procedure. Do not rely on a large fixed alpha as proof of
feasibility. The prior audit recorded QAOA-file commit
`1703b781e8ffaa3f930a1a4100f49c42ddd9ddd1` dated June 2026; verify it against
the actual checkout and pin the version used in our experiment.

Interpretation: FO preserves the central recursive-cut idea, but preprocessing,
normalization, optimization safeguards, and baseline construction can materially
change behavior. The audit supports testing those differences, not a blanket
claim that the published work is invalid.

## 6. What is implemented and verified

`qphylo_core.py` is a new, standard-library-only implementation. Use Python
3.10 or later. It does not import Qiskit, NumPy, BioPython, or D-Wave.

Its graph convention is symmetric, finite, nonnegative weights and a zero
diagonal. Zero diagonal is an explicit starter convention; verify the original
self-loop/diagonal handling before claiming exact biological reproduction.
Zero-volume groups raise an error. A policy for isolated vertices or missing
BLAST connectivity is still needed before processing real affinity matrices.

Implemented functions:

- `validate_weights` and `validate_bits`: validate the graph and solver outputs.
- `cut_weight` and `normalized_cut`: score actual partitions.
- `qubo_energy` and `qubo_coefficients`: include the constant offset.
- `exact_fixed_size_cut`: exhaustive minimum cut for a specified c.
- `choose_split`: compare fixed-size candidates by current-graph Ncut.
- `reconstruct`: recurse, return Newick topology and a split trace.
- `self_test`: mathematical and regression checks.

Solver adapter contract:

```python
solve(weights, c) -> sequence_of_binary_values
```

Return one bit per current-graph vertex, with exactly c ones. Supply the adapter
to `choose_split` or `reconstruct`. The core scores its returned bits itself.
It does not silently replace a failed quantum solver with the exact solver.
The core trace records taxa, groups, target size, cut, and Ncut. The prepared
adapter adds per-evaluation weights, parameters, penalty, seed, feasible rate,
counts, energy expectation, and backend/job metadata, plus attempt outcomes.
The pilot prints a JSON report; protein dataset/provenance logging is implemented.

Verified commands:

```bash
python3 qphylo_core.py --self-test
python3 qphylo_core.py
```

Observed check result:

```text
PASS: local normalized cut, QUBO expansion, feasibility, recursive scoring.
```

The four-taxon demo returns `(('A','B'),('C','D'));` for a graph with high
within-pair similarities and low between-pair similarities.

The checks include a six-vertex regression graph whose induced four-vertex
subgraph chooses a different split when outside degrees are incorrectly included.
They also check QUBO expansion on every four-bit assignment, complement cut
symmetry, invalid cardinality rejection, recursive local scoring, and invalid
affinity rejection. These are mathematical checks, not biological validation.

The returned Newick root represents the recursion. No evolutionary root or
branch lengths are inferred. Evaluate topology with unrooted metrics.

## 7. Immediate implementation sequence

This was the initial implementation outline. Steps 1–5 and the model-NJ/metric
parts of step7 are complete for the pilot. The next choice is repeat protein
seeds or a controlled benchmark with verified generating-tree provenance.
DESIGN_DECISIONS.md records the current acceptance/reinforcement policy and
supersedes any abbreviated wording in this original outline.

1. **Inspect the project folder and Git state.** Locate this handoff, the core,
   papers, upstream clone/fork, and any existing Qollab notebook. Record the
   upstream commit rather than assuming current main matches the paper.
2. **Obtain the working Qollab tutorial/backend cell.** Check installed Qiskit
   and IonQ integration versions. Use platform-managed credentials/environment
   variables; never commit keys. Verify one tiny supported circuit execution.
3. **Implement a Qiskit fixed-size cut adapter.** Build the QUBO/Hamiltonian,
   check bit order and energy/offset against the exact core, run shallow QAOA,
   and select/validate feasible measured samples. If none are feasible, surface
   the failure and apply a documented retry policy. Log all retries/resources.
4. **Run one cut on the required backend.** Compare its measured cut with the
   exact minimum for the same c and graph. Only then attempt an eight-taxon tree.
5. **Implement biological input/preprocessing.** Parse matched FASTA/Newick
   labels. Establish BLAST-style affinity handling, symmetrization if needed,
   missing hits, diagonals, and zero-volume behavior. These choices must be
   documented and shared by exact-cut and QAOA runs.
6. **Create the small controlled divergence dataset.** Choose an available
   protein simulator (Pyvolve or INDELible were discussed), generating model,
   sequence length, branch settings, and fixed seeds. Original-style notes
   mentioned WAG, five gamma categories, 1% invariant sites, 500 residues, and
   indels; verify details and decide which can be reproduced in time.
7. **Add model-based NJ and comparison/logging.** Use the same simulated
   sequences and label mapping. Save unrooted topology metrics and optimization
   diagnostics. Keep every method's provenance and execution settings explicit.
8. **Prepare a Qollab notebook and honest demo.** Show the working pipeline,
   exact-versus-QAOA results, divergence response if enough runs finish, platform
   execution evidence, limitations, and required attribution.

The main comparisons answer different questions:

- Exact-cut versus QAOA, same affinities: optimization/quantum approximation.
- Exact-cut versus NJ: biological reconstruction behavior of the methods.
- Original-style versus FO-style preprocessing/scoring: algorithmic differences.

Do not change all these factors at once and attribute the outcome solely to QAOA.

## 8. Reporting and open decisions

Recommended biological logging (dataset-specific fields not implemented): dataset identifier, generating tree,
simulation model/seed/divergence setting, method, software/backend versions,
depth, shots, optimizer budget, penalty/retries, feasible sample rate, cut
optimality gap, Ncut, reconstructed tree, unrooted RF distance, and wall time.
Synthetic optimizer/trace/resource logging is implemented; normalized unrooted
RF divides by 2(n-3). Define biological-reference provenance separately.
Separate backend queue time from
execution/optimization time where the platform exposes it.

Still open:

- Local quantum-package setup and additional Qollab package versions (fork and
  working branch are confirmed above; cloud Qiskit reports 2.5.2).
- Qollab notebook/backend interface, allowed budgets, and final submission rules.
- Simulator/model, exact divergence settings, and feasible replicate count.
- BLAST diagonal, symmetry, disconnected graph, and tie-handling conventions.
- Revisit penalty, sample selection and retry choices in D04/D14 if evidence warrants.
- NJ alignment/tool-version fidelity beyond the implemented D23 baseline.
- Whether time permits a separate long-branch experiment or preprocessing ablation.

Smaller Ncut does not guarantee a more accurate evolutionary tree. QAOA sample
concentration is not biological branch support. Original work already studied
high-divergence behavior and graph perturbations; distinguish optimizer stability
from evolutionary confidence. No speedup, improved accuracy, or quantum advantage
has yet been demonstrated here.

## 9. Keeping this handoff useful

Place this file and `qphylo_core.py` at the repository root. Ask the new chat to
read this file, inspect the working tree, and continue with the current milestone.
Put lasting coding conventions into `AGENTS.md` if useful; this handoff is
ordinary project documentation and need not replace existing agent guidance.

Update the milestone, verified commands, decisions, and next step as work
progresses. Keep experimental outputs in reproducible files rather than relying
on conversation history. Preserve attribution when incorporating upstream code.

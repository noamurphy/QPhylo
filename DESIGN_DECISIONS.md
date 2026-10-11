# Design decisions

Recorded: 10 October 2026. This is the running record of what we chose, why,
what evidence informed it, and what remains undecided. PROJECT_HANDOFF.md is
the broader project snapshot; qollab_smoke/README.md contains run instructions.

Add an entry whenever a choice changes the scientific method, implementation,
execution workflow, or experiment scope. Keep earlier decisions and mark them
superseded when evidence changes them. Separate implemented choices from
proposals and backend measurements from local mathematical checks.

## D01 — Focus on methodological alignment

**Status:** Active; replaces divergence as the defining research objective.

**Choice:** Adapt FO's Qiskit QAOA approach to more closely follow Onodera's
constraint reinforcement, recursive normalization, and model-based NJ baseline.

**Why:** Differences in the NMcut implementations may matter to the reported
performance difference. An aligned comparison is more informative than changing
sequence divergence alone. Those differences have not been shown to cause the
published results. Divergence and long-branch experiments remain possible later
evaluation axes. Onodera used a quantum-inspired Digital Annealer, not gate-based
quantum hardware or a quantum annealing QPU.

**Evidence:** User's revised objective; paper and FO source review summarized in
[the handoff](PROJECT_HANDOFF.md).

## D02 — Keep a small exact reference and a replaceable cut solver

**Status:** Implemented.

**Choice:** The core calls `solve(weights, c)` and requires exactly c selected
vertices. An exact solver provides a reference; QAOA implements the same interface.

**Why:** This separates errors in recursion/scoring from errors in approximate
optimization and lets both methods receive the same graph. Exact comparison is
tractable on the small pilot, not a scalability claim or a hidden QAOA fallback.

**Implementation:** [qphylo_core.py](qphylo_core.py),
[qphylo_qaoa.py](qphylo_qaoa.py).

## D03 — Score each recursive split on its current graph

**Status:** Implemented for the graph core; biological fidelity remains pending.

**Choice:** Evaluate candidates for c=1..floor(n/2), compare normalized cut on
the current induced graph, then recurse. Use finite symmetric nonnegative
weights and a zero diagonal. Stop directly on one or two taxa.

**Why:** Local normalization follows the intended recursive NMcut structure;
using denominators from the entire original graph can change split selection.
Symmetry and diagonal conventions make the objective unambiguous. These graph
conventions do not settle how raw biological similarities should be constructed.
The result is a topology; branch lengths and biological rooting are not inferred.

## D04 — Enforce size through a penalized objective and reinforcement

**Status:** Implemented with a bounded pilot budget.

**Choice:** Minimize `cut + alpha*(sum(bits)-c)^2`. In a fresh final sample,
reinforce alpha if the lowest observed penalized-energy sample is infeasible,
even when another sample has the correct size. Among accepted feasible samples,
select the minimum actual cut. Fail explicitly when the retry budget is exhausted.

**Why:** This preserves the intended reinforcement trigger rather than silently
repairing or replacing an infeasible answer. It mirrors an aspect of Onodera's
method, not Digital Annealer hardware behavior. Exact solving, projecting bits,
and merging complementary samples are not substitutes for the QAOA result.

**Pilot:** alpha=100, increment=500, at most one retry. These are recorded pilot
settings, not established optimal settings for biological graphs.

## D05 — Explicit Hamiltonian scaling and circuit conventions

**Status:** Implemented and checked mathematically; relevant gates tested on Qollab.

**Choice:** Map binary variables using x=(1-Z)/2. Divide phases and optimizer
objective by the largest absolute nonconstant Ising coefficient (minimum 1).
Keep the identity offset when scoring energy; omit its global circuit phase.
Use H initialization, a transverse-X mixer, RZ fields and CX-RZ-CX interactions.

**Why:** Scaling makes pilot angle steps usable with large penalties while
preserving each fixed-penalty minimizer. It changes gamma units and must be
recorded. Explicit basic gates avoid depending on opaque RZZ support.

**Convention:** Parameters are all gamma values followed by all beta values.
Reverse Qiskit measurement strings so core tuple position i means vertex/qubit i.
Asymmetric q0/q2 tests, RX/RZ/ZZ checks, and independent energy checks support
these choices; symmetric Bell counts alone cannot establish bit ordering.

## D06 — Start with a bounded coordinate-search optimizer

**Status:** Original pilot choice retained as the default; deferral of COBYLA
preparation is superseded by D14. Live COBYLA is now verified (D16).

**Choice:** Seeded coordinate search, p=1, 1,000 shots per evaluation, at most
nine evaluations per penalty attempt including a fresh final measurement.

**Why:** Initially the browser's SciPy and optional algorithm packages were
unverified. A small dependency-light optimizer made execution and scoring
reviewable first. FO used COBYLA, so this is not an optimizer-matched reproduction.
The supplied QAOA tutorial now demonstrates SciPy, QAOA and COBYLA imports in
Qollab. Package availability is no longer the sole reason to retain this choice;
keeping the pilot stable isolates execution/recursion checks before comparing
optimizers. The initialization seed does not seed backend shots or noise.

## D07 — Isolate local dependencies and keep the fork as reference

**Status:** Active.

**Choice:** Work locally in the qphylo Conda environment; install dependencies
only as needed. Keep the new core/adapter separate from the cloned FO repository.

**Why:** Browser packages and local Python packages are distinct environments.
Isolation avoids affecting unrelated projects. Separate modules make deviations
from FO inspectable. Preserve upstream license notices if upstream code is reused.
No local quantum packages were installed during these preparation steps.

## D08 — Use top-level Qollab scripts with an injected backend

**Status:** Active; the attempted main-function workaround was abandoned.

**Choice:** Python/Qiskit project shape, paste-ready bundled script, top-level
execution, and the `backend` supplied by Run Experiment. No provider credentials
or guessed `main()` contract.

**Why:** This follows the working tutorial. Early X tests finished their jobs
but then reported `main` undefined. Wrapping the code did not work for the user;
a later Bell run was clean in top-level style. The original runner error's cause
remains unknown. Bundling avoids relying on untested multi-file imports.

## D09 — Checkpoint IonQ evaluations between code runs

**Status:** Implemented; scope narrowed to restricted backends by D10.

**Choice:** Save measured counts and requests as explicit JSON; replay previous
optimizer decisions and submit only one missing measurement per code run.

**Why:** An IonQ optimizer run produced: "only one job is allowed per code run".
The original advice to run the whole optimizer in one invocation was incompatible
with that execution path. Checkpoints preserve the sampled optimizer, fresh final
measurement and reinforcement without assuming files/globals survive in Qollab.
Errors/timeouts stop work rather than triggering an automatic new cloud job.

**Implementation:** [qphylo_run_steps.py](qphylo_run_steps.py). Rebuild using
`prepare_qollab.py --resume <saved-output>`; no old jobs are resubmitted.
Continuation has local test coverage; actual IonQ continuation remains pending.

## D10 — Use unrestricted local execution for development

**Status:** Active; supersedes treating the one-job limit as applying to all Qollab runs.

**Choice:** Keep the selected-backend adapter and expose `execution_mode`:
`multi_job` completes optimization/recursion on local simulators; `one_job`
checkpoints for the observed IonQ restriction. Default to `one_job` until the
user explicitly selects the unrestricted path. The graph method and optimizer
remain identical across these execution modes.

**Why:** The user's original optimizer completed nine jobs in a single built-in
simulator run. We can proceed without replacing the sampler. A StatevectorSampler
option is unnecessary for this milestone and would bypass Run Experiment's
selected noise backend. The tutorial's local StatevectorSampler is another valid
ideal-simulation path, not evidence of multi-job IonQ access. A future local-angle
optimization followed by one IonQ measurement is a distinct experiment from
tuning angles with noisy IonQ measurements; do not conflate them.

**Evidence:** [Built-in optimized-cut output](qollab_smoke/noamurphy-qaoa-phylogenetic-tree-reconst-output-completed-2026-10-10T19-32-54.722Z.html).
Nine DONE jobs, BasicSimulator, Qiskit 2.5.2, eight training evaluations plus a
fresh final sample. Independently rescored every count dictionary. Selected
AB|CD, cut=4, exact gap=0, local Ncut=4/11. Final feasible rate=70.2%, optimal
cut shots=305/1,000, mean penalized energy=86.166 versus 110.16 initially.
One penalty attempt, alpha=100, no reinforcement needed. This is one toy ideal
run, not biological evidence or demonstrated quantum advantage.

**Metadata discrepancy:** The user-entered label says "IonQ Aria 1 (Noise
Model)", while every job identifies BasicSimulator and no noise model. Preserve
the original output; classify this run as local ideal and correct the label on
future runs. A manually entered label does not select or prove the backend.

## D11 — Keep ideal, noise-model and hardware results distinct

**Status:** Active; fuller comparisons pending.

**Choice:** Develop on the local ideal simulator; later compare selected noise
models and, if appropriate, real hardware as separately labeled experiments.

**Why:** Local simulations allow rapid multi-evaluation checks. IonQ noise models
represent different devices from the local IBM models. Simulators remain classical
computations and their noise models do not guarantee physical-device performance.
Hardware access and any required hackathon backend remain separate choices.

## D12 — Verify a toy cut before recursion and biological comparisons

**Status:** Cut, four-taxon tree and deeper eight-taxon recursion passed (D16).

**Choice:** Gate/ordering checks, fixed-angle cut probe, optimized cut, then the
four-taxon recursive tree. Build biological preprocessing and NJ after this path
is verified. Use exact calculations as diagnostics, never to replace a failed
approximate split.

**Why:** Each stage isolates a different source of error. The fixed probe proved
execution/scoring; its optimum among many samples alone did not prove optimization.
The built-in runs verify optimization execution and the root split through the
recursive interface. Biological performance requires matched data and independent
baselines. Save full reports/counts rather than only the best answer.

**New evidence:** [Four-taxon toy-tree output](qollab_smoke/qaoa_toyTestResults_noamurphy-qaoa-phylogenetic-tree-reconst-output-completed-2026-10-10T19-58-08.397Z.html).
BasicSimulator, correct ideal-backend label, multi_job mode, Qiskit 2.5.2.
All 18 jobs DONE, 1,000 shots each. Both c=1/c=2 attempts accepted alpha=100
without reinforcement. Final feasible rates: 47.1% for c=1, 46.2% for c=2.
Candidate c=1 cut=11/Ncut=4/3; c=2 cut=4/Ncut=4/11. QAOA and exact trees both
`(('A','B'),('C','D'));`. Independently rescored every evaluation and replayed
all optimizer requests/decisions and attempts from saved counts without new jobs.
The trace contains only one split: the children are pairs and need no further
optimization. This does not yet validate multi-level recursive splitting on Qollab
or live penalty reinforcement; those have local test coverage.

## D13 — First working phylogenetic draft workflow

**Recorded:** 10 October 2026.
**Status:** Workflow accepted by the user; implementation incomplete. Eight-taxon
synthetic validation is specified in D14; biological dataset selection is pending.

**Choice:** Use a single small six-to-eight-taxon protein dataset for an initial
reproducible comparison. First exercise deeper recursion and add a selectable
COBYLA optimizer, retaining the existing coordinate-search result as a reference.
Then connect original-style BLASTP affinities to exact and QAOA reconstruction,
and implement model-based NJ from an MSA of the same sequences. Export trees,
split/optimizer logs, topology comparisons and runtime/resource counts.

**Why:** The quantum execution and four-taxon integration now work. Child-graph
splitting is the next validation gap; restoring FO's optimizer addresses the
documented optimizer mismatch. Biological preprocessing and NJ are essential
to the revised goal. One small dataset keeps integration reviewable before
expanding to a divergence study or noise/hardware comparisons.

**Draft completion criterion:** One documented dataset produces QAOA and exact
NMcut trees from the same affinity matrix, a model-based NJ tree from the same
taxa, and a saved comparison report. Failures remain explicit. Establish dataset
and reference-tree provenance before claiming biological accuracy; an empirical
inferred reference is not a known generating tree. Sequence simulation, larger
replicate studies and hardware benchmarking are separate later milestones.

**Revisit:** Select the biological data source and toolchain before fixing an
implementation schedule. This acceptance does not select a protein dataset,
generating model or final biological benchmark settings.

## D14 — Eight-taxon recursion, then matched coordinate/COBYLA runs

**Recorded:** 10 October 2026.
**Status:** Implemented; live eight-taxon smoke and five COBYLA seeds completed.
Five-seed paired comparison is incomplete: four coordinate seeds need reruns
as documented in D16.

**Choice:** The user chose eight taxa as the next validation step. Use a
synthetic hierarchical graph: affinities 9 within AB/CD/EF/GH, 3 between pairs
in each four-taxon half and 1 across ABCD|EFGH. Run `eight_smoke.py` first on
the built-in ideal simulator with nine evaluations per penalty attempt. Then
use `eight_coordinate.py` and `eight_cobyla.py` with matching settings and
seeds 42–46: at most 30 training evaluations plus one fresh final sample per
attempt. Preserve alpha=100, increment=500, one retry, p=1, 1,000 shots,
Hamiltonian scaling and the selected-backend sampler. Seed local shot sampling
separately with 20261010 + initialization seed.

**Why:** The graph requires root and child optimization, unlike the completed
four-taxon test. Equal evaluation budgets and initial angles isolate optimizer
behavior inside our adapted NMcut. Five paired initialization seeds give a
small descriptive pilot; they are not biological replicates. Exact-topology
no-retry runs need 72 smoke jobs or up to 248 comparison jobs. Strict arbitrary
recursion/retry ceilings are 270 and 930, respectively.

**Implementation:** COBYLA calls SciPy directly, the algorithm underneath FO's
Qiskit wrapper. Record SciPy version, termination and measured angles. Use
rhobeg=1.0, tol=1e-4, clip trial angles into the same gamma/beta box as coordinate
search, and validate the best measured training parameters with fresh shots.
An explicit guard prevents training submissions beyond budget. Budget exhaustion
is distinct from cut infeasibility. COBYLA currently requires multi_job;
coordinate-search IonQ checkpoints remain supported.

**Evidence:** User request selecting eight taxa and a COBYLA comparison;
[run guide](qollab_smoke/EIGHT_TAXA_GUIDE.md), mathematical fixture checks and
explicit backend/minimizer test doubles. FO's inspected `qaoa_functions.py`
uses COBYLA(maxiter=300). [SciPy documentation](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-cobyla.html)
defines its evaluation budget and notes the implementation change in 1.16.

**Limits / revisit:** This uses FO's optimizer algorithm, not all FO settings or
the full FO pipeline. Our 30-evaluation pilot, initialization, bounds, scaling
and penalties differ. At preparation, actual COBYLA optimization needed its
Qollab run; D16 records its subsequent verification. Local SciPy/Qiskit were
not installed. The hand-built graph is not evolutionary data.
This prepared comparison is a concrete initial protocol, revisable with evidence.

## D15 — Validate saved counts and export a reproducible descriptive figure

**Recorded:** 10 October 2026.
**Status:** Implemented; checked using the actual four-taxon output and the valid
eight-taxon seed-42 pair. Full five-seed comparative results remain pending.

**Choice:** Save full HTML outputs, including failures. The collector independently
rescores raw counts, validates traces, checks paired settings/backend/version
metadata and exports runs, evaluations and final-cut CSV files with a manifest.
Reject duplicate optimizer/seed entries and mixed smoke/comparison budgets.
Plot root c=n/2, attempt-0 training curves, fresh-final optimal-shot frequencies,
unrooted RF distance to the exact affinity-tree reference and total job counts.
Retain child-graph diagnostics in the tables. Use already-installed base R for
PNG/PDF; no local quantum or plotting packages were installed.

**Why:** Comparing scores across different child graphs or penalties can obscure
what improved. The shared root problem supports a direct optimizer comparison.
RF measures agreement with the exact graph-cut tree, not evolutionary accuracy.
Raw counts and original outputs make percentages and failures independently
reviewable. Exact calculations diagnose results after sampling; they never
substitute for failed QAOA splits.

**Evidence:** `collect_qaoa_results.py`, `plot_qaoa_results.R`, report validation
tests and [real four-taxon preview](results/toy4_validation/optimizer_comparison.png).
The preview contains coordinate search only, not a fabricated eight-taxon or
COBYLA result.

**Limits / revisit:** Uniform-shot references are analytic, not measured random
solver runs. Recorded time sums backend job elapsed times, including queue time,
not whole-program wall time. Failures retain missing topology metrics and must
be counted. Five seeds on one graph support description, not significance,
biological accuracy, speedup or quantum advantage. NJ remains a later milestone.

## D16 — Preserve mislabeled seed outputs and rerun missing initializations

**Recorded:** 10 October 2026.
**Status:** Diagnosed and recovery scripts prepared; coordinate seeds 43–46 pending.
Retention and replacement-folder choice superseded by D17.

**Evidence:** The collector rejected duplicate coordinate/seed 42. All five
`comparison/coordinate_seed42.html` through `coordinate_seed46.html` record
config.seed=42 and simulator_seed=20261052, with identical parameters and counts
in every evaluation. Filenames cannot establish which seed ran. The COBYLA
reports record distinct seeds 42–46, SciPy 1.18.0 and 169/220/222/194/214 jobs.
Each report and the 72-job smoke output was independently rescored and its trace
validated. All recover the exact topology with three split nodes, RF=0 and no
reinforcement. Actual COBYLA and deeper recursion are therefore now verified on
the built-in ideal simulator. The exact reason run settings stayed at 42 is unknown.

**Choice / why:** Keep originals and the duplicate guard. Count coordinate seed
42 once; do not relabel saved report seeds or treat identical seeded measurements
as independent initialization replicates. Prepare separate paste-ready files for
coordinate seeds 43–46 in `qollab_smoke/seed_runs/`, and save replacement runs in
`eight_taxa_results/reruns/`. Reuse the five valid COBYLA runs. Add `--inspect`
to show stored seeds/settings, identify all conflicting files in errors, and
allow `prepare_qollab.py --output` to retain separate scripts per seed.

**Available result at this stage:** `results/eight_taxa_seed42_preview/optimizer_comparison.png`
(the preview was subsequently superseded/removed after D17/D18)
and accompanying CSV/manifest use only the genuine pair. Both trees match;
coordinate search uses 248 jobs, COBYLA 169. Root c=4 fresh-final optimum frequency
is 9/1,000 for both. This is one descriptive pair, not the completed five-seed
study or evidence of superiority. Full source inventory is in
`results/eight_taxa_run_audit/stored_run_settings.csv`.

**Next:** Run the four prepared coordinate scripts, verify the printed/stored
seed, and collect explicitly selected original and replacement outputs using
the recovery command in the eight-taxon guide. Live reinforcement and biological
input/NJ remain pending.

## D17 — Replace repeated seed outputs rather than archive them

**Recorded:** 10 October 2026.
**Status:** User-selected workflow; corrected coordinate runs pending.

**Choice:** The user will replace the four mislabeled coordinate seed43–seed46
outputs directly in `eight_taxa_results/comparison/`. Do not require archiving
the repeated seed-42 outputs or a separate reruns folder. Retain valid coordinate
seed 42 and the five valid COBYLA outputs. Use the ordinary folder collector
command after replacement; the duplicate guard and stored-seed checks remain.

**Why / evidence:** The user explicitly declined retaining these repeats and
suspects the wrong seed variable was edited. The repeats have identical measured
angles/counts and do not supply additional initialization replicates. Their
diagnosis is already recorded in D16 and the audit. This supersedes D16's
retention advice for these four files; it does not change how genuine failed
runs or other distinct experimental results are counted.

## D18 — Complete the paired synthetic optimizer pilot

**Recorded:** 10 October 2026.
**Status:** Completed after the user replaced the four mislabeled outputs.

**Evidence:** `results/eight_taxa_comparison/manifest.json` records ten runs,
five paired seeds 42–46, zero failures, 2,290 measured evaluations and 81 final
samples. Both optimizers recover the exact unrooted topology in 5/5 runs.
Mean tree jobs are COBYLA 203.8 and coordinate search 254.2; balanced-root
optimal-shot frequency averages 1.34% for both. Balanced-root feasibility
averages 36.8% and 59.32%, respectively. Coordinate seed 43 misses the exact
balanced cut but reaches the same unrooted topology through different splits.

**Interpretation / next:** Keep COBYLA available as the FO-aligned optimizer
and coordinate search as a comparison. Early termination contributes to fewer
COBYLA jobs; this does not establish better optimization or biological accuracy.
Uniform bitstrings have an analytic 99.96% probability of including an optimal
balanced root cut among 1,000 shots, so tree recovery alone is weak evidence of
optimizer benefit. The user chose protein preprocessing as the next integration
step; measured uniform baselines and harder graph studies remain useful later.

## D19 — Audited protein FASTA and BLASTP affinity bridge

**Recorded:** 10 October 2026.
**Status:** Implemented, tested and run successfully on the eight-protein input.

**Choice:** `qphylo_proteins.py` reads first-token FASTA labels in input order,
rejects duplicate/empty IDs and invalid protein records, uppercases residues
and removes alignment gaps '-' and '.'. Preserve full descriptions and gap
counts. Use safe internal BLAST identifiers with a saved label map. Initial
integration input is FO's eight-taxon `assembled_sequences_1822.fasta` (389
residues each), not a newly generated accuracy benchmark.

**Scoring:** Run local BLASTP with BLOSUM62, E-value <=10, gap costs 11/1,
SEG off, composition adjustment 2, one best HSP per directed pair, all n
targets allowed and one thread. Require positive self bit scores. Normalize
each directed score as `100*b_ij/mean(b_ii,b_jj)`, then average both directions
to obtain symmetry; a missing direction contributes zero. Keep raw bit scores,
directed normalized scores and symmetric similarities for audit. No clipping
above 100 and no summed HSP scores.

**Why / limits:** The paper specifies BLASTP bit scores and self-score
normalization rather than FO's global raw scores. Bidirectional averaging is
an explicit adaptation: the papers inspected do not settle that convention.
Retain normalized self similarities of 100 in the audit matrix, then set the
graph diagonal to zero for the existing core's convention. Self-edge handling
affects Ncut volumes and should be examined before claiming complete fidelity.
BLAST version and explicit settings are recorded; the paper used BLAST+2.12.0,
so a newer executable is not an exact version reproduction.

**Exports / checks:** Export graph JSON, labeled matrices, raw hits, input
hashes, settings/version and connectivity diagnostics. Disconnected graphs are
saved for review without artificial edges or automatic component reconstruction.
Connected graphs of <=12 taxa also get an exact fixed-size-cut reference tree
and trace; this is not ground truth or a QAOA result. No NJ or biological
reference accuracy is computed yet. Standard-library tests check parsing,
normalization, HSP choice, missing hits, failure reporting and core integration;
the source FASTA has been validated. Imported hits have unverified invocation
settings and must use the exported input/IDs with a supplied BLAST version.

**Sources:** Onodera paper section 2.1 and companion paper section 3.3;
[NCBI BLASTP options](https://www.ncbi.nlm.nih.gov/sites/books/NBK279684/table/appendices.T.blastp_application_options/),
[NCBI output/HSP options](https://www.ncbi.nlm.nih.gov/sites/books/NBK279684/).
**Live evidence:** BLAST+2.17.0 and its dependencies were installed in the
existing qphylo Conda environment with filesystem approval; Python3.11.17 and
JupyterLab4.6.4 remain unchanged. The first real run in `results/protein1822/`
has 64 qualifying directed hits, no missing cross hits, no directional score
asymmetry, 28/28 positive edges and one connected component. No duplicate
sequences. The exported exact affinity tree is
`(('taxon8',('taxon1',('taxon6','taxon7'))),(('taxon2','taxon3'),('taxon4','taxon5')));`.
Every exported graph weight was independently recomputed from raw BLAST scores;
input/hit hashes and the tree export match the manifest. This verifies the
protein-to-graph-to-exact-core bridge, not biological accuracy. Forty local
tests pass; one real-Qiskit test remains skipped. No local Qiskit was installed.

## D20 — Freeze the protein graph into a standalone Qollab pilot

**Recorded:** 10 October 2026.
**Status:** Prepared and locally tested; live protein QAOA run pending.

**Choice:** Generate `qollab_smoke/protein1822_cobyla.py` from the existing
core/solver and audited `results/protein1822/graph.json`. Embed the weights,
taxon order, graph hash and preprocessing manifest. Require a connected,
graph-ready export whose labels, source/hit hashes and weights agree with its
adjacent manifest and affinity matrix. No BLAST, FASTA or local module imports
are needed inside Qollab. Use top-level injected-backend execution on the
built-in ideal simulator, multi_job, COBYLA, p=1, 1,000 shots, seed42,
simulator_seed=20261010+seed, budget31 and alpha100/+500 with one retry.

**Why:** Single-file bundling continues the verified Qollab workflow while
isolating protein preprocessing from browser execution. Keeping pilot settings
makes the input graph the principal change from the synthetic COBYLA study;
COBYLA is selected for continuity with FO, not proven superiority in D18.
The new `protein_tree` report mode distinguishes real protein affinities from
synthetic fixtures. The existing collector/figures accept it, preserve source
provenance and label protein-derived runs correctly.

**Evidence / limits:** Forty-two local tests pass, one real-Qiskit test is
skipped. New integration coverage executes the embedded eight-taxon graph with
explicit circuit/backend/minimizer doubles and verifies recursive scoring,
reference RF and collector provenance. Generation rejects altered matrices,
hash mismatches and exports requiring review. These are integration checks,
not measured QAOA outcomes. Exact affinity-tree diagnostics never replace a
failed sampled split. Retain failure reports; reinforcement may be needed with
these larger weights. Protein COBYLA remains multi_job only, and NJ/biological
reference comparisons are not implemented by this step.

## D21 — Retain direct script replacement for the remaining pilot

**Recorded:** 10 October 2026.
**Status:** User-selected workflow; launcher experiment abandoned.

**Choice / why:** Continue pasting each standalone experiment into Qollab's
`main.py`. The user withdrew the swappable-file launcher request to prioritize
the previously verified execution method within the remaining time. The local
launcher passed basic execution checks but was not tested on Qollab; its changes
were reverted. D20's protein bundle and settings are unchanged. The user will
run it on the built-in simulator and save the full output; live results remain
pending.

## D22 — Verify the first protein QAOA run and allow numerical roundoff in audits

**Recorded:** 10 October 2026.
**Status:** First live protein reconstruction verified; model-based NJ remains next.

**Evidence:** `qollab_smoke/protein1822_cobyla_run1.html` completes 403 DONE
BasicSimulator jobs / 403,000 shots, Qiskit2.5.2, SciPy1.18.0, seed42. Its four
recursive splits exactly match the exported affinity-tree reference (RF0).
Nine fixed-size problems are accepted, each with zero cut gap within roundoff.
Five require the one permitted reinforcement from alpha100 to600: root c1/c2/c3,
the three-taxon child c1 and the 2/3/4/5 child c1. All then accept. Exhaustive
QUBO checks show these five alpha100 objectives have infeasible global minima.
This verifies live reinforcement, previously unexercised on the synthetic graphs.
Root c4 wins Ncut1.044440 over c3 Ncut1.047917; its optimum appears31/1000 times.
The accepted root c2 optimum appears only3/1000 times, so one successful seed
does not establish reliability. No biological accuracy or quantum advantage claim.

**Audit choice / why:** Independent local rescoring differs from the browser
by at most2.73e-12 on calculated energies/cuts. Replace literal equality for
floating-point summary fields with relative/absolute tolerances1e-12/1e-12;
keep discrete sample fields exact. Reject meaningful numeric changes and changed
shot counts. Preserve original outputs. All403 scores, trace scores, DONE records,
best-measured final angles, retries and graph/provenance match after this change.
Forty-two tests pass, one real-Qiskit test skipped; checks cover tolerance and
continued rejection of corrupted scores/counts. Derived tables and a readable
summary are saved in `results/protein1822_cobyla_run1/`.

**Limits / next:** Nine optimizer attempts reached the training budget; that is
not a circuit failure. Reusing the simulator seed for each job means fresh-final
execution is not statistically independent sampling at identical angles. Keep
that limitation explicit. Prioritize the agreed model-based NJ integration and
three-method/reference comparison for the first draft; more QAOA seeds can follow.

## D23 — Model-based NJ from the same supplied protein alignment

**Recorded:** 10 October 2026.
**Status:** Implemented and verified on protein1822; three-method comparison complete.

**Choice:** Use local R/phangorn `dist.ml(model="WAG", k=5, shape=0.5,
exclude="pairwise")`, then ape `nj`. Use WAG model frequencies, no separate
invariant-site component, 100 site-bootstrap replicates, seed42, randomized
taxon order per replicate, single-process execution. Match Onodera section2.8's
distance model, gamma parameters, NJ algorithm and bootstrap count. Reuse the
supplied FO eight-protein alignment (389 columns, canonical residues, no gaps)
rather than perform a new MUSCLE3.8.31 alignment. Record this adaptation and
current package versions; do not claim full version/alignment reproduction.

**Environment:** The user explicitly chose installation of R/ape/phangorn
into qphylo. A proposed separate qphylo-nj environment was not created. Inspect
Conda's plan with existing packages frozen before installation; preserve the
current Python3.11.17, JupyterLab4.6.4 and ipykernel7.4.0 versions.
The reviewed plan added92 packages without replacing any existing packages.
Installed R4.5.3, ape5.8.1 and phangorn2.12.1; all three protected Python/Jupyter
versions were rechecked unchanged after installation.

**Implementation / validation:** `qphylo_nj.py` checks alignment columns,
taxon order and source hash against the verified protein QAOA report, then calls
`build_model_nj.R`. Export distances, plain/support-annotated NJ trees, all
bootstrap trees, exact/QAOA/reference trees, canonical splits, unrooted RF,
versions, logs, manifest and a topology figure. Reject invalid model distances
or mismatched taxa without repairing them. Negative NJ branches, if any, are
retained and counted. Preserve source files; failed reruns clear stale owned
outputs and keep a failure manifest. The live integration identified two format
conventions: ape retains our Newick quote delimiters, which are decoded before
exact taxon validation, and phangorn2.12.1 returns FBP proportions, which are
converted explicitly to percentages. Access bootstrap trees through multiPhylo
`[[` to restore compressed tip labels. Reroot support as edge labels for display.
Forty-seven tests pass, one real-Qiskit test skipped, including a real R regression
on the protein data. These local classical checks submit no quantum jobs.

**Measured evidence:** `results/protein1822_comparison/` contains completed
NJ distances, all100 bootstrap trees, support, RF tables, hashes and a visually
checked PNG. All distances are finite, symmetric and nonnegative with zero
diagonal; off-diagonal range0.029124–3.622905 substitutions/site. Independently
parsed Newick splits confirm all six RF comparisons and the Python QAOA trace.
QAOA/exact RF0; NJ/provided-reference RF0; either NMcut method versus NJ/reference
RF2, normalized0.2. Four of five unrooted internal splits agree. The differing
NMcut split is {taxon1,taxon6,taxon7}; NJ/reference use {taxon1,taxon8} instead.
Independently counting saved bootstrap splits confirms NJ support: 100% for
2/3,4/5,6/7 pairs; 89% for the2/3/4/5-versus1/6/7/8 split; 94% for1/8.
Two negative NJ terminal branch lengths are retained and counted; the figure
and RF analysis compare topology only. Thus the observed NMcut/NJ difference
is present even with exact graph solving, rather than QAOA topology error here.
NJ+bootstrap took about0.3 seconds locally; this is not a matched runtime study.

**Scope / sources:** The provided `tree_best_1822.newick` is a labeled external
reference; generating-tree provenance remains unverified. The display anchor
is not a biological outgroup. Bootstrap percentages measure site-resampling
support, not the probability of correctness. Source: Onodera paper section2.8,
[phangorn distance documentation](https://klausvigo.github.io/phangorn/reference/dist.hamming.html),
[bootstrap documentation](https://klausvigo.github.io/phangorn/reference/bootstrap.pml.html).

## D24 — Verified submission requirements and tonight's completion priorities

**Recorded:** 10 October 2026, following the user's deadline clarification.
**Status:** Rules verified; submission preparation and final backend check pending.

**Constraint / proposed scope:** The user will travel/visit family tomorrow and
wants almost everything completed tonight. Prioritize a runnable final entry
point, existing validated results, attribution, repository packaging, Project
Card and demo. Additional protein seeds or controlled benchmarks are optional
and should not delay required submission materials. This narrows the earlier
suggested larger evaluation; no new accuracy or quantum-advantage claim follows.

**Verified requirements:** Official terms give October11,19:00 Eastern
(20:00 Halifax) as the deadline. Publish a runnable Qiskit project on Qollab,
release our code/documentation and external project repo under MIT, disclose
starting material and contributions, and include the exact credit line:
"This effort is supported by Qollab & IonQ". Video must be2–3 minutes and
publicly viewable without sign-in. Projects remain available on Qollab12 months.
Judges use free local simulation; IonQ results do not change scoring. The FAQ
confirms an IonQ run must fit one job and should be tested before submission.
The captain submits through the registration page and can edit until deadline.

**Implementation consequence:** The verified full protein COBYLA workflow
uses multiple jobs and will exceed IonQ's per-run limit. Preserve that full
local-simulator workflow. Prepare and clearly label a one-job demonstration
within the final entry point; this is not completed IonQ tree optimization.
Its final live run still needs the user's Qollab execution. Existing IonQ
gate/cut evidence alone does not verify the final submitted code.

**Evidence / limits:** Read the official terms and expanded FAQ in the browser:
[terms](https://qollab.xyz/programs/hackathon/terms),
[submission/video/job-limit FAQ](https://qollab.xyz/programs/hackathon/how-it-works#faq).
No submission, publication, cloud job or external message was performed.
Retain third-party notices and distinguish supplied references from known
generating trees. Reference provenance remains unverified (D23).

## D25 — One self-contained submission entry, selected-backend modes

**Recorded:** 10 October 2026.
**Status:** Implemented and locally checked; final live Qollab execution pending.

**Choice:** `prepare_submission.py` builds `qollab_smoke/submission_main.py`
from the maintained core, QAOA, benchmark and `submission_entry.py` sources.
`RUN_MODE='auto'` runs fresh full-tree COBYLA only on the recognized ideal
BasicSimulator. IonQ/noisy/unknown backends receive one eight-qubit balanced
root-cut circuit. `single_cut` is an explicit override; unsupported forced
`full_tree` fails before jobs. Preserve the D21 direct-paste workflow and the
selected injected backend; no StatevectorSampler substitution or main wrapper.

**Why:** The full local reconstruction is verified and useful for judging;
IonQ's one-job limit prevents completing its optimizer loop in one invocation.
Reuse the D22 accepted root-c4 angles for a clearly labeled cut demonstration,
rather than claim a new full tree on that restricted path. Existing coordinate
checkpoint code remains separate; this does not replace it or claim IonQ
COBYLA checkpoint optimization.

**Evidence:** Freeze parameters `[1.1552531532226773,1.1696075648591386]`,
alpha100, Ising scale35.00307314074985, from the verified protein report's
call3/attempt0 final measurement. The builder revalidates training counts,
accepted attempt, best training point, graph/provenance and NJ output hashes.
The previous source HTML SHA-256 is
`c007843de0901f051dae13ae3c10937d8ca0281fab4a993298156a447268c242`.
Full reconstruction keeps the original budgets/reinforcement/scoring/audit
contract. Seven new tests use explicit backend/circuit/minimizer doubles to
check routing, current/historical job separation, collector compatibility,
preflight rejection, failure/no-feasible-sample behavior and altered snapshot
rejection. Full local suite:55 tests,54 pass,one real-Qiskit skip.

**Consequences / limits:** Saved tree/NJ/bootstrap comparisons are labeled
historical; only current `jobs`/`counts` constitute new backend execution.
The single-cut mode performs no optimization, retry, projection, exact fallback
or new tree. Exact cut diagnostics only score the actual measurement. Noisy
sampling can return no feasible bitstring or a nonoptimal cut; a completed
probe signifies execution, not tree-solving success. Fresh Qollab checks on
built-in and IonQ are required; earlier small IonQ checks do not establish
this final eight-qubit circuit's compatibility. Conservative routing means
other unrestricted simulators also receive one cut until explicitly verified.

## D26 — Finish documentation and a curated release before extra experiments

**Recorded:** 10 October 2026.
**Status:** Local documentation/release prepared; application/publication/video pending.

**Choice:** The user accepted D24 priorities, explicitly deferred extra
experiments until submission is ready, and left video recording near last.
While away, the user authorized finishing other documentation if waiting for
approval or live execution. Prepare README, copy-ready Project Card, reproduction
guide, final checks, submission checklist, MIT license, credits and a timed
2:45 demo outline. Package reusable modules/tests, selected inputs, audited
results and original output evidence into a local curated release; exclude
the full fork, Git internals, research PDFs, scratch outputs and environments.

**Why / evidence:** Required deliverables and existing measured results have
priority under the user's travel deadline. The selected input/reference remain
unchanged, and research/starting-code attribution distinguishes method credit,
code origin, data source and the hackathon's additions. The supplied FO PDF
verified its authors and DOI; local upstream MIT license was preserved. The
RAxMLGrove repository explicitly states ODbL1.0/Database Contents License1.0
and gives Höhler et al.(2022), DOI10.1093/bioinformatics/btab863. Record those
source notices in CREDITS rather than relicense upstream data as original MIT
software. Conda recipe records the existing local versions; no dependencies
were installed during this step. Release files are copied with a hash inventory.

**Verification:** The curated `submission_release/` contains115 inventoried
files with matching hashes and valid local Markdown links/Python syntax. Its
own test run completed55 tests:54 pass,one real-Qiskit skip. Obsolete recovery
links were corrected without erasing D16 history; the current guide identifies
the completed D18 comparison. No source measurements were regenerated or edited.

**Consequences / limits:** Local artifacts are ready for review and copying;
they do not mean the Qollab Project Card was applied, GitHub/Qollab published,
video recorded/uploaded or captain submission completed. Exact selected-data
origin/generating history remains unverified despite the paper's general
RAxMLGrove source statement. A fresh-machine Conda solve is untested. Additional
experiments can be appended later without changing present claims. See
`docs/SUBMISSION_CHECKLIST.md` for remaining external actions.

**Sources:** [RAxMLGrove notice](https://github.com/angtft/RAxMLGrove#readme),
[FO repository](https://github.com/Nicofero/PhyloTree_ReConst), user priority
and away-time instructions; [CREDITS.md](CREDITS.md) preserves the notices.

## D27 — Separate license/credits files in the Qollab upload

**Recorded:** 10 October 2026.
**Status:** Implemented at the user's request; final live execution still pending.

**Choice:** Replace the full commented MIT texts in `submission_main.py` with
a short copyright/SPDX header pointing to `LICENSE` and `CREDITS.md`. Upload
three files at the Qollab project root. `python prepare_submission.py` now
regenerates the bundled source and prepares `qollab_upload/main.py`, `LICENSE`
and `CREDITS.md`; release packaging also includes this upload set.

**Why:** The earlier embedded notices kept all attribution with a standalone
copy. The user confirmed Qollab can accept separate files and preferred the
usual separate-license layout. This supersedes the full-notice header choice
in the initial D26 packaging. The copyright, complete MIT license, upstream
MIT notice, method/data credits and source-specific terms remain preserved.

**Evidence / limits:** Generated executable syntax is compared against the
previous bundle; only comments change. Upload copies are byte-checked against
the generated code and original license/credits. Code still contains all graph
data and imports no local Python modules or documentation files at runtime.
This changes packaging, not the scientific algorithm, settings, backend routing
or measured results. The user's upload and final live runs remain pending.

## D28 — Final auto-mode live checks verified on both simulator paths

**Recorded:** 10 October 2026.
**Status:** Both final-entry execution paths independently verified; publication/
Project Card/video/submission remain to finish.

**Choice / why:** Accept the user's two saved final-entry runs as satisfying
the outstanding D25/D27 compatibility checks. Keep the tested executable code
and scientific settings unchanged. Mark live execution complete in submission
documentation and include these original outputs plus a separate audit in the
curated release. Treat the ideal repeat as a compatibility rerun, not another
independent protein initialization; classify the cloud run by actual noise metadata.

**Evidence:** `qollab_smoke/submission_builtin_run1.html` selected `full_tree`
from `auto`, completed 403 DONE BasicSimulator jobs / 403,000 shots, four split
nodes, nine accepted cuts and five alpha100->600 reinforcements. Its tree has
unrooted RF0 to exact NMcut and RF2 to the saved model NJ/provided reference.
All stored count summaries, angle bounds, seeded initializations, budgets,
Hamiltonian scales, final best-training-point selections, job/evaluation metadata
and penalty decisions were checked. Replaying the accepted measured partitions
through the core reproduced the reported tree without exact fallback. New job
IDs are disjoint from the original pilot; parameters/counts match that same-seed
run, consistent with deterministic ideal simulation.

`submission_ionq_run1.html` selected `single_cut`, status `probe_completed`:
exactly one DONE job / 1,000 shots on IonQSimulatorBackend, actual
`noise_model='aria-1'`. Independent rescoring confirms 218 feasible shots (21.8%),
27 minimum-cut shots (2.7%), best cut795.4323749543944 / gap0, local
Ncut1.04444015003035, partition {1,6,7,8}|{2,3,4,5}. Frozen angles/source and
gate counts match the prepared bundle. There are no current optimizer attempts,
evaluation history or inferred tree on this path. No post-run errors were found.
Both reports retain an unchanged, explicitly historical comparison snapshot.

**Provenance / limits:** [Audit summary](results/submission_validation/summary.md),
[audit JSON/source hashes](results/submission_validation/audit.json) and
[rescored cloud counts](results/submission_validation/ionq_counts.csv).
Built-in HTML SHA-256:
`823375d4ebde4fa3f78ddcb6e37dcd5780f0ef44773e2d0b4517717b43a2b52f`;
IonQ HTML SHA-256:
`f65f4152e8f20094485dd265e4e325de72cff0484e59a7020056464c4616b932`.
Qiskit2.5.2 on both; SciPy1.18.0 for local-simulator COBYLA. One noisy fixed-angle
cut is not full-tree optimization, a controlled noise comparison or hardware
performance evidence. Prior biological/provenance/advantage limits remain.
IonQ coordinate-checkpoint continuation is still a separate unverified feature.

## D29 — Use the upload README as Project Card, with original results visuals

**Recorded:** 11 October 2026.
**Status:** Implemented locally; Qollab rendering/upload remains pending.
**Choice:** `qollab_upload/README.md` is the maintained public Project Card,
replacing the D26 field-copy draft. Upload README and `assets/` beside the three
D27 code/notice files. Keep the root README as the local/GitHub source guide.
Preserve the former draft in `docs/archive/PROJECT_CARD_D26.md`. Recommend the
original two-tree `assets/project_header.png` for the separate project header;
show full tree and optimizer comparisons inside Markdown. Build custom upload
folders with the same documents/assets and include everything in the release.

**Why:** The user reports README/Project Card synchronization and requests
finished Markdown with images. A simpler header makes the project's central
finding readable at the top without the density of the four-panel results plot.
Original data-derived artwork represents the work and avoids a stock-image or
AI-image rights question. Figures are original project plots, under the existing
project terms, with upstream data/method credit retained.

**Evidence:** User request; saved Newick trees and `splits.csv` in
`results/protein1822_comparison/`; RF0 QAOA/exact, RF2 QAOA/NJ, 94% NJ 1/8
bootstrap support. `prepare_project_header.R` asserts those saved values,
highlights the different unrooted split in each tree, and copies full results
PNGs unchanged. The header was visually checked for readable labels and margins.
Default/custom upload builds and a fresh 129-file release were verified, all
local Markdown links resolve, and both copied results images match the originals.

**Consequences / limits:** No new optimization, NJ/bootstrap run, hardware result
or biological claim. Layout shows topology only, with a display anchor rather
than an inferred evolutionary root. Original raw results/figures are unchanged.
The main.py SHA-256 remains the D27/D28 hash. Relative image rendering and
header cropping must be checked in Qollab after upload; add verified repository
and video URLs when available. The builder maintains the upload README/assets
as source artifacts, rather than regenerating editorial text from a template.

## D30 — Use decorative protein/quantum art for the shallow project header

**Recorded:** 11 October 2026.
**Status:** Alternative generated and prepared locally; final selection/crop in
Qollab remains pending. Supersedes D29's header recommendation only.
**Choice:** Recommend `qollab_upload/assets/project_banner_v2.png` as the header:
AI-generated folded-protein forms and a horizontal branching/quantum motif on
a dark background, with the recognizable subjects kept near the vertical center.
Keep the full results diagrams inside the README and retain the original header.

**Why / evidence:** The user's screenshot shows a roughly 5.3:1 header crop,
which hides the topology image's titles, captions and parts of its taxa. A
decorative banner conveys the topic without depending on scientific labels.
The built-in image tool generated a new asset without a stock-image input;
the prompt is preserved in `docs/HEADER_ART.md`. The generated image was visually
inspected; no empirical claim relies on it.

**Consequences / limits:** This is conceptual art, not a computed protein
structure, tree or circuit specification. Disclose image-generation assistance
in credits and the public README. Do not describe the image as results or claim
exclusive rights. No new experiment or change to tested main.py. The user still
needs to check the actual crop in Qollab. The upload set and curated release
include the alternative; D29's original topology plot remains intact.

## D31 — Accept the Qollab presentation and proceed to publication/video

**Recorded:** 11 October 2026.
**Status:** Presentation complete by user report; publication/submission pending.
**Choice:** Accept the current header rendering and stop image iterations. Treat
Qollab files, title, about, README and header as applied and reviewed. Prioritize
public GitHub publication of the curated release, public-access/link checks,
the planned 2–3 minute demo, then captain submission by the deadline. Additional
experiments remain deferred under D24 until the entry is ready.

**Why / evidence:** The user explicitly says Qollab looks ready and accepts
the remaining platform rendering. The runnable entry and both backend paths
were already verified (D28). A further art iteration is unnecessary for their
accepted presentation and would compete with the remaining submission work.

**Consequences / limits:** This is the user's presentation review, not an
independent platform inspection. No verified public project/repository/video URL
or submission receipt is available yet. Only local status/continuity documents
change; no Qollab code, README, credits, license or image reupload is required.
Refresh matching release copies and preserve original scientific evidence.
Exclude incidental `.Rhistory` session files from the curated inventory, as with
editor checkpoints; keep those files in the working directory.

## D32 — Publish a new public QPhylo repository; creation awaiting user

**Recorded:** 11 October 2026.
**Status:** Authorized by user; repository creation blocked by browser access.
**Choice:** Publish the curated release as a new public GitHub repository named
`QPhylo`. The authenticated connector profile is `noamurphy`. Do not publish the
entire working directory. The user must create the empty repository (with an
initial README commit) and provide its URL before connector publication resumes.

**Why / evidence:** The user confirms there is no existing repository and asks
for `QPhylo`. The connector supports existing-repository operations but has no
repository-creation tool. Opening GitHub's new-repository page was rejected
because browser permission was denied. No alternate browser or raw credential
workaround was attempted. The prepared release remains ready and hash-verified.

**Consequences / limits:** No repository, remote commit or public link was
created. An initial README allows subsequent publication against an existing
commit using connector tools. Preserve the supplied upstream licenses, selected
data notices and provenance limits. Project/video/public-access checks and
captain submission remain pending. These local tracking changes do not require
Qollab file reuploads.

## D33 — Existing public repository verified; upload awaiting GitHub access

**Recorded:** 11 October 2026.
**Status:** Repository exists; content publication blocked by integration access.
**Choice:** Publish the existing curated 131-file inventory plus release manifest
to `https://github.com/noamurphy/QPhylo`, preserving the user's initial commit.
Request completion/reconnection of the GitHub app's access to this repository
before retrying writes. Do not broaden access or obtain alternate credentials.

**Evidence:** User supplied the repository URL. GitHub metadata confirms
`private=false`, default branch `main`; head
`0b3816ec9708137fbb7bea53af54c45fecd70dc8` contains only an eight-byte README.
The first blob upload returned HTTP403, `Resource not accessible by integration`.
Repository listing for `noamurphy` and managed installation listing are empty;
plugin discovery confirms GitHub is installed. These support a missing resource
authorization/installation as the likely cause, not a malformed release.
All inventoried local hashes were verified. No project file or commit was uploaded.

**Consequences / limits:** The public repository is currently a stub, not the
published project source. Do not claim publication, add a misleading source
link to the Project Card, or retry through another credential/browser path.
User must finish repository access for the installed GitHub connection; request
only QPhylo access. Once access changes, retry the authorized write and verify
the committed tree's paths/content hashes before marking publication complete.
Only local tracking files change; no Qollab reupload is needed yet.

**Follow-up:** The user reported reconnecting. A fresh write retry still
returned the same403; repository and installation listings remained empty,
and `main` still pointed to the initial README commit. Reconnection alone did
not fix the access failure. Next inspect the connection's installed GitHub
App repository access and add QPhylo if absent. If no matching app is installed,
the installation/authorization flow still needs completion. Do not keep
retrying writes without another access change. GitHub's official guidance:
[review/modify installed apps](https://docs.github.com/en/apps/using-github-apps/reviewing-and-modifying-installed-github-apps).

## D34 — Publish the curated release through the authorized GitHub connection

**Recorded:** 11 October 2026.
**Status:** Write access verified; publication in progress, pending tree audit.
**Choice:** Publish the 131 inventoried release files plus `RELEASE_MANIFEST.json`
to the user's public `noamurphy/QPhylo` repository. Preserve the initial commit
as parent, use content-addressed Git blobs/tree, then update `main` only after
checking every uploaded path and blob hash. Add the actual source/reproduction
links to the Qollab upload README and refresh matching release copies.

**Why / evidence:** After the user's further access change, repository listing
now includes QPhylo and the license-blob write succeeds. Its Git blob SHA is
`d031c9fb2968e8faf5ac32e0d0003a89ad0e6859`, matching the prepared bytes. This
resolves the D33 access failure without changing credentials or broadening
permissions through agent actions. The parent is the user's initial README
commit, `0b3816ec9708137fbb7bea53af54c45fecd70dc8`.

**Consequences / limits:** Publish only selected source, notices and evidence;
the full fork, research PDFs, environments and incidental session files remain
excluded. No change to tested main.py or scientific results. Only the Qollab
README needs replacement to include its repository link. Public Qollab/video
access and captain submission remain pending. Record successful publication
and verification after the remote branch audit; do not count object uploads
alone as completed publication.

## Open choices

- Protein-derived Qollab reconstruction is verified for one seed (D22); assess
  repeatability later without claiming a full reproduction of FO's pipeline.
- Model-based NJ with shared taxa is verified (D23); verify supplied-reference
  provenance and revisit affinity symmetry/diagonal conventions before claiming
  complete original-method fidelity or biological accuracy.
- Choose biological benchmark sizes, replicates, generating trees, sequence
  simulator and divergence conditions. D14 specifies synthetic validation only.
- Verify local IBM noise-backend compilation before relying on that path.
- IonQ coordinate checkpoint continuation remains unverified. D25 implements
  local-trained angles for a separate one-job cut demonstration; its final
  Aria-1 noise-model validation passed (D28).
- Qollab presentation is applied and accepted by the user (D31). Public GitHub
  publication, public access checks, video and submission remain pending;
  final entry compatibility passed (D28).
- The supplied QAOA tutorial prints probabilities summing to 0.4617. Their square
  roots sum to 1, suggesting a conversion issue. Use raw integer counts in our
  adapter; do not import those percentages as validated probabilities.

## Template for future entries

### Dxx — Concrete decision title

**Recorded:** Date.
**Status:** Proposed / active / implemented / superseded by Dxx.
**Choice:** What we will do.
**Why:** Problem, tradeoff, and reason for this choice.
**Evidence:** Link to source, result, or user decision; distinguish inference.
**Consequences / revisit:** Limits and what evidence would change the choice.

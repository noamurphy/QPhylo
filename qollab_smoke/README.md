# Ready for the next Qollab session

**Submission presentation (D31):** The user confirms files, title, about, README
and header are uploaded and ready. Remaining header rendering is accepted.
Next actions are GitHub publication, public-access checks, demo video and captain
submission; see the [finish-line checklist](../docs/SUBMISSION_CHECKLIST.md).
The instructions below are retained for rebuilding or repeating experiments.
GitHub publication is authorized for a new public `QPhylo` repository (D32).
The user created [QPhylo](https://github.com/noamurphy/QPhylo), but connector
upload is blocked by repository authorization (D33). The public repository
currently contains only its initial README; the curated release is ready.

**Final entry point (D25/D27/D29):** Upload the contents of `../qollab_upload/`
to the Qollab project root: [main.py](../qollab_upload/main.py),
[LICENSE](../qollab_upload/LICENSE), [CREDITS.md](../qollab_upload/CREDITS.md),
[README.md](../qollab_upload/README.md), and its `assets/` folder. The upload
README is the Project Card; select `assets/project_banner_v2.png` for the separate
header image (D30). This is conceptual AI artwork designed for the shallow crop,
not a results diagram. The original results images are copied unchanged, and
the earlier topology header remains available (D29).
The code matches generated [submission_main.py](submission_main.py), leaving
`RUN_MODE = 'auto'`. The recognized ideal built-in simulator performs a fresh
full protein-tree optimization; IonQ/other backends submit exactly one root-cut
circuit using previously trained ideal-run angles. The graph, provenance and
historical comparison are embedded, with saved versus current results explicitly
labeled. No supporting code/data uploads or StatevectorSampler are needed.
Both final live checks are now verified (D28): built-in 403 DONE jobs / RF0 to exact;
IonQ Aria-1 noise model 1 DONE job, 218 feasible and 27 optimal shots/1,000, best
cut gap 0. See the [audit summary](../results/submission_validation/summary.md).
For repeat runs, follow the
[two-run check](../docs/QOLLAB_FINAL_CHECK.md), saving outputs as
`submission_builtin_run1.html` and `submission_ionq_run1.html`.
Rebuild the code and project upload set with `python prepare_submission.py`;
edit `submission_entry.py` for lasting
runner changes. The earlier scripts below remain experimental/reproduction tools.

Rules are verified in D24. Public demo video must be 2–3 minutes; deadline is
October11,20:00 Halifax (19:00 Eastern). The finished Project Card is `../qollab_upload/README.md` (D29); credits and
reproduction materials are also prepared (D26).

**Current step:** The first protein run passed: 403 DONE jobs, RF0 against
the exact affinity tree and five successful penalty reinforcements (D22).
See the [readable summary](../results/protein1822_cobyla_run1/summary.md).
Model-based NJ and the three-method comparison are now complete (D23).
See the [tree comparison](../results/protein1822_comparison/tree_comparison.png)
and [NJ summary](../results/protein1822_comparison/summary.md): QAOA/exact RF0,
NJ/provided-reference RF0, and NMcut versus NJ RF2. The reference's generating
provenance is still unverified. Repeat the classical comparison locally with
`python qphylo_nj.py` from the project root after activating `qphylo`.

For repeat runs, paste the whole [protein1822_cobyla.py](protein1822_cobyla.py)
into Qollab's `main.py` and select **Built-in simulator (24q)**. It embeds the
verified eight-protein graph and preprocessing provenance; no supporting uploads
or BLAST installation in Qollab are required. Settings: multi_job, COBYLA, p=1,
1,000 shots, seed 42, 30 training evaluations plus one fresh final per attempt,
alpha=100, increment=500, at most one retry. Save the complete output, including
REPORT_JSON and any failure, separately from the synthetic comparison results.
Review status, QAOA/exact affinity trees, unrooted RF, attempts and job counts.
RF=0 means agreement with exact graph-cut recursion, not evolutionary accuracy.
The first live protein QAOA run and local model-based NJ are verified.

Rebuild after regenerating the local protein graph with:
`python prepare_qollab.py --graph results/protein1822/graph.json --output qollab_smoke/protein1822_cobyla.py`
(run from the project root). The builder checks the adjacent manifest and
affinity matrix before embedding them. `qphylo_proteins.py` performs local
preprocessing; D19 records its conventions and verified BLAST execution.

**Current execution choice:** built-in local simulation has completed nine jobs
in one run. Use `execution_mode = "multi_job"` for local simulators. The observed
IonQ path allows one job per code run; use `execution_mode = "one_job"` there
and continue from saved coordinate-search checkpoints. The generic probe defaults
to `one_job`; the new protein COBYLA file explicitly uses `multi_job`.

Decisions and reasons are maintained in [DESIGN_DECISIONS.md](../DESIGN_DECISIONS.md).
The optimized cut, four-taxon tree and eight-taxon recursion passed on BasicSimulator.
The synthetic scripts remain available for repeat experiments. The
[eight-taxon guide](EIGHT_TAXA_GUIDE.md) gives pass criteria, then paired
coordinate/COBYLA scripts, saved-output folders and CSV/PNG/PDF collection.
The corrected five-seed paired comparison is complete; its CSV/PNG/PDF outputs
are in `../results/eight_taxa_comparison/` (D18).
The source uses the backend selected in Run Experiment. It does not instantiate
StatevectorSampler or silently replace a selected noise model with an ideal one.

The one-X and two-X jobs returned the expected measurements. Their post-run
`main` error remains unexplained. The Bell test completed without that error:
00=499, 01=3, 10=6, 11=492, 1,000 shots, 12.1 seconds, Qiskit 2.5.2,
IonQSimulatorBackend. Its exact simulator/noise selection was not recorded.

All prepared scripts use the tutorial's top-level execution style. Paste the
selected script into the **main Python file** of the existing Python/Qiskit
project. A source filename on disk does not imply a required `main()` function.
Do not paste two experiment scripts into the same main file at once.

## 1. Check ordering and rotations

Both bit-order checks have now completed cleanly on **IonQ Aria 1 (Noise
Model)**: qubit 0 gave 001=999 and 000=1; qubit 2 gave 100=1000, each with 1,000
shots. The rotation checks also completed cleanly on the same selection:
`rx_pi` gave 1=1000, `rz_phase` gave 1=1000, and `zz_phase` gave 01=1000.
Outputs are in `Next Check Results/`. The fixed cut probe has also passed below.

Paste `next_checks.py`. Fill `backend_label` with the exact simulator and noise
selection from Run Experiment. Change `test_name` for each run; each run submits
one job with 1,000 shots. Save each output HTML in this folder.

| `test_name` | Circuit purpose | Ideal measured string |
| --- | --- | --- |
| `bit_order_q0` | Flip qubit 0 among three qubits | `001` |
| `bit_order_q2` | Flip qubit 2 among three qubits | `100` |
| `rx_pi` | Bind an RX parameter to pi | `1` |
| `rz_phase` | H, bound RZ(pi), H converts phase into a measurable flip | `1` |
| `zz_phase` | Test the CX-RZ(pi)-CX cost interaction with interference | `01` |

These are expected ideal outcomes, not hard thresholds for a noisy simulator.
Record off-pattern counts. If a bit-order test favors the opposite string,
investigate that mapping before running cut scoring. If parameter binding or a
rotation fails, preserve the error and do not proceed to the cut probe.

Qiskit count strings display the highest bit index on the left. The adapter
reverses these strings so tuple position i means vertex/qubit i. The earlier
`00`/`11` Bell result could not distinguish a reversed mapping.
Reference: [Qiskit bit ordering](https://quantum.cloud.ibm.com/docs/en/guides/bit-ordering).

## 2. Run one four-qubit cut probe

The first fixed probe has completed cleanly. Independent rescoring confirms
1,000 shots, 196 feasible shots (19.6%), 51 optimal-cut shots (5.1%), selected
AB|CD, cut=4, gap=0 and Ncut=4/11. Backend metadata reports noise_model=aria-1.
The output is saved as `qaoa_cut_standaloneResults_...html` in this folder.
Angles were fixed, so this verifies execution/decoding/scoring rather than
optimization. The single optimized cut has also passed; section 3 covers the next
toy-tree check.

Paste **`qaoa_cut_standalone.py`**, not `cut_pilot.py`. The standalone file embeds
the core and adapter, so the first run does not depend on Qollab importing other
uploaded project files. Its configuration is near the bottom of the file.
Leave `mode = "fixed_probe"` and enter `backend_label`.

The toy graph has strong A-B and C-D affinities (9) and weak cross-pair
affinities (1). For target size c=2, the exact cut is **4**, attained by AB|CD.
The exact current-graph Ncut is **4/11**, about 0.363636. No biological input
is used in this check.

The probe runs p=1, gamma=0.37, beta=0.41, alpha=100, 1,000 shots. It submits
**one job**, reports feasible sample rate, the best observed feasible cut, its
gap from exact, and the backend/job information. It does not optimize angles
and does not reinforce the penalty. No feasible sample is a recorded probe
outcome, not permission to substitute an exact or repaired partition.

An independent mathematical calculation predicts about 19.0% feasible shots
and 5.0% optimal-cut shots for these particular angles. These are numerical
reference predictions, not Qiskit or IonQ measurements. This probe validates
execution and scoring; obtaining a good cut among many samples alone does not
demonstrate that QAOA improves on uniform random sampling.

Save the complete output, especially the final `REPORT_JSON` line. Compare the
partition using vertex-order tuples, not the printed measurement-string order.

## 3. Optimize one cut, then try the toy tree

The original multi-job optimized cut has now passed on the built-in simulator:
eight training jobs and one fresh final job, each 1,000 shots, all DONE on
BasicSimulator (Qiskit 2.5.2). Independent rescoring of every evaluation confirms
AB|CD, cut=4, gap=0, local Ncut=4/11; final feasible rate=70.2%, optimal-cut
shots=305/1,000. Mean penalized energy was 110.16 initially and 86.166 in the
fresh final sample. Alpha=100 was accepted without reinforcement. This is a
single ideal toy run, not a biological comparison or an advantage claim.

Output: `noamurphy-qaoa-phylogenetic-tree-reconst-output-completed-2026-10-10T19-32-54.722Z.html`.
Its manual backend_label still says IonQ Aria 1; job metadata establishes that
this run used BasicSimulator. Preserve that output and fix the label on new runs.

### Four-taxon tree passed; instructions for repeating it

Saved output: `qaoa_toyTestResults_noamurphy-qaoa-phylogenetic-tree-reconst-output-completed-2026-10-10T19-58-08.397Z.html`.
All 18 jobs DONE on BasicSimulator, 1,000 shots each. Both candidate-size attempts
accepted alpha=100 without retry; final feasible rates were c=1:47.1%, c=2:46.2%.
The current-graph Ncuts were 4/3 and 4/11, so AB|CD won. QAOA and exact trees
match. All measured counts and optimizer decisions were independently replayed.
This exercises the root selection; pair children terminate without further jobs.
For deeper recursion, run the prepared eight-taxon case before declaring multi-level
Qollab behavior verified. The first-draft workflow is D13; D14/D15 specify this
synthetic optimizer pilot.

Paste the updated **qaoa_cut_standalone.py** and set these values near the bottom:

```python
mode = "toy_tree"
execution_mode = "multi_job"
backend_label = "Built-in simulator (ideal)"
```

Select the built-in simulator in Run Experiment too; changing backend_label only
records a description and does not select a backend. Leave resume_json empty.
The same sampled solver runs all required recursive cuts in this one invocation.
Expect status `completed`, QAOA and exact toy topologies with the AB|CD split
(`(('A','B'),('C','D'));`), and traces showing the root c=1/c=2 candidates. Compare
splits rather than treating child order as biologically significant. Inspect
every attempt's feasibility/penalty action and preserve the full REPORT_JSON.
If the run fails, save its report rather than accepting the exact diagnostic tree
as the QAOA result. Total requested jobs remain bounded by the pilot budgets below.

To repeat the single cut, use `mode = "optimized_cut"` with the same multi_job
setting. No optional algorithm package or StatevectorSampler is required.

### IonQ continuation and shared pilot settings

For IonQ, choose `execution_mode = "one_job"`, enter the actual backend label
and select that backend in Run Experiment. The checkpoint path uses one new
measurement per code run. Both execution modes use the same pilot settings:

- p=1, 1,000 shots per evaluation, optimizer initialization seed 42.
- Seeded coordinate search, with at most **9 jobs per penalty attempt**, including
  a fresh final sample at the best parameters. These are separate code runs only
  in one_job mode; multi_job mode completes them in one invocation.
- Initial penalty alpha=100; increment 500 on a feasibility failure.
- At most one retry: **18 jobs maximum for one cut**, 18,000 requested shots total.
  On the restricted IonQ path that means at most 18 separate code runs.
- Final-sample acceptance requires the lowest observed penalized-energy sample
  to have exactly c ones. Finding some other feasible samples does not override
  this reinforcement trigger. Among feasible samples, choose minimum actual
  cut, breaking ties by the first selected-index tuple.
- No exact fallback, bit projection, or complement merging. Backend errors and
  timeouts stop execution rather than consuming penalty retries.

This is a deliberate pilot optimizer choice. The library panel now shows
optional Qiskit Algorithms 0.4.0 and Qiskit Optimization 0.7.0. The supplied QAOA
tutorial ran QAOA/COBYLA and SciPy imports successfully in Qollab. Our default
pilot still uses coordinate search; selectable COBYLA now calls SciPy directly
and is prepared in `eight_cobyla.py`. Neither optional wrapper package is
required by our adapter. Enabling packages does not remove the one-job limit
or provide automatic persistence between runs. COBYLA currently requires local
multi_job mode; the existing coordinate checkpoint path remains available.

FO used COBYLA; this pilot is not an optimizer-matched reproduction. Its QAOA
formulation, measured expectation objective, and recursive solver interface are
prepared. The eight-taxon comparison uses matched 30-evaluation training budgets
for both methods; FO used COBYLA(maxiter=300). Actual COBYLA execution in our
adapter has completed for five seeds. Settings and differences are recorded
in D14 and the eight-taxon guide. Reference:
[Qiskit Optimization COBYLA](https://qiskit-community.github.io/qiskit-optimization/stubs/qiskit_optimization.optimizers.COBYLA.html).

The seed controls optimizer initialization, not IonQ shot/noise randomness.
Gamma is searched in [0,2*pi] and beta in [0,pi]. Cost phases and objective are
divided by the largest absolute nonconstant Ising coefficient (minimum scale
1). This preserves the fixed-penalty minimizer, changes the angle units, and
must remain in the recorded configuration. The identity offset is omitted only
from circuit phases; it is retained in energy and expectation comparisons.

After each run, look for `Optimization status: awaiting_next_run` and save the
complete output. This status is a planned pause, not a failure. The
`RESUME_STATE_JSON` line contains measured counts and exact optimizer requests.
The next script reconstructs past decisions from those counts and submits only
the next missing measurement. It does not submit old jobs again.

From the project root, prepare the next script from the saved output:

```bash
python prepare_qollab.py --resume "qollab_smoke/path-to-saved-output.html"
```

This updates `qaoa_cut_standalone.py` with the checkpoint and its original
configuration. Paste the regenerated file into Qollab and run once. Repeat
until status is `completed` or `failed`. Codex can prepare the next file when
each saved output is added to the folder. Re-running the initial script without
its checkpoint starts over; use the regenerated script to continue.

The helper can also recover the first completed measurement from the old failed
multi-job output if its `REPORT_JSON` includes that measurement and the known
one-job-limit error. Save that output before starting again. Other backend
errors/timeouts are not resumed automatically; inspect the job first. Do not
change graph, shots, optimizer settings or the backend label while resuming.
Keep the same actual Run Experiment selection too.

For `mode = "toy_tree"`, the four-taxon root requires c=1 and c=2. Depending on
the returned split, a three-taxon child can require one additional cut. This
means at most **54 jobs / 54,000 shots** under these pilot limits: one invocation
locally in multi_job mode, or separate runs on IonQ in one_job mode. Queue and execution time are combined in elapsed time
because the supplied backend interface has not exposed separate timing fields.

## Local work completed and remaining verification

The core is unchanged. New reusable adapter: `../qphylo_qaoa.py`.

Run from the project root:

```bash
conda activate qphylo
python qphylo_core.py --self-test
python -m unittest discover -s tests -v
python prepare_qollab.py
```

The last command without `--resume` regenerates the initial standalone file
without a checkpoint. It combines the core, adapter, `qphylo_run_steps.py`,
`qphylo_benchmark.py`,
and `cut_pilot.py`; edit those source files rather than editing the generated
copy when keeping lasting changes. Small run settings may be changed directly
in Qollab, but record them with the output.

No new packages were installed for final-entry preparation. Earlier protein/NJ
work installed BLAST+2.17.0, R4.5.3, ape5.8.1 and phangorn2.12.1 in qphylo with
user authorization. Python3.11.17, JupyterLab4.6.4 and ipykernel7.4.0 remain.
Local Qiskit, NumPy, SciPy, Aer and IonQ provider are absent. The current suite
has 55 tests: 54 pass and one real-Qiskit Statevector test is skipped until
Qiskit is available. Gate-plan checks use a separate mathematical reference;
backend and COBYLA protocol tests use explicit test doubles, not cloud jobs.
The collector/figures have been checked on saved four/eight-taxon and protein
outputs; the NJ suite includes a real R regression when dependencies are present.

The Qiskit builder uses H, RX, RZ, CX and measurement. Parameter binding is
completed before backend submission. No IBM Runtime Sampler, D-Wave packages,
credentials, or provider instantiation is used. References:
[Qiskit circuits](https://quantum.cloud.ibm.com/docs/en/api/qiskit/2.3/qiskit.circuit.QuantumCircuit),
[IonQ provider gate serialization](https://qiskit-community.github.io/qiskit-ionq/_modules/qiskit_ionq/helpers.html).
These upstream references do not establish Qollab's installed provider behavior;
the prepared Qollab checks are still required.

## Biological experiment: completed pilot and remaining limits

The protein bridge, live QAOA run and model-NJ comparison are complete (D19,
D22, D23). Both exact and QAOA paths received the same affinity matrix. Explicit
symmetry averaging, zero diagonal, reused supplied alignment and a one-retry
limit are documented adaptations. They do not establish full original-method
fidelity. A controlled benchmark with verified generating trees and additional
protein seeds remains future work; it is deferred until submission is ready.

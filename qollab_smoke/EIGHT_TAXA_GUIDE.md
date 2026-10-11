# Eight-taxon recursion and optimizer pilot

The eight-taxon smoke test passed on Qollab: 72 DONE jobs, three split nodes and
RF=0. All five COBYLA seeds completed and recovered the exact graph-cut tree.
The user corrected the coordinate outputs and completed all five pairs (D18).
All ten trees match exact topology, with no failed runs. Current results are in
`results/eight_taxa_comparison/`. The recovery notes below document the earlier
seed issue and how to prepare another run; recovery is no longer pending.
Use the **built-in ideal simulator** throughout. No local quantum packages need
to be installed.

## Recovery for the saved seed mismatch

The four original seed-specific preparation files are no longer present.
For a new run, regenerate a separate script using `--seed` and `--output`, as
shown at the end of this guide, then paste the entire generated script.

Use the same built-in simulator. Check the printed initialization seed before
saving, and the stored `config.seed` afterward. The top-level `seed = ...`
in the run settings controls the experiment. Changing a filename or only the
embedded QAOAConfig class's default seed does not change those run settings.
The cause of the unchanged settings cannot be determined from output alone.

Save new outputs as `coordinate_seed43.html` through `coordinate_seed46.html`
in `eight_taxa_results/comparison/`, replacing the four mislabeled outputs.
Keep the valid coordinate seed 42 and all five valid COBYLA outputs.
Do not edit stored report seeds: the circuits were actually sampled with seed 42.

To inspect filenames against recorded settings without generating figures:

```bash
python collect_qaoa_results.py qollab_smoke/eight_taxa_results/comparison --inspect
```

After the four replacements, the folder contains the ten intended runs. Use
the original collection command:

```bash
python collect_qaoa_results.py qollab_smoke/eight_taxa_results/comparison --output results/eight_taxa_comparison --plot
```

The earlier `results/eight_taxa_seed42_preview/` was superseded and is no longer
present. Use the completed five-pair results in `results/eight_taxa_comparison/`.

## 1. Run the recursion smoke test first

Paste all of [eight_smoke.py](eight_smoke.py) into Qollab's main Python file.
Select **Built-in simulator (24q)** in Run Experiment. Settings are already:
hierarchical8, toy_tree, multi_job, coordinate search, p=1, 1,000 shots, seed 42,
nine jobs per penalty attempt and at most one reinforcement retry. Backend shot
sampling is seeded with `20261010 + seed` via the local backend's seed_simulator
option; initialization and sampling seeds are separately recorded.

The synthetic graph has affinity 9 within each pair AB/CD/EF/GH, 3 between
different pairs in the same four-taxon half, and 1 across ABCD|EFGH. It is a
hand-built graph, not a biological sequence dataset or a known evolutionary tree.

The exact reference is:

```text
((('A','B'),('C','D')),(('E','F'),('G','H')));
```

Successful validation means:

- Status `completed`, with three recursive split nodes, including child graphs.
- Unrooted RF distance to exact = 0 (child order need not match the printed string).
- Root split ABCD|EFGH: cut 16, local Ncut 8/19 (~0.421053).
- Child splits AB|CD and EF|GH: cut 12, local Ncut 0.8 each. These denominators
  use the four-taxon child graph, not the original eight-taxon degrees.
- All final answers accepted with their recorded size; inspect reinforcement
  outcomes, failures and all candidate cut gaps rather than just the topology.

An exact-topology run solves eight candidate cuts: root c=1..4, then c=1..2
in each child. With nine evaluations and no retries, expect 72 jobs/72,000
shots. A different split shape or retries changes the cost. The strict worst
case over any eight-taxon binary recursion is 270 jobs/270,000 shots. The script
prints this upper bound before submitting. Save complete output even if it fails.

Save the output under `eight_taxa_results/smoke/`, e.g. `coordinate_seed42.html`.
Keep these nine-evaluation smoke results separate from the comparison below.

## 2. Run a matched optimizer comparison after the smoke check

Use [eight_coordinate.py](eight_coordinate.py) and
[eight_cobyla.py](eight_cobyla.py). Both use the same fixture, backend, shots,
QAOA depth, penalties, angle scaling, initial angles and fresh-final rule.
Each has at most **30 training jobs + one fresh final job per penalty attempt**.
Both submit angles in gamma [0,2*pi], beta [0,pi]. COBYLA proposals are clipped
to that domain before circuit submission; both select their best measured
training angles for fresh final validation.

Run each optimizer once for each `seed` in **42, 43, 44, 45, 46**: ten runs total.
Change only the seed assignment in each script's settings near the bottom.
The simulator seed expression changes with it, and matches between optimizers
for each pair. This controls the local sampling RNG but does not make different
circuits produce identical samples. These are five initialization replicates on
one graph, not five biological datasets or proof of general performance.

COBYLA imports SciPy inside Qollab; the supplied QAOA tutorial established its
availability. It does not require the optional Qiskit Algorithms/Optimization
wrappers. Actual COBYLA execution in our adapter has now completed for five seeds;
local protocol checks use an explicit minimizer test double because SciPy is
not installed in the local project environment.

FO's inspected Qiskit path uses `COBYLA(maxiter=300)`. Its wrapper calls SciPy's
COBYLA; our direct SciPy call uses the same algorithm, with `rhobeg=1.0` and
`tol=1e-4`. The first study uses 30 training evaluations for **both** methods,
not 300 for COBYLA versus 30 for coordinate search. Bound clipping, explicit
initialization, scaled Hamiltonian and adaptive penalty are recorded choices,
so this is **an optimizer comparison inside our adapted NMcut**, not a complete
FO reproduction. SciPy versions matter (the implementation changed in 1.16);
the report records the version and minimizer status/message.
References: [SciPy COBYLA](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-cobyla.html),
[Qiskit COBYLA wrapper](https://qiskit-community.github.io/qiskit-optimization/stubs/qiskit_optimization.optimizers.COBYLA.html),
[inspected FO source](https://github.com/Nicofero/PhyloTree_ReConst/blob/main/qaoa_functions.py).

COBYLA budget exhaustion is an optimizer termination, not necessarily a failed
cut: feasibility is decided from the fresh final sample. Backend failures still
stop the run. COBYLA checkpointing is not implemented; use local multi_job mode.
The coordinate-search IonQ checkpoint path is retained.

The expected exact-topology/no-retry cost is up to 248 jobs/248,000 shots per
comparison run, 2,480 jobs/2.48 million shots across ten runs. Early optimizer
convergence can reduce it. Worst-case arbitrary split shape plus retries is
930 jobs per run. These are local simulator jobs, not requested paid hardware runs.

Save every output, including failures, in `eight_taxa_results/comparison/`:

```text
coordinate_seed42.html    cobyla_seed42.html
coordinate_seed43.html    cobyla_seed43.html
coordinate_seed44.html    cobyla_seed44.html
coordinate_seed45.html    cobyla_seed45.html
coordinate_seed46.html    cobyla_seed46.html
```

Don't overwrite earlier results. If repeating a seed, keep it elsewhere and
explicitly choose one output for the first paired analysis. Changing only
optimizer or seed is intended; changing budget/backend/graph requires a new study.

## 3. Collect data and create figures

From the project folder:

```bash
python collect_qaoa_results.py qollab_smoke/eight_taxa_results/comparison --output results/eight_taxa_comparison --plot
```

The helper reads HTML without executing it, checks stored scores against raw
counts, recomputes current-graph cut diagnostics, and rejects unmatched graph,
budget or backend settings and duplicate optimizer/seed replicates. It currently
requires the built-in BasicSimulator for this ideal study. Failed runs remain
in the summary with missing topology metrics; they aren't converted to successes.

Outputs:

- `runs.csv`: trees, unrooted RF to the exact reference, jobs, requested shots,
  recorded job elapsed time, retries and run failures.
- `evaluations.csv`: every sampled objective value, best training value so far,
  parameters, initialization, penalty, feasible rate, job ID and SciPy version.
- `final_cuts.csv`: each candidate's final feasible rate, optimal-shot rate,
  actual cut gap, local Ncut and analytic uniform-bitstring reference rates.
- `manifest.json`: original files, paired seeds, graph identity and interpretation.
- `optimizer_comparison.png` and `.pdf`: training curves for the shared root
  c=4 problem, paired final optimal-shot frequencies, whole-tree RF distances
  and whole-tree job costs. PNG is for slides; PDF is for sharp exported figures.

CSV export uses standard Python. Figures use already-installed base R here;
no matplotlib installation is required. Codex can run the collector after
outputs are saved. The pipeline has been checked on the actual four-taxon run
and the eight-taxon seed-42 pair:
`../results/toy4_validation/` contains a **single-method preview**, not an
eight-taxon/COBYLA result.

We plot the same root cut rather than averaging incompatible child-graph
energies. Child-graph diagnostics remain available in the CSV. Uniform rates
are analytic baselines, not separately measured random-solver runs. If both
optimizers recover the tree, frequency, convergence and resources can still
differ; if not, show the mismatches and failures. Five seeds support a small
descriptive pilot, not a significance or quantum-advantage claim.

## Rebuild the prepared files

```bash
python prepare_qollab.py --preset eight_smoke
python prepare_qollab.py --preset eight_coordinate
python prepare_qollab.py --preset eight_cobyla
```

Optional `--seed 43` embeds a different seed. Presets generate separate files
and do not overwrite qaoa_cut_standalone.py. Lasting implementation changes go
in the source modules, then rebuild. The default four-taxon bundle and IonQ
coordinate checkpoints remain available.

Use `--output` to retain a separate file for each seed, for example:

```bash
python prepare_qollab.py --preset eight_coordinate --seed 43 --output qollab_smoke/seed_runs/eight_coordinate_seed43.py
```

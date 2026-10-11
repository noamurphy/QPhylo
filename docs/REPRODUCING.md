# Reproduce the pilot

Run local commands from the project root. Saved results remain evidence;
reproductions below use new output folders. Qollab quantum runs are started
manually through Run Experiment, not by these local commands.

## Environment

Use the existing environment with `conda activate qphylo`. A new installation can
use `conda env create -f environment.yml`, then activate `qphylo`. The recipe
records the versions used on macOS/Apple Silicon; installation on a fresh
machine has not been verified. Dependencies may resolve differently elsewhere.

Local recorded versions: Python 3.11.17, BLAST+ 2.17.0, R 4.5.3, ape 5.8.1,
phangorn 2.12.1, JupyterLab 4.6.4, ipykernel 7.4.0. Qollab's verified protein
run recorded Qiskit 2.5.2 and SciPy 1.18.0. Browser libraries are managed by
Qollab; local Qiskit/NumPy/SciPy installation is not needed for these audit steps.
The broad upstream requirements file is not our installation recipe.

JupyterLab is installed in qphylo. Launch `jupyter lab` from the project directory
to browse local files using that environment. For another Jupyter installation,
register the kernel with
`python -m ipykernel install --user --name qphylo --display-name 'Python (qphylo)'`.

## Core and execution checks

```bash
python qphylo_core.py --self-test
python -m unittest discover -s tests
python prepare_submission.py
```

The suite currently has 55 tests: 54 pass and one real-Qiskit statevector test
skips in the local environment. Mocked backend/minimizer checks verify control
flow and limits; they are not measured quantum experiments. With R dependencies
present, the NJ suite includes a real local regression. The builder verifies
training counts, provenance and comparison hashes, then embeds audited data and
frozen cut angles. The command also prepares `qollab_upload/main.py`, `LICENSE`
and `CREDITS.md`, preserving the maintained project `README.md` and `assets/`
for uploading together at the Qollab project root. Run
`Rscript prepare_project_header.R` separately to regenerate the header from saved
topologies and copy the unchanged results plots. This does not run new experiments. See
[final Qollab check](QOLLAB_FINAL_CHECK.md) for live runs.

## Protein preprocessing

The selected alignment and supplied reference are copied from the credited
upstream repository; their relative paths are retained in the curated release.

```bash
python qphylo_proteins.py --validate-only --output results/reproduce_input_check
python qphylo_proteins.py --output results/reproduce_protein1822
```

The second command runs BLASTP locally with the documented defaults. Inspect
the new manifest, graph, affinity matrix and connectivity before using it.
An exact small-graph reference is exported for diagnosis. Do not infer
evolutionary truth from that graph-reference topology.

To create a new paste-ready protein runner:

```bash
python prepare_qollab.py --graph results/reproduce_protein1822/graph.json --output qollab_smoke/reproduce_protein1822.py
```

Paste that script into Qollab `main.py`, select the built-in ideal simulator and
save its complete output separately. This experiment-only runner optimizes a
full tree; the final `submission_main.py` is the entry with backend routing.

## Audit and visualize saved QAOA results

```bash
python collect_qaoa_results.py qollab_smoke/protein1822_cobyla_run1.html --output results/reproduce_protein_audit
python collect_qaoa_results.py qollab_smoke/eight_taxa_results/comparison --output results/reproduce_eight_taxa --plot
```

The collector independently rescores all counts and validates optimizer decisions,
feasibility, traces, seeds and backend metadata. Only one report per optimizer/seed
is allowed; choose repeats explicitly. PNG/PDF plotting calls Rscript. Protein
topology plotting is produced by the separate NJ comparison below.

## Model-based NJ and comparison

```bash
python qphylo_nj.py --output results/reproduce_nj
```

Defaults use the original supplied alignment, verified protein QAOA HTML,
reference file, 100 site bootstraps and seed42. The alignment's taxa/order and
SHA-256 must match the QAOA provenance. Distances use WAG, five Gamma categories,
shape0.5, model frequencies, pairwise gap exclusion and no separate invariant
component. The existing alignment is reused rather than regenerated with MUSCLE.

Outputs include distances, all bootstrap trees, supports, four comparison
topologies, RF tables, PNG, versions and hashes. Negative NJ branches are retained;
topology plots and unrooted RF ignore branch lengths. Bootstrap support is
resampling stability, not a biological ground-truth label.

## Extend without changing the question accidentally

Use the same alignment and graph for exact/QAOA comparison. Keep model-NJ based
on the same taxa and document preprocessing differences. Record new choices in
`DESIGN_DECISIONS.md`. Increase seeds, controlled divergence or circuit depth
only after the small-input audit remains valid. Scaling beyond these eight taxa
has not been established, and exact enumeration is exponential.

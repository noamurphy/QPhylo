# First protein-derived QAOA reconstruction

The saved Qollab run completed successfully on the built-in ideal simulator.
QAOA reconstructed the same eight-taxon topology as exact fixed-size-cut
recursion on the BLASTP-derived affinity graph. Unrooted RF distance is 0.

| Check | Result |
| --- | --- |
| Circuit jobs | 403, all DONE |
| Measurements | 403,000; 1,000 per job |
| Recursive split nodes | 4 |
| Fixed-size cut problems | 9, all accepted |
| Penalty reinforcement | 5 problems retried from alpha 100 to 600 |
| Accepted candidate cut gaps | Zero within floating-point roundoff, for all 9 problems |
| Backend | BasicSimulator, no noise model; Qiskit 2.5.2 |
| Optimizer | SciPy 1.18.0 COBYLA, p=1, initialization seed 42 |

The root divides taxa 1/6/7/8 from 2/3/4/5. Recursion isolates taxon8,
then taxon1, retaining the pairs 6/7, 2/3 and 4/5. This describes the
algorithm's recursion, not an inferred evolutionary root or branch lengths.

## What reinforcement accomplished

The three root candidates selecting 1, 2 or 3 taxa, the three-taxon child
problem, and the 2/3/4/5 child's one-taxon candidate needed reinforcement.
At alpha100 their lowest observed penalized-energy answers had the wrong
number of selected taxa. Each was reoptimized at alpha600 and accepted on
its only allowed retry. Exhaustive checks also show these five alpha100
objectives have an infeasible global minimum. This is substantive live evidence
that the reinforcement path operates as designed on the larger protein weights.

Reinforcement does not guarantee every measurement has the correct size or
improve every sampling percentage. Acceptance checks the lowest observed
penalized-energy sample, then chooses the minimum actual cut among feasible
samples. The root two-taxon candidate had only 3 optimal measurements in 1,000
on its accepted attempt, so repeatability still needs more seeds.

## Root candidates

The algorithm compares the accepted minimum-cut candidates using normalized
cut on the current graph. Smaller Ncut wins.

| Selected taxa count | Accepted alpha | Ncut | Optimal shots / 1,000 |
| --- | --- | --- | --- |
| 1 | 600 | 1.078915 | 63 |
| 2 | 600 | 1.107349 | 3 |
| 3 | 600 | 1.047917 | 4 |
| 4 | 100 | **1.044440** | 31 |

The four-taxon candidate wins by about 0.003477 over the three-taxon candidate.
Its optimal-shot frequency is 3.1%, versus the analytic uniform-bitstring
expectation of 0.78125%; that comparison alone is not evidence of optimizer
superiority. Uniform sampling has a high chance of hitting an optimum somewhere
among 1,000 shots, and this is one seed on one graph.

## Audit and limits

The local audit independently rescored all 403 count dictionaries, validated
the recursive partitions and scores, and matched the graph and complete
preprocessing provenance to the local exports. Each penalty attempt uses one
final job at the best measured training angles; all five retries follow the
recorded increment. Exact calculations are diagnostics, not substituted cuts.
The simulator seed is reused per job: a fresh final job is a separate execution,
not statistically independent sampling when angles and seed are identical.

Nine COBYLA attempts reached their 30-training-evaluation limit; five stopped
when the trust-region radius reached its lower bound. A budget-limit optimizer
status is not a failed circuit job or an infeasible accepted split.

Tiny local/browser score differences (maximum 2.73e-12) required the collector
to allow floating-point roundoff: relative tolerance 1e-12 and absolute tolerance
1e-12 for calculated floating-point fields. Discrete sample fields remain exact;
changed shot totals and meaningful score changes are still rejected. The original
HTML is unchanged. Forty-two local tests pass; the real-Qiskit test remains skipped.

The recorded job elapsed times sum to about 5.39 seconds; this excludes package
loading and other program work. No biological accuracy, noise robustness or
quantum advantage is established by matching the exact affinity tree.

## Next milestone

Prepare model-based neighbor joining from the same eight proteins, then compare
QAOA, exact affinity NMcut and NJ alongside the supplied reference tree, with
reference provenance stated. This completes the missing biological baseline for
the first draft. Additional QAOA initialization/simulator seeds can then assess
repeatability without delaying that integration step.

Source: [original Qollab output](../../qollab_smoke/protein1822_cobyla_run1.html).
Validated exports: `runs.csv`, `evaluations.csv`, `final_cuts.csv`, `manifest.json`.

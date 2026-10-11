# Final Qollab entry: both execution paths verified

The user ran `main.py` with `RUN_MODE='auto'` on both backends. Independent audits
of the saved outputs passed. No post-run traceback or execution error was found.
Original HTML outputs are preserved unchanged.

| Actual backend | Auto selected | Completed jobs / shots | Verified outcome |
| --- | --- | --- | --- |
| BasicSimulator, ideal | Full tree | 403 / 403,000 | Four recursive split nodes; nine accepted cuts; five reinforcements; RF0 to exact NMcut |
| IonQSimulatorBackend, Aria-1 noise model | Single cut | 1 / 1,000 | 218 feasible shots (21.8%); 27 optimal-cut shots (2.7%); best observed cut gap 0 |

The IonQ partition was **{taxon1,taxon6,taxon7,taxon8} | {taxon2,taxon3,taxon4,taxon5}**,
cut 795.4323749543944, normalized cut 1.04444015003035. Angles were previously
trained on the ideal backend. This run sampled one eight-qubit circuit and
performed no new angle optimization or full-tree inference.

Every built-in count dictionary was rescored. The audit checked graph/config/
provenance, actual job metadata, angle bounds, initialization seeds, evaluation
budgets, Hamiltonian scales, final best-training-point selection and penalty
actions. Accepted measured cuts were replayed through the recursion, reproducing
the reported tree without using an exact fallback. RF against exact NMcut is 0;
RF against the saved NJ and provided reference is 2, as in the original pilot.
IonQ's raw counts, shot total, feasibility, expectation, optimum rate, partition,
normalized cut and gate counts were independently checked.

The ideal rerun has fresh job IDs and the same parameters/counts as the original
seeded run. This verifies repeat execution, not an independent reliability
replicate. The single noisy cut verifies compatibility, not hardware performance,
biological accuracy, a controlled noise study or quantum advantage. Historical
tree/NJ information remains separate from the current measurements in both reports.

Evidence:

- [Built-in original output](../../qollab_smoke/submission_builtin_run1.html)
- [IonQ original output](../../qollab_smoke/submission_ionq_run1.html)
- [Audit, source hashes and recorded versions](audit.json)
- [Rescored IonQ counts](ionq_counts.csv)

Recorded versions: Qiskit 2.5.2 on both paths; SciPy 1.18.0 for built-in COBYLA.

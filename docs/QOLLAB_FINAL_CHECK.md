# Final Qollab check: same file, two backends

**Completed and independently verified (D28):** the saved built-in run completed
403 jobs and matched the exact NMcut tree. The IonQ Aria-1 noise run completed
one job with 218 feasible and 27 optimal-cut shots/1,000, best cut gap 0. No post-run
errors were found. See [the audit summary](../results/submission_validation/summary.md).
The instructions below remain available for repeating the checks.

Upload [main.py](../qollab_upload/main.py), [LICENSE](../qollab_upload/LICENSE),
[CREDITS.md](../qollab_upload/CREDITS.md), [README.md](../qollab_upload/README.md),
and the `assets/` folder from `qollab_upload/` to the
Python/Qiskit project root, replacing its existing `main.py`. The code is identical
to [submission_main.py](../qollab_smoke/submission_main.py). Leave `RUN_MODE = 'auto'`.
Choose `assets/project_banner_v2.png` for the separate project header. Check its
crop and ensure that
both inline README figures render after upload. README/images and license/credits
are documentation; the program does not read them at runtime.
A function named `main` and supporting Python-file imports are unnecessary.

## 1. Built-in simulator (24q)

Select it in **Run Experiment**, then run once. Check the **CURRENT RUN** section:

- Backend `BasicSimulator`, mode `full_tree`, status `completed`.
- A fresh recursive QAOA grouping and Newick tree are printed.
- RF to exact NMcut is reported; the historical verified run had RF=0. Preserve
  any new discrepancy rather than overwriting it with the expected result.
- Current jobs are DONE; each circuit has 1,000 shots. The prior run used 403
  jobs; different package behavior may change optimizer termination and counts.
- `REPORT_JSON` contains current evaluations, attempts, trace and actual metadata.

Save the complete output as `qollab_smoke/submission_builtin_run1.html`.

## 2. IonQ simulator

Keep the same code and choose an IonQ free simulator in Run Experiment. Record
the exact ideal/noise selection with the saved file. Check:

- Mode `single_cut`, status `probe_completed`.
- **Exactly one** current circuit job, DONE, with 1,000 measured shots.
- Angles are explicitly identified as previously trained on BasicSimulator.
- Feasible shots, optimal-cut shots and best observed feasible partition are
  reported from this measurement. There is no fresh tree or live optimizer.
- No feasible sample is a valid recorded sampling outcome, not a reason to
  substitute an exact cut. Backend errors indicate a failed execution.

Save as `qollab_smoke/submission_ionq_run1.html`. A noisy result need not match
ideal counts or achieve the exact cut. Earlier gate/four-qubit tests and the
final eight-qubit Aria-1 noise-model execution have passed.

The **SAVED COMPARISON** table in both modes is historical context. Its 403 jobs,
trees and NJ support must not be counted as this run's work.

## If a run fails

Keep the entire output and backend selection. Do not switch to a sampler that
ignores the chosen backend, or accept the saved/exact tree as a new QAOA result.
Local source changes belong in `submission_entry.py`/the reusable modules; run
`python prepare_submission.py` to regenerate the upload set afterward.

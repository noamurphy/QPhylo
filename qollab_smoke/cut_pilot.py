"""Source for the standalone Qollab bundle; imports work locally with Qiskit.

The default is ONE fixed-parameter QAOA probe, not an optimization run.
Switch mode only after the probe and gate/ordering checks pass on Qollab.
"""

from dataclasses import asdict
import json

from qphylo_core import cut_weight, exact_fixed_size_cut, normalized_cut, reconstruct
from qphylo_qaoa import (
    QAOAConfig, QAOACutSolver, QiskitBackendSampler,
    summarize_counts,
)
from qphylo_run_steps import run_one_job_step
from qphylo_benchmark import benchmark_metrics, job_upper_bound, toy_fixture

mode = "fixed_probe"  # fixed_probe, optimized_cut, toy_tree, or protein_tree.
execution_mode = "one_job"  # one_job for IonQ; multi_job for local simulators.
resume_json = ""  # prepare_qollab.py --resume embeds the saved checkpoint here.
backend_label = "UNRECORDED: copy exact simulator/noise selection here"
fixture = "toy4"  # hierarchical8 for the eight-taxon test.
optimizer = "seeded_coordinate_search"  # or COBYLA (local multi_job only).
seed = 42  # Use paired seeds 42,43,44,45,46 for the comparison.
evaluation_budget = 9  # Comparison: 31 total = 30 training + fresh final.
simulator_seed = None  # Comparison: integer seed for supported local backends.
dataset_provenance = None  # --graph embeds the audited preprocessing manifest.
weights, labels = toy_fixture(fixture)
target_size = 2
parameters = (0.37, 0.41)  # One gamma and one beta, in scaled-Hamiltonian units.
config = QAOAConfig(p=1, shots=1000, seed=seed, max_evaluations=evaluation_budget,
                    initial_penalty=100, penalty_step=500, max_retries=1,
                    optimizer=optimizer)

sampler = QiskitBackendSampler(backend, job_limit=1 if execution_mode == "one_job" else None)
solver = None
report = {"mode": mode, "execution_mode": execution_mode,
          "backend_label": backend_label, "weights": weights,
          "labels": labels, "config": asdict(config)}

try:
    if execution_mode not in ("one_job", "multi_job"):
        raise ValueError("execution_mode must be one_job or multi_job.")
    if resume_json and execution_mode != "one_job":
        raise ValueError("A checkpoint must continue in one_job mode.")
    if mode == "protein_tree" and execution_mode != "multi_job":
        raise ValueError("This protein pilot requires multi_job on the built-in simulator.")
    if simulator_seed is not None:
        if execution_mode != "multi_job":
            raise ValueError("simulator_seed is for local multi_job runs only.")
        backend.set_options(seed_simulator=simulator_seed)
    print("Fixture:", fixture, "optimizer:", config.optimizer,
          "taxa:", len(labels), "initialization seed:", config.seed)
    print("Maximum jobs for any tree shape:", job_upper_bound(len(labels), config))
    if mode == "fixed_probe":
        batch = sampler(weights, target_size, config.initial_penalty, parameters,
                        {"shots": config.shots, "phase": "fixed_probe"})
        stats = summarize_counts(weights, target_size, config.initial_penalty,
                                 batch.counts, config.shots)
        exact = exact_fixed_size_cut(weights, target_size)
        report.update({"target_size": target_size, "parameters": parameters,
                       "counts": dict(batch.counts), "summary": stats,
                       "exact_bits": exact, "exact_cut": cut_weight(weights, exact)})
        print("Exact fixed-size minimum cut:", report["exact_cut"])
        print(f"Feasible measurement rate: {stats['feasible_rate']:.1%}")
        sample = stats["best_feasible_sample"]
        if sample:
            report["observed_cut_gap"] = sample["cut"] - report["exact_cut"]
            report["observed_ncut"] = normalized_cut(weights, sample["bits"])
            print("Best observed feasible partition:", sample["bits"])
            print("Its cut:", sample["cut"], "gap from exact:", report["observed_cut_gap"])
        else:
            print("No feasible sample in this probe. This is a recorded result, not a fallback.")
        report["status"] = "probe_completed"
    elif mode in ("optimized_cut", "toy_tree", "protein_tree"):
        if execution_mode == "one_job":
            state = json.loads(resume_json) if resume_json else None
            report = run_one_job_step(weights, labels, target_size, mode, config,
                                      sampler, backend_label, state)
            report["execution_mode"] = execution_mode
        else:
            solver = QAOACutSolver(sampler, config)
            if mode == "optimized_cut":
                bits = solver(weights, target_size)
                exact = exact_fixed_size_cut(weights, target_size)
                report.update({"target_size": target_size, "bits": bits,
                    "cut": cut_weight(weights, bits),
                    "ncut": normalized_cut(weights, bits),
                    "exact_cut": cut_weight(weights, exact),
                    "cut_gap": cut_weight(weights, bits) - cut_weight(weights, exact)})
            else:
                tree, trace = reconstruct(weights, labels, solver)
                exact_tree, exact_trace = reconstruct(weights, labels)
                report.update({"tree": tree, "trace": trace,
                    "exact_tree": exact_tree, "exact_trace": exact_trace})
            report.update({"status": "completed", "evaluations": solver.evaluations,
                "attempts": solver.attempts,
                "completed_measurements": len(solver.evaluations)})
        print("Optimization status:", report["status"])
        print("Completed measurements:", report["completed_measurements"])
        if report["evaluations"]:
            latest = report["evaluations"][-1]
            print("Last measurement phase:", latest["phase"])
            print("Current penalty:", latest["penalty"])
            print("Measured mean penalized energy:", latest["summary"]["expectation"])
            print(f"Feasible measurement rate: {latest['summary']['feasible_rate']:.1%}")
        if report["status"] == "awaiting_next_run":
            print("Checkpoint saved below. Save this output; prepare the next run from it.")
            print("Next measurement phase:", report["next_request"]["context"]["phase"])
        elif report["status"] == "completed":
            if mode == "optimized_cut":
                print("Selected partition:", report["bits"], "cut gap:", report["cut_gap"])
            else:
                print("QAOA topology:", report["tree"])
                print("Exact affinity topology:", report["exact_tree"])
        else:
            print("Optimization failed:", report["error"])
        if "resume_state" in report:
            print("RESUME_STATE_JSON:", json.dumps(report["resume_state"], sort_keys=True))
    else:
        raise ValueError(f"Unknown mode: {mode}")
except Exception as exc:
    report["status"] = "failed"
    report["error"] = f"{type(exc).__name__}: {exc}"
    raise
finally:
    if solver is not None:
        report["evaluations"] = solver.evaluations
        report["attempts"] = solver.attempts
    report["jobs"] = sampler.jobs
    report["fixture"] = fixture
    report["simulator_seed"] = simulator_seed
    report["dataset_provenance"] = dataset_provenance
    # Diagnostics are computed after solving; they never pick or repair bits.
    benchmark_metrics(report)
    if report.get("topology_metrics"):
        print("Unrooted RF distance to exact:", report["topology_metrics"]["unrooted_rf_to_exact"])
        print("Recursive split nodes:", report["topology_metrics"]["split_nodes"])
    print("REPORT_JSON:", json.dumps(report, sort_keys=True))

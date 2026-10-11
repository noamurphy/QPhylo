"""Final Qollab runner, bundled by prepare_submission.py.

Qollab supplies backend. The generated file supplies RUN_MODE and SUBMISSION_DATA.
Auto runs a fresh tree on BasicSimulator; other backends get one frozen-angle cut.
"""

from dataclasses import asdict
from itertools import combinations
import json
from math import isclose

from qphylo_core import cut_weight, exact_fixed_size_cut, normalized_cut, reconstruct
from qphylo_qaoa import (QAOAConfig, QAOACutSolver, QiskitBackendSampler,
                         summarize_counts, decode_bitstring, hamiltonian_scale)
from qphylo_benchmark import benchmark_metrics, graph_signature, job_upper_bound, topology_splits


def submission_mode(backend, requested):
    if requested not in ("auto", "full_tree", "single_cut"):
        raise ValueError("RUN_MODE must be auto, full_tree, or single_cut.")
    ideal_basic = (type(backend).__name__ == "BasicSimulator" and
                   getattr(getattr(backend, "options", None), "noise_model", None) is None)
    if requested == "full_tree" and not ideal_basic:
        raise ValueError("Full-tree optimization requires Built-in simulator (BasicSimulator). "
                         "Use auto or single_cut for IonQ/other backends.")
    if requested == "single_cut":
        return "single_cut", "Explicit single-circuit mode."
    if ideal_basic:
        return "full_tree", "Built-in ideal simulator supports multiple jobs."
    return "single_cut", "Conservative one-job mode for the selected backend."


def display_trace_tree(trace, labels):
    """Render the measured recursive grouping as readable plain text."""
    rows = iter(trace)
    def group(taxa):
        if len(taxa) == 1:
            return taxa[0]
        if len(taxa) == 2:
            return tuple(taxa)
        row = next(rows)
        if list(row["taxa"]) != list(taxa):
            raise ValueError("Display trace taxa mismatch.")
        return (group(row["left"]), group(row["right"]))
    tree = group(labels)
    def show(node, indent="", last=True):
        if isinstance(node, str):
            print(indent + ("`-- " if last else "|-- ") + node)
        else:
            print(indent + ("`-- " if last else "|-- ") + "group")
            child_indent = indent + ("    " if last else "|   ")
            for i, child in enumerate(node):
                show(child, child_indent, i == len(node)-1)
    show(tree)


def display_saved_comparison(saved):
    print("\nSAVED COMPARISON - previous runs, not new measurements in this execution")
    print(f"{saved['n_taxa']} proteins; BLASTP affinities; model NJ from the same {saved['columns']}-column alignment.")
    print("Comparison                              Unrooted RF")
    for row in saved["pairwise_comparisons"]:
        print(f"{row['method1'] + ' vs ' + row['method2']:<40}{float(row['rf']):g}")
    print(f"Saved protein QAOA: {saved['submitted_jobs']} jobs, {saved['reinforcements']} penalty reinforcements.")
    model = saved["model"]
    print(f"NJ used {model['name']} + Gamma{model['gamma_categories']} (shape {model['gamma_shape']}), "
          f"{saved['bootstrap']['replicates']} site bootstraps.")
    for row in saved["nj_support"]:
        print(f"Saved NJ split {{{', '.join(row['side'])}}}: {row['percent']:g}% bootstrap support.")
    print("Reference generating-tree provenance is unverified; no accuracy/quantum-advantage claim.")


def run_submission(backend, requested, data):
    weights, labels = data["weights"], data["labels"]
    config = QAOAConfig(**data["config"])
    report = {"entry_point": "qphylo-submission-v1", "requested_run_mode": requested,
              "weights": weights, "labels": labels, "config": asdict(config),
              "fixture": "protein1822", "dataset_provenance": data["dataset_provenance"],
              "backend_label": type(backend).__name__ + " (actual class; see job metadata)",
              "saved_comparison": data["saved_comparison"], "status": "preparing",
              "graph_signature": data["graph_signature"],
              "simulator_seed": None, "evaluations": [], "attempts": [], "jobs": []}
    sampler = solver = None
    try:
        selected, reason = submission_mode(backend, requested)
        report.update(selected_run_mode=selected, mode_reason=reason,
                      mode="protein_tree" if selected == "full_tree" else "fixed_probe",
                      execution_mode="multi_job" if selected == "full_tree" else "one_job")
        if graph_signature(weights, labels) != data["graph_signature"]:
            raise ValueError("Embedded graph/provenance disagree.")
        print("QPhylo: protein phylogeny with QAOA normalized minimum cut")
        print("This effort is supported by Qollab & IonQ")
        print("Backend:", type(backend).__name__, "| Mode:", selected)
        print(reason)
        print("Input:", len(labels), "protein taxa | QAOA depth:", config.p, "| Shots/circuit:", config.shots)
        display_saved_comparison(data["saved_comparison"])
        print("\nCURRENT RUN - measurements below use the backend selected in Run Experiment")
        job_limit = job_upper_bound(len(labels), config) if selected == "full_tree" else 1
        report["current_job_limit"] = job_limit
        sampler = QiskitBackendSampler(backend, job_limit=job_limit)
        if selected == "full_tree":
            print(f"Fresh {config.optimizer} tree reconstruction: up to {config.max_evaluations-1} training "
                  "+ 1 final circuit per attempt.")
            print(f"Penalty starts at {config.initial_penalty:g}; step {config.penalty_step:g}, "
                  f"at most {config.max_retries} reinforcement per cut.")
            print("No measured failure is replaced by the exact solver. Maximum jobs:", job_limit)
            simulator_seed = data["training_simulator_seed"]
            backend.set_options(seed_simulator=simulator_seed)
            report["simulator_seed"] = simulator_seed
            solver = QAOACutSolver(sampler, config)
            tree, trace = reconstruct(weights, labels, solver)
            report.update(tree=tree, trace=trace, status="completed",
                          completed_measurements=len(solver.evaluations))
            report["evaluations"], report["attempts"] = solver.evaluations, solver.attempts
            report["jobs"] = sampler.jobs
            benchmark_metrics(report)
            actual = topology_splits(trace, labels)
            nj = {tuple(side) for side in data["saved_comparison"]["nj_splits"]}
            reference = {tuple(side) for side in data["saved_comparison"]["reference_splits"]}
            report["baseline_metrics"] = {
                "unrooted_rf_to_saved_model_nj": len(actual ^ nj),
                "unrooted_rf_to_provided_reference": len(actual ^ reference),
                "baseline_origin": "Previously computed model NJ and supplied reference for this same alignment"}
            print("\nFresh QAOA grouping (branch lengths and evolutionary root are not inferred):")
            display_trace_tree(trace, labels)
            print("Newick:", tree)
            print("RF to exact NMcut:", report["topology_metrics"]["unrooted_rf_to_exact"])
            print("RF to previously computed model NJ:", report["baseline_metrics"]["unrooted_rf_to_saved_model_nj"])
            print("RF to supplied reference:", report["baseline_metrics"]["unrooted_rf_to_provided_reference"])
            print("Penalty reinforcements:", sum(a["action"] == "reinforce" for a in solver.attempts))
            print("Identical simulator seeds make identical-angle repeats deterministic.")
        else:
            frozen = data["frozen_cut"]
            c, alpha, parameters = frozen["target_size"], frozen["penalty"], frozen["parameters"]
            if frozen["graph_signature"] != data["graph_signature"]:
                raise ValueError("Frozen-angle graph mismatch.")
            if not isclose(hamiltonian_scale(weights, c, alpha), frozen["scale"], rel_tol=1e-12):
                raise ValueError("Frozen-angle Hamiltonian scale mismatch.")
            print(f"Single circuit: root cut with exactly {c} of {len(labels)} taxa, using previously trained angles.")
            print("This samples one cut; it does not optimize angles or reconstruct a new tree here.")
            print("Angle source:", frozen["source_file"], "| seed", frozen["training_config"]["seed"],
                  "COBYLA, built-in ideal simulator.")
            batch = sampler(weights, c, alpha, parameters,
                            {"shots": config.shots, "phase": "frozen_protein_cut"})
            stats = summarize_counts(weights, c, alpha, batch.counts, config.shots)
            exact = exact_fixed_size_cut(weights, c)
            optimum = cut_weight(weights, exact)
            optimal_shots = sum(int(count) for key, count in batch.counts.items()
                if sum(decode_bitstring(key, len(labels))) == c and
                isclose(cut_weight(weights, decode_bitstring(key, len(labels))), optimum, abs_tol=1e-9))
            optimum_states = sum(isclose(cut_weight(weights, tuple(int(i in side) for i in range(len(labels)))),
                                               optimum, abs_tol=1e-9)
                                 for side in combinations(range(len(labels)), c))
            report.update(target_size=c, penalty=alpha, parameters=parameters, scale=frozen["scale"],
                          angle_source=frozen, counts=dict(batch.counts), summary=stats,
                          exact_bits=exact, exact_cut=optimum, optimal_shots=optimal_shots,
                          optimal_shot_rate=optimal_shots/config.shots,
                          uniform_optimal_rate=optimum_states/(2**len(labels)), status="probe_completed")
            print(f"Job completed | Feasible shots: {stats['feasible_shots']}/{config.shots} ({stats['feasible_rate']:.1%})")
            print(f"Exact minimum-cut shots: {optimal_shots}/{config.shots} ({optimal_shots/config.shots:.1%})")
            best = stats["best_feasible_sample"]
            if best:
                left = [label for label, bit in zip(labels, best["bits"]) if bit]
                right = [label for label, bit in zip(labels, best["bits"]) if not bit]
                report.update(observed_cut_gap=best["cut"]-optimum,
                              observed_ncut=normalized_cut(weights, best["bits"]))
                print("Best observed feasible cut:", ", ".join(left), "|", ", ".join(right))
                print(f"Cut: {best['cut']:.6f} | Gap to exact: {report['observed_cut_gap']:.6f}")
            else:
                print("No feasible sample observed. Recorded as measured; no retry or fallback.")
            print("Noisy/stochastic cut results need not match the saved ideal run.")
        return report
    except Exception as exc:
        report["status"] = "failed"
        report["error"] = f"{type(exc).__name__}: {exc}"
        print("Run failed:", report["error"])
        raise
    finally:
        if solver is not None:
            report["evaluations"], report["attempts"] = solver.evaluations, solver.attempts
        if sampler is not None:
            report["jobs"] = sampler.jobs
        report["graph_signature"] = graph_signature(weights, labels)
        if report.get("mode") == "protein_tree":
            benchmark_metrics(report)
        else:
            report["resource_metrics"] = {"submitted_jobs": len(report["jobs"]),
                "requested_shots": sum(j.get("shots", 0) for j in report["jobs"]),
                "job_elapsed_seconds": sum(j.get("elapsed_seconds", 0) for j in report["jobs"])}
        print("Current-run circuit jobs:", report["resource_metrics"]["submitted_jobs"])
        print("Current-run status:", report["status"])
        print("\nCredits: Onodera et al. (NMcut/Digital Annealer method); Fernandez-Otero et al.")
        print("(QAOA phylogeny work and supplied dataset); Qiskit, BLAST+, ape and phangorn.")
        print("Hackathon additions: BLAST affinity bridge, local recursive normalization,")
        print("bounded penalty reinforcement, explicit sampled-cut validation and NJ comparison.")
        print("Audit report follows; saved_comparison is historical, jobs/counts are this run.")
        print("REPORT_JSON:", json.dumps(report, sort_keys=True))


report = run_submission(backend, RUN_MODE, SUBMISSION_DATA)

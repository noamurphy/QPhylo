"""Resume the same sampled optimizer across one-job Qollab code runs.

Replay saved counts to reconstruct deterministic optimizer decisions; only the
next missing measurement is submitted. Checkpoints are explicit JSON, not an
assumption that Qollab preserves local files or variables between runs.
"""

from dataclasses import asdict
import json

from qphylo_core import cut_weight, exact_fixed_size_cut, normalized_cut, reconstruct
from qphylo_qaoa import QAOACutSolver, SampleBatch, summarize_counts


class NextRunRequired(Exception):
    def __init__(self, request):
        self.request = request
        super().__init__("The next measurement requires a new Qollab code run.")


def canonical(value):
    return json.dumps(value, sort_keys=True, allow_nan=False)


def sample_request(weights, c, penalty, parameters, context):
    return {"weights": weights, "target_size": c, "penalty": penalty,
            "parameters": parameters, "context": dict(context)}


class OneJobReplaySampler:
    """Read old measurements, submit at most one new job, then pause."""
    def __init__(self, sample, history):
        self.sample = sample
        self.history = json.loads(canonical(history))
        self.position = 0
        self.submitted = False

    def __call__(self, weights, c, penalty, parameters, context):
        request = sample_request(weights, c, penalty, parameters, context)
        if self.position < len(self.history):
            row = self.history[self.position]
            if canonical(row["request"]) != canonical(request):
                raise ValueError("Checkpoint does not match the next optimizer request.")
            summarize_counts(weights, c, penalty, row["counts"], context["shots"])
            self.position += 1
            return SampleBatch(row["counts"], row["backend"])
        if self.submitted:
            raise NextRunRequired(request)
        # Set this BEFORE submitting, so backend/protocol errors cannot trigger
        # a second submission in this invocation.
        self.submitted = True
        batch = self.sample(weights, c, penalty, parameters, context)
        summarize_counts(weights, c, penalty, batch.counts, context["shots"])
        self.history.append(json.loads(canonical({"request": request,
            "counts": dict(batch.counts), "backend": batch.metadata})))
        self.position += 1
        return batch


def run_one_job_step(weights, labels, target_size, mode, config, sample,
                     backend_label, resume_state=None):
    """Return a report and resume_state after zero or one new submissions.

    A final sample remains distinct from training, even at the same angles.
    Penalty reinforcement and recursive calls continue across run boundaries.
    Changing graph/config/backend label invalidates a checkpoint before any job.
    """
    if mode not in ("optimized_cut", "toy_tree"):
        raise ValueError("Only optimized_cut and toy_tree use checkpoints.")
    if config.optimizer == "COBYLA":
        raise ValueError("COBYLA is currently supported in multi_job local runs only; "
                         "use coordinate search for IonQ checkpoints.")
    contract = json.loads(canonical({"version": 1, "weights": weights,
        "labels": labels, "target_size": target_size, "mode": mode,
        "config": asdict(config), "backend_label": backend_label}))
    if resume_state is not None:
        saved_contract = dict(resume_state.get("contract", {}))
        # Preserve existing coordinate-search checkpoints from before the new
        # optimizer settings were added to QAOAConfig.
        saved_config = dict(saved_contract.get("config", {}))
        for key in ("optimizer", "cobyla_rhobeg", "cobyla_tol"):
            saved_config.setdefault(key, contract["config"][key])
        saved_contract["config"] = saved_config
        if canonical(saved_contract) != canonical(contract):
            raise ValueError("Graph, settings or backend label changed since this checkpoint.")
        if resume_state.get("status") == "failed":
            raise ValueError("The previous run failed. Inspect its job/error before resuming.")
        history = resume_state["history"]
    else:
        history = []
    replay = OneJobReplaySampler(sample, history)
    solver = QAOACutSolver(replay, config)
    report = {**contract, "status": "running"}
    try:
        if mode == "optimized_cut":
            bits = solver(weights, target_size)
            exact = exact_fixed_size_cut(weights, target_size)
            report.update({"bits": bits, "cut": cut_weight(weights, bits),
                           "ncut": normalized_cut(weights, bits),
                           "exact_cut": cut_weight(weights, exact),
                           "cut_gap": cut_weight(weights, bits) - cut_weight(weights, exact)})
        else:
            tree, trace = reconstruct(weights, labels, solver)
            exact_tree, exact_trace = reconstruct(weights, labels)
            report.update({"tree": tree, "trace": trace, "exact_tree": exact_tree,
                           "exact_trace": exact_trace})
        if replay.position != len(replay.history):
            raise ValueError("Checkpoint includes unused measurements.")
        report["status"] = "completed"
    except NextRunRequired as pause:
        report["status"] = "awaiting_next_run"
        report["next_request"] = pause.request
    except Exception as exc:
        report["status"] = "failed"
        report["error"] = f"{type(exc).__name__}: {exc}"
    report["evaluations"] = solver.evaluations
    report["attempts"] = solver.attempts
    report["new_submission_attempts"] = int(replay.submitted)
    report["completed_measurements"] = len(replay.history)
    report["resume_state"] = {"contract": contract, "history": replay.history,
                              "status": report["status"]}
    return report


def checkpoint_from_report(report):
    """Accept a current report, or recover successful jobs from the old run.

    The old multi-job report can contain a failed second submission; only
    completed evaluation batches are reused. No failed request is replayed.
    """
    if "resume_state" in report:
        return report["resume_state"]
    if report.get("mode") not in ("optimized_cut", "toy_tree"):
        raise ValueError("This output is not an optimizer checkpoint.")
    evaluations = report.get("evaluations", [])
    if not evaluations:
        raise ValueError("No completed optimizer measurements were saved.")
    target_size = report.get("target_size", evaluations[0]["target_size"])
    contract = {"version": 1, "weights": report["weights"],
                "labels": report["labels"], "target_size": target_size,
                "mode": report["mode"], "config": report["config"],
                "backend_label": report["backend_label"]}
    history = []
    for row in evaluations:
        if row["backend"].get("status") != "DONE":
            raise ValueError("An evaluation has no completed backend job.")
        context = {key: row[key] for key in (
            "call", "attempt", "evaluation", "phase", "shots", "seed")}
        history.append({"request": sample_request(row["weights"], row["target_size"],
            row["penalty"], row["parameters"], context), "counts": row["counts"],
            "backend": row["backend"]})
    error = report.get("error", "")
    if "only one job is allowed per code run" not in error:
        raise ValueError("Only the known one-job-limit error is safe to recover automatically.")
    return json.loads(canonical({"contract": contract, "history": history,
                               "status": "awaiting_next_run"}))

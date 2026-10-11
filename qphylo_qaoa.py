"""Small QAOA adapter for the graph-cut core, with no D-Wave dependencies.

Only build_qaoa and QiskitBackendSampler require Qiskit. The math, sample
scoring, and default coordinate search use standard Python; COBYLA needs SciPy.
This is a new adapter following the FO QAOA formulation, not copied FO code.
The default optimizer is seeded coordinate search; selectable COBYLA uses the
same SciPy algorithm underneath FO's Qiskit wrapper, with recorded pilot settings.
"""

from dataclasses import dataclass, field
from math import isfinite, pi
from numbers import Integral
import random
import time
from typing import Mapping

from qphylo_core import (
    cut_weight, normalized_cut, qubo_coefficients, qubo_energy,
    validate_bits, validate_weights,
)


def ising_coefficients(weights, c, penalty):
    """Convert x_i=(1-Z_i)/2; return offset, Z fields and ZZ couplings."""
    constant, linear, quadratic = qubo_coefficients(weights, c, penalty)
    offset = constant + sum(linear) / 2 + sum(quadratic.values()) / 4
    fields = [-value / 2 for value in linear]
    couplings = {}
    for (i, j), value in quadratic.items():
        fields[i] -= value / 4
        fields[j] -= value / 4
        couplings[i, j] = value / 4
    return offset, tuple(fields), couplings


def hamiltonian_scale(weights, c, penalty):
    """Rescale phases/objective by max absolute nonconstant coefficient."""
    _, fields, couplings = ising_coefficients(weights, c, penalty)
    return max(1.0, *(abs(v) for v in fields),
               *(abs(v) for v in couplings.values()))


def gate_plan(weights, c, penalty, parameters):
    """QAOA H|0>, then exp(-i gamma H/scale), exp(-i beta sum X).

    Parameters are [gamma_0,...,gamma_(p-1), beta_0,...,beta_(p-1)].
    Omit only the identity term (a global phase). Use CX-RZ-CX for ZZ so
    Qollab need not support an opaque RZZ instruction. No barriers/layout swaps.
    """
    w = validate_weights(weights)
    parameters = tuple(float(v) for v in parameters)
    if not parameters or len(parameters) % 2 or not all(map(isfinite, parameters)):
        raise ValueError("Provide finite gamma and beta values for every layer.")
    p = len(parameters) // 2
    _, fields, couplings = ising_coefficients(w, c, penalty)
    scale = hamiltonian_scale(w, c, penalty)
    plan = [("h", None, (i,)) for i in range(len(w))]
    for layer in range(p):
        gamma, beta = parameters[layer], parameters[p + layer]
        for i, coefficient in enumerate(fields):
            if coefficient:
                plan.append(("rz", 2 * gamma * coefficient / scale, (i,)))
        for (i, j), coefficient in couplings.items():
            if coefficient:
                plan.extend([("cx", None, (i, j)),
                             ("rz", 2 * gamma * coefficient / scale, (j,)),
                             ("cx", None, (i, j))])
        plan.extend(("rx", 2 * beta, (i,)) for i in range(len(w)))
    return plan


def build_qaoa(weights, c, penalty, parameters, measure=True):
    from qiskit import QuantumCircuit
    n = len(validate_weights(weights))
    qc = QuantumCircuit(n, n) if measure else QuantumCircuit(n)
    for gate, angle, qubits in gate_plan(weights, c, penalty, parameters):
        args = qubits if angle is None else (angle, *qubits)
        getattr(qc, gate)(*args)
    if measure:
        qc.measure(list(range(n)), list(range(n)))
    return qc


def decode_bitstring(bitstring, n):
    """get_counts strings are c_(n-1)..c_0; core bits are x_0..x_(n-1)."""
    if not isinstance(bitstring, str):
        raise ValueError("Expected a binary string from get_counts().")
    key = bitstring.replace(" ", "")
    if len(key) != n or any(char not in "01" for char in key):
        raise ValueError(f"Expected exactly {n} binary digits, got {bitstring!r}.")
    return tuple(int(char) for char in reversed(key))


def summarize_counts(weights, c, penalty, counts: Mapping, expected_shots=None):
    """Classically score observed samples; no projection or invented samples."""
    w = validate_weights(weights)
    qubo_energy(w, [0] * len(w), c, penalty)
    rows, total, feasible_shots = [], 0, 0
    for key, count in counts.items():
        if isinstance(count, bool) or not isinstance(count, Integral) or count < 0:
            raise ValueError("Counts must be nonnegative integer shot counts.")
        bits = decode_bitstring(key, len(w))
        if count == 0:
            continue
        count = int(count)
        selected = tuple(i for i, bit in enumerate(bits) if bit)
        feasible = len(selected) == c
        total += count
        feasible_shots += count if feasible else 0
        rows.append({"bits": bits, "bitstring": key, "count": count,
                     "selected": selected, "feasible": feasible,
                     "cut": cut_weight(w, bits),
                     "energy": qubo_energy(w, bits, c, penalty)})
    if total <= 0:
        raise ValueError("No measurements were returned.")
    if expected_shots is not None and total != expected_shots:
        raise ValueError(f"Expected {expected_shots} shots, got {total}.")
    # When total energies tie, prefer feasible, then first selected-index tuple.
    best = min(rows, key=lambda row: (row["energy"], not row["feasible"],
                                      row["selected"]))
    feasible_rows = [row for row in rows if row["feasible"]]
    best_feasible = min(feasible_rows,
                        key=lambda row: (row["cut"], row["selected"])) if feasible_rows else None
    return {"shots": total, "feasible_shots": feasible_shots,
            "feasible_rate": feasible_shots / total,
            "expectation": sum(row["count"] * row["energy"] for row in rows) / total,
            "best_energy_sample": best, "best_feasible_sample": best_feasible}


@dataclass(frozen=True)
class SampleBatch:
    counts: Mapping
    metadata: dict = field(default_factory=dict)


class QiskitBackendSampler:
    """Use the backend supplied by Qollab; never instantiate a provider/key."""
    def __init__(self, backend, poll_seconds=5, timeout_seconds=600, progress=print,
                 job_limit=None):
        if poll_seconds < 0 or timeout_seconds <= 0:
            raise ValueError("Polling interval and timeout must be valid.")
        self.backend = backend
        self.poll_seconds = poll_seconds
        self.timeout_seconds = timeout_seconds
        self.progress = progress
        if job_limit is not None and (isinstance(job_limit, bool) or
                not isinstance(job_limit, int) or job_limit < 1):
            raise ValueError("job_limit must be a positive integer or None.")
        self.job_limit = job_limit
        self.jobs = []

    def __call__(self, weights, c, penalty, parameters, context):
        if self.job_limit is not None and len(self.jobs) >= self.job_limit:
            raise RuntimeError("This run's job limit is reached; resume in a new code run.")
        import qiskit
        qc = build_qaoa(weights, c, penalty, parameters)
        shots = context["shots"]
        metadata = {"qiskit_version": qiskit.__version__,
                    "backend_type": type(self.backend).__name__,
                    "circuit_depth": qc.depth(),
                    "gate_counts": dict(qc.count_ops()), **context}
        options = getattr(self.backend, "options", None)
        # Store only explicit simulator settings, not arbitrary provider options.
        for name in ("noise_model", "seed_simulator"):
            value = getattr(options, name, None)
            if isinstance(value, (str, int, float, bool)) or value is None:
                metadata[name] = value
        started = time.monotonic()
        self.jobs.append(metadata)
        try:
            job = self.backend.run(qc, shots=shots)
            metadata["job_id"] = str(job.job_id())
            if self.progress:
                self.progress(f"Job {metadata['job_id']} ({context.get('phase', 'probe')})")
            last_status = None
            while True:
                status = job.status()
                name = getattr(status, "name", str(status).split(".")[-1])
                metadata["status"] = name
                if name != last_status and self.progress:
                    self.progress(f"  {name}")
                last_status = name
                if name in ("DONE", "ERROR", "CANCELLED"):
                    break
                if time.monotonic() - started >= self.timeout_seconds:
                    raise TimeoutError(
                        f"Stopped waiting for job {metadata['job_id']}; it may still "
                        "be running. Check Qollab before submitting again.")
                time.sleep(self.poll_seconds)
            if name != "DONE":
                raise RuntimeError(f"Job {metadata['job_id']} ended with {name}.")
            counts = job.result().get_counts()
            # Validate before returning; protocol errors are not penalty retries.
            summarize_counts(weights, c, penalty, counts, shots)
            return SampleBatch(counts, metadata)
        except Exception as exc:
            metadata["error"] = f"{type(exc).__name__}: {exc}"
            raise
        finally:
            metadata["elapsed_seconds"] = time.monotonic() - started


@dataclass(frozen=True)
class QAOAConfig:
    p: int = 1
    shots: int = 1000
    seed: int = 42
    max_evaluations: int = 9  # Includes one fresh final sample, per attempt.
    initial_penalty: float = 100.0
    penalty_step: float = 500.0
    max_retries: int = 1  # At most two attempts; no unbounded retries.
    coordinate_step: float = 0.6
    optimizer: str = "seeded_coordinate_search"
    cobyla_rhobeg: float = 1.0
    cobyla_tol: float = 1e-4

    def __post_init__(self):
        for name in ("p", "shots", "max_evaluations"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer.")
        if self.max_evaluations < 2:
            raise ValueError("Reserve at least one training and one final sample.")
        if isinstance(self.max_retries, bool) or not isinstance(self.max_retries, int) or self.max_retries < 0:
            raise ValueError("max_retries must be a nonnegative integer.")
        if isinstance(self.seed, bool) or not isinstance(self.seed, int):
            raise ValueError("seed must be an integer.")
        if self.optimizer not in ("seeded_coordinate_search", "COBYLA"):
            raise ValueError("optimizer must be seeded_coordinate_search or COBYLA.")
        for name in ("initial_penalty", "penalty_step", "coordinate_step",
                     "cobyla_rhobeg", "cobyla_tol"):
            value = getattr(self, name)
            if not isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive.")


class FeasibilityError(RuntimeError):
    pass


class EvaluationBudgetReached(RuntimeError):
    """Internal stop if an optimizer attempts to exceed its training budget."""


class QAOACutSolver:
    """Callable solve(weights,c), using an injected QAOA sample function.

    Optimize measured penalized expectation using coordinate search or COBYLA.
    Select from a fresh final sample at the best measured parameters. Reinforce
    alpha when its lowest observed QUBO-energy sample has the wrong cardinality,
    even if other feasible samples exist. Retry budget is explicit. This mirrors
    Onodera's feasibility reinforcement logic, not the Digital Annealer itself.
    """
    def __init__(self, sample, config=None):
        self.sample = sample
        self.config = config or QAOAConfig()
        self.evaluations = []
        self.attempts = []
        self.calls = 0

    def __call__(self, weights, c):
        w = validate_weights(weights)
        cfg = self.config
        minimize, scipy_version = None, None
        if cfg.optimizer == "COBYLA":
            # Lazy import: coordinate search still needs no SciPy installation.
            import scipy
            from scipy.optimize import minimize
            scipy_version = scipy.__version__
        qubo_energy(w, [0] * len(w), c, cfg.initial_penalty)
        call = self.calls
        self.calls += 1
        for attempt in range(cfg.max_retries + 1):
            penalty = cfg.initial_penalty + cfg.penalty_step * attempt
            scale = hamiltonian_scale(w, c, penalty)
            seed = cfg.seed + 1009 * call + attempt
            rng = random.Random(seed)
            parameters = tuple(rng.uniform(0, 2 * pi) for _ in range(cfg.p)) + tuple(
                rng.uniform(0, pi) for _ in range(cfg.p))
            evaluation = 0
            training_budget = cfg.max_evaluations - 1
            initial_parameters = parameters
            best_parameters, best_value = parameters, float("inf")
            optimizer_result = {}

            def evaluate(theta, phase):
                nonlocal evaluation, best_parameters, best_value
                theta = tuple(float(v) for v in theta)
                if phase == "training" and evaluation >= training_budget:
                    raise EvaluationBudgetReached("Training budget exhausted before submission.")
                context = {"call": call, "attempt": attempt,
                           "evaluation": evaluation, "phase": phase,
                           "shots": cfg.shots, "seed": seed}
                evaluation += 1
                batch = self.sample(w, c, penalty, theta, context)
                stats = summarize_counts(w, c, penalty, batch.counts, cfg.shots)
                self.evaluations.append({**context, "weights": w, "target_size": c,
                    "penalty": penalty, "scale": scale, "parameters": theta,
                    "optimizer": cfg.optimizer, "scipy_version": scipy_version,
                    "initial_parameters": initial_parameters,
                    "counts": dict(batch.counts), "summary": stats,
                    "backend": batch.metadata})
                value = stats["expectation"] / scale
                if phase == "training" and value < best_value:
                    best_value, best_parameters = value, theta
                return stats

            if cfg.optimizer == "COBYLA":
                bounds = [(0.0, 2 * pi)] * cfg.p + [(0.0, pi)] * cfg.p
                def cobyla_objective(theta):
                    # Older/newer SciPy versions can evaluate infeasible trial
                    # points while enforcing bounds. Submit the clipped angles
                    # so both optimizers actually sample the same angle domain.
                    bounded = tuple(min(upper, max(lower, float(value)))
                                    for value, (lower, upper) in zip(theta, bounds))
                    return evaluate(bounded, "training")["expectation"] / scale
                try:
                    result = minimize(
                        cobyla_objective,
                        initial_parameters, method="COBYLA", bounds=bounds,
                        options={"maxiter": training_budget, "rhobeg": cfg.cobyla_rhobeg,
                                 "tol": cfg.cobyla_tol})
                    optimizer_result = {"success": bool(result.success),
                        "status": int(result.status), "message": str(result.message),
                        "nfev": int(result.nfev), "returned_parameters": list(map(float, result.x))}
                except EvaluationBudgetReached:
                    optimizer_result = {"success": False, "message": "Explicit evaluation budget reached."}
                parameters = best_parameters
            else:
                evaluate(parameters, "training")
                step = cfg.coordinate_step
                while evaluation < training_budget:
                    improved, before = False, evaluation
                    for axis in range(2 * cfg.p):
                        for direction in (1, -1):
                            if evaluation >= training_budget:
                                break
                            proposal = list(parameters)
                            limit = 2 * pi if axis < cfg.p else pi
                            proposal[axis] = min(limit, max(0, proposal[axis] + direction * step))
                            proposal = tuple(proposal)
                            if proposal == parameters:
                                continue
                            previous_best = best_value
                            evaluate(proposal, "training")
                            if best_value < previous_best:
                                parameters, improved = best_parameters, True
                    if not improved:
                        step /= 2
                    if evaluation == before or step < 1e-6:
                        break

            final = evaluate(parameters, "final")
            best = final["best_energy_sample"]
            feasible = best["feasible"]
            self.attempts.append({"call": call, "attempt": attempt, "penalty": penalty,
                "target_size": c, "seed": seed, "parameters": parameters,
                "initial_parameters": initial_parameters, "optimizer": cfg.optimizer,
                "optimizer_result": optimizer_result, "scipy_version": scipy_version,
                "evaluations": evaluation, "feasible_rate": final["feasible_rate"],
                "best_energy_sample": best,
                "action": "accept" if feasible else (
                    "reinforce" if attempt < cfg.max_retries else "fail")})
            if feasible:
                bits = validate_bits(final["best_feasible_sample"]["bits"], len(w), c)
                return bits
        raise FeasibilityError(
            f"No acceptable final result for c={c} after {cfg.max_retries + 1} "
            "attempts. See solver.attempts and solver.evaluations; no exact fallback.")

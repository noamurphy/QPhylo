"""Mathematical/adapter checks, independent of quantum package installation."""

import cmath
from contextlib import redirect_stdout
from io import StringIO
import json
from tempfile import TemporaryDirectory
import shutil
from itertools import product
from math import cos, sin, sqrt
import importlib.util
from pathlib import Path
import sys
from types import ModuleType
import unittest
from unittest.mock import Mock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from qphylo_core import exact_fixed_size_cut, qubo_energy, reconstruct
from qphylo_qaoa import (
    FeasibilityError, QAOAConfig, QAOACutSolver, QiskitBackendSampler,
    SampleBatch, decode_bitstring, gate_plan, hamiltonian_scale,
    ising_coefficients, summarize_counts,
)

W = ((0, 9, 1, 2), (9, 0, 3, 1), (1, 3, 0, 8), (2, 1, 8, 0))


def apply_gate_plan(n, plan):
    """Tiny test-only gate simulator; does not replace Qiskit execution."""
    state = [0j] * (1 << n)
    state[0] = 1
    for name, angle, qubits in plan:
        if name == "cx":
            control, target = qubits
            new = [0j] * len(state)
            for index, amplitude in enumerate(state):
                destination = index ^ (1 << target) if index & (1 << control) else index
                new[destination] = amplitude
            state = new
            continue
        qubit = qubits[0]
        if name == "h":
            a, b, c, d = 1 / sqrt(2), 1 / sqrt(2), 1 / sqrt(2), -1 / sqrt(2)
        elif name == "rz":
            a, b, c, d = cmath.exp(-0.5j * angle), 0, 0, cmath.exp(0.5j * angle)
        elif name == "rx":
            a, b, c, d = cos(angle / 2), -1j * sin(angle / 2), -1j * sin(angle / 2), cos(angle / 2)
        else:
            raise AssertionError(name)
        new = state[:]
        for zero in range(len(state)):
            if zero & (1 << qubit):
                continue
            one = zero | (1 << qubit)
            new[zero] = a * state[zero] + b * state[one]
            new[one] = c * state[zero] + d * state[one]
        state = new
    return state


def direct_qaoa_state(weights, c, penalty, parameters):
    """Independent diagonal-QUBO phases and full tensor-product mixer."""
    n, p = len(weights), len(parameters) // 2
    scale = hamiltonian_scale(weights, c, penalty)
    offset = ising_coefficients(weights, c, penalty)[0]
    energies = [qubo_energy(weights, tuple((index >> i) & 1 for i in range(n)),
                            c, penalty) for index in range(1 << n)]
    state = [1 / sqrt(1 << n)] * (1 << n)
    for layer in range(p):
        gamma, beta = parameters[layer], parameters[p + layer]
        state = [amplitude * cmath.exp(-1j * gamma * (energy - offset) / scale)
                 for amplitude, energy in zip(state, energies)]
        state = [sum(state[source] * cos(beta) ** (n - (source ^ target).bit_count())
                     * (-1j * sin(beta)) ** (source ^ target).bit_count()
                     for source in range(1 << n)) for target in range(1 << n)]
    return state


class MathTests(unittest.TestCase):
    def test_ising_preserves_all_assignment_energies(self):
        for c in (1, 2, 3):
            for penalty in (0.125, 100, 600):
                offset, fields, pairs = ising_coefficients(W, c, penalty)
                for bits in product((0, 1), repeat=4):
                    spins = [1 - 2 * bit for bit in bits]
                    energy = offset + sum(a * z for a, z in zip(fields, spins))
                    energy += sum(value * spins[i] * spins[j] for (i, j), value in pairs.items())
                    self.assertAlmostEqual(energy, qubo_energy(W, bits, c, penalty))

    def test_gate_plan_matches_independent_diagonal_and_mixer(self):
        for c in (1, 2):
            for parameters in ((0.37, 0.41), (0.37, 1.13, 0.41, 0.19)):
                expected = direct_qaoa_state(W, c, 100, parameters)
                actual = apply_gate_plan(4, gate_plan(W, c, 100, parameters))
                for a, b in zip(actual, expected):
                    self.assertAlmostEqual(abs(a - b), 0, places=12)
                self.assertAlmostEqual(sum(abs(a) ** 2 for a in actual), 1)

    def test_asymmetric_count_keys_use_qubit_index_order(self):
        self.assertEqual(decode_bitstring("001", 3), (1, 0, 0))
        self.assertEqual(decode_bitstring("100", 3), (0, 0, 1))
        stats = summarize_counts(W, 1, 100, {"1000": 800, "0001": 200}, 1000)
        self.assertEqual(stats["best_feasible_sample"]["bits"], (0, 0, 0, 1))
        self.assertAlmostEqual(stats["expectation"], 11.2)

    def test_malformed_counts_are_rejected(self):
        for counts in ({}, {"000": 1}, {"xxxx": 1}, {"0000": -1},
                       {"0000": 0.5}, {"0000": True}, {0: 1}):
            with self.subTest(counts=counts), self.assertRaises(ValueError):
                summarize_counts(W, 2, 100, counts)
        with self.assertRaises(ValueError):
            summarize_counts(W, 2, 100, {"0011": 999}, 1000)

    def test_scoring_keeps_cut_separate_from_penalty(self):
        stats = summarize_counts(W, 2, 100, {"0000": 900, "0011": 100})
        self.assertAlmostEqual(stats["feasible_rate"], 0.1)
        self.assertEqual(stats["best_feasible_sample"]["cut"], 7)
        self.assertEqual(stats["best_energy_sample"]["bits"], (1, 1, 0, 0))


class SolverTests(unittest.TestCase):
    def test_reinforce_even_when_some_feasible_samples_exist(self):
        def sample(w, c, alpha, theta, context):
            return SampleBatch({"00": 400, "01": 600})
        solver = QAOACutSolver(sample, QAOAConfig(
            max_evaluations=2, initial_penalty=0.1, max_retries=1))
        self.assertEqual(solver(((0, 2), (2, 0)), 1), (1, 0))
        self.assertEqual([r["action"] for r in solver.attempts], ["reinforce", "accept"])
        self.assertEqual([r["penalty"] for r in solver.attempts], [0.1, 500.1])
        self.assertEqual(len(solver.evaluations), 4)
        self.assertEqual([r["phase"] for r in solver.evaluations],
                         ["training", "final", "training", "final"])

    def test_fail_at_retry_limit_without_fallback(self):
        sampler = Mock(return_value=SampleBatch({"0000": 1000}))
        solver = QAOACutSolver(sampler, QAOAConfig(max_evaluations=5, max_retries=1))
        with self.assertRaises(FeasibilityError):
            solver(W, 2)
        self.assertLessEqual(sampler.call_count, 10)
        self.assertEqual(solver.attempts[-1]["action"], "fail")

    def test_final_measurement_is_required(self):
        def sample(w, c, alpha, theta, context):
            key = "0011" if context["phase"] == "training" else "0000"
            return SampleBatch({key: 1000})
        solver = QAOACutSolver(sample, QAOAConfig(max_evaluations=2, max_retries=0))
        with self.assertRaises(FeasibilityError):
            solver(W, 2)

    def test_recursive_contract_with_explicit_classical_test_double(self):
        def exact_test_double(w, c, alpha, theta, context):
            bits = exact_fixed_size_cut(w, c)
            return SampleBatch({"".join(str(bit) for bit in reversed(bits)): 1000})
        solver = QAOACutSolver(exact_test_double, QAOAConfig(max_evaluations=2))
        self.assertEqual(reconstruct(W, list("ABCD"), solver), reconstruct(W, list("ABCD")))

    def test_optimizer_budget_and_seed_are_reproducible(self):
        def sample(w, c, alpha, theta, context):
            return SampleBatch({"0011": 1000})
        solvers = [QAOACutSolver(sample, QAOAConfig(max_evaluations=9, seed=42))
                   for _ in range(2)]
        for solver in solvers:
            solver(W, 2)
            self.assertLessEqual(len(solver.evaluations), 9)
            self.assertEqual(solver.evaluations[-1]["phase"], "final")
        self.assertEqual(solvers[0].evaluations, solvers[1].evaluations)


class BackendProtocolTests(unittest.TestCase):
    def test_done_and_error_jobs(self):
        qiskit_stub = ModuleType("qiskit")
        qiskit_stub.__version__ = "test-double"
        backend = Mock()
        backend.options = None
        job = backend.run.return_value
        job.job_id.return_value = "test-job"
        job.status.return_value = type("Status", (), {"name": "DONE"})()
        job.result.return_value.get_counts.return_value = {"0011": 1000}
        circuit = Mock()
        circuit.depth.return_value = 3
        circuit.count_ops.return_value = {"cx": 2}
        sampler = QiskitBackendSampler(backend, progress=None)
        with patch.dict(sys.modules, {"qiskit": qiskit_stub}), patch(
                "qphylo_qaoa.build_qaoa", return_value=circuit):
            batch = sampler(W, 2, 100, (0.3, 0.4), {"shots": 1000})
            self.assertEqual(batch.counts, {"0011": 1000})
            self.assertEqual(batch.metadata["job_id"], "test-job")
            job.reset_mock()
            job.status.return_value = type("Status", (), {"name": "ERROR"})()
            with self.assertRaises(RuntimeError):
                sampler(W, 2, 100, (0.3, 0.4), {"shots": 1000})
            job.result.assert_not_called()
            self.assertIn("error", sampler.jobs[-1])

    def test_explicit_job_limit_stops_second_submission_locally(self):
        sampler = QiskitBackendSampler(Mock(), job_limit=1, progress=None)
        sampler.jobs.append({"status": "DONE"})
        with self.assertRaises(RuntimeError):
            sampler(W, 2, 100, (0.3, 0.4), {"shots": 1000})
        sampler.backend.run.assert_not_called()

    def test_timeout_retains_job_identity(self):
        qiskit_stub = ModuleType("qiskit")
        qiskit_stub.__version__ = "test-double"
        backend, circuit = Mock(), Mock()
        backend.options = None
        backend.run.return_value.job_id.return_value = "still-running"
        backend.run.return_value.status.return_value = type("Status", (), {"name": "RUNNING"})()
        circuit.count_ops.return_value = {}
        circuit.depth.return_value = 1
        sampler = QiskitBackendSampler(backend, timeout_seconds=10, progress=None)
        with patch.dict(sys.modules, {"qiskit": qiskit_stub}), patch(
                "qphylo_qaoa.build_qaoa", return_value=circuit), patch(
                "qphylo_qaoa.time.monotonic", side_effect=[0, 11, 12]):
            with self.assertRaises(TimeoutError):
                sampler(W, 2, 100, (0.3, 0.4), {"shots": 1000})
        self.assertEqual(sampler.jobs[-1]["job_id"], "still-running")
        backend.run.return_value.result.assert_not_called()


class QollabPackagingTests(unittest.TestCase):
    def test_multi_job_bundle_completes_cut_and_tree_on_selected_backend(self):
        source = (Path(__file__).resolve().parents[1] /
                  "qollab_smoke/qaoa_cut_standalone.py").read_text()
        qiskit_stub = ModuleType("qiskit")
        qiskit_stub.__version__ = "test-double"
        circuit = Mock()
        circuit.depth.return_value = 18
        circuit.count_ops.return_value = {"h": 4}
        qiskit_stub.QuantumCircuit = Mock(return_value=circuit)
        for mode in ("optimized_cut", "toy_tree"):
            with self.subTest(mode=mode):
                backend = Mock()
                backend.options = None
                backend.run.return_value.job_id.return_value = "local-test"
                backend.run.return_value.status.return_value = type(
                    "Status", (), {"name": "DONE"})()
                # Include feasible outcomes for both root target sizes. The
                # energy scorer must choose the appropriate one for each cut.
                backend.run.return_value.result.return_value.get_counts.return_value = {
                    "0001": 500, "0011": 500}
                configured = source.replace('mode = "fixed_probe"', f'mode = "{mode}"').replace(
                    'execution_mode = "one_job"', 'execution_mode = "multi_job"')
                namespace = {"backend": backend, "__name__": "__main__"}
                with patch.dict(sys.modules, {"qiskit": qiskit_stub}), redirect_stdout(StringIO()):
                    exec(compile(configured, "qaoa_cut_standalone.py", "exec"), namespace)
                report = namespace["report"]
                self.assertEqual(report["status"], "completed")
                self.assertEqual(report["execution_mode"], "multi_job")
                self.assertGreater(backend.run.call_count, 1)
                self.assertEqual(backend.run.call_count, len(report["evaluations"]))
                self.assertNotIn("resume_state", report)
                if mode == "optimized_cut":
                    self.assertEqual(report["cut_gap"], 0)
                    self.assertEqual(report["evaluations"][-1]["phase"], "final")
                else:
                    self.assertEqual(report["tree"], "(('A','B'),('C','D'));")
                    self.assertEqual(report["tree"], report["exact_tree"])
                    self.assertEqual(len(report["attempts"]), 2)

    def test_multi_job_mode_rejects_checkpoint_before_submission(self):
        source = (Path(__file__).resolve().parents[1] /
                  "qollab_smoke/qaoa_cut_standalone.py").read_text().replace(
                    'execution_mode = "one_job"', 'execution_mode = "multi_job"').replace(
                    'resume_json = ""', 'resume_json = "{}"')
        backend = Mock()
        with redirect_stdout(StringIO()), self.assertRaisesRegex(ValueError, "checkpoint"):
            exec(compile(source, "qaoa_cut_standalone.py", "exec"),
                 {"backend": backend, "__name__": "__main__"})
        backend.run.assert_not_called()

    def test_standalone_default_submits_one_probe(self):
        # This checks packaging/control flow with explicit test doubles. Gate
        # mathematics is checked separately above; this is not a cloud run.
        source = (Path(__file__).resolve().parents[1] /
                  "qollab_smoke/qaoa_cut_standalone.py").read_text()
        qiskit_stub = ModuleType("qiskit")
        qiskit_stub.__version__ = "test-double"
        circuit = Mock()
        circuit.depth.return_value = 42
        circuit.count_ops.return_value = {"h": 4, "cx": 12, "rz": 6, "rx": 4, "measure": 4}
        qiskit_stub.QuantumCircuit = Mock(return_value=circuit)
        backend = Mock()
        backend.options = None
        backend.run.return_value.job_id.return_value = "packaging-test"
        backend.run.return_value.status.return_value = type("Status", (), {"name": "DONE"})()
        backend.run.return_value.result.return_value.get_counts.return_value = {
            "0011": 500, "1100": 500}
        namespace = {"backend": backend, "__name__": "__main__"}
        output = StringIO()
        with patch.dict(sys.modules, {"qiskit": qiskit_stub}), redirect_stdout(output):
            exec(compile(source, "qaoa_cut_standalone.py", "exec"), namespace)
        backend.run.assert_called_once()
        report_line = next(line for line in output.getvalue().splitlines()
                           if line.startswith("REPORT_JSON: "))
        report = json.loads(report_line.removeprefix("REPORT_JSON: "))
        self.assertEqual(report["mode"], "fixed_probe")
        self.assertEqual(report["status"], "probe_completed")
        self.assertEqual(report["exact_cut"], 4)
        self.assertEqual(report["observed_cut_gap"], 0)
        self.assertEqual(len(report["jobs"]), 1)

    def test_optimizer_bundle_resumes_from_saved_output_with_one_new_job(self):
        import prepare_qollab
        root = Path(__file__).resolve().parents[1]
        qiskit_stub = ModuleType("qiskit")
        qiskit_stub.__version__ = "test-double"
        circuit = Mock()
        circuit.depth.return_value = 1
        circuit.count_ops.return_value = {"h": 4}
        qiskit_stub.QuantumCircuit = Mock(return_value=circuit)

        def execute(source):
            backend = Mock()
            backend.options = None
            backend.run.return_value.job_id.return_value = "packaging-step"
            backend.run.return_value.status.return_value = type("Status", (), {"name": "DONE"})()
            backend.run.return_value.result.return_value.get_counts.return_value = {"0011": 1000}
            namespace = {"backend": backend, "__name__": "__main__"}
            output = StringIO()
            with patch.dict(sys.modules, {"qiskit": qiskit_stub}), redirect_stdout(output):
                exec(compile(source, "qaoa_cut_standalone.py", "exec"), namespace)
            backend.run.assert_called_once()
            return namespace["report"], output.getvalue()

        with TemporaryDirectory() as folder:
            temporary_root = Path(folder)
            (temporary_root / "qollab_smoke").mkdir()
            for relative in ("qphylo_core.py", "qphylo_qaoa.py", "qphylo_run_steps.py",
                             "qphylo_benchmark.py",
                             "qollab_smoke/cut_pilot.py"):
                shutil.copy2(root / relative, temporary_root / relative)
            with patch.object(prepare_qollab, "ROOT", temporary_root):
                initial = prepare_qollab.build().read_text().replace(
                    'mode = "fixed_probe"', 'mode = "optimized_cut"')
                first, output = execute(initial)
                self.assertEqual(first["status"], "awaiting_next_run")
                self.assertEqual(first["completed_measurements"], 1)
                saved = temporary_root / "first.txt"
                saved.write_text(output)
                resumed = prepare_qollab.build(saved).read_text()
                second, _ = execute(resumed)
                self.assertEqual(second["status"], "awaiting_next_run")
                self.assertEqual(second["completed_measurements"], 2)
                self.assertEqual(second["resume_state"]["history"][0],
                                 first["resume_state"]["history"][0])


@unittest.skipUnless(importlib.util.find_spec("qiskit"), "Qiskit is not installed locally")
class ActualQiskitTests(unittest.TestCase):
    def test_statevector_matches_independent_formula(self):
        from qiskit.quantum_info import Statevector
        from qphylo_qaoa import build_qaoa
        for c in (1, 2):
            parameters = (0.37, 1.13, 0.41, 0.19)
            actual = Statevector.from_instruction(build_qaoa(W, c, 100, parameters, measure=False))
            for a, b in zip(actual.data, direct_qaoa_state(W, c, 100, parameters)):
                self.assertAlmostEqual(abs(a - b), 0, places=12)


if __name__ == "__main__":
    unittest.main()

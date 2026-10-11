"""Submission scheduler tests use explicit circuit/backend/optimizer doubles.

These validate integration and job limits, not new measured quantum results.
"""

from contextlib import redirect_stdout
from io import StringIO
import json
from pathlib import Path
import shutil
import sys
from tempfile import TemporaryDirectory
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import Mock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from prepare_submission import build_submission, ROOT
from qphylo_core import exact_fixed_size_cut
from collect_qaoa_results import validated_report


class BasicSimulator:
    def __init__(self, counts=None, terminal="DONE"):
        self.calls = 0
        self.options = SimpleNamespace(noise_model=None, seed_simulator=None)
        self.counts, self.terminal = counts, terminal

    def set_options(self, **options):
        self.options.__dict__.update(options)

    def run(self, qc, shots):
        self.calls += 1
        job = Mock()
        job.job_id.return_value = str(self.calls)
        job.status.return_value = SimpleNamespace(name=self.terminal)
        job.result.return_value.get_counts.return_value = self.counts or {qc.test_bits: shots}
        return job


class IonQSimulatorBackend(BasicSimulator):
    pass


class UnknownSimulator(BasicSimulator):
    pass


def dependencies():
    qiskit, scipy, optimize = (ModuleType(name) for name in ("qiskit", "scipy", "scipy.optimize"))
    qiskit.__version__ = scipy.__version__ = "test-double"
    def circuit(n, m):
        qc = Mock()
        qc.depth.return_value = 18
        qc.count_ops.return_value = {"h": n}
        return qc
    qiskit.QuantumCircuit = circuit
    def minimize(fun, x0, method, bounds, options):
        for _ in range(options["maxiter"]):
            fun(x0)
        return SimpleNamespace(success=False, status=2, message="budget", nfev=options["maxiter"], x=x0)
    optimize.minimize = minimize
    return {"qiskit": qiskit, "scipy": scipy, "scipy.optimize": optimize}


class SubmissionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.folder = TemporaryDirectory()
        cls.bundle = build_submission(output=Path(cls.folder.name)/"main.py")
        cls.source = cls.bundle.read_text()

    @classmethod
    def tearDownClass(cls):
        cls.folder.cleanup()

    def execute(self, backend, mode="auto", change_data=None):
        modules = dependencies()
        output, namespace = StringIO(), {"backend": backend, "__name__": "__main__"}
        prefix, entry = self.source.split("# ---- Source: qollab_smoke/submission_entry.py ----", 1)
        with patch.dict(sys.modules, modules), redirect_stdout(output):
            exec(compile(prefix, "definitions", "exec"), namespace)
            namespace["RUN_MODE"] = mode
            if change_data:
                change_data(namespace["SUBMISSION_DATA"])
            builder = namespace["build_qaoa"]
            def test_builder(w, c, alpha, parameters):
                qc = builder(w, c, alpha, parameters)
                qc.test_bits = "".join(map(str, reversed(exact_fixed_size_cut(w,c))))
                return qc
            namespace["build_qaoa"] = test_builder
            try:
                exec(compile(entry, "submission", "exec"), namespace)
                error = None
            except Exception as exc:
                error = exc
        text = output.getvalue()
        line = next(line for line in text.splitlines() if line.startswith("REPORT_JSON:"))
        return json.loads(line.split(":",1)[1]), text, error

    def test_auto_single_job_on_IonQ_and_unknown_backend(self):
        noisy_basic = BasicSimulator()
        noisy_basic.options.noise_model = "test-noise"
        for backend in (IonQSimulatorBackend(), UnknownSimulator(), noisy_basic):
            report, text, error = self.execute(backend)
            self.assertIsNone(error)
            self.assertEqual(backend.calls, 1)
            self.assertEqual(report["status"], "probe_completed")
            self.assertEqual(report["selected_run_mode"], "single_cut")
            self.assertEqual(report["current_job_limit"], 1)
            self.assertEqual(report["summary"]["shots"], 1000)
            self.assertNotIn("tree", report)
            self.assertEqual(report["angle_source"]["source_file"], "protein1822_cobyla_run1.html")
            self.assertIn("SAVED COMPARISON - previous runs", text)
            self.assertIn("does not optimize angles or reconstruct a new tree", text)
            # Its one current job is not conflated with 403 historical jobs.
            self.assertEqual(report["resource_metrics"]["submitted_jobs"], 1)
            self.assertEqual(report["saved_comparison"]["submitted_jobs"], 403)

    def test_auto_full_tree_is_collector_compatible(self):
        backend = BasicSimulator()
        report, text, error = self.execute(backend)
        self.assertIsNone(error)
        self.assertEqual(report["status"], "completed")
        self.assertEqual(report["mode"], "protein_tree")
        self.assertEqual(backend.calls, 279)  # Nine accepted cuts, 31 calls each, no retry in this double.
        self.assertEqual(report["topology_metrics"]["unrooted_rf_to_exact"], 0)
        self.assertEqual(report["baseline_metrics"]["unrooted_rf_to_saved_model_nj"], 2)
        self.assertEqual(report["config"]["optimizer"], "COBYLA")
        self.assertEqual(backend.options.seed_simulator, 20261052)
        self.assertIn("Fresh QAOA grouping", text)
        with TemporaryDirectory() as folder:
            path = Path(folder)/"double_report.json"
            path.write_text(json.dumps(report))
            self.assertEqual(validated_report(path)["resource_metrics"]["submitted_jobs"], 279)

    def test_forced_full_tree_on_IonQ_rejected_before_submission(self):
        backend = IonQSimulatorBackend()
        report, text, error = self.execute(backend, "full_tree")
        self.assertIsInstance(error, ValueError)
        self.assertEqual(backend.calls, 0)
        self.assertEqual(report["status"], "failed")

    def test_single_cut_no_feasible_sample_kept_without_retry_or_tree(self):
        backend = IonQSimulatorBackend(counts={"00000000":1000})
        report, text, error = self.execute(backend)
        self.assertIsNone(error)
        self.assertEqual(backend.calls, 1)
        self.assertEqual(report["summary"]["feasible_rate"], 0)
        self.assertIsNone(report["summary"]["best_feasible_sample"])
        self.assertNotIn("tree", report)
        self.assertNotIn("observed_cut_gap", report)
        self.assertIn("No feasible sample observed", text)

    def test_backend_failure_preserves_one_failed_job_and_stops(self):
        backend = IonQSimulatorBackend(terminal="ERROR")
        report, text, error = self.execute(backend)
        self.assertIsInstance(error, RuntimeError)
        self.assertEqual(backend.calls, 1)
        self.assertEqual(report["status"], "failed")
        self.assertEqual(report["jobs"][0]["status"], "ERROR")
        self.assertNotIn("counts", report)

    def test_frozen_angle_graph_or_scale_mismatch_rejected_before_job(self):
        for field, value in (("graph_signature", "changed"), ("scale", 99)):
            backend = IonQSimulatorBackend()
            report, _, error = self.execute(backend, change_data=lambda d:d["frozen_cut"].update({field:value}))
            self.assertIsInstance(error, ValueError)
            self.assertEqual(backend.calls, 0)

    def test_builder_rejects_altered_comparison_snapshot(self):
        with TemporaryDirectory() as folder:
            comparison = Path(folder)/"comparison"
            shutil.copytree(ROOT/"results/protein1822_comparison", comparison)
            with (comparison/"splits.csv").open("a") as stream:
                stream.write("\nchanged\n")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                build_submission(output=Path(folder)/"main.py", comparison=comparison)


if __name__ == "__main__":
    unittest.main()

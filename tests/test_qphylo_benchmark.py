"""Exact fixture, optimizer-budget and saved-report comparison checks."""

from contextlib import redirect_stdout
from copy import deepcopy
from dataclasses import asdict
from io import StringIO
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import Mock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from qphylo_benchmark import (benchmark_metrics, job_upper_bound, topology_splits,
                              toy_fixture)
from qphylo_core import exact_fixed_size_cut, reconstruct
from qphylo_qaoa import QAOACutSolver, QAOAConfig, SampleBatch
from qphylo_run_steps import run_one_job_step
from collect_qaoa_results import comparison_tables, read_report, run_inventory, export_tables
from prepare_qollab import build, protein_graph


def scipy_stub(overrun=False):
    scipy, optimize = ModuleType("scipy"), ModuleType("scipy.optimize")
    scipy.__version__ = "test-double"
    def minimize(fun, x0, method, bounds, options):
        assert method == "COBYLA"
        budget = options["maxiter"]
        for step in range(budget + int(overrun)):
            theta = [v + step * 0.1 for v in x0]
            fun(theta)
        return SimpleNamespace(success=False, status=2, message="budget", nfev=budget,
                               x=[upper + 10 for lower, upper in bounds])
    optimize.minimize = minimize
    return {"scipy": scipy, "scipy.optimize": optimize}


class BenchmarkTests(unittest.TestCase):
    def test_eight_taxon_exact_tree_and_resource_bound(self):
        weights, labels = toy_fixture("hierarchical8")
        tree, trace = reconstruct(weights, labels)
        self.assertEqual(tree, "((('A','B'),('C','D')),(('E','F'),('G','H')));")
        self.assertEqual([len(row["taxa"]) for row in trace], [8, 4, 4])
        self.assertEqual([row["cut"] for row in trace], [16, 12, 12])
        self.assertEqual(len(topology_splits(trace, labels)), 5)
        self.assertEqual(job_upper_bound(8, QAOAConfig()), 270)
        self.assertEqual(job_upper_bound(8, QAOAConfig(max_evaluations=31)), 930)

    def test_unrooted_splits_ignore_child_order_and_detect_wrong_split(self):
        _, labels = toy_fixture("toy4")
        correct = [{"left": ["A", "B"], "right": ["C", "D"]}]
        reversed_order = [{"right": ["A", "B"], "left": ["C", "D"]}]
        wrong = [{"left": ["A", "C"], "right": ["B", "D"]}]
        self.assertEqual(topology_splits(correct, labels), topology_splits(reversed_order, labels))
        self.assertEqual(len(topology_splits(correct, labels) ^ topology_splits(wrong, labels)), 2)

    def test_cobyla_budget_guard_and_fresh_best_measured_sample(self):
        weights, _ = toy_fixture("toy4")
        for overrun in (False, True):
            def sample(w, c, alpha, theta, context):
                # Later measured evaluations improve; fake minimizer returns an
                # unmeasured/out-of-bounds point, which must not be selected.
                feasible = min(1000, 100 * (context["evaluation"] + 1))
                return SampleBatch({"0011": feasible, "0000": 1000-feasible})
            sampler = Mock(side_effect=sample)
            cfg = QAOAConfig(optimizer="COBYLA", max_evaluations=6, max_retries=0)
            with patch.dict(sys.modules, scipy_stub(overrun)):
                solver = QAOACutSolver(sampler, cfg)
                solver(weights, 2)
            self.assertEqual(sampler.call_count, 6)
            self.assertEqual(solver.evaluations[-1]["phase"], "final")
            self.assertEqual(solver.evaluations[-1]["parameters"], solver.evaluations[-2]["parameters"])
            self.assertEqual(solver.attempts[0]["scipy_version"], "test-double")
            for row in solver.evaluations:
                self.assertTrue(0 <= row["parameters"][0] <= 2*3.141592653589793)
                self.assertTrue(0 <= row["parameters"][1] <= 3.141592653589793)

    def test_cobyla_checkpoint_rejected_before_sampling(self):
        w, labels = toy_fixture("toy4")
        sampler = Mock()
        with self.assertRaisesRegex(ValueError, "multi_job"):
            run_one_job_step(w, labels, 2, "toy_tree", QAOAConfig(optimizer="COBYLA"),
                             sampler, "test")
        sampler.assert_not_called()

    def test_presets_execute_eight_taxon_recursion(self):
        root = Path(__file__).resolve().parents[1]
        qiskit = ModuleType("qiskit")
        qiskit.__version__ = "test-double"
        def circuit(n, m):
            qc = Mock()
            qc.num_qubits = n
            qc.depth.return_value = 18
            qc.count_ops.return_value = {"h": n}
            return qc
        qiskit.QuantumCircuit = circuit
        for preset in ("eight_smoke", "eight_cobyla"):
            backend = Mock()
            backend.options = None
            def run(qc, shots):
                keys = ("00000001", "00000011", "00000111", "00001111") if qc.num_qubits == 8 else ("0001", "0011")
                job = Mock()
                job.job_id.return_value = str(backend.run.call_count)
                job.status.return_value = SimpleNamespace(name="DONE")
                job.result.return_value.get_counts.return_value = {k: shots//len(keys) for k in keys}
                return job
            backend.run.side_effect = run
            ns = {"backend": backend, "__name__": "__main__"}
            with patch.dict(sys.modules, {"qiskit": qiskit, **scipy_stub()}), redirect_stdout(StringIO()):
                exec(compile((root/"qollab_smoke"/f"{preset}.py").read_text(), preset, "exec"), ns)
            report = ns["report"]
            self.assertEqual(report["status"], "completed")
            self.assertEqual(report["topology_metrics"]["unrooted_rf_to_exact"], 0)
            self.assertEqual(report["topology_metrics"]["split_nodes"], 3)
            self.assertEqual(len(report["attempts"]), 8)
            self.assertEqual([row["taxa"] for row in report["final_cut_metrics"]][-1], tuple("EFGH"))
            self.assertEqual(backend.run.call_count, 8*report["config"]["max_evaluations"])

    def test_saved_report_validation_and_matching(self):
        root = Path(__file__).resolve().parents[1]
        actual = root/"qollab_smoke/qaoa_toyTestResults_noamurphy-qaoa-phylogenetic-tree-reconst-output-completed-2026-10-10T19-58-08.397Z.html"
        report = read_report(actual)
        with TemporaryDirectory() as folder:
            path = Path(folder)/"copy.json"
            path.write_text(json.dumps(report))
            reports, tables, pairs = comparison_tables([path])
            self.assertEqual(tables["runs"][0]["rf_to_exact"], 0)
            self.assertEqual(tables["runs"][0]["submitted_jobs"], 18)
            self.assertEqual(len(tables["evaluations"]), 18)
            self.assertEqual(len(tables["final_cuts"]), 2)
            self.assertEqual(reports[0]["final_cut_metrics"][1]["optimal_shot_rate"], 0.088)
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                comparison_tables([path, path])
            misnamed = Path(folder)/"coordinate_seed43.json"
            misnamed.write_text(json.dumps(report))
            self.assertEqual(run_inventory([misnamed])[0]["stored_seed"], 42)
            with self.assertRaises(ValueError) as exc:
                comparison_tables([path, misnamed])
            self.assertIn("coordinate_seed43.json", str(exc.exception))
            self.assertIn("copy.json", str(exc.exception))
            self.assertIn("not filenames", str(exc.exception))
            mismatched_seed = deepcopy(report)
            mismatched_seed["simulator_seed"] = 1234
            path.write_text(json.dumps(mismatched_seed))
            with self.assertRaisesRegex(ValueError, "simulator seed"):
                comparison_tables([path])
            corrupted = deepcopy(report)
            # Independent runtimes may differ in the last floating-point bits.
            corrupted["evaluations"][0]["summary"]["expectation"] += 1e-13
            path.write_text(json.dumps(corrupted))
            comparison_tables([path])
            corrupted = deepcopy(report)
            corrupted["evaluations"][0]["summary"]["shots"] += 1
            path.write_text(json.dumps(corrupted))
            with self.assertRaisesRegex(ValueError, "disagree"):
                comparison_tables([path])
            corrupted = deepcopy(report)
            corrupted["evaluations"][0]["summary"]["expectation"] += 1
            path.write_text(json.dumps(corrupted))
            with self.assertRaisesRegex(ValueError, "disagree"):
                comparison_tables([path])

    def test_separate_seed_scripts_override_dataclass_default(self):
        import ast
        with TemporaryDirectory() as folder:
            paths = [build(preset="eight_coordinate", seed=seed,
                           output_path=Path(folder)/f"coordinate_seed{seed}.py")
                     for seed in (43, 44)]
            for seed, path in zip((43, 44), paths):
                # Inspect actual top-level settings, not the unrelated default
                # seed in the embedded QAOAConfig class.
                source = path.read_text()
                node = next(node for node in ast.parse(source).body
                            if isinstance(node, ast.Assign) and
                            isinstance(node.targets[0], ast.Name) and node.targets[0].id == "seed")
                self.assertEqual(ast.literal_eval(node.value), seed)
                self.assertIn(f"initialization seed: {seed}.", source)
            self.assertTrue(all(path.exists() for path in paths))

    def test_protein_bundle_executes_embedded_graph_and_exports_provenance(self):
        # Explicit circuit/backend/minimizer doubles test integration, not QAOA
        # accuracy. No cloud jobs or invented measured result files are created.
        root = Path(__file__).resolve().parents[1]
        graph_path = root / "results/protein1822/graph.json"
        graph = json.loads(graph_path.read_text())
        qiskit = ModuleType("qiskit")
        qiskit.__version__ = "test-double"
        def circuit(n, m):
            qc = Mock()
            qc.depth.return_value = 18
            qc.count_ops.return_value = {"h": n}
            return qc
        qiskit.QuantumCircuit = circuit
        class BasicSimulator:
            def __init__(self):
                self.options = SimpleNamespace(noise_model=None, seed_simulator=None)
                self.calls = 0
            def set_options(self, **kwargs):
                self.options.__dict__.update(kwargs)
            def run(self, qc, shots):
                self.calls += 1
                job = Mock()
                job.job_id.return_value = str(self.calls)
                job.status.return_value = SimpleNamespace(name="DONE")
                job.result.return_value.get_counts.return_value = {qc.test_bits: shots}
                return job
        with TemporaryDirectory() as folder:
            path = build(graph_path=graph_path, output_path=Path(folder)/"main.py")
            definitions, runner = path.read_text().split("# ---- Source: qollab_smoke/cut_pilot.py ----", 1)
            ns = {"backend": BasicSimulator(), "__name__": "__main__"}
            with patch.dict(sys.modules, {"qiskit": qiskit, **scipy_stub()}), redirect_stdout(StringIO()):
                exec(compile(definitions, str(path), "exec"), ns)
                actual_builder = ns["build_qaoa"]
                def test_builder(w, c, alpha, theta):
                    qc = actual_builder(w, c, alpha, theta)
                    qc.test_bits = "".join(map(str, reversed(exact_fixed_size_cut(w, c))))
                    return qc
                ns["build_qaoa"] = test_builder
                exec(compile(runner, str(path), "exec"), ns)
            report = ns["report"]
            self.assertEqual(report["weights"], graph["weights"])
            self.assertEqual(report["labels"], graph["labels"])
            self.assertEqual(report["mode"], "protein_tree")
            self.assertEqual(report["status"], "completed")
            self.assertEqual(report["topology_metrics"]["unrooted_rf_to_exact"], 0)
            self.assertEqual(report["topology_metrics"]["split_nodes"], 4)
            self.assertEqual(len(report["attempts"]), 9)
            self.assertEqual(ns["backend"].calls, 279)
            saved = Path(folder)/"test_double_report.json"
            saved.write_text(json.dumps(report))
            manifest = export_tables([saved], Path(folder)/"tables")
            self.assertTrue(manifest["scope"].startswith("Protein-derived"))
            self.assertEqual(manifest["dataset_provenance"]["preprocessing"]["source_sha256"],
                             graph["provenance"]["source_sha256"])
            self.assertEqual(read_report(saved)["dataset_provenance"]["kind"], "protein_derived")

    def test_protein_bundle_rejects_unreviewed_or_tampered_export(self):
        import shutil
        root = Path(__file__).resolve().parents[1]
        source = root/"results/protein1822"
        with TemporaryDirectory() as folder:
            folder = Path(folder)
            for name in ("graph.json", "manifest.json", "affinity_matrix.csv"):
                shutil.copyfile(source/name, folder/name)
            path = folder/"graph.json"
            original = json.loads(path.read_text())
            changed = deepcopy(original)
            changed["weights"][0][1] += 1
            changed["weights"][1][0] += 1
            path.write_text(json.dumps(changed))
            with self.assertRaisesRegex(ValueError, "audited affinity"):
                protein_graph(path)
            changed = deepcopy(original)
            changed["provenance"]["source_sha256"] = "mismatched"
            path.write_text(json.dumps(changed))
            with self.assertRaisesRegex(ValueError, "source_sha256 disagree"):
                protein_graph(path)
            path.write_text(json.dumps(original))
            manifest_path = folder/"manifest.json"
            manifest = json.loads(manifest_path.read_text())
            manifest["status"] = "requires_review"
            manifest_path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "requires review"):
                protein_graph(path)

    def test_failed_run_is_retained_without_invented_topology(self):
        w, labels = toy_fixture("toy4")
        report = {"weights": w, "labels": labels, "mode": "toy_tree", "status": "failed",
                  "execution_mode": "multi_job", "config": asdict(QAOAConfig()),
                  "jobs": [], "evaluations": [], "attempts": [], "error": "backend error"}
        with TemporaryDirectory() as folder:
            path = Path(folder)/"failed.json"
            path.write_text(json.dumps(report))
            reports, tables, pairs = comparison_tables([path])
            self.assertEqual(tables["runs"][0]["status"], "failed")
            self.assertIsNone(tables["runs"][0]["rf_to_exact"])
            self.assertIsNone(reports[0]["topology_metrics"])


if __name__ == "__main__":
    unittest.main()

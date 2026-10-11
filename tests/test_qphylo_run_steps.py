"""Verify checkpoint execution against uninterrupted optimization."""

from dataclasses import replace
from html import escape
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from qphylo_qaoa import QAOAConfig, QAOACutSolver, SampleBatch
from qphylo_run_steps import canonical, checkpoint_from_report, run_one_job_step
from prepare_qollab import checkpoint_from_file

W = ((0, 9, 1, 1), (9, 0, 1, 1), (1, 1, 0, 9), (1, 1, 9, 0))
LABELS = list("ABCD")


class ResumeTests(unittest.TestCase):
    def step(self, sampler, config, state=None, weights=W, labels=LABELS, c=2,
             mode="optimized_cut", backend_label="test-backend"):
        return run_one_job_step(weights, labels, c, mode, config,
                                sampler, backend_label, state)

    def test_each_run_submits_once_and_replays_same_optimizer(self):
        def data(w, c, alpha, theta, context):
            # Depend on the candidate to exercise real acceptance decisions.
            optimum = min(1000, max(0, int(1000 * theta[0] / 7)))
            return SampleBatch({"0011": optimum, "0000": 1000 - optimum},
                               {"status": "DONE", "job_id": str(context["evaluation"])})
        uninterrupted = QAOACutSolver(data, QAOAConfig(max_evaluations=9, max_retries=0))
        expected_bits = uninterrupted(W, 2)
        sampler = Mock(side_effect=data)
        state = None
        for _ in range(9):
            before = sampler.call_count
            report = self.step(sampler, uninterrupted.config, state)
            self.assertLessEqual(sampler.call_count - before, 1)
            self.assertNotEqual(report["status"], "failed", report.get("error"))
            state = json.loads(canonical(report["resume_state"]))
            if report["status"] == "completed":
                break
        else:
            self.fail("Optimizer did not complete within its evaluation budget.")
        self.assertEqual(report["bits"], expected_bits)
        self.assertEqual(canonical(report["evaluations"]), canonical(uninterrupted.evaluations))
        self.assertEqual(sampler.call_count, len(uninterrupted.evaluations))
        self.assertEqual(report["evaluations"][-1]["phase"], "final")
        before = sampler.call_count
        repeated = self.step(sampler, uninterrupted.config, state)
        self.assertEqual(repeated["status"], "completed")
        self.assertEqual(sampler.call_count, before)

    def test_reinforcement_crosses_run_boundaries(self):
        sampler = Mock(return_value=SampleBatch({"00": 400, "01": 600}))
        cfg = QAOAConfig(max_evaluations=2, initial_penalty=0.1, max_retries=1)
        state = None
        for run in range(4):
            before = sampler.call_count
            report = self.step(sampler, cfg, state, weights=((0, 2), (2, 0)),
                               labels=["A", "B"], c=1)
            self.assertEqual(sampler.call_count - before, 1)
            state = json.loads(canonical(report["resume_state"]))
        self.assertEqual(report["status"], "completed")
        self.assertEqual([row["action"] for row in report["attempts"]], ["reinforce", "accept"])
        self.assertEqual([row["penalty"] for row in report["attempts"]], [0.1, 500.1])

    def test_config_change_is_rejected_before_any_submission(self):
        sampler = Mock(return_value=SampleBatch({"0011": 1000}))
        cfg = QAOAConfig(max_evaluations=2)
        state = self.step(sampler, cfg)["resume_state"]
        before = sampler.call_count
        with self.assertRaises(ValueError):
            self.step(sampler, replace(cfg, shots=2000), state)
        with self.assertRaises(ValueError):
            self.step(sampler, cfg, state, backend_label="another-backend")
        self.assertEqual(sampler.call_count, before)

    def test_failed_submission_does_not_retry_or_silently_resume(self):
        sampler = Mock(side_effect=RuntimeError("backend failed"))
        cfg = QAOAConfig(max_evaluations=9)
        report = self.step(sampler, cfg)
        self.assertEqual(report["status"], "failed")
        self.assertEqual(sampler.call_count, 1)
        self.assertEqual(report["completed_measurements"], 0)
        with self.assertRaises(ValueError):
            self.step(sampler, cfg, report["resume_state"])
        self.assertEqual(sampler.call_count, 1)

    def test_recursive_calls_resume_without_repeat_jobs(self):
        def data(w, c, alpha, theta, context):
            key = "0001" if c == 1 else "0011"
            return SampleBatch({key: 1000})
        sampler = Mock(side_effect=data)
        cfg = QAOAConfig(max_evaluations=2, max_retries=0)
        state = None
        for _ in range(6):
            before = sampler.call_count
            report = self.step(sampler, cfg, state, mode="toy_tree")
            self.assertLessEqual(sampler.call_count - before, 1)
            self.assertNotEqual(report["status"], "failed", report.get("error"))
            state = report["resume_state"]
            if report["status"] == "completed":
                break
        self.assertEqual(report["status"], "completed")
        self.assertEqual(report["tree"], "(('A','B'),('C','D'));")
        self.assertEqual(sampler.call_count, 4)

    def test_saved_html_checkpoint_is_read_without_execution(self):
        sampler = Mock(return_value=SampleBatch({"0011": 1000}))
        report = self.step(sampler, QAOAConfig(max_evaluations=2))
        with TemporaryDirectory() as folder:
            path = Path(folder) / "result.html"
            path.write_text("<div>RESUME_STATE_JSON: " + escape(
                json.dumps(report["resume_state"])) + "</div>")
            self.assertEqual(checkpoint_from_file(path), report["resume_state"])

    def test_recover_first_measurement_from_known_old_limit_error(self):
        sampler = Mock(return_value=SampleBatch({"0011": 1000}, {"status": "DONE"}))
        report = self.step(sampler, QAOAConfig(max_evaluations=2))
        legacy = {key: report[key] for key in ("weights", "labels", "mode",
                                             "config", "backend_label", "evaluations")}
        legacy["error"] = "RuntimeError: The code has already filed one circuit execution job; only one job is allowed per code run"
        recovered = checkpoint_from_report(legacy)
        self.assertEqual(recovered["history"], report["resume_state"]["history"])
        result = self.step(sampler, QAOAConfig(max_evaluations=2), recovered)
        self.assertEqual(result["status"], "completed")
        self.assertEqual(sampler.call_count, 2)


if __name__ == "__main__":
    unittest.main()

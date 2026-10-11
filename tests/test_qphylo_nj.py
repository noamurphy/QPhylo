"""Input checks, failure double and optional real R protein-baseline regression."""

import csv
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from qphylo_nj import read_alignment, run, find_rscript, DEFAULT_FASTA, DEFAULT_REPORT, DEFAULT_REFERENCE


class NJPreparationTests(unittest.TestCase):
    def test_alignment_preserves_gap_columns_and_first_token_labels(self):
        with TemporaryDirectory() as folder:
            path = Path(folder)/"aligned.fasta"
            path.write_text("\n>A description\na-c.\n\n>B\nac--\n>C\na--c\n")
            labels, seqs = read_alignment(path)
            self.assertEqual(labels, ["A", "B", "C"])
            self.assertEqual(seqs, ["A-C-", "AC--", "A--C"])

    def test_rejects_unaligned_ambiguous_duplicate_or_too_few_records(self):
        examples = [("\n>A\nACD\n>B\nACDE\n>C\nACE\n", "equal column"),
                    (">A\nACX\n>B\nACD\n>C\nACE\n", "canonical"),
                    (">A\nACD\n>A\nACE\n>C\nACF\n", "Duplicate"),
                    (">A\nACD\n>B\nACE\n", "three taxa")]
        with TemporaryDirectory() as folder:
            path = Path(folder)/"input.fasta"
            for text, error in examples:
                path.write_text(text)
                with self.assertRaisesRegex(ValueError, error):
                    read_alignment(path)

    def test_mismatched_source_rejected_before_R_or_outputs(self):
        with TemporaryDirectory() as folder:
            path = Path(folder)/"changed.fasta"
            path.write_text(DEFAULT_FASTA.read_text() + "\n")
            output = Path(folder)/"out"
            with patch("qphylo_nj.subprocess.run") as execute:
                with self.assertRaisesRegex(ValueError, "source hash"):
                    run(fasta=path, output=output)
            execute.assert_not_called()
            self.assertFalse(output.exists())

    def test_R_failure_keeps_manifest_and_clears_stale_derived_tree(self):
        with TemporaryDirectory() as folder:
            output = Path(folder)/"out"
            output.mkdir()
            (output/"nj_tree.newick").write_text("stale")
            failed = SimpleNamespace(returncode=1, stdout="", stderr="Missing R package: phangorn")
            with patch("qphylo_nj.find_rscript", return_value="Rscript"), \
                    patch("qphylo_nj.subprocess.run", return_value=failed) as execute:
                with self.assertRaisesRegex(RuntimeError, "Missing R package"):
                    run(output=output)
            self.assertEqual(execute.call_count, 1)
            self.assertFalse((output/"nj_tree.newick").exists())
            manifest = json.loads((output/"manifest.json").read_text())
            self.assertEqual(manifest["status"], "failed")
            self.assertEqual(manifest["model"]["gamma_categories"], 5)
            self.assertEqual(manifest["model"]["gamma_shape"], 0.5)
            self.assertEqual(manifest["alignment"]["columns"], 389)
            self.assertEqual(manifest["alignment"]["gap_characters"], 0)
            self.assertEqual(manifest["bootstrap"]["replicates"], 100)


class NJIntegrationTests(unittest.TestCase):
    def test_real_R_protein_topology_and_bootstrap_percentage_regression(self):
        try:
            executable = find_rscript()
        except FileNotFoundError:
            self.skipTest("Rscript unavailable")
        probe = subprocess.run([executable, "-e",
            'quit(status=if (requireNamespace("ape", quietly=TRUE) && requireNamespace("phangorn", quietly=TRUE)) 0 else 1)'],
            capture_output=True, text=True)
        if probe.returncode:
            self.skipTest("ape/phangorn unavailable in selected R")
        with TemporaryDirectory() as folder:
            output = Path(folder)
            manifest = run(output=output, rscript=executable)
            self.assertEqual(manifest["status"], "completed")
            pairs = {(row["method1"], row["method2"]): int(row["rf"])
                     for row in manifest["pairwise_comparisons"]}
            self.assertEqual(pairs[("QAOA", "Exact_NMcut")], 0)
            self.assertEqual(pairs[("QAOA", "Model_NJ")], 2)
            self.assertEqual(pairs[("Model_NJ", "Provided_reference")], 0)
            bootstrap_lines = (output/"nj_bootstrap_trees.newick").read_text().splitlines()
            self.assertEqual(len(bootstrap_lines), 100)
            with (output/"splits.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            # Known protein1822 split support, in percent rather than proportions.
            pair = next(row for row in rows if row["method"] == "Model_NJ" and row["side1"] == "taxon1,taxon8")
            self.assertEqual(float(pair["bootstrap_percent"]), 94)


if __name__ == "__main__":
    unittest.main()

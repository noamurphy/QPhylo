"""Protein input, BLAST scoring conventions and core integration checks."""

import csv
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from qphylo_core import validate_weights
from qphylo_proteins import (Protein, blast_command, graph_diagnostics,
                            normalized_affinities, parse_blast_scores,
                            prepare_protein_graph, read_protein_fasta)


def hit(i, j, score, evalue="0", end=4):
    return f"qphylo_{i:06d}\tqphylo_{j:06d}\t{score}\t{evalue}\t1\t{end}\t1\t{end}\t{end}\n"


class ProteinTests(unittest.TestCase):
    def test_fasta_labels_gaps_case_and_real_input(self):
        with TemporaryDirectory() as folder:
            fasta = Path(folder)/"input.fa"
            fasta.write_text(">A first protein\nac-d\n.e f\n\n>B\nWQRT\n")
            proteins = read_protein_fasta(fasta)
            self.assertEqual(proteins[0], Protein("A", "A first protein", "ACDEF", 7, 2))
            self.assertEqual(proteins[1].sequence, "WQRT")
        real = read_protein_fasta(Path(__file__).resolve().parents[1]/
                                  "PhyloTree_ReConst/alignments/8/assembled_sequences_1822.fasta")
        self.assertEqual(len(real), 8)
        self.assertEqual([r.label for r in real], [f"taxon{i}" for i in range(1,9)])
        self.assertEqual({len(r.sequence) for r in real}, {389})

    def test_invalid_or_ambiguous_identifiers_fail(self):
        cases = ["ACDE\n", ">\nACDE\n", ">A\n---\n", ">A\nACDE*\n", ">A\nAC3E\n",
                 ">A one\nACDE\n>A two\nWQRT\n", ""]
        with TemporaryDirectory() as folder:
            path = Path(folder)/"input.fa"
            for content in cases:
                path.write_text(content)
                with self.subTest(content=content), self.assertRaises(ValueError):
                    read_protein_fasta(path)

    def test_hsp_choice_evalue_missing_direction_and_normalization(self):
        records = (Protein("A", "A", "ACDE", 4, 0), Protein("B", "B", "WQRT", 4, 0))
        with TemporaryDirectory() as folder:
            path = Path(folder)/"hits.tsv"
            path.write_text(hit(1,1,100)+hit(2,2,200)+hit(1,2,30)+hit(1,2,60,"10")+
                            hit(2,1,90,"10.1"))
            scores, info = parse_blast_scores(path, records)
            self.assertEqual(scores, ((100.0,60.0),(0.0,200.0)))
            self.assertEqual(info["excluded_above_evalue"], 1)
            directed, similarity, weights = normalized_affinities(scores)
            self.assertEqual(directed, ((100.0,40.0),(0.0,100.0)))
            self.assertEqual(similarity, ((100.0,20.0),(20.0,100.0)))
            self.assertEqual(weights, ((0.0,20.0),(20.0,0.0)))
            self.assertEqual(validate_weights(weights), weights)

    def test_bad_hits_and_missing_self_are_rejected(self):
        records = (Protein("A", "A", "ACDE", 4, 0), Protein("B", "B", "WQRT", 4, 0))
        malformed = [hit(1,1,100), hit(1,1,100)+hit(3,2,30),
                     hit(1,1,100)+hit(2,2,"nan"), hit(1,1,100)+hit(2,2,-1),
                     hit(1,1,100)+hit(2,2,100,end=10), "wrong\tcolumns\n"]
        with TemporaryDirectory() as folder:
            path = Path(folder)/"hits.tsv"
            for content in malformed:
                path.write_text(content)
                with self.subTest(content=content), self.assertRaises(ValueError):
                    parse_blast_scores(path, records)

    def test_symmetry_keeps_values_above_100_and_detects_disconnection(self):
        _, similarity, weights = normalized_affinities(((10,20,0),(16,10,0),(0,0,10)))
        self.assertEqual(weights[0][1], 180)
        diag = graph_diagnostics(weights, ("A","B","C"))
        self.assertEqual(diag["connected_components"], [["A","B"],["C"]])
        self.assertEqual(diag["zero_degree_taxa"], ["C"])
        self.assertFalse(diag["connected"])

    def test_live_runner_contract_exports_graph_and_exact_core_tree(self):
        # Explicit subprocess double checks orchestration; not a BLAST measurement.
        with TemporaryDirectory() as folder:
            fasta, out = Path(folder)/"input.fa", Path(folder)/"out"
            fasta.write_text(">A\nACDE\n>B\nWQRT\n>C\nPLKV\n>D\nMHSN\n")
            def run(command, **kwargs):
                if "-version" in command:
                    return SimpleNamespace(returncode=0,stdout="blastp: test-double\n",stderr="")
                target = Path(command[command.index("-out")+1])
                target.write_text("".join(hit(i,j,100 if i==j else 90 if (i-1)//2==(j-1)//2 else 10)
                                          for i in range(1,5) for j in range(1,5)))
                self.assertEqual(Path(command[command.index("-subject")+1]).resolve(), (out/"ungapped.fasta").resolve())
                self.assertEqual(command[command.index("-max_target_seqs")+1], "4")
                return SimpleNamespace(returncode=0,stdout="",stderr="")
            with patch("qphylo_proteins.shutil.which",return_value="/test/blastp"), \
                 patch("qphylo_proteins.subprocess.run",side_effect=run):
                manifest=prepare_protein_graph(fasta,out)
            self.assertEqual(manifest["status"],"graph_ready")
            self.assertEqual(manifest["exact_reference"]["tree"],"(('A','B'),('C','D'));")
            graph=json.loads((out/"graph.json").read_text())
            self.assertEqual(graph["labels"],["A","B","C","D"])
            self.assertEqual(graph["weights"][0],[0,90,10,10])
            self.assertEqual(len(manifest["blast_hits_sha256"]),64)
            self.assertTrue((out/"directed_bit_scores.csv").exists())

    def test_failed_execution_has_manifest_and_no_successful_graph(self):
        with TemporaryDirectory() as folder:
            fasta, out = Path(folder)/"input.fa", Path(folder)/"out"
            fasta.write_text(">A\nACDE\n>B\nWQRT\n")
            out.mkdir();(out/"graph.json").write_text("stale")
            def run(command,**kwargs):
                return SimpleNamespace(returncode=0 if "-version" in command else 1,
                                       stdout="test-double",stderr="bad BLAST input")
            with patch("qphylo_proteins.shutil.which",return_value="/test/blastp"), \
                 patch("qphylo_proteins.subprocess.run",side_effect=run), \
                 self.assertRaisesRegex(RuntimeError,"BLASTP exited"):
                prepare_protein_graph(fasta,out)
            self.assertEqual(json.loads((out/"manifest.json").read_text())["status"],"failed")
            self.assertFalse((out/"graph.json").exists())

    def test_disconnected_import_exported_for_review_without_fake_tree(self):
        with TemporaryDirectory() as folder:
            fasta, hits, out=Path(folder)/"input.fa",Path(folder)/"hits.tsv",Path(folder)/"out"
            fasta.write_text(">A\nACDE\n>B\nWQRT\n")
            hits.write_text(hit(1,1,100)+hit(2,2,100))
            manifest=prepare_protein_graph(fasta,out,imported_hits=hits,imported_version="test-double")
            self.assertEqual(manifest["status"],"graph_requires_review")
            self.assertTrue((out/"graph.json").exists())
            self.assertFalse((out/"exact_tree.newick").exists())
            with patch("qphylo_proteins.subprocess.run") as runner:
                prepare_protein_graph(fasta,out,validate_only=True)
            runner.assert_not_called()
            self.assertFalse((out/"graph.json").exists())


if __name__ == "__main__":
    unittest.main()

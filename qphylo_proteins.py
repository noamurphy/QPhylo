"""Protein FASTA -> local BLASTP bit scores -> symmetric graph-cut affinities.

Standard Python only. BLAST+ is an external program, not a Python dependency.
Never substitutes global alignments, rescales to another range, or repairs edges.
The exact tree is an optional small-graph diagnostic, not a quantum result.

From the project folder after activating qphylo:
    python qphylo_proteins.py                 # real eight-protein BLAST graph
    python qphylo_proteins.py --validate-only # FASTA check without BLAST
Defaults: existing FO protein1822 input; outputs in results/protein1822/.
Validation-only uses results/protein1822_input_check/ by default.
--fasta and --output select a different input and a separate run directory.
"""

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
from math import isfinite
from pathlib import Path
import shutil
import subprocess
import time

from qphylo_core import reconstruct, validate_weights


PROTEIN_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO")
BLAST_FIELDS = "qseqid sseqid bitscore evalue qstart qend sstart send length"


@dataclass(frozen=True)
class Protein:
    label: str
    description: str
    sequence: str
    original_length: int
    removed_gaps: int


def read_protein_fasta(path):
    """Keep first-token labels/order; uppercase sequences and strip '-' / '.'."""
    records, seen = [], set()
    description, chunks = None, []

    def finish():
        if description is None:
            return
        label = description.split()[0]
        raw = "".join(chunks).upper()
        sequence = raw.replace("-", "").replace(".", "")
        if not sequence:
            raise ValueError(f"FASTA record {label!r} has no amino acids.")
        invalid = sorted(set(sequence) - PROTEIN_ALPHABET)
        if invalid:
            raise ValueError(f"FASTA record {label!r} contains unsupported residues {invalid}; "
                             "provide protein sequences without stop symbols or digits.")
        records.append(Protein(label, description, sequence, len(raw), len(raw)-len(sequence)))

    for line_number, line in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            description = line[1:].strip()
            if not description:
                raise ValueError(f"Empty FASTA header at line {line_number}.")
            label = description.split()[0]
            if label in seen:
                raise ValueError(f"Duplicate FASTA identifier {label!r}; first header tokens must be unique.")
            seen.add(label)
            chunks = []
        else:
            if description is None:
                raise ValueError(f"Sequence before first FASTA header at line {line_number}.")
            chunks.append("".join(line.split()))
    finish()
    if not records:
        raise ValueError("FASTA contains no sequences.")
    return tuple(records)


def blast_ids(records):
    # Safe internal IDs avoid BLAST rewriting descriptive/user-supplied labels.
    return tuple(f"qphylo_{i+1:06d}" for i in range(len(records)))


def write_blast_fasta(records, path):
    with Path(path).open("w") as stream:
        for identifier, record in zip(blast_ids(records), records):
            stream.write(f">{identifier}\n")
            for start in range(0, len(record.sequence), 80):
                stream.write(record.sequence[start:start+80] + "\n")


def blast_command(executable, fasta, output, n, evalue=10.0):
    if not isfinite(evalue) or evalue <= 0:
        raise ValueError("E-value threshold must be finite and positive.")
    return [str(executable), "-task", "blastp", "-query", str(fasta),
            "-subject", str(fasta), "-out", str(output), "-outfmt", "6 " + BLAST_FIELDS,
            "-evalue", str(evalue), "-matrix", "BLOSUM62", "-gapopen", "11",
            "-gapextend", "1", "-seg", "no", "-comp_based_stats", "2",
            "-max_hsps", "1", "-max_target_seqs", str(n), "-num_threads", "1"]


def parse_blast_scores(path, records, evalue=10.0):
    """Maximum qualifying single-HSP bit score per directed pair; absent=0.

    Self hits are required and positive. Score sums across HSPs are not used.
    Strictly validate internal IDs, numeric values and alignment coordinates.
    """
    if not isfinite(evalue) or evalue <= 0:
        raise ValueError("E-value threshold must be finite and positive.")
    identifiers = {identifier: i for i, identifier in enumerate(blast_ids(records))}
    n = len(records)
    scores = [[0.0] * n for _ in range(n)]
    rows, excluded = 0, 0
    for line_number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            raise ValueError(f"BLAST line {line_number}: expected nine tab-separated columns ({BLAST_FIELDS}).")
        query, subject = fields[:2]
        if query not in identifiers or subject not in identifiers:
            raise ValueError(f"BLAST line {line_number}: unknown query/subject ID; use the exported ungapped.fasta.")
        i, j = identifiers[query], identifiers[subject]
        try:
            score, expect = map(float, fields[2:4])
            qs, qe, ss, se, length = map(int, fields[4:])
        except ValueError as exc:
            raise ValueError(f"BLAST line {line_number}: invalid numeric value.") from exc
        if not all(map(isfinite, (score, expect))) or score <= 0 or expect < 0:
            raise ValueError(f"BLAST line {line_number}: require positive finite bit score and nonnegative finite E-value.")
        if not (1 <= qs <= qe <= len(records[i].sequence) and
                1 <= ss <= se <= len(records[j].sequence) and length > 0):
            raise ValueError(f"BLAST line {line_number}: invalid protein alignment coordinates.")
        rows += 1
        if expect > evalue:
            excluded += 1
            continue
        scores[i][j] = max(scores[i][j], score)
    missing = [records[i].label for i in range(n) if scores[i][i] <= 0]
    if missing:
        raise ValueError(f"Missing positive qualifying self scores for {missing}; cannot normalize. "
                         "Check BLAST output/settings and protein lengths.")
    return tuple(map(tuple, scores)), {"parsed_hsps": rows, "excluded_above_evalue": excluded}


def normalized_affinities(scores):
    """Return directed normalized scores, symmetric similarity, zero-diag graph.

    b_ij is the best directed BLAST bit score, s_i=b_ii.
    d_ij=100*b_ij/((s_i+s_j)/2); a_ij=(d_ij+d_ji)/2.
    A missing direction contributes zero. Keep a_ii=100 for audit, w_ii=0
    for the existing core's no-self-edge convention. Never clip values to 100.
    """
    n = len(scores)
    if not n or any(len(row) != n for row in scores):
        raise ValueError("Directed scores must be a nonempty square matrix.")
    if any(not isfinite(v) or v < 0 for row in scores for v in row):
        raise ValueError("Directed scores must be finite and nonnegative.")
    if any(scores[i][i] <= 0 for i in range(n)):
        raise ValueError("Self scores must be positive.")
    directed = tuple(tuple(100.0*scores[i][j]/((scores[i][i]+scores[j][j])/2)
                           for j in range(n)) for i in range(n))
    symmetric = tuple(tuple((directed[i][j]+directed[j][i])/2 for j in range(n))
                      for i in range(n))
    weights = validate_weights(tuple(tuple(0.0 if i == j else symmetric[i][j]
                                          for j in range(n)) for i in range(n)))
    return directed, symmetric, weights


def graph_diagnostics(weights, labels):
    remaining, components = set(range(len(labels))), []
    while remaining:
        frontier, component = [min(remaining)], set()
        while frontier:
            i = frontier.pop()
            if i in component:
                continue
            component.add(i)
            frontier.extend(j for j, value in enumerate(weights[i]) if value > 0 and j not in component)
        remaining -= component
        components.append([labels[i] for i in sorted(component)])
    return {"connected_components": components,
            "connected": len(components) == 1,
            "zero_degree_taxa": [labels[i] for i, row in enumerate(weights) if sum(row) == 0],
            "positive_edges": sum(weights[i][j] > 0 for i in range(len(labels)) for j in range(i)),
            "possible_edges": len(labels)*(len(labels)-1)//2,
            "max_off_diagonal_affinity": max((v for i, row in enumerate(weights)
                                                for j, v in enumerate(row) if i != j), default=0)}


def write_matrix(path, labels, matrix):
    with Path(path).open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["taxon", *labels])
        writer.writerows([label, *row] for label, row in zip(labels, matrix))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare_protein_graph(fasta_path, output, blastp="blastp", evalue=10.0,
                          imported_hits=None, imported_version=None, validate_only=False):
    """Export audited graph and small exact tree. No quantum jobs or installs."""
    source, folder = Path(fasta_path).resolve(), Path(output).resolve()
    records = read_protein_fasta(source)
    if len(records) < 2:
        raise ValueError("Graph reconstruction needs at least two protein records.")
    if not isfinite(evalue) or evalue <= 0:
        raise ValueError("E-value threshold must be finite and positive.")
    if source == folder / "ungapped.fasta":
        raise ValueError("Choose a separate output folder; the source FASTA must not be overwritten.")
    # Resolve dependencies before producing a run folder, unless explicitly
    # exporting a validated input only. This avoids a false 'completed' graph.
    executable, version = None, imported_version
    if not validate_only and imported_hits is None:
        executable = shutil.which(str(blastp))
        if not executable:
            raise FileNotFoundError("BLASTP is not available. Activate qphylo and install the blast package, "
                                    "or supply --blastp /path/to/blastp. Use --validate-only to check FASTA first.")
        version = subprocess.run([executable, "-version"], check=True, capture_output=True,
                                 text=True, timeout=30).stdout.strip()
    if imported_hits is not None and not imported_version:
        raise ValueError("Importing hits requires --blast-version; use the matching exported ungapped.fasta.")
    folder.mkdir(parents=True, exist_ok=True)
    # A new attempt owns its derived products. Never leave a stale completed
    # graph after a failed/validation-only rerun in this same output directory.
    for filename in ("graph.json", "affinity_matrix.csv", "directed_bit_scores.csv",
                     "directed_normalized_scores.csv", "normalized_similarity.csv",
                     "exact_tree.newick", "exact_trace.json", "blast.stderr.txt"):
        (folder / filename).unlink(missing_ok=True)
    ungapped, hits = folder / "ungapped.fasta", folder / "blast_hits.tsv"
    write_blast_fasta(records, ungapped)
    with (folder / "label_map.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["blast_id", "taxon", "description", "original_length", "ungapped_length", "removed_gaps"])
        for identifier, record in zip(blast_ids(records), records):
            writer.writerow([identifier, record.label, record.description, record.original_length,
                             len(record.sequence), record.removed_gaps])
    groups = {}
    for record in records:
        groups.setdefault(record.sequence, []).append(record.label)
    manifest = {"schema": "qphylo-protein-graph-v1", "status": "input_validated",
        "source_fasta": str(source), "source_sha256": sha256(source),
        "ungapped_sha256": sha256(ungapped), "n_taxa": len(records),
        "labels": [record.label for record in records],
        "duplicate_sequence_groups": [group for group in groups.values() if len(group) > 1],
        "blast_version": version,
        "execution": "imported_hits" if imported_hits is not None else "local_blastp",
        "imported_settings_verified": False if imported_hits is not None else None,
        "settings": {"task": "blastp", "matrix": "BLOSUM62", "evalue": evalue,
            "gapopen": 11, "gapextend": 1, "seg": "no", "comp_based_stats": 2,
            "max_hsps": 1, "max_target_seqs": len(records), "num_threads": 1},
        "normalization": "100 * best directed HSP bit score / mean of the two positive self bit scores",
        "symmetry": "arithmetic mean of both directed normalized scores; missing direction contributes zero",
        "graph_diagonal": 0, "normalized_similarity_diagonal": 100,
        "clipping": "none", "reference_tree": None,
        "scope": "Protein-derived affinity integration; not biological accuracy or a complete Onodera reproduction."}
    manifest_path = folder / "manifest.json"
    def save():
        manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    save()
    if validate_only:
        return manifest
    started = time.monotonic()
    try:
        if imported_hits is not None:
            imported = Path(imported_hits).resolve()
            if imported != hits:
                shutil.copyfile(imported, hits)
            manifest["imported_hits"] = str(imported)
        else:
            command = blast_command(executable, ungapped, hits, len(records), evalue)
            manifest["command"] = command
            manifest["status"] = "running_blast"
            save()
            result = subprocess.run(command, capture_output=True, text=True, timeout=600)
            (folder / "blast.stderr.txt").write_text(result.stderr)
            if result.returncode:
                raise RuntimeError(f"BLASTP exited {result.returncode}: {result.stderr.strip()}")
        scores, stats = parse_blast_scores(hits, records, evalue)
        directed, symmetric, weights = normalized_affinities(scores)
        labels = tuple(record.label for record in records)
        manifest.update(stats)
        manifest["blast_hits_sha256"] = sha256(hits)
        manifest["directed_asymmetric_pairs"] = [
            {"taxa": [labels[i], labels[j]], "bits_i_to_j": scores[i][j], "bits_j_to_i": scores[j][i]}
            for i in range(len(labels)) for j in range(i) if scores[i][j] != scores[j][i]]
        manifest["missing_directed_off_diagonal_hits"] = sum(
            scores[i][j] == 0 for i in range(len(labels)) for j in range(len(labels)) if i != j)
        manifest["diagnostics"] = graph_diagnostics(weights, labels)
        write_matrix(folder / "directed_bit_scores.csv", labels, scores)
        write_matrix(folder / "directed_normalized_scores.csv", labels, directed)
        write_matrix(folder / "normalized_similarity.csv", labels, symmetric)
        write_matrix(folder / "affinity_matrix.csv", labels, weights)
        graph = {"schema": manifest["schema"], "labels": labels, "weights": weights,
                 "provenance": {"manifest": str(manifest_path), "source_sha256": manifest["source_sha256"],
                                "blast_hits_sha256": manifest["blast_hits_sha256"]}}
        (folder / "graph.json").write_text(json.dumps(graph, indent=2, allow_nan=False) + "\n")
        if not manifest["diagnostics"]["connected"]:
            manifest["status"] = "graph_requires_review"
            manifest["exact_reference"] = {"status": "not_run", "reason": "Disconnected graph; no artificial edges or automatic component trees."}
        else:
            manifest["status"] = "graph_ready"
            if len(labels) <= 12:
                tree, trace = reconstruct(weights, labels)
                (folder / "exact_tree.newick").write_text(tree + "\n")
                (folder / "exact_trace.json").write_text(json.dumps(trace, indent=2) + "\n")
                manifest["exact_reference"] = {"status": "completed", "tree": tree,
                                                "meaning": "Exact fixed-size-cut recursion on this affinity graph, not evolutionary ground truth."}
            else:
                manifest["exact_reference"] = {"status": "not_run", "reason": "Exponential diagnostic limited to at most 12 taxa."}
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        manifest["elapsed_seconds"] = time.monotonic() - started
        save()
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", default="PhyloTree_ReConst/alignments/8/assembled_sequences_1822.fasta")
    parser.add_argument("--output", help="Default: results/protein1822, or protein1822_input_check for validation only")
    parser.add_argument("--blastp", default="blastp", help="Executable name/path in the active environment")
    parser.add_argument("--evalue", type=float, default=10.0)
    parser.add_argument("--validate-only", action="store_true", help="Validate FASTA and export ungapped input; no BLAST run")
    parser.add_argument("--hits", help="Import nine-column BLAST output made from this input's internal IDs")
    parser.add_argument("--blast-version", help="Required provenance when importing --hits")
    args = parser.parse_args()
    if args.output is None:
        args.output = "results/protein1822_input_check" if args.validate_only else "results/protein1822"
    try:
        result = prepare_protein_graph(args.fasta, args.output, args.blastp, args.evalue,
                                       args.hits, args.blast_version, args.validate_only)
    except (ValueError, FileNotFoundError, RuntimeError, subprocess.SubprocessError) as exc:
        parser.exit(2, f"Protein preprocessing stopped: {exc}\n")
    print("Status:", result["status"], "taxa:", result["n_taxa"])
    print("Output folder:", Path(args.output).resolve())
    if "diagnostics" in result:
        print("Connected components:", result["diagnostics"]["connected_components"])
        print("Positive edges:", result["diagnostics"]["positive_edges"], "/", result["diagnostics"]["possible_edges"])
    if result.get("exact_reference", {}).get("status") == "completed":
        print("Exact affinity-tree reference:", result["exact_reference"]["tree"])
    if result["status"] == "graph_requires_review":
        parser.exit(2, "Graph saved for inspection; resolve disconnected components before reconstruction.\n")

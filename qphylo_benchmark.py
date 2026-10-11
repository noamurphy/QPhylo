"""Synthetic fixtures and exact diagnostics; never used to pick QAOA answers."""

from itertools import combinations
from math import isclose
import hashlib
import json

from qphylo_core import cut_weight, exact_fixed_size_cut, normalized_cut, reconstruct
from qphylo_qaoa import decode_bitstring


def toy_fixture(name):
    if name == "toy4":
        return ((0, 9, 1, 1), (9, 0, 1, 1), (1, 1, 0, 9), (1, 1, 9, 0)), tuple("ABCD")
    if name != "hierarchical8":
        raise ValueError("fixture must be toy4 or hierarchical8.")
    # Pair affinity 9, different pairs in the same half 3, other half 1.
    weights = tuple(tuple(0 if i == j else 9 if i // 2 == j // 2 else
                          3 if i // 4 == j // 4 else 1
                          for j in range(8)) for i in range(8))
    return weights, tuple("ABCDEFGH")


def graph_signature(weights, labels):
    payload = {"weights": [[float(v) for v in row] for row in weights],
               "labels": list(labels)}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def job_upper_bound(n, config):
    # The longest possible recursion peels off one taxon at every split.
    calls = sum(k // 2 for k in range(3, n + 1))
    return calls * config.max_evaluations * (config.max_retries + 1)


def topology_splits(trace, labels):
    """Canonical unrooted nontrivial bipartitions, independent of child order."""
    all_taxa = frozenset(labels)
    splits = set()
    for row in trace:
        for side in (row["left"], row["right"]):
            group = frozenset(side)
            other = all_taxa - group
            if 2 <= len(group) and 2 <= len(other):
                a, b = tuple(sorted(group)), tuple(sorted(other))
                splits.add(min((a, b), key=lambda x: (len(x), x)))
    return splits


def candidate_taxa(trace, labels):
    """Match solver call indices to taxa using the completed recursive trace."""
    nodes, call = iter(trace), 0
    result = {}
    def visit(taxa):
        nonlocal call
        if len(taxa) <= 2:
            return
        row = next(nodes)
        if list(row["taxa"]) != list(taxa):
            raise ValueError("Trace is not in recursive visit order.")
        for _ in range(len(taxa) // 2):
            result[call] = tuple(taxa)
            call += 1
        visit(row["left"])
        visit(row["right"])
    visit(labels)
    if next(nodes, None) is not None:
        raise ValueError("Trace contains extra nodes.")
    return result


def tree_from_trace(trace, labels):
    nodes = iter(trace)
    def visit(taxa):
        if len(taxa) <= 2:
            quoted = ["'" + s.replace("'", "''") + "'" for s in taxa]
            return quoted[0] if len(taxa) == 1 else "(" + ",".join(quoted) + ")"
        row = next(nodes)
        if (list(row["taxa"]) != list(taxa) or
                set(row["left"]) & set(row["right"]) or
                sorted(row["left"] + row["right"]) != sorted(taxa)):
            raise ValueError("Invalid recursive trace partition.")
        return "(" + visit(row["left"]) + "," + visit(row["right"]) + ")"
    tree = visit(labels) + ";"
    if next(nodes, None) is not None:
        raise ValueError("Unused trace nodes.")
    return tree


def benchmark_metrics(report):
    """Add diagnostics from actual samples, preserving failures and penalties."""
    labels = report["labels"]
    report["graph_signature"] = graph_signature(report["weights"], labels)
    if report.get("mode") not in ("toy_tree", "protein_tree"):
        return report
    exact_tree, exact_trace = reconstruct(report["weights"], labels)
    report["exact_tree"], report["exact_trace"] = exact_tree, exact_trace
    completed = report.get("status") == "completed"
    taxa = candidate_taxa(report["trace"], labels) if completed else {}
    if completed:
        if tree_from_trace(report["trace"], labels) != report["tree"]:
            raise ValueError("Stored tree disagrees with recursive trace.")
        for row in report["trace"]:
            indices = [labels.index(label) for label in row["taxa"]]
            w = tuple(tuple(report["weights"][i][j] for j in indices) for i in indices)
            bits = tuple(int(label in row["left"]) for label in row["taxa"])
            if (sum(bits) != row["target_size"] or
                    not isclose(cut_weight(w, bits), row["cut"]) or
                    not isclose(normalized_cut(w, bits), row["ncut"])):
                raise ValueError("Stored trace scores disagree with its induced graph.")
        actual = topology_splits(report["trace"], labels)
        expected = topology_splits(exact_trace, labels)
        rf = len(actual ^ expected)
        report["topology_metrics"] = {"unrooted_rf_to_exact": rf,
            "normalized_rf_to_exact": rf / (2 * (len(labels) - 3)) if len(labels) > 3 else 0.0,
            "matches_exact_splits": actual == expected,
            "split_nodes": len(report["trace"])}
    else:
        report["topology_metrics"] = None
    finals = []
    for row in report.get("evaluations", []):
        if row["phase"] != "final":
            continue
        w, c = row["weights"], row["target_size"]
        exact = exact_fixed_size_cut(w, c)
        optimum = cut_weight(w, exact)
        optimal_shots = 0
        for key, count in row["counts"].items():
            bits = decode_bitstring(key, len(w))
            if sum(bits) == c and isclose(cut_weight(w, bits), optimum, abs_tol=1e-9):
                optimal_shots += count
        optimum_states = 0
        for selected in combinations(range(len(w)), c):
            bits = tuple(int(i in selected) for i in range(len(w)))
            optimum_states += isclose(cut_weight(w, bits), optimum, abs_tol=1e-9)
        best = row["summary"]["best_feasible_sample"]
        finals.append({"call": row["call"], "attempt": row["attempt"],
            "taxa": taxa.get(row["call"]), "n": len(w), "target_size": c,
            "penalty": row["penalty"], "exact_cut": optimum,
            "best_feasible_cut_gap": best["cut"] - optimum if best else None,
            "best_feasible_ncut": normalized_cut(w, best["bits"]) if best else None,
            "feasible_rate": row["summary"]["feasible_rate"],
            "optimal_shot_rate": optimal_shots / row["shots"],
            "uniform_feasible_rate": len(list(combinations(range(len(w)), c))) / (2 ** len(w)),
            "uniform_optimal_rate": optimum_states / (2 ** len(w)),
            "expectation": row["summary"]["expectation"],
            "scaled_expectation": row["summary"]["expectation"] / row["scale"]})
    report["final_cut_metrics"] = finals
    report["resource_metrics"] = {"submitted_jobs": len(report.get("jobs", [])),
        "requested_shots": sum(j.get("shots", 0) for j in report.get("jobs", [])),
        "job_elapsed_seconds": sum(j.get("elapsed_seconds", 0) for j in report.get("jobs", []))}
    return report

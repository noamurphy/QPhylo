"""Small, dependency-free reference core for the Qollab phylogeny experiment.

Run: python3 qphylo_core.py --self-test
Demo: python3 qphylo_core.py

Algorithm: minimize cut weight separately for each size c=1..floor(n/2),
compare those candidates by normalized cut, then recurse on induced graphs.
This is NOT a global exact normalized-cut solver or a biological benchmark.
The matrix convention here is symmetric, nonnegative, zero diagonal.
Affinity construction (including BLAST), sequence simulation, NJ, Qiskit,
backend execution, and biological accuracy measurements are not implemented.

Solver interface: solve(weights, c) returns n binary values, exactly c ones.
QAOA must return a feasible measured bitstring; its expected/penalized energy
is not used as the cut score. No silent exact fallback is performed.
On exact ties this reference uses the first lexicographic subset.
"""

from dataclasses import dataclass
from itertools import combinations, product
from math import isclose, isfinite
from typing import Callable, Sequence
import argparse

Matrix = tuple[tuple[float, ...], ...]
Bits = tuple[int, ...]
Solver = Callable[[Matrix, int], Sequence[int]]


def validate_weights(weights: Sequence[Sequence[float]]) -> Matrix:
    """Validate the explicit graph convention; never silently clip affinities."""
    w = tuple(tuple(float(v) for v in row) for row in weights)
    n = len(w)
    if n == 0 or any(len(row) != n for row in w):
        raise ValueError("Weights must be a nonempty square matrix.")
    for i in range(n):
        if w[i][i] != 0:
            raise ValueError("This reference requires zero diagonal weights.")
        for j in range(n):
            if not isfinite(w[i][j]) or w[i][j] < 0:
                raise ValueError("Affinities must be finite and nonnegative.")
            if not isclose(w[i][j], w[j][i], rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError("Affinities must be symmetric.")
    return w


def validate_bits(bits: Sequence[int], n: int, c: int | None = None) -> Bits:
    x = tuple(bits)
    if len(x) != n or any(v not in (0, 1) for v in x):
        raise ValueError("The solver must return one binary value per vertex.")
    if c is not None and sum(x) != c:
        raise ValueError(f"Infeasible solver result: need {c} ones, got {sum(x)}.")
    return tuple(int(v) for v in x)


def cut_weight(weights: Sequence[Sequence[float]], bits: Sequence[int]) -> float:
    w = validate_weights(weights)
    x = validate_bits(bits, len(w))
    return sum(w[i][j] for i in range(len(w)) for j in range(i + 1, len(w))
               if x[i] != x[j])


def normalized_cut(weights: Sequence[Sequence[float]], bits: Sequence[int]) -> float:
    """Use ONLY the current graph's degrees in both volume denominators."""
    w = validate_weights(weights)
    x = validate_bits(bits, len(w))
    if sum(x) in (0, len(w)):
        raise ValueError("A cut must contain two nonempty groups.")
    volumes = [sum(sum(w[i]) for i in range(len(w)) if x[i] == side)
               for side in (0, 1)]
    if min(volumes) <= 0:
        raise ValueError("Normalized cut is undefined for a zero-volume group.")
    return cut_weight(w, x) * (1 / volumes[0] + 1 / volumes[1])


def qubo_energy(weights: Sequence[Sequence[float]], bits: Sequence[int],
                c: int, penalty: float) -> float:
    """Cut(x) + penalty*(sum(x)-c)^2; includes the constant penalty*c*c."""
    w = validate_weights(weights)
    if not 1 <= c < len(w):
        raise ValueError("The target size must lie between 1 and n-1.")
    if not isfinite(penalty) or penalty <= 0:
        raise ValueError("The penalty must be finite and positive.")
    x = validate_bits(bits, len(w))
    return cut_weight(w, x) + penalty * (sum(x) - c) ** 2


def qubo_coefficients(weights: Sequence[Sequence[float]], c: int,
                      penalty: float) -> tuple[float, tuple[float, ...], dict]:
    """Return offset, linear[i], quadratic[(i,j)] for a Qiskit adapter.

    Evaluate as offset + sum(linear[i]*x[i]) + sum(q[i,j]*x[i]*x[j]).
    For binary variables x_i^2=x_i. Preserve the offset when comparing energies.
    """
    w = validate_weights(weights)
    qubo_energy(w, [0] * len(w), c, penalty)  # Validate c and penalty.
    linear = tuple(sum(row) + penalty * (1 - 2 * c) for row in w)
    quadratic = {(i, j): 2 * (penalty - w[i][j])
                 for i in range(len(w)) for j in range(i + 1, len(w))}
    return penalty * c * c, linear, quadratic


def exact_fixed_size_cut(weights: Matrix, c: int) -> Bits:
    """Exponential reference solver, intended only for small taxon sets."""
    w = validate_weights(weights)
    n = len(w)
    if not 1 <= c < n:
        raise ValueError("The target size must lie between 1 and n-1.")
    best_x, best_value = None, float("inf")
    for selected in combinations(range(n), c):
        selected = set(selected)
        x = tuple(int(i in selected) for i in range(n))
        value = cut_weight(w, x)
        if value < best_value:
            best_x, best_value = x, value
    assert best_x is not None
    return best_x


@dataclass(frozen=True)
class Split:
    bits: Bits
    target_size: int
    cut: float
    ncut: float


def choose_split(weights: Sequence[Sequence[float]],
                 solve: Solver = exact_fixed_size_cut) -> Split:
    """Rescore each feasible fixed-size candidate using the current graph."""
    w = validate_weights(weights)
    if len(w) < 2:
        raise ValueError("A split needs at least two vertices.")
    candidates = []
    for c in range(1, len(w) // 2 + 1):
        x = validate_bits(solve(w, c), len(w), c)
        candidates.append(Split(x, c, cut_weight(w, x), normalized_cut(w, x)))
    return min(candidates, key=lambda s: (s.ncut, s.target_size, s.bits))


def reconstruct(weights: Sequence[Sequence[float]], labels: Sequence[str],
                solve: Solver = exact_fixed_size_cut) -> tuple[str, list[dict]]:
    """Return a Newick topology and trace; no branch lengths are estimated.

    The recursive root is a representation of the splitting process, not an
    inferred evolutionary root. Use unrooted metrics for biological evaluation.
    """
    w = validate_weights(weights)
    names = tuple(labels)
    if len(names) != len(w) or any(not isinstance(s, str) or not s for s in names):
        raise ValueError("Provide one nonempty string label per vertex.")
    if len(set(names)) != len(names):
        raise ValueError("Taxon labels must be unique.")
    trace = []

    def visit(current: Matrix, taxa: tuple[str, ...]) -> str:
        if len(taxa) <= 2:
            quoted = ["'" + s.replace("'", "''") + "'" for s in taxa]
            return quoted[0] if len(taxa) == 1 else "(" + ",".join(quoted) + ")"
        split = choose_split(current, solve)
        groups = [tuple(i for i, bit in enumerate(split.bits) if bit == side)
                  for side in (1, 0)]
        trace.append({"taxa": taxa, "target_size": split.target_size,
                      "left": tuple(taxa[i] for i in groups[0]),
                      "right": tuple(taxa[i] for i in groups[1]),
                      "cut": split.cut, "ncut": split.ncut})
        children = []
        for group in groups:
            subgraph = tuple(tuple(current[i][j] for j in group) for i in group)
            children.append(visit(subgraph, tuple(taxa[i] for i in group)))
        return "(" + ",".join(children) + ")"

    return visit(w, names) + ";", trace


def self_test() -> None:
    """Regression checks for scientific scoring and the future solver adapter."""
    full = validate_weights([
        [0, 19, 13, 3, 19, 14], [19, 0, 13, 20, 2, 9],
        [13, 13, 0, 17, 14, 18], [3, 20, 17, 0, 9, 11],
        [19, 2, 14, 9, 0, 18], [14, 9, 18, 11, 18, 0]])
    local = tuple(tuple(full[i][j] for j in range(4)) for i in range(4))
    x1, x2 = exact_fixed_size_cut(local, 1), exact_fixed_size_cut(local, 2)
    assert x1 == (1, 0, 0, 0) and x2 == (1, 1, 0, 0)
    assert isclose(normalized_cut(local, x1), 1.2592592592592593)
    assert isclose(normalized_cut(local, x2), 49 * (1 / 87 + 1 / 83))
    assert choose_split(local).target_size == 2
    # External degrees reverse the ranking, detecting accidental global scoring.
    def wrong_score(x):
        va = sum(sum(full[i]) for i in range(4) if x[i])
        vb = sum(sum(full[i]) for i in range(4) if not x[i])
        return cut_weight(local, x) * (1 / va + 1 / vb)
    assert wrong_score(x1) < wrong_score(x2)

    offset, linear, quadratic = qubo_coefficients(local, 2, 100)
    for x in product((0, 1), repeat=4):
        expanded = offset + sum(a * b for a, b in zip(linear, x))
        expanded += sum(v * x[i] * x[j] for (i, j), v in quadratic.items())
        assert isclose(expanded, qubo_energy(local, x, 2, 100))
        assert isclose(cut_weight(local, x), cut_weight(local, [1 - b for b in x]))
    try:
        choose_split(local, lambda w, c: [0] * len(w))
    except ValueError as exc:
        assert "Infeasible" in str(exc)
    else:
        raise AssertionError("Infeasible quantum samples must be rejected.")

    observed = []
    def recording_solver(w, c):
        observed.append((len(w), c))
        return exact_fixed_size_cut(w, c)
    tree, trace = reconstruct(full, list("ABCDEF"), recording_solver)
    assert tree.endswith(";") and len(trace) >= 2
    assert any(n < 6 for n, c in observed)
    for row in trace:
        group = ["ABCDEF".index(label) for label in row["taxa"]]
        sub = tuple(tuple(full[i][j] for j in group) for i in group)
        x = tuple(int(label in row["left"]) for label in row["taxa"])
        assert isclose(row["ncut"], normalized_cut(sub, x))
    for invalid in ([[0, -1], [-1, 0]], [[0, 1], [2, 0]], [[1, 2], [2, 0]]):
        try:
            validate_weights(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid affinities must be rejected.")
    print("PASS: local normalized cut, QUBO expansion, feasibility, recursive scoring.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    else:
        demo = [[0, 9, 1, 1], [9, 0, 1, 1], [1, 1, 0, 9], [1, 1, 9, 0]]
        tree, trace = reconstruct(demo, ["A", "B", "C", "D"])
        print(tree)
        for row in trace:
            print(row)

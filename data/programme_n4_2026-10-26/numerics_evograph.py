"""N4 exact numerics: fixation probabilities of the birth-death Moran process on weighted digraphs.

Independent of the Lean development. The configuration chain on all 2^n mutant sets is built from its
definition (node i chosen with probability mu_i f_i / sum mu f, then j with probability W_ij, offspring
of i replaces j) and its absorption probabilities into the all-mutant state are obtained by solving the
first-step linear system exactly over the rationals (Fraction Gaussian elimination).

Positive controls: well-mixed W_H = J/n, the undirected 4-cycle, a random rational doubly stochastic
matrix (Birkhoff mixture of permutation matrices); a stationary-selection (non-uniform mu) example.
Negative controls: the star graphs n = 3, 4 (not isothermal), and a random row-stochastic matrix that
is not doubly stochastic. Output: N4/numerics_evograph.json.
"""

from __future__ import annotations

import itertools
import json
import random
from fractions import Fraction as Fr
from pathlib import Path


def update(s: frozenset[int], i: int, j: int) -> frozenset[int]:
    return s | {j} if i in s else s - {j}


def kernel(n: int, w: list[list[Fr]], mu: list[Fr], r: Fr, s: frozenset[int]) -> dict[frozenset[int], Fr]:
    fit = [r if i in s else Fr(1) for i in range(n)]
    tot = sum(mu[i] * fit[i] for i in range(n))
    out: dict[frozenset[int], Fr] = {}
    for i in range(n):
        for j in range(n):
            c = mu[i] * fit[i] / tot * w[i][j]
            if c:
                t = update(s, i, j)
                out[t] = out.get(t, Fr(0)) + c
    return out


def solve(a: list[list[Fr]], b: list[Fr]) -> list[Fr]:
    m = len(b)
    aug = [row[:] + [b[k]] for k, row in enumerate(a)]
    for col in range(m):
        piv = next(k for k in range(col, m) if aug[k][col] != 0)
        aug[col], aug[piv] = aug[piv], aug[col]
        pv = aug[col][col]
        aug[col] = [x / pv for x in aug[col]]
        for k in range(m):
            if k != col and aug[k][col] != 0:
                f = aug[k][col]
                aug[k] = [x - f * y for x, y in zip(aug[k], aug[col])]
    return [aug[k][m] for k in range(m)]


def fixation(n: int, w: list[list[Fr]], mu: list[Fr], r: Fr) -> dict[frozenset[int], Fr]:
    states = [frozenset(c) for k in range(1, n) for c in itertools.combinations(range(n), k)]
    idx = {s: k for k, s in enumerate(states)}
    full = frozenset(range(n))
    a = [[Fr(0)] * len(states) for _ in states]
    b = [Fr(0)] * len(states)
    for s in states:
        k = idx[s]
        a[k][k] += 1
        for t, p in kernel(n, w, mu, r, s).items():
            if t == full:
                b[k] += p
            elif t in idx:
                a[k][idx[t]] -= p
    x = solve(a, b)
    res = {s: x[idx[s]] for s in states}
    res[frozenset()] = Fr(0)
    res[full] = Fr(1)
    return res


def moran(n: int, r: Fr, k: int) -> Fr:
    if r == 1:
        return Fr(k, n)
    return (1 - r ** -k) / (1 - r ** -n)


def forward_bias(n: int, w: list[list[Fr]], mu: list[Fr], r: Fr, s: frozenset[int]) -> tuple[Fr, Fr]:
    up = dn = Fr(0)
    for t, p in kernel(n, w, mu, r, s).items():
        if len(t) == len(s) + 1:
            up += p
        elif len(t) == len(s) - 1:
            dn += p
    return up, dn


def report(name: str, n: int, w: list[list[Fr]], mu: list[Fr], r: Fr) -> dict[str, object]:
    fx = fixation(n, w, mu, r)
    col = [sum(w[i][j] for i in range(n)) for j in range(n)]
    stationary = all(sum(mu[i] * w[i][j] for i in range(n)) == mu[j] for j in range(n))
    singles = {str(v): str(fx[frozenset({v})]) for v in range(n)}
    avg = sum(fx[frozenset({v})] for v in range(n)) / n
    all_moran = all(fx[s] == moran(n, r, len(s)) for s in fx)
    ratios_r = True
    for s in fx:
        if 0 < len(s) < n:
            up, dn = forward_bias(n, w, mu, r, s)
            if up != r * dn:
                ratios_r = False
    return {"name": name, "n": n, "r": str(r), "col_sums": [str(c) for c in col],
            "mu_stationary": stationary, "single_mutant_fixation": singles,
            "uniform_single_average": str(avg), "moran_value_k1": str(moran(n, r, 1)),
            "every_configuration_equals_moran": all_moran, "forward_bias_r_everywhere": ratios_r}


def star(n: int) -> list[list[Fr]]:
    w = [[Fr(0)] * n for _ in range(n)]
    for j in range(1, n):
        w[0][j] = Fr(1, n - 1)
        w[j][0] = Fr(1)
    return w


def main() -> None:
    r = Fr(2)
    out: list[dict[str, object]] = []
    for n in (3, 4):
        out.append(report(f"complete_with_loops_WH_n{n}", n, [[Fr(1, n)] * n for _ in range(n)], [Fr(1)] * n, r))
    cyc = [[Fr(0)] * 4 for _ in range(4)]
    for i in range(4):
        cyc[i][(i + 1) % 4] = Fr(1, 2)
        cyc[i][(i - 1) % 4] = Fr(1, 2)
    out.append(report("cycle4", 4, cyc, [Fr(1)] * 4, r))
    rng = random.Random(20261026)
    n = 4
    perms = list(itertools.permutations(range(n)))
    wts = [Fr(rng.randint(1, 9)) for _ in range(5)]
    chosen = rng.sample(perms, 5)
    ds = [[Fr(0)] * n for _ in range(n)]
    for c, p in zip(wts, chosen):
        for i in range(n):
            ds[i][p[i]] += c / sum(wts)
    out.append(report("random_doubly_stochastic_n4", n, ds, [Fr(1)] * n, r))
    out.append(report("star_n3_NEG", 3, star(3), [Fr(1)] * 3, r))
    out.append(report("star_n4_NEG", 4, star(4), [Fr(1)] * 4, r))
    rs = [[Fr(rng.randint(1, 9)) for _ in range(n)] for _ in range(n)]
    rs = [[x / sum(row) for x in row] for row in rs]
    out.append(report("random_row_stochastic_n4_NEG", n, rs, [Fr(1)] * n, r))
    # stationary (non-uniform) selection: star with mu = stationary distribution of the star walk
    mu3 = [Fr(1, 2), Fr(1, 4), Fr(1, 4)]
    out.append(report("star_n3_stationary_selection", 3, star(3), mu3, r))
    # star n = 3 full table of configuration fixation probabilities (for the Lean negative control)
    fx = fixation(3, star(3), [Fr(1)] * 3, r)
    table = {",".join(map(str, sorted(s))) or "empty": str(v) for s, v in sorted(fx.items(), key=lambda kv: (len(kv[0]), sorted(kv[0])))}
    cyc_bias = {name: [str(x) for x in forward_bias(4, cyc, [Fr(1)] * 4, r, frozenset(s))]
                for name, s in (("0,1", {0, 1}), ("0,2", {0, 2}))}
    path = Path(__file__).with_name("numerics_evograph.json")
    path.write_text(json.dumps({"r": "2", "cases": out, "star_n3_table": table,
                                "cycle4_up_down_two_mutants": cyc_bias,
                                "star_n3_center_bias": [str(x) for x in forward_bias(3, star(3), [Fr(1)] * 3, r, frozenset({0}))]}, indent=2) + "\n")
    print(path.read_text())


if __name__ == "__main__":
    main()

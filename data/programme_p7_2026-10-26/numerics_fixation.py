"""Tier X cross-check for lane P7 (never supports a claim; the claims are the Lean theorems).

1. Exact rational (fractions.Fraction) positive controls matching the Lean controls:
   WF N=2: P^2(1,2) = 3/8; Moran N=3, r=2: solve the first-step linear system exactly -> 4/7, 6/7.
   Exact solve of the Moran first-step system for N=10, r=3/2 compared with (1-r^-i)/(1-r^-N).
2. mpmath (60 digits) matrix powers P^(2^k) for WF and Moran (r in {3/2, 4/5, 1}), N=10:
   |P^n(i,N) - formula| should shrink to ~0 while the interior mass shrinks.
3. Negative controls: the exponent-flipped formula (1-r^i)/(1-r^N) and the neutral i/N under
   selection must NOT match the matrix-power limit.
Writes numerics_fixation.json next to this script.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

import mpmath as mp

mp.mp.dps = 60
OUT = Path(__file__).resolve().parent / "numerics_fixation.json"


def wf_matrix(n_pop: int, exact: bool) -> list[list]:
    one = Fraction(1) if exact else mp.mpf(1)
    rows = []
    for i in range(n_pop + 1):
        x = (Fraction(i, n_pop) if exact else mp.mpf(i) / n_pop)
        rows.append([comb(n_pop, j) * x ** j * (one - x) ** (n_pop - j) for j in range(n_pop + 1)])
    return rows


def moran_rates(n_pop: int, r, i: int, exact: bool) -> tuple:
    if exact:
        ii, nn = Fraction(i), Fraction(n_pop)
    else:
        ii, nn = mp.mpf(i), mp.mpf(n_pop)
    den = r * ii + (nn - ii)
    up = r * ii / den * ((nn - ii) / nn)
    down = (nn - ii) / den * (ii / nn)
    return up, down


def moran_matrix(n_pop: int, r, exact: bool) -> list[list]:
    zero = Fraction(0) if exact else mp.mpf(0)
    one = Fraction(1) if exact else mp.mpf(1)
    rows = []
    for i in range(n_pop + 1):
        up, down = moran_rates(n_pop, r, i, exact)
        row = [zero] * (n_pop + 1)
        if i + 1 <= n_pop:
            row[i + 1] += up
        if i - 1 >= 0:
            row[i - 1] += down
        row[i] += one - up - down
        rows.append(row)
    return rows


def matmul(a: list[list], b: list[list]) -> list[list]:
    n = len(a)
    return [[sum(a[i][k] * b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def solve_first_step_exact(p: list[list]) -> list[Fraction]:
    """Solve h(0)=0, h(N)=1, h(i)=sum_j P(i,j)h(j) for 0<i<N by Gaussian elimination over Q."""
    n = len(p) - 1
    m = n - 1
    a = [[Fraction(0)] * m for _ in range(m)]
    b = [Fraction(0)] * m
    for row in range(m):
        i = row + 1
        for j in range(1, n):
            a[row][j - 1] = (Fraction(1) if i == j else Fraction(0)) - p[i][j]
        b[row] = p[i][n]
    for col in range(m):
        piv = next(rr for rr in range(col, m) if a[rr][col] != 0)
        a[col], a[piv] = a[piv], a[col]
        b[col], b[piv] = b[piv], b[col]
        for rr in range(m):
            if rr != col and a[rr][col] != 0:
                f = a[rr][col] / a[col][col]
                a[rr] = [x - f * y for x, y in zip(a[rr], a[col])]
                b[rr] -= f * b[col]
    return [Fraction(0)] + [b[k] / a[k][k] for k in range(m)] + [Fraction(1)]


def closed(n_pop: int, r, i: int):
    return (1 - r ** (-i)) / (1 - r ** (-n_pop))


def flipped(n_pop: int, r, i: int):
    return (1 - r ** i) / (1 - r ** n_pop)


def main() -> int:
    res: dict[str, object] = {"tier": "X", "note": "cross-check only; claims rest on the Lean kernel"}
    # 1. exact positive controls
    p2 = wf_matrix(2, exact=True)
    res["wf_N2_P2_1_2"] = str(matmul(p2, p2)[1][2])
    h3 = solve_first_step_exact(moran_matrix(3, Fraction(2), exact=True))
    res["moran_N3_r2_exact"] = [str(x) for x in h3]
    r = Fraction(3, 2)
    h10 = solve_first_step_exact(moran_matrix(10, r, exact=True))
    res["moran_N10_r3_2_exact_equals_closed"] = all(h10[i] == closed(10, r, i) for i in range(11))
    res["moran_N10_r3_2_exact_equals_flipped"] = all(h10[i] == flipped(10, r, i) for i in range(11))
    # 2. matrix powers
    runs = []
    cases = [("WF", None), ("Moran", mp.mpf(3) / 2), ("Moran", mp.mpf(4) / 5), ("Moran", mp.mpf(1))]
    for model, rr in cases:
        n_pop = 10
        mat = wf_matrix(n_pop, exact=False) if model == "WF" else moran_matrix(n_pop, rr, exact=False)
        power, steps = mat, 1
        for _ in range(12):  # up to n = 2^12 = 4096 steps
            power, steps = matmul(power, power), steps * 2
        for i in (1, 3, 7):
            target = mp.mpf(i) / n_pop if (model == "WF" or rr == 1) else closed(n_pop, rr, i)
            interior = sum(power[i][k] for k in range(1, n_pop))
            row = {"model": model, "r": None if rr is None else float(rr), "N": n_pop, "i": i,
                   "n_steps": steps, "PnN": mp.nstr(power[i][n_pop], 20),
                   "formula": mp.nstr(target, 20),
                   "abs_err": mp.nstr(abs(power[i][n_pop] - target), 5),
                   "interior_mass": mp.nstr(interior, 5)}
            if model == "Moran" and rr != 1:
                row["neg_flipped_abs_diff"] = mp.nstr(abs(power[i][n_pop] - flipped(n_pop, rr, i)), 5)
                row["neg_neutral_abs_diff"] = mp.nstr(abs(power[i][n_pop] - mp.mpf(i) / n_pop), 5)
            runs.append(row)
    res["matrix_power_runs"] = runs
    OUT.write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Tier X cross-check for the P7 Price / Fisher Lean file (never supports a claim).

Exact-rational (fractions.Fraction) evaluation of both sides of the Price equation and of the
haploid selection increment on (a) the explicit instances used as Lean controls and (b) seeded
random instances; plus a float64 run to show the size of rounding error. Writes price_numerics.json.
"""

from __future__ import annotations

import json
import random
from fractions import Fraction
from pathlib import Path
from typing import Sequence

Num = Fraction


def wmean(n: Sequence[Num], f: Sequence[Num]) -> Num:
    return sum((a * b for a, b in zip(n, f)), Num(0)) / sum(n, Num(0))


def wcov(n: Sequence[Num], f: Sequence[Num], g: Sequence[Num]) -> Num:
    mf, mg = wmean(n, f), wmean(n, g)
    return sum((a * (b - mf) * (c - mg) for a, b, c in zip(n, f, g)), Num(0)) / sum(n, Num(0))


def zbar_next(n: Sequence[Num], w: Sequence[Num], zp: Sequence[Num]) -> Num:
    return wmean([a * b for a, b in zip(n, w)], zp)


def price_sides(n: Sequence[Num], w: Sequence[Num], z: Sequence[Num], zp: Sequence[Num]) -> tuple[Num, Num]:
    lhs = wmean(n, w) * (zbar_next(n, w, zp) - wmean(n, z))
    rhs = wcov(n, w, z) + wmean(n, [b * (d - c) for b, c, d in zip(w, z, zp)])
    return lhs, rhs


def fisher_sides(n: Sequence[Num], w: Sequence[Num], wp: Sequence[Num]) -> tuple[Num, Num]:
    """Delta wbar (offspring fitness wp) vs Var(w)/wbar."""
    return zbar_next(n, w, wp) - wmean(n, w), wcov(n, w, w) / wmean(n, w)


def F(xs: Sequence[int]) -> list[Num]:
    return [Num(x) for x in xs]


def main() -> None:
    out: dict[str, object] = {"tier": "X", "note": "never supports a claim"}
    # (a) Lean control instances
    lhs, rhs = price_sides(F([1, 2, 3]), F([1, 2, 4]), F([0, 1, 3]), F([1, 1, 2]))
    out["pos_three_types"] = {"lhs": str(lhs), "rhs": str(rhs), "equal": lhs == rhs}
    n2, w2 = F([1, 1]), F([1, 2])
    d_bias, v_bias = fisher_sides(n2, w2, F([1, 1]))
    pl, pr = price_sides(n2, w2, w2, F([1, 1]))
    out["neg_transmission_bias"] = {
        "delta_wbar": str(d_bias), "var_over_wbar": str(v_bias), "equal": d_bias == v_bias,
        "price_lhs": str(pl), "price_rhs": str(pr), "price_equal": pl == pr,
    }
    d0, v0 = fisher_sides(n2, w2, w2)
    out["neg_missing_wbar"] = {
        "delta_wbar": str(d0), "var_over_wbar": str(v0), "var": str(wcov(n2, w2, w2)),
        "delta_equals_var_over_wbar": d0 == v0, "delta_equals_var": d0 == wcov(n2, w2, w2),
    }
    # (b) seeded random exact instances
    rng = random.Random(20261025)
    ok_price = ok_fisher = 0
    trials = 500
    for _ in range(trials):
        k = rng.randint(1, 8)
        n = [Num(rng.randint(1, 50), rng.randint(1, 9)) for _ in range(k)]
        w = [Num(rng.randint(1, 50), rng.randint(1, 9)) for _ in range(k)]
        z = [Num(rng.randint(-30, 30), rng.randint(1, 9)) for _ in range(k)]
        zp = [Num(rng.randint(-30, 30), rng.randint(1, 9)) for _ in range(k)]
        a, b = price_sides(n, w, z, zp)
        ok_price += a == b
        c, d = fisher_sides(n, w, w)
        ok_fisher += c == d
    out["random_exact"] = {"trials": trials, "price_equal": ok_price, "fisher_equal": ok_fisher}
    # (c) float64 rounding size on the same seed
    rng = random.Random(7)
    worst = 0.0
    for _ in range(trials):
        k = rng.randint(1, 8)
        nf = [rng.uniform(0.1, 10) for _ in range(k)]
        wf = [rng.uniform(0.1, 10) for _ in range(k)]
        zf = [rng.uniform(-5, 5) for _ in range(k)]
        zpf = [rng.uniform(-5, 5) for _ in range(k)]
        N = sum(nf)
        A = sum(a * b for a, b in zip(nf, wf))
        mw, mz = A / N, sum(a * c for a, c in zip(nf, zf)) / N
        lhs_f = mw * (sum(a * b * d for a, b, d in zip(nf, wf, zpf)) / A - mz)
        rhs_f = sum(a * (b - mw) * (c - mz) for a, b, c in zip(nf, wf, zf)) / N + sum(
            a * b * (d - c) for a, b, c, d in zip(nf, wf, zf, zpf)) / N
        worst = max(worst, abs(lhs_f - rhs_f))
    out["float64_worst_abs_residual"] = worst
    path = Path(__file__).with_name("price_numerics.json")
    path.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

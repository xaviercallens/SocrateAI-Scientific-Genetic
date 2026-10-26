"""P7 auditor's independent number trace and spot re-computation (Tier X, never supports a claim).

1. Recompute every gen_numbers.tex macro and table_cvode.tex row from data/results.jsonl.
2. Exact Fraction values S_i/S_N for every accepted Moran row (and i/N for Wright-Fisher): compare with p_fix.
3. Data-certificate literals of FixationDataCert.lean: find them verbatim in results.jsonl, recompute errors.
4. Independent matrix power: exact Fraction P^n for Moran N=3, r=2 (n up to 400) and mpmath 50-digit P^2048 for
   Moran N=10, r=3/2, i=1 (own kernel constructor). Positive control: limit 4/7; negative controls: 1/3, 1/7.
5. Price positive control -13/36 and both negative controls, exact.
"""
import json
import re
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp

P = Path("/home/callensxavier_gmail_com/AutoevolveAI/.claude/worktrees/openai-math-discovery/papers/prog_population_genetics_fixation")
rows = [json.loads(l) for l in (P / "data/results.jsonl").read_text().splitlines()]
acc = [r for r in rows if r["verdict"] == "accepted"]
ref = [r for r in rows if r["verdict"] == "refused"]
out = {}
out["macros"] = {
    "nRows": len(rows), "nAcc": len(acc), "nRef": len(ref),
    "maxErr": f"{max(r['abs_err'] for r in acc):.1e}",
    "maxInterior": f"{max(r['interior'] for r in acc):.1e}",
    "maxDrift": f"{max(r['mass_drift'] for r in acc):.1e}",
    "nMoranRuns": sum(r["kind"] == "moran" for r in acc),
    "nWFRuns": sum(r["kind"] == "wright_fisher" for r in acc),
}
out["gen_numbers_tex"] = (P / "gen/gen_numbers.tex").read_text()
by = {}
for r in acc:
    by.setdefault(r["label"], []).append(r)
out["table_rows"] = {k: {"N_values": sorted({x["N"] for x in v}), "r": v[0]["r"], "starts": len(v),
                         "max_err": f"{max(x['abs_err'] for x in v):.1e}"} for k, v in by.items()}


def rat(x: float) -> F:
    return F(x).limit_denominator(10**6)


def moran_exact(n: int, r: F, i: int) -> F:
    s = lambda m: sum((1 / r) ** k for k in range(m))
    return s(i) / s(n)


worst = 0.0
for r in acc:
    ex = F(r["i"], r["N"]) if r["kind"] == "wright_fisher" else moran_exact(r["N"], rat(r["r"]), r["i"])
    worst = max(worst, abs(float(F(r["p_fix"]) - ex)))
out["max_abs_pfix_minus_exact_fraction"] = worst

src = (P.parents[1] / "formal_cert/Programmes/PopGen/FixationDataCert.lean").read_text()
lits = re.findall(r"\|\(([0-9.]+) : ℝ\) - moranFix (\d+) \(?([0-9/ ]+?)\)? (\d+)\|", src)
cert = []
for lit, n, r, i in lits:
    rr = F(r.replace(" ", ""))
    ex = moran_exact(int(n), rr, int(i))
    found = [x for x in acc if x["kind"] == "moran" and x["N"] == int(n) and rat(x["r"]) == rr and x["i"] == int(i)]
    cert.append({"literal": lit, "N": n, "r": str(rr), "i": i, "in_results": any(repr(x["p_fix"]) == lit or
                 f"{x['p_fix']!r}" == lit or str(x["p_fix"]) == lit for x in found),
                 "p_fix_in_results": [x["p_fix"] for x in found], "err_vs_exact": f"{abs(float(F(lit) - ex)):.2e}"})
out["data_cert"] = cert


def moran_kernel_fr(n: int, r: F) -> list:
    k = [[F(0)] * (n + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        up = r * i / (r * i + n - i) * F(n - i, n)
        dn = F(n - i) / (r * i + n - i) * F(i, n)
        if i < n:
            k[i][i + 1] = up
        if i > 0:
            k[i][i - 1] = dn
        k[i][i] = 1 - up - dn
    return k


K = moran_kernel_fr(3, F(2))
v = [F(0), F(1), F(0), F(0)]
traj = {}
for step in range(1, 401):
    v = [sum(v[k] * K[k][j] for k in range(4)) for j in range(4)]
    if step in (1, 2, 10, 50, 100, 400):
        traj[step] = float(v[3])
out["moran_N3_r2_i1_exact_Pn_1_3"] = traj
out["moran_N3_r2_i1_gap_to_4_7_at_400"] = float(F(4, 7) - v[3])
out["moran_N3_r2_monotone_check_note"] = "Pn(1,3) values listed above should be non-decreasing and approach 4/7=0.571428..."

mp.mp.dps = 50
n = 10
Kf = moran_kernel_fr(n, F(3, 2))
M = mp.matrix([[mp.mpf(x.numerator) / x.denominator for x in row] for row in Kf])
Mp = M ** 2048
pn = Mp[1, n]
ex = moran_exact(10, F(3, 2), 1)
out["moran_N10_r1p5_i1_P2048"] = mp.nstr(pn, 25)
out["moran_N10_r1p5_i1_exact"] = mp.nstr(mp.mpf(ex.numerator) / ex.denominator, 25)
out["moran_N10_r1p5_i1_abs_err"] = mp.nstr(abs(pn - mp.mpf(ex.numerator) / ex.denominator), 5)
out["neg_neutral_gap"] = mp.nstr(abs(pn - mp.mpf(1) / 10), 5)
out["neg_flipped_gap"] = mp.nstr(abs(pn - (1 - mp.mpf(1.5)) / (1 - mp.mpf(1.5) ** 10)), 5)


def wmean(nn, f):
    return sum(a * b for a, b in zip(nn, f)) / sum(nn)


def wcov(nn, f, g):
    mf, mg = wmean(nn, f), wmean(nn, g)
    return sum(a * (b - mf) * (c - mg) for a, b, c in zip(nn, f, g)) / sum(nn)


def zbar_next(nn, w, zp):
    return wmean([a * b for a, b in zip(nn, w)], zp)


nn, w, z, zp = [F(1), F(2), F(3)], [F(1), F(2), F(4)], [F(0), F(1), F(3)], [F(1), F(1), F(2)]
lhs = wmean(nn, w) * (zbar_next(nn, w, zp) - wmean(nn, z))
rhs = wcov(nn, w, z) + wmean(nn, [a * (b - c) for a, b, c in zip(w, zp, z)])
out["price_pos"] = [str(lhs), str(rhs)]
nn, w = [F(1), F(1)], [F(1), F(2)]
out["price_neg_transmission"] = {"delta": str(zbar_next(nn, w, [F(1), F(1)]) - wmean(nn, w)),
                                 "var_over_wbar": str(wcov(nn, w, w) / wmean(nn, w))}
out["price_neg_missing_wbar"] = {"delta": str(zbar_next(nn, w, w) - wmean(nn, w)), "var": str(wcov(nn, w, w))}
txt = json.dumps(out, indent=1)
Path("/home/callensxavier_gmail_com/programmes-2026-10-25/P7/audit/audit_numbers.json").write_text(txt)
print(txt)

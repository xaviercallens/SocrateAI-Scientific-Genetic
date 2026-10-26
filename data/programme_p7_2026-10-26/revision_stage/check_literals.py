"""Check C3/C7 facts: raw JSONL text of p_fix vs Lean literals; exact errors of decimal literals; r=1 mpmath error."""
from __future__ import annotations

import json
import re
from fractions import Fraction as F
from pathlib import Path

W = Path("/home/callensxavier_gmail_com/AutoevolveAI/.claude/worktrees/openai-math-discovery")
P = W / "papers/prog_population_genetics_fixation/data"
raw = (P / "results.jsonl").read_text().splitlines()
lits = {(2, 2.0, 1): "0.6666666666666601", (3, 2.0, 1): "0.5714285714284759", (3, 2.0, 2): "0.8571428571427792",
        (10, 1.5, 1): "0.33921585523416337", (10, 1.5, 5): "0.8836363636355173", (10, 1.0, 3): "0.29999999999999977"}


def exact(n: int, r: float, i: int) -> F:
    if r == 1.0:
        return F(i, n)
    rr = F(r).limit_denominator(100)
    s = lambda m: sum(rr ** (-k) for k in range(m))
    return s(i) / s(n)


out = {}
for line in raw:
    d = json.loads(line)
    if d["verdict"] != "accepted" or d["kind"] != "moran":
        continue
    key = (d["N"], float(d["r"]), d["i"])
    if key in lits:
        m = re.search(r'"p_fix":([^,]+),', line)
        lit = lits[key]
        assert repr(d["p_fix"]) == lit, (key, repr(d["p_fix"]), lit)
        ex = exact(*key)
        out[str(key)] = {"raw_jsonl_p_fix_text": m.group(1), "python_repr": repr(d["p_fix"]), "lean_literal": lit,
                         "decimal_literal_minus_exact": float(abs(F(lit) - ex)),
                         "float64_p_fix_minus_float64_exact": abs(d["p_fix"] - float(ex)),
                         "jsonl_abs_err": d["abs_err"]}
print(json.dumps(out, indent=1))
nm = json.loads((P / "numerics_fixation.json").read_text())
vals = [r for r in nm["matrix_power_runs"] if r.get("r") == 1.0 and r["model"] == "Moran"]
print([(v["N"], v["i"], v["abs_err"]) for v in vals])
(Path("/home/callensxavier_gmail_com/programmes-2026-10-25/P7/revise/literals.json")).write_text(json.dumps(out, indent=1))

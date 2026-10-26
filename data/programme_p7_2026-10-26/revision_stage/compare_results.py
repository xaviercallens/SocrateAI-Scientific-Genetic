"""Compare the re-run results.jsonl (with min_component, identity r null) against the v1 file field by field."""
from __future__ import annotations

import json
from pathlib import Path

W = Path("/home/callensxavier_gmail_com/AutoevolveAI/.claude/worktrees/openai-math-discovery")
old = [json.loads(x) for x in (W / "papers/prog_population_genetics_fixation/data/results.jsonl").read_text().splitlines()]
new = [json.loads(x) for x in Path("/home/callensxavier_gmail_com/programmes-2026-10-25/P7/revise/cvode_out/results.jsonl").read_text().splitlines()]
assert len(old) == len(new) == 115
diffs = []
mins = []
for a, b in zip(old, new):
    for k in set(a) | set(b):
        if a.get(k, "<absent>") != b.get(k, "<absent>"):
            diffs.append((a["label"], a["i"], k, a.get(k, "<absent>") if k != "trace" else "trace", b.get(k, "<absent>") if k != "trace" else "trace"))
    if "min_component" in b:
        mins.append(b["min_component"])
from collections import Counter

c = Counter((d[2], str(d[3]), str(d[4])) for d in diffs if d[2] != "min_component")
print("differing fields other than min_component:", dict(c))
print("rows with min_component:", sum(1 for d in diffs if d[2] == "min_component"), "of", len(new))
print("min over all accepted runs of min_component:", min(mins), "max:", max(mins))
print("identity row new r:", [b["r"] for b in new if b["label"] == "identity_kernel_no_absorption"])
out = {"rows": len(new), "fields_differing_other_than_min_component": {f"{k[0]}|{k[1]}->{k[2]}": v for k, v in c.items()},
       "min_component_min": min(mins), "min_component_max": max(mins), "n_with_min_component": len(mins)}
Path("/home/callensxavier_gmail_com/programmes-2026-10-25/P7/revise/compare_results.json").write_text(json.dumps(out, indent=1))

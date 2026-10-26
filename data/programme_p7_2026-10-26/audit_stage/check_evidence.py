"""P7 audit: compare every sha256 recorded inside ledger evidence blobs with the current file on disk."""
import hashlib
import json
from pathlib import Path

E = Path("/home/callensxavier_gmail_com/AutoevolveAI/.claude/worktrees/openai-math-discovery/papers/prog_population_genetics_fixation/ledger/evidence")
ok = bad = missing = 0


def walk(o, out):
    if isinstance(o, dict):
        if "path" in o and "sha256" in o:
            out.append((o["path"], o["sha256"]))
        for v in o.values():
            walk(v, out)
    elif isinstance(o, list):
        for v in o:
            walk(v, out)


for f in sorted(E.glob("*.json")):
    refs = []
    walk(json.loads(f.read_text()), refs)
    for p, h in refs:
        q = Path(p)
        if not q.exists():
            missing += 1
            print("MISSING", f.name[:12], p)
            continue
        cur = hashlib.sha256(q.read_bytes()).hexdigest()
        if cur == h.replace("sha256:", ""):
            ok += 1
        else:
            bad += 1
            print("MISMATCH", f.name[:12], p, h[:16], cur[:16])
print(json.dumps({"refs_ok": ok, "mismatch": bad, "missing": missing}))

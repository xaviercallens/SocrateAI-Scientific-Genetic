"""P7 auditor's independent Lean recompile (producer != verifier).

Builds every P7 module and its H11/H12 dependencies FROM SOURCE into the auditor's own olean directory
(the producer's formal_cert/.olean_out is NOT on LEAN_PATH), then compiles the import-all, controls and
data-certificate files, and the auditor's own planted / mutation files. Records rc, output and sha256.
"""
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

W = Path("/home/callensxavier_gmail_com/AutoevolveAI/.claude/worktrees/openai-math-discovery")
FC = W / "formal_cert"
PINNED = Path("/mnt/disks/disk-socrateai-local-1/callensxavier_home_data/openai_math_pinned/lean")
BIN = Path("/mnt/disks/disk-socrateai-local-1/callensxavier_home_data/elan/toolchains/leanprover--lean4---v4.34.1/bin")
A = Path("/home/callensxavier_gmail_com/programmes-2026-10-25/P7/audit")
MYOUT = A / "olean"
OUTD = A / "out"

env0 = dict(os.environ, PATH=f"{BIN}:{os.environ.get('PATH', '')}")
lp = subprocess.run(["lake", "env", "printenv", "LEAN_PATH"], cwd=PINNED, env=env0, capture_output=True, text=True,
                    check=True).stdout.strip()
assert ".olean_out" not in lp, lp
env = dict(env0, LEAN_PATH=lp + ":" + str(MYOUT))


def run(src: Path, root: Path, module: str | None, tag: str) -> dict:
    cmd = [str(BIN / "lean"), "-DautoImplicit=false", f"--root={root}"]
    if module:
        ol = MYOUT / (module.replace(".", "/") + ".olean")
        ol.parent.mkdir(parents=True, exist_ok=True)
        cmd += ["-o", str(ol), "-i", str(ol.with_suffix(".ilean"))]
    cmd.append(str(src))
    t0 = time.time()
    r = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=3000)
    res = {"tag": tag, "file": str(src), "module": module, "rc": r.returncode, "seconds": round(time.time() - t0, 1),
           "sha256": hashlib.sha256(src.read_bytes()).hexdigest(), "output": r.stdout + r.stderr}
    (OUTD / f"{tag}.json").write_text(json.dumps(res, indent=1))
    ax = [l for l in res["output"].splitlines() if "depends on axioms" in l]
    errs = [l for l in res["output"].splitlines() if ": error" in l]
    print(tag, res["rc"], res["seconds"], "axiom_lines", len(ax), "errors", len(errs), flush=True)
    return res


builds = [
    ("Autoresearch/PopGenWrightFisher.lean", "Autoresearch.PopGenWrightFisher"),
    ("Autoresearch/PopGenFisher.lean", "Autoresearch.PopGenFisher"),
    ("Programmes/PopGen/Absorption.lean", "Programmes.PopGen.Absorption"),
    ("Programmes/PopGen/MoranModel.lean", "Programmes.PopGen.MoranModel"),
    ("Programmes/PopGen/WrightFisherFixation.lean", "Programmes.PopGen.WrightFisherFixation"),
    ("Programmes/PopGen/MoranFixation.lean", "Programmes.PopGen.MoranFixation"),
    ("Programmes/PopGen/Price.lean", "Programmes.PopGen.Price"),
]
for rel, mod in builds:
    r = run(FC / rel, FC, mod, "build_" + mod.split(".")[-1])
    if r["rc"] != 0:
        print("BUILD FAILED", mod, r["output"][-3000:])
        sys.exit(1)
for f in ["PopGenProgImportAll", "PriceImportAll", "FixationDataCert", "PopGenProgControls", "PriceControls",
          "FixationDataCertControls"]:
    run(FC / "Programmes/PopGen" / f"{f}.lean", FC, None, "compile_" + f)
for f in ["AuditAxioms", "AuditPlanted", "AuditMutation"]:
    run(A / "src" / f"{f}.lean", A / "src", None, "audit_" + f)

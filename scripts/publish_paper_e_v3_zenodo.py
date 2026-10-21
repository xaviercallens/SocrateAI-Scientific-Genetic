"""Deposit Paper E edition 3 as a new version of Zenodo record 23525664.

Edition 3 adds a machine-checked no-go theorem for proofreading-stage-removal models and its
refutation by published yeast data. The scripted "Opus 5-round peer review" JSON is not, and
will not be, part of any bundle.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import urllib.request
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[2]
ZENODO = "https://zenodo.org/api"
PREVIOUS = 23525664
RECEIPT = HERE / "data" / "zenodo_paper_e_v3_receipt.json"

STEM = "PaperE_Genetics_DNAReplication_ART_LeanFlow_RustySundials_Lean4"
FIGURES = [
    "fig1_hopfield_art_proofreading_and_kjma.pdf",
    "fig2_calugareanu_supercoiling_rrk_and_hecke.pdf",
    "fig3_c1_c2_torque_haar_art_convergence.pdf",
    "fig4_e8_dicodon_chimeric_art_exo_humanitarian_538x.pdf",
    "figE5_stage_removal_feasible_surface.pdf",
    "figE6_variant_panel_vs_ceiling.pdf",
    "figE7_in_vitro_in_vivo_reversal.pdf",
]
DATA = [
    "paper_e_telemetry.json",
    "tcga_pole_empirical_check.json",
    "pole_proofreading_literature.json",
    "proofreading_stage_removal_test.json",
]
SCRIPTS = [
    "run_tcga_pole_empirical_check.py",
    "run_proofreading_stage_removal_test.py",
    "make_proofreading_figures.py",
]
LEAN = [
    "ProofreadingStageRemoval.lean",
    "Breakthroughs3DomainsSpec.lean",
    "Breakthroughs3Domains.lean",
    "GeneticAutoresearchSpec.lean",
    "GeneticAutoresearch.lean",
    "StochasticThermodynamicsSpec.lean",
    "StochasticThermodynamics.lean",
    "LeanMasterAutoTactics.lean",
]

TITLE = (
    "Non-Equilibrium Hopfield-Ninio Proofreading, Chimeric ART-Exo Fidelity, Gosset E8 Dicodon "
    "Quantization, and Topological DNA Supercoiling in LeanFlow, rusty-SUNDIALS, and Lean 4"
)
DESCRIPTION = (
    "<p>Edition 3. A Lean 4 and exact-rational model of kinetic proofreading, torque-coupled "
    "fidelity, a chimeric reverse transcriptase with an exonuclease stage, E8 dicodon quantization "
    "and DNA supercoiling.</p>"
    "<p><b>New in this edition: a no-go theorem, and the refutation of this paper's own mechanism "
    "for POLE P286R.</b> If a replicative error rate is an insertion/extension error multiplied by "
    "one escape factor in [0,1] per proofreading stage, then no mutant that merely removes stages "
    "can exceed the stage-free (exonuclease-null) error rate, whatever the stage count and whatever "
    "the per-stage factors (theorem PSR_2, formal/ANSE/ProofreadingStageRemoval.lean; 7 declarations, "
    "every one carrying a #print axioms line, no sorryAx, axioms within {propext, Classical.choice, "
    "Quot.sound}). Kane and Shcherbakova (Cancer Res 2014, Table 2) measured the yeast analogue of "
    "POLE P286R, pol2-P301R, and the exonuclease-null allele pol2-4 in the same table and the same "
    "isogenic mismatch-repair-proficient background: 2300 versus 44 (x10^-8, CAN1), i.e. 52.3-fold "
    "above the ceiling, at least 35.3-fold on the published 95% confidence intervals. Three further "
    "cancer-associated alleles (P286H, F367S, S459F) break the same ceiling. Theorem GAR.H15's "
    "treatment of P286R as the loss of two of three proofreading stages is therefore refuted as a "
    "mechanism, not merely mis-parameterised. Negative controls behave as required: V411L and six "
    "functionally neutral alleles are not flagged.</p>"
    "<p><b>Mechanistic consequence.</b> In the purified system without mismatch repair, P301R is "
    "more accurate than the exonuclease-null enzyme (59 versus 155 x10^-6 base substitutions, Xing "
    "et al. 2019), yet far worse in the cell: the ordering reverses by a factor of 137, so intrinsic "
    "error production cannot be the explanatory variable and a mismatch-repair interaction term is "
    "required. A tension is reported rather than smoothed over: pol2-P301R msh6-delta is synthetic "
    "lethal while pol2-4 msh6-delta is viable, which argues against complete evasion. The decisive "
    "experiment, a rate measurement in a hypomorphic rather than null mismatch-repair background, is "
    "identified but has not been performed.</p>"
    "<p>Also in this edition: all 24 declarations of StochasticThermodynamics.lean now verify "
    "individually (edition 2 left 18 helper lemmas without their own #print axioms line).</p>"
    "<p><b>Carried over from edition 2.</b> The description of earlier editions as peer reviewed is "
    "withdrawn: the bundled 'Opus 5-round peer review' files were produced by a script that calls no "
    "reviewer, and they are not distributed. The constant 0.915 'epistemic score' is removed. The "
    "538x ANSE ratio and the clinical figures are labelled model-internal. The ART discovery is cited "
    "as Yoon et al. 2026, correcting a non-existent bioRxiv reference. <b>This paper has not been "
    "peer reviewed.</b></p>"
)
KEYWORDS = [
    "kinetic proofreading", "DNA polymerase epsilon", "POLE P286R", "Lean 4",
    "formal verification", "no-go theorem", "falsification", "mismatch repair",
]


def token() -> str:
    spec = importlib.util.spec_from_file_location(
        "p", REPO / "scripts/three_breakthroughs/publish_3breakthroughs_to_zenodo.py")
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return str(mod.load_zenodo_token())


def req(method: str, url: str, tok: str, data: bytes | None = None,
        ctype: str = "application/json") -> dict[str, Any]:
    r = urllib.request.Request(
        url, data=data, method=method,
        headers={"Authorization": f"Bearer {tok}", "Content-Type": ctype})
    with urllib.request.urlopen(r, timeout=600) as resp:
        body = resp.read()
        return json.loads(body) if body else {}


def bundle() -> list[Path]:
    files = [HERE / "papers" / f"{STEM}.pdf", HERE / "papers" / f"{STEM}.tex"]
    files += [HERE / "papers" / f for f in FIGURES]
    files += [HERE / "data" / f for f in DATA]
    files += [HERE / "scripts" / f for f in SCRIPTS]
    files += [HERE / "formal" / "ANSE" / f for f in LEAN]
    files += [
        HERE / "docs" / "AUDIT_2026-10-21.md",
        HERE / "tests" / "test_proofreading_stage_removal.py",
    ]
    missing = [str(f) for f in files if not f.exists()]
    assert not missing, missing
    forbidden = [f for f in files if "peer_review" in f.name]
    assert not forbidden, forbidden
    return files


def main() -> int:
    tok = token()
    prev = req("GET", f"{ZENODO}/deposit/depositions/{PREVIOUS}", tok)
    assert prev["state"] == "done", prev["state"]

    nv = req("POST", f"{ZENODO}/deposit/depositions/{PREVIOUS}/actions/newversion", tok)
    dep = req("GET", nv["links"]["latest_draft"], tok)
    for old in dep.get("files", []):
        req("DELETE", f"{ZENODO}/deposit/depositions/{dep['id']}/files/{old['id']}", tok)

    sums: dict[str, str] = {}
    for f in bundle():
        payload = f.read_bytes()
        req("PUT", f"{dep['links']['bucket']}/{f.name}", tok, payload, "application/octet-stream")
        sums[f.name] = hashlib.sha256(payload).hexdigest()

    meta = {"metadata": {
        "title": TITLE,
        "upload_type": "publication",
        "publication_type": "preprint",
        "description": DESCRIPTION,
        "creators": [{"name": "Callens, Xavier", "affiliation": "SocrateAI Lab"}],
        "keywords": KEYWORDS,
        "license": "cc-by-4.0",
        "version": "3",
        "related_identifiers": [
            {"identifier": "https://github.com/xaviercallens/SocrateAI-Scientific-Genetic",
             "relation": "isSupplementTo", "scheme": "url"},
            {"identifier": f"10.5281/zenodo.{PREVIOUS}", "relation": "isNewVersionOf", "scheme": "doi"},
        ],
    }}
    req("PUT", f"{ZENODO}/deposit/depositions/{dep['id']}", tok, json.dumps(meta).encode())
    done = req("POST", f"{ZENODO}/deposit/depositions/{dep['id']}/actions/publish", tok)

    receipt = {
        "previous": PREVIOUS,
        "id": done.get("id"),
        "doi": done.get("doi"),
        "conceptrecid": done.get("conceptrecid"),
        "state": done.get("state"),
        "files": sorted(sums),
        "sha256": sums,
    }
    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(receipt["id"], receipt["doi"], receipt["state"], f"({len(sums)} files)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""
Harness for `SocrateAI-Scientific-Agora-LeanFlow` (`leanflow_mcp`),
`SocrateAI-Scientific-Mensura`, and `SocrateAI-Scientific-Elenchus` across the 3 Breakthrough Domains.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

LEANFLOW_ROOT = Path("/home/xavkal/xdev/SocrateAI-Scientific-Agora-LeanFlow")
if str(LEANFLOW_ROOT) not in sys.path:
    sys.path.insert(0, str(LEANFLOW_ROOT))

from leanflow_mcp import (  # noqa: E402
    invariants_check,
    lean_certify,
    ode_solve,
    ode_solve_frontier,
)


def run_leanflow_3domains_benchmark() -> dict[str, Any]:
    """
    Runs the LeanFlow MCP solver and formal certification suite across:
    - Stiff Robertson ODE (`scipy_radau` + `rusty_cvode` analytic Jacobian)
    - Domain 3A (`3A_fusion_stellarator` + Lean 4.34.1 certificate)
    - Domain 3B (`3B_navier_stokes_mensura_cap` + Lean 4.34.1 certificate)
    - Domain 3C (`3C_flt_pqc_quantum` + Lean 4.34.1 certificate)
    - Fail-closed negative controls (`sorry` rejection, false Euler refutation)
    """
    robertson_cvode = ode_solve("robertson", t_span=(0.0, 4.0e10), solver="rusty_cvode", use_analytic_jac=True)
    inv_robertson = invariants_check(
        robertson_cvode["results"]["times"],
        robertson_cvode["results"]["states"],
        [
            {"type": "linear_sum", "weights": [1.0, 1.0, 1.0], "target": 1.0, "tol": 1e-8, "name": "mass"},
            {"type": "non_negative", "tol": -1e-12, "name": "positivity"},
        ],
    )

    dom_3a = ode_solve_frontier("3A_fusion_stellarator")
    cert_3a = lean_certify({"type": "domain_3a_littlejohn"})

    dom_3b = ode_solve_frontier("3B_navier_stokes_mensura_cap")
    cert_3b = lean_certify({"type": "domain_3b_radii_poly_mensura"})

    dom_3c = ode_solve_frontier("3C_flt_pqc_quantum")
    cert_3c = lean_certify({"type": "domain_3c_flt_pqc"})

    dom_3e = ode_solve_frontier("3E_genetics_dna_replication_art")
    cert_3e = lean_certify({"type": "domain_3e_genetics_dna_art"})
    xgudhi_betti_cert = lean_certify({"type": "betti_euler", "betti": [1, 6, 0], "expected_chi": -5})

    sorry_reject = lean_certify({"type": "raw_lean_snippet", "source": "theorem bad : 1 = 2 := by sorry"})
    false_refute = lean_certify({"type": "betti_euler", "betti": [1, 2, 1], "expected_chi": 2})

    gate_passed = bool(
        robertson_cvode["status"] == "completed"
        and inv_robertson["results"]["all_satisfied"]
        and dom_3a["status"] == "verified"
        and cert_3a["status"] == "verified"
        and dom_3b["status"] == "verified"
        and cert_3b["status"] == "verified"
        and dom_3c["status"] == "verified"
        and cert_3c["status"] == "verified"
        and dom_3e["status"] == "verified"
        and cert_3e["status"] == "verified"
        and xgudhi_betti_cert["status"] == "verified"
        and sorry_reject["status"] == "failed"
        and sorry_reject["results"]["reason"] == "contains_sorry"
        and false_refute["status"] == "refuted"
    )

    return {
        "subproject": "socrateai_scientific_agora_leanflow",
        "robertson_cvode_steps": robertson_cvode["results"]["n_steps"],
        "robertson_cvode_rhs_evals": robertson_cvode["results"]["n_rhs_evals"],
        "robertson_invariants_satisfied": inv_robertson["results"]["all_satisfied"],
        "domain_3a": dom_3a["results"],
        "domain_3a_lean_tier": cert_3a["tier"],
        "domain_3b": dom_3b["results"],
        "domain_3b_lean_tier": cert_3b["tier"],
        "domain_3c_deuring_norm": dom_3c["results"]["flt_modular_pqc"]["deuring_response_norm"],
        "domain_3c_hellegouarch_c4": dom_3c["results"]["pqc_isogeny"]["c4"],
        "domain_3c_lean_tier": cert_3c["tier"],
        "domain_3e": dom_3e["results"],
        "domain_3e_lean_tier": cert_3e["tier"],
        "xgudhi_hecke_hairpin_betti_certified": xgudhi_betti_cert["status"] == "verified",
        "fail_closed_sorry_rejected": sorry_reject["results"]["reason"] == "contains_sorry",
        "fail_closed_false_refuted": false_refute["status"] == "refuted",
        "gate_passed": gate_passed,
    }


def get_subproject_manifest() -> dict[str, Any]:
    """Returns the certified subproject manifest for SocrateAI-Scientific-Agora-LeanFlow."""
    bench = run_leanflow_3domains_benchmark()
    return {
        "subproject": {
            "slug": "socrateai_scientific_agora_leanflow",
            "frontier_id": "LEANFLOW",
            "github_url": "https://github.com/xaviercallens/SocrateAI-Scientific-Agora-LeanFlow",
            "robertson_cvode_steps": bench["robertson_cvode_steps"],
            "robertson_cvode_rhs_evals": bench["robertson_cvode_rhs_evals"],
            "robertson_invariants_satisfied": bench["robertson_invariants_satisfied"],
            "domain_3a_lean_tier": bench["domain_3a_lean_tier"],
            "domain_3b_lean_tier": bench["domain_3b_lean_tier"],
            "domain_3c_lean_tier": bench["domain_3c_lean_tier"],
            "domain_3c_deuring_norm": bench["domain_3c_deuring_norm"],
            "fail_closed_sorry_rejected": bench["fail_closed_sorry_rejected"],
            "fail_closed_false_refuted": bench["fail_closed_false_refuted"],
            "gate_passed": bench["gate_passed"],
        },
        "benchmark": bench,
    }


if __name__ == "__main__":
    res = get_subproject_manifest()
    print(json.dumps(res, indent=2))
    if not res["subproject"]["gate_passed"]:
        sys.exit(1)

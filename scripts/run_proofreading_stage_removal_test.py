"""Test the proofreading-stage-removal model class against published Pol-epsilon data.

The model class (Paper E, GAR.H15, and multi-stage Hopfield-Ninio models generally) says a
mutator allele acts by removing proofreading stages. Theorem PSR_2 of
``formal/ANSE/ProofreadingStageRemoval.lean`` proves that such an allele can never be more
error-prone than the stage-free (exonuclease-null) enzyme, because each stage contributes an
escape factor in [0, 1]. That is a falsifiable ceiling, and this script tests it.

Decision rule: a variant breaks the ceiling when the lower bound of its 95% confidence
interval exceeds the upper bound of the exonuclease-null interval, both expressed as fold
change over the wild type of the same experiment. Confidence intervals are the published
ones; no distribution is assumed and nothing is fitted.

Run: .venv/bin/python scripts/run_proofreading_stage_removal_test.py
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data" / "pole_proofreading_literature.json"
OUT = HERE / "data" / "proofreading_stage_removal_test.json"


def fold(rate: float, ci: list[float], wt_rate: float, wt_ci: list[float]) -> dict[str, float]:
    """Fold change over wild type, with a conservative interval (worst case over both intervals)."""
    return {
        "fold": rate / wt_rate,
        "fold_low": ci[0] / wt_ci[1],
        "fold_high": ci[1] / wt_ci[0],
    }


def classify(variant: dict[str, float], ceiling: dict[str, float]) -> str:
    """Compare a variant's fold change with the exonuclease-null ceiling."""
    if variant["fold_low"] > ceiling["fold_high"]:
        return "breaks_ceiling"
    if variant["fold_high"] < ceiling["fold_low"]:
        return "below_ceiling"
    return "indistinguishable"


def run() -> dict[str, Any]:
    lit = json.loads(DATA.read_text(encoding="utf-8"))
    kane = lit["in_vivo_can1"]["kane2014"]["strains"]
    barbari = lit["in_vivo_can1"]["barbari2018"]["strains"]

    k_wt = kane["POL2_wild_type"]
    ceiling = fold(
        kane["pol2-4_exonuclease_null"]["rate"], kane["pol2-4_exonuclease_null"]["ci"],
        k_wt["rate"], k_wt["ci"],
    )

    # --- Primary test: both alleles measured side by side in Kane 2014 Table 2 ---
    p301r = fold(kane["pol2-P301R"]["rate"], kane["pol2-P301R"]["ci"], k_wt["rate"], k_wt["ci"])
    primary = {
        "comparison": "pol2-P301R vs pol2-4, same table, same isogenic MMR-proficient background",
        "exonuclease_null_ceiling_fold": ceiling,
        "P301R_fold": p301r,
        "ratio_P301R_over_exonuclease_null": kane["pol2-P301R"]["rate"] / kane["pol2-4_exonuclease_null"]["rate"],
        "ratio_lower_bound": kane["pol2-P301R"]["ci"][0] / kane["pol2-4_exonuclease_null"]["ci"][1],
        "verdict": classify(p301r, ceiling),
        "PSR_3_implied_insertion_error_multiple": kane["pol2-P301R"]["rate"] / kane["pol2-4_exonuclease_null"]["rate"],
    }

    # --- Secondary: the 13-variant panel, normalised to its own wild type (cross-paper) ---
    b_wt = barbari["POL2_wild_type"]
    panel = {}
    for name, rec in barbari.items():
        if name == "POL2_wild_type":
            continue
        f = fold(rec["rate"], rec["ci"], b_wt["rate"], b_wt["ci"])
        panel[name] = {
            "human_equivalent": rec.get("human_equivalent"),
            **f,
            "verdict": classify(f, ceiling),
        }

    # --- Controls ---
    exo_self = classify(ceiling, ceiling)
    controls = {
        "positive_control_P301R_must_break": {
            "expected": "breaks_ceiling",
            "observed": primary["verdict"],
            "passed": primary["verdict"] == "breaks_ceiling",
            "rationale": "Kane 2014 states the effect is ~50-fold above pol2-4; a test that cannot detect this is not sensitive enough to run.",
        },
        "negative_control_exonuclease_null_vs_itself": {
            "expected": "indistinguishable",
            "observed": exo_self,
            "passed": exo_self == "indistinguishable",
            "rationale": "The stage-free reference is the ceiling, so it must never be reported as breaking it.",
        },
        "negative_control_neutral_variants_must_not_break": {
            "expected": "no variant among the functionally neutral alleles breaks the ceiling",
            "variants": {n: panel[n]["verdict"] for n in
                         ("pol2-F139L", "pol2-R252H", "pol2-N351S", "pol2-R778W", "pol2-A979V", "pol2-D1757N")},
            "passed": all(panel[n]["verdict"] != "breaks_ceiling" for n in
                          ("pol2-F139L", "pol2-R252H", "pol2-N351S", "pol2-R778W", "pol2-A979V", "pol2-D1757N")),
            "rationale": "These alleles are reported as not significantly mutator; flagging any of them would show the rule fires spuriously.",
        },
        "negative_control_V411L_must_not_break": {
            "expected": "below_ceiling or indistinguishable",
            "observed": panel["pol2-V426L"]["verdict"],
            "passed": panel["pol2-V426L"]["verdict"] != "breaks_ceiling",
            "rationale": "Barbari: pol2-V426L 'showed no significant increase in Can-R mutation'. A stage-removal account is adequate for it.",
        },
    }

    # --- In vitro versus in vivo discordance ---
    iv = lit["in_vitro_fidelity"]
    in_vitro_ratio = iv["base_substitution_error_rate"]["P301R"] / iv["base_substitution_error_rate"]["exonuclease_null"]
    in_vitro_ratio_mf = (
        (iv["mutant_frequency"]["P301R"] - iv["mutant_frequency"]["background"])
        / (iv["mutant_frequency"]["exonuclease_null"] - iv["mutant_frequency"]["background"])
    )
    discordance = {
        "in_vitro_ratio_P301R_over_exonuclease_null_base_substitutions": in_vitro_ratio,
        "in_vitro_ratio_background_corrected_mutant_frequency": in_vitro_ratio_mf,
        "in_vivo_ratio_P301R_over_exonuclease_null": primary["ratio_P301R_over_exonuclease_null"],
        "discordance_factor": primary["ratio_P301R_over_exonuclease_null"] / in_vitro_ratio,
        "meaning": (
            "In the purified system, with no mismatch repair, P301R is more accurate than the "
            "exonuclease-null enzyme; in cells, with mismatch repair present, it is far worse. "
            "The ordering reverses, so the difference cannot lie in the polymerase's intrinsic "
            "error production alone: it requires a term describing how errors interact with "
            "mismatch repair. The factor below bounds the size of that term."
        ),
        "caveat": (
            "Different systems: per-nucleotide lacZ error rates with purified enzyme versus "
            "per-locus CAN1 rates per generation in vivo, where both Pol delta and Pol epsilon "
            "replicate. The ratio of ratios is more defensible than any absolute value, and the "
            "direction of the reversal does not depend on the normalisation."
        ),
        "constraint_on_mmr_escape": (
            "Writing the in vivo rate as (insertion error) x (proofreading escape) x (mismatch "
            "repair escape), the first two factors are what the purified assay measures. The "
            "discordance therefore bounds the ratio of mismatch-repair escape between the two "
            "alleles. Because escape cannot exceed 1, the exonuclease-null allele's escape is at "
            "most the reciprocal of the discordance factor."
        ),
        "counter_evidence": lit["mmr_interaction"]["finding"] + " " + lit["mmr_interaction"]["implication"],
    }
    discordance["max_exonuclease_null_mmr_escape"] = 1.0 / discordance["discordance_factor"]

    breaking = [n for n, v in panel.items() if v["verdict"] == "breaks_ceiling"]
    return {
        "question": (
            "Can a mutator allele of DNA polymerase epsilon be explained as the loss of "
            "proofreading stages?"
        ),
        "theorem_tested": "ANSE.ProofreadingStageRemoval.PSR_2_bounded_by_exo_null (formal/ANSE/ProofreadingStageRemoval.lean)",
        "falsification_criterion": (
            "The model class is refuted if any allele's mutation rate significantly exceeds the "
            "exonuclease-null rate, since PSR_2 proves that is impossible for any allele that "
            "only removes stages, whatever the number of stages or the per-stage factors."
        ),
        "primary_result": primary,
        "variant_panel_cross_paper": panel,
        "variants_breaking_ceiling": breaking,
        "controls": controls,
        "in_vitro_in_vivo_discordance": discordance,
        "conclusion": (
            "Refuted. The P301R allele, which models human POLE P286R, exceeds the "
            "exonuclease-null ceiling by a factor of "
            f"{primary['ratio_P301R_over_exonuclease_null']:.1f} (at least "
            f"{primary['ratio_lower_bound']:.1f} on published confidence intervals), and "
            f"{len(breaking)} further alleles break the same ceiling. By PSR_3 the insertion or "
            "extension step itself must be at least that much more error-prone, so no choice of "
            "stage count or per-stage factor can rescue the model class. Paper E's own GAR.H15 "
            "treatment of P286R as the loss of two of three stages is therefore refuted as a "
            "mechanism, not merely mis-parameterised."
        ),
    }


if __name__ == "__main__":
    result = run()
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    ok = all(c["passed"] for c in result["controls"].values())
    print(json.dumps(result, indent=2))
    print("\nCONTROLS PASSED:", ok)
    raise SystemExit(0 if ok else 1)

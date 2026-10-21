"""Compare the Hopfield-Ninio POLE fold-change predictions of Paper E with published
TCGA endometrial exome mutation rates (Kandoth et al., Nature 497:67-73, 2013).

Observed rates are copied from the paper's Results text (read 2026-10-21 at
https://pmc.ncbi.nlm.nih.gov/articles/PMC3704730/); the source prints the unit as
"x 10^-6 mutations per Mb", which is read here as mutations per Mb (the only reading
consistent with the >100 mut/Mb "ultramutated" convention). Nothing is fitted.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "data" / "tcga_pole_empirical_check.json"

OBSERVED_MUT_PER_MB: dict[str, float] = {
    "POLE_ultramutated": 232.0,
    "MSI_hypermutated": 18.0,
    "copy_number_low_mostly_MSS_endometrioid": 2.9,
    "copy_number_high_serous_like": 2.3,
}
GAMMA = Fraction(3, 250)


def model_fold(m_wild_type: int, m_mutant: int) -> Fraction:
    """Fold increase in error fraction when proofreading stages drop m_wild_type -> m_mutant."""
    return GAMMA ** m_mutant / GAMMA ** m_wild_type


def run() -> dict[str, object]:
    pole = OBSERVED_MUT_PER_MB["POLE_ultramutated"]
    baselines = {k: v for k, v in OBSERVED_MUT_PER_MB.items() if k != "POLE_ultramutated"}
    observed_fold = {k: pole / v for k, v in baselines.items()}
    predictions = {
        "m3_to_m1_as_in_paper_E_GAR_H15": float(model_fold(3, 1)),
        "m3_to_m2_single_stage_loss": float(model_fold(3, 2)),
        "m3_to_m3_negative_control_no_loss": float(model_fold(3, 3)),
    }
    ref = observed_fold["copy_number_low_mostly_MSS_endometrioid"]
    return {
        "source": "Kandoth et al. 2013 Nature 497:67-73, exome analysis section; 17 POLE-ultramutated tumours of 248 sequenced; "
                  "hotspots P286R/V411L in 13 of 17",
        "observed_mut_per_Mb": OBSERVED_MUT_PER_MB,
        "observed_fold_POLE_over": observed_fold,
        "model_predicted_fold": predictions,
        "paper_E_absolute_TMB_P286R_mut_per_Mb": 12500 / 9,
        "absolute_overshoot_vs_observed_POLE": (12500 / 9) / pole,
        "fold_overshoot_m3_to_m1_vs_CNlow": predictions["m3_to_m1_as_in_paper_E_GAR_H15"] / ref,
        "fold_ratio_m3_to_m2_vs_CNlow": predictions["m3_to_m2_single_stage_loss"] / ref,
        "interpretation": "The m=3->1 prediction (6944x) overshoots the observed ~80x by ~87x: refuted as a quantitative claim. "
                          "A single-stage loss (1/gamma = 83.3x) is within 4% of the observed ratio; this is a post-hoc observation "
                          "on one data point, not a confirmation.",
    }


if __name__ == "__main__":
    result = run()
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))

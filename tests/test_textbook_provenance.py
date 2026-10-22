"""The textbook's quoted rates must come from the provenance record.

The first edition of Section 4 carried four hard-coded variant ratios (21, 12.5,
9.1, 0.55) that could not be reproduced from
`data/pole_proofreading_literature.json` by any normalisation, and a wild-type
rate that disagreed with it by a factor of ten. These tests pin the numbers to
the record so that a disagreement fails here rather than reaching a reader.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[1]
LITERATURE = ROOT / "data" / "pole_proofreading_literature.json"
VALUES = ROOT / "textbook" / "data" / "figure_values.json"


def _load_figures_module() -> Any:
    """Import make_figures.py by path; it is a script, not an installed module."""
    spec = importlib.util.spec_from_file_location(
        "make_figures", ROOT / "scripts" / "textbook" / "make_figures.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def mf() -> Any:
    return _load_figures_module()


@pytest.fixture(scope="module")
def literature() -> dict[str, Any]:
    return json.loads(LITERATURE.read_text())


def test_kane_rates_match_the_provenance_record(mf: Any,
                                                literature: dict[str, Any]) -> None:
    """The three rates carrying Theorem 4.1 are read, not retyped."""
    strains = literature["in_vivo_can1"]["kane2014"]["strains"]
    assert mf.KANE_P301R == strains["pol2-P301R"]["rate"]
    assert mf.KANE_EXO_NULL == strains["pol2-4_exonuclease_null"]["rate"]
    assert mf.KANE_WT == strains["POL2_wild_type"]["rate"]
    # The decimal slip that the internal review caught: 1.5 would be wrong.
    assert mf.KANE_WT == 15.0


def test_wild_type_reproduces_the_sources_own_fold_statement(mf: Any) -> None:
    """Kane & Shcherbakova report pol2-P301R as ~150-fold over wild type.

    That statement is the independent check on the wild-type rate: 2300/15 is
    153, while 2300/1.5 would be 1533.
    """
    fold = mf.KANE_P301R / mf.KANE_WT
    assert 145.0 < fold < 160.0, f"P301R/WT = {fold}, not the published ~150-fold"


def test_rho_is_past_the_ceiling_even_on_the_confidence_intervals(mf: Any) -> None:
    """PSR_2 bounds rho by 1; the worst case on the published CIs must still fail."""
    rho = mf.KANE_P301R / mf.KANE_EXO_NULL
    worst = mf.KANE_P301R_CI[0] / mf.KANE_EXO_NULL_CI[1]
    assert rho > 1.0 and worst > 1.0
    assert 52.0 < rho < 52.5
    assert 35.0 < worst < 35.5


def test_barbari_rho_is_derived_and_not_hard_coded(mf: Any,
                                                   literature: dict[str, Any]) -> None:
    """Each panel ratio equals the cross-paper normalisation, recomputed here."""
    strains = literature["in_vivo_can1"]["barbari2018"]["strains"]
    wt = strains["POL2_wild_type"]["rate"]
    exo_null_fold = mf.KANE_EXO_NULL / mf.KANE_WT
    for allele in ("pol2-S474F", "pol2-F382S", "pol2-P301H", "pol2-D290V",
                   "pol2-V426L"):
        expected = (strains[allele]["rate"] / wt) / exo_null_fold
        assert mf.barbari_rho(allele) == pytest.approx(expected)
    # The retired literals must not reappear.
    retired = {21.0, 12.5, 9.1, 0.55}
    computed = {round(mf.barbari_rho(a), 2)
                for a in ("pol2-S474F", "pol2-F382S", "pol2-P301H", "pol2-V426L")}
    assert not (computed & retired)


def test_the_panel_splits_cancer_alleles_from_the_rest(mf: Any) -> None:
    """Control 1 needs alleles above the ceiling; control 2 needs alleles below."""
    above = [a for a in ("pol2-S474F", "pol2-F382S", "pol2-P301H")
             if mf.barbari_rho(a) > 1.0]
    below = [a for a in ("pol2-D290V", "pol2-V426L")
             if mf.barbari_rho(a) < 1.0]
    assert len(above) == 3, "control 1 fails: no further allele breaks the ceiling"
    assert len(below) == 2, "control 2 fails: every allele is flagged"


def test_recorded_figure_values_agree_with_the_record(mf: Any) -> None:
    """figure_values.json is what the prose quotes; it must match the module."""
    sec04 = json.loads(VALUES.read_text())["sec04"]
    assert sec04["kane_wild_type_e8"] == mf.KANE_WT
    assert sec04["rho_ci_floor"] == pytest.approx(
        mf.KANE_P301R_CI[0] / mf.KANE_EXO_NULL_CI[1])
    by_allele = {row["allele"]: row for row in sec04["panel"]}
    for allele in ("pol2-S474F", "pol2-F382S", "pol2-P301H", "pol2-V426L"):
        row = by_allele[allele]
        assert row["rho"] == pytest.approx(mf.barbari_rho(allele))
        # "Breaks the ceiling" is decided by the number, never asserted beside it.
        assert row["breaks_ceiling"] is (row["rho"] > 1.0)


def test_luria_delbruck_variance_is_not_reported_as_a_stable_statistic() -> None:
    """Section 6 must carry the reseed spread, not one draw.

    The sample variance-to-mean ratio is biased far below the closed form, so
    recording only the single-seed value is the defect the review removed.
    """
    sec06 = json.loads(VALUES.read_text())["sec06"]
    assert sec06["reseeds"] >= 10
    lo, hi = sec06["reseed_variance_over_mean_min"], sec06["reseed_variance_over_mean_max"]
    assert hi > 10.0 * lo, "the point of the figure is that the estimator is unstable"
    # Every reseed should fall far short of the closed form N/g.
    theory = sec06["theory_variance_over_mean"]
    assert theory == pytest.approx(
        sec06["final_population"] / sec06["generations"])
    assert hi < theory / 10.0, "sample variance should underestimate N/g by orders"

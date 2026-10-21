"""Regression tests for the proofreading-stage-removal falsification (Paper E, edition 3).

These lock three things that a future edit could silently break: the published literature
values the result rests on, the ceiling test's controls, and the refutation itself.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from run_proofreading_stage_removal_test import classify, fold, run  # noqa: E402


@pytest.fixture(scope="module")
def result() -> dict:
    return run()


def test_published_literature_values_are_unaltered() -> None:
    """The refutation is only as good as the numbers it quotes from the primary sources."""
    lit = json.loads((REPO_ROOT / "data" / "pole_proofreading_literature.json").read_text())
    kane = lit["in_vivo_can1"]["kane2014"]["strains"]
    assert kane["POL2_wild_type"]["rate"] == 15
    assert kane["pol2-P301R"]["rate"] == 2300
    assert kane["pol2-P301R"]["ci"] == [1800, 2800]
    # The stage-free reference must be the double catalytic-carboxylate substitution.
    assert kane["pol2-4_exonuclease_null"]["rate"] == 44
    assert "D290A" in kane["pol2-4_exonuclease_null"]["genotype"]
    # Xing's in vitro ordering is the opposite of the in vivo one; both must survive edits.
    iv = lit["in_vitro_fidelity"]["base_substitution_error_rate"]
    assert iv["P301R"] < iv["exonuclease_null"]
    assert (iv["exonuclease_null"], iv["P301R"]) == (155, 59)


def test_ceiling_is_broken_by_p301r_with_non_overlapping_intervals(result: dict) -> None:
    """The primary claim: P301R exceeds the exonuclease-null ceiling, interval-disjointly."""
    primary = result["primary_result"]
    assert primary["verdict"] == "breaks_ceiling"
    assert primary["ratio_P301R_over_exonuclease_null"] == pytest.approx(2300 / 44, rel=1e-9)
    # Even the most conservative reading of the published intervals keeps the excess large.
    assert primary["ratio_lower_bound"] > 30.0
    assert primary["P301R_fold"]["fold_low"] > primary["exonuclease_null_ceiling_fold"]["fold_high"]


def test_controls_pass_and_are_capable_of_failing(result: dict) -> None:
    """A control that cannot fail proves nothing, so check both the verdicts and the rule."""
    assert all(c["passed"] for c in result["controls"].values())
    ceiling = result["primary_result"]["exonuclease_null_ceiling_fold"]
    # The classifier must not flag the reference against itself...
    assert classify(ceiling, ceiling) == "indistinguishable"
    # ...but must flag something an order of magnitude above it.
    far_above = {"fold": 100.0, "fold_low": 90.0, "fold_high": 110.0}
    assert classify(far_above, ceiling) == "breaks_ceiling"
    # ...and must place a clearly smaller rate below it.
    far_below = {"fold": 0.2, "fold_low": 0.15, "fold_high": 0.25}
    assert classify(far_below, ceiling) == "below_ceiling"


def test_v411l_and_neutral_alleles_are_not_flagged(result: dict) -> None:
    """Negative controls: a non-mutator allele must stay consistent with stage removal."""
    panel = result["variant_panel_cross_paper"]
    assert panel["pol2-V426L"]["verdict"] != "breaks_ceiling"
    for neutral in ("pol2-F139L", "pol2-R252H", "pol2-N351S", "pol2-R778W", "pol2-D1757N"):
        assert panel[neutral]["verdict"] != "breaks_ceiling", neutral
    # The alleles that do break it are the known strong mutators, not an arbitrary set.
    assert set(result["variants_breaking_ceiling"]) == {"pol2-P301H", "pol2-F382S", "pol2-S474F"}


def test_in_vitro_in_vivo_ordering_actually_reverses(result: dict) -> None:
    """The discordance is the paper's mechanistic claim; its sign must not silently flip."""
    d = result["in_vitro_in_vivo_discordance"]
    assert d["in_vitro_ratio_P301R_over_exonuclease_null_base_substitutions"] < 1.0
    assert d["in_vivo_ratio_P301R_over_exonuclease_null"] > 1.0
    assert d["discordance_factor"] > 100.0
    # Escape fractions are probabilities, so the implied bound must be a usable one.
    assert 0.0 < d["max_exonuclease_null_mmr_escape"] < 1.0


def test_fold_helper_uses_conservative_interval_arithmetic() -> None:
    """fold() must widen, never narrow, the published uncertainty."""
    f = fold(100.0, [80.0, 120.0], 10.0, [8.0, 12.0])
    assert f["fold"] == pytest.approx(10.0)
    assert f["fold_low"] == pytest.approx(80.0 / 12.0)
    assert f["fold_high"] == pytest.approx(120.0 / 8.0)
    assert f["fold_low"] < f["fold"] < f["fold_high"]


def test_psr2_ceiling_holds_for_every_stage_assignment() -> None:
    """The arithmetic content of PSR_2: no product of factors in [0,1] exceeds the empty product.

    The Lean proof in formal/ANSE/ProofreadingStageRemoval.lean is the real guarantee; this is
    an independent exact-rational spot check that the statement being tested is the right one.
    """
    eps = Fraction(1)
    for q in (Fraction(1, 100), Fraction(3, 250), Fraction(7, 10), Fraction(1)):
        exo_null = eps  # empty product
        for kept in range(4):
            assert eps * q**kept <= exo_null
        # Removing a stage never lowers the error rate (PSR_1).
        for kept in range(3):
            assert eps * q ** (kept + 1) <= eps * q**kept

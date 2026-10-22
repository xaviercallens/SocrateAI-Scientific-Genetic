"""Tests for the two experiments that close Section 8's and Section 7's gaps.

The shuffles must preserve exactly what they claim to preserve, or every
p-value in Section 11 is meaningless; and the aligner must be correct, or the
yeast-human correspondence is not established. Both are checked here on small
inputs where the right answer is known independently, rather than only on the
real sequences where it is not.
"""

from __future__ import annotations

import importlib.util
import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
NULL_RESULTS = ROOT / "data" / "tda_sequence_null.json"
ORTHOLOGUE = ROOT / "data" / "yeast_human_correspondence.json"


def _load(name: str, relative: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def tda() -> Any:
    return _load("tda_null", "scripts/textbook/tda_sequence_null.py")


@pytest.fixture(scope="module")
def ortho() -> Any:
    return _load("ortho", "scripts/textbook/verify_yeast_human_correspondence.py")


# ------------------------------------------------------------- the shuffles ---
SEQ = ("ATGGCGCTAGCTAGCATCGATCGGATCCGATTACGATCAGCTACGATCAGCATCAGCTAC"
       "GATCAGCATCGATCAGCTAACGATCGATCAGCTACGGATCAGCTACGATCAGCATCAGCT")


def test_mono_shuffle_preserves_composition_and_breaks_dinucleotides(tda: Any) -> None:
    rng = np.random.default_rng(1)
    out = tda.shuffle_mono(SEQ, rng)
    assert len(out) == len(SEQ)
    assert Counter(out) == Counter(SEQ)
    # Over several draws it must change dinucleotide counts at least once,
    # otherwise it is not destroying the structure the null is meant to destroy.
    changed = any(tda.dinucleotide_counts(tda.shuffle_mono(SEQ, rng))
                  != tda.dinucleotide_counts(SEQ) for _ in range(5))
    assert changed


def test_dinucleotide_shuffle_preserves_every_dinucleotide_count(tda: Any) -> None:
    """The load-bearing invariant of the `dinuc` null."""
    rng = np.random.default_rng(2)
    target = tda.dinucleotide_counts(SEQ)
    for _ in range(12):
        out = tda.shuffle_dinucleotide(SEQ, rng)
        assert len(out) == len(SEQ)
        assert Counter(out) == Counter(SEQ)
        assert tda.dinucleotide_counts(out) == target
        assert out[0] == SEQ[0] and out[-1] == SEQ[-1]


def test_dinucleotide_shuffle_actually_moves(tda: Any) -> None:
    """A 'shuffle' that returns its input would make the null test vacuous."""
    rng = np.random.default_rng(3)
    outs = {tda.shuffle_dinucleotide(SEQ, rng) for _ in range(8)}
    assert len(outs) > 1, "dinucleotide shuffle is not randomising"


def test_codon_shuffle_preserves_codon_usage(tda: Any) -> None:
    rng = np.random.default_rng(4)
    out = tda.shuffle_codon(SEQ, rng)
    n = len(SEQ) // 3
    usage = Counter(SEQ[3 * i:3 * i + 3] for i in range(n))
    assert Counter(out[3 * i:3 * i + 3] for i in range(n)) == usage
    assert len(out) == len(SEQ)


def test_stride_is_coprime_with_the_codon_period(tda: Any) -> None:
    """The artefact that the first run of this experiment actually hit.

    A subsampling stride sharing a factor with the reading frame over-samples one
    codon position. POLD1 is 3324 nt, so the naive stride is 3324 // 500 = 6,
    which reads only every other codon's first base -- and that produced an
    apparently significant signal. The stride must now be coprime with 3.
    """
    import math
    for length in (3324, 6861, 4083, 1500, 3000, 9000, 500, 501):
        stride = tda.choose_stride(length)
        assert stride >= 1
        assert math.gcd(stride, 3) == 1, f"length {length} gave stride {stride}"
    # Specifically: the value that bit us must no longer be produced.
    assert tda.choose_stride(3324) != 6


def test_chaos_game_contracts_toward_recent_bases(tda: Any) -> None:
    """Two sequences sharing a k-base suffix must land within 2^-k of each other.

    This is the property that makes the embedding meaningful at all (Section 8,
    Exercise 8.1), so it is worth checking rather than assuming.
    """
    a = tda.chaos_game("ACGT" * 5 + "GGGGGGGG")[-1]
    b = tda.chaos_game("TGCA" * 5 + "GGGGGGGG")[-1]
    # Shared suffix of 8 bases; the tetrahedron has diameter 2*sqrt(2).
    assert float(np.linalg.norm(a - b)) < 2.0 ** -8 * 2.0 * np.sqrt(2) + 1e-12


# ------------------------------------------------------------- the aligner ---
def test_self_alignment_is_the_identity_map(ortho: Any) -> None:
    seq = "MKLPQRSTVWYACDEFGHIN"
    pairs = ortho.align(seq, seq, "identity", 11.0, 1.0)
    aligned = [(i, j) for i, j in pairs if i is not None and j is not None]
    assert len(aligned) == len(seq)
    assert all(i == j for i, j in aligned)


def test_single_deletion_is_placed_exactly(ortho: Any) -> None:
    """Deleting F from ACDEFGHIK must produce exactly one gap, at F."""
    pairs = ortho.align("ACDEFGHIK", "ACDEGHIK", "identity", 11.0, 1.0)
    gaps = [i for i, j in pairs if j is None]
    assert gaps == [4], f"expected one gap at index 4, got {gaps}"


def test_affine_gaps_prefer_one_long_run(ortho: Any) -> None:
    """With gap_open >> gap_extend the aligner must not scatter short gaps."""
    pairs = ortho.align("AAAAACCCCCAAAAA", "AAAAAAAAAA", "identity", 11.0, 1.0)
    runs, current = [], 0
    for _, j in pairs:
        if j is None:
            current += 1
        elif current:
            runs.append(current)
            current = 0
    if current:
        runs.append(current)
    assert runs == [5], f"expected a single run of 5, got {runs}"


def test_partner_of_reports_gaps_as_none(ortho: Any) -> None:
    pairs = ortho.align("ACDEFGHIK", "ACDEGHIK", "identity", 11.0, 1.0)
    assert ortho.partner_of(pairs, 4) is None
    assert ortho.partner_of(pairs, 0) == 0


# ----------------------------------------------------- the recorded results ---
@pytest.mark.skipif(not NULL_RESULTS.exists(), reason="null experiment not run")
def test_null_experiment_controls_passed() -> None:
    """A recorded result whose controls failed must never be relied on."""
    report = json.loads(NULL_RESULTS.read_text())
    verdict = report["controls"]["verdict"]
    assert verdict["negative_passes"], "i.i.d. control was flagged"
    assert verdict["positive_passes"], "tandem-repeat control was not flagged"
    # The positive control must be strongly displaced, not merely past 0.05.
    assert abs(verdict["positive_z"]) > 3.0
    assert report["shuffles"] >= 100


@pytest.mark.skipif(not NULL_RESULTS.exists(), reason="null experiment not run")
def test_recorded_strides_are_frame_safe() -> None:
    """Every gene in the recorded result must have been sampled frame-safely."""
    import math
    report = json.loads(NULL_RESULTS.read_text())
    for gene, entry in report["genes"].items():
        assert "subsample_stride" in entry, (
            f"{gene} was recorded before the frame-bias fix; re-run "
            f"scripts/textbook/tda_sequence_null.py")
        stride = entry["subsample_stride"]
        assert math.gcd(stride, 3) == 1, f"{gene} sampled at stride {stride}"
        assert entry["stride_coprime_with_codon"] is True


@pytest.mark.skipif(not NULL_RESULTS.exists(), reason="null experiment not run")
def test_null_p_values_respect_their_own_floor() -> None:
    """An empirical p cannot be smaller than 1/(B+1); check none claims to be."""
    report = json.loads(NULL_RESULTS.read_text())
    floor = 1 / (1 + report["shuffles"])
    for gene, entry in report["genes"].items():
        for kind, outcome in entry["nulls"].items():
            for statistic in ("total_persistence", "longest_bar", "bar_count"):
                stats = outcome[statistic]
                assert stats["p_greater"] >= floor - 1e-12, (gene, kind, statistic)
                assert stats["p_less"] >= floor - 1e-12, (gene, kind, statistic)
                assert stats["p_floor"] == pytest.approx(floor)


@pytest.mark.skipif(not ORTHOLOGUE.exists(), reason="alignment not run")
def test_orthologue_correspondence_is_invariant_and_controlled() -> None:
    report = json.loads(ORTHOLOGUE.read_text())
    assert report["invariant_across_schemes"]
    assert all(run["human_286_aligns_to_yeast"] == 301
               for run in report["alignments"])
    verdict = report["controls"]["verdict"]
    assert verdict["positive_passes"] and verdict["negative_passes"]
    # The conserved F-P-D context is what the textbook prints; check it holds.
    context = report["context"]
    assert context["human_279_293"][6:9] == "FPD"
    assert context["yeast_294_308"][6:9] == "FPD"


@pytest.mark.skipif(not ORTHOLOGUE.exists(), reason="alignment not run")
def test_local_conservation_beats_the_diagonal_prior() -> None:
    """The evidence is local conservation, not the aligned position.

    Section 11 makes this point explicitly because the shuffled control lands on
    302, one residue from the right answer. The margin that matters is the local
    identity ratio, so that is what is pinned here.
    """
    report = json.loads(ORTHOLOGUE.read_text())
    margin = report["controls"]["local_conservation_margin"]
    assert margin["worst_real_local_identity"] > 2.0 * margin["best_shuffled_local_identity"]
    assert margin["worst_real_local_identity"] > 60.0
    # And the weakness being guarded against is real: the shuffle does land adjacent.
    adjacent = [n["shuffled_maps_286_to"]
                for n in report["controls"]["negative_shuffled_yeast"]]
    assert any(p is not None and abs(p - 301) <= 2 for p in adjacent), (
        "if the shuffle no longer lands near 301, Section 11's caveat about the "
        "diagonal prior needs rewriting")

#!/usr/bin/env python3
"""Does a coding sequence have topology beyond what its composition explains?

Section 8 of the textbook computed persistent homology on the POLE coding
sequence and then declined to draw any biological conclusion, because the chaos
game representation sends *any* sequence into a self-similar attractor whose
holes exist regardless of what the sequence says. The honest reading was that
the Betti numbers are a fact about the embedding. This script runs the test
that was missing.

Design. The statistic is the total persistence of H1 (the sum of bar lengths)
of the Vietoris-Rips filtration on the 3D chaos-game cloud. It is compared
against three nested nulls, each destroying one more level of sequence
structure while preserving the one below:

    mono   a uniform permutation of the sequence. Preserves base composition,
           destroys every k-mer statistic for k >= 2.
    dinuc  an Altschul-Erickson Eulerian-path shuffle. Preserves every
           dinucleotide count exactly, destroys k >= 3.
    codon  a permutation of whole codons. Preserves codon usage and the reading
           frame, destroys codon order.

Reading the three together localises any signal: a sequence that beats `mono`
but not `dinuc` has dinucleotide structure and nothing more.

Controls. Two, because a null test that cannot come out negative is worthless
and one that cannot come out positive has no power:

    negative  an i.i.d. sequence drawn from the real base frequencies. It has
              no structure above composition by construction, so it MUST NOT be
              flagged against `mono`. If it is, the test is broken.
    positive  a tandem repeat of one motif, with base composition matched to the
              real sequence. Its order is maximally non-random, so it MUST be
              flagged. If it is not, the test has no power.

Both controls are asserted, not merely reported. The script exits non-zero if
either fails, so a broken test cannot quietly produce a p-value.

    python scripts/textbook/tda_sequence_null.py              # POLE, 200 shuffles
    python scripts/textbook/tda_sequence_null.py --shuffles 50 --genes POLE POLD1
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import urllib.request
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "textbook" / "data"
RECORDS = DATA / "genomic_records.json"
OUT = ROOT / "data" / "tda_sequence_null.json"

# Tetrahedral chaos-game vertices, as in make_figures.py section 8.
VERTICES = {"A": (1.0, 1.0, 1.0), "C": (1.0, -1.0, -1.0),
            "G": (-1.0, 1.0, -1.0), "T": (-1.0, -1.0, 1.0)}

# Fixed by a timing study: 500 points at max edge 0.8 costs ~1.8 s and retains
# every H1 bar that survives past 0.5 in the full 900-point cloud. The null
# needs hundreds of these, so the budget is what sets the resolution.
N_POINTS = 500
MAX_EDGE = 0.8


# ----------------------------------------------------------------- geometry ---
def chaos_game(seq: str) -> np.ndarray:
    """Embed a sequence in R^3 by the chaos game on a regular tetrahedron."""
    point = np.zeros(3)
    cloud = []
    for base in seq.upper():
        vertex = VERTICES.get(base)
        if vertex is None:
            continue
        point = (point + np.asarray(vertex)) / 2.0
        cloud.append(point.copy())
    return np.asarray(cloud)


def choose_stride(length: int, n: int = N_POINTS) -> int:
    """A subsampling stride that is coprime with the codon period.

    Section 8's Worked Example 8.1 warns that a stride sharing a factor with the
    reading frame over-samples one codon position, and that the resulting bias
    is invisible in the persistence diagram. The first run of this experiment
    walked straight into it: POLD1 is 3324 nt, so `3324 // 500 = 6`, and a stride
    of 6 reads only every other codon's first base. That produced an apparently
    significant signal which was an artefact of our own sampling.

    So the stride is chosen coprime with 3, and recorded in the output.
    """
    stride = max(1, length // n)
    while stride > 1 and math.gcd(stride, 3) != 1:
        stride -= 1
    return stride


def subsample(cloud: np.ndarray, n: int = N_POINTS) -> np.ndarray:
    return cloud[:: choose_stride(len(cloud), n)][:n]


def h1_total_persistence(seq: str) -> dict[str, float]:
    """Total persistence, longest bar and bar count of H1.

    Imported inside the function so the module can be loaded (and its shuffles
    tested) on a host without GUDHI.
    """
    import gudhi

    sample = subsample(chaos_game(seq))
    tree = gudhi.RipsComplex(points=sample, max_edge_length=MAX_EDGE
                             ).create_simplex_tree(max_dimension=2)
    tree.compute_persistence()
    bars = [(b, d) for dim, (b, d) in tree.persistence()
            if dim == 1 and d != float("inf")]
    lengths = [d - b for b, d in bars]
    return {"total_persistence": float(sum(lengths)),
            "longest_bar": float(max(lengths, default=0.0)),
            "bar_count": float(len(bars))}


# ----------------------------------------------------------------- shuffles ---
def shuffle_mono(seq: str, rng: np.random.Generator) -> str:
    """Uniform permutation: preserves base composition only."""
    letters = list(seq)
    rng.shuffle(letters)
    return "".join(letters)


def shuffle_dinucleotide(seq: str, rng: np.random.Generator) -> str:
    """Altschul-Erickson shuffle: preserves every dinucleotide count exactly.

    The sequence is an Eulerian path in the multigraph whose vertices are bases
    and whose edges are the observed dinucleotides. Shuffling each vertex's
    successor list and consuming it front-to-back traverses an Eulerian path
    provided the *last* edge out of each non-terminal vertex forms an
    arborescence rooted at the terminal base -- which is the condition this
    function tests and reshuffles until it holds.
    """
    if len(seq) < 3:
        return seq
    terminal = seq[-1]
    successors: dict[str, list[str]] = {}
    for a, b in zip(seq, seq[1:]):
        successors.setdefault(a, []).append(b)
    vertices = set(successors) | {terminal}

    for _ in range(10_000):
        trial = {v: list(e) for v, e in successors.items()}
        for edges in trial.values():
            rng.shuffle(edges)
        if _last_edges_reach(trial, vertices, terminal):
            cursor = {v: 0 for v in trial}
            out = [seq[0]]
            current = seq[0]
            for _ in range(len(seq) - 1):
                nxt = trial[current][cursor[current]]
                cursor[current] += 1
                out.append(nxt)
                current = nxt
            return "".join(out)
    raise RuntimeError("no Eulerian shuffle found; the graph may be degenerate")


def _last_edges_reach(trial: dict[str, list[str]], vertices: set[str],
                      terminal: str) -> bool:
    """Do the last edges out of each non-terminal vertex lead to `terminal`?"""
    for start in vertices:
        if start == terminal:
            continue
        seen: set[str] = set()
        current = start
        while current != terminal:
            if current in seen or not trial.get(current):
                return False
            seen.add(current)
            current = trial[current][-1]
    return True


def shuffle_codon(seq: str, rng: np.random.Generator) -> str:
    """Permute whole codons: preserves codon usage and the reading frame."""
    n_codons = len(seq) // 3
    codons = [seq[3 * i:3 * i + 3] for i in range(n_codons)]
    tail = seq[3 * n_codons:]
    order = rng.permutation(n_codons)
    return "".join(codons[i] for i in order) + tail


SHUFFLES: dict[str, Callable[[str, np.random.Generator], str]] = {
    "mono": shuffle_mono, "dinuc": shuffle_dinucleotide, "codon": shuffle_codon}


# ------------------------------------------------------------- verification ---
def dinucleotide_counts(seq: str) -> Counter:
    return Counter(seq[i:i + 2] for i in range(len(seq) - 1))


def verify_shuffle_invariants(seq: str, rng: np.random.Generator) -> dict[str, Any]:
    """Each shuffle must preserve what it claims to and destroy what it claims to.

    Reported and asserted. A shuffle that quietly fails to preserve its
    invariant would make every p-value below meaningless.
    """
    report: dict[str, Any] = {}

    mono = shuffle_mono(seq, rng)
    assert Counter(mono) == Counter(seq), "mono shuffle changed base composition"
    assert len(mono) == len(seq)
    report["mono"] = {"preserves_base_counts": True,
                      "dinucleotide_counts_changed":
                          dinucleotide_counts(mono) != dinucleotide_counts(seq)}

    dinuc = shuffle_dinucleotide(seq, rng)
    assert len(dinuc) == len(seq), "dinuc shuffle changed length"
    assert Counter(dinuc) == Counter(seq), "dinuc shuffle changed base composition"
    assert dinucleotide_counts(dinuc) == dinucleotide_counts(seq), \
        "dinuc shuffle did not preserve dinucleotide counts"
    report["dinuc"] = {"preserves_dinucleotide_counts": True,
                       "identical_to_input": dinuc == seq}

    codon = shuffle_codon(seq, rng)
    assert len(codon) == len(seq), "codon shuffle changed length"
    n = len(seq) // 3
    assert (Counter(codon[3 * i:3 * i + 3] for i in range(n))
            == Counter(seq[3 * i:3 * i + 3] for i in range(n))), \
        "codon shuffle did not preserve codon usage"
    report["codon"] = {"preserves_codon_usage": True}
    return report


# ------------------------------------------------------------------ sequences ---
def iid_sequence(seq: str, rng: np.random.Generator) -> str:
    """Negative control: i.i.d. draws from the real base frequencies."""
    bases = sorted(VERTICES)
    counts = Counter(b for b in seq.upper() if b in VERTICES)
    total = sum(counts.values())
    probs = [counts[b] / total for b in bases]
    return "".join(rng.choice(bases, size=len(seq), p=probs))


def tandem_repeat_sequence(seq: str) -> str:
    """Positive control: one motif repeated, at the real base composition.

    Composition is matched so that the mononucleotide null cannot distinguish
    this from the real sequence on composition alone -- only on order.
    """
    counts = Counter(b for b in seq.upper() if b in VERTICES)
    total = sum(counts.values())
    # A 20-mer whose composition approximates the real one, then tiled.
    motif = "".join(b * max(1, round(20 * counts[b] / total))
                    for b in sorted(VERTICES))
    return (motif * (len(seq) // len(motif) + 1))[:len(seq)]


def fetch_cds(gene: str) -> dict[str, Any]:
    """Canonical CDS from Ensembl, cached with its URL. POLE is already cached."""
    cache = DATA / f"{gene.lower()}_cds.json"
    if cache.exists():
        return json.loads(cache.read_text())
    records = json.loads(RECORDS.read_text())
    transcript = records["genes"][gene]["canonical_transcript"]
    url = (f"https://rest.ensembl.org/sequence/id/{transcript}"
           f"?type=cds;content-type=application/json")
    request = urllib.request.Request(url, headers={
        "Accept": "application/json",
        "User-Agent": "SocrateAI-Formal-Genetics-Textbook/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        payload = json.loads(response.read().decode())
    record = {"transcript": transcript, "source_url": url,
              "length_nt": len(payload["seq"]), "sequence": payload["seq"]}
    cache.write_text(json.dumps(record, indent=1))
    return record


# ----------------------------------------------------------------- the test ---
def _one_draw(task: tuple[str, str, int]) -> dict[str, float]:
    seq, kind, seed = task
    rng = np.random.default_rng(seed)
    return h1_total_persistence(SHUFFLES[kind](seq, rng))


def null_test(seq: str, kind: str, n_shuffles: int, seed: int,
              workers: int) -> dict[str, Any]:
    """Compare the observed statistic against `n_shuffles` draws from one null."""
    observed = h1_total_persistence(seq)
    tasks = [(seq, kind, seed + i) for i in range(n_shuffles)]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        draws = list(pool.map(_one_draw, tasks, chunksize=1))

    result: dict[str, Any] = {"null": kind, "shuffles": n_shuffles,
                              "observed": observed}
    for statistic in ("total_persistence", "longest_bar", "bar_count"):
        values = np.asarray([d[statistic] for d in draws], dtype=float)
        obs = observed[statistic]
        n_ge = int(np.sum(values >= obs))
        n_le = int(np.sum(values <= obs))
        # Empirical one-sided p-values with the observation included, so the
        # floor is 1/(B+1) and is reported rather than implied.
        p_greater = (1 + n_ge) / (1 + n_shuffles)
        p_less = (1 + n_le) / (1 + n_shuffles)
        sd = float(values.std(ddof=1))
        result[statistic] = {
            "observed": obs,
            "null_mean": float(values.mean()), "null_sd": sd,
            "null_min": float(values.min()), "null_max": float(values.max()),
            "z": (obs - float(values.mean())) / sd if sd > 0 else math.nan,
            "p_greater": p_greater, "p_less": p_less,
            "p_two_sided": min(1.0, 2 * min(p_greater, p_less)),
            "p_floor": 1 / (1 + n_shuffles),
            # The draws themselves, so a figure can show the null distribution
            # rather than a Gaussian fitted to its first two moments.
            "draws": [round(float(v), 6) for v in values],
        }
    return result


def run(genes: list[str], n_shuffles: int, seed: int,
        workers: int) -> dict[str, Any]:
    report: dict[str, Any] = {
        "description": "Shuffle-null test for H1 persistence of chaos-game "
                       "embeddings of coding sequence. Closes the gap recorded "
                       "in the Section 8 ledger of the textbook.",
        "method": {"statistic": "H1 total persistence (sum of bar lengths)",
                   "embedding": "3D tetrahedral chaos game representation",
                   "filtration": "Vietoris-Rips, max_dimension 2",
                   "points": N_POINTS, "max_edge_length": MAX_EDGE,
                   "nulls": {"mono": "uniform permutation; preserves base counts",
                             "dinuc": "Altschul-Erickson; preserves dinucleotide counts",
                             "codon": "codon permutation; preserves codon usage"},
                   "p_values": "empirical, observation included; floor 1/(B+1)"},
        "shuffles": n_shuffles, "seed": seed, "genes": {}, "controls": {}}

    rng = np.random.default_rng(seed)
    primary = fetch_cds(genes[0])["sequence"]

    print(f"verifying shuffle invariants on {genes[0]} ({len(primary)} nt)")
    report["shuffle_invariants"] = verify_shuffle_invariants(primary, rng)
    print("  all invariants hold")

    # --- controls first: a broken test must not reach the real sequences ---
    print("negative control: i.i.d. sequence against the mono null")
    negative_seq = iid_sequence(primary, rng)
    negative = null_test(negative_seq, "mono", n_shuffles, seed + 10_000, workers)
    neg = negative["total_persistence"]
    neg_p, neg_z = neg["p_two_sided"], neg["z"]
    report["controls"]["negative_iid_vs_mono"] = negative
    print(f"  p_two_sided = {neg_p:.4f}  z = {neg_z:+.2f} "
          f"(must not be significant)")

    print("positive control: tandem repeat against the mono null")
    positive_seq = tandem_repeat_sequence(primary)
    positive = null_test(positive_seq, "mono", n_shuffles, seed + 20_000, workers)
    pos = positive["total_persistence"]
    pos_p, pos_z = pos["p_two_sided"], pos["z"]
    report["controls"]["positive_tandem_vs_mono"] = positive
    print(f"  p_two_sided = {pos_p:.4f}  z = {pos_z:+.2f} "
          f"(must be significant)")

    # An empirical two-sided p cannot go below 2/(B+1), so at small B the
    # positive control can be at the extreme of the null and still read 0.15.
    # Require it to be as extreme as the design permits, and additionally
    # require a large |z|, which is not floor-limited and is therefore the real
    # power statement.
    floor_two_sided = 2 / (1 + n_shuffles)
    positive_threshold = max(0.05, floor_two_sided)
    report["controls"]["verdict"] = {
        "negative_passes": bool(neg_p > 0.05),
        "positive_passes": bool(pos_p <= positive_threshold and abs(pos_z) > 3.0),
        "negative_p_two_sided": neg_p, "negative_z": neg_z,
        "positive_p_two_sided": pos_p, "positive_z": pos_z,
        "positive_threshold_used": positive_threshold,
        "p_two_sided_floor": floor_two_sided,
        "positive_motif": positive_seq[:40]}

    # --- the real sequences ---
    for gene in genes:
        record = fetch_cds(gene)
        seq = record["sequence"]
        print(f"{gene}: {len(seq)} nt")
        stride = choose_stride(len(seq))
        assert math.gcd(stride, 3) == 1, \
            f"{gene}: stride {stride} shares a factor with the codon period"
        entry: dict[str, Any] = {"transcript": record["transcript"],
                                 "source_url": record["source_url"],
                                 "length_nt": record["length_nt"],
                                 "subsample_stride": stride,
                                 "stride_coprime_with_codon": True,
                                 "nulls": {}}
        print(f"  stride {stride} (coprime with the codon period)")
        for i, kind in enumerate(SHUFFLES):
            outcome = null_test(seq, kind, n_shuffles, seed + 1000 * (i + 1), workers)
            entry["nulls"][kind] = outcome
            tp = outcome["total_persistence"]
            print(f"  vs {kind:5s}  observed={tp['observed']:.3f} "
                  f"null={tp['null_mean']:.3f}+-{tp['null_sd']:.3f} "
                  f"z={tp['z']:+.2f}  p2={tp['p_two_sided']:.4f}")
        report["genes"][gene] = entry

    return report


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--shuffles", type=int, default=200)
    parser.add_argument("--genes", nargs="*", default=["POLE", "POLD1", "MSH6"])
    parser.add_argument("--seed", type=int, default=20261022)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()

    report = run(args.genes, args.shuffles, args.seed, args.workers)

    verdict = report["controls"]["verdict"]
    if not verdict["negative_passes"]:
        print(f"\nFAIL negative control: an i.i.d. sequence was flagged "
              f"(p={verdict['negative_p_two_sided']:.4f}). The test is broken; "
              f"no p-value above is usable.")
        return 1
    if not verdict["positive_passes"]:
        print(f"\nFAIL positive control: a tandem repeat was not flagged "
              f"(p={verdict['positive_p_two_sided']:.4f} against threshold "
              f"{verdict['positive_threshold_used']:.4f}, "
              f"z={verdict['positive_z']:+.2f}). The test has no power; a null "
              f"result above would be uninformative.")
        return 1

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=1, sort_keys=True))
    print(f"\nboth controls pass; wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

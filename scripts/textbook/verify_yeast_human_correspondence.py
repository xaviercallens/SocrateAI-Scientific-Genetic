#!/usr/bin/env python3
"""Does yeast pol2-P301R actually model human POLE P286R?

The textbook's refutation runs through a yeast allele. Section 7's ledger and the
reproducibility appendix both record the correspondence human P286 <-> yeast P301
as "assumed throughout this book and established nowhere". This script
establishes it, or fails trying.

Method. Global pairwise alignment (Gotoh, affine gaps) of the full-length human
POLE catalytic subunit (UniProt Q07864, 2286 aa) against the full-length
Saccharomyces cerevisiae POL2 (P21951, 2222 aa), retrieved live with their URLs
and a UTC timestamp. We then ask which yeast residue human position 286 is
aligned to.

Why not BLOSUM62. We have no verified substitution matrix on this host, and
transcribing 400 numbers from memory is exactly the kind of unchecked assertion
this book argues against. Instead the alignment is run under several explicit
scoring schemes and gap costs, and the result is reported as *invariant* across
them or not at all. A correspondence that survives a change of scoring is better
evidence than one produced by a single matrix nobody re-derives; a correspondence
that does not survive should not be relied on, and we would rather discover that
here than assert it in a textbook.

Controls.
    positive  aligning the human sequence to itself must map 286 -> 286 under
              every scheme. If it does not, the aligner is broken.
    negative  aligning human against a shuffled yeast sequence (same amino-acid
              composition, order destroyed) must NOT recover 286 -> 301. If it
              does, the correspondence is an artefact of composition rather than
              of homology.

Both are asserted. The script writes nothing and exits non-zero if either fails.

    python scripts/textbook/verify_yeast_human_correspondence.py
    python scripts/textbook/verify_yeast_human_correspondence.py --offline
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "textbook" / "data" / "orthologue_sequences.json"
OUT = ROOT / "data" / "yeast_human_correspondence.json"

HUMAN = {"accession": "Q07864", "gene": "POLE", "residue": 286, "expect": "P"}
YEAST = {"accession": "P21951", "gene": "POL2", "residue": 301, "expect": "P"}

# Amino-acid classes used by the "class" scoring scheme. Standard
# physicochemical groupings; the point is that they are stated, not that they
# are uniquely correct.
CLASSES = {
    "hydrophobic": "AVLIMC", "aromatic": "FWY", "polar": "STNQ",
    "acidic": "DE", "basic": "KRH", "special": "GP",
}
CLASS_OF = {aa: name for name, members in CLASSES.items() for aa in members}


# -------------------------------------------------------------- retrieval ---
def fetch_sequence(accession: str) -> dict[str, Any]:
    url = (f"https://rest.uniprot.org/uniprotkb/{accession}.json"
           f"?fields=accession,id,protein_name,sequence,length")
    request = urllib.request.Request(url, headers={
        "Accept": "application/json",
        "User-Agent": "SocrateAI-Formal-Genetics-Textbook/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        payload = json.loads(response.read().decode())
    return {
        "accession": payload["primaryAccession"],
        "uniprot_id": payload.get("uniProtkbId"),
        "length_aa": payload["sequence"]["length"],
        "sequence": payload["sequence"]["value"],
        "source_url": url,
        "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    }


def load_sequences(offline: bool) -> dict[str, Any]:
    if offline:
        if not CACHE.exists():
            raise SystemExit(f"--offline but no cache at {CACHE}")
        return json.loads(CACHE.read_text())
    records = {"human": fetch_sequence(HUMAN["accession"]),
               "yeast": fetch_sequence(YEAST["accession"])}
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(records, indent=1))
    return records


# --------------------------------------------------------------- scoring ---
def score_matrix(scheme: str) -> dict[tuple[str, str], float]:
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    table: dict[tuple[str, str], float] = {}
    for a in alphabet:
        for b in alphabet:
            if a == b:
                table[(a, b)] = 4.0
            elif scheme == "identity":
                table[(a, b)] = -2.0
            elif scheme == "class":
                same = CLASS_OF.get(a) == CLASS_OF.get(b)
                table[(a, b)] = 1.0 if same else -2.0
            else:
                raise ValueError(f"unknown scheme {scheme}")
    return table


# -------------------------------------------------------------- alignment ---
def align(seq_a: str, seq_b: str, scheme: str, gap_open: float,
          gap_extend: float) -> list[tuple[int | None, int | None]]:
    """Gotoh global alignment. Returns (i, j) pairs of 0-based indices.

    Three state matrices -- M (match/mismatch), X (gap in B), Y (gap in A) --
    with affine gap costs, and an explicit traceback. Vectorised over `seq_b`.
    """
    table = score_matrix(scheme)
    n, m = len(seq_a), len(seq_b)
    neg = -1e18

    b_arr = np.frombuffer(seq_b.encode(), dtype="S1").astype("U1")

    prev_m = np.full(m + 1, neg); prev_x = np.full(m + 1, neg)
    prev_y = np.full(m + 1, neg)
    prev_m[0] = 0.0
    for j in range(1, m + 1):
        prev_y[j] = -gap_open - (j - 1) * gap_extend

    # Traceback pointers: which state we came from, per cell, per state.
    ptr_m = np.zeros((n + 1, m + 1), dtype=np.int8)
    ptr_x = np.zeros((n + 1, m + 1), dtype=np.int8)
    ptr_y = np.zeros((n + 1, m + 1), dtype=np.int8)

    # Substitution scores as a lookup indexed by byte value, so each row's
    # scores are one fancy-index rather than a Python comprehension.
    lut = np.full((128, 128), -2.0)
    for (a, b), v in table.items():
        lut[ord(a), ord(b)] = v
    a_codes = np.frombuffer(seq_a.encode(), dtype=np.uint8).astype(np.intp)
    b_codes = np.frombuffer(seq_b.encode(), dtype=np.uint8).astype(np.intp)

    arange = np.arange(m + 1, dtype=float)

    for i in range(1, n + 1):
        sub = lut[a_codes[i - 1], b_codes]

        cur_m = np.full(m + 1, neg)
        cur_y = np.full(m + 1, neg)

        # X: gap in B (consume seq_a) -- depends only on the previous row.
        open_x = prev_m - gap_open
        extend_x = prev_x - gap_extend
        cur_x = np.maximum(open_x, extend_x)
        ptr_x[i, :] = np.where(open_x >= extend_x, 0, 1)
        cur_x[0] = -gap_open - (i - 1) * gap_extend

        # M: depends only on the previous row, so it vectorises outright.
        stacked = np.stack((prev_m[:-1], prev_x[:-1], prev_y[:-1]))
        best_prev = stacked.max(axis=0)
        cur_m[1:] = best_prev + sub
        ptr_m[i, 1:] = stacked.argmax(axis=0).astype(np.int8)

        # Y: gap in A, sequential in j. A gap run of length L costs
        #   gap_open + (L-1) * gap_extend,
        # so  cur_y[j] = max_{k<j} ( cur_m[k] - gap_open - (j-k-1)*gap_extend )
        #              = -gap_open - (j-1)*gap_extend
        #                + max_{k<=j-1} ( cur_m[k] + k*gap_extend ),
        # and the running maximum is a prefix scan rather than a loop.
        shifted = cur_m[:-1] + arange[:-1] * gap_extend
        running = np.maximum.accumulate(shifted)
        cur_y[1:] = running - gap_open - (arange[1:] - 1) * gap_extend
        # Opened here if the best predecessor is the immediately preceding M.
        opened = cur_m[:-1] - gap_open
        ptr_y[i, 1:] = np.where(cur_y[1:] <= opened + 1e-9, 0, 2).astype(np.int8)

        prev_m, prev_x, prev_y = cur_m, cur_x, cur_y

    # Traceback from the best end state.
    i, j = n, m
    finals = [prev_m[m], prev_x[m], prev_y[m]]
    state = int(np.argmax(finals))
    pairs: list[tuple[int | None, int | None]] = []
    while i > 0 or j > 0:
        if state == 0 and i > 0 and j > 0:
            pairs.append((i - 1, j - 1))
            state = int(ptr_m[i, j]); i -= 1; j -= 1
        elif state == 1 and i > 0:
            pairs.append((i - 1, None))
            state = int(ptr_x[i, j]); i -= 1
        elif state == 2 and j > 0:
            pairs.append((None, j - 1))
            state = int(ptr_y[i, j]); j -= 1
        elif i > 0:
            pairs.append((i - 1, None)); i -= 1
        else:
            pairs.append((None, j - 1)); j -= 1
    pairs.reverse()
    return pairs


def partner_of(pairs: list[tuple[int | None, int | None]],
               index_a: int) -> int | None:
    """Which index in B is aligned to index_a in A (None if a gap)?"""
    for ia, ib in pairs:
        if ia == index_a:
            return ib
    return None


def identity_fraction(pairs: list[tuple[int | None, int | None]],
                      seq_a: str, seq_b: str) -> float:
    aligned = [(ia, ib) for ia, ib in pairs if ia is not None and ib is not None]
    if not aligned:
        return 0.0
    same = sum(1 for ia, ib in aligned if seq_a[ia] == seq_b[ib])
    return same / len(aligned)


def local_identity(pairs: list[tuple[int | None, int | None]], seq_a: str,
                   seq_b: str, centre_a: int, half_width: int = 20
                   ) -> tuple[float, int]:
    """Identity in a window of `seq_a` around `centre_a`, over aligned pairs.

    This, not the aligned position itself, is the statistic that discriminates
    homology from chance here. A global alignment of two sequences of similar
    length maps position i to roughly position i whatever the sequences say, so
    "286 lands on 301" is weak evidence on its own -- the shuffled control lands
    on 302, one residue away. Local conservation has no such diagonal prior.
    """
    window = [(ia, ib) for ia, ib in pairs
              if ia is not None and ib is not None
              and abs(ia - centre_a) <= half_width]
    if not window:
        return 0.0, 0
    same = sum(1 for ia, ib in window if seq_a[ia] == seq_b[ib])
    return 100 * same / len(window), len(window)


# ------------------------------------------------------------------- main ---
SCHEMES = [("identity", 11.0, 1.0), ("identity", 8.0, 2.0),
           ("class", 11.0, 1.0), ("class", 8.0, 2.0)]


def run(offline: bool) -> dict[str, Any]:
    records = load_sequences(offline)
    human, yeast = records["human"]["sequence"], records["yeast"]["sequence"]

    report: dict[str, Any] = {
        "description": "Verification that yeast POL2 P301 is the orthologous "
                       "position of human POLE P286. Closes the gap recorded in "
                       "the Section 7 ledger of the textbook.",
        "method": {
            "alignment": "Gotoh global alignment, affine gaps, full length",
            "scoring": "no verified substitution matrix available on this host; "
                       "the correspondence is required to be invariant across "
                       "the schemes and gap costs listed",
            "schemes": [{"scheme": s, "gap_open": go, "gap_extend": ge}
                        for s, go, ge in SCHEMES],
        },
        "sequences": {k: {kk: vv for kk, vv in v.items() if kk != "sequence"}
                      for k, v in records.items()},
    }

    # The residues must be what the literature says before anything else runs.
    human_aa = human[HUMAN["residue"] - 1]
    yeast_aa = yeast[YEAST["residue"] - 1]
    assert human_aa == HUMAN["expect"], \
        f"human {HUMAN['accession']} residue {HUMAN['residue']} is {human_aa}"
    assert yeast_aa == YEAST["expect"], \
        f"yeast {YEAST['accession']} residue {YEAST['residue']} is {yeast_aa}"
    report["residue_identity_check"] = {
        "human_286": human_aa, "yeast_301": yeast_aa, "both_proline": True}
    print(f"human {HUMAN['accession']} 286 = {human_aa}; "
          f"yeast {YEAST['accession']} 301 = {yeast_aa}")

    # --- the alignments ---
    print("aligning full-length sequences under each scheme")
    runs = []
    for scheme, gap_open, gap_extend in SCHEMES:
        pairs = align(human, yeast, scheme, gap_open, gap_extend)
        partner = partner_of(pairs, HUMAN["residue"] - 1)
        partner_1based = None if partner is None else partner + 1
        local, window = local_identity(pairs, human, yeast, HUMAN["residue"] - 1)
        entry = {
            "scheme": scheme, "gap_open": gap_open, "gap_extend": gap_extend,
            "human_286_aligns_to_yeast": partner_1based,
            "aligned_yeast_residue": None if partner is None else yeast[partner],
            "matches_P301": partner_1based == YEAST["residue"],
            "percent_identity": round(100 * identity_fraction(pairs, human, yeast), 2),
            "local_percent_identity": round(local, 2),
            "local_window_pairs": window,
        }
        runs.append(entry)
        print(f"  {scheme:8s} open={gap_open:<5} ext={gap_extend:<4} "
              f"286 -> {partner_1based} ({entry['aligned_yeast_residue']}) "
              f"global={entry['percent_identity']}% "
              f"local={entry['local_percent_identity']}%")
    report["alignments"] = runs
    report["invariant_across_schemes"] = all(r["matches_P301"] for r in runs)

    # --- controls ---
    print("positive control: human against itself")
    positive = []
    for scheme, gap_open, gap_extend in SCHEMES:
        pairs = align(human, human, scheme, gap_open, gap_extend)
        partner = partner_of(pairs, HUMAN["residue"] - 1)
        positive.append({"scheme": scheme, "gap_open": gap_open,
                         "self_maps_286_to": None if partner is None else partner + 1})
    report["controls"] = {"positive_self_alignment": positive}
    positive_ok = all(p["self_maps_286_to"] == HUMAN["residue"] for p in positive)
    print(f"  286 -> {[p['self_maps_286_to'] for p in positive]} "
          f"({'pass' if positive_ok else 'FAIL'})")

    print("negative control: human against a composition-matched shuffle of yeast")
    rng = np.random.default_rng(20261022)
    letters = list(yeast)
    rng.shuffle(letters)
    shuffled = "".join(letters)
    assert Counter(shuffled) == Counter(yeast), "shuffle changed composition"
    negative = []
    for scheme, gap_open, gap_extend in SCHEMES:
        pairs = align(human, shuffled, scheme, gap_open, gap_extend)
        partner = partner_of(pairs, HUMAN["residue"] - 1)
        local, window = local_identity(pairs, human, shuffled, HUMAN["residue"] - 1)
        negative.append({
            "scheme": scheme, "gap_open": gap_open,
            "shuffled_maps_286_to": None if partner is None else partner + 1,
            "percent_identity": round(100 * identity_fraction(pairs, human, shuffled), 2),
            "local_percent_identity": round(local, 2),
            "local_window_pairs": window})
    report["controls"]["negative_shuffled_yeast"] = negative

    # The discriminating statistic is local conservation, not the aligned
    # position: the shuffle lands on 302, one residue off, purely because a
    # global alignment of similar-length sequences has a strong diagonal prior.
    # Require the real local identity to beat the shuffled one by 2x under
    # every scheme, which is a claim about homology rather than about position.
    worst_real = min(r["local_percent_identity"] for r in runs)
    best_shuffled = max(n["local_percent_identity"] for n in negative)
    negative_ok = bool(worst_real > 2.0 * best_shuffled)
    report["controls"]["local_conservation_margin"] = {
        "worst_real_local_identity": worst_real,
        "best_shuffled_local_identity": best_shuffled,
        "ratio": round(worst_real / best_shuffled, 2) if best_shuffled else None,
        "criterion": "real local identity > 2x shuffled, under every scheme",
        "note": "the aligned position alone is weak evidence: the shuffled "
                "control lands on 302, adjacent to 301, because global "
                "alignment of similar-length sequences follows the diagonal",
    }
    print(f"  286 -> {[n['shuffled_maps_286_to'] for n in negative]} "
          f"global={[n['percent_identity'] for n in negative]}% "
          f"local={[n['local_percent_identity'] for n in negative]}%")
    print(f"  local conservation: real {worst_real}% vs shuffled "
          f"{best_shuffled}% ({'pass' if negative_ok else 'FAIL'})")

    # Sequence context, for the textbook to print.
    report["context"] = {
        "human_279_293": human[278:293],
        "yeast_294_308": yeast[293:308],
        "note": "the ExoI region; human numbering 279-293, yeast 294-308",
    }
    report["controls"]["verdict"] = {
        "positive_passes": bool(positive_ok),
        "negative_passes": bool(negative_ok),
    }
    return report


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--offline", action="store_true",
                        help="use the cached sequences instead of fetching")
    args = parser.parse_args()

    report = run(args.offline)
    verdict = report["controls"]["verdict"]

    if not verdict["positive_passes"]:
        print("\nFAIL positive control: self-alignment does not fix residue 286. "
              "The aligner is broken and no result above is usable.")
        return 1
    if not verdict["negative_passes"]:
        margin = report["controls"]["local_conservation_margin"]
        print(f"\nFAIL negative control: local conservation around the site is "
              f"{margin['worst_real_local_identity']}% for the real pair against "
              f"{margin['best_shuffled_local_identity']}% for a shuffle. The "
              f"correspondence would be an artefact of composition and length, "
              f"not homology.")
        return 1
    if not report["invariant_across_schemes"]:
        print("\nThe correspondence is NOT invariant across scoring schemes. "
              "Writing the result, because this is a finding and not an error: "
              "the textbook must not rely on it.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=1, sort_keys=True))
    print(f"\ncontrols pass; wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

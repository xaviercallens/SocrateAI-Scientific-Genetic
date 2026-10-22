#!/usr/bin/env python3
"""Extract homologous viral populations from the Hugging Face Virus-Host-Genomes set.

The Reticulata test needs what Phase 2 did not have: a *population* of real,
homologous, comparably long sequences, with an independent expectation about
whether its genealogy is tree-like. This script builds such populations from

    https://huggingface.co/datasets/hiyata/Virus-Host-Genomes

(58 046 viral genomes with taxonomy; `test.parquet` is 50 MB and enough).
Download it once:

    curl -L -o test.parquet \\
      https://huggingface.co/datasets/hiyata/Virus-Host-Genomes/resolve/main/test.parquet

then

    python scripts/reticulata_hf_populations.py --parquet test.parquet \\
      --out data/reticulata_phase3/populations

Selection is pre-registered and deterministic, because the choice of which
sequences enter a population is exactly the kind of decision that can
manufacture a result. Two strategies are written for every taxon and both are
reported:

* ``accession`` -- sort the surviving sequences by accession and take the first
  *n*. Arbitrary with respect to phylogeny, and chosen before any margin was
  computed.
* ``spread``    -- greedy max-min 6-mer diversity, seeded from the
  accession-first sequence. Deliberately harder: more divergence means more
  alignment error and more saturation.

Nothing here computes a distance on unaligned sequences; the 6-mer profile is
used only to *choose* sequences, and the FASTA it writes is the input to the
aligner, not to the test.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import pandas as pd

# Populations, with the expectation that was fixed before any test was run.
# `expectation` is a biological prior from the literature, not a label derived
# from this data, and the test is allowed to contradict it.
POPULATIONS: Dict[str, Dict[str, str]] = {
    "hepatitis_b": {
        "virus_name": "Hepatitis B virus",
        "expectation": "reticulate",
        "why": (
            "Inter-genotype recombination in HBV is textbook: genotype I and "
            "genotype C/D recombinants are recognised taxonomic entities, and "
            "the 3.2 kb genome has overlapping reading frames throughout."
        ),
    },
    "coxsackievirus_a6": {
        "virus_name": "Coxsackievirus A6",
        "expectation": "reticulate",
        "why": (
            "Enterovirus A recombination is concentrated in the non-structural "
            "3' half; CVA6 outbreak lineages are named by their recombinant "
            "form (RF-A ... RF-E)."
        ),
    },
    "rabies": {
        "virus_name": "Lyssavirus rabies",
        "expectation": "tree-like",
        "why": (
            "Rabies virus is a non-segmented negative-strand RNA virus; "
            "recombination is reported only exceptionally and the species is "
            "normally analysed under a strictly bifurcating phylogeny."
        ),
    },
    "hpv52": {
        "virus_name": "human papillomavirus 52",
        "expectation": "tree-like",
        "why": (
            "Papillomaviruses co-diverge with their hosts and recombination "
            "within a type is considered rare; HPV52 sublineages are defined "
            "by a tree."
        ),
    },
    "jc_polyomavirus": {
        "virus_name": "JC polyomavirus",
        "expectation": "tree-like",
        "why": (
            "JCV genotypes track human migration and are analysed as a tree; "
            "recombination is reported mainly in the regulatory region, not "
            "genome-wide."
        ),
    },
}

LENGTH_TOLERANCE: float = 0.02
MAX_AMBIGUOUS_FRACTION: float = 0.001
KMER: int = 6


def _ambiguous_fraction(sequence: str) -> float:
    """Fraction of residues that are not an unambiguous A, C, G or T."""
    upper = sequence.upper()
    good = sum(upper.count(base) for base in "ACGT")
    return 1.0 - good / len(upper) if upper else 1.0


def _kmer_profile(sequence: str, k: int = KMER) -> Counter:
    """Multiset of k-mers, used only to choose which sequences to include."""
    upper = sequence.upper()
    return Counter(upper[i : i + k] for i in range(len(upper) - k + 1))


def _profile_distance(left: Counter, right: Counter) -> float:
    """Normalised L1 distance between two k-mer profiles, in [0, 1]."""
    total_left = sum(left.values())
    total_right = sum(right.values())
    if total_left == 0 or total_right == 0:
        return 1.0
    keys = set(left) | set(right)
    l1 = sum(
        abs(left.get(key, 0) / total_left - right.get(key, 0) / total_right)
        for key in keys
    )
    return l1 / 2.0


def filter_population(frame: pd.DataFrame, virus_name: str) -> pd.DataFrame:
    """Apply the pre-registered quality filters to one taxon.

    Keeps near-median-length, unambiguous, deduplicated sequences. The
    deduplication is by exact sequence: two identical genomes carry no
    information about tree-likeness and would only pad the quartet count.
    """
    subset = frame[frame["virus_name"] == virus_name].copy()
    if subset.empty:
        return subset
    subset["length"] = subset["sequence"].str.len()
    median = float(subset["length"].median())
    low, high = median * (1.0 - LENGTH_TOLERANCE), median * (1.0 + LENGTH_TOLERANCE)
    subset = subset[(subset["length"] >= low) & (subset["length"] <= high)]
    subset = subset[
        subset["sequence"].map(_ambiguous_fraction) <= MAX_AMBIGUOUS_FRACTION
    ]
    subset = subset.drop_duplicates(subset=["sequence"])
    return subset.sort_values("accession", kind="mergesort").reset_index(drop=True)


def select_accession(subset: pd.DataFrame, n: int) -> List[int]:
    """The first *n* rows in accession order."""
    return list(range(min(n, len(subset))))


def select_spread(subset: pd.DataFrame, n: int) -> List[int]:
    """Greedy max-min 6-mer diversity, seeded from the accession-first row.

    Deterministic: ties are broken by the (already accession-sorted) row index.
    """
    if len(subset) <= n:
        return list(range(len(subset)))
    profiles = [_kmer_profile(seq) for seq in subset["sequence"]]
    chosen: List[int] = [0]
    best: List[float] = [_profile_distance(profiles[0], profile) for profile in profiles]
    while len(chosen) < n:
        candidate = -1
        candidate_score = -1.0
        for index in range(len(profiles)):
            if index in chosen:
                continue
            if best[index] > candidate_score:
                candidate_score = best[index]
                candidate = index
        chosen.append(candidate)
        for index in range(len(profiles)):
            distance = _profile_distance(profiles[candidate], profiles[index])
            if distance < best[index]:
                best[index] = distance
    return sorted(chosen)


def write_fasta(path: Path, labels: Sequence[str], sequences: Sequence[str]) -> str:
    """Write a FASTA file and return the SHA-256 of its bytes."""
    lines: List[str] = []
    for label, sequence in zip(labels, sequences):
        lines.append(f">{label}")
        upper = sequence.upper()
        for start in range(0, len(upper), 80):
            lines.append(upper[start : start + 80])
    payload = ("\n".join(lines) + "\n").encode("ascii")
    path.write_bytes(payload)
    return hashlib.sha256(payload).hexdigest()


def pairwise_kmer_summary(sequences: Iterable[str]) -> Tuple[float, float]:
    """Mean and maximum 6-mer profile distance within a population.

    A coarse pre-alignment diversity read, reported so that a reader can see
    which population the aligner was handed before any distance is computed on
    columns.
    """
    profiles = [_kmer_profile(sequence) for sequence in sequences]
    values: List[float] = []
    for i in range(len(profiles)):
        for j in range(i + 1, len(profiles)):
            values.append(_profile_distance(profiles[i], profiles[j]))
    if not values:
        return 0.0, 0.0
    return sum(values) / len(values), max(values)


def build(parquet: Path, out_dir: Path, n: int) -> Dict[str, object]:
    """Build every population at both strategies and return the provenance record."""
    frame = pd.read_parquet(parquet)
    out_dir.mkdir(parents=True, exist_ok=True)
    record: Dict[str, object] = {
        "source": "huggingface:hiyata/Virus-Host-Genomes",
        "source_file": parquet.name,
        "source_sha256": hashlib.sha256(parquet.read_bytes()).hexdigest(),
        "rows_in_file": int(len(frame)),
        "filters": {
            "length_tolerance_about_median": LENGTH_TOLERANCE,
            "max_ambiguous_fraction": MAX_AMBIGUOUS_FRACTION,
            "deduplicated": "exact sequence",
        },
        "taxa_per_population": n,
        "populations": [],
    }
    populations: List[Dict[str, object]] = []
    for slug, spec in POPULATIONS.items():
        subset = filter_population(frame, spec["virus_name"])
        entry: Dict[str, object] = {
            "slug": slug,
            "virus_name": spec["virus_name"],
            "expectation": spec["expectation"],
            "why": spec["why"],
            "surviving_sequences": int(len(subset)),
            "median_length": int(subset["length"].median()) if len(subset) else 0,
            "files": {},
        }
        for strategy, selector in (
            ("accession", select_accession),
            ("spread", select_spread),
        ):
            picked = selector(subset, n)
            rows = subset.iloc[picked]
            labels = [str(accession) for accession in rows["accession"]]
            sequences = [str(sequence) for sequence in rows["sequence"]]
            path = out_dir / f"{slug}_{strategy}.fasta"
            digest = write_fasta(path, labels, sequences)
            mean_kmer, max_kmer = pairwise_kmer_summary(sequences)
            entry["files"][strategy] = {
                "file": path.name,
                "sha256": digest,
                "taxa": len(labels),
                "accessions": labels,
                "lengths": [len(sequence) for sequence in sequences],
                "mean_kmer_distance": round(mean_kmer, 6),
                "max_kmer_distance": round(max_kmer, 6),
            }
        populations.append(entry)
    record["populations"] = populations
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parquet", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--taxa", type=int, default=16)
    args = parser.parse_args()
    record = build(args.parquet, args.out, args.taxa)
    manifest = args.out / "provenance.json"
    manifest.write_text(json.dumps(record, indent=2) + "\n")
    for entry in record["populations"]:  # type: ignore[index]
        files = entry["files"]  # type: ignore[index]
        print(
            f"{entry['slug']:18s} {entry['virus_name']:28s} "  # type: ignore[index]
            f"len={entry['median_length']:6d} pool={entry['surviving_sequences']:4d} "  # type: ignore[index]
            f"accession_kmer={files['accession']['mean_kmer_distance']:.4f} "
            f"spread_kmer={files['spread']['mean_kmer_distance']:.4f}"
        )
    print(f"wrote {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

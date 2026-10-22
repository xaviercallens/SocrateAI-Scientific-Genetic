#!/usr/bin/env python3
"""Build SARS-CoV-2 populations whose recombination has a *published breakpoint*.

This is the dataset the Reticulata programme was missing. Phase 3's ten viral
populations carried only a literature *expectation* about tree-likeness, and
they had to pass through an aligner we wrote ourselves, whose error we then had
to measure and gate on. Here both problems change shape:

* **The recombinants are named and their breakpoints are published.** The Pango
  designation list records recombinant lineages as ``X*`` with their parent
  lineages and, for many, an explicit numeric breakpoint interval in reference
  coordinates -- each traceable to a reviewed GitHub issue. That is ground truth
  of a different quality from "the literature expects this family to recombine".
* **The data arrives already aligned**, from Nextstrain's open (non-GISAID)
  release, aligned by Nextclade to the MN908947 reference. So the breakpoint
  coordinates and the alignment share one coordinate system, and *our* aligner
  is not involved at all.

The second point needs an honest caveat, stated here rather than buried. Using a
published alignment removes **our** aligner from the trust surface and replaces
it with **theirs**, which is reference-anchored too -- Nextclade projects onto
MN908947 and strips insertions, the same construction and therefore the same
class of error that `reticulation align-error` measures. What improves is
provenance: the alignment is produced independently of this project, by a
standard tool, and is citable. What does not improve is coverage of the
measurement: we have not measured Nextclade's error, and the admissibility gate
calibrated for our own aligner does not transfer to it. A population built here
should be read with that trade in mind.

Usage:

    python scripts/reticulata_sarscov2_breakpoint_populations.py \\
      --out data/reticulata_phase4/populations

Streams about 2.2 GB and stores only what it selects. Nothing is cached.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Sequence

METADATA_URL = "https://data.nextstrain.org/files/ncov/open/metadata.tsv.zst"
SEQUENCES_URL = "https://data.nextstrain.org/files/ncov/open/aligned.fasta.zst"
LINEAGE_NOTES_URL = (
    "https://raw.githubusercontent.com/cov-lineages/pango-designation/"
    "master/lineage_notes.txt"
)
REFERENCE = "MN908947"

# Triads fixed before any margin was computed. `breakpoint` is quoted verbatim
# from the Pango lineage notes; `issue` is the reviewed proposal behind it.
TRIADS: Dict[str, Dict[str, object]] = {
    "xcm": {
        "recombinant": "XCM",
        "parents": ("XBB.2.3", "DV.7.1"),
        "breakpoint": "breakpoint between 22017-22031",
        "issue": "sars-cov-2-variants/lineage-proposals#674",
        "why": (
            "The tightest published breakpoint interval we found -- 14 "
            "nucleotides -- with all three lineages well sampled."
        ),
    },
    "xcp": {
        "recombinant": "XCP",
        "parents": ("EG.5.1.10", "JD.1.1"),
        "breakpoint": "breakpoint between 14190-14633",
        "issue": "sars-cov-2-variants/lineage-proposals#772",
        "why": (
            "A breakpoint in ORF1b rather than in S, so the two inherited "
            "segments are of very different lengths. One parent lineage is "
            "thinly sampled."
        ),
    },
}

# Quality filters, fixed in advance.
MIN_COVERAGE = 0.97
MAX_MISSING = 1000

PER_PARENT = 5
PER_RECOMBINANT = 6


def stream_metadata(lineages: Sequence[str]) -> List[Dict[str, str]]:
    """Stream the metadata once, keeping only rows in `lineages`.

    The file is ~950 MB compressed and several gigabytes expanded, so it is
    filtered in the pipe and never stored.
    """
    wanted = set(lineages)
    script = r"""
BEGIN { FS="\t"; OFS="\t" }
NR==1 {
  for (i = 1; i <= NF; i++) col[$i] = i
  next
}
{
  lin = $(col["Nextclade_pango"])
  if (lin == "" || lin == "?") lin = $(col["pango_lineage"])
  if (lin in want)
    print $(col["strain"]), lin, $(col["length"]), \
          $(col["missing_data"]), $(col["coverage"]), $(col["nonACGTN"]), \
          $(col["genbank_accession"]), $(col["date"])
}
"""
    assignments = "".join(f'want["{name}"]=1;' for name in sorted(wanted))
    program = script.replace("BEGIN { FS=\"\\t\"; OFS=\"\\t\" }",
                             f'BEGIN {{ FS="\\t"; OFS="\\t"; {assignments} }}')
    pipeline = (
        f"curl -sL {METADATA_URL} | zstd -dc 2>/dev/null | awk '{program}'"
    )
    done = subprocess.run(
        ["bash", "-o", "pipefail", "-c", pipeline],
        capture_output=True,
        text=True,
        timeout=3600,
    )
    if done.returncode != 0:
        raise RuntimeError(f"metadata stream failed: {done.stderr[:400]}")
    rows: List[Dict[str, str]] = []
    for line in done.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) != 8:
            continue
        rows.append(
            {
                "strain": parts[0],
                "lineage": parts[1],
                "length": parts[2],
                "missing_data": parts[3],
                "coverage": parts[4],
                "nonACGTN": parts[5],
                "genbank_accession": parts[6],
                "date": parts[7],
            }
        )
    return rows


def passes_quality(row: Dict[str, str]) -> bool:
    """The pre-registered quality filter."""
    try:
        coverage = float(row["coverage"])
        missing = float(row["missing_data"])
    except ValueError:
        return False
    return coverage >= MIN_COVERAGE and missing <= MAX_MISSING


def select(rows: Sequence[Dict[str, str]], lineage: str, count: int) -> List[Dict[str, str]]:
    """The first `count` qualifying strains in strain-name order.

    Deterministic and arbitrary with respect to phylogeny, which is the point:
    the selection rule was fixed before any margin was computed.
    """
    pool = sorted(
        (row for row in rows if row["lineage"] == lineage and passes_quality(row)),
        key=lambda row: row["strain"],
    )
    return pool[:count]


def stream_sequences(strains: Sequence[str], out: Path) -> int:
    """Stream the aligned release once, writing only the selected records."""
    listing = out.with_suffix(".strains")
    listing.write_text("\n".join(strains) + "\n")
    program = r"""
BEGIN { while ((getline line < WANT) > 0) want[line] = 1; keep = 0 }
/^>/ { name = substr($0, 2); sub(/[ \t].*/, "", name); keep = (name in want);
       if (keep) print ">" name; next }
keep { print }
"""
    pipeline = (
        f"curl -sL {SEQUENCES_URL} | zstd -dc 2>/dev/null | "
        f"awk -v WANT={listing} '{program}' > {out}"
    )
    done = subprocess.run(
        ["bash", "-o", "pipefail", "-c", pipeline],
        capture_output=True,
        text=True,
        timeout=7200,
    )
    if done.returncode != 0:
        raise RuntimeError(f"sequence stream failed: {done.stderr[:400]}")
    listing.unlink()
    return out.read_text().count(">")


def lineage_note(recombinant: str) -> str:
    """The verbatim Pango note for a recombinant, fetched rather than recalled."""
    done = subprocess.run(
        ["curl", "-sL", LINEAGE_NOTES_URL],
        capture_output=True,
        text=True,
        timeout=600,
    )
    if done.returncode != 0:
        raise RuntimeError("could not fetch the lineage notes")
    for line in done.stdout.splitlines():
        if re.match(rf"^{re.escape(recombinant)}\t", line):
            return line.split("\t", 1)[1].strip()
    raise RuntimeError(f"{recombinant} is not in the lineage notes")


def build(out_dir: Path, which: Sequence[str]) -> Dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    lineages: List[str] = []
    for key in which:
        spec = TRIADS[key]
        lineages.append(str(spec["recombinant"]))
        lineages.extend(spec["parents"])  # type: ignore[arg-type]

    print(f"streaming metadata for {len(lineages)} lineages ...")
    rows = stream_metadata(lineages)
    print(f"  {len(rows)} rows matched")

    record: Dict[str, object] = {
        "source": {
            "sequences": SEQUENCES_URL,
            "metadata": METADATA_URL,
            "lineage_notes": LINEAGE_NOTES_URL,
            "alignment": (
                f"already aligned by Nextclade to {REFERENCE}; this project's "
                "aligner is not involved"
            ),
            "caveat": (
                "Using a published alignment removes OUR aligner from the trust "
                "surface and replaces it with Nextclade's, which is "
                "reference-anchored too and whose error we have not measured. "
                "Provenance improves; measurement coverage does not."
            ),
        },
        "filters": {
            "min_coverage": MIN_COVERAGE,
            "max_missing_data": MAX_MISSING,
            "selection": "first N qualifying strains in strain-name order",
        },
        "populations": [],
    }

    populations: List[Dict[str, object]] = []
    for key in which:
        spec = TRIADS[key]
        recombinant = str(spec["recombinant"])
        parent_a, parent_b = spec["parents"]  # type: ignore[misc]
        note = lineage_note(recombinant)
        picked = (
            select(rows, parent_a, PER_PARENT)
            + select(rows, parent_b, PER_PARENT)
            + select(rows, recombinant, PER_RECOMBINANT)
        )
        shortfall = {
            parent_a: PER_PARENT - len(select(rows, parent_a, PER_PARENT)),
            parent_b: PER_PARENT - len(select(rows, parent_b, PER_PARENT)),
            recombinant: PER_RECOMBINANT
            - len(select(rows, recombinant, PER_RECOMBINANT)),
        }
        path = out_dir / f"{key}_{recombinant}_triad.fasta"
        print(f"streaming sequences for {recombinant} ({len(picked)} strains) ...")
        written = stream_sequences([row["strain"] for row in picked], path)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        populations.append(
            {
                "key": key,
                "file": path.name,
                "sha256": digest,
                "recombinant_lineage": recombinant,
                "parent_lineages": [parent_a, parent_b],
                "published_breakpoint": spec["breakpoint"],
                "pango_note_verbatim": note,
                "proposal_issue": spec["issue"],
                "why_this_triad": spec["why"],
                "records_requested": len(picked),
                "records_written": written,
                "shortfall_by_lineage": shortfall,
                "members": [
                    {
                        "strain": row["strain"],
                        "lineage": row["lineage"],
                        "genbank_accession": row["genbank_accession"],
                        "date": row["date"],
                        "coverage": row["coverage"],
                        "missing_data": row["missing_data"],
                    }
                    for row in picked
                ],
                "prediction": (
                    "Pre-registered: the quartet certificate's dominant pairing "
                    "should separate the two PARENT lineages, with recombinant "
                    "members on the other side -- the planted-recombinant "
                    "signature. A dominant pairing that instead separates two "
                    "members of the same lineage would contradict this."
                ),
            }
        )
    record["populations"] = populations
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--triads",
        nargs="+",
        default=["xcm"],
        choices=sorted(TRIADS),
        help="which published-breakpoint triads to build",
    )
    args = parser.parse_args()
    record = build(args.out, args.triads)
    manifest = args.out / "provenance.json"
    manifest.write_text(json.dumps(record, indent=2) + "\n")
    for entry in record["populations"]:  # type: ignore[index]
        print(
            f"{entry['file']}: {entry['records_written']} records, "  # type: ignore[index]
            f"breakpoint {entry['published_breakpoint']!r}"  # type: ignore[index]
        )
    print(f"wrote {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

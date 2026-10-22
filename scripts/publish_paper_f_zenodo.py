"""Deposit Paper F as a new Zenodo record.

Paper F reports the tree-level null, the measured cost of the alignment step,
and the Gromov hinge of Buneman's converse. First deposit, so a new record;
later editions must go through `/actions/newversion` on the concept DOI this run
returns.

Two rules this script enforces rather than trusts:

* **No peer-review artefact enters the bundle.** The repository once described
  work as "5-round Opus peer reviewed" on the strength of a function that made
  no model call. Any file whose name mentions peer review aborts the run, and
  the description states plainly that the book has not been peer reviewed.
* **Publication is verified, never assumed.** The Zenodo API has been observed
  returning 5xx for writes that in fact succeeded, so a failed publish is
  followed by a state read, and every uploaded file's MD5 is compared against
  the locally computed one before publishing. A blind retry would risk a
  duplicate record.

Usage:

    python scripts/publish_paper_f_zenodo.py --dry-run   # bundle + checks only
    python scripts/publish_paper_f_zenodo.py             # deposit and publish
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[2]
ZENODO = "https://zenodo.org/api"
RECEIPT = HERE / "data" / "zenodo_paper_f_edition2_receipt.json"
STEM = "PaperF_Reticulate_Evolution_TreeLevel_Null_Lean4"

TITLE = (
    "Proof-Carrying Detection of Reticulate Evolution: a Tree-Level Null, the "
    "Measured Cost of the Alignment Step, and the Gromov Hinge of Buneman's "
    "Converse"
)

DESCRIPTION = """
<p><strong>Not peer reviewed.</strong> This repository once described work as
&ldquo;5-round Opus peer reviewed&rdquo; on the basis of a function that made no
model call; that description was withdrawn and no equivalent claim is made here.</p>

<p>Tests for recombination are normally calibrated against a permutation
ensemble. This paper develops one calibrated against a <em>theorem</em>, reports
the two ways its own first version was wrong, and gives the measurement we
believe is new.</p>

<p><strong>The null is a theorem.</strong> A strict molecular clock on a clonal
genealogy makes expected p-distances ultrametric, and an ultrametric admits no
unfilled Vietoris&ndash;Rips cycle at any scale, so three sequences violating the
strong triangle inequality are a <em>certificate</em> rather than a surprise
relative to shuffles. The inference rule is formalised in Lean&nbsp;4, and emitted
certificates embed the aligned rows so the kernel recomputes every count itself.</p>

<p><strong>That statistic is nonetheless the wrong one.</strong> On a clonal tree
with a relaxed clock and <em>no reticulation whatever</em>, divergence-matched so
only the clock differs, the ultrametric margin exceeds every single-recombinant
margin; at 32 leaves min(no reticulation) = 99/4000 &gt; 79/4000 = max(one
recombinant) at every seed. The four-point condition encodes <em>tree</em> without
<em>clock</em>; we prove in Lean that every symmetric ultrametric satisfies it, and
the quartet margin on the same data does not move.</p>

<p><strong>The main empirical contribution is about the alignment step.</strong>
Running the test on real data required writing an aligner, which puts a step on
the trust surface that no downstream proof covers. Measured on simulated
populations with no reticulation and a known true alignment, reference-anchored
alignment manufactures 0 of 9100 four-point violations at a 0.4% gap load and
1600 of 9100 at 8.2%. We found no published benchmark that varies alignment
quality as an independent variable, so we believe this is the first direct
measurement of a mechanism the field already asserts &mdash; confirmatory but
novel. Read as an admissibility threshold it refuses 7 of 10 real viral
populations; the 3 that pass agree with a literature prior fixed in advance, and
the two certificates clearing their matched floor are kernel-accepted with both
corrupted twins rejected.</p>

<p><strong>Two corrections to our own apparatus are reported as results.</strong>
The band-edge test by which the aligner had been certified is not a soundness
certificate and had really failed &mdash; a pair scoring 16555 against an unbanded
optimum of 17275 with the band edge never touched. Certifying every alignment
against the unbanded optimum dropped the figure above from 4876 to 1600: two
thirds of the headline was an uncertified step inside the measuring apparatus,
inflating the result in the favourable direction. Separately, the reported
statistic is a maximum over subsets where the established treelikeness tests are
averages; run seed-paired against the hard null, the maximum separates 16.57x
against 5.07x for mean delta_q and 1.22x for the proportion of violating
quadruples.</p>

<p><strong>New in edition 2.</strong> The head-to-head edition 1 listed as its
chief limitation has been run, against the one incumbent whose null also assumes
no clock. On the five populations where both tests can be measured <em>they agree
exactly</em>: both specific on the strict-clock and the relaxed-clock null, both
detecting all three reticulate populations at every seed. The incumbent is not
beaten on power, so the paper's claims are revised &mdash; the case rests on an
exhibitable witness and on the absence of a permutation null, not on power. Two
asymmetries are reported, one each way: the incumbent needs no recalibration when
divergence changes, and it returns no measurement at all on gap-blocked data.</p>

<p><strong>A new theorem.</strong> The hinge of Buneman's converse, which this
programme had only cited, is proved: a symmetric dissimilarity satisfies the
four-point condition if and only if its Gromov product satisfies the max
inequality at every root. We found no formalisation of tree metrics, splits or
Buneman's theorem in Lean, Isabelle or Coq. All Lean declarations compile with
zero sorry and no axioms beyond {propext, Classical.choice, Quot.sound}.</p>

<p><strong>Bundle.</strong> The PDF and LaTeX source, the Lean&nbsp;4 modules, and
the raw JSON output and data provenance behind every number, so the claims can be
checked or disputed without trusting this description. Repository:
<a href="https://github.com/xaviercallens/SocrateAI-Scientific-Genetic">SocrateAI-Scientific-Genetic</a>.</p>
""".strip()

KEYWORDS = [
    "reticulate evolution",
    "recombination detection",
    "four-point condition",
    "tree metric",
    "ultrametric",
    "Buneman theorem",
    "Gromov product",
    "Lean 4",
    "formal verification",
    "proof-carrying certificate",
    "sequence alignment",
    "alignment error",
    "molecular clock",
    "phylogenetics",
    "topological data analysis",
    "reproducibility",
]

CREATORS = [
    {
        "name": "Callens, Xavier",
        "affiliation": (
            "SocrateAI Scientific Research & Formal Verification Laboratory"
        ),
    }
]


def token() -> str:
    """Load the Zenodo token through the repository's existing loader."""
    path = REPO / "scripts/three_breakthroughs/publish_3breakthroughs_to_zenodo.py"
    spec = importlib.util.spec_from_file_location("zenodo_token_loader", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load the token helper at {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return str(module.load_zenodo_token())


def req(
    method: str,
    url: str,
    tok: str,
    data: bytes | None = None,
    ctype: str = "application/json",
) -> Dict[str, Any]:
    """One Zenodo API call. Raises on HTTP error; the caller decides what to do."""
    request = urllib.request.Request(
        url,
        data=data,
        method=method,
        headers={"Authorization": f"Bearer {tok}", "Content-Type": ctype},
    )
    with urllib.request.urlopen(request, timeout=900) as response:
        body = response.read()
        return json.loads(body) if body else {}


def md5_of(path: Path) -> str:
    """MD5 of a file, for comparison against what Zenodo reports it stored."""
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bundle() -> List[Path]:
    """The files to deposit, with the checks that must pass before any upload."""
    files: List[Path] = [
        HERE / "papers" / f"{STEM}.pdf",
        HERE / "papers" / f"{STEM}.tex",
        HERE / "formal" / "ANSE" / "GeneticaReticulata.lean",
        HERE / "formal" / "ANSE" / "GeneticaReticulataTree.lean",
        HERE / "formal" / "ANSE" / "ReticulationCertificate.lean",
        HERE / "docs" / "GENETICA_RETICULATA.md",
        HERE / "data" / "reticulata_phase3" / "populations" / "provenance.json",
        HERE / "data" / "reticulata_phase3" / "README.md",
    ]
    files.extend(sorted((HERE / "data" / "reticulata_phase2").glob("*.json")))
    files.extend(sorted((HERE / "data" / "reticulata_phase3").glob("*.json")))
    files.extend(
        sorted((HERE / "data" / "reticulata_phase3" / "emitted").glob("*.lean"))
    )

    missing = [str(f) for f in files if not f.exists()]
    if missing:
        raise RuntimeError(f"bundle is missing files: {missing}")
    forbidden = [str(f) for f in files if "peer_review" in f.name.lower()]
    if forbidden:
        raise RuntimeError(
            f"refusing to deposit a peer-review artefact: {forbidden}"
        )
    names = [f.name for f in files]
    collisions = {n for n in names if names.count(n) > 1}
    if collisions:
        raise RuntimeError(f"duplicate upload names in the bundle: {collisions}")
    return files


def new_version(previous: int, tok: str, files, metadata):
    """Deposit as a NEW VERSION of an existing record, never as a fresh record.

    The project rule, learned the hard way: a corrected edition published as a
    new record leaves the old one as a separate citable artefact with no link to
    its correction. `/actions/newversion` keeps one concept DOI resolving to the
    latest edition.
    """
    draft = req(
        "POST",
        f"{ZENODO}/deposit/depositions/{previous}/actions/newversion",
        tok,
    )
    latest = draft["links"]["latest_draft"]
    deposition = req("GET", latest, tok)
    deposition_id = deposition["id"]
    print(f"new version draft {deposition_id} of record {previous}")

    # A new version inherits the previous version's files; drop them so the
    # bundle is exactly what this edition says it is.
    for old in deposition.get("files", []):
        req(
            "DELETE",
            f"{ZENODO}/deposit/depositions/{deposition_id}/files/{old['id']}",
            tok,
        )
    print(f"  cleared {len(deposition.get('files', []))} inherited files")

    req(
        "PUT",
        f"{ZENODO}/deposit/depositions/{deposition_id}",
        tok,
        json.dumps(metadata).encode(),
    )
    return deposition_id, deposition["links"]["bucket"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="assemble and check the bundle, contact nothing",
    )
    parser.add_argument(
        "--new-version-of",
        type=int,
        default=23574921,
        help=(
            "deposit as a new version of this published record id, which is the "
            "only correct way to publish a corrected edition; 0 creates a fresh "
            "record and should only be used for a first deposit"
        ),
    )
    args = parser.parse_args()

    files = bundle()
    print(f"bundle: {len(files)} files")
    for path in files:
        print(f"  {path.relative_to(HERE)}  {path.stat().st_size:>9,} B")
    if args.dry_run:
        print("\ndry run: nothing was uploaded and no record was created")
        return 0

    tok = token()
    metadata = {
        "metadata": {
            "upload_type": "publication",
            "publication_type": "preprint",
            "title": TITLE,
            "description": DESCRIPTION,
            "creators": CREATORS,
            "keywords": KEYWORDS,
            "license": "cc-by-4.0",
            "access_right": "open",
            "version": "edition-2",
            "language": "eng",
            "related_identifiers": [
                {
                    "identifier": (
                        "https://github.com/xaviercallens/"
                        "SocrateAI-Scientific-Genetic/releases/tag/v1.8.0"
                    ),
                    "relation": "isSupplementTo",
                    "scheme": "url",
                }
            ],
        }
    }

    if args.new_version_of:
        deposition_id, bucket = new_version(args.new_version_of, tok, files, metadata)
    else:
        deposition = req(
            "POST",
            f"{ZENODO}/deposit/depositions",
            tok,
            json.dumps(metadata).encode(),
        )
        deposition_id = deposition["id"]
        bucket = deposition["links"]["bucket"]
        print(f"created deposition {deposition_id}")

    uploaded: Dict[str, str] = {}
    for path in files:
        local = md5_of(path)
        result = req(
            "PUT",
            f"{bucket}/{path.name}",
            tok,
            path.read_bytes(),
            "application/octet-stream",
        )
        remote = str(result.get("checksum", "")).removeprefix("md5:")
        if remote != local:
            raise RuntimeError(
                f"{path.name}: Zenodo stored md5 {remote!r}, local file is "
                f"{local!r}. Not publishing a bundle whose contents do not match."
            )
        uploaded[path.name] = local
        print(f"  uploaded {path.name} (md5 verified)")

    try:
        published = req(
            "POST",
            f"{ZENODO}/deposit/depositions/{deposition_id}/actions/publish",
            tok,
        )
    except urllib.error.HTTPError as error:
        # The API has been seen reporting 5xx for writes that succeeded. Read the
        # state back rather than retrying, which could create a second record.
        print(f"publish returned HTTP {error.code}; reading the state back")
        published = req(
            "GET", f"{ZENODO}/deposit/depositions/{deposition_id}", tok
        )
        if published.get("state") != "done" or not published.get("submitted"):
            raise RuntimeError(
                f"publish failed and the deposition is still "
                f"{published.get('state')!r}; inspect deposition "
                f"{deposition_id} by hand before retrying"
            ) from error
        print("the deposition is published despite the error response")

    receipt = {
        "record_id": published.get("id"),
        "doi": published.get("doi"),
        "concept_doi": published.get("conceptdoi"),
        "concept_record_id": published.get("conceptrecid"),
        "state": published.get("state"),
        "submitted": published.get("submitted"),
        "title": TITLE,
        "version": "edition-2",
        "files_uploaded": uploaded,
        "peer_reviewed": False,
        "note": (
            "First Zenodo deposit of this paper, so a new record. Later editions "
            "must use /actions/newversion on the concept record above, never a "
            "fresh record. This paper has not been peer reviewed."
        ),
    }
    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"\nDOI        : {receipt['doi']}")
    print(f"concept DOI: {receipt['concept_doi']}")
    print(f"receipt    : {RECEIPT.relative_to(HERE)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

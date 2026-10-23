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
RECEIPT = HERE / "data" / "zenodo_paper_g_edition1_receipt.json"
STEM = "PaperG_Split_Level_Certificates_NonTransfer_Lean4"

TITLE = (
    "A Proved Lift That Does Not Transfer: Split-Level Certificates of "
    "Non-Tree-Likeness, and Three Mechanisms of Failure"
)

DESCRIPTION = """
<p>A quartet certificate refutes the four-point condition on four sequences with
an exact integer comparison a proof kernel can recheck. It names four sequences;
what a biologist reads is a <em>split</em>, one group of sequences against
another. This paper promotes the certificate to that level and reports both
halves of the result.</p>

<p><strong>The formal half succeeds.</strong> Because complete deletion gives
every quadruple its own denominator, the Bandelt&ndash;Dress isolation index
&mdash; a <em>minimum</em> over straddling quadruples &mdash; is unavailable
without an incommensurable comparison, so the split certificate is defined as a
<em>logical conjunction</em> of per-quadruple integer comparisons. We then prove
that two certified <em>incompatible</em> splits force a four-point violation,
with the witness quadruple <em>constructed</em> from the four side-intersections
rather than searched for. Four theorems; the Lean module is 27/27 declarations
with no <code>sorry</code> and no axioms beyond {propext, Classical.choice,
Quot.sound}. One of them proves a claim earlier editions asserted in prose.</p>

<p><strong>The empirical half fails, and the failures are the
contribution.</strong> Reaching a working simulated demonstration (20/20
specificity, 15/15 sensitivity over seven populations at five seeds) required
three negative controls to fail first: the criterion had no noise floor and fired
on a two-count margin in two thousand where an exact tree forces equality; no
constant floor can separate, because saturation retains exactly n&minus;3 splits
at a floor 64&times; the one that erases every other dataset; and the matched
null was never matched, because the relaxed-clock simulator's divergence at fixed
mutation rate has a coefficient of variation near 50% spanning
10&ndash;28&times;, so matching it <em>in expectation</em> is ill-posed however
carefully implemented. Rejection sampling on each replicate's <em>realised</em>
divergence repairs it.</p>

<p><strong>The demonstration then does not transfer.</strong> At sixteen taxa on
real viral populations the strict conjunction certifies 0, 2 and 4 splits where a
tree has 13: every population reads as a tree, so its single agreement with a
tree-like prior is vacuous. The mechanism is measured, not inferred &mdash;
conjuncts for a balanced split number C(n/2,2)^2, which is 36 at n=8 but 784 at
n=16, and the closest-missing splits failed by <em>one</em> conjunct out of 91.
An exception budget calibrated against the matched null recovers power and does
not separate: all three populations, including the only tree-like control, stop
reading as trees at a smaller budget than a genuine tree does (&minus;21,
&minus;15, &minus;11). Any threshold between 13 and 24 would classify all three
correctly, which on three populations with one negative control is fitting a
cutoff to the answer, so none is chosen.</p>

<p><strong>The binding constraint is not the statistic</strong> but the alignment
admissibility gate, which leaves three of ten populations usable with a single
negative control &mdash; and one negative control cannot separate &ldquo;this
method false-positives&rdquo; from &ldquo;that population&rsquo;s prior is
wrong&rdquo;. The next step for this line is more admissible negative controls,
not another statistic.</p>

<p>A methodological finding has scope beyond this paper: a matched null
calibrated on the <em>mean</em> of a nuisance parameter is not matched when the
nuisance has a heavy-tailed realised distribution, and averaging more replicates
does not repair it &mdash; averaging shrinks the standard error of the tuning,
not the spread of the replicates the statistic is computed on. The diagnostic is
cheap: measure the across-realisation coefficient of variation of the quantity
being matched.</p>

<p>All figures are produced by committed code and stored artifacts; none is
hand-entered. <strong>Not peer reviewed.</strong></p>
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
    files.extend(sorted((HERE / "data" / "reticulata_phase5").glob("*.json")))
    files.append(HERE / "rust" / "reticulation" / "src" / "splits.rs")
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


def adopt_draft(deposition_id: int, tok: str, metadata):
    """Finish an EXISTING unsubmitted draft instead of creating another one.

    Zenodo's API returns 504 Gateway Time-out on writes that in fact succeeded,
    so `/actions/newversion` can leave a draft behind while the caller sees a
    failure. Re-running the whole script then creates a SECOND draft and the
    first is stranded. This adopts the draft that already exists: clear its
    inherited files, set this edition's metadata, and hand back the bucket so
    the upload and publish proceed normally.
    """
    deposition = req("GET", f"{ZENODO}/deposit/depositions/{deposition_id}", tok)
    if deposition.get("submitted"):
        raise RuntimeError(
            f"deposition {deposition_id} is already submitted; adopting it would "
            "modify a published record"
        )
    inherited = deposition.get("files", [])
    for old in inherited:
        req(
            "DELETE",
            f"{ZENODO}/deposit/depositions/{deposition_id}/files/{old['id']}",
            tok,
        )
    print(f"adopted draft {deposition_id}; cleared {len(inherited)} inherited files")
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
        default=0,
        help=(
            "deposit as a new version of this published record id, which is the "
            "only correct way to publish a corrected edition; 0 creates a fresh "
            "record and should only be used for a first deposit"
        ),
    )
    parser.add_argument(
        "--use-draft",
        type=int,
        default=0,
        help=(
            "finish this existing unsubmitted draft instead of creating a new "
            "version; use after a 504 that left a draft behind, so a second "
            "draft is not created and the first stranded"
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
            "version": "edition-1",
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

    if args.use_draft:
        deposition_id, bucket = adopt_draft(args.use_draft, tok, metadata)
    elif args.new_version_of:
        deposition_id, bucket = new_version(
            args.new_version_of, tok, files, metadata
        )
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
        "version": "edition-7",
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

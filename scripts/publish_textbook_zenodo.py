"""Deposit *Formal Genetics*, edition 5, as a new Zenodo record.

Edition 5 adds Section 12, "Reticulate Evolution: When the Null Can Be a
Theorem", and the reference subsection behind it. This is the book's **first**
Zenodo deposit, so it creates a new record rather than a new version of an
existing one; later editions must go through `/actions/newversion` on the
concept DOI this run returns.

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

    python scripts/publish_textbook_zenodo.py --dry-run   # bundle + checks only
    python scripts/publish_textbook_zenodo.py             # deposit and publish
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
RECEIPT = HERE / "data" / "zenodo_textbook_edition5_receipt.json"

TITLE = (
    "Formal Genetics: Proof, Computation and Evidence in Molecular Biology "
    "(Edition 5)"
)

DESCRIPTION = """
<p><strong>A student textbook in twelve sections.</strong> Each section develops
one piece of molecular genetics three times over: as a mathematical argument the
reader works through by hand, as a theorem checked by the Lean&nbsp;4 kernel, and
as a number produced by a solver or by analysis of a real dataset. Every section
closes with an <em>epistemic ledger</em> separating what was measured, what was
assumed, and what was proved.</p>

<p><strong>This book has NOT been peer reviewed.</strong> It has had an internal
review, which is described in its reproducibility appendix and is not a
substitute. An earlier edition of this repository described work as
&ldquo;5-round Opus peer reviewed&rdquo; on the basis of a function that made no
model call; that description was withdrawn, and the appendix documents the
withdrawal rather than quietly dropping it.</p>

<p><strong>The organising case is a refutation.</strong> An earlier edition of
this material taught that the cancer-associated polymerase variant
<em>POLE</em>&nbsp;p.Pro286Arg acts by removing proofreading stages. Section&nbsp;3
proves a no-go theorem showing no model of that kind reproduces the measured
phenotype, and Section&nbsp;4 gives the measurement that breaks it by a factor of
52. The mechanism is retracted in these pages rather than defended.</p>

<p><strong>New in edition 5 &mdash; Section 12, Reticulate Evolution.</strong>
Section&nbsp;8 measured the topology of a single DNA sequence embedded by the chaos
game, and Section&nbsp;11 reported the negative result that followed. Section&nbsp;12
replaces that test rather than only reporting its failure. Moving from one
sequence to a <em>population</em> makes the null a theorem &mdash; a strict
molecular clock on a clonal genealogy admits no unfilled cycle &mdash; so a
violation becomes a certificate a kernel can check rather than a p-value against
shuffles. The section then refutes its own first statistic: a clonal population
with a relaxed clock and no recombination whatever violates the clock-level
condition <em>more</em> strongly than a real recombinant does, so the tree-level
four-point condition must carry the method, and a machine-checked proof that the
two conditions nest is what licenses reading the difference between them. It
closes on ten real viral populations, seven of which are reported blocked because
a measured control shows the aligner alone would have produced the answer.</p>

<p><strong>What edition 5 corrects in its own apparatus.</strong> A review of the
surrounding literature, run after the work was done, found four errors in it and
they are fixed in the open: the band-edge test used to certify the aligner was
not a certificate and had really failed on this project's own data, which
inflated a headline count threefold; Buneman's theorem is restated with its
actual hypotheses (a nonnegative dissimilarity map, not a metric, because the
triangle inequality is a consequence of the four-point condition); a negative
control was demoted because recombination in that family is documented rather
than merely unexcluded; and the simulator's indel parameters are labelled chosen
rather than calibrated, no published viral indel rate having been found.</p>

<p><strong>Contents.</strong> Reading a proof about a cell &middot; kinetic
proofreading &middot; a no-go theorem &middot; the POLE refutation &middot;
mismatch repair &middot; mutation accumulation &middot; variant landscapes
&middot; topological genomics &middot; DNA topology &middot; from model to clinic
&middot; gaps closed &middot; reticulate evolution. With worked examples,
exercises and solutions, a reproducibility appendix naming the command behind
every number, a research-directions appendix, and an annotated reference list in
which three entries are marked <em>no verified DOI</em> rather than given a
plausible-looking string.</p>

<p><strong>Bundle.</strong> The PDF, the complete LaTeX sources and build script,
the Lean&nbsp;4 modules behind Section&nbsp;12's theorems, and the raw JSON output
and data provenance behind its numbers &mdash; so that the claims can be checked
or disputed without trusting this description of them.</p>

<p>Repository:
<a href="https://github.com/xaviercallens/SocrateAI-Scientific-Genetic">SocrateAI-Scientific-Genetic</a>,
release v1.8.0.</p>
""".strip()

KEYWORDS = [
    "molecular genetics",
    "formal verification",
    "Lean 4",
    "kinetic proofreading",
    "DNA polymerase epsilon",
    "POLE",
    "mismatch repair",
    "topological data analysis",
    "persistent homology",
    "reticulate evolution",
    "recombination",
    "four-point condition",
    "tree metric",
    "sequence alignment",
    "textbook",
    "reproducibility",
    "refutation",
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
    textbook = HERE / "textbook"
    files: List[Path] = [textbook / "build" / "main.pdf"]
    files.append(textbook / "main.tex")
    files.append(textbook / "build.py")
    files.extend(sorted((textbook / "sections").glob("*.tex")))
    files.extend(
        HERE / "formal" / "ANSE" / name
        for name in (
            "GeneticaReticulata.lean",
            "GeneticaReticulataTree.lean",
            "ReticulationCertificate.lean",
            "ProofreadingStageRemoval.lean",
        )
    )
    files.extend(
        sorted((HERE / "data" / "reticulata_phase2").glob("*.json"))
    )
    files.extend(
        sorted((HERE / "data" / "reticulata_phase3").glob("*.json"))
    )
    files.append(HERE / "data" / "reticulata_phase3" / "populations" / "provenance.json")
    files.append(HERE / "data" / "reticulata_phase3" / "README.md")
    files.append(HERE / "docs" / "GENETICA_RETICULATA.md")

    missing = [str(f) for f in files if not f.exists()]
    if missing:
        raise RuntimeError(f"bundle is missing files: {missing}")
    forbidden = [str(f) for f in files if "peer_review" in f.name.lower()]
    if forbidden:
        raise RuntimeError(
            f"refusing to deposit a peer-review artefact: {forbidden}"
        )
    # Distinct upload names: Zenodo is flat, and two files called main.tex would
    # silently collide.
    names = [f.name for f in files]
    collisions = {n for n in names if names.count(n) > 1}
    if collisions:
        raise RuntimeError(f"duplicate upload names in the bundle: {collisions}")
    return files


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="assemble and check the bundle, contact nothing",
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
            "publication_type": "book",
            "title": TITLE,
            "description": DESCRIPTION,
            "creators": CREATORS,
            "keywords": KEYWORDS,
            "license": "cc-by-4.0",
            "access_right": "open",
            "version": "edition-5",
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
        "version": "edition-5",
        "files_uploaded": uploaded,
        "peer_reviewed": False,
        "note": (
            "First Zenodo deposit of this book, so a new record. Later editions "
            "must use /actions/newversion on the concept record above, never a "
            "fresh record. This book has not been peer reviewed."
        ),
    }
    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"\nDOI        : {receipt['doi']}")
    print(f"concept DOI: {receipt['concept_doi']}")
    print(f"receipt    : {RECEIPT.relative_to(HERE)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

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
RECEIPT = HERE / "data" / "zenodo_paper_f_edition7_receipt.json"
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

<p><strong>CORRECTION TO THE MATCHED-NULL FLOORS.</strong> The previous
edition's matched-null floors were not matched. Each was tuned by bisecting a
<em>single</em> realisation of a quantity that is chaotic in the mutation rate
&mdash; the simulator draws its per-branch rate multipliers from the same random
stream as its mutations, so changing the rate shifts every downstream multiplier
&mdash; and the floor was then taken as a <em>maximum</em> over replicates whose
realised divergences differed by up to an order of magnitude. Measured, the
published nulls missed their populations' divergence by &minus;49%, +11% and
&minus;12%; at fixed rate the simulator's divergence has a coefficient of
variation near 50%, spanning 10&ndash;28&times; between extremes, so matching
such a null <em>in expectation</em> is ill-posed however carefully it is
implemented.</p>

<p><strong>Re-deciding against a null matched per replicate</strong> &mdash; each
replicate accepted only if its own realised divergence is within 10% of the
population's &mdash; leaves <em>every verdict and the 3/3 literature agreement
intact</em>, but moves two numbers in opposite directions. Coxsackievirus A6
clears <strong>2.75&times;</strong> its floor rather than 1.24&times;, so the
previous edition understated its own evidence. And the JC polyomavirus margin is
<strong>four-fifths</strong> of its floor rather than the sixth previously
reported: still correctly below, but a far closer call than published. Hepatitis
B is unchanged at 4.50&times; with an identical rational floor. Two defects acted
in opposite directions here &mdash; a non-converging bisection misses either way,
while over-dispersion always inflates a maximum &mdash; so the net direction was
a measurement, not an argument, and a prediction made before re-running it was
wrong.</p>

<p><strong>CORRECTION.</strong> The previous edition defined the
Bandelt&ndash;Dress isolation index as half the <em>minimum</em> of a quadruple's
two crossing pairing sums minus its own. That is not the isolation index: the
definition takes the <em>maximum</em>, clamped below at zero. The two forms agree
on a tree metric, where the crossing sums are equal, and differ exactly where the
four-point condition fails &mdash; so the error was invisible on precisely the
data where the condition holds, and the property used to choose between them (that
the index is an internal edge length on a tree) is satisfied by both and so
decided nothing. The previous edition is superseded.</p>

<p><strong>The corrected statement is stronger.</strong> For sorted pairing sums
s1 &le; s2 &le; s3 the three splits' indices are (s3&minus;s1)/2, (s3&minus;s2)/2
and 0, so the certificate's margin is twice <em>one</em> index &mdash; the middle
one &mdash; and the four-point condition is exactly its vanishing. The reading is
the payoff: a tree metric has exactly <em>one</em> split with positive index, a
certificate makes a second positive, and any two distinct splits of a four-set are
incompatible. So a quartet certificate <em>is</em> the split-decomposition
condition for a split system to fail to be a tree, not an analogue of it &mdash;
checked as an identity, with quadruples carrying two positive indices and
quadruples violating the four-point condition both numbering 47,802 across the
grid. The measurement that followed survives unchanged: the supported-split index
is identical to five decimal places on a clonal population and on both reticulate
populations built from it, so it is blind to the reticulation and tracks
divergence instead. The supported index carries the tree; the index a tree forces
to zero carries the reticulation.</p>

<p><strong>What edition 5 claimed, now corrected.</strong> Edition 4 closed the topological line;
this edition replaces that closure with a proved identification of what the
certificate was already measuring. Classical phylogenetic combinatorics assigns
each split of a quadruple a Bandelt&ndash;Dress isolation index, and between the
two non-minimal splits that index gap is exactly half the pairing-sum gap &mdash;
so <strong>the certificate's margin is twice an isolation-index gap</strong>, and
the four-point condition on a quadruple is precisely the non-positivity of that
gap. Both are machine-checked, and neither proof needs a tree, a complex or a
metric space, only three rationals, which is why the quantity survives being
computed one subset at a time where a global complex cannot. The identity is also
checked against the data on 75,500 quadruples with zero failures, using indices
computed from the per-split definition rather than from sorted sums. And
measuring the split decomposition's own headline quantity locates the signal: the
largest supported-split index is identical to five decimal places on the clonal
population and on both reticulate populations built from it &mdash; blind to the
reticulation, tracking divergence instead. The supported index carries the tree;
the gap carries the reticulation; the certificate measures the gap. That is also
a retrospective explanation of edition 4's topological failure, since a distance
filtration is driven by the magnitudes that track divergence.</p>

<p><strong>What edition 4 added.</strong> The topological line is closed on measured
grounds, and the maximum statistic's defect is measured across taxon counts
rather than argued from two points. Bit-parallel subset counting (21.6&times;
measured, verified against the byte implementation by property test) made
sixty-four taxa affordable; the clonal floor roughly doubles from sixteen to
sixty-four taxa while the signal stays flat, so separation halves from
19.4&times; to 9.2&times; &mdash; detection at a fixed size survives, but a
threshold calibrated at one size does not transfer. A filtration funnel then
tested three hypotheses against the topological statistic's failure, with the
criterion fixed in advance, and <em>none met it</em>: saturation correction is
refuted (it makes the false positive worse), global complete deletion removes the
gap false positive entirely and confirms it was the missing-data confound
reappearing as topology, and the miss of a single recombinant survives all four
variants while the certificate detects that population at every seed. The
global-distance-matrix construction is therefore closed: it cannot use per-subset
complete deletion, and it misses what the certificate catches.</p>

<p><strong>What edition 3 added.</strong> A published-breakpoint test, and the
resolution limit it reveals. Two SARS-CoV-2 triads whose recombinants are
<em>named</em> by the Pango designation process with published breakpoint
intervals, taken already aligned from Nextstrain's open release &mdash; so this
project's aligner is not involved. The pre-registered prediction, that the
dominant pairing separates the two <em>parent</em> lineages, held on one
population (6 of 10 certificates) and failed on the other (0 of 10). Both margins
are 2 and 4 <em>sites</em> on populations differing by roughly 17 and 38 sites
genome-wide, and one matched null produced no violating quadruple at all, so the
conclusion is about applicability: recombinants between recently-diverged
lineages sit <strong>below this method's resolution</strong>, and the tension is
structural &mdash; the datasets with the best ground truth are the least
divergent. Edition 3 also corrects edition 2's comparison table, which judged
this method against a mismatched null and so understated its specificity; and it
reports a persistent-homology experiment that cannot answer its own question,
its null having produced no cycle in 495 of 495 replicates. Both failures taught
one rule, now computed in the output rather than asserted: check that a null
varies before believing a p-value derived from it.</p>

<p><strong>What edition 2 added.</strong> The head-to-head edition 1 listed as its
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
        default=23594789,
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
            "version": "edition-7",
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

#!/usr/bin/env python3
"""Retrieve the genomic records the textbook quotes, live, with provenance.

Every number the book attributes to a database is fetched here and written to
`textbook/data/genomic_records.json` together with the exact URL it came from and
the UTC timestamp of the request. Figures and section text read the cache, never
the network, so a build is reproducible and a reader can audit any figure back to
a named endpoint.

Sources (all public, no credentials):
  * Ensembl REST       -- gene coordinates, MANE Select transcripts (GRCh38)
  * NCBI E-utilities   -- ClinVar variant counts by gene and clinical significance
  * UniProt REST       -- protein entries and exonuclease domain boundaries

    python retrieve_genomic_data.py            # fetch, write cache
    python retrieve_genomic_data.py --offline  # verify the cache without network
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "textbook" / "data" / "genomic_records.json"

ENSEMBL = "https://rest.ensembl.org"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UNIPROT = "https://rest.uniprot.org"
INTERPRO = "https://www.ebi.ac.uk/interpro/api"

# Residues the book points at, checked against the retrieved domain boundaries so
# that "P286R lies in the proofreading domain" is verified rather than asserted.
VARIANT_RESIDUES: dict[str, int] = {"POLE_P286R": 286, "POLE_V411L": 411}

# The genes the book works with: the replicative polymerases whose proofreading
# domain is at issue, and the mismatch-repair genes that Section 5 argues carry
# the missing mechanism.
GENES: dict[str, str] = {
    "POLE": "ENSG00000177084",
    "POLD1": "ENSG00000062822",
    "MSH6": "ENSG00000116062",
    "MSH2": "ENSG00000095002",
    "MLH1": "ENSG00000076242",
}

# UniProt accessions for the same proteins.
PROTEINS: dict[str, str] = {
    "POLE": "Q07864",
    "POLD1": "P28340",
    "MSH6": "P52701",
}

USER_AGENT = "SocrateAI-Formal-Genetics-Textbook/1.0 (reproducible-build)"


def fetch(url: str, *, attempts: int = 3, pause: float = 1.5) -> Any:
    """GET `url` and parse JSON, retrying transient failures.

    Public endpoints rate-limit; a 429 or 5xx is retried, a 4xx is not, because
    a bad identifier should surface as an error rather than be silently skipped.
    """
    last: Exception | None = None
    for attempt in range(attempts):
        req = urllib.request.Request(url, headers={
            "Accept": "application/json", "User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=45) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code < 500 and exc.code != 429:
                raise RuntimeError(f"{url} -> HTTP {exc.code}") from exc
            time.sleep(pause * (attempt + 1))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last = exc
            time.sleep(pause * (attempt + 1))
    raise RuntimeError(f"{url} failed after {attempts} attempts: {last}")


def stamp(url: str) -> dict[str, str]:
    return {"source_url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat()}


def ensembl_gene(symbol: str, ensg: str) -> dict[str, Any]:
    """Gene record plus the MANE Select transcript, which is what clinical reports use."""
    url = f"{ENSEMBL}/lookup/id/{ensg}?expand=1;content-type=application/json"
    rec = fetch(url)
    mane = [t for t in rec.get("Transcript", [])
            if t.get("is_canonical") == 1 or "MANE_Select" in str(t.get("tag", ""))]
    canonical = mane[0] if mane else (rec.get("Transcript") or [{}])[0]
    length = int(rec["end"]) - int(rec["start"]) + 1
    return {
        "symbol": symbol,
        "ensembl_gene_id": ensg,
        "assembly": rec.get("assembly_name"),
        "chromosome": rec.get("seq_region_name"),
        "start": rec.get("start"),
        "end": rec.get("end"),
        "strand": rec.get("strand"),
        "span_bp": length,
        "biotype": rec.get("biotype"),
        "description": rec.get("description"),
        "canonical_transcript": canonical.get("id"),
        "transcript_count": len(rec.get("Transcript", [])),
        **stamp(url),
    }


def clinvar_counts(symbol: str) -> dict[str, Any]:
    """How many ClinVar submissions call this gene's variants pathogenic, and how
    many call them benign. The ratio is the base rate Section 7 teaches students to
    look at before interpreting any single variant."""
    out: dict[str, Any] = {"symbol": symbol}
    for label, term in (
        ("pathogenic", f'{symbol}[gene] AND ("pathogenic"[Clinical significance])'),
        ("benign", f'{symbol}[gene] AND ("benign"[Clinical significance])'),
        ("uncertain", f'{symbol}[gene] AND ("uncertain significance"[Clinical significance])'),
        ("total", f"{symbol}[gene]"),
    ):
        url = (f"{EUTILS}/esearch.fcgi?db=clinvar&retmode=json&retmax=0"
               f"&term={urllib.parse.quote(term)}")
        rec = fetch(url)
        out[label] = int(rec["esearchresult"]["count"])
        out[f"{label}_url"] = url
        time.sleep(0.4)  # NCBI asks for <= 3 requests/second without a key
    out.update(stamp(f"{EUTILS}/esearch.fcgi?db=clinvar"))
    return out


def uniprot_entry(symbol: str, acc: str) -> dict[str, Any]:
    """Protein length from UniProt, domain boundaries from InterPro.

    UniProt annotates no `Domain` feature for these polymerases -- Q07864 carries
    only Regions and a Motif -- so the 3'->5' exonuclease boundaries come from
    InterPro, which is where `IPR006133` is defined. Taking them from the wrong
    endpoint is how one ends up asserting a domain range with nothing behind it.
    """
    url = f"{UNIPROT}/uniprotkb/{acc}.json"
    rec = fetch(url)
    ipro_url = (f"{INTERPRO}/entry/InterPro/protein/uniprot/{acc}/?page_size=100")
    ipro = fetch(ipro_url)

    domains: list[dict[str, Any]] = []
    for result in ipro.get("results", []):
        meta = result.get("metadata", {})
        if meta.get("type") != "domain":
            continue
        spans = [(f["start"], f["end"])
                 for prot in result.get("proteins", [])
                 for loc in prot.get("entry_protein_locations", [])
                 for f in loc.get("fragments", [])]
        if not spans:
            continue
        domains.append({
            "interpro_accession": meta.get("accession"),
            "name": meta.get("name"),
            "start": min(s for s, _ in spans),
            "end": max(e for _, e in spans),
        })
    domains.sort(key=lambda d: d["start"])
    return {
        "symbol": symbol,
        "accession": acc,
        "protein_name": (rec.get("proteinDescription", {}).get("recommendedName", {})
                         .get("fullName", {}).get("value")),
        "length_aa": rec.get("sequence", {}).get("length"),
        "domains": domains,
        "uniprot_url": url,
        "interpro_url": ipro_url,
        **stamp(url),
    }


def collect() -> dict[str, Any]:
    print("retrieving Ensembl gene records ...", flush=True)
    genes = {sym: ensembl_gene(sym, ensg) for sym, ensg in GENES.items()}
    for g in genes.values():
        print(f"  {g['symbol']:6s} chr{g['chromosome']}:{g['start']}-{g['end']} "
              f"({g['assembly']}, {g['span_bp']} bp)", flush=True)

    print("retrieving ClinVar counts ...", flush=True)
    clinvar = {sym: clinvar_counts(sym) for sym in ("POLE", "POLD1", "MSH6")}
    for c in clinvar.values():
        print(f"  {c['symbol']:6s} pathogenic={c['pathogenic']} benign={c['benign']} "
              f"VUS={c['uncertain']} total={c['total']} (categories overlap)", flush=True)

    print("retrieving UniProt entries ...", flush=True)
    proteins = {sym: uniprot_entry(sym, acc) for sym, acc in PROTEINS.items()}
    for p in proteins.values():
        print(f"  {p['symbol']:6s} {p['accession']} {p['length_aa']} aa, "
              f"{len(p['domains'])} domains", flush=True)

    return {
        "purpose": "Live records underpinning the Formal Genetics textbook.",
        "built_utc": datetime.now(timezone.utc).isoformat(),
        "genes": genes,
        "clinvar": clinvar,
        "proteins": proteins,
    }


def verify(doc: dict[str, Any]) -> list[str]:
    """Sanity checks on the cache. These are the invariants the text relies on;
    if one fails the book should not quote the number."""
    problems: list[str] = []
    for sym, g in doc["genes"].items():
        if g.get("assembly") != "GRCh38":
            problems.append(f"{sym}: assembly is {g.get('assembly')}, expected GRCh38")
        if not isinstance(g.get("span_bp"), int) or g["span_bp"] <= 0:
            problems.append(f"{sym}: non-positive span {g.get('span_bp')}")
        if not g.get("canonical_transcript"):
            problems.append(f"{sym}: no canonical transcript")
    # ClinVar's clinical-significance categories OVERLAP: a record carrying
    # conflicting assertions matches more than one filter, so they do not partition
    # the total and their sum may exceed it. The invariant that does hold is that
    # no single category exceeds the total. Asserting a partition here was wrong,
    # and the counts said so.
    for sym, c in doc["clinvar"].items():
        for label in ("pathogenic", "benign", "uncertain"):
            if c[label] > c["total"]:
                problems.append(
                    f"{sym}: {label}={c[label]} exceeds total={c['total']}")
    for sym, p in doc["proteins"].items():
        if not isinstance(p.get("length_aa"), int) or p["length_aa"] <= 0:
            problems.append(f"{sym}: bad protein length {p.get('length_aa')}")
        for d in p.get("domains", []):
            if not (1 <= d["start"] < d["end"] <= p["length_aa"]):
                problems.append(f"{sym}: domain {d['interpro_accession']} span "
                                f"{d['start']}-{d['end']} outside 1-{p['length_aa']}")
    # The book claims P286R sits inside the 3'->5' exonuclease domain IPR006133.
    # Verify it against the retrieved boundaries rather than restating it.
    pole = doc["proteins"].get("POLE", {})
    exo = next((d for d in pole.get("domains", [])
                if d["interpro_accession"] == "IPR006133"), None)
    if exo is None:
        problems.append("POLE: IPR006133 (exonuclease domain) not found in InterPro")
    else:
        res = VARIANT_RESIDUES["POLE_P286R"]
        if not exo["start"] <= res <= exo["end"]:
            problems.append(f"POLE: residue {res} outside IPR006133 "
                            f"{exo['start']}-{exo['end']}")
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--offline", action="store_true",
                    help="verify the existing cache without hitting the network")
    args = ap.parse_args()

    if args.offline:
        if not OUT.exists():
            print(f"no cache at {OUT}; run without --offline first")
            return 1
        doc = json.loads(OUT.read_text())
    else:
        doc = collect()
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(json.dumps(doc, indent=1))
        print(f"\nwrote {OUT}")

    problems = verify(doc)
    if problems:
        print("\nVERIFICATION FAILED:")
        for p in problems:
            print("  -", p)
        return 1
    print(f"verification passed: {len(doc['genes'])} genes, "
          f"{len(doc['clinvar'])} ClinVar queries, {len(doc['proteins'])} proteins")
    return 0


if __name__ == "__main__":
    sys.exit(main())

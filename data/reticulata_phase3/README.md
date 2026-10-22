# Reticulata Phase 3 — real sequence data, and the aligner that had to be built for it

Phase 2 closed the missing-data gap and refuted its own Phase 1 statistic, but it
could not be run on real data at all, for a mundane reason recorded in
[`docs/GENETICA_RETICULATA.md`](../../docs/GENETICA_RETICULATA.md) §10.7: the
test consumes an *alignment*, and this machine has no `mafft`, no `muscle` and no
`clustalo`. Phase 3 removes that block by implementing the alignment step
(`rust/reticulation/src/align.rs`) and then runs the four-point test on real
viral populations — but only on the ones where a measured control says the
aligner is not the thing producing the answer.

Nothing here is hand-edited. The interpretation is in §11 of the same document.

## Regenerating everything

```bash
# 1. the data: 58 046 viral genomes with taxonomy, 50 MB
curl -L -o test.parquet \
  https://huggingface.co/datasets/hiyata/Virus-Host-Genomes/resolve/main/test.parquet

# 2. the populations: pre-registered filters, two selection strategies
python scripts/reticulata_hf_populations.py --parquet test.parquet \
  --out data/reticulata_phase3/populations --taxa 16

# 3. the aligner's control, and the admissibility threshold it measures
cd rust/reticulation && cargo build --release
./target/release/reticulation align-error --length 4000 --deletion-rate 0.3 \
  --mean-deletion 30 > ../../data/reticulata_phase3/align_error.json

# 4. the test on real data, with Lean certificates for what survives
./target/release/reticulation real-data \
  --populations ../../data/reticulata_phase3/populations \
  --emit-lean ../../data/reticulata_phase3/emitted \
  > ../../data/reticulata_phase3/real_data.json
```

Steps 3 and 4 are deterministic given their seeds. Step 2 is deterministic given
the parquet, whose SHA-256 is recorded in `populations/provenance.json`.

| file | what it is |
| :--- | :--- |
| `populations/*.fasta` | Ten populations of 16 real viral genomes: five taxa × two selection strategies (`accession` order, and greedy 6-mer `spread`). Unaligned — these are the aligner's input, not the test's. |
| `populations/provenance.json` | The source digest, the filters, and for every population its accessions, lengths, pre-alignment diversity, and the **literature expectation** about whether its genealogy is tree-like. The expectation was written before any margin was computed. |
| `align_error.json` | The control the Lean kernel cannot supply: simulated clonal, relaxed-clock populations with deletion events, where the true alignment is known. Realign, compare per-quadruple *verdicts*, and count the four-point violations that alignment error alone manufactures. Includes the admissibility ladder. |
| `real_data.json` | The test on the ten real populations, each run twice on two different alignment references, each margin compared against a floor from a simulated tree at matched width, matched divergence and matched gap load. |
| `emitted/*.lean` | **The evidence for the real-data findings.** Each file embeds the four aligned rows it is about; the kernel recomputes every count by `rfl`. The `_bad` twins must fail. |

## Why most of the populations are reported BLOCKED

`align_error.json`'s ladder is the reason, and it is the most important number
Phase 3 produced. On populations with **no reticulation in them at all**,
reference-anchored alignment starts calling quadruples four-point violations as
soon as the indel load rises: at the heaviest rung, 4876 of 9100 quadruples over
five seeds. A margin computed on such an alignment is a statement about the
aligner.

So `real-data` applies a threshold taken from the ladder — the heaviest rung at
which *nothing* was manufactured — and reports anything above it as BLOCKED,
with its counts, rather than interpreting it. Both thresholds are rates per
aligned cell, not counts: an absolute insertion count is not comparable between
a 3 kb and a 12 kb genome.

This is deliberately costly. It blocks populations whose numbers look
interesting. That is the point: the alternative is a published claim about a
virus that is really a claim about a dynamic-programming tie-break.

## Re-checking a real-data certificate without rerunning anything

```bash
cd <repo with formal/ANSE built>
lake env lean data/reticulata_phase3/emitted/reticulata_real_<population>.lean      # accepted
lake env lean data/reticulata_phase3/emitted/reticulata_real_<population>_bad.lean  # must FAIL
```

The kernel computes the six pair-counts and the common column count from the
embedded rows itself, so a transposed pair or a miscopied count cannot survive.
What remains trusted is narrower than in Phase 2 by nothing and wider by one
item: that the embedded rows are what the aligner produced from the FASTA. The
encoding is printed in each file's header (`A=0 C=1 G=2 T=3`, anything else `4`
and not comparable), the FASTA is in `populations/`, and the aligner is 900 lines
with 14 tests — so that step is checkable too, which is why `align_error.json`
exists.

## Kernel verdicts actually obtained

| file | columns | kernel | time |
| :--- | --: | :--: | --: |
| `reticulata_real_hepatitis_b_accession.lean` | 3 182 | accepted, no axioms beyond the whitelist | 67.8 s |
| `reticulata_real_hepatitis_b_accession_bad.lean` | 3 182 | rejected | 70.6 s |
| `reticulata_real_coxsackievirus_a6_spread.lean` | 7 434 | accepted, no axioms beyond the whitelist | 130.5 s |
| `reticulata_real_coxsackievirus_a6_spread_bad.lean` | 7 434 | rejected | 157.5 s |

The first attempt reported all four rejected in 0.3 s for want of a file: `lake`
runs from the formal directory and the emitted paths were relative. Note that
the *rejection* control passed at the same time, both directions failing for the
same irrelevant reason — a control that could not fail. Fixed by canonicalising
the path, and the verdicts above are from the fixed run.

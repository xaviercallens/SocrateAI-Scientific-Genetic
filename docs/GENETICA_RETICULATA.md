# GENETICA RETICULATA

### Certified detection of reticulate evolution

> **The tree is the hypothesis. The loop is the refutation. The triple is the proof.**

---

## 1. The name

**Genetica Reticulata** — "the netted genealogy".

| Why | |
| :--- | :--- |
| **Accurate** | *Reticulate* evolution is the exact technical term for what the method detects: recombination, reassortment, horizontal transfer — genealogy that is a net rather than a tree. |
| **House style** | Matches the project's Latin scholarly naming (`Mechanica Fluidorum`, `Physique Cinétique`, `SocrateAI Mathesis`) rather than inventing an acronym. |
| **It names the object, not the technique** | "Persistent homology of genetic distance" names a tool. *Reticulata* names the thing we are looking for, which is the right way round: when the tool changes, the name should not. |
| **Easter egg** | *reticulata* is a real taxonomic epithet (*Poecilia reticulata*, the guppy). A biologist reads it as a species name and smiles before they read the subtitle. |

**Short form:** `Reticulata`.
**Rust crate:** `reticulata`.
**Lean namespace:** `ANSE.GeneticaReticulata`.
**CLI:** `reticulata`.

Runner-up, recorded because it was close: **Elenchus Genetica**. The method *is* a Socratic elenchus — it refutes the clonal hypothesis from the hypothesis's own commitments — and the project already uses `elenchus` for its refutation machinery. Rejected only because that name is taken project-wide and reusing it would blur two distinct things.

---

## 2. What it is, in one page

The clonal hypothesis has a mathematical commitment: if descent is a tree and the clock is strict, the pairwise genetic distances form an **ultrametric**, meaning the two largest distances in every triple are equal.

> **No-go theorem.** If `d` is an ultrametric then `d(·,·) ≤ ε` is transitive at every scale. So every Vietoris–Rips component is a clique, every cycle carries all its chords, and `β₁ = 0` at **every scale simultaneously**.

Therefore three sequences with `d(x,y) ≤ s`, `d(y,z) ≤ s` and `d(x,z) > s` cannot exist under the hypothesis. Such a triple is not *evidence against* the hypothesis in a statistical sense — it is a **certificate** that refutes it outright.

```
         y                         d(x,y) ≤ s        ✓ close
        / \                        d(y,z) ≤ s        ✓ close
       /   \                       d(x,z) >  s       ✗ impossible if ultrametric
      x     z     ← far apart
```

This is the whole idea, and everything below is engineering around it.

---

## 3. Why Lean 4 is the foundation and not the decoration

In the usual arrangement, a program computes a statistic and a paper asserts it. The reader trusts the program. Nobody reads the program.

Here the arrangement is inverted. **Lean defines what counts as evidence; Rust merely searches for it.** The trust boundary is the certificate, and it is three numbers wide.

```
┌─────────────────────────────────────────────────────────────┐
│  LEAN 4  —  the foundation (trusted)                        │
│                                                             │
│   • IsUltrametric            what the null commits to       │
│   • Certificate              what counts as a refutation    │
│   • certificate_refutes_…    the inference rule             │
│   • ultrametric_admits_no_…  the control: nothing findable   │
│                              under the null                 │
└─────────────────────────────────────────────────────────────┘
                        ▲
                        │  emits candidate certificates
                        │  (3 labels + integer site counts)
┌───────────────────────┴─────────────────────────────────────┐
│  RUST  —  the search (UNTRUSTED, and that is the point)     │
│                                                             │
│   • distances, O(n²), rayon                                 │
│   • triple search, O(n³), exact                             │
│   • Vietoris–Rips front-end → lophat → representative cycles│
└─────────────────────────────────────────────────────────────┘
```

> ### ⚠ Corrected after an adversarial audit (2026-10-22)
>
> An earlier version of this section claimed that "a Rust bug cannot produce a false conclusion that survives checking". **A four-lens adversarial audit falsified that, and it was right.** The integer-only certificates are checked *given the counts*: the kernel proves an arithmetic fact about six integers, and nothing ties those integers to an alignment. A mislabelled triple or a transposed pair would still typecheck. The counting was trusted and we had written otherwise.
>
> The fix is implemented, not just documented: `render_with_sequences` embeds the three sequences so the **kernel computes the counts itself** by `rfl`, then instantiates `sequences_refute_from_counts` against them. Measured at 4000 sites: ~50 s, no axioms beyond the whitelist. The accurate claim is now:
>
> | | |
> | :--- | :--- |
> | **Excluded by the kernel** | a false inequality, a wrong count, a transposed pair, a miscomputed margin — everything downstream of the sequences |
> | **Still trusted** | that the embedded lists *are* the dataset's sequences; checkable by diffing against the alignment, far narrower than trusting a search, but not nothing |
> | **Not a kernel matter** | whether non-ultrametricity means recombination (see §4) |

Compare §11.3 of the textbook, where a wrong subsampling stride silently produced a significant result that only a lucky nested-null design caught. Under the sequence-level form that class of failure *is* excluded — but it took an audit to get the claim right, and the first version of this document overstated it in exactly the way the project keeps finding elsewhere.

---

## 4. The keystone: exact rationals, no floating point in the trust path

This is the foundational change the plan is built on, and it has been verified before being proposed.

A p-distance is **(differing sites) / (comparable sites)** — a rational number with an integer numerator and denominator. There is no reason for a float to appear anywhere in the evidence chain. So a certificate should carry **integers**, and Lean should check it over `ℚ`:

```lean
structure QCertificate where
  d_xy_num : ℤ        -- differing sites between x and y
  d_yz_num : ℤ
  d_xz_num : ℤ
  sites    : ℤ        -- comparable sites (the common denominator)
  sites_pos : 0 < sites
  close_xy : (d_xy_num : ℚ) / sites ≤ (max d_xy_num d_yz_num : ℚ) / sites
  close_yz : (d_yz_num : ℚ) / sites ≤ (max d_xy_num d_yz_num : ℚ) / sites
  far_xz   : (max d_xy_num d_yz_num : ℚ) / sites < (d_xz_num : ℚ) / sites
```

**Feasibility confirmed.** A certificate instantiated with real counts from our simulated recombinant — 4000 sites, 69 / 70 / 131 differing — typechecks with every field discharged by `norm_num`, and the margin comes out as the exact rational `61/4000`. A forged inequality (`70/4000 < 69/4000`) is refutable by `norm_num`, so no forged certificate can be constructed at all. Checked against the project's Lean 4.33 toolchain; script retained at
`scratchpad/q_cert_feasibility.lean`.

The consequence is worth stating plainly: **Rust emits a `.lean` file, and the Lean kernel validates each certificate.** Not a Rust `assert`, not a float comparison with an epsilon — the kernel. That is what "Lean 4 as foundation" should mean, and the current prototype does not yet do it (it validates in Rust, over `f64`). Closing that gap is Phase 1.

---

## 5. Implementation plan

Each phase states its deliverable, its Lean obligation, its Rust obligation, an **acceptance criterion that can fail**, and what would make the phase not worth finishing. The last column is the one usually omitted.

### Phase 0 — Foundations · ✅ complete

| | |
| :--- | :--- |
| **Lean** | `ReticulationCertificate.lean`: no-go, inference rule, non-vacuity witness, control. 10/10 clean, no `sorryAx`. |
| **Rust** | `reticulation` crate: distances, exact triple search, VR front-end over lophat with representative cycles, clonal/recombinant simulators. 25 tests. |
| **Accepted on** | Negative control (clonal) defect 0.00450 vs typical distance 0.030; positive control 0.01525 = 3.39× the clonal floor; second-ranked certificate names the true parents exactly. |

### Phase 1 — Move the trust boundary into the kernel · the critical path

Everything the architecture claims depends on this phase, so it comes first.

| | |
| :--- | :--- |
| **Lean** | Generalise `Certificate` from `ℝ` to a `LinearOrderedField`, add `QCertificate` over `ℚ`, and prove `qcertificate_refutes_ultrametric`. Add `margin_exact : margin = (d_xz_num - max …) / sites`. |
| **Rust** | Replace `f64` distances with exact `(num, den)` integer pairs throughout the evidence path. Emit `certificates.lean` — a file of `QCertificate` instances with integer literals. Keep `f64` only for ranking and display, never for a reported inequality. |
| **Tooling** | `reticulata verify` runs `lake env lean certificates.lean` and reports the kernel's verdict. A certificate that the kernel rejects is a **bug report**, logged, never silently dropped. |
| **Accepts if** | Every certificate from the positive control is kernel-accepted; a deliberately corrupted certificate is kernel-**rejected**. Both directions are required — a verifier that accepts everything is worthless. |
| **Abandon if** | Kernel-checking a realistic certificate set (say 10⁴ certificates) is too slow to run in CI. Mitigation already identified: check only the top-`k` by margin and say so. |

### Phase 2 — Real data, where the answer is already known

| | |
| :--- | :--- |
| **Data** | An aligned viral population with **published** recombination structure. HIV-1 circulating recombinant forms (CRFs) are the best target: the recombinants are named, their parental subtypes are documented, and the breakpoints are mapped. Influenza reassortment is the second target. |
| **Rust** | FASTA input already exists; add a sequence-ID → known-recombinant map so recovery can be scored rather than eyeballed. |
| **Accepts if** | Certificate midpoints recover documented CRFs at a rate better than chance, with the comparison pre-registered: state the expected recovery rate and the null before running. |
| **Abandon if** | Recovery is no better than chance *and* the clonal-margin calibration of Phase 3 shows the margins are within the noise floor. Then the method is a correct theorem without empirical purchase, and that is the finding — publishable, and the end of the line. |
| **Risk, named** | HIV is a hard first case: hypermutation, strong rate heterogeneity, and saturation at the p-distance level. A negative result here would be ambiguous between "method fails" and "HIV violates the clock so badly that the ultrametric null is uninformative". Running a **low-divergence** set first (e.g. a single outbreak) is the de-risking move. |

### Phase 3 — The margin's null distribution · the piece that removes the simulation

The honest limitation of the whole approach is that finite sequences are never exactly ultrametric: **539 of 560 clonal triples violate the null** in our control. Every such certificate is valid; the theorem simply proves something weaker than the interesting claim. Only the margin's *magnitude* discriminates, and today that needs simulation.

| | |
| :--- | :--- |
| **Theory** | Characterise the distribution of `d_max − d_mid` for three sequences related by a clock tree, at `L` sites and divergence `θ`. The three pairwise counts are multinomial with dependencies induced by the tree; the margin is a difference of order statistics. A normal approximation with an explicit variance is almost certainly enough, since only a threshold is needed. |
| **Lean** | Modest and worth doing: formalise the *monotonicity* claims the threshold argument rests on, so the calibration's structure is checked even if the distribution is not. |
| **Rust** | Implement the closed-form threshold; keep the simulator as the cross-check and **assert they agree**. |
| **Accepts if** | The closed-form threshold matches simulation within Monte-Carlo error across a grid of `(L, θ, n)`. |
| **Abandon if** | The dependency structure makes a usable closed form intractable. Then simulation-based calibration is the permanent answer and should be documented as such rather than apologised for. |

### Phase 4 — The tree-level no-go, proved rather than cited

| | |
| :--- | :--- |
| **Lean** | Prove that a four-point (tree) metric has acyclic Vietoris–Rips filtrations. Currently cited to Bauer & Roll. Needs 0-hyperbolicity machinery; the tractable route is the finite combinatorial statement (a tree metric's Rips graph at any scale is a disjoint union of *chordal* graphs) rather than the geometric one. |
| **Payoff** | Lifts the four-point statistic from *suggestive* to *proved*, making the two nested nulls symmetric in status. Separates **rate variation** (ultrametric fails, four-point holds) from **reticulation** (both fail) with equal rigour on both. |
| **Accepts if** | Gate-clean, with a non-vacuity witness and a control, like every other module here. |
| **Abandon if** | It balloons past a few hundred lines. The citation is honest and the code already labels the asymmetry; this is an improvement, not a prerequisite. |

### Phase 5 — Stable summaries

| | |
| :--- | :--- |
| **Context** | Raw Betti numbers are unstable summaries — which is exactly why Lesnick and collaborators introduced the **novelty profile**, and why a UC Davis replication found persistent homology "variable in its detection of recombination". We should not ignore a warning from the literature we rely on. |
| **Rust** | Implement the novelty profile alongside the certificate count; compare their variance under reseeding, the way §11.3 compares the overdispersion estimator across seeds. |
| **Accepts if** | The novelty profile's variance under reseeding is materially lower than the raw `β₁` count's. **If it is not, report that** — a negative methodological result is still a result, and this project has published several. |

### Phase 6 — Upstream and publish

| | |
| :--- | :--- |
| **Upstream** | Offer the Vietoris–Rips front-end to `tomchaplin/lophat` as a PR. It fills a gap the maintainer has documented in the README, and review by someone who knows the reduction is worth more than our own tests. |
| **Cross-validate first** | Diff our `H₁` intervals against GUDHI (already cloned at `/home/xavkal/xdev/gudhi-devel`) and Ripser on random clouds. **Any disagreement is our bug.** This must precede the PR, not follow it. |
| **Publish** | The paper's contribution is the *architecture*, not the biology: a null that is a theorem, and a trust boundary three numbers wide. Deposit with the certificates as data, so a referee can re-check the evidence with a kernel rather than rerun the pipeline. |

---

## 6. Module layout

```
formal/ANSE/
  GeneticaReticulata/
    Ultrametric.lean      IsUltrametric, margin, defect ↔ ultrametric
    Certificate.lean      Certificate, QCertificate, the inference rule
    Controls.lean         non-vacuity witness, the null-admits-nothing control
    TreeMetric.lean       four-point condition  (Phase 4)

rust/reticulata/
  src/
    exact.rs              Rational(num, den) over i64/i128  (Phase 1)
    distance.rs           alignment → exact distance matrix
    certificate.rs        triple search, exact defect, ranking
    emit_lean.rs          certificates → QCertificate instances  (Phase 1)
    rips.rs               VR front-end → lophat  (the upstreamable part)
    novelty.rs            novelty profile  (Phase 5)
    calibrate.rs          clonal threshold, closed form + simulator  (Phase 3)
    simulate.rs           clonal / recombinant generators
    bin/reticulata.rs     controls | fasta | verify | calibrate
```

**Dependency direction is one-way and deliberate: Rust depends on the Lean definitions, never the reverse.** If the two disagree about what a certificate is, Lean is right by construction.

---

## 7. What would make the whole programme wrong

Stated because a plan with no failure mode is a wish.

1. **The null is too weak to be interesting.** Already partly true: finite sampling alone refutes ultrametricity, so "not an ultrametric" is a much weaker claim than "recombinant". Phase 3 decides whether the margin threshold rescues it. If it does not, the honest outcome is a correct theorem with no empirical purchase.
2. **Rate heterogeneity swamps reticulation.** Real clocks are not strict. If biological rate variation produces margins comparable to recombination's, the ultrametric null is the wrong null and the four-point null (Phase 4) must carry the method alone — with the tree-level theorem then becoming a prerequisite rather than an improvement.
3. **The certificate is not the bottleneck anyone cared about.** Existing recombination detection (RDP, GARD, 3SEQ) is mature and well-calibrated. If this adds only "and it is machine-checked", that is a contribution to *epistemics* rather than to *detection*, and the paper should say so rather than claim an advance in power it has not demonstrated.

On (3) we have a view: a method whose evidence a referee can verify with a kernel instead of rerunning a pipeline is worth having even at equal statistical power. But that is a claim about scientific practice, not about recombination, and it should be argued as one.

---

## 9. Adversarial audit of the proof of concept (2026-10-22)

Phase 0 and Phase 1 were built, then attacked by a four-lens adversarial
workflow — forgery, floating point on the evidence path, wrong counts, and
controls that cannot fail — with every claimed finding independently re-verified
by a second auditor instructed to default to REFUTED. 20 agents, 10 findings
confirmed. The audit is the most useful thing that has happened to this design,
so it is recorded in full rather than summarised away.

### Confirmed, and what was done

| # | Finding | Severity | Action |
| :-- | :--- | :--- | :--- |
| 1 | **Emitted file was not about the data** — the kernel proved `0 < margin cert0`, a closed fact about six integers; the inference rule and its dataset hypotheses were never instantiated. A mislabelled triple would typecheck. | high | **Fixed.** `render_with_sequences` embeds the sequences; the kernel computes the counts by `rfl` and instantiates `sequences_refute_from_counts`. Claim corrected everywhere. |
| 2 | **Shadowing defeats the checker** — a hand-written file can redefine `CountCertificate` locally with its obligations as `True := trivial`, and `kernel_check` never verified the type came from the import. | minor–moderate | **Fixed.** All emitted names fully qualified (`ANSE.GeneticaReticulata.…`). |
| 3, 10 | **`kernel_check` did not enforce the axiom whitelist** — it grepped for `sorryAx` only, so a file could declare its own axiom and pass. Also a repo-rule violation. | minor–moderate | **Fixed.** The verifier now parses every `depends on axioms:` line and rejects anything outside `{propext, Classical.choice, Quot.sound}`. |
| 4 | **The PASS/FAIL discrimination gate was decided in `f64`** (`ratio > 2.0`). A verdict depending on rounding. | minor | **Fixed.** Gate is now an exact integer cross-product; `f64` survives only in the displayed ratio. |
| 5 | **The overflow test could not fail** — the product of two `i64` always fits in `i128`, so `checked_cmp` never returns `None`. A forced test, which is this project's signature defect. | cosmetic | **Fixed.** Renamed to state the invariant honestly, and added a real test of `sub_exact`, which *can* overflow because it narrows back to `i64`. |
| 6 | **Pairwise deletion of missing data can manufacture a certificate on an exactly ultrametric genealogy.** Different pairs get different comparable-site counts, and that alone can break the ultrametric inequality. | medium | **Documented** as a first-class caveat below. Not fixable by code: it is a property of pairwise deletion. |
| 7 | **RNA silently analysed** — every `U` was dropped as non-`ACGT`, producing confident output on input the model does not fit. | medium | **Fixed.** `validate_alphabet` refuses RNA with an explicit message. |
| 8 | **Protein alignments silently analysed** over whatever A/C/G/T columns they contained. | minor | **Fixed.** Refused, along with alignments more than half gaps. |
| 9 | `CountCertificate` does not bound its integers, so impossible counts (`n > d`) are kernel-accepted. | minor | **Superseded** by #1's fix: in the sequence-level form the counts are computed, so impossible ones cannot arise. The integer-only form keeps the looser type and is now documented as the weaker of the two. |

Also found and fixed, from the forgery lens as an untested suspicion that we then
tested: **sequence labels came from FASTA headers and were interpolated into a
Lean doc comment unescaped**, so a header containing `-/` closes the comment. We
confirmed a label `a-/axiom` reaches the emitter. Full code injection was blocked
only because the FASTA reader keeps the first whitespace token, which is a
fragile accident rather than a defence; labels are now sanitised and a test
asserts hostile headers cannot escape.

### The caveat that cannot be coded away

Finding 6 deserves more than a table row, because it limits the biological
reading rather than the implementation.

Under pairwise deletion, two pairs of sequences are compared over *different
sets of columns*. Two distances computed over different columns are not
commensurable, and their ratio can violate the ultrametric inequality even when
the underlying genealogy is exactly a clock-like tree. So a kernel-accepted
certificate on gappy data may be a true statement about the *distance matrix*
and a false statement about the *genealogy*.

The kernel cannot help here: it faithfully checks the inequality it is given. The
honest mitigations are to report the per-pair denominators alongside every
certificate — which the emitted files do — and to prefer complete-deletion within
the triple when gap density is non-trivial.

> **Resolved in Phase 2 (§10).** Complete deletion within the triple is now
> implemented in Rust (`subset.rs`), proved in Lean
> (`complete_triple_refute_from_counts`), and *measured*: on strictly clonal,
> gap-blocked data, pairwise deletion manufactures **33–52 violating triples out
> of 560** that complete deletion denies, at all five seeds. The kernel itself
> witnesses one such triple in
> `pairwise_deletion_manufactures_what_complete_deletion_denies`. The gap on
> gappy real data is therefore closed; §10 describes what replaced it, and what
> it costs.

### What the auditors tried that did *not* break it

Recorded because a list of findings without this is not an audit result.

- Constructing a certificate with `n > d` that the kernel accepts: possible, but
  it cannot arise from the sequence-level form, where counts are computed.
- Finding floating point in the distance computation, the triple search, or the
  certificate ranking: none. All exact integer cross-products.
- Making `emit_lean::corrupt` produce a still-valid certificate: impossible. The
  swap inverts the strict inequality `far_xy` by construction, so the rejection
  control cannot silently pass.
- Getting a Rust bug to emit the shadowing text: impossible. The emitter writes a
  fixed preamble and integer literals.

---

## 10. Phase 2 — the tree-level null, and the experiment that nearly sank the method

Phase 1 moved the trust boundary into the kernel. It left two things open, both
named by Phase 1's own audit rather than discovered later: the ultrametric null
conflates *tree* with *clock*, and pairwise deletion can manufacture
certificates. Phase 2 addresses both and then runs the experiment that decides
whether the first one matters. It does.

### 10.1 What was built

| | |
| :--- | :--- |
| **Lean** | `GeneticaReticulataTree.lean`, 16/16 clean, no `sorryAx`, no custom axioms. `IsFourPoint`; `fourPoint_of_ultrametric` (**the nulls are nested — proved, not cited**); `quartet_violation_refutes_ultrametric`; complete-deletion counting (`comparable3`, `diff3xy/yz/xz`) and quadruple counting (`comparable4`, six `diff4*`); the two emitted inference rules; and the finding-6 witness. |
| **Rust** | `subset.rs`: `Alignment` with subset-complete counting, the complete-deletion triple search, the exhaustive four-point quartet search, and `compare_deletion_policies`. `simulate.rs`: a divergence-matched relaxed clock, reciprocal mosaics, gap blocks. `emit_lean.rs`: two new renderers and two new corruptions. 70 tests. |
| **CLI** | `reticulation experiment` (the seven-dataset grid), `reticulation saturation` (the length-scaling diagnostic), `reticulation verify-tree` (kernel acceptance and rejection for both new forms). |

### 10.2 The idea: a null that does not assume a clock

`IsUltrametric d` says `d x z ≤ max (d x y) (d y z)` — "clonal genealogy **and**
strict molecular clock". Real clocks are not strict, so a violation is ambiguous.

The **four-point condition** says the same kind of thing one level up, about
quadruples and with sums in place of single distances:

```lean
def IsFourPoint (d : α → α → ℚ) : Prop :=
  ∀ a b c e, d a b + d c e ≤ max (d a c + d b e) (d a e + d b c)
```

A weighted tree always leaves the two largest pairing sums equal, **whatever its
branch lengths**. So the four-point condition encodes *tree* without *clock*. That
a metric satisfying it on every quadruple is realised by a weighted tree is
Buneman's theorem — **cited, not proved**, exactly the status "a strict clock
gives an ultrametric" has for the Phase 1 null.

What *is* proved here is the nesting:

```lean
theorem fourPoint_of_ultrametric (hsym : ∀ x y, d x y = d y x)
    (hu : IsUltrametric d) : IsFourPoint d
```

and therefore its contrapositive, `quartet_violation_refutes_ultrametric`. One
quartet certificate refutes both nulls. The converse fails, and that asymmetry is
the whole experiment: **large triple margins with floor-level quartet margins
means rate variation on a tree, not reticulation.** Without the nesting theorem
that reading would be an assumption; with it, it is an inference.

The proof is a 256-branch case analysis over the four instances of the strong
triangle inequality available in a quadruple. Unlovely, but it is a theorem.

### 10.3 Complete deletion makes the certificate smaller, not bigger

Counting only the columns where *every* member of the subset carries a base gives
all distances in a certificate **one denominator**. The obligations then stop
being cross-products and become comparisons between integers:

```lean
theorem complete_triple_forces_violation (n_xy n_yz n_xz D : ℕ) (hD : 0 < D)
    (far_xy : n_xy < n_xz) (far_yz : n_yz < n_xz) :
    max (ratioN n_xy D) (ratioN n_yz D) < ratioN n_xz D
```

Compare Phase 1's `count_certificate_forces_violation`, which needs
`n_xy * d_xz < n_xz * d_xy` because the three pairs have three denominators. The
trust surface shrinks and the confound disappears together, which is the kind of
coincidence that suggests the formulation is right.

### 10.4 The experiment: seven datasets, five seeds, predictions registered first

`reticulation experiment` builds seven populations at 4000 sites and 16–18
sequences, five seeds each, and decides every comparison by an exact integer
cross-product against **D1 at the same seed**. Seed-paired, because a margin from
one seed is not a result — a lesson this project has had to learn more than once.

The adversarial dataset is **D2: a clonal tree with a relaxed clock**. Its rate
multipliers are normalised to mean 1, so it has the *same expected divergence* as
D1 and only the clock differs. Without that normalisation a relaxed clock is also
a faster clock, and a violation could not be attributed to either cause — the
first version of this experiment made exactly that mistake, and it changed the
answer.

Largest exact margins, `/4000`, five seeds:

| dataset | reticulation? | clock | complete-deletion **triple** margin | **quartet** margin |
| :--- | :--: | :--- | :--- | :--- |
| **D1** clonal | no | strict | 18, 24, 19, 39, 25 | 7, 5, 6, 7, 7 |
| **D2** clonal | **no** | **relaxed** | **106, 105, 72, 64, 15** | **7, 5, 2, 6, 1** |
| **D3** one recombinant | yes | strict | 61, 47, 61, 67, 61 | 51, 44, 63, 62, 52 |
| **D4** one recombinant | yes | relaxed | 106, 106, 72, 64, 31 | 84, 121, 34, 48, 20 |
| **D5** reciprocal mosaics | yes | strict | 61, 57, 61, 67, 61 | **116, 94, 122, 134, 122** |
| **D7** clonal, 10× divergence | no | strict | 74, 50, 42, 58, 79 | 46, 54, 59, 46, 46 |

Read the D2 row twice. **On data with no reticulation whatsoever, the ultrametric
test produces margins up to 106/4000 — larger than every single-recombinant
margin in D3.** It exceeds twice D1's floor at 3 of 5 seeds. The quartet column on
the same data is indistinguishable from the floor, and is sometimes *below* it: 0
of 5 seeds exceed twice the floor.

Worst-seed separation of positive control from its own floor:

| statistic | positive control | weakest seed |
| :--- | :--- | --: |
| ultrametric triple | D3, one recombinant | **1.72×** |
| four-point quartet | D5, reciprocal mosaics | **16.57×** |
| ultrametric triple | D2, *no reticulation* | 0.60× … up to 5.9× |
| four-point quartet | D2, *no reticulation* | **0.14×** |

**The conclusion is not subtle.** The clock-level statistic is confounded by rate
variation badly enough that its certificates cannot be read as evidence of
reticulation. The tree-level statistic, on the same data, is not moved at all,
and separates its positive control by an order of magnitude more. Risk 2 of §7 —
"rate heterogeneity swamps reticulation" — is **realised for the Phase 1 null and
refuted for the Phase 2 null**. Phase 4 was listed as an improvement; it is a
prerequisite, and this is the evidence.

### 10.5 Finding 6, measured

D6 is D1 with a contiguous gap block in 15 of 16 sequences, and the statistic is
the per-triple **verdict**, not the margin:

| seed | violating under both | **pairwise deletion only** | complete deletion only |
| --: | --: | --: | --: |
| 2026 | 521 | **33** | 6 |
| 2027 | 519 | **36** | 4 |
| 2028 | 501 | **52** | 7 |
| 2029 | 500 | **50** | 9 |
| 2030 | 499 | **49** | 11 |

Each number in the middle column is a triple that pairwise deletion calls a
violation — a certificate the kernel would accept — on a genealogy that is
*exactly* a clock-like tree. 6–9% of all triples. The control that gives this
meaning: on D1, which has no missing data, the two policies agree on **all 560
triples at all five seeds**, as they must, since without gaps they are the same
policy.

The price is stated too: complete deletion leaves roughly 1000–1700 of 4000
columns common to a triple at this gap load. That is a real loss of power, and it
is why complete deletion should be the default only when there is missing data to
delete.

### 10.6 Two of my own criteria were wrong, and both corrections are recorded

Not footnotes — the same error class in three guises, and worth naming because it
is the error this project keeps making.

1. **Comparing margins across different denominators.** The first version of the
   finding-6 test compared the largest *margin* under each deletion policy.
   Complete deletion shrinks the denominator, so one differing site is worth
   more, and its margins come out **larger** even where it is strictly more
   conservative. The two numbers are no more commensurable than the two pairwise
   distances finding 6 is about. Replaced by the per-triple verdict comparison.

2. **A threshold calibrated on one seed.** The power gate required the positive
   control to exceed **twice** the clonal floor — v1.6.0's criterion. But the
   floor is a maximum over 560 triples, an extreme-value statistic, and it varies
   from 18/4000 to 39/4000 across five seeds. v1.6.0's headline **3.389× was the
   best of five seeds**; the same comparison holds at only 3 of 5. The gate is now
   the seed-paired strict inequality, and the 2× result is kept as an ungated
   finding rather than dropped, because its failure *is* the result.

3. **Comparing a margin against a floor from a different divergence.** D7's
   quartet margin exceeds twice D1's at all five seeds, which looked like
   saturation breaking additivity. It is not. `reticulation saturation`
   quadruples the alignment length and watches the margin: sampling noise scales
   as `L^(-1/2)` and a deterministic bias does not shrink at all.

   | divergence | margin at L=1000 | 4000 | 16000 | ratio per 4× L |
   | :--- | --: | --: | --: | :--- |
   | mean 0.026 | 0.00380 | 0.00160 | 0.00071 | 0.42, 0.45 |
   | mean 0.223 | 0.02620 | 0.01255 | 0.00665 | **0.48, 0.53** |

   Both shrink like `L^(-1/2)`. There is **no saturation bias** at ten times the
   divergence — the quartet noise *floor* simply grows with divergence, so D7 was
   being compared against a floor from a population ten times less divergent.
   The practical consequence is the opposite of a divergence ceiling: high
   divergence is usable, provided the null is calibrated at **matched divergence
   and matched length**. The divergence sweep in that command's output gives the
   floor as a function of mean p-distance from 0.006 to 0.33.

   This matters for Phase 2's target data. Had the bias been real there would
   have been no fix available: a Jukes–Cantor-corrected distance is a logarithm,
   hence irrational, and the exact-rational certificate cannot represent it.

### 10.7 What is still not established

- **Buneman's theorem is cited.** We prove that a quartet certificate refutes
  `IsFourPoint` and that `IsFourPoint` is implied by `IsUltrametric`. That
  `IsFourPoint` on every quadruple implies a weighted tree realises the metric is
  not proved here. The honest chain is: certificate ⟹ not four-point ⟹ (Buneman)
  not tree-realisable.
- **Four-point certificates still exist under the null**, by sampling noise, as
  §10.4's D1 column shows. The margin's magnitude still carries the inference and
  the calibration is still simulation-based. Phase 3 is unaffected by this work.
- **No real data yet.** *(Closed in Phase 3 — see §11. The aligner was built,
  ten real viral populations were tested, and seven were refused by a measured
  admissibility threshold.)* The gate that blocked it (finding 6) is closed, but this
  machine has no multiple-sequence aligner and no aligned homologous population
  with documented recombination. The genomic-benchmark sets cached locally are
  unaligned, non-homologous promoter and enhancer sequences; running a
  phylogenetic test on them would be precisely the category error Phase 1's audit
  punished. Phase 2 on real data needs an alignment first, and that is a separate
  step with its own failure modes.

### 10.8 Both new forms, kernel-checked in both directions

`reticulation verify-tree` emits four files from the D5 positive control at seed
2026 and runs the project's Lean toolchain on each. The kernel computes every
count from the embedded sequences; nothing is handed to it but the sequences and
the claimed counts.

| file | what it asserts | kernel | time |
| :--- | :--- | :--: | --: |
| `reticulata_complete_triple.lean` | `complete_triple_refute_from_counts`, 69 / 70 / 131 differing over **4000 columns common to all three** | **accepted** | 42.4 s |
| `reticulata_quartet.lean` | `quartet_refute_from_counts`, pairing sums **256** vs 116 and 140 over 4000 common columns | **accepted** | 71.1 s |
| `..._complete_triple_bad.lean` | far and near numerators swapped | **rejected**, `unsolved goals` | — |
| `..._quartet_bad.lean` | dominant pairing exchanged with a dominated one | **rejected**, `unsolved goals` | — |

Two details worth keeping:

- The complete-deletion triple reports **exactly** Phase 1's counts, 69/70/131
  over 4000. It must: D5 has no missing data, so complete deletion within the
  triple *is* pairwise deletion there. The agreement is a cross-check of the new
  counting path against the old one, not a coincidence.
- The quartet certificate's dominant pairing is `clonal_003` against
  `clonal_013` — **the two true parents** — with both planted mosaics in the
  non-dominant pair. The method did not merely fire; it named the right four
  sequences and put them on the right sides. The sums are symmetric
  (128 + 128 against 58 + 58 and 70 + 70), which is the reciprocal-mosaic
  signature of `reciprocal_recombinants_witness` at 4000 sites instead of 4.

Both corruptions are invalid in Rust before they are emitted, and invalid by
construction rather than by luck: swapping the far and near numerators inverts a
strict inequality, and exchanging the dominant pairing makes the new first
conjunct of `fourPointViolates` the old second-largest sum compared against the
old largest.

### 10.9 Where this leaves the plan

| phase | status after this work |
| :--- | :--- |
| **0, 1** | complete, and Phase 1's audit findings all closed or documented |
| **2** (real data) | **unblocked** — finding 6 is closed and measured. Still needs an aligned homologous population; no aligner on this machine. |
| **3** (the margin's null) | unchanged, and now clearly the critical path: §10.4's D1 and §10.6's length scaling are the first two data points of exactly the calibration Phase 3 is supposed to provide. |
| **4** (tree-level no-go) | **promoted from improvement to prerequisite, and half-done.** The nesting is proved and the quartet certificate works; Buneman's direction remains cited. |
| **5, 6** | unchanged |

The recommendation that follows from §10.4 is blunt: **report the quartet
statistic, not the triple statistic, as the method's primary evidence.** Keep the
triple certificate — it is cheaper, it is what the nesting theorem is stated in
terms of, and a triple violation with no quartet violation is itself an
informative signal about rate variation. But a triple certificate alone should
not be published as evidence of reticulation, and v1.6.0's headline overstated
what a 3.389× margin ratio means.

### 10.10 Scale check at 32 taxa, where the result gets worse for the triple test

§10.4 used 16 leaves. Doubling to 32 (`--generations 5`, three seeds) means 4960
triples and **35 960–46 376 quadruples per dataset, searched exhaustively** — the
quartet search refuses to sample, so the ceiling has to be raised deliberately
rather than quietly. Largest exact margins, `/4000`:

| dataset | reticulation? | clock | **triple** margin | **quartet** margin |
| :--- | :--: | :--- | :--- | :--- |
| **D1** clonal | no | strict | 24, 30, 29 | 10, 9, 10 |
| **D2** clonal | **no** | **relaxed** | **155, 157, 99** | **11, 18, 8** |
| **D3** one recombinant | yes | strict | 79, 77, 75 | 79, 70, 69 |
| **D5** reciprocal mosaics | yes | strict | 79, 77, 75 | 158, 154, 126 |

The D2 and D3 triple margins no longer merely overlap. **At every seed, the
clonal population with a relaxed clock has a larger ultrametric margin than the
population containing an actual recombinant** — `min(D2) = 99 > 79 = max(D3)`. A
larger triple margin is, on this data, weak evidence *against* reticulation. The
statistic is not conservative, not noisy, but ordered the wrong way round.

The quartet column on the same populations separates with no overlap at all:
`{8, 11, 18}` with no reticulation against `{69, 70, 79}` for one recombinant and
`{126, 154, 158}` for reciprocal mosaics.

Adding taxa makes the confound worse for the triple test, which is what the
mechanism predicts: more leaves means more chances for a pair of lineages with
extreme rates, and the statistic is a maximum over triples. The four-point floor
grows far more slowly (9–10 against 5–7 at 16 taxa) because a tree metric
constrains quadruples however the rates vary.

The gap-manufacturing rate holds at scale too: 286–366 manufactured triples out
of 4960, the same 6–7% as at 16 taxa, with zero disagreement on gapless data.

Raw output: [`data/reticulata_phase2/experiment_n32.json`](../data/reticulata_phase2/experiment_n32.json).

## 11. Phase 3 — real data, and the aligner that had to be built to get there

§10.7 recorded the one thing Phase 2 could not do: *"No real data yet … this
machine has no multiple-sequence aligner and no aligned homologous population
with documented recombination."* Phase 3 removes both halves of that block. It
ends with three real viral populations tested and seven refused, and with a
number that is more important than any of them: how many four-point violations
reference-anchored alignment manufactures on a genealogy that is exactly a tree.

### 11.1 What was built

| piece | where | size |
| :--- | :--- | :--- |
| Reference-anchored alignment: affine-gap global DP (Needleman–Wunsch–Gotoh), banded, plus the **unbanded optimum** every banded path is certified against (§11.4), centroid reference selection, FASTA | `rust/reticulation/src/align.rs` | 14 unit tests |
| The identity regression test for the gap-pair tie | `rust/reticulation/tests/aligner_identity.rs` | 2 tests |
| Clonal, relaxed-clock descent **with deletion events and its own true alignment** | `simulate::clonal_tree_relaxed_with_deletions` | 2 tests |
| `reticulation align-error` — the aligner's control and the admissibility ladder | `src/main.rs` | — |
| `reticulation real-data` — the test on real populations, two references each, matched floors, Lean emission | `src/main.rs` | — |
| Population extraction from Hugging Face under pre-registered filters | `scripts/reticulata_hf_populations.py` | — |

### 11.2 The data: ten real populations, with expectations fixed in advance

`hiyata/Virus-Host-Genomes` (58 046 viral genomes with taxonomy, MIT licence) is
reachable and carries enough near-complete genomes per species to build
populations. Five taxa were chosen for one reason only: each has an *independent*
literature expectation about whether its genealogy is tree-like, so the test can
be wrong in a visible way.

| taxon | genome | pool after filters | expectation | why |
| :--- | --: | --: | :--- | :--- |
| Hepatitis B virus | 3 182 nt | 862 | **reticulate** | inter-genotype recombinants are recognised taxonomic entities |
| Coxsackievirus A6 | 7 434 nt | 559 | **reticulate** | enterovirus A lineages are named by recombinant form |
| Lyssavirus rabies | 11 924 nt | 326 | tree-like | non-segmented negative-strand RNA; recombination reported only exceptionally |
| Human papillomavirus 52 | 7 937 nt | 161 | tree-like (weak) | co-divergence with the host; within-type recombination is rare but **documented** — a literature audit after the run flagged this as the weakest prior in the table. A clonal bacterial lineage such as *M. tuberculosis* would be the stronger negative control. No verdict depends on it: both HPV52 populations were refused by the gate. |
| JC polyomavirus | 5 121 nt | 89 | tree-like | genotypes track human migration and are analysed as a tree |

Filters, fixed before any margin was computed: within 2% of the taxon's median
length, at most 0.1% non-ACGT, exact duplicates removed, 16 sequences per
population. Two selection strategies are built for every taxon and **both are
reported**: `accession` order (arbitrary with respect to phylogeny) and greedy
6-mer `spread` (deliberately harder — more divergence means more alignment error
and more saturation). Ten populations in all. Provenance, accessions and the
source digest: [`data/reticulata_phase3/populations/provenance.json`](../data/reticulata_phase3/populations/provenance.json).

### 11.3 The alignment is part of the trust surface, and this is how big it is

A certificate proves an arithmetic fact about counts the kernel takes from the
sequences the file embeds. If the alignment put non-homologous residues in the
same column, those counts are computed faithfully from the wrong data and the
kernel accepts them. Real data therefore adds exactly one item to the trust
surface — the aligner — and the only honest way to size it is to align data whose
true alignment is known.

`clonal_tree_relaxed_with_deletions` provides that. Deletions are frame
coordinates set to `-`, inherited by every descendant, so the founder frame is a
genuine alignment rather than bookkeeping. **There is no reticulation in these
populations at all.** Strip the gaps, hand the ragged sequences to the aligner,
realign, and compare per-quadruple *verdicts* — not margins, because the two
alignments have different widths and §10.6 is the record of what comparing
magnitudes across denominators does.

16 taxa, 4000 nt, five seeds, 1820 quadruples each (9100 per rung):

| mean deletion length | inferred gap fraction | worst column recovery | quadruples manufactured | seeds with any |
| --: | --: | --: | --: | :--: |
| 30 | 0.00000 | 1.00000 | 0 | 0/5 |
| 30 | 0.00398 | 0.99997 | 0 | 0/5 |
| 30 | 0.00549 | 0.99880 | 251 | 2/5 |
| 30 | 0.00870 | 0.99712 | 108 | 1/5 |
| 30 | 0.01507 | 0.99906 | 110 | 3/5 |
| 30 | 0.02012 | 0.99831 | 149 | 3/5 |
| 120 | 0.00000 | 1.00000 | 0 | 0/5 |
| 120 | 0.01322 | 0.99901 | 57 | 2/5 |
| 120 | 0.01627 | 0.99901 | 54 | 2/5 |
| 120 | 0.03731 | 0.99467 | 324 | 3/5 |
| 120 | 0.05254 | 0.97501 | 253 | 3/5 |
| 120 | 0.08233 | 0.96032 | **1600** | 5/5 |

At the heaviest rung the aligner invents 1600 of 9100 quadruple verdicts. **Past
a threshold, alignment error is not a correction to this test; it is the
result.** The curve is not monotone (0.00870 manufactures fewer than 0.00549),
which is what an extreme-value statistic over 1820 quadruples at five seeds
looks like — one more reason to gate on *whether anything was manufactured*
rather than on how much.

**These are the certified numbers, and they are much smaller than the first
run's.** Before §11.4's certificate was in place the same ladder reported 4876
manufactured at the heaviest rung with column recovery down to 0.8681. Roughly
two thirds of that was not reference-anchored alignment error at all: it was the
band silently returning sub-optimal paths. With every alignment proved optimal,
column recovery never falls below 0.960 and the worst rung manufactures 1600.
The qualitative conclusion is unchanged — above a threshold the aligner invents
violations — but its magnitude was overstated threefold, and the overstatement
came from our own unchecked heuristic rather than from the method being
measured.

So the threshold is read off the ladder instead of chosen: **gap fraction
≤ 0.00398 and insertions dropped ≤ 8.06 × 10⁻⁴ per aligned cell**, the heaviest
rung at which nothing was manufactured at any seed. Both are rates, not counts:
an absolute insertion count is not comparable between a 3 kb and a 12 kb genome,
and treating it as if it were would be the fourth instance of the incommensurable
comparison this project keeps making.

Three predictions were registered before the run. **A1** (the aligner is the
identity on gapless input, and its optimal score equals the diagonal score) is
the gate and holds. **A2** (worst-seed column recovery > 0.99) holds at 0.99726.
**A3** (nothing manufactured) holds at the registered load and *fails* higher up
the ladder — which is the finding, not a failed run.

The thresholds themselves did **not** move under certification: the first rung
that manufactures anything is the same one, so the gate of §11.6 and its seven
refusals stand on the certified run as they did on the uncertified one.

Raw output: [`data/reticulata_phase3/align_error.json`](../data/reticulata_phase3/align_error.json).

### 11.4 The band-edge test is not a certificate, and an audit of the literature caught it

The aligner's dynamic programme is banded, and §11.1 originally claimed it
"never returns a path the band constrained without saying so" on the strength of
a band-edge test. **That claim was wrong.** Gibrat (BMC Bioinformatics 19:226,
2018) measured both failure directions on real data: **2.2% of alignments were
optimal yet reached the band edge, and 0.5% were sub-optimal without ever
reaching it.** Not touching the edge is evidence of optimality, not proof of it.

The repair needs no theorem. `optimal_global_score` computes the unbanded
optimum by full Gotoh in `O(nm)` time and `O(m)` memory — no band, no heuristic,
no traceback — and `align_pair_certified` widens the band until the banded path
attains it, refusing to return one that does not. At this scale that is
affordable (15 pairs per population, one extra `O(nm)` pass each), and both
`align-error` and `real-data` now use the certified form, so every number in
§11.3 and §11.6 rests on an alignment proved optimal rather than one that passed
a proxy.

**And the proxy had in fact failed, on our own data.** The first certified run
exited non-zero immediately:

```
error: aligning 'indel_000' to reference 'indel_012': the banded path scores
16555 against an unbanded optimum of 17275: the band of half-width 64 changed
the answer, and the band-edge test did not notice.
```

A deficit of 720 score points, with the traceback never touching an edge — so
edge-driven widening never fired, and the pre-certification numbers in this
section were computed partly on alignments the band had quietly degraded. This is
Gibrat's 0.5% class, observed here rather than cited. The widening loop is now
driven by the **score** rather than by the edge test: widen on either signal,
and return only a path that attains the optimum or an error. The band-edge test
survives as a cheap early trigger, with no soundness role.

Two related findings from the same audit, acted on:

- **Gotoh's published initialisation is buggy** (Flouri et al.,
  `doi:10.1101/031500`): it omits one of two gap-opens and can return a score
  that is not the cost of any alignment, and the error is reproduced in standard
  textbooks and in several widely used implementations. Our initialisation turns
  out to be the corrected one, but that was luck until it was checked. It is now
  checked against *exhaustive enumeration* of every alignment of short
  sequences — `brute_force_agrees_with_the_dp_on_small_inputs` — rather than
  against another DP.
- **The centre-star 2-approximation does not apply here.** Gusfield's `2 − 2/k`
  bound requires a *metric* cost; with affine gaps it is void. The general ratio
  is `2 − l/k` (Bafna–Lawler–Pevzner 1997), and the citation for NP-hardness
  under fixed metrics is Elias 2006, not Wang & Jiang 1994, whose scheme is
  explicitly non-metric. Nothing in this repo claimed the bound, and nothing
  should.

Also recorded as an honest limitation rather than repaired: **the deletion
simulator's indel rate and mean length are chosen, not calibrated.** No
published viral indel length distribution or per-site viral indel rate could be
found, so the ladder varies a load rather than matching one to biology.

### 11.5 A tie-break that was nearly a result

A1 failed at 1 of 5 seeds on the first implementation, and the cause was not a
bug in the DP. With the affine gap at `-10` to open and `-1` per residue, three
consecutive mismatches score `3×(-4) + 5 = -7`, and deleting one reference
residue, aligning the rest, then inserting one query residue scores
`-11 + 15 - 11 = -7`. The identity was one of *two* optimal alignments and the
implementation returned the other, because `max_by_key` in Rust returns the last
maximum.

The repair was to state the preference — M over X over Y, a substitution over a
gap pair — and to strengthen A1 so it checks the score as well as the rows. The
lesson generalises past this crate: an aligner has free parameters that no amount
of downstream proof constrains, and a tie-break is one of them. Regression test:
`rust/reticulation/tests/aligner_identity.rs`.

### 11.6 The real-data result

Every population is aligned twice, on the two most central sequences by 6-mer
distance, and the report says whether the verdict moved. Every margin is compared
against a floor from a simulated clonal **relaxed-clock** tree at matched
alignment width, matched mean p-distance (tuned by bisection on `mu`) and matched
gap load — because a margin on its own is not a statistic, and a floor from the
wrong divergence measures the floor.

**Seven of ten populations are reported BLOCKED** by the §11.3 gate, with their
counts, rather than interpreted:

| population | gap fraction | insertions / cell | why blocked |
| :--- | --: | --: | :--- |
| `coxsackievirus_a6_accession` | 0.0025 | 3.5 × 10⁻³ | insertion load |
| `hepatitis_b_spread` | 0.0588 | 5.8 × 10⁻² | both, by an order of magnitude |
| `hpv52_accession` | 0.0005 | 1.7 × 10⁻³ | insertion load |
| `hpv52_spread` | 0.0004 | 3.6 × 10⁻³ | insertion load |
| `jc_polyomavirus_spread` | 0.0624 | 6.2 × 10⁻² | both |
| `rabies_accession` | 0.0013 | 1.4 × 10⁻³ | insertion load |
| `rabies_spread` | 0.0020 | 1.8 × 10⁻³ | insertion load |

This is expensive on purpose. `hepatitis_b_spread` has a margin of 78/1737 and
`jc_polyomavirus_spread` 141/3363; both are discarded. The alternative is a
published claim about a virus that is really a claim about a dynamic-programming
tie-break. Note also which strategy survives: `spread` maximises diversity, and
diversity is exactly what breaks reference-anchored alignment.

The three admissible populations:

| population | expectation | mean p | quartet margin | matched floor | floor's `mu` | verdict | references agree |
| :--- | :--- | --: | --: | --: | --: | :--- | :--: |
| `hepatitis_b_accession` | **reticulate** | 0.0088 | **9/3182** | 2/3182 | 0.00121 | exceeds **2×** the floor | yes |
| `coxsackievirus_a6_spread` | **reticulate** | 0.1096 | 329/7312 | 268/7390 | 0.03256 | exceeds the floor, not 2× | yes |
| `jc_polyomavirus_accession` | tree-like | 0.0096 | 4/5108 | 24/5121 | 0.00588 | **below** the floor | yes |

Three populations, three agreements with a prior fixed in advance, and the two
alignment references agree on the verdict in all three. The JC polyomavirus row
is the one that could most easily have embarrassed the method and did not: its
margin is a sixth of the floor a tree produces at the same width and divergence.

Two cautions that belong next to that table and not in a footnote:

- **Three is not a validation.** It is the first time this test has run on data
  that was not simulated, with the aligner's contribution bounded by a
  measurement instead of assumed away.
- **The Coxsackievirus row is weak and at high divergence** (mean p = 0.11). It
  clears 1× and not 2×, and §10.6's length-scaling diagnostic says the floor
  itself grows with divergence. Read it as "not inconsistent with the prior",
  not as a detection.

Both populations that clear their floor carry a Lean file embedding the four
aligned rows, checked in both directions by the project's toolchain:

| file | columns | kernel | time |
| :--- | --: | :--: | --: |
| `reticulata_real_hepatitis_b_accession.lean` | 3 182 | **accepted**, no axioms beyond the whitelist | 52.7 s |
| `reticulata_real_hepatitis_b_accession_bad.lean` | 3 182 | **rejected** | 59.5 s |
| `reticulata_real_coxsackievirus_a6_spread.lean` | 7 434 | **accepted**, no axioms beyond the whitelist | 117.6 s |
| `reticulata_real_coxsackievirus_a6_spread_bad.lean` | 7 434 | **rejected** | 135.0 s |

The kernel recomputes all six pair counts and the common column count from the
embedded rows by `rfl`, so these are the first *real-data* findings in this
programme that a reader can re-verify without running any of our code. Two
process notes belong here, both instances of the same failure shape. The first
attempt reported all four files rejected in 0.3 s for want of a file —
`kernel_check` runs `lake` from the formal directory and the paths were
relative — and the rejection control came out vacuously true at the same time,
both directions failing for the same irrelevant reason. Fixed by canonicalising
the path. Then a later run reported the Hepatitis B file rejected in 4.4 s while
the Coxsackievirus file was accepted in 132.8 s: the cause was a concurrent
`lake build` of `GeneticaReticulataTree` in another shell, which left the import
briefly broken. Run directly, that file exits 0 with only whitelisted axioms.
The verdicts tabulated above are from a clean run with nothing else touching
`lake`, and the lesson is the one already in this document twice: **a fast
rejection is a suspicious rejection**, because the cheapest way to fail is not
to have run at all.

Raw output: [`data/reticulata_phase3/real_data.json`](../data/reticulata_phase3/real_data.json),
certificates in [`data/reticulata_phase3/emitted/`](../data/reticulata_phase3/emitted/).

### 11.7 Where this leaves the plan

| phase | status after Phase 3 |
| :--- | :--- |
| **0, 1** | complete |
| **2** (real data) | **done, with a gate.** Ten populations built, three tested, seven refused by a measured admissibility threshold. The missing-aligner block is gone. |
| **3** (the margin's null) | still the critical path, and now sharper: the matched-floor machinery of §11.6 is a simulation-based stand-in for exactly the distribution Phase 3 owes. |
| **4** (tree-level no-go) | unchanged: nesting proved, Buneman cited. |
| **5** (stable summaries) | unchanged. |
| **6** (upstream, publish) | unchanged. |
| **new** | **a progressive or consistency-based aligner.** Reference-anchored alignment is what blocked seven of ten populations. A literature audit (§11.8) softens this: reference-anchored projection *is* standard surveillance practice and is the final step of MUSCLE5's Super5 mode, and its one real advantage is that without a guide tree there is no error-propagation failure mode. It remains the binding constraint on coverage, but the cheaper route is to use data that is **already aligned**. |

### 11.8 What a literature audit changed, and the directions it opened

A seven-thread review of the surrounding literature was run after the work above
was committed. It corrected four things in it — the band-edge claim (§11.4), the
hypotheses of Buneman's theorem, the HPV52 prior (§11.2), and the status of the
simulator's indel parameters — and produced five directions that are better than
anything in the original plan:

1. **A real positive control that bypasses our aligner entirely.** The
   SARS-CoV-2 Pango designation list carries **1 107 recombinant `X*` lineages,
   111 of them with explicit numeric breakpoint intervals** in reference
   coordinates, each traceable to a reviewed GitHub issue; and Nextstrain
   publishes an **open, non-GISAID, already-aligned** `aligned.fasta.zst`
   (Nextclade-aligned to MN908947), so the breakpoints and the alignment share
   one coordinate system. That removes the aligner from the trust surface and
   the admissibility gate with it. Licence: pango-designation is CC BY-NC 4.0,
   non-commercial.
2. **The head-to-head has now been run, and the incumbent is not beaten.**
   `incompat::phi_test` implements a restricted pairwise-homoplasy index — the
   four-gamete test on biallelic, parsimony-informative, gap-free columns, with
   the site-order permutation null of Bruen, Philippe and Bryant. Both tests run
   on the same alignments under the same complete-deletion policy. Because Φw is
   a count of compatible neighbouring pairs over a count of pairs, and the pair
   count is invariant under permutation, the permutation test compares integers:
   exact, like the certificate.

   | dataset | reticulation | Φw at 0.05 | ours vs **D1's** floor | ours vs its **own matched** floor | calibrated mean δq |
   | :--- | :--: | :--: | :--: | :--: | :--: |
   | D1 clonal, strict clock | no | 0/5 | 0/5 | **0/5** | 0/5 |
   | D2 clonal, **relaxed clock** | no | **0/5** | 0/5 | **0/5** | 0/5 |
   | D3 one recombinant | yes | 5/5 | 5/5 | **5/5** | 5/5 |
   | D4 relaxed clock + one recombinant | yes | 5/5 | 5/5 | **5/5** | 5/5 |
   | D5 reciprocal mosaics | yes | 5/5 | 5/5 | **5/5** | 5/5 |
   | D6 gap blocks | no | **no measurement** | 5/5 | **0/5** | 0/5 |
   | D7 saturated | no | **0/5** | 5/5 | **0/5** | 0/5 |

   **The fourth column was an unfair comparison and the fifth replaces it.** An
   earlier version of this table reported only "ours vs D1's floor" and so
   recorded the certificate as firing 5/5 on the two no-reticulation datasets D6
   and D7. That was the incommensurability recorded three times in this document,
   committed a fourth time against ourselves: D6 has gaps, so a different
   denominator, and D7 has ten times the divergence, so a different noise floor.
   Given a floor matched to its *own* width, divergence and gap load — the same
   construction `real-data` uses, and the same the calibrated averaged statistic
   gets — **the certificate is specific on all four no-reticulation datasets and
   detects all three reticulate ones at every seed.** Nothing about the method
   changed; the comparison was corrected.

   A third detector is in the table for the same reason. Mean δq, calibrated as a
   *test* by parametric bootstrap under a matched tree null rather than reported
   as a score, is the fairest possible competitor: same quadruples, same
   alignment, same matched-null discipline, no permutation. It too is 0/5 and
   5/5 throughout.

   **Once every detector is judged against a matched null, all three agree on
   all seven datasets** — perfect specificity on the four without reticulation,
   detection at every seed on the three with it — except that Φw cannot be
   measured on the gap-blocked data at all. The incumbent is not beaten on power,
   and **this method's case cannot rest on power.** That is a correction to how
   Paper F edition 1 framed the missing comparison: it was framed as the
   validation we lacked, and it turns out to be a three-way tie.

   Two asymmetries survive the correction, one each way. The incumbent's null is
   a **permutation of site order**, so it needs no matched simulation at all;
   both the certificate and the calibrated average need one, and that is a real
   operational cost the incumbent does not pay. Conversely the incumbent returns
   **no measurement** on the gap-blocked dataset, having too few biallelic
   gap-free informative columns once complete deletion is applied, where both of
   the others correctly return "no reticulation". And the permutation null is not
   free of consequence either: it is precisely the null class the RDP5 manual
   reports as blind to misalignment artefacts.

   What this method has that the incumbent does not is not power: it is an
   **exhibitable witness** a kernel recomputes, and the **absence of a
   permutation null**. The second is load-bearing given §11.3 — the RDP5 manual
   states misalignment artefacts are detected with significant p-values in the
   real alignment yet are undetectable in the shuffled alignments, so a
   site-order permutation null, which is precisely Φw's null, is blind by
   construction to the artefact class our ladder measures. Neither claim is
   established by this grid, and both are now stated as the actual case rather
   than power.
3. **Buneman's converse is formalisable, and the first step is now done.** The
   constructive proof factors through the Gromov product rather than attacking
   the tree directly, and that reduction is pure algebra. It is now proved:

   ```
   theorem fourPoint_iff_gromov_strongTriangle {d : α → α → ℚ}
       (hsym : ∀ x y, d x y = d y x) :
       IsFourPoint d ↔ ∀ r, StrongTriangle (gromov d r)
   ```

   with `gromov d r x y = (d x y − d r x − d r y) / 2` and `StrongTriangle f`
   the max inequality `∀ x y z, f x z ≤ max (f x y) (f y z)`. Both directions are
   relabelling plus symmetry — instantiate at `(a,b,c,e) := (x,z,r,y)` one way and
   root at `c` with the order `(a,e,b)` the other. `#print axioms` gives
   `[propext, Classical.choice, Quot.sound]`, no `sorryAx`.

   Two details worth keeping. The condition is stated as `StrongTriangle` and
   **not** as `IsUltrametric`, because the Gromov product takes negative values:
   calling it an ultrametric would import a positivity hypothesis it does not
   satisfy, and the literature's phrase "the Gromov product is an ultrametric"
   means only the max inequality. And the kernel earned its place here — the
   first attempt instantiated the reverse direction at the wrong permutation,
   `linarith` failed on two goals, and the axiom gate reported `sorryAx` rather
   than letting a plausible-looking proof through.

   What remains cited is the *tree construction*: the equidistant
   cherry-picking recursion and the re-rooting that turns it back into a tree
   metric. So the chain is now: certificate ⟹ not four-point ⟹ (proved) some
   root's Gromov product fails the max inequality ⟹ (cited) not tree-realisable.
   Still no formalisation of tree metrics, splits or Buneman's theorem was found
   in Lean, Isabelle or Coq, so finishing it would be first-of-kind.
4. **Our statistic is a maximum; the tests that work are averages — so the
   comparison was run.** A maximum over subsets is an extreme-value statistic
   whose floor grows with taxon count (the quartet floor moves 5–7 → 9–10 from 16
   to 32 taxa), so margins are not comparable across population sizes — the
   incommensurability trap again, now latent in the primary statistic. The
   published treelikeness tests avoid it by averaging. `subset::quartet_distribution`
   now computes the averaged statistics (mean δq, Q-residual, proportion of
   violating quadruples) and `reticulation experiment` reports a seed-paired
   head-to-head against the hard null, D2. **The result did not go the way the
   critique predicted:**

   | statistic | what it is | beats the relaxed-clock null | worst-seed separation |
   | :--- | :--- | :--: | --: |
   | quartet max margin | the certificate's own, exact, with a witness | 5/5 | **16.57×** |
   | mean δq | the averaged statistic the literature favours | 5/5 | 5.07× |
   | proportion violating | the variant the literature reports as failing | 5/5 | 1.22× |

   All three rank both reticulate populations above the relaxed-clock null at
   every seed, so **none** of them inherits the confound that destroyed the
   triple statistic. But at fixed taxon count the maximum separates more than
   three times as widely as the mean, and `proportion_violating` barely separates
   at all — which is consistent with the published finding that that variant
   fails as a test, obtained here independently.

   What this does *not* settle: separation at one taxon count is not the max
   statistic's weakness. Its weakness is incomparability *across* sizes, and this
   grid is all at one size. The honest position is that the maximum is better
   than assumed and still has the structural defect, so both are reported and
   neither is retired. Note also that the averaged statistics are **floating
   point and are not certificates** — only the max-margin row is exact and
   kernel-checkable, because a mean has no witness to exhibit. That is the real
   trade, and it is now priced.
5. **The calibration argument generalises beyond this method.** RDP5's manual
   states that misalignment-induced false positives *"will be detected as
   recombination events with significant p-values"* but *"will likely be
   undetectable in the shuffled alignments"*. Column-shuffling nulls are blind to
   this artefact class by construction, which is an argument about how every
   method in this family should be calibrated, independent of the four-point
   statistic.

Two framings to keep honest. No recombination-detection benchmark varies
alignment quality as an independent variable — so §11.3's measurement is
**confirmatory but novel**, the first direct measurement of a mechanism the field
already asserts, not a surprise. And four-point *satisfaction* does not certify
the absence of reticulation, nor does a violation pin down the tree
(arXiv:1405.2965), which bounds what a certificate means even granting Buneman.

Full report: [`reports/Reticulata next research directions.md`](../../../../reports/Reticulata%20next%20research%20directions.md) in the workspace root, with
seven note files behind it.

## 12. Phase 4 — published breakpoints, and the resolution limit they reveal

§11's real-data result rested on two weak things: a literature *expectation* that
a family recombines, and an alignment from our own aligner, whose error we
therefore had to measure and gate on. Phase 4 replaces both — and finds the
method's resolution limit by doing so.

### 12.1 Ground truth that is named rather than expected

The Pango designation process records SARS-CoV-2 recombinants as `X*` lineages
with their parent lineages and, for many, an explicit numeric breakpoint interval
in reference coordinates, each traceable to a reviewed proposal. Nextstrain
publishes an open, non-GISAID release **already aligned** by Nextclade to
MN908947 — the same coordinate frame the breakpoints are quoted in. So:

| population | recombinant | parents | published breakpoint |
| :--- | :--- | :--- | :--- |
| `xcm_XCM_triad` | XCM | XBB.2.3 × DV.7.1 | between 22017–22031 |
| `xcp_XCP_triad` | XCP | EG.5.1.10 × JD.1.1 | between 14190–14633 |

Sixteen genomes each — 5 + 5 from the parent lineages, 6 from the recombinant —
selected by `scripts/reticulata_sarscov2_breakpoint_populations.py` under filters
fixed in advance (coverage ≥ 0.97, missing data ≤ 1000, first *n* in strain-name
order), with no shortfall in any lineage. The script streams ~2.2 GB and stores
only what it selects.

**The caveat, stated rather than buried.** This removes *our* aligner from the
trust surface and replaces it with Nextclade's, which is reference-anchored too —
the same construction, and therefore the same class of error §11.3 measures.
Provenance improves: the alignment is produced independently, by a standard tool,
and is citable. Measurement coverage does not: we have not measured Nextclade's
error, so the admissibility gate calibrated for our aligner does not transfer and
is not applied. The gap load is reported instead (0.003 on both).

### 12.2 A prediction about *which* sequences, not about a margin

If the certificate detects the recombination the designation process named, its
dominant pairing should separate the two **parent** lineages with a recombinant
on the other side — the planted-recombinant shape. That constrains which four
sequences the method picks, which a margin cannot, and it was registered before
the test ran.

| population | prediction | certificates separating the parents | margin | matched floor |
| :--- | :--: | :--: | --: | :--- |
| XCP | **held** | 6 of 10 | 2/28920 | **none** (degenerate) |
| XCM | **failed** | 0 of 10 | 4/28511 | 3/29531 |

On XCP the method picked the two true parents, EG.5.1.10 and JD.1.1, as the far
pair in six of its ten strongest certificates. On XCM it did not: every one of
its ten certificates paired a parent with a recombinant on each side.

One of two is a split result, and it is reported as one.

### 12.3 The finding that matters is a resolution limit

Both margins are at the edge of what the statistic can resolve, and the artifacts
now say so themselves:

| population | mean p-distance | implied differing sites genome-wide | margin, in sites |
| :--- | --: | --: | --: |
| XCP | 0.00057 | ~17 | **2** |
| XCM | 0.00127 | ~38 | **4** |

A margin is a count of differing sites over a quadruple's common columns. When a
population differs by a few tens of sites across 29 903 columns, a margin of two
or four sites cannot carry a verdict however the floor comparison comes out. And
XCP's floor comparison is itself degenerate: the matched null produced **no
violating quadruple in any replicate**, so `exceeds_floor` reduces to "the
observed margin is positive" and carries no calibration. Both the degeneracy flag
and the resolution warning are now computed in the output rather than left to
prose, the same discipline §12.4 applies to the topology experiment.

**So the honest conclusion is about applicability, and it cuts against the
direction we pursued:** SARS-CoV-2 recombinants between recently-diverged
lineages sit *below* this method's resolution. There is a real tension here worth
naming, because it will recur — **the datasets with the best ground truth are the
least divergent**, since recombinants get named between lineages that are close
enough to be confidently assigned. Phase 3's populations, at mean p-distance
0.0088 (Hepatitis B) and 0.11 (Coxsackievirus A6), were far better suited to the
statistic and far worse supplied with ground truth.

The one encouraging signal, and it is labelled as suggestive rather than
evidence: on XCP the *pairing* was right in 6 of 10 certificates at a margin of
two sites. The method pointed at the correct two parents in a regime where its
magnitude means nothing.

### 12.4 Persistent homology: the experiment is inconclusive, and the artifact says so

`reticulation topology` tests something the literature review found missing — no
peer-reviewed independent calibration of persistent-homology specificity for this
application appears to exist, and the one unpublished evaluation reports the
cycle rate tracking point-mutation rate rather than recombination rate. Our grid
has the ground truth to settle it: the saturated dataset differs from the
baseline only in having ten times the substitution rate.

The experiment **cannot answer the question as set up**, and the output reports
`"this_experiment_is_interpretable": false`. The matched null produced no cycle
in **495 of 495 replicates** on five of the seven datasets, so the p-value
collapses to two values — the floor 0.010 whenever the observed value is
non-zero, and 1.000 otherwise. The bootstrap calibrates nothing; the test reduces
to "did a cycle appear at all". Consequently the 2/5 counts on the three
reticulate datasets *and* on the saturated one are that binary, and must not be
read as power or specificity. The concern is **not** reproduced here; neither is
it refuted.

Diagnosed cause: at 16 taxa with a Rips threshold at the median off-diagonal
distance, the complex is too small for H₁ to populate. The fix is more taxa, so
that cycles occur in the null often enough to calibrate against, and/or a larger
threshold — which reorders the plan in
[`RETICULATA_NUMERICS_TOPOLOGY_ML.md`](./RETICULATA_NUMERICS_TOPOLOGY_ML.md):
the bit-parallel counting of its §2.1 has to come *before* the topology test, not
after it.

One observation appeared to survive the degeneracy: when cycles did appear, the
reciprocal-mosaic dataset's total persistence was an order of magnitude larger
(0.0107, 0.0095) than the single-recombinant or saturated ones (0.0003–0.0013),
suggesting the *magnitude* might separate where the binary verdict could not.
**§13.3 refutes that at sixty-four taxa**, where the mosaic persistence sits
within a factor of two of both no-reticulation datasets. It is recorded here as
a hypothesis that was tested and failed rather than deleted.

### 12.5 The four-way comparison, corrected

§11.8's table reported the certificate against the strict-clock dataset's floor
and the competitors against matched nulls. That was the incommensurability this
document has recorded three times, committed a fourth time against ourselves. The
corrected table, with every detector judged against a null matched to each
dataset's own width, divergence and gap load, is in §11.8 — and with that
correction the certificate is specific on all four no-reticulation datasets and
detects all three reticulate ones at every seed, as do Φw and the calibrated
averaged statistic. Φw alone cannot be measured on gap-blocked data.

Raw output: [`data/reticulata_phase4/`](../data/reticulata_phase4/) and
[`data/reticulata_phase2/topology.json`](../data/reticulata_phase2/topology.json).

## 13. The max statistic's defect, measured across taxon counts

§11.8's head-to-head concluded that the certificate's real weakness is not power
but **incomparability across taxon counts**, and it argued that from two data
points: the quartet floor moving from 5–7 at sixteen taxa to 9–10 at thirty-two.
Arguing from two points is how this project has twice talked itself into an
incommensurable comparison, so the defect is now measured across a range.

### 13.1 Bit-parallel counting, because the measurement was unaffordable

Sixty-four taxa is 635 376 quadruples, and the byte path visits 4000 columns per
quadruple — 2.5e9 column visits per population. `bitcount::BitAlignment` packs
each sequence into three bitvectors over columns (the base's two bits under
`A=00, C=01, G=10, T=11`, plus a mask of definite bases), so a quadruple's seven
counts come from one pass over `L/64` words: intersect four masks for the common
columns, then `popcount` the differing bits for each of the six pairs inside it.

It is a **second** implementation, not a replacement. The byte path of
`Alignment` remains the specification;
`bit_and_byte_counting_agree_on_random_alignments_with_gaps` asserts identical
counts on random alignments containing gaps, ambiguity codes and mixed case at
widths 1, 7, 63, 64, 65, 130 and 257 — the partial-final-word cases a packing bug
would hide in — and a second test does the same on gap-blocked simulated data.
An optimisation that disagrees with the obvious implementation is a silent wrong
answer, and these counts are what the whole method rests on.

Measured, not asserted: **21.6× faster** at thirty-two taxa (0.24 s against
5.14 s over 35 960 quadruples), with the two implementations returning the same
exact margin. The `scale` command re-checks that agreement on every run and exits
non-zero if it ever fails.

### 13.2 The floor grows, the signal does not

Three seeds per rung, 4000 nt, exhaustive over every quadruple. Margins are exact
rationals with denominator 4000; the table gives numerators.

| taxa | quadruples | violating fraction, clonal | strict-clock floor | reciprocal-mosaic margin | worst-seed separation |
| --: | --: | --: | :--- | :--- | --: |
| 16 | 1 820 | 0.656 | 7, 5, 6 | 136, 112, 124 | **19.43×** |
| 32 | 35 960 | 0.755 | 10, 9, 10 | 148, 154, 108 | **10.80×** |
| 64 | 635 376 | 0.821 | 17, 16, 14 | 156, 148, 166 | **9.18×** |

The floor roughly **doubles** from sixteen to sixty-four taxa (mean 6 → 15.7)
while the mosaic margin stays essentially flat (124 → 157). So the separation
**halves**, 19.4× → 9.2×, over a fourfold increase in population size. That is
the defect, measured: more subsets, more chances at an extreme pair of rates, a
higher null level, and a signal that does not grow to match.

Two readings, and both matter:

- **Detection at a fixed taxon count is not lost.** Nine-fold separation at
  sixty-four taxa is still wide. The growth shifts both numbers rather than
  drowning the signal.
- **A threshold calibrated at one size does not transfer to another.** "Twice the
  floor" means something different at sixteen taxa than at sixty-four, which is
  exactly why §11.6's real-data verdicts use a floor matched to each population
  and why the raw-margin column of §11.8 was the wrong comparison. The defect is
  real, bounded, and already handled by the matched-null discipline — but it must
  not be forgotten when a margin is quoted.

A third observation is worth recording because it looks alarming and is not. The
fraction of quadruples that violate the four-point condition on a **clonal,
strict-clock** population rises from 0.656 to 0.821 as taxa grow: at sixty-four
taxa, 82% of all quadruples violate. Finite sequences are never exactly
additive, so this is expected — and it is precisely why the *margin* rather than
the *count* carries the inference, and why `proportion_violating` was the variant
the literature reported as failing as a test and that §11.8 measured separating
by only 1.22×.

Raw output: [`data/reticulata_phase2/scale.json`](../data/reticulata_phase2/scale.json).

### 13.3 Topology at sixty-four taxa: where the calibration becomes valid, the statistic fires on the wrong things

§12.4's topology experiment was uninterpretable because its matched null produced
no cycle in 495 of 495 replicates. The diagnosed fix was more taxa, which §13.1
made affordable. Re-run at sixty-four taxa, three seeds, 49 bootstrap replicates
(p-value floor 0.02):

| dataset | reticulation | mean H₁ total | mean cycles | null replicates with no cycle | significant |
| :--- | :--: | --: | --: | :--- | :--: |
| D1 clonal, strict clock | no | 0.00000 | 0.0 | 147/147 | 0/3 |
| D2 clonal, relaxed clock | no | 0.00000 | 0.0 | 147/147 | 0/3 |
| D3 one recombinant | **yes** | 0.00000 | **0.0** | 147/147 | **0/3** |
| D4 relaxed clock + one recombinant | yes | 0.00008 | 0.3 | 147/147 | 1/3 |
| D5 reciprocal mosaics | yes | 0.01475 | 1.7 | 147/147 | 3/3 |
| D6 gap blocks | **no** | 0.01072 | 13.3 | **66/147** | **3/3** |
| D7 saturated | **no** | 0.00867 | 7.0 | **138/147** | **3/3** |

The output still reports `"this_experiment_is_interpretable": false`, because
five of the seven nulls remain degenerate. But **two nulls are now
non-degenerate, and they are the two datasets containing no reticulation at
all** — and the statistic fires on both at 3/3:

- **D7, saturated.** No reticulation; differs from D1 only in having ten times
  the substitution rate. This is the concern from the one unpublished evaluation
  we found, now observed on data with a known answer in a cell where the
  calibration is valid: the statistic responds to **mutation rate**.
- **D6, gap blocks.** No reticulation; differs from D1 only in missing data. The
  statistic responds to **gaps**, and 13.3 cycles per population is by far the
  largest count in the table.

Against that, **D3 — a genuine single recombinant — produces no cycles at all**,
and is missed at 0/3. And the magnitude hypothesis of §12.4 fails: the mosaic
persistence (0.0147) is within a factor of 1.4 to 1.7 of the two
no-reticulation datasets (0.0107, 0.0087), so magnitude does not separate them
either.

**Scope of this claim, stated carefully.** This is a result about the statistic
*as implemented and calibrated here* — Vietoris–Rips H₁ on a p-distance matrix
with the threshold at the median off-diagonal distance, against a matched
relaxed-clock null — at three seeds. It is not a refutation of the published
method, which used different data, thresholds and emphasis, and we have not
reproduced that work's own analysis. What it does establish is that this
construction, on this grid, cannot be relied on to distinguish reticulation from
divergence or from missing data.

### 13.4 Why the gap result is structural, not just empirical

D6 deserves more than a row. The audit's finding 6 was that pairwise deletion
gives different pairs different denominators, and that this alone manufactures
violations on a genealogy that is exactly a clock tree; §10.3's answer was
**complete deletion within the subset**, which makes every distance in a
certificate share one denominator.

**The Rips construction cannot do that.** Persistent homology needs a full
distance matrix over the population, and a matrix has one distance per pair —
there is no "within this subset" to restrict to, because the complex is built
from all pairs at once. So the topological route necessarily uses pairwise
deletion and necessarily inherits the confound the four-point certificate was
redesigned to avoid. D6's 13.3 cycles per population on gap-blocked clonal data
is that confound showing up as topology.

This is an argument from construction rather than from these three seeds, and it
is the strongest thing in this section: on gappy data, a statistic defined over
a single global distance matrix has no analogue of complete deletion available
to it.

Raw output: [`data/reticulata_phase2/topology_n64.json`](../data/reticulata_phase2/topology_n64.json).

## 14. The filtration funnel, and closing the topology line

§13.3 left one instruction: change the filtration or stop. Five hypotheses were
written down for *why* Vietoris–Rips H₁ fires on gap-blocked and saturated
populations while missing a single recombinant, three were selected on
**cost-to-decide** — each cheap and each decisive about a different cause — and
the criterion was fixed before the run: a variant succeeds if the
reciprocal-mosaic population is significant at every seed **and** both
no-reticulation populations at none.

Sixty-four taxa, three seeds, 49 bootstrap replicates, matched null per cell.

| variant | hypothesis | mosaics | single recombinant | gap blocks | saturated | criterion |
| :--- | :--- | :--: | :--: | :--: | :--: | :--: |
| baseline, p-distance total | — | 3/3 | **0/3** | 2/3 | 3/3 | no |
| **F1** Jukes–Cantor total | saturation: p-distance concavity | 3/3 | **0/3** | 2/3 | 3/3 | no |
| **F2** global complete deletion | denominators: finding 6 as topology | 3/3 | **0/3** | **0/3** | 3/3 | no |
| **F3** longest bar, not total | many weak cycles conflated with few strong | 3/3 | **0/3** | **0/3** | 3/3 | no |

**No variant meets the criterion.** But three things are now settled.

**F1 is refuted.** Jukes–Cantor correction leaves the saturated false positive
at 3/3 and *raises* its statistic (0.00867 → 0.01421). The sensitivity is not
p-distance concavity.

**F2 is confirmed, and confirms the structural argument of §13.4.** Global
complete deletion takes the gap-blocked population from 13.3 cycles and 2/3
significant to **0.0 cycles and 0/3**. The gap false positive *was* the audit's
finding 6 reappearing as topology, exactly as the construction argument
predicted. This is the round's most satisfying result: the same confound, in a
second guise, removed by the same idea. F3 removes it too, by a different route
— distinguishing many weak cycles from few strong ones.

**The miss is not a filtration artefact.** A single recombinant produces
**zero cycles in all four variants**. That is not something a filtration choice
fixes, and it is what closes the line, because the certificate does not have the
problem. Measured on the same populations at the same taxon counts:

| taxa | clonal floor | single-recombinant margin | above floor | mosaic margin |
| --: | :--- | :--- | :--: | :--- |
| 16 | 7, 5, 6 | 60, 52, 53 | **3/3** | 136, 112, 124 |
| 32 | 10, 9, 10 | 60, 67, 98 | **3/3** | 148, 154, 108 |
| 64 | 17, 16, 14 | 89, 82, 83 | **3/3** | 156, 148, 166 |

At sixty-four taxa the four-point certificate detects the single recombinant at
every seed with a roughly fivefold margin over the floor, on the population
where the topological statistic finds no cycle at all.

### 14.1 Decision: the topology line is closed

Not because topology is the wrong idea. The warrant — *a tree metric admits no
unfilled cycle* — is the same fact the certificate rests on, and it is the reason
any of this works. What is closed is the **global-distance-matrix
implementation** of it, for two reasons that the funnel established rather than
assumed:

1. It **cannot use per-subset complete deletion.** A complex is built from all
   pairs at once, so there is no "within this subset" to restrict to; F2 shows
   that the available substitute, one global denominator, does remove the gap
   artefact but at the cost of every column that any taxon is missing.
2. It **misses the signal the certificate catches.** Zero cycles on a single
   recombinant in four variants, against 3/3 detection by the certificate at
   matched taxon count.

So the programme's topological content is the certificate itself: the same
theorem, evaluated combinatorially per quadruple instead of through a global
complex. Evaluating it per subset is precisely what makes complete deletion
within the subset possible, which is what makes the gap confound avoidable rather
than merely reducible.

The two expensive hypotheses held back from the funnel are closed with it. The
residual filtration (F4) addresses a baseline problem that is not the failure
here. The quartet-margin complex (F5) would use per-subset complete deletion and
so would not inherit finding 6 — but it adds information only about how
violations *overlap*, and the certificate already captures the quadruple-level
signal that the global construction was missing. F5 is therefore demoted from
"the justified pivot" to "optional and speculative", and §3.3 of the directions
document is marked accordingly.

One honest limit on all of this: three seeds, one generative model, 64 taxa, and
the saturated cell's null is degenerate in every variant, so that column carries
no calibration and its 3/3 is the binary "a cycle appeared". The two conclusions
that do not depend on that cell are F1's refutation and the single-recombinant
miss, and those are the two the decision rests on.

Raw output: [`data/reticulata_phase2/filtration.json`](../data/reticulata_phase2/filtration.json)
and [`scale.json`](../data/reticulata_phase2/scale.json).

## 15. The geometry the certificate was already doing

§14 closed the Vietoris–Rips route and concluded that this programme's
topological content is the certificate itself. That is a slogan unless it is made
precise. It is precise, it is a theorem — and the first version of this section
got the key definition wrong, so the withdrawal comes first.

### 15.1 A correction: the isolation index takes the maximum, not the minimum

The first version of this section, published as **Paper F edition 5
(10.5281/zenodo.23589159)** and **textbook edition 9
(10.5281/zenodo.23589161)**, defined

```
WRONG:  α(S) = ½ ( min of the other two pairing sums − S's own sum )
```

and called it the Bandelt–Dress isolation index. It is not. The definition takes
the **maximum** of the two crossing sums, clamped below at zero:

```
α(A|B) = ½ · min over a,a' ∈ A, b,b' ∈ B of
             max{ 0, cross₁ − own, cross₂ − own }
```

which on a single quadruple is `½ max(0, max(cross₁, cross₂) − own)`. The
degenerate case confirms it: for `A = {x}`, `B = {y,z}` the formula gives
`½(d(x,y) + d(x,z) − d(y,z))`, the Gromov product — consistent with
§11.8's `fourPoint_iff_gromov_strongTriangle`.

**Why it slipped through.** The two forms **agree exactly on a tree metric**,
where the crossing sums are equal, and differ only where the four-point condition
fails. So the mistake was invisible on precisely the data where the condition
holds. Worse, the check used to pick between them — "the index should equal the
internal edge length on a tree" — is satisfied by *both* forms, so it was not a
discriminating test and was treated as one. This project's rule is to verify every
non-classic citation against a primary record before deposit; the rule was applied
to the DOIs in the reference list and not to a definition assumed known.

The theorems proved about the wrong function were true *about that function*; only
the name was wrong. And the corrected statement is **stronger**, which is luck
rather than vindication.

### 15.2 The margin is twice an isolation index

Write a quadruple's three pairing sums sorted as `s₁ ≤ s₂ ≤ s₃`. Under the
correct definition the three splits' indices are

| split | index |
| :--- | :--- |
| smallest sum (the supported one) | `½(s₃ − s₁)` |
| middle sum | `½(s₃ − s₂)` |
| largest sum | `0` |

so the **middle** split's index is exactly half the certificate's margin
`s₃ − s₂` — one index, not a gap between two. Five theorems in
`GeneticaReticulataTree.lean`, 0 `sorry`, no axioms beyond the whitelist:

| theorem | content |
| :--- | :--- |
| `buneman_middle_eq_half_margin` | the middle split's index is `(s₃ − s₂)/2` |
| `buneman_largest_eq_zero` | the dominant split's index vanishes |
| `buneman_supported_eq_half_gap` | the smallest split's index is `(s₃ − s₁)/2` |
| `fourPoint_iff_middle_buneman_zero` | **the four-point condition is the vanishing of the middle index** |
| `buneman_supported_is_edge_when_fourPoint` | where four-point holds, the supported index is `(s₂ − s₁)/2`, the internal edge's length |

### 15.3 The reading, which is why the correction was worth making

On a tree metric exactly **one** of the three splits carries a positive index —
the one the tree supports. A certificate says the middle index is positive too,
and **any two distinct splits of a four-set are incompatible**. So a quartet
certificate is the statement that *two incompatible splits both carry positive
isolation index*, which is the textbook condition for a split system to fail to
be a tree.

The certificate is therefore not *analogous* to the split-decomposition notion of
non-tree-likeness. At the quadruple level it **is** that notion. The code checks
this as an identity rather than asserting it: across all seven datasets at five
seeds, quadruples with two positive indices number **47 802**, and four-point
violations number **47 802** — equal, as they must be.

Note what the proofs do not use: no tree, no complex, no metric space, only three
rationals. That is exactly why the quantity survives being computed one subset at
a time, where a complex over a global distance matrix cannot, and so why §14
closed an implementation rather than an idea.

### 15.4 Checked against the data, not only proved

`subset::buneman_indices` computes each split's index from the per-split
definition — the larger of the *other two* sums minus its own, clamped at zero —
rather than from the sorted sums, so a mistake surfaces as a mismatch instead of
cancelling out. `quartet_distribution` then checks on every quadruple that the
middle index equals the margin numerator and the dominant index vanishes.

Across all seven datasets at five seeds: **the identities hold on 75 500
quadruples and fail on 0.**

### 15.5 Reticulation lives in the index that should be zero

Largest supported-split index per population, first three seeds:

| dataset | reticulation | max supported isolation |
| :--- | :--: | :--- |
| D1 clonal, strict clock | no | 0.02863, 0.02637, 0.02750 |
| D3 one recombinant | **yes** | 0.02863, 0.02637, 0.02750 |
| D5 reciprocal mosaics | **yes** | 0.02863, 0.02637, 0.02750 |
| D2 clonal, relaxed clock | no | 0.03587, 0.03963, 0.02312 |
| D6 gap blocks | no | 0.03388, 0.02981, 0.02980 |
| D7 saturated | no | 0.19350, 0.19225, 0.19850 |

**The supported-split index is identical, to five decimal places, on the clonal
tree and on both reticulate populations built from it.** It is blind to the
reticulation — as it should be, since the strongest supported split is a deep
split of the underlying tree and adding recombinants does not disturb it. What it
tracks is divergence: the saturated population's index is seven times the others'.
This observation survives the correction of §15.1 unchanged.

So the statement the pivot arrives at:

> **The supported split's index carries the tree. The index that a tree forces to
> zero — the middle one — carries the reticulation. The certificate measures that
> index.**

That also explains §13.3's failure after the fact. A Rips filtration is driven by
distance *magnitudes*, which is the quantity the supported index tracks and hence
the quantity that follows divergence; the signal lives in an index that a tree
sets to zero, and magnitudes do not expose it. A statistic that fires on saturated
data and misses a single recombinant is what the geometry predicts.

Raw output: the `the_margin_is_twice_an_isolation_index` block and the per-seed
`max_supported_isolation` and `two_positive_indices` in
[`data/reticulata_phase2/experiment.json`](../data/reticulata_phase2/experiment.json).

## 16. The null's floor, predicted instead of simulated

Phase 3 of §5's plan — the margin's null distribution — has been the oldest open
item in this programme, and the literature review found no analytic sampling
distribution for it either. Every calibrated answer here costs a matched
simulation: `real-data` spends a bisection plus five tree simulations per
population, each calibrated statistic spends 99 replicates per dataset-seed, and
§12.4 and §12.3 both reported counts against nulls that could not vary because
there was no cheap way to know in advance where the null is degenerate.

`reticulation null-sweep` and `scripts/reticulata_null_surrogate.py` close that,
with a restriction.

### 16.1 Why a fitted model is admissible here at all

**Only the floor is predicted.** The certificate's margin remains an exact
integer comparison the Lean kernel can check, so a wrong surrogate costs power or
wastes a query — it cannot make a certificate wrong. That is the admissibility
condition of
[`RETICULATA_NUMERICS_TOPOLOGY_ML.md`](./RETICULATA_NUMERICS_TOPOLOGY_ML.md) §1,
and it is the only reason this sits inside a method whose whole claim is
checkability.

The sweep records, for each of **288 parameter cells** at five seeds, the maximum
quartet margin on a clonal population with no reticulation — the same floor
`real-data` computes. No bisection: `mu` is swept directly and the *achieved*
divergence recorded, which is the covariate a surrogate wants and removes the
40-step tuning from every cell. Bit-parallel counting (§13.1) is what makes 32
taxa across 288 cells affordable.

### 16.2 The falsifiable check, fixed before the fit

The model is a scaling law, not a black box:

```
log(floor) ~ intercept + b·log(L) + c·log(mean p-distance)
             + d·log(quadruples) + e·gap_load + f·clock_spread
```

and `b` has an **independent prediction**. The margin is a difference of counts
over a denominator of order `L`; the numerator carries multinomial noise scaling
as `√L`, so the floor should scale as `√L / L = L^(-1/2)`. §10.6's length-scaling
diagnostic measured exactly that signature independently (0.42–0.53 across two
divergences). So a tolerance of ±0.15 around **−0.5** was fixed before fitting.

**Fitted: `b = −0.556`.** The scaling law is confirmed rather than merely
interpolated — the floor is noise-dominated, and the surrogate rests on a
mechanism instead of on a curve through a table.

### 16.3 Held-out regions, and the honest verdict

Tested by **region**, not by random split, because extrapolation is the failure
mode that matters:

| held out | train rows | test rows | median relative error | p90 |
| :--- | --: | --: | --: | --: |
| longest alignments (`L = 8000`) | 903 | 343 | **20.5%** | 55.5% |
| largest taxon count (32) | 781 | 465 | **28.7%** | 66.5% |

In-sample median error is 23.3%, so the model is not overfitting — it is simply
this accurate.

§4.1 of the directions document fixed the abandonment criterion in advance:
*abandon if held-out error at the quantiles that matter exceeds the margin
differences the method is trying to resolve.* Applying it honestly:

- **Wide-margin verdicts survive.** The Hepatitis B population's margin is 4.5×
  its floor; a 30% error in the floor does not move that verdict.
- **Borderline verdicts do not.** The Coxsackievirus population clears its floor
  by 1.23×; a 25% floor error flips it.

So the criterion is **partially met, and the outcome is a restricted deployment
rather than adoption**: the surrogate is adequate for screening parameter regions
and for predicting degeneracy, and it must not decide a verdict near the floor.
Those still get a simulation. Stating that as a restriction rather than reporting
the 20% median as a success is the whole point of having fixed the criterion
first.

### 16.4 What it buys immediately, demonstrated rather than claimed

**Degeneracy is now predictable.** Ten of the 288 cells are degenerate — no
replicate produced a violating quadruple — and they cluster exactly where
expected: the smallest taxon count (8) at the lowest substitution rate, i.e. few
subsets and little divergence. These are reported, never imputed: a floor of zero
has no logarithm and is not a floor any model should invent.

And the surrogate reproduces §12.3's hard-won conclusion in advance. Its
calibrated divergence range is **0.00239 to 0.51429**. The XCP
published-breakpoint population sits at mean p-distance **0.00057** — **4.2×
below the lower bound**, so the surrogate refuses to predict for it and flags it
as outside the calibrated regime. That is precisely the verdict §12.3 reached
only after a 2.2 GB download, an alignment, a test, and a degenerate floor.
Having the table first would have saved all of it.

That is the honest case for this piece of machinery: not that it predicts floors
well, because at 20–29% it predicts them only adequately, but that **it knows
where it cannot be used, and that boundary is itself the finding** the programme
kept paying to rediscover.

Raw output: [`data/reticulata_phase3/null_sweep.json`](../data/reticulata_phase3/null_sweep.json)
and [`null_surrogate.json`](../data/reticulata_phase3/null_surrogate.json).

## 17. The certificate promoted to a split, and three controls that failed first

§15 established that a quartet certificate *is* the split-decomposition
condition for four taxa. A quartet is not what a biologist reads, though; a
**split** is — this group of sequences against that group. This section promotes
the certificate to that level. The formal half worked first time. The empirical
half failed its negative control three times, and each failure found something
that mattered more than the original question.

### 17.1 A conjunction, not a minimum

Bandelt and Dress define a split's isolation index as a **minimum over the
quadruples straddling it**. Taken literally that is unusable here: under
complete deletion within the subset every quadruple has its own denominator, and
a minimum over quantities with different denominators is an incommensurable
comparison — the defect this project has caught in its own apparatus more than
once. So the split certificate is a **logical conjunction** of per-quadruple
integer comparisons:

> `A|B` is certified iff for **every** quadruple `{a,a',b,b'}` with `a,a' ∈ A`
> and `b,b' ∈ B`, the split's own pairing sum is strictly below the larger of
> the two crossing sums — each comparison exact, over that quadruple's own
> denominator.

Nothing is averaged, nothing needs a common metric, and every conjunct is
kernel-checkable. The split-level statement inherits exactness from the
quadruple level rather than paying for it.

### 17.2 The lift, proved

Two splits are **incompatible** when all four of `A∩C, A∩D, B∩C, B∩D` are
non-empty. Take one taxon from each: on that quadruple the two splits induce two
*different* pairings. If both splits are certified, both of those pairings are
non-maximal there, so the third is strictly the largest — exactly a four-point
violation.

| theorem | content |
| :--- | :--- |
| `third_pairing_strictly_largest_of_two_certified` | two certified pairings force the third strictly largest |
| `fourPoint_violated_of_two_certified` | the same, stated as a four-point violation |
| `fourSet_splits_pairwise_incompatible` | no two of a four-set's three splits are compatible — **the claim editions 5/9 asserted without proof** |
| `splitsCompatible_self` | a split never crosses itself, so one certified split yields no spurious witness |

`GeneticaReticulataTree.lean` is **27/27 declarations, 0 `sorry`, no axioms
beyond the whitelist**. So a population-level claim *forces* a quadruple-level
certificate, with the witness **constructed** from the four intersections rather
than searched for. The code builds that witness and checks it really violates
four-point: true on **350 of 350** dataset-seed-floor combinations.

### 17.3 Failure 1 — the criterion had no noise floor

With the definition as stated, the negative control failed **0 of 5**: a clonal
strict-clock population with no reticulation certified 8 splits where its tree
has 5, and 6 of the resulting pairs were incompatible.

The exhibit showed why immediately. Pairing sums `[66, 110, 112]` over a
denominator of 2000 — and on an *exact* tree metric the two crossing sums are
**equal**. Their observed difference is sampling noise of order `√L`, and
`own < max(cross)` reads the sign of that noise. The split with 110 was
certified on a 2-count margin in 2000.

**And one of the predictions was rigged.** P4 asked whether a planted
recombinant was implicated, accepting it either in a witness quadruple or
anywhere in the incompatible pair's splits. A split's two sides *cover every
taxon*, so the second branch is true whenever any incompatible pair exists. It
was a forced test and it recorded a forced pass; only the witness quadruple,
four taxa of nine, carries information.

### 17.4 Failure 2 — no constant floor can work

Requiring each conjunct to clear a rate `θ` of its denominator, swept over ten
values, **no separating floor exists**. The per-dataset breakdown says why:

| floor rate | D1 | D2 | D6 | D7 | D3* | D4* | D5* |
| ---: | :-: | :-: | :-: | :-: | :-: | :-: | :-: |
| 0.00100 | 5/5 | 5/5 | 1/5 | 1/5 | 0/5 | 0/5 | 0/5 |
| 0.00400 | 5/5 | 5/5 | 5/5 | 3/5 | 2/5 | 2/5 | 0/5 |
| 0.01600 | 5/5 | 5/5 | 5/5 | 5/5 | 5/5 | 5/5 | 5/5 |

(tree-readings per 5 seeds; `*` = reticulate, where a tree-reading is a miss)

The saturated population **keeps exactly 5 certified splits at a floor 64×
higher than the one that erases every other dataset**, because saturation
inflates margins *in rate terms*. A single global rate is the wrong functional
form — which is what §16's surrogate already implied when its fit put a large
coefficient on divergence.

### 17.5 Failure 3 — the matched null was never matched

Calibrating the floor per population against a clonal null matched on divergence
fixed specificity but reported a null whose divergence missed the population's
by **77%**, and only on the two relaxed-clock datasets. That exposed a bug in
machinery shared with the published results:

> `clonal_tree_relaxed` draws its per-branch rate multipliers from the **same RNG
> stream** as `mutate`, and `mutate` consumes a number of draws that depends on
> `mu`. Changing `mu` shifts the stream and changes every downstream rate
> multiplier, so a single realisation's divergence is a **chaotic, non-monotone**
> function of `mu`. `tune_mu_for_divergence` bisected exactly that, with
> `spread = 1.0` in every caller. The strict-clock simulator has nothing
> downstream of its mutations, which is why only relaxed-clock nulls were hit.

Averaging the objective barely helped (worst error 147.5% → 143.7%), so the
assumption underneath was measured instead. At **fixed** `mu`, over 200
realisations:

| `mu` | mean divergence | CV | min | max | max/min |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.001 | 0.00658 | 57.6% | 0.00088 | 0.02495 | 28.4× |
| 0.004 | 0.02677 | 49.1% | 0.00521 | 0.07064 | 13.6× |
| 0.016 | 0.09299 | 48.9% | 0.02272 | 0.22849 | 10.1× |

Each lineage's rate is its parent's times `exp(spread·u)`, so rate multipliers
**compound** across generations and their log-variance grows with depth. The
distributions at `mu = 0.001` and `0.004` overlap heavily. So matching a
relaxed-clock null on divergence **in expectation is ill-posed**: the floor is a
max over replicates whose divergences differ by up to 28×, dominated by
whichever replicate drew a fast-rate tree. The null was inflated — conservative,
but not matched.

### 17.6 What works: rejection sampling on the realised nuisance

A replicate counts only if **its own realised divergence** is within 10% of the
population's. With that, all four pre-registered predictions hold at full count:

| prediction | result |
| :--- | :--- |
| **M1** every non-reticulate population reads as a tree | **20/20** |
| **M2** every reticulate population reads as a non-tree | **15/15** |
| **M3** reticulate populations still certify something | **15/15** |
| **M4** a planted recombinant appears in a witness quadruple | **15/15** |

| dataset | reticulate | floor | certified | `n-3` | trees | incompatible pairs | null acceptance |
| :--- | :--: | ---: | ---: | ---: | :--: | ---: | ---: |
| D1 clonal, strict clock | no | 0.00170 | **5.0** | 5 | 5/5 | 0.0 | 72.5% |
| D2 clonal, relaxed clock | no | 0.00120 | 4.4 | 5 | 5/5 | 0.0 | 17.1% |
| D6 gap blocks | no | 0.00339 | 4.0 | 5 | 5/5 | 0.0 | 52.5% |
| D7 saturated | no | **0.04000** | **5.0** | 5 | 5/5 | 0.0 | 100.0% |
| D3 one recombinant | **yes** | 0.00160 | 9.0 | 6 | 0/5 | 7.0 | 54.7% |
| D4 relaxed + recombinant | **yes** | 0.00190 | 7.4 | 6 | 0/5 | 4.6 | 15.5% |
| D5 reciprocal mosaics | **yes** | 0.00170 | 10.6 | 7 | 0/5 | 18.8 | 51.5% |

D1 and D7 recover **exactly `n-3` = 5** certified splits, their true internal
split count. D7's floor is **23× D1's**, which is the divergence dependence no
constant rate could express. Acceptance rates are reported because they matter:
at 15–17% for the relaxed-clock populations the sampler needs ~6 draws per
accepted replicate, and a population that could not be matched at all would get
**no floor and no verdict** rather than a floor from an unmatched null.

### 17.7 What this does and does not establish

It establishes that the split-level certificate separates reticulation from four
distinct tree-like confounders — relaxed clocks, gap blocks, and saturation —
**perfectly on this grid**, 35 of 35 dataset-seeds, with a biologist-readable
exhibit and a kernel-checkable witness beneath it.

It does **not** establish that on real data. This is simulated data with planted
truth at one parameter setting (2000 columns, 8–10 taxa). And the result arrived
only after three negative-control failures, which is the honest summary of its
maturity: the formal content is proved, the empirical content is one grid.

**The methodological finding is the transferable part.** Matching a null in
expectation is not matching it, whenever the nuisance has a heavy-tailed
realised distribution; and a quantity whose across-realisation CV is ~50%
cannot calibrate anything through its mean.

Raw output: [`split_system.json`](../data/reticulata_phase5/split_system.json),
[`tuning_bias.json`](../data/reticulata_phase5/tuning_bias.json).

## 18. Re-deciding the real populations against a null that is actually matched

§17.5 found that `tune_mu_for_divergence` bisected a single realisation of a
chaotic objective, and that matching a relaxed-clock null on *expected*
divergence is ill-posed when the realised divergence has CV ≈ 50%. That function
set the matched floors behind **the published real-data verdicts**, so those
verdicts had to be re-decided rather than reasoned about.

### 18.1 How badly the published nulls missed

The published artifact recorded each floor's null divergence, so the miss is
recoverable from it directly:

| population | observed divergence | published null | error |
| :--- | ---: | ---: | ---: |
| coxsackievirus A6 | 0.10957 | 0.05588 | **−49.0%** |
| hepatitis B | 0.00882 | 0.00977 | +10.8% |
| JC polyomavirus | 0.00961 | 0.00842 | −12.3% |

The coxsackievirus null sat at **half** its population's divergence. The seven
refused populations are unaffected: they fail the aligner *admissibility* gate,
which never consults a floor.

### 18.2 Two defects pulling opposite ways, so the net had to be measured

Before re-running I claimed the bug's direction was "conservative", and then
predicted the coxsackievirus floor would **rise** once its null was matched at
the correct, higher divergence. **Both claims were wrong as stated**, and the
second was wrong on the record because it ignored a mechanism §17.5 had already
established. Two defects act in opposite directions:

* A **non-converging bisection** misses its divergence target in *either*
  direction, so the null's central divergence can be too low or too high.
* The floor is a **max over replicates**, and over-dispersion always *inflates*
  a max — the extreme draw sets it, regardless of the nominal centre.

So the net is not deducible from the divergence error alone. Measured, the
inflation term dominated in all three populations: **every corrected floor is
lower or equal**, even where the correctly-matched null sits at twice the
divergence.

### 18.3 The corrected verdicts

Each null replicate is now accepted only if **its own realised divergence** is
within 10% of the population's:

| population | margin | published floor | ratio | corrected floor | ratio | acceptance | literature |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| coxsackievirus A6 | 329/7312 | 268/7390 | 1.24× | **121/7390** | **2.75×** | 5/30 | reticulate |
| hepatitis B | 9/3182 | 2/3182 | 4.50× | 2/3182 | 4.50× | 5/23 | reticulate |
| JC polyomavirus | 4/5108 | 24/5121 | 0.17× | **5/5121** | **0.80×** | 5/16 | tree-like |

**No verdict is overturned**, and all three still agree with their literature
expectation (3/3). But two published numbers move, in opposite directions:

* **Coxsackievirus A6 strengthens.** 1.24× → 2.75×, so it now clears *twice* its
  floor as well. The published result understated its own evidence, and the
  summary's `exceeding_twice_matched_floor` gains it.
* **JC polyomavirus's negative is far less comfortable than published.** At
  0.17× it read as decisively tree-like; against a properly matched floor it is
  a **near-miss at 0.80×**. The verdict is unchanged and still correct, but the
  published figure overstated how far from the line it sat. This is the honest
  cost of the correction and it is not visible anywhere in the verdict itself.

Hepatitis B is untouched: its floor is identical rational (`2/3182`), which is
the stability a population well inside its regime should show.

### 18.4 What this does and does not change

The method's published *conclusions* survive — three usable populations, seven
refused on admissibility, 3/3 agreement with the literature. What changes is
that the floors those conclusions rest on are now **matched per replicate rather
than in expectation**, with acceptance rates reported (17%, 22%, 31%) and a
population that cannot be matched receiving no floor and no comparison instead
of a floor from an unmatched null.

The transferable point, restated because getting it wrong twice is the evidence
for it: **when two defects in an apparatus act in opposite directions, the net
direction is not an argument, it is a measurement.**

Raw output: [`real_data.json`](../data/reticulata_phase3/real_data.json),
regenerated — never hand-edited — and
[`tuning_bias.json`](../data/reticulata_phase5/tuning_bias.json).

## 19. The split-level certificate does not transfer to the real populations

§17 ended with the split-level certificate separating reticulation from relaxed
clocks, gap blocks and saturation on **35 of 35** simulated dataset-seeds, and
with the scope stated as "one parameter grid, 8–10 taxa". That scope statement
named the wrong axis as harmless. Taxon count is not a detail of the grid; it is
the variable the method's central requirement scales with. On the real
populations, at 16 taxa, the split-level certificate **has no power at all**, and
the extension that recovers power **does not separate**.

### 19.1 Making it computable: one table instead of 39.7 million scans

Certifying every split of 16 taxa means evaluating about **39.7 million**
straddling quadruples. There are only `C(16,4) = 1820` distinct ones, each reused
about **21 800** times, so the quadruple counts are precomputed once
(`QuadrupleTable`, ~1.3 × 10⁷ column comparisons) and the rest are lookups. Done
naively it is ~2.8 × 10¹¹.

The table identifies which of the three canonical pairing sums a split induces by
asking whom the smallest index is paired with — one piece of index arithmetic
that is easy to get subtly wrong and invisible in aggregate output. It is checked
exhaustively against directly computed counts for **every quadruple and every
2+2 split**, and the whole refactor is confirmed behaviour-preserving: the
simulated experiment reproduces its committed artifact exactly, all 35 matched
rows and all 350 floor-sweep rows identical.

### 19.2 Zero power at sixteen taxa

A 16-taxon tree has `n-3 = 13` internal splits. Against a matched-null floor,
with the strict conjunction:

| population | prior | certified splits | reads as |
| :--- | :--- | ---: | :--- |
| coxsackievirus A6 | reticulate | **0** / 32751 | tree |
| hepatitis B | reticulate | **2** / 32751 | tree |
| JC polyomavirus | tree-like | **4** / 32751 | tree |

Every population reads as a tree, so the single "agreement" with the literature
is vacuous — a method that always answers *tree* agrees with every tree-like
prior for no reason.

**The mechanism, measured rather than inferred.** The conjunction requires every
straddling quadruple to clear the floor, and the number of conjuncts for a
balanced split is `C(n/2,2)²`: **36 at n = 8, but 784 at n = 16**. The
near-miss histogram shows the consequence precisely — the closest-missing splits
failed by **1 conjunct out of 91** (coxsackievirus), **2 of 91** (hepatitis B)
and **11 of 234** (JC polyomavirus). The rule is not wrong about which splits
exist. It is too strict to certify one at that conjunction length, and it got
away with it on the simulated grid only because 8 taxa made the requirement
roughly twenty times weaker.

### 19.3 An exception budget, calibrated and not chosen

Allowing at most `k` dissenting conjuncts keeps the condition **exact and
countable** — "at most `k` of these integer comparisons fail" is as checkable as
"none does" — so it is a weaker claim, not a vaguer one. What it must never be is
a *chosen* number. `k` is calibrated as the largest value at which **every
matched null replicate still reads as a tree**: one step beyond it, a genuine
clonal tree at this population's own divergence would be called reticulate.

Calibrated `k = 23` for all three populations, with the null first breaking at
24. That recovers power, and produces a false positive:

| population | prior | certified at `k=23` | reads as | verdict |
| :--- | :--- | ---: | :--- | :--- |
| coxsackievirus A6 | reticulate | 8 | non-tree | correct |
| hepatitis B | reticulate | 18 | non-tree | correct |
| JC polyomavirus | **tree-like** | 6 | **non-tree** | **wrong** |

### 19.4 The comparison that replaces the threshold, and why it still fails

A threshold invites fitting, so the better statistic is a **comparison**: the
budget at which the population stops reading as a tree, against the budget at
which its own matched null does. No cutoff is chosen anywhere.

| population | prior | population breaks at | null breaks at | difference |
| :--- | :--- | ---: | ---: | ---: |
| hepatitis B | reticulate | 3 | 24 | −21 |
| coxsackievirus A6 | reticulate | 9 | 24 | −15 |
| JC polyomavirus | **tree-like** | 13 | 24 | **−11** |

**All three break earlier than a tree, including the negative control.** The
ordering is in the right direction — the two reticulate populations break
earliest, the tree-like one is intermediate, a genuine tree latest — but the
statistic does not separate. Any threshold between 13 and 24 would classify all
three correctly, and with three populations and **one** negative control,
choosing such a threshold is fitting it to the answer. That is the forced test
this document has had to label in its own apparatus more than once, and it is not
being done here.

Both alignment references agree on every verdict (3/3), so none of this is an
alignment artefact.

### 19.5 The honest verdict, and the binding constraint

The split-level certificate is **proved** (§17.2, four theorems, module 27/27
with 0 `sorry`) and **unvalidated**. It works on one simulated grid at 8–10 taxa
and does not transfer to 16 real taxa, where the strict form has no power and the
calibrated form cannot be shown to separate.

What limits this is not the statistic. It is that **the aligner admissibility
gate leaves only 3 of 10 populations usable, and only one of them is a negative
control.** No amount of further statistical work on three populations can
establish specificity; one negative control cannot distinguish "this method has a
false-positive problem" from "this population's tree-like prior is wrong" — and
§12's own caution already records that one of the negative-control families has
documented recombination, so these priors are priors, not labels.

So the next move for this line is **more admissible negative controls**, by
better alignment or by pre-aligned data of known provenance, and not another
statistic. Until then the quartet-level certificate remains the programme's only
validated instrument, and the split-level lift remains a proved theorem with a
simulated demonstration.

Raw output: [`split_real.json`](../data/reticulata_phase5/split_real.json).

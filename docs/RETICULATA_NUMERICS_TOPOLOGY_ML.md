# Reticulata: numerics, topology and where machine learning is allowed

Status: proposal, 2026-10-22. Written after the Phase 2/3 work and the
seven-thread literature review, and against what that review found. Nothing here
is a result; §1 states the rule the rest of the document obeys, §2–§4 are ranked
directions, and each carries a first step and a criterion for abandoning it.

The companion documents are [`GENETICA_RETICULATA.md`](./GENETICA_RETICULATA.md)
(what exists) and `reports/Reticulata next research directions.md` in the
workspace root (what the literature says).

---

## 1. The rule: where machine learning may enter, and where it may not

This method's entire claim is that its conclusions are *checkable*. The null is a
theorem rather than an ensemble; the certificate is a witness a kernel
recomputes; the aligner's contribution is measured rather than assumed. A
learned component can be added to that without damage **only if a wrong model
costs time or power, never correctness**. That gives a sharp admissibility test,
and it rules more out than in.

**Admissible.** A learned component whose output is either (a) validated against
a cheaper ground truth that we can generate at will, or (b) incapable of changing
a verdict even when wrong.

**Inadmissible, and why:**

| proposal | why not |
| :--- | :--- |
| A learned aligner, or learned alignment refinement | The aligner is already the largest item on the trust surface, and §11.3 of the companion document had to measure its error to gate on it. A learned aligner's error is measurable in the same way, but it would be *harder* to characterise across the parameter space, and the gate would then rest on a model's generalisation. This makes the worst-covered step worse-covered. |
| A learned statistic reported as the primary detector | It has no exhibitable witness. The one thing this method has that the incumbents do not, per the head-to-head, is that a reader can re-verify a finding without running our code. A learned score discards exactly that. |
| ML-chosen thresholds without held-out validation | The admissibility threshold is currently *read off* a measured dose–response curve. Fitting it is fine; fitting it without a held-out parameter region is how a threshold stops being a measurement. |
| Learned imputation of gappy columns | Complete deletion within the subset exists because differing denominators manufacture certificates. Imputation would reintroduce that confound in a form no longer visible in the counts. |

Everything in §4 is designed to pass the admissibility test rather than to argue
with it.

---

## 2. Numeric Rust

### 2.1 Bit-parallel subset counting · **done, 21.6× measured**

The certificate's statistic is a maximum over subsets, and the head-to-head
showed its real defect is not power but **incomparability across taxon counts**:
the quartet floor moves from 5–7 at 16 taxa to 9–10 at 32. That defect can only
be characterised by measuring it at several taxon counts, and the current counting
is byte-at-a-time over columns — `Alignment::quadruple_counts` walks the
alignment once per quadruple, so the work is `O(n⁴ · L)` and 32 taxa already
costs 35 960 quadruples × 4000 columns.

`incompat.rs` shows the fix, because it already does this: pack one bit per
sequence per column and the four-gamete test becomes four bitwise ops. For
distances the equivalent is a **2-bit-per-base transposed layout** — for each
column, two `u64` words holding the base's two bits across up to 64 sequences.
Then the differing-site count for a pair is a `popcount` over `(a0^b0) | (a1^b1)`
masked by the comparable-site mask, and a quadruple's six counts come from six
such reductions over `L/64` words. That is a 20–60× reduction in inner-loop work
on plausible parameters, which is what makes 64 and 128 taxa measurable.

**Implemented as `bitcount::BitAlignment`.** Measured **21.6×** at thirty-two
taxa (0.24 s against 5.14 s over 35 960 quadruples) with both implementations
returning the same exact margin, which the `scale` command re-checks on every run
and exits non-zero if it ever fails. The byte path remains the specification and
the agreement is a property test across widths 1, 7, 63, 64, 65, 130 and 257 plus
gap-blocked simulated data.

What it bought, immediately: the taxon-count sweep of §13 in the companion
document. The floor roughly doubles from sixteen to sixty-four taxa while the
mosaic margin stays flat, so separation halves from 19.4× to 9.2× — the defect
measured across a range rather than argued from two points. Detection at a fixed
size survives; a threshold calibrated at one size does not transfer.
* **Note.** `popcount` is a single instruction on anything this will run on;
  explicit SIMD is not needed and would add an unsafe surface for little gain.

### 2.2 Exactness guards at scale

Margins are exact rationals and comparisons are integer cross-products. At 128
taxa with `L = 30 000` the products stay far inside `i64`, but the
cross-multiplication in `exceeds_k_times` already uses `i128` and the subset
counting does not need to. **First step:** a property test asserting no overflow
at the largest parameters the search will accept, so the bound is checked rather
than reasoned about. Cheap, and it is the kind of thing that fails silently.

### 2.3 The aligner, if it is kept at all

The literature review's verdict was that reference-anchored projection is
standard surveillance practice and the final step of MUSCLE5's `Super5`, so the
construction is defensible; and the Phase 4 direction (§3.2) removes our aligner
from the critical path for SARS-CoV-2 entirely. If it is kept for other families:

* The cheap optimality certificate of Gibrat — attaining an upper bound on any
  path leaving the band certifies optimality in constant time — would replace the
  current full unbanded recomputation. His derivation assumes a band about the
  main diagonal while ours is about the proportional diagonal, so **the
  derivation must be redone, not transplanted.** That is a small piece of real
  mathematics and it is worth doing: it turns one extra `O(nm)` pass per pair
  into `O(1)`.
* `poasta` (BSD-3, pure Rust, peer-reviewed 2025) implements exact gap-affine
  partial-order alignment with a proved-admissible A\* heuristic. That is the
  route to a non-reference-anchored aligner without writing a progressive one,
  and it is the single biggest available reduction in the trust surface for
  gappy, divergent families. **Abandon if** its divergence range turns out not to
  cover the families of interest.

---

## 3. Topology

### 3.1 Is persistent homology specific, or does it track mutation? · **built, not yet run**

The construction this programme started from reads an unfilled 1-cycle in the
Vietoris–Rips complex of a genetic-distance point cloud as reticulate evolution.
The theoretical warrant is the same fact the four-point null rests on: a tree
metric admits no unfilled cycle. What the literature review did **not** find is
any peer-reviewed independent calibration of that method's specificity; the one
unpublished evaluation it found reports the irreducible-cycle rate rising with
the point-mutation rate rather than with the recombination rate.

That is testable here with a known answer, because the simulation grid has the
ground truth: three datasets contain reticulation and four do not, and the
saturated dataset differs from the baseline *only* in having ten times the
substitution rate. `reticulation topology` is implemented and calibrates
`h1_total_persistence` by parametric bootstrap under a matched tree null — the
same calibration the other three detectors get, which is the point: Section 8 of
the textbook reported a persistence value against a shuffle null and the result
was a negative finding obtained the hard way.

**Run, and inconclusive — for a diagnosable reason.** The output reports
`"this_experiment_is_interpretable": false`. The matched null produced no cycle
in **495 of 495 replicates** on five of seven datasets, so the p-value collapses
to two values: the floor whenever the observed value is non-zero, 1.000
otherwise. The bootstrap calibrates nothing and the test reduces to "did a cycle
appear at all", so the 2/5 counts on the reticulate datasets *and* on the
saturated one are that binary and cannot be read as power or specificity. The
concern is neither reproduced nor refuted. See §12.4 of the companion document.

* **Diagnosed cause.** At 16 taxa with a Rips threshold at the median
  off-diagonal distance, the complex is too small for H₁ to populate.
* **Re-run at 64 taxa, and the answer arrived.** §2.1 made it affordable. Five
  of seven nulls are still degenerate, so the run is still formally flagged
  uninterpretable — **but the two that are not degenerate are the two datasets
  containing no reticulation**, and the statistic fires on both at 3/3:
  saturated data (ten times the divergence) and gap-blocked data. Meanwhile a
  genuine single recombinant produces *no cycles at all* and is missed 0/3. The
  magnitude hypothesis below also fails: the mosaic persistence is within a
  factor of 1.7 of both false positives. See §13.3 of the companion document for
  the scope caveats, which matter — this is a result about this construction and
  calibration at three seeds, not a refutation of the published method.
* **And one structural argument, which is stronger than the counts.** The Rips
  construction needs a *global* distance matrix, so it cannot use complete
  deletion within a subset — there is no "within this subset" when the complex is
  built from all pairs at once. The topological route therefore necessarily uses
  pairwise deletion and necessarily inherits the audit's finding 6, the confound
  the four-point certificate was redesigned to avoid. The gap-blocked dataset's
  13.3 cycles per population is that confound appearing as topology. This follows
  from the construction, not from three seeds.
* **Hypothesis tested and failed, kept rather than deleted.** At 16 taxa the
  mosaic persistence looked an order of magnitude larger than the
  single-recombinant and saturated ones, suggesting magnitude might separate
  where the binary verdict could not. At 64 taxa it does not.

### 3.2 A four-way benchmark on one ground-truth grid

The field, per the review, has no benchmark that puts these detectors on the same
simulated data with the same calibration. We are two statistics away from one:

| detector | status | calibration |
| :--- | :--- | :--- |
| quartet certificate (max margin) | exists | matched tree null |
| pairwise homoplasy index Φ_w | exists | site-order permutation |
| mean δ_q, calibrated | exists | matched tree null, parametric bootstrap |
| H₁ total persistence | implemented, §3.1 | matched tree null, parametric bootstrap |

With the matched-null column added to all four, the table is commensurable and
reports specificity on four no-reticulation datasets and power on three. The
honest expectation, after the head-to-head, is that **the certificate will not
win on power** and its contribution will remain the exhibitable witness. Finding
that out precisely is worth more than assuming it either way.

### 3.3 A filtration built from quadruple violations · **demoted; see §14 of the companion document**

The filtration funnel closed the Rips line, and closed this with it. A complex
over violating subsets would use per-subset complete deletion and so would not
inherit finding 6 — but it adds information only about how violations *overlap*,
and the certificate already captures the quadruple-level signal the global
construction was missing. Optional and speculative, not the pivot.

Every statistic above summarises the quartet distribution by one number. The
quadruples themselves carry structure: *which* quadruples violate, and how they
overlap. A complex whose simplices are the violating subsets, filtered by margin,
would have homology that says something about how the violations are arranged —
one reticulation event should produce violations sharing the two parents, where
population structure should produce a different pattern.

This is a guess, and it is flagged as one. **First step:** build the complex on
simulated data with one planted recombinant versus two, and check whether any
homological invariant separates them; if a simple count of shared parents does
the same job, prefer the count. **Abandon if** the invariant needs the margin
threshold tuned per dataset to show anything, which would make it a parameter
search dressed as topology.

---

## 4. Machine learning, under the §1 rule

### 4.1 A learned surrogate for the null floor · **built; restricted deployment**

Phase 3 of the original plan — the margin's null distribution — remains open, and
the review confirmed that no analytic sampling distribution of the four-point
margin exists in the literature either. Current practice here is to *simulate* a
matched null per query: the `real-data` command spends a 40-step bisection plus
five tree simulations per population, and the calibrated statistics above spend
99 replicates per dataset-seed. That is thousands of simulations to answer one
question, and it makes large sweeps impractical.

The floor is a smooth function of a handful of parameters: alignment width,
mean p-distance, taxon count, gap load, clock spread. So fit it. Simulate across
that space once, regress the null margin's upper quantiles on the parameters, and
calibration becomes a lookup.

**Why this passes the §1 rule.** The surrogate predicts a quantity we can always
generate the truth for, so it is validated by held-out simulation at measured
error — and, critically, the certificate is untouched: it stays an exact integer
comparison, and only the *floor* it is compared against is predicted. A
deployment can even keep the surrogate honest at runtime by simulating a handful
of replicates and checking they fall where the surrogate says.

* **First step.** A parameter sweep writing `(parameters → empirical quantiles)`
  to a JSON table; then the simplest model that fits — a monotone regression or a
  small gradient-boosted tree — with a held-out *region* of parameter space, not
  a random split, because the failure mode that matters is extrapolation.
* **Report** predicted-versus-simulated floor error as a distribution, and refuse
  the surrogate outside the convex hull of the training parameters rather than
  extrapolating. A surrogate that silently extrapolates is exactly the
  unexamined step this project keeps finding in its own apparatus.
* **Abandon if** held-out error at the quantiles that matter exceeds the margin
  differences the method is trying to resolve, in which case it cannot support a
  verdict and simulation stays.

**Built, and the criterion is partially met** (§16 of the companion document).
The falsifiable check passed: the fitted length exponent is **−0.556** against a
prediction of −0.5 fixed in advance, so the floor is noise-dominated and the law
is confirmed rather than interpolated. Held-out *regions* give median relative
errors of **20.5%** (longest alignments) and **28.7%** (largest taxon count),
p90 up to 66%. By the criterion above that is enough for wide-margin verdicts
(Hepatitis B clears its floor by 4.5×) and **not** enough for borderline ones
(Coxsackievirus clears by 1.23×, which a 25% floor error flips). So: adequate for
screening and for predicting degeneracy, not for deciding near the floor, which
still gets a simulation.

What it actually buys is the boundary rather than the prediction. Its calibrated
divergence range is 0.00239–0.51429, and the XCP published-breakpoint population
sits at 0.00057 — 4.2× below the lower bound, so the surrogate refuses it and
flags it as out of regime. That is the conclusion §12.3 reached only after a
2.2 GB download and a degenerate floor.

### 4.2 A learned search order · **admissible, cannot corrupt the answer**

The quartet search is exhaustive and refuses to sample, which is correct and is
why 32 taxa needs an explicit ceiling. A model that *ranks* quadruples by
predicted margin would let the search find its strongest witness early and report
a bound on what remains unexamined — without ever changing which quadruples are
admissible.

This is the ideal place for a learned component: **a wrong model costs time, not
correctness.** The certificate is still exact, the enumeration is still complete
unless a ceiling is hit, and the model's only effect is the order of discovery.

* **First step.** Cheap features per quadruple — pairwise distance ranks, shared
  informative sites — and measure *rank correlation* with the true margin on
  simulated data. If a hand-written heuristic on pairwise distances gets most of
  the achievable ordering, use that and drop the model.
* **Abandon if** the ordering gain does not change the taxon count that is
  practically reachable, since that is the only reason to want it.

### 4.3 A learned detector as a ceiling, not a component · **admissible, diagnostic**

Train a classifier to separate reticulate from non-reticulate populations on the
simulated grid, with held-out seeds and held-out parameter regions. Do not deploy
it. Use it to answer one question the four hand-built statistics cannot answer
about themselves: **how much detectable signal are they leaving on the table?**

* If a trained model substantially beats all four on held-out simulations, that
  bounds the gap and tells us the hand-built statistics are lossy — useful, and it
  motivates looking for what feature the model is using.
* If it does not, that is equally informative: it suggests the statistics are
  near the achievable limit on this generative model, and that further work
  belongs on calibration and provenance rather than on detection.
* **The trap to avoid**, and it is the project's own recurring one: a model
  trained and tested on the same generative process measures the process, not the
  biology. Any such number must be reported as "on this simulator", and the real
  test is held-out *parameter regions*, not held-out seeds.

### 4.4 What a learned component must never be allowed to do here

Produce the verdict without a witness. The head-to-head established that this
method does not win on power; its distinctive contribution is that a reader can
re-verify a finding from four sequences and six integers. Every direction above
preserves that, and any direction that does not should be refused on those
grounds rather than evaluated on accuracy.

---

## 5. Ranking

Corrected after running §3.1 and Phase 4, both of which hit the same wall: a
**degenerate null**. The topology bootstrap produced no cycle in 495 of 495
replicates, and the matched floor for one SARS-CoV-2 population produced no
violating quadruple at all. In both cases the calibration collapsed to "is the
observed value non-zero", and in both cases the cause is that the regime is too
small or too similar for the statistic to have anything to measure. That moves
scale from an optimisation to a precondition.

| # | direction | effort | why first |
| :-: | :--- | :--- | :--- |
| 1 | §2.1 bit-parallel counting | days | **promoted.** Larger taxon counts are the precondition for §3.1 to be interpretable at all, and the only way to characterise the max statistic's real weakness |
| 2 | ~~§4.1 learned null floor surrogate~~ **built** | — | Done, with a restriction: the falsifiable exponent check passed (−0.556 against −0.5), held-out error is 20–29% median, so it screens and predicts degeneracy but does not decide borderline verdicts. It does reproduce §12.3's conclusion in advance by refusing the out-of-regime population |
| 3 | ~~topology~~ **closed** | — | **The line is closed (§14 of the companion document).** No cheap filtration change met the criterion; F1 (saturation) was refuted, F2 confirmed that the gap false positive is finding 6 reappearing as topology, and a single recombinant produces zero cycles in all four variants while the certificate detects it 3/3 at matched taxon count. Earlier note retained: **§3.1 is answered as far as this construction goes.** Where its calibration is valid the statistic tracks divergence and missing data rather than reticulation, and the structural argument above says the Rips route cannot adopt complete deletion. Further topology work should either change the filtration (§3.3) or stop; re-running the same construction at more taxa will not change the structural point |
| 4 | §3.2 complete the four-way benchmark | hours after 3 | the field has no such table |
| 5 | §4.3 learned ceiling | days | tells us whether to keep working on detection at all |
| 6 | §2.3 `poasta`, or Gibrat's bound rederived | days | reduces the largest trust-surface item for gappy families |
| 7 | §3.3 violation-complex homology | unknown | speculative, and may reduce to counting shared parents |

**A lesson worth carrying into every item above.** Two experiments this round
produced a confident-looking count from a null that could not vary. Both now
compute and report a degeneracy flag, and any new calibrated statistic should do
the same before its numbers are read: *check that the null varies before
believing a p-value derived from it.* A null that is identically zero makes every
positive observation significant at the floor, which is indistinguishable from a
powerful test until you look.

Two things this list deliberately does not contain: a learned aligner, and a
learned primary statistic. §1 says why.

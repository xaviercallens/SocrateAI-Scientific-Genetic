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

**What a Rust bug can do:** miss a certificate, or emit a triple that is not one.
**What a Rust bug cannot do:** produce a false conclusion that survives checking.

Compare §11.3 of the textbook, where a wrong subsampling stride silently produced a significant result that only a lucky nested-null design caught. Under this architecture that class of failure is not merely unlikely; it is structurally excluded, because the conclusion is carried by the certificate and not by the pipeline that found it.

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

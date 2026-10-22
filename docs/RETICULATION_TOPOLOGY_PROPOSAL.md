# An alternative to "TDA doesn't work on DNA"

**Status:** implemented, controls passing, not yet run on real data.
**Lean:** `formal/ANSE/ReticulationCertificate.lean` — 10/10 declarations clean, no `sorryAx`.
**Rust:** `rust/reticulation/` — 25 tests passing, extends [lophat](https://github.com/tomchaplin/lophat).
**Controls:** `data/reticulation_controls.json`.

---

## 1. What actually failed, and why that matters

Section 8 of the textbook embedded **one** DNA sequence by the chaos game and measured the topology of the resulting point cloud. Section 11 then ran the shuffle null and found nothing: three nested nulls, 200 shuffles, three genes, no signal beyond base composition.

It is tempting to summarise that as "TDA doesn't work on DNA". That is the wrong lesson, and Worked Example 11.1 says so. The failure was specific and diagnosable:

| What went wrong | Consequence |
| :--- | :--- |
| The chaos game contracts *any* sequence onto a self-similar attractor | The holes measured were the attractor's, not the sequence's |
| One sequence, so no population structure to measure | Nothing for a loop to *mean* |
| No mathematical null available | Had to calibrate against shuffles |
| The shuffle null had no proved semantics | p-value floor, frame-stride artefact (§11.3), statistic with no interpretation |

Every one of those is a property of the *construction*, not of persistent homology. So the question is not whether to abandon topology but whether there is a construction in which a loop has a meaning. There is, it is well established, and we had simply picked the wrong one.

## 2. Literature: where topology on sequence data genuinely works

The canonical positive result is **[Chan, Carlsson & Rabadan, PNAS 110(46):18566–18571 (2013)](https://doi.org/10.1073/pnas.1313480110)**, "Topology of viral evolution". The move that makes it work is to compute persistence on a **population** of sequences under a genetic-distance metric, not on an embedding of a single sequence. Then:

- `H₀` ↔ clustering / population structure.
- **`H₁` ≠ 0 ↔ non-tree-like (reticulate) evolution** — recombination, reassortment, horizontal transfer.

The reason this is not another heuristic is structural: *a tree has no loops*. Phylogenetic methods commit to one tree; persistence instead asks, across all scales at once, whether the data is consistent with **any** tree. The authors state the dichotomy directly — clonal evolution can be summarised as a tree, while **reticulate evolution exhibits nontrivial topology of dimension greater than zero**. They applied it to dengue, HIV and influenza, and bounded the rate of reticulate events at roughly 20 reassortments per year in avian influenza.

Follow-up theory matters for honesty about what the statistic measures. Lesnick and collaborators introduced the **novelty profile**, a stable statistic that "not only counts recombination events but also quantifies how recombination creates genetic diversity" — which exists because raw Betti numbers are *not* stable summaries. And a student replication study ([UC Davis REU](https://reu.math.ucdavis.edu/application/files/8416/4826/7705/geisler-report.pdf)) found persistent homology "variable in its detection of recombination", with some homology measures preferable to others. We take both as warnings rather than discouragement: the construction is right, the summary statistic needs care.

On the mathematical side, the relevant fact is about **tree-likeness**. A tree metric is 0-hyperbolic, and for finite tree metrics the Vietoris–Rips complexes collapse onto their corresponding subforests ([Bauer & Roll, arXiv:2112.06781](https://arxiv.org/abs/2112.06781)) — hence acyclic. We searched for a citable statement of the sharper claim we wanted (that an **ultrametric** has disjoint-clique Rips complexes at every scale) and did not find one, so we proved it ourselves rather than assert it. That turned out to be the most useful decision in the whole design.

Supporting references: **Gromov hyperbolicity and Rips collapses** (above); **Altschul & Erickson (1985)**, [doi:10.1093/oxfordjournals.molbev.a040370](https://doi.org/10.1093/oxfordjournals.molbev.a040370), for sequence nulls that state what they preserve; **Gotoh (1982)**, [doi:10.1016/0022-2836(82)90398-9](https://doi.org/10.1016/0022-2836(82)90398-9), for the alignment layer.

## 3. The proposal: a null that is a theorem, not an ensemble

Here is the part that is new, and it is the direct answer to what went wrong in §2 of this document.

An **ultrametric** satisfies `d(x,z) ≤ max(d(x,y), d(y,z))`. Biologically that is exactly a *strict molecular clock on a clonal genealogy*: all distances are set by divergence times, so in every triple the two largest distances agree.

> **No-go theorem.** If `d` is an ultrametric then for every scale `ε` the relation `d(·,·) ≤ ε` is **transitive**. Hence every Vietoris–Rips component is a clique, every cycle has all its chords, and `β₁ = 0` at every scale simultaneously.

The proof is three lines and it is in Lean (`ultrametric_no_induced_path`). The consequence is the whole point:

**An unfilled triangle is not "surprising relative to shuffles". It is a certificate that the data is not ultrametric.** No shuffle ensemble, no p-value floor, no stride artefact.

### The certificate, and why it changes the trust model

For an unordered triple, sort the distances `d_min ≤ d_mid ≤ d_max`. The violation is

```
margin = d_max − d_mid > 0
```

and the matrix is an ultrametric exactly when no triple has a positive margin (`ultrametric_iff_no_positive_margin`). The statistic is **exact and parameter-free** — no scale, no threshold, nothing fitted.

A violating triple maps precisely onto the Lean `Certificate` structure: `x, z` realise `d_max`, `y` is the third point, `scale = d_mid`. Its three fields are the three hypotheses of the theorem.

This gives an architecture worth stating explicitly:

```
Rust searches  →  emits a certificate (3 labels, 3 distances)
Lean proves    →  certificate ⟹ ¬ultrametric   (certificate_refutes_ultrametric)
Human rechecks →  three inequalities, by hand
```

**The scientific conclusion does not depend on the correctness of the Rust.** A bug can fail to find a certificate, or emit a triple that is not one — and the latter is caught by `Certificate::validate`, which re-derives all three inequalities before anything is reported. What a bug *cannot* do is produce a false conclusion that survives checking. Compare the chaos-game pipeline, where a wrong subsampling stride silently produced a significant result.

### Two nested nulls, both mathematical

| Null | Biology | Statistic | Status |
| :--- | :--- | :--- | :--- |
| **ultrametric** | strict clock + clonal | `max(d_max − d_mid)` over triples | **proved** in Lean |
| **four-point / tree** | clonal, any branch rates | `max(S_max − S_mid)` over quadruples | **cited**, not proved here |

Reading them together separates *rate variation* (ultrametric fails, four-point holds) from *reticulation* (both fail). This is the nested-null idea of §11.2 — but with proved nulls instead of shuffles. The asymmetry in status is deliberate and the code enforces it: four-point defects are labelled `suggestive`, never `proved`.

## 4. The honest limitation, stated first rather than last

**Finite sequences are never exactly ultrametric.** Distances are stochastic, so clonal data violates the null too. In our negative control:

> **539 of 560 triples "violate" ultrametricity in purely clonal data.**

Every one of those certificates is *valid* — the data really is not an ultrametric. The theorem is correct. It simply proves something **weaker than the interesting claim**: it proves "not a strict clock on a clonal genealogy", and multinomial sampling at finite sequence length is enough to break that.

So the existence of a certificate is not evidence of recombination. Only the **margin magnitude** discriminates, and that still requires calibration.

What changed is the *kind* of calibration available: a **mechanistic generative model** of clonal evolution at matched sequence length and divergence, rather than a permutation of the observed data whose meaning was unclear. That is a real improvement in interpretability. It is not an escape from calibration and must not be described as one.

This is the same shape as the book's organising lesson: a true theorem whose hypotheses are not quite the biological claim you wanted.

## 5. Implementation

### Repository selected to clone and extend

**[tomchaplin/lophat](https://github.com/tomchaplin/lophat)** (MIT, v0.14.1, edition 2024) — "Lockfree Persistent Homology Algorithm Toolbox", implementing the parallel lock-free reduction of Morozov & Nigmetov, with `R = DV` decomposition and **representative cycles**.

Why this one, having looked at the alternatives:

| Candidate | Verdict |
| :--- | :--- |
| **lophat** | **Selected.** Maintained, edition 2024, rayon, exposes `diagram_with_reps` |
| `teia` | Stale (~7 years); its README lists Vietoris–Rips as future work |
| `rust-ripser` | Only surfaced as a Python wheel; no crate or repo found |
| Ripser (C++) | Reference implementation, but not Rust; kept as a validation oracle |
| GUDHI (C++/Python) | Already cloned at `/home/xavkal/xdev/gudhi-devel`; kept as the second oracle |

The gap to fill was identified from lophat's own README: *"LoPHAT is **not** specialised to compute PH of common filtrations or even filtered simplicial complexes."* It consumes a filtered boundary matrix. So the extension is exactly the **Vietoris–Rips front-end over a genetic-distance matrix** — which is upstreamable as a PR rather than requiring a fork. The clone sits at `/home/xavkal/xdev/temp/oss-clones/lophat` for reference; the crate depends on the published `lophat = "0.14"` so builds stay reproducible.

### What was built — `rust/reticulation/`

| Module | Contents |
| :--- | :--- |
| `distance.rs` | p-distance and Jukes–Cantor over an alignment, rayon-parallel; pairwise gap/ambiguity skipping; **saturation flagged, never silently reported as identity** |
| `certificate.rs` | exhaustive `O(n³)` triple search, exact defect, `validate()` re-deriving every inequality; sampled four-point defect labelled non-exhaustive |
| `rips.rs` | **the lophat extension**: filtered VR boundary matrix to dimension 2, faces before cofaces, → lophat → `H₁` intervals *with representative cycles mapped back to sequence labels* |
| `simulate.rs` | SplitMix64 (written out, so a seed reproduces a run with no crate-version dependence); clonal-tree and recombinant generators for the controls |
| `main.rs` | `controls` and `fasta` subcommands; **exits non-zero if either control fails** |

Leverage of the existing machine: conventions, edition and rayon/faer usage follow `/home/xavkal/xdev/rusty-SUNDIALS` (a cargo workspace at edition 2024); GUDHI at `/home/xavkal/xdev/gudhi-devel` is the cross-validation oracle. The crate is deliberately *not* added to the SUNDIALS workspace — a genetics/TDA crate is a poor architectural fit inside an ODE-solver workspace, and the science belongs with the textbook.

### Verification actually run

| Check | Result |
| :--- | :--- |
| `cargo test --release` | **25 passed** |
| Lean gate on `ReticulationCertificate.lean` | **10/10 clean**, no `sorryAx`, axioms ⊆ `{propext, Classical.choice, Quot.sound}` |
| `reticulation controls` | exit 0, both controls pass |
| `ultrametric_has_no_finite_h1_class` | computational counterpart of the no-go: **no finite `H₁`** |
| `square_cycle_has_a_finite_h1_class` | birth 1.0, death 2.0 — exactly as predicted analytically |
| `validate_rejects_a_forged_certificate` | falsifiability probe: a forged certificate **is** rejected |

### Control results (`data/reticulation_controls.json`)

16 clonal leaves, 4000 nt, μ = 0.004/generation, seed 2026.

| | Negative (clonal) | Positive (+1 recombinant) |
| :--- | ---: | ---: |
| ultrametric defect | 0.00450 | **0.01525** |
| ratio to clonal floor | — | **3.39×** |
| violating triples | 539 / 560 | — |
| recombinant is a certificate midpoint | — | **yes** |
| all certificates validated | — | yes |

The detail worth more than the ratio: the **second-ranked certificate names the true parents exactly** — `clonal_003 — recombinant — clonal_013`, where the recombinant was constructed from `clonal_003` and `clonal_013`. The method does not merely flag that something non-clonal happened; it identifies the sequences involved, which is what the representative cycle is for and what no global p-value can do.

## 6. What has *not* been done

- **No real data.** The controls are simulated. The method needs an aligned *population*; the textbook's `POLE`/`POLD1`/`MSH6` are one sequence each and not mutually homologous, so they are not a valid input. This is the immediate next step, not a completed result.
- **No GUDHI cross-validation run yet.** The oracle is on the machine and the test cases are analytic; a diff against GUDHI on random point clouds is the obvious next check.
- **The tree-level no-go is cited, not proved.** Only the ultrametric case is in Lean.
- **The margin's null distribution is not characterised analytically.** We calibrate by simulation. A closed form for the clonal margin distribution as a function of sequence length and divergence would remove the simulation entirely, and is the most valuable piece of theory left.
- **No stability analysis.** Betti numbers are unstable summaries; the novelty profile exists for that reason. We have not implemented it.

## 7. Next steps, ranked

1. **Run it on an aligned viral population** where recombination is independently known — the HIV and influenza datasets of Chan et al. are the natural benchmark, since the answer is published. Concrete first action: obtain an aligned HIV-1 *env* set and compare recovered certificate midpoints against known circulating recombinant forms.
2. **Cross-validate `rips.rs` against GUDHI** on random clouds; any disagreement in `H₁` intervals is a bug in the new front-end.
3. **Characterise the clonal margin distribution** in sequence length and divergence. This is the piece that would turn calibration-by-simulation into a closed-form threshold, and it is a tractable probability problem.
4. **Upstream the VR front-end to lophat** as a PR. It fills a gap the maintainer has documented, and the project benefits from review by someone who knows the reduction.
5. **Prove the tree-metric no-go in Lean**, lifting the four-point statistic from *suggestive* to *proved* and making the nested nulls symmetric in status.
6. **Implement the novelty profile** as a stable alternative to raw Betti numbers.

## 8. Why this is a better bet than the approach it replaces

Not because topology is fashionable, but because each specific failure has a specific fix:

| §8's failure | Fixed by |
| :--- | :--- |
| topology of the embedding, not the data | population under a genetic-distance metric |
| loops had no meaning | a loop means non-tree-like evolution |
| null was a shuffle ensemble | null is a theorem |
| statistic had no interpretation | exact, parameter-free defect with a proved characterisation |
| output was a global p-value | output is a certificate naming specific sequences |
| a pipeline bug produced a significant result | certificate re-validated; conclusion independent of the search code |

And one property the old approach could not have had at all: if the Rust is wrong, the conclusion still cannot be. That is the only form of software trust that scales to code nobody will read.

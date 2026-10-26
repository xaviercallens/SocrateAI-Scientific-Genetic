# Contribution note for SocrateAI-Scientific-Genetic: Price equation and multi-type haploid selection increment

Target repo: github.com/xaviercallens/SocrateAI-Scientific-Genetic. The parent session opens the PR; this note only
prepares it. Source: AutoevolveAI programme P7 WORK2, 2026-10-25. Details are in `PRICE.md` next to this note.

## What would be added

The proposed file is `formal/ANSE/ProgPopGenPrice.lean`, a copy of
`AutoevolveAI/formal_cert/Programmes/PopGen/Price.lean` (sha256
`1309c57a7818da57c9e955764df28fb26fb83c930661bd7c8439994d4d07438f`). It imports only `Mathlib.Algebra.BigOperators.Fin`
and `Mathlib.Tactic`, and contains:

- **Price equation** `price_equation`. For a finite population over any field, `w̄ (z̄' − z̄) = Cov_n(w,z) +
  E_n(w Δz)`, assuming `∑ n ≠ 0` and `∑ n w ≠ 0`. `price_equation_pos` gives the same identity under `n, w > 0` and a
  nonempty type set.
- **Haploid selection increment** `fisher_haploid_div` / `fisher_haploid_pos`: `Δw̄ = Var_n(w)/w̄` and `w̄ ≤ w̄'`, for
  any number of types. It assumes constant type fitnesses and no transmission bias.
- Controls:
  - one positive control;
  - two negative controls: transmission bias breaks the variance formula while Price still holds, and the `1/w̄`
    factor is needed;
  - non-vacuity.
  - The planted `sorry`/`axiom` controls stay in a separate file, `PriceControls.lean`, which is never imported.

## Overlap with `formal/ANSE/AlixPopulationGenetics.lean`

- The Alix file has H12 for **two alleles** over ℚ:
  - `aurelia_mean_fitness_increment`: `w̄' − w̄ = p(1−p)(w₁−w₂)²/w̄`;
  - `aurelia_mean_fitness_nondecreasing`;
  - `kal_increment_is_variance`, which says the numerator is the fitness variance.
- The new file **generalises the increment and the monotonicity to any finite number of types**, deriving the
  increment from Price's equation. Monotonicity is generalised over ℝ (`fisher_haploid_pos`); it was not re-derived
  over ℚ in the Alix form, so `aurelia_mean_fitness_nondecreasing` is not in the crosswalk.
- Namespace `ProgPopGen.Price.Alix` holds verbatim copies of `meanFit` and `selectP`, and restates two of the three
  Alix H12 theorems with the same statement text (`s = Fin 2`, `n = ![p, 1−p]`, `w = ![w₁, w₂]`):
  - `aurelia_mean_fitness_increment` is derived from the general `fisher_haploid_div`;
  - `kal_increment_is_variance` is the general `wvar` definition evaluated at `Fin 2` (`wvar_two`). It is not a
    consequence of Price.

  The match is textual; no Lean term links the two namespaces.
- Suggested follow-up inside the Genetic repo, not done here: state `ANSE.Alix.aurelia_mean_fitness_increment` as a
  corollary of the general theorem, or add an `example` linking the two. This would make the two-allele result a
  specialisation, not a separate proof.
- No overlap with H11 (Wright–Fisher harmonic uniqueness) in that file.
- Other P7 work, also suitable for a PR: the WORK1 fixation-limit modules close the gap H11 notes as
  "NOT formalised here", namely that `i/N` is the fixation probability. They are `Absorption.lean`,
  `WrightFisherFixation.lean`, `MoranModel.lean` and `MoranFixation.lean`; see `FIXATION.md`.
- A whole-repo grep of the Genetic clone (`*.lean`, `*.md`, `*.tex`, `*.py`) found no Price-equation or covariance
  decomposition result. The only matches for "price" are unrelated English uses. "covariance" appears only as a phylogenetic
  tree covariance matrix (`scripts/run_research_directions_rd1_rd5.py:182`) and once in the textbook source
  (`papers/University_Textbook_Mathematical_Genetics_Biophysics_Lean4.tex:491`, the covariance operator of a
  reticulate coalescent tree, `M_ret = M_tree + U Vᵀ`). Neither is a Price-type decomposition.

## Build compatibility (stated honestly)

- Everything was compiled under the AutoevolveAI pin: openai/math Lean 4.34.1 + full Mathlib.
- The Genetic repo README says Lean 4.34.1, but this clone has no `lean-toolchain` or `lakefile` at its root, so
  **the file has not been compiled inside the Genetic repo's own build**.
- As a reference, the existing `AlixPopulationGenetics.lean` compiles under the same pin: rc 0, 11 declarations, all
  `[propext, Classical.choice, Quot.sound]`; see `compile_genetic_alix_reference.json`.
- Before merging, the PR should be compiled in whatever build that repo actually uses.

## Wording for the PR / README (tier-honest)

> Adds a kernel-checked (Lean 4.34.1 + Mathlib; axioms propext, Classical.choice, Quot.sound) formalization of the
> classical Price equation for a finite population (Price 1970). Adds the multi-type haploid selection increment
> Δw̄ = Var(w)/w̄ under constant type fitnesses and no transmission bias, which generalises the existing two-allele
> H12. This is not Fisher's fundamental theorem in the partial-change / additive-variance sense, and not the diploid
> case. Known results; no novelty claimed. Producer-side check; independent recheck pending.

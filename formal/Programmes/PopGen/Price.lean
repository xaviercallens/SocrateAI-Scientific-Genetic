import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Tactic

/-!
# The Price equation for a finite population, and the haploid selection increment (P7, 2026-10-25)

Formalization of a KNOWN classical identity (G. R. Price, Nature 227 (1970) 520–21; modern statement in
T. LaGatta, arXiv:2202.10289v4, eq. (1.6)). Nothing here is new mathematics; the contribution is a
kernel-checked statement with explicit hypotheses, over an arbitrary field `K` and an arbitrary finite index
set `s`.

## Model (the modelling assumption is in the definitions)

A parent population is split into types `i ∈ s` with abundances (weights) `n i`, fitness `w i`
(offspring per individual), parental trait value `z i` and **mean offspring trait** `z' i`.
Type `i` contributes `n i * w i` offspring whose mean trait is `z' i`. Hence

* `wmean s n f = (∑ n f) / (∑ n)` — weighted (population) mean,
* `wcov s n f g = (∑ n (f - f̄)(g - ḡ)) / (∑ n)` — weighted (population, not sample) covariance,
* `zbarNext s n w z' = (∑ n w z') / (∑ n w)` — mean trait of the offspring generation.

With these definitions the Price equation `w̄ Δz̄ = Cov(w,z) + E(w Δz)` is an exact algebraic identity
(`price_equation`); it says nothing about which `w`, `z'` arise from a biological mechanism.

## Fisher-type corollary (scope, stated honestly)

Taking the trait to be fitness itself, `z = w`, and **no transmission bias** for that trait, `z' = w`
(each type's offspring have the parental type fitness: haploid/asexual types, constant type fitnesses,
no mutation), the transmission term vanishes and one gets `Δw̄ = Var_n(w) / w̄` (`fisher_haploid_div`)
and `w̄ ≤ w̄'` (`fisher_haploid_nondecreasing`). This is the selection-only increment for haploid
(asexual, clonal) types. It is NOT Fisher's fundamental theorem in the Price (1972) / Ewens
"partial change" reading, it involves the total, not the additive genetic, variance, and it is NOT the
diploid multi-allele case (where mean-fitness monotonicity is an inequality, not this equality).
-/

namespace ProgPopGen.Price

open Finset

section Field

variable {ι K : Type*} [Field K]

/-- Weighted mean `(∑ n f) / (∑ n)` over `s`. -/
def wmean (s : Finset ι) (n f : ι → K) : K :=
  (∑ i ∈ s, n i * f i) / ∑ i ∈ s, n i

/-- Weighted (population) covariance `(∑ n (f - f̄)(g - ḡ)) / (∑ n)` over `s`. -/
def wcov (s : Finset ι) (n f g : ι → K) : K :=
  (∑ i ∈ s, n i * (f i - wmean s n f) * (g i - wmean s n g)) / ∑ i ∈ s, n i

/-- Weighted (population) variance. -/
def wvar (s : Finset ι) (n f : ι → K) : K :=
  wcov s n f f

/-- Mean trait of the offspring generation: type `i` leaves `n i * w i` offspring of mean trait
`z' i`. -/
def zbarNext (s : Finset ι) (n w z' : ι → K) : K :=
  wmean s (fun i => n i * w i) z'

/-- Mean fitness of the offspring generation when type fitnesses are transmitted unchanged
(no transmission bias for fitness): `zbarNext` with `z' = w`. -/
def wbarNext (s : Finset ι) (n w : ι → K) : K :=
  zbarNext s n w w

/-- Expansion of the weighted covariance numerator into raw sums. -/
theorem cov_num_expand (s : Finset ι) (n f g : ι → K) (a b : K) :
    ∑ i ∈ s, n i * (f i - a) * (g i - b)
      = ∑ i ∈ s, n i * f i * g i - a * ∑ i ∈ s, n i * g i - b * ∑ i ∈ s, n i * f i
          + a * b * ∑ i ∈ s, n i := by
  have h : ∀ i ∈ s, n i * (f i - a) * (g i - b)
      = n i * f i * g i - a * (n i * g i) - b * (n i * f i) + a * b * n i := by
    intro i _; ring
  rw [Finset.sum_congr rfl h]
  simp only [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum]

/-- Expansion of the transmission-term numerator into raw sums. -/
theorem trans_num_expand (s : Finset ι) (n w z z' : ι → K) :
    ∑ i ∈ s, n i * (w i * (z' i - z i))
      = ∑ i ∈ s, n i * w i * z' i - ∑ i ∈ s, n i * w i * z i := by
  rw [← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  ring

/-- Scalar core of the Price equation (five sums treated as atoms). -/
theorem price_scalar (N A B C D : K) (hN : N ≠ 0) (hA : A ≠ 0) :
    A / N * (D / A - B / N)
      = (C - A / N * B - B / N * A + A / N * (B / N) * N) / N + (D - C) / N := by
  field_simp
  ring

/-- **Price equation** (finite population, exact identity):
`w̄ (z̄' - z̄) = Cov_n(w, z) + E_n(w Δz)`, with `Δz i = z' i - z i`.
Hypotheses: total abundance `∑ n ≠ 0` and total offspring number `∑ n w ≠ 0`. -/
theorem price_equation (s : Finset ι) (n w z z' : ι → K)
    (hN : ∑ i ∈ s, n i ≠ 0) (hW : ∑ i ∈ s, n i * w i ≠ 0) :
    wmean s n w * (zbarNext s n w z' - wmean s n z)
      = wcov s n w z + wmean s n (fun i => w i * (z' i - z i)) := by
  unfold wcov zbarNext
  rw [cov_num_expand]
  unfold wmean
  rw [trans_num_expand]
  have hz : ∀ i ∈ s, n i * w i * z i = n i * z i * w i := by intro i _; ring
  rw [Finset.sum_congr rfl hz]
  exact price_scalar _ _ _ _ _ hN hW

/-- The covariance in "raw moment" form: `Cov_n(f,g) = E_n(f g) - E_n(f) E_n(g)`. -/
theorem wcov_eq_moment (s : Finset ι) (n f g : ι → K) (hN : ∑ i ∈ s, n i ≠ 0) :
    wcov s n f g = wmean s n (fun i => f i * g i) - wmean s n f * wmean s n g := by
  unfold wcov
  rw [cov_num_expand]
  unfold wmean
  have h : ∀ i ∈ s, n i * (f i * g i) = n i * f i * g i := by intro i _; ring
  rw [Finset.sum_congr rfl h]
  field_simp
  ring

/-- Mean fitness is nonzero when both totals are. -/
theorem wmean_ne_zero (s : Finset ι) (n w : ι → K)
    (hN : ∑ i ∈ s, n i ≠ 0) (hW : ∑ i ∈ s, n i * w i ≠ 0) : wmean s n w ≠ 0 :=
  div_ne_zero hW hN

/-- **Haploid selection increment, from Price** (`z = w`, no transmission bias `z' = w`):
`w̄ (w̄' - w̄) = Var_n(w)`. The transmission term is identically `0`. -/
theorem fisher_haploid (s : Finset ι) (n w : ι → K)
    (hN : ∑ i ∈ s, n i ≠ 0) (hW : ∑ i ∈ s, n i * w i ≠ 0) :
    wmean s n w * (wbarNext s n w - wmean s n w) = wvar s n w := by
  have h := price_equation s n w w w hN hW
  simp only [sub_self, mul_zero] at h
  unfold wbarNext wvar
  rw [h]
  simp [wmean]

/-- `Δw̄ = Var_n(w) / w̄` (haploid/asexual types, constant type fitnesses, no transmission bias). -/
theorem fisher_haploid_div (s : Finset ι) (n w : ι → K)
    (hN : ∑ i ∈ s, n i ≠ 0) (hW : ∑ i ∈ s, n i * w i ≠ 0) :
    wbarNext s n w - wmean s n w = wvar s n w / wmean s n w := by
  rw [eq_div_iff (wmean_ne_zero s n w hN hW), mul_comm]
  exact fisher_haploid s n w hN hW

end Field

/-! ## Order consequences over `ℝ` -/

section Real

variable {ι : Type*}

theorem wvar_nonneg (s : Finset ι) (n f : ι → ℝ) (hn : ∀ i ∈ s, 0 ≤ n i) :
    0 ≤ wvar s n f := by
  unfold wvar wcov
  apply div_nonneg
  · apply Finset.sum_nonneg
    intro i hi
    rw [mul_assoc]
    exact mul_nonneg (hn i hi) (mul_self_nonneg _)
  · exact Finset.sum_nonneg hn

/-- Mean fitness does not decrease under haploid selection with constant type fitnesses
(nonnegative abundances, positive total abundance, positive mean fitness). -/
theorem fisher_haploid_nondecreasing (s : Finset ι) (n w : ι → ℝ)
    (hn : ∀ i ∈ s, 0 ≤ n i) (hN : 0 < ∑ i ∈ s, n i) (hW : 0 < ∑ i ∈ s, n i * w i) :
    wmean s n w ≤ wbarNext s n w := by
  have h := fisher_haploid_div s n w hN.ne' hW.ne'
  have hm : 0 < wmean s n w := div_pos hW hN
  have : 0 ≤ wvar s n w / wmean s n w := div_nonneg (wvar_nonneg s n w hn) hm.le
  linarith

/-- Price equation under the biological positivity hypotheses (`s` nonempty, `n i > 0`, `w i > 0`). -/
theorem price_equation_pos (s : Finset ι) (n w z z' : ι → ℝ) (hs : s.Nonempty)
    (hn : ∀ i ∈ s, 0 < n i) (hw : ∀ i ∈ s, 0 < w i) :
    wmean s n w * (zbarNext s n w z' - wmean s n z)
      = wcov s n w z + wmean s n (fun i => w i * (z' i - z i)) :=
  price_equation s n w z z' (Finset.sum_pos hn hs).ne'
    (Finset.sum_pos (fun i hi => mul_pos (hn i hi) (hw i hi)) hs).ne'

/-- Haploid selection increment under the biological positivity hypotheses. -/
theorem fisher_haploid_pos (s : Finset ι) (n w : ι → ℝ) (hs : s.Nonempty)
    (hn : ∀ i ∈ s, 0 < n i) (hw : ∀ i ∈ s, 0 < w i) :
    wbarNext s n w - wmean s n w = wvar s n w / wmean s n w ∧ wmean s n w ≤ wbarNext s n w := by
  have hN := Finset.sum_pos hn hs
  have hW := Finset.sum_pos (fun i hi => mul_pos (hn i hi) (hw i hi)) hs
  exact ⟨fisher_haploid_div s n w hN.ne' hW.ne',
    fisher_haploid_nondecreasing s n w (fun i hi => (hn i hi).le) hN hW⟩

end Real

/-! ## Crosswalk: the two-allele haploid statements of `ANSE.Alix` (SocrateAI-Scientific-Genetic,
`formal/ANSE/AlixPopulationGenetics.lean`), re-derived from `fisher_haploid_div`.

`meanFit` and `selectP` below are copied verbatim from that file (namespace changed); the two theorems
restate `aurelia_mean_fitness_increment` and `kal_increment_is_variance` with the same statements.
The proofs go through the general multi-type result with `s = univ : Finset (Fin 2)`,
`n = ![p, 1 - p]`, `w = ![w₁, w₂]`. The Genetic repo file itself is not imported. -/

namespace Alix

/-- Mean fitness with allele frequency `p` and allele fitnesses `w₁`, `w₂` (copied). -/
def meanFit (w₁ w₂ p : ℚ) : ℚ := p * w₁ + (1 - p) * w₂

/-- Allele frequency after one round of haploid selection (copied). -/
def selectP (w₁ w₂ p : ℚ) : ℚ := p * w₁ / meanFit w₁ w₂ p

theorem wmean_two (w₁ w₂ p : ℚ) :
    wmean (Finset.univ : Finset (Fin 2)) ![p, 1 - p] ![w₁, w₂] = meanFit w₁ w₂ p := by
  simp [wmean, meanFit, Fin.sum_univ_two]

theorem wbarNext_two (w₁ w₂ p : ℚ) (hm : meanFit w₁ w₂ p ≠ 0) :
    wbarNext (Finset.univ : Finset (Fin 2)) ![p, 1 - p] ![w₁, w₂]
      = meanFit w₁ w₂ (selectP w₁ w₂ p) := by
  unfold meanFit at hm
  simp only [wbarNext, zbarNext, wmean, selectP, meanFit, Fin.sum_univ_two, Matrix.cons_val_zero,
    Matrix.cons_val_one, Matrix.cons_val_fin_one]
  field_simp
  ring

theorem wvar_two (w₁ w₂ p : ℚ) :
    wvar (Finset.univ : Finset (Fin 2)) ![p, 1 - p] ![w₁, w₂] = p * (1 - p) * (w₁ - w₂) ^ 2 := by
  simp only [wvar, wcov, wmean, Fin.sum_univ_two, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_fin_one]
  ring

/-- Same statement as `ANSE.Alix.aurelia_mean_fitness_increment`, derived from `fisher_haploid_div`. -/
theorem aurelia_mean_fitness_increment (w₁ w₂ p : ℚ) (hm : meanFit w₁ w₂ p ≠ 0) :
    meanFit w₁ w₂ (selectP w₁ w₂ p) - meanFit w₁ w₂ p
      = p * (1 - p) * (w₁ - w₂) ^ 2 / meanFit w₁ w₂ p := by
  have hN : ∑ i : Fin 2, (![p, 1 - p] : Fin 2 → ℚ) i ≠ 0 := by simp [Fin.sum_univ_two]
  have hW : ∑ i : Fin 2, (![p, 1 - p] : Fin 2 → ℚ) i * (![w₁, w₂] : Fin 2 → ℚ) i ≠ 0 := by
    simpa [Fin.sum_univ_two, meanFit] using hm
  have h := fisher_haploid_div (Finset.univ : Finset (Fin 2)) ![p, 1 - p] ![w₁, w₂] hN hW
  rwa [wbarNext_two w₁ w₂ p hm, wmean_two, wvar_two] at h

/-- Same statement as `ANSE.Alix.kal_increment_is_variance`, derived from `wvar_two` and the
general variance definition. -/
theorem kal_increment_is_variance (w₁ w₂ p : ℚ) :
    p * (1 - p) * (w₁ - w₂) ^ 2
      = p * (w₁ - meanFit w₁ w₂ p) ^ 2 + (1 - p) * (w₂ - meanFit w₁ w₂ p) ^ 2 := by
  rw [← wvar_two]
  simp only [wvar, wcov, Fin.sum_univ_two, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_fin_one, wmean_two]
  ring

end Alix

/-! ## Controls (claimed, kernel-checked; each side evaluated by `norm_num`, not by the theorem) -/

/-- Positive control: three types, `n = (1,2,3)`, `w = (1,2,4)`, `z = (0,1,3)`, `z' = (1,1,2)`.
Both sides of the Price equation evaluate to `-13/36` (computed directly). -/
theorem ctl_pos_three_types :
    wmean (Finset.univ : Finset (Fin 3)) ![1, 2, 3] ![1, 2, 4]
        * (zbarNext Finset.univ ![1, 2, 3] ![1, 2, 4] ![1, 1, 2]
          - wmean Finset.univ ![1, 2, 3] ![0, 1, 3]) = (-13 / 36 : ℚ) ∧
    wcov (Finset.univ : Finset (Fin 3)) ![1, 2, 3] ![1, 2, 4] ![0, 1, 3]
        + wmean Finset.univ ![1, 2, 3]
          (fun i => (![1, 2, 4] : Fin 3 → ℚ) i * ((![1, 1, 2] : Fin 3 → ℚ) i - (![0, 1, 3] : Fin 3 → ℚ) i))
        = (-13 / 36 : ℚ) := by
  constructor <;>
  simp [wmean, wcov, zbarNext, Fin.sum_univ_three] <;> norm_num

/-- Negative control (transmission bias): `n = (1,1)`, `w = z = (1,2)` but offspring fitness
`z' = (1,1)`. Then `Δw̄ = -1/2` while `Var(w)/w̄ = 1/6`: the variance formula FAILS. The Price
equation still holds on the same instance: both sides equal `-3/4`. -/
theorem ctl_neg_transmission_bias :
    zbarNext (Finset.univ : Finset (Fin 2)) ![1, 1] ![1, 2] ![1, 1]
        - wmean Finset.univ ![1, 1] ![1, 2] = (-1 / 2 : ℚ) ∧
    wvar (Finset.univ : Finset (Fin 2)) ![1, 1] ![1, 2] / wmean Finset.univ ![1, 1] ![1, 2]
        = (1 / 6 : ℚ) ∧
    wmean (Finset.univ : Finset (Fin 2)) ![1, 1] ![1, 2]
        * (zbarNext Finset.univ ![1, 1] ![1, 2] ![1, 1] - wmean Finset.univ ![1, 1] ![1, 2])
        = (-3 / 4 : ℚ) ∧
    wcov (Finset.univ : Finset (Fin 2)) ![1, 1] ![1, 2] ![1, 2]
        + wmean Finset.univ ![1, 1]
          (fun i => (![1, 2] : Fin 2 → ℚ) i * ((![1, 1] : Fin 2 → ℚ) i - (![1, 2] : Fin 2 → ℚ) i))
        = (-3 / 4 : ℚ) := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;>
  simp [wmean, wcov, wvar, zbarNext, Fin.sum_univ_two] <;> norm_num

/-- Negative control (the `1/w̄` factor matters): `n = (1,1)`, `w = (1,2)`, no transmission bias.
`Δw̄ = 1/6 = Var(w)/w̄`, but `Var(w) = 1/4 ≠ Δw̄`. -/
theorem ctl_neg_missing_wbar :
    wbarNext (Finset.univ : Finset (Fin 2)) ![1, 1] ![1, 2] - wmean Finset.univ ![1, 1] ![1, 2]
        = (1 / 6 : ℚ) ∧
    wvar (Finset.univ : Finset (Fin 2)) ![1, 1] (![1, 2] : Fin 2 → ℚ) = 1 / 4 ∧
    wvar (Finset.univ : Finset (Fin 2)) ![1, 1] (![1, 2] : Fin 2 → ℚ)
      ≠ wbarNext Finset.univ ![1, 1] ![1, 2] - wmean Finset.univ ![1, 1] ![1, 2] := by
  refine ⟨?_, ?_, ?_⟩ <;>
  simp [wmean, wcov, wvar, wbarNext, zbarNext, Fin.sum_univ_two] <;> norm_num

/-- Non-vacuity: the positivity hypotheses of `price_equation_pos` hold for a concrete population. -/
theorem nonvacuous_pos :
    (Finset.univ : Finset (Fin 3)).Nonempty ∧ (∀ i ∈ (Finset.univ : Finset (Fin 3)),
      0 < (![1, 2, 3] : Fin 3 → ℝ) i) ∧ (∀ i ∈ (Finset.univ : Finset (Fin 3)),
      0 < (![1, 2, 4] : Fin 3 → ℝ) i) := by
  refine ⟨Finset.univ_nonempty, ?_, ?_⟩ <;> intro i _ <;> fin_cases i <;> norm_num

end ProgPopGen.Price

#print axioms ProgPopGen.Price.cov_num_expand
#print axioms ProgPopGen.Price.trans_num_expand
#print axioms ProgPopGen.Price.price_scalar
#print axioms ProgPopGen.Price.price_equation
#print axioms ProgPopGen.Price.wcov_eq_moment
#print axioms ProgPopGen.Price.wmean_ne_zero
#print axioms ProgPopGen.Price.fisher_haploid
#print axioms ProgPopGen.Price.fisher_haploid_div
#print axioms ProgPopGen.Price.wvar_nonneg
#print axioms ProgPopGen.Price.fisher_haploid_nondecreasing
#print axioms ProgPopGen.Price.price_equation_pos
#print axioms ProgPopGen.Price.fisher_haploid_pos
#print axioms ProgPopGen.Price.Alix.wmean_two
#print axioms ProgPopGen.Price.Alix.wbarNext_two
#print axioms ProgPopGen.Price.Alix.wvar_two
#print axioms ProgPopGen.Price.Alix.aurelia_mean_fitness_increment
#print axioms ProgPopGen.Price.Alix.kal_increment_is_variance
#print axioms ProgPopGen.Price.ctl_pos_three_types
#print axioms ProgPopGen.Price.ctl_neg_transmission_bias
#print axioms ProgPopGen.Price.ctl_neg_missing_wbar
#print axioms ProgPopGen.Price.nonvacuous_pos

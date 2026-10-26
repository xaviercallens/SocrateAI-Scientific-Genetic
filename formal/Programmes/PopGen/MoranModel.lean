import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Algebra.Field.GeomSum
import Mathlib.Tactic

/-!
# Moran model with selection: rates and the fixation function (P7, programmes 2026-10-25)

Population of `N` individuals, `i` mutants of relative fitness `r > 0`. One step of the Moran process:
an individual is chosen to reproduce with probability proportional to fitness, and an individual
chosen uniformly at random dies. The number of mutants is a birth–death chain on `{0,…,N}` with

* `p i = (r i / (r i + N - i)) · ((N - i) / N)` (up),
* `q i = ((N - i) / (r i + N - i)) · (i / N)` (down),
* stay with probability `1 - p i - q i`.

Fixation function: `moranFix N r i = S_i / S_N` with `S_n = ∑_{k<n} r⁻ᵏ`. This file proves, at the level
of ℕ-indexed first-step equations: the rates are probabilities, `p i = r · q i`, the first-step
(harmonic) equation for `moranFix` at every interior state, positivity of the one-step variance, the
closed form `(1 - r⁻ⁱ)/(1 - r⁻ᴺ)` for `r ≠ 1`, the neutral identity `moranFix N 1 i = i/N`, and the
limit `r → 1`. The chain-level statements (absorption limit, uniqueness) are in `MoranFixation.lean`.

Known classical result (Moran 1958 and textbooks); novelty is not claimed.
-/

namespace ProgPopGen

open Finset Filter Topology

/-- Birth–death kernel on `{0,…,N}`: up with probability `p i`, down with `q i`, stay otherwise. -/
noncomputable def bdP (N : ℕ) (p q : ℕ → ℝ) (i j : Fin (N + 1)) : ℝ :=
  (if (j : ℕ) = (i : ℕ) + 1 then p i else 0) + (if (j : ℕ) + 1 = (i : ℕ) then q i else 0) +
    (if (j : ℕ) = (i : ℕ) then 1 - p i - q i else 0)

theorem sum_fin_ite_val {N : ℕ} (c : ℕ) (a : ℝ) (f : ℕ → ℝ) :
    ∑ j : Fin (N + 1), (if (j : ℕ) = c then a else 0) * f j
      = if c < N + 1 then a * f c else 0 := by
  rw [Fin.sum_univ_eq_sum_range (fun m => (if m = c then a else 0) * f m) (N + 1)]
  simp [ite_mul, Finset.mem_range]

/-- First-step expectation for a birth–death kernel (`p N = 0`, `q 0 = 0`). -/
theorem bdP_sum {N : ℕ} (p q : ℕ → ℝ) (hpN : p N = 0) (hq0 : q 0 = 0) (i : Fin (N + 1))
    (f : ℕ → ℝ) :
    ∑ j, bdP N p q i j * f j
      = p i * f (i + 1) + q i * f (i - 1) + (1 - p i - q i) * f i := by
  have hi := i.isLt
  have hdown : ∑ j : Fin (N + 1), (if (j : ℕ) + 1 = (i : ℕ) then q i else 0) * f j
      = q i * f (i - 1) := by
    rcases Nat.eq_zero_or_pos (i : ℕ) with h0 | hpos
    · rw [h0, hq0]; simp
    · obtain ⟨m, hm⟩ : ∃ m, (i : ℕ) = m + 1 := ⟨(i : ℕ) - 1, by omega⟩
      have e : ∀ j : Fin (N + 1), (if (j : ℕ) + 1 = (i : ℕ) then q i else 0)
          = (if (j : ℕ) = m then q i else 0) := by
        intro j
        rw [hm]
        simp only [Nat.add_right_cancel_iff]
      simp only [e]
      rw [sum_fin_ite_val, if_pos (by omega), hm]
      simp
  have e : ∀ j : Fin (N + 1), bdP N p q i j * f j =
      (if (j : ℕ) = (i : ℕ) + 1 then p i else 0) * f j
        + (if (j : ℕ) + 1 = (i : ℕ) then q i else 0) * f j
        + (if (j : ℕ) = (i : ℕ) then 1 - p i - q i else 0) * f j := fun j => by
    unfold bdP; ring
  rw [Finset.sum_congr rfl (fun j _ => e j), Finset.sum_add_distrib, Finset.sum_add_distrib,
    hdown, sum_fin_ite_val, sum_fin_ite_val, if_pos hi]
  by_cases hiN : (i : ℕ) + 1 < N + 1
  · rw [if_pos hiN]
  · rw [if_neg hiN]
    have : (i : ℕ) = N := by omega
    rw [this, hpN]
    ring

/-- Up-rate of the Moran chain. -/
noncomputable def moranUp (N : ℕ) (r : ℝ) (i : ℕ) : ℝ :=
  r * i / (r * i + ((N : ℝ) - i)) * (((N : ℝ) - i) / N)

/-- Down-rate of the Moran chain. -/
noncomputable def moranDown (N : ℕ) (r : ℝ) (i : ℕ) : ℝ :=
  ((N : ℝ) - i) / (r * i + ((N : ℝ) - i)) * ((i : ℝ) / N)

/-- Moran transition kernel on `Fin (N+1)`. -/
noncomputable def moranP (N : ℕ) (r : ℝ) : Fin (N + 1) → Fin (N + 1) → ℝ :=
  bdP N (moranUp N r) (moranDown N r)

theorem moranUp_last (N : ℕ) (r : ℝ) : moranUp N r N = 0 := by simp [moranUp]
theorem moranUp_zero (N : ℕ) (r : ℝ) : moranUp N r 0 = 0 := by simp [moranUp]
theorem moranDown_zero (N : ℕ) (r : ℝ) : moranDown N r 0 = 0 := by simp [moranDown]
theorem moranDown_last (N : ℕ) (r : ℝ) : moranDown N r N = 0 := by simp [moranDown]

/-- `p i = r · q i`: the ratio of up to down rates is the fitness. -/
theorem moranUp_eq (N : ℕ) (r : ℝ) (i : ℕ) : moranUp N r i = r * moranDown N r i := by
  unfold moranUp moranDown
  ring

theorem moran_den_pos {N : ℕ} {r : ℝ} (hN : 0 < N) (hr : 0 < r) {i : ℕ} (hi : i ≤ N) :
    0 < r * i + ((N : ℝ) - i) := by
  have hiN : (i : ℝ) ≤ N := by exact_mod_cast hi
  rcases Nat.eq_zero_or_pos i with h | h
  · subst h
    simpa using (show (0 : ℝ) < N by exact_mod_cast hN)
  · have : (0 : ℝ) < r * i := mul_pos hr (by exact_mod_cast h)
    linarith

theorem moranDown_nonneg {N : ℕ} {r : ℝ} (hN : 0 < N) (hr : 0 < r) {i : ℕ} (hi : i ≤ N) :
    0 ≤ moranDown N r i := by
  have hiN : (0 : ℝ) ≤ (N : ℝ) - i := by
    have : (i : ℝ) ≤ N := by exact_mod_cast hi
    linarith
  have hD := moran_den_pos hN hr hi
  unfold moranDown
  exact mul_nonneg (div_nonneg hiN hD.le) (div_nonneg (Nat.cast_nonneg i) (Nat.cast_nonneg N))

theorem moranUp_nonneg {N : ℕ} {r : ℝ} (hN : 0 < N) (hr : 0 < r) {i : ℕ} (hi : i ≤ N) :
    0 ≤ moranUp N r i := by
  rw [moranUp_eq]
  exact mul_nonneg hr.le (moranDown_nonneg hN hr hi)

/-- Holding probability: `1 - p - q = (r i² + (N - i)²) / ((r i + N - i) N)`. -/
theorem moran_stay_eq {N : ℕ} {r : ℝ} (hN : 0 < N) (hr : 0 < r) {i : ℕ} (hi : i ≤ N) :
    1 - moranUp N r i - moranDown N r i
      = (r * (i : ℝ) ^ 2 + ((N : ℝ) - i) ^ 2) / ((r * i + ((N : ℝ) - i)) * N) := by
  have hD := (moran_den_pos hN hr hi).ne'
  have hNr : (N : ℝ) ≠ 0 := by exact_mod_cast hN.ne'
  unfold moranUp moranDown
  field_simp
  ring

theorem moran_stay_nonneg {N : ℕ} {r : ℝ} (hN : 0 < N) (hr : 0 < r) {i : ℕ} (hi : i ≤ N) :
    0 ≤ 1 - moranUp N r i - moranDown N r i := by
  rw [moran_stay_eq hN hr hi]
  apply div_nonneg
  · exact add_nonneg (mul_nonneg hr.le (sq_nonneg _)) (sq_nonneg _)
  · exact mul_nonneg (moran_den_pos hN hr hi).le (Nat.cast_nonneg N)

theorem moranDown_pos {N : ℕ} {r : ℝ} (hN : 0 < N) (hr : 0 < r) {i : ℕ} (h1 : 0 < i)
    (h2 : i < N) : 0 < moranDown N r i := by
  have hiN : (0 : ℝ) < (N : ℝ) - i := by
    have : (i : ℝ) < N := by exact_mod_cast h2
    linarith
  have hD := moran_den_pos hN hr h2.le
  unfold moranDown
  exact mul_pos (div_pos hiN hD) (div_pos (by exact_mod_cast h1) (by exact_mod_cast hN))

theorem moranUp_pos {N : ℕ} {r : ℝ} (hN : 0 < N) (hr : 0 < r) {i : ℕ} (h1 : 0 < i)
    (h2 : i < N) : 0 < moranUp N r i := by
  rw [moranUp_eq]
  exact mul_pos hr (moranDown_pos hN hr h1 h2)

/-! ## The fixation function -/

/-- `S_n = ∑_{k<n} r⁻ᵏ`. -/
noncomputable def geomS (r : ℝ) (n : ℕ) : ℝ := ∑ k ∈ range n, r⁻¹ ^ k

/-- Moran fixation function `S_i / S_N`. -/
noncomputable def moranFix (N : ℕ) (r : ℝ) (i : ℕ) : ℝ := geomS r i / geomS r N

theorem geomS_succ (r : ℝ) (n : ℕ) : geomS r (n + 1) = geomS r n + r⁻¹ ^ n := by
  unfold geomS
  rw [Finset.sum_range_succ]

theorem geomS_pos {r : ℝ} (hr : 0 < r) {n : ℕ} (hn : 0 < n) : 0 < geomS r n := by
  unfold geomS
  exact Finset.sum_pos (fun k _ => pow_pos (inv_pos.mpr hr) k)
    (Finset.nonempty_range_iff.mpr hn.ne')

theorem moranFix_zero (N : ℕ) (r : ℝ) : moranFix N r 0 = 0 := by
  simp [moranFix, geomS]

theorem moranFix_last {N : ℕ} {r : ℝ} (hN : 0 < N) (hr : 0 < r) : moranFix N r N = 1 :=
  div_self (geomS_pos hr hN).ne'

/-- First-step (harmonic) equation of the Moran chain for `moranFix`, at interior states. -/
theorem moranFix_harmonic_nat {N : ℕ} {r : ℝ} (hr : 0 < r) {i : ℕ} (h1 : 0 < i) :
    moranUp N r i * moranFix N r (i + 1) + moranDown N r i * moranFix N r (i - 1)
      + (1 - moranUp N r i - moranDown N r i) * moranFix N r i = moranFix N r i := by
  obtain ⟨m, rfl⟩ : ∃ m, i = m + 1 := ⟨i - 1, by omega⟩
  have hpq := moranUp_eq N r (m + 1)
  have hrr : r * r⁻¹ = 1 := mul_inv_cancel₀ hr.ne'
  simp only [moranFix, Nat.add_sub_cancel, geomS_succ]
  linear_combination (r⁻¹ ^ m * r⁻¹ / geomS r N) * hpq
    + (moranDown N r (m + 1) * r⁻¹ ^ m / geomS r N) * hrr

theorem moranFix_step (N : ℕ) (r : ℝ) (i : ℕ) :
    moranFix N r (i + 1) - moranFix N r i = r⁻¹ ^ i / geomS r N := by
  unfold moranFix
  rw [geomS_succ, ← sub_div, add_sub_cancel_left]

/-- One-step variance of `moranFix` is positive at interior states. -/
theorem moran_var_pos_nat {N : ℕ} {r : ℝ} (hN : 0 < N) (hr : 0 < r) {i : ℕ} (h1 : 0 < i)
    (h2 : i < N) :
    0 < moranUp N r i * (moranFix N r (i + 1) - moranFix N r i) ^ 2
      + moranDown N r i * (moranFix N r (i - 1) - moranFix N r i) ^ 2
      + (1 - moranUp N r i - moranDown N r i) * (moranFix N r i - moranFix N r i) ^ 2 := by
  rw [sub_self, moranFix_step, zero_pow two_ne_zero, mul_zero, add_zero]
  have hS := geomS_pos hr hN
  have hx : 0 < (r⁻¹ ^ i / geomS r N) ^ 2 := pow_pos (div_pos (pow_pos (inv_pos.mpr hr) i) hS) 2
  have hu := mul_pos (moranUp_pos hN hr h1 h2) hx
  have hd := mul_nonneg (moranDown_nonneg hN hr h2.le)
    (sq_nonneg (moranFix N r (i - 1) - moranFix N r i))
  linarith

/-- Closed form for `r ≠ 1`: `moranFix N r i = (1 - r⁻ⁱ) / (1 - r⁻ᴺ)`. -/
theorem moranFix_closed {N : ℕ} {r : ℝ} (hr1 : r ≠ 1) (i : ℕ) :
    moranFix N r i = (1 - r⁻¹ ^ i) / (1 - r⁻¹ ^ N) := by
  have hρ : r⁻¹ ≠ 1 := fun h => hr1 (inv_eq_one.mp h)
  have hρ' : r⁻¹ - 1 ≠ 0 := sub_ne_zero.mpr hρ
  unfold moranFix geomS
  rw [geom_sum_eq hρ, geom_sum_eq hρ, div_div_div_cancel_right₀ hρ', ← neg_div_neg_eq, neg_sub,
    neg_sub]

theorem geomS_one (n : ℕ) : geomS 1 n = n := by
  simp [geomS]

/-- Neutral case `r = 1`: the fixation function is exactly `i/N`. -/
theorem moranFix_neutral (N i : ℕ) : moranFix N 1 i = (i : ℝ) / N := by
  simp [moranFix, geomS_one]

theorem geomS_tendsto_one (n : ℕ) : Tendsto (fun r : ℝ => geomS r n) (𝓝 1) (𝓝 (n : ℝ)) := by
  have hinv : Tendsto (fun r : ℝ => r⁻¹) (𝓝 1) (𝓝 1) := by
    simpa using (tendsto_inv₀ (one_ne_zero : (1 : ℝ) ≠ 0))
  have := tendsto_finsetSum (range n) (fun k _ => hinv.pow k)
  simpa [geomS] using this

/-- `moranFix N r i → i/N` as `r → 1`. -/
theorem moranFix_tendsto_one {N : ℕ} (hN : 0 < N) (i : ℕ) :
    Tendsto (fun r : ℝ => moranFix N r i) (𝓝 1) (𝓝 ((i : ℝ) / N)) := by
  have hNr : (N : ℝ) ≠ 0 := by exact_mod_cast hN.ne'
  exact (geomS_tendsto_one i).div (geomS_tendsto_one N) hNr

/-- Neutral limit of the closed form: `(1 - r⁻ⁱ)/(1 - r⁻ᴺ) → i/N` as `r → 1`, `r ≠ 1`. -/
theorem moran_closed_tendsto_neutral {N : ℕ} (hN : 0 < N) (i : ℕ) :
    Tendsto (fun r : ℝ => (1 - r⁻¹ ^ i) / (1 - r⁻¹ ^ N)) (𝓝[≠] 1) (𝓝 ((i : ℝ) / N)) := by
  refine ((moranFix_tendsto_one hN i).mono_left nhdsWithin_le_nhds).congr' ?_
  filter_upwards [self_mem_nhdsWithin] with r hr
  exact moranFix_closed hr i

/-! ## Controls (real theorems; planted sorry/axiom controls live in `PopGenProgControls.lean`) -/

/-- Positive control: rates for `N = 2`, `r = 2`, `i = 1`. -/
theorem ctl_pos_moran_rates : moranUp 2 2 1 = 1 / 3 ∧ moranDown 2 2 1 = 1 / 6 := by
  constructor <;> norm_num [moranUp, moranDown]

/-- Positive control: `N = 2`, `r = 2`, `i = 1` gives `2/3`; equals the closed form. -/
theorem ctl_pos_moranFix_N2 : moranFix 2 2 1 = 2 / 3 ∧ (1 - (2 : ℝ)⁻¹ ^ 1) / (1 - (2 : ℝ)⁻¹ ^ 2) = 2 / 3 := by
  constructor <;> norm_num [moranFix, geomS, Finset.sum_range_succ]

/-- Positive control: `N = 3`, `r = 2` gives `4/7` and `6/7`, and the first-step equation holds at
`i = 1` with these numbers. -/
theorem ctl_pos_moranFix_N3 :
    moranFix 3 2 1 = 4 / 7 ∧ moranFix 3 2 2 = 6 / 7 ∧
    moranUp 3 2 1 * (6 / 7) + moranDown 3 2 1 * 0 + (1 - moranUp 3 2 1 - moranDown 3 2 1) * (4 / 7)
      = (4 / 7 : ℝ) := by
  refine ⟨?_, ?_, ?_⟩ <;> norm_num [moranFix, geomS, Finset.sum_range_succ, moranUp, moranDown]

/-- Exponent-flipped (wrong) candidate `(1 - rⁱ)/(1 - rᴺ)`. -/
noncomputable def flipFix (N : ℕ) (r : ℝ) (i : ℕ) : ℝ := (1 - r ^ i) / (1 - r ^ N)

/-- Negative control: the flipped formula gives `1/7 ≠ 4/7` at `N = 3`, `r = 2`, `i = 1`. -/
theorem ctl_neg_flipped_value : flipFix 3 2 1 = 1 / 7 ∧ flipFix 3 2 1 ≠ moranFix 3 2 1 := by
  constructor
  · norm_num [flipFix]
  · norm_num [flipFix, moranFix, geomS, Finset.sum_range_succ]

/-- Negative control: the flipped formula violates the first-step equation at `N = 3`, `r = 2`,
`i = 1` (`3/14 ≠ 1/7`). -/
theorem ctl_neg_flipped_not_harmonic :
    moranUp 3 2 1 * flipFix 3 2 2 + moranDown 3 2 1 * flipFix 3 2 0
      + (1 - moranUp 3 2 1 - moranDown 3 2 1) * flipFix 3 2 1 ≠ flipFix 3 2 1 := by
  norm_num [moranUp, moranDown, flipFix]

/-- Negative control: the neutral formula `i/N` violates the first-step equation under selection
(`N = 3`, `r = 2`, `i = 1`: `7/18 ≠ 1/3`). -/
theorem ctl_neg_neutral_not_harmonic :
    moranUp 3 2 1 * (2 / 3) + moranDown 3 2 1 * 0
      + (1 - moranUp 3 2 1 - moranDown 3 2 1) * (1 / 3) ≠ (1 / 3 : ℝ) := by
  norm_num [moranUp, moranDown]

end ProgPopGen

#print axioms ProgPopGen.bdP_sum
#print axioms ProgPopGen.moranUp_eq
#print axioms ProgPopGen.moran_stay_eq
#print axioms ProgPopGen.moran_stay_nonneg
#print axioms ProgPopGen.moranFix_harmonic_nat
#print axioms ProgPopGen.moran_var_pos_nat
#print axioms ProgPopGen.moranFix_closed
#print axioms ProgPopGen.moranFix_neutral
#print axioms ProgPopGen.moranFix_tendsto_one
#print axioms ProgPopGen.moran_closed_tendsto_neutral
#print axioms ProgPopGen.ctl_pos_moran_rates
#print axioms ProgPopGen.ctl_pos_moranFix_N2
#print axioms ProgPopGen.ctl_pos_moranFix_N3
#print axioms ProgPopGen.ctl_neg_flipped_value
#print axioms ProgPopGen.ctl_neg_flipped_not_harmonic
#print axioms ProgPopGen.ctl_neg_neutral_not_harmonic

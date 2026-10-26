import Programmes.PopGen.Absorption
import Programmes.PopGen.MoranModel

/-!
# Moran model with selection: fixation probability (P7, programmes 2026-10-25)

Chain-level statements for the Moran birth–death chain `moranP N r` (fitness `r > 0`) on `{0,…,N}`:

* `moran_absorbing`: stochastic, `0` and `N` absorbing.
* `moran_interiorHarmonic`: `moranFix N r` satisfies the first-step equation at interior states.
* `moran_fixation_unique`: it is the unique such function with boundary values `0`, `1`.
* `moran_fixation_limit`: `Pⁿ(i,N) → moranFix N r i` and `Pⁿ(i,0) → 1 - moranFix N r i`.
* `moran_fixation_limit_closed` (`r ≠ 1`): the limit is `(1 - r⁻ⁱ)/(1 - r⁻ᴺ)`.
* `moran_fixation_limit_neutral` (`r = 1`): the limit is `i/N`.

Uniqueness is obtained from the absorption limit (two harmonic functions with the same boundary values
are both limits of the same sequence), not from a separate maximum principle. Known classical result;
novelty not claimed.
-/

namespace ProgPopGen

open Finset Filter Topology

theorem moran_absorbing {N : ℕ} {r : ℝ} (hN : 0 < N) (hr : 0 < r) :
    AbsorbingChain (moranP N r) where
  nonneg i j := by
    have hi : (i : ℕ) ≤ N := Nat.lt_succ_iff.mp i.isLt
    have h1 := moranUp_nonneg hN hr hi
    have h2 := moranDown_nonneg hN hr hi
    have h3 := moran_stay_nonneg hN hr hi
    unfold moranP bdP
    refine add_nonneg (add_nonneg ?_ ?_) ?_ <;> split_ifs <;> first | assumption | norm_num
  rowSum i := by
    have := bdP_sum (moranUp N r) (moranDown N r) (moranUp_last N r) (moranDown_zero N r) i
      (fun _ => 1)
    simp only [mul_one] at this
    unfold moranP
    rw [this]
    ring
  absorb0 j := by
    unfold moranP bdP
    simp only [Fin.val_zero, moranUp_zero, moranDown_zero, sub_zero, ite_self, zero_add]
    by_cases hj : j = 0
    · subst hj; simp
    · have hj' : (j : ℕ) ≠ 0 := fun h => hj (Fin.ext h)
      simp [hj, hj']
  absorbN j := by
    unfold moranP bdP
    simp only [Fin.val_last, moranUp_last, moranDown_last, sub_zero, ite_self, zero_add]
    by_cases hj : j = Fin.last N
    · subst hj; simp
    · have hj' : (j : ℕ) ≠ N := fun h => hj (Fin.ext h)
      simp [hj, hj']

theorem moran_interiorHarmonic {N : ℕ} {r : ℝ} (hr : 0 < r) :
    InteriorHarmonic (moranP N r) (fun j => moranFix N r j) := by
  intro i h1 _
  show ∑ j, bdP N (moranUp N r) (moranDown N r) i j * moranFix N r j = moranFix N r i
  rw [bdP_sum _ _ (moranUp_last N r) (moranDown_zero N r)]
  exact moranFix_harmonic_nat hr h1

theorem moran_var_pos {N : ℕ} {r : ℝ} (hN : 0 < N) (hr : 0 < r) (k : Fin (N + 1))
    (h1 : 0 < (k : ℕ)) (h2 : (k : ℕ) < N) :
    0 < ∑ j, moranP N r k j * (moranFix N r j - moranFix N r k) ^ 2 := by
  have e := bdP_sum (moranUp N r) (moranDown N r) (moranUp_last N r) (moranDown_zero N r) k
    (fun m : ℕ => (moranFix N r m - moranFix N r k) ^ 2)
  refine lt_of_lt_of_eq (moran_var_pos_nat hN hr h1 h2) ?_
  exact e.symm

/-- Absorption-limit form of the Moran fixation probability, for every fitness `r > 0`. -/
theorem moran_fixation_limit {N : ℕ} {r : ℝ} (hN : 0 < N) (hr : 0 < r) (i : Fin (N + 1)) :
    Tendsto (fun n => kPow (moranP N r) n i (Fin.last N)) atTop (𝓝 (moranFix N r i)) ∧
    Tendsto (fun n => kPow (moranP N r) n i 0) atTop (𝓝 (1 - moranFix N r i)) :=
  fixation_limit (moran_absorbing hN hr) hN (moran_interiorHarmonic hr)
    (by simp [moranFix_zero]) (by simpa using moranFix_last hN hr) (moran_var_pos hN hr) i

/-- Selection, `r ≠ 1`: `Pⁿ(i,N) → (1 - r⁻ⁱ)/(1 - r⁻ᴺ)`. -/
theorem moran_fixation_limit_closed {N : ℕ} {r : ℝ} (hN : 0 < N) (hr : 0 < r) (hr1 : r ≠ 1)
    (i : Fin (N + 1)) :
    Tendsto (fun n => kPow (moranP N r) n i (Fin.last N)) atTop
      (𝓝 ((1 - r⁻¹ ^ (i : ℕ)) / (1 - r⁻¹ ^ N))) := by
  rw [← moranFix_closed hr1]
  exact (moran_fixation_limit hN hr i).1

/-- Neutral Moran chain, `r = 1`: `Pⁿ(i,N) → i/N`. -/
theorem moran_fixation_limit_neutral {N : ℕ} (hN : 0 < N) (i : Fin (N + 1)) :
    Tendsto (fun n => kPow (moranP N 1) n i (Fin.last N)) atTop (𝓝 (((i : ℕ) : ℝ) / N)) := by
  rw [← moranFix_neutral]
  exact (moran_fixation_limit hN one_pos i).1

/-- Uniqueness: any function harmonic at the interior states of the Moran chain with boundary values
`0` and `1` is `moranFix N r`. -/
theorem moran_fixation_unique {N : ℕ} {r : ℝ} (hN : 0 < N) (hr : 0 < r) {g : Fin (N + 1) → ℝ}
    (hg : InteriorHarmonic (moranP N r) g) (g0 : g 0 = 0) (g1 : g (Fin.last N) = 1) :
    ∀ i, g i = moranFix N r i :=
  harmonic_unique_of_var (moran_absorbing hN hr) hN (moran_interiorHarmonic hr)
    (by simp [moranFix_zero]) (by simpa using moranFix_last hN hr) (moran_var_pos hN hr) hg g0 g1

/-- Negative control at chain level: `Pⁿ(1,3)` for `N = 3`, `r = 2` does NOT tend to the neutral
value `1/3`, nor to the exponent-flipped value `1/7`; it tends to `4/7`. -/
theorem ctl_neg_wrong_limits :
    ¬ Tendsto (fun n => kPow (moranP 3 2) n 1 (Fin.last 3)) atTop (𝓝 (1 / 3)) ∧
    ¬ Tendsto (fun n => kPow (moranP 3 2) n 1 (Fin.last 3)) atTop (𝓝 (1 / 7)) := by
  have h2 := (moran_fixation_limit (N := 3) (r := 2) (by norm_num) (by norm_num) 1).1
  have hv : moranFix 3 2 ((1 : Fin 4) : ℕ) = 4 / 7 := by
    rw [show ((1 : Fin 4) : ℕ) = 1 from rfl]
    norm_num [moranFix, geomS, Finset.sum_range_succ]
  rw [hv] at h2
  constructor
  · intro h
    have := tendsto_nhds_unique h h2
    norm_num at this
  · intro h
    have := tendsto_nhds_unique h h2
    norm_num at this

end ProgPopGen

#print axioms ProgPopGen.moran_absorbing
#print axioms ProgPopGen.moran_interiorHarmonic
#print axioms ProgPopGen.moran_var_pos
#print axioms ProgPopGen.moran_fixation_limit
#print axioms ProgPopGen.moran_fixation_limit_closed
#print axioms ProgPopGen.moran_fixation_limit_neutral
#print axioms ProgPopGen.moran_fixation_unique
#print axioms ProgPopGen.ctl_neg_wrong_limits

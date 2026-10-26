import Programmes.PopGen.Absorption
import Autoresearch.PopGenWrightFisher

/-!
# Neutral Wright–Fisher: absorption-limit form of the fixation probability (P7, closes the H11 gap)

H11 (`Autoresearch/PopGenWrightFisher.lean`) kernel-checked the harmonic characterisation: `i/N` is the
unique function harmonic at interior states with boundary values `0`, `1`. It did not connect this to
an absorption probability. Here, with the general theorem of `Programmes/PopGen/Absorption.lean`:

* `wright_fisher_fixation_limit`: `Pⁿ(i,N) → i/N` and `Pⁿ(i,0) → 1 - i/N` as `n → ∞`, where `Pⁿ` is
  the `n`-step kernel of the neutral Wright–Fisher chain (`kPow (wfP N) n`).
* `wright_fisher_fixation_mono`: `Pⁿ(i,N)` is non-decreasing in `n`.

`Pⁿ(i,N)` is the probability of being at `N` at time `n` (by the definition of the `n`-step law,
Chapman–Kolmogorov); `N` is absorbing, so it is the probability of having fixed by time `n`, and its
limit is the fixation probability. That identification with a path measure is the standard definition
and is not itself formalised (no path space is built). Known result; novelty not claimed.
-/

namespace ProgPopGen

open Finset Filter Topology AutoresearchPopGen

theorem wf_absorbing {N : ℕ} (hN : 0 < N) : AbsorbingChain (wfP N) :=
  ⟨wfP_nonneg hN, wfP_sum_one, wfP_zero, wfP_last hN⟩

theorem wf_interiorHarmonic {N : ℕ} (hN : 0 < N) :
    InteriorHarmonic (wfP N) (fun j => ((j : ℕ) : ℝ) / N) :=
  fun i _ _ => frac_harmonic hN i

/-- One-step variance of `j/N` is positive at interior states: the jump to `0` alone contributes
`P(k,0) (k/N)² > 0`. -/
theorem wf_var_pos {N : ℕ} (hN : 0 < N) (k : Fin (N + 1)) (hk1 : 0 < (k : ℕ)) (hk2 : (k : ℕ) < N) :
    0 < ∑ j, wfP N k j * (((j : ℕ) : ℝ) / N - ((k : ℕ) : ℝ) / N) ^ 2 := by
  have hpos := wfP_interior_zero_pos k hk2
  have hk : (0 : ℝ) < ((k : ℕ) : ℝ) / N :=
    div_pos (by exact_mod_cast hk1) (by exact_mod_cast hN)
  have hterm : 0 < wfP N k 0 * (((((0 : Fin (N + 1)) : ℕ) : ℝ)) / N - ((k : ℕ) : ℝ) / N) ^ 2 := by
    apply mul_pos hpos
    rw [show ((((0 : Fin (N + 1)) : ℕ) : ℝ)) = 0 by simp, zero_div, zero_sub, neg_sq]
    exact pow_pos hk 2
  exact lt_of_lt_of_le hterm (Finset.single_le_sum
    (f := fun j => wfP N k j * (((j : ℕ) : ℝ) / N - ((k : ℕ) : ℝ) / N) ^ 2)
    (fun j _ => mul_nonneg (wfP_nonneg hN k j) (sq_nonneg _)) (Finset.mem_univ 0))

/-- **Closes the H11 gap.** For the neutral Wright–Fisher chain on `{0,…,N}`, the `n`-step probability
of being at `N` tends to `i/N`, and of being at `0` tends to `1 - i/N`. -/
theorem wright_fisher_fixation_limit {N : ℕ} (hN : 0 < N) (i : Fin (N + 1)) :
    Tendsto (fun n => kPow (wfP N) n i (Fin.last N)) atTop (𝓝 (((i : ℕ) : ℝ) / N)) ∧
    Tendsto (fun n => kPow (wfP N) n i 0) atTop (𝓝 (1 - ((i : ℕ) : ℝ) / N)) := by
  have hNr : (N : ℝ) ≠ 0 := by exact_mod_cast hN.ne'
  exact fixation_limit (wf_absorbing hN) hN (wf_interiorHarmonic hN) (by simp)
    (by simp [hNr]) (wf_var_pos hN) i

theorem wright_fisher_fixation_mono {N : ℕ} (hN : 0 < N) (n : ℕ) (i : Fin (N + 1)) :
    kPow (wfP N) n i (Fin.last N) ≤ kPow (wfP N) (n + 1) i (Fin.last N) :=
  kPow_last_mono (wf_absorbing hN) n i

/-- Positive control: after two generations from `i = 1`, `N = 2`, the probability of having fixed is
`3/8` (not yet the limit `1/2`): `P²(1,2) = P(1,1)P(1,2) + P(1,2) = 1/8 + 1/4`. -/
theorem ctl_pos_wf_two_steps : kPow (wfP 2) 2 1 (Fin.last 2) = 3 / 8 := by
  simp [kPow, Fin.sum_univ_three, wfP]
  norm_num

end ProgPopGen

#print axioms ProgPopGen.wf_absorbing
#print axioms ProgPopGen.wf_var_pos
#print axioms ProgPopGen.wright_fisher_fixation_limit
#print axioms ProgPopGen.wright_fisher_fixation_mono
#print axioms ProgPopGen.ctl_pos_wf_two_steps

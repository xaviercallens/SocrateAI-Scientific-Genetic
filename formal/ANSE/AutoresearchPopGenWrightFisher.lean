import Mathlib.Data.Nat.Choose.Sum
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Tactic

/-!
# Neutral Wright–Fisher fixation probability (H11, autoresearch 2026-10-25)

First formalization attempt of a KNOWN classical result (Fisher 1930, Wright 1931; e.g. Ewens,
*Mathematical Population Genetics*, §3): in the neutral Wright–Fisher chain on `{0,…,N}` the
probability of fixation (absorption at `N`) from state `i` is `i/N`.

Formulation achieved here (harmonic-function characterisation): the transition kernel
`P(i,j) = C(N,j) (i/N)^j (1-i/N)^(N-j)` is stochastic, `0` and `N` are absorbing, the state is a
martingale (`∑ⱼ P(i,j) j = i`), `h(i) = i/N` is harmonic, and it is the UNIQUE function that is
harmonic at interior states with `h 0 = 0`, `h N = 1` (maximum principle; every interior state
reaches `0` in one step with positive probability).
-/

namespace AutoresearchPopGen

open Finset

/-- Wright–Fisher transition probability `P(i → j)` on `Fin (N+1)`. -/
noncomputable def wfP (N : ℕ) (i j : Fin (N + 1)) : ℝ :=
  (N.choose (j : ℕ) : ℝ) * ((i : ℕ) / N : ℝ) ^ (j : ℕ) * (1 - ((i : ℕ) / N : ℝ)) ^ (N - (j : ℕ))

theorem frac_mem {N : ℕ} (hN : 0 < N) (i : Fin (N + 1)) :
    0 ≤ ((i : ℕ) / N : ℝ) ∧ ((i : ℕ) / N : ℝ) ≤ 1 := by
  have hNr : (0 : ℝ) < N := by exact_mod_cast hN
  have hi : ((i : ℕ) : ℝ) ≤ N := by exact_mod_cast Nat.lt_succ_iff.mp i.isLt
  exact ⟨by positivity, (div_le_one hNr).mpr hi⟩

theorem wfP_nonneg {N : ℕ} (hN : 0 < N) (i j : Fin (N + 1)) : 0 ≤ wfP N i j := by
  obtain ⟨h0, h1⟩ := frac_mem hN i
  have : 0 ≤ 1 - ((i : ℕ) / N : ℝ) := by linarith
  unfold wfP
  positivity

/-- Each row of the kernel sums to one. -/
theorem wfP_sum_one {N : ℕ} (i : Fin (N + 1)) : ∑ j, wfP N i j = 1 := by
  set x : ℝ := ((i : ℕ) / N : ℝ)
  have h := add_pow x (1 - x) N
  rw [show x + (1 - x) = 1 by ring, one_pow] at h
  rw [h, ← Fin.sum_univ_eq_sum_range (fun k => x ^ k * (1 - x) ^ (N - k) * (N.choose k : ℝ))]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  unfold wfP
  ring

/-- Binomial mean identity. -/
theorem binom_mean (x y : ℝ) (n : ℕ) :
    ∑ k ∈ range (n + 2), ((n + 1).choose k : ℝ) * x ^ k * y ^ (n + 1 - k) * (k : ℝ)
      = (n + 1) * x * (x + y) ^ n := by
  rw [Finset.sum_range_succ', add_pow, Finset.mul_sum]
  simp only [Nat.cast_zero, mul_zero, add_zero]
  refine Finset.sum_congr rfl (fun k _ => ?_)
  have h := Nat.add_one_mul_choose_eq n k
  have h' : ((n + 1).choose (k + 1) : ℝ) * ((k : ℝ) + 1) = ((n : ℝ) + 1) * (n.choose k : ℝ) := by
    exact_mod_cast h.symm
  rw [show n + 1 - (k + 1) = n - k by omega, pow_succ]
  push_cast
  linear_combination (x ^ k * y ^ (n - k) * x) * h'

/-- Martingale property: `∑ⱼ P(i,j) j = i` (expected next state equals current state). -/
theorem wfP_mean {N : ℕ} (hN : 0 < N) (i : Fin (N + 1)) :
    ∑ j, wfP N i j * ((j : ℕ) : ℝ) = ((i : ℕ) : ℝ) := by
  obtain ⟨n, rfl⟩ := Nat.exists_eq_succ_of_ne_zero hN.ne'
  set x : ℝ := ((i : ℕ) / (n + 1 : ℕ) : ℝ) with hx
  have hm := binom_mean x (1 - x) n
  rw [show x + (1 - x) = 1 by ring, one_pow, mul_one] at hm
  have e : ∑ j : Fin (n + 1 + 1), wfP (n + 1) i j * ((j : ℕ) : ℝ)
      = ∑ k ∈ range (n + 2), ((n + 1).choose k : ℝ) * x ^ k * (1 - x) ^ (n + 1 - k) * (k : ℝ) := by
    rw [← Fin.sum_univ_eq_sum_range
      (fun k => ((n + 1).choose k : ℝ) * x ^ k * (1 - x) ^ (n + 1 - k) * (k : ℝ))]
    rfl
  rw [e, hm, hx]
  have : ((n + 1 : ℕ) : ℝ) ≠ 0 := by positivity
  push_cast at this ⊢
  field_simp

/-- State `0` is absorbing. -/
theorem wfP_zero {N : ℕ} (j : Fin (N + 1)) : wfP N 0 j = if j = 0 then 1 else 0 := by
  unfold wfP
  by_cases hj : j = 0
  · subst hj; simp
  · have : (j : ℕ) ≠ 0 := fun h => hj (Fin.ext h)
    simp [hj, this]

/-- State `N` is absorbing. -/
theorem wfP_last {N : ℕ} (hN : 0 < N) (j : Fin (N + 1)) :
    wfP N (Fin.last N) j = if j = Fin.last N then 1 else 0 := by
  have hNr : (N : ℝ) ≠ 0 := by exact_mod_cast hN.ne'
  unfold wfP
  simp only [Fin.val_last, div_self hNr, sub_self, one_pow, mul_one]
  by_cases hj : j = Fin.last N
  · subst hj; simp
  · have hlt : (j : ℕ) < N := Fin.val_lt_last hj
    have : N - (j : ℕ) ≠ 0 := by omega
    simp [hj, zero_pow this]

/-- From an interior state the chain jumps to `0` with positive probability. -/
theorem wfP_interior_zero_pos {N : ℕ} (i : Fin (N + 1)) (hi : (i : ℕ) < N) :
    0 < wfP N i 0 := by
  have hN : 0 < N := lt_of_le_of_lt (Nat.zero_le _) hi
  have hNr : (0 : ℝ) < N := by exact_mod_cast hN
  have hir : ((i : ℕ) : ℝ) < N := by exact_mod_cast hi
  have : 0 < 1 - ((i : ℕ) / N : ℝ) := by
    rw [sub_pos, div_lt_one hNr]; exact hir
  unfold wfP
  simp only [Fin.val_zero, Nat.choose_zero_right, Nat.cast_one, pow_zero, one_mul, Nat.sub_zero]
  positivity

/-- `h(i) = i/N` is harmonic for the kernel (at every state). -/
theorem frac_harmonic {N : ℕ} (hN : 0 < N) (i : Fin (N + 1)) :
    ∑ j, wfP N i j * (((j : ℕ) : ℝ) / N) = ((i : ℕ) : ℝ) / N := by
  rw [← wfP_mean hN i, Finset.sum_div]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  ring

/-- Maximum principle: a function harmonic at interior states and vanishing at `0` and `N`
is `≤ 0` everywhere. -/
theorem max_principle {N : ℕ} (hN : 0 < N) (d : Fin (N + 1) → ℝ) (h0 : d 0 = 0)
    (hN' : d (Fin.last N) = 0)
    (hharm : ∀ i : Fin (N + 1), 0 < (i : ℕ) → (i : ℕ) < N → d i = ∑ j, wfP N i j * d j) :
    ∀ i, d i ≤ 0 := by
  obtain ⟨m, -, hm⟩ := Finset.exists_max_image Finset.univ d Finset.univ_nonempty
  suffices hM : d m ≤ 0 from fun i => le_trans (hm i (Finset.mem_univ _)) hM
  by_cases hm0 : (m : ℕ) = 0
  · have : m = 0 := Fin.ext hm0
    rw [this, h0]
  by_cases hmN : (m : ℕ) = N
  · have : m = Fin.last N := Fin.ext hmN
    rw [this, hN']
  have hpos : 0 < (m : ℕ) := Nat.pos_of_ne_zero hm0
  have hlt : (m : ℕ) < N := lt_of_le_of_ne (Nat.lt_succ_iff.mp m.isLt) hmN
  have hsum : ∑ j, wfP N m j * (d m - d j) = 0 := by
    have e1 : ∑ j, wfP N m j * (d m - d j) = d m * ∑ j, wfP N m j - ∑ j, wfP N m j * d j := by
      rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
      refine Finset.sum_congr rfl (fun j _ => ?_)
      ring
    rw [e1, wfP_sum_one, ← hharm m hpos hlt]
    ring
  have hnn : ∀ j ∈ (Finset.univ : Finset (Fin (N + 1))), 0 ≤ wfP N m j * (d m - d j) :=
    fun j _ => mul_nonneg (wfP_nonneg hN m j) (sub_nonneg.mpr (hm j (Finset.mem_univ _)))
  have hz := (Finset.sum_eq_zero_iff_of_nonneg hnn).mp hsum 0 (Finset.mem_univ _)
  have hp := wfP_interior_zero_pos m hlt
  rcases mul_eq_zero.mp hz with h | h
  · exact absurd h hp.ne'
  · rw [h0] at h; linarith

/-- Uniqueness: `i/N` is the only function harmonic at interior states with boundary values
`0` at `0` and `1` at `N`. Hence the fixation probability (which is such a function by first-step
analysis) equals `i/N`. -/
theorem harmonic_unique {N : ℕ} (hN : 0 < N) (g : Fin (N + 1) → ℝ) (h0 : g 0 = 0)
    (h1 : g (Fin.last N) = 1)
    (hharm : ∀ i : Fin (N + 1), 0 < (i : ℕ) → (i : ℕ) < N → g i = ∑ j, wfP N i j * g j) :
    ∀ i, g i = ((i : ℕ) : ℝ) / N := by
  have hNr : (N : ℝ) ≠ 0 := by exact_mod_cast hN.ne'
  set f : Fin (N + 1) → ℝ := fun i => ((i : ℕ) : ℝ) / N with hf
  have hfh : ∀ i, ∑ j, wfP N i j * f j = f i := fun i => frac_harmonic hN i
  have hd : ∀ s : ℝ, s = 1 ∨ s = -1 → ∀ i, s * (g i - f i) ≤ 0 := by
    intro s hs
    apply max_principle hN (fun i => s * (g i - f i))
    · simp [h0, hf]
    · simp [h1, hf, hNr]
    · intro i hi0 hiN
      have e : ∑ j, wfP N i j * (s * (g j - f j))
          = s * ((∑ j, wfP N i j * g j) - ∑ j, wfP N i j * f j) := by
        rw [mul_sub, Finset.mul_sum, Finset.mul_sum, ← Finset.sum_sub_distrib]
        refine Finset.sum_congr rfl (fun j _ => ?_)
        ring
      show s * (g i - f i) = ∑ j, wfP N i j * (s * (g j - f j))
      rw [e, hfh i, ← hharm i hi0 hiN]
  intro i
  have a := hd 1 (Or.inl rfl) i
  have b := hd (-1) (Or.inr rfl) i
  show g i = f i
  linarith

/-- Combined statement (H11, formulation (a)). -/
theorem wright_fisher_fixation_harmonic {N : ℕ} (hN : 0 < N) :
    (∀ i, ∑ j, wfP N i j * (((j : ℕ) : ℝ) / N) = ((i : ℕ) : ℝ) / N) ∧
    (∀ g : Fin (N + 1) → ℝ, g 0 = 0 → g (Fin.last N) = 1 →
      (∀ i : Fin (N + 1), 0 < (i : ℕ) → (i : ℕ) < N → g i = ∑ j, wfP N i j * g j) →
      ∀ i, g i = ((i : ℕ) : ℝ) / N) :=
  ⟨frac_harmonic hN, fun g h0 h1 hh => harmonic_unique hN g h0 h1 hh⟩

/-! ## Controls on the model (informative)

Positive control (definition fidelity): `wfP 2 1 1 = 1/2` — real division, `C(2,1)(1/2)(1/2)`.

Negative control (neutrality is essential): the Wright–Fisher kernel with selection samples with
probability `ψ ≠ i/N`. For `N = 2`, `i = 1`, `ψ = 2/3` the expected next frequency is
`2/3 ≠ 1/2`, so `h(i) = i/N` is NOT harmonic for the selection kernel and the fixation
probability is not `i/N` there. The harmonicity proved above therefore really uses neutrality. -/

/-- Wright–Fisher kernel with a general sampling probability `ψ` (selection: `ψ ≠ i/N`). -/
noncomputable def wfPsel (N : ℕ) (ψ : ℝ) (j : Fin (N + 1)) : ℝ :=
  (N.choose (j : ℕ) : ℝ) * ψ ^ (j : ℕ) * (1 - ψ) ^ (N - (j : ℕ))

theorem ctl_positive_wfP_value : wfP 2 1 1 = 1 / 2 := by
  norm_num [wfP]

theorem ctl_negative_selection_not_harmonic :
    ∑ j : Fin 3, wfPsel 2 (2 / 3) j * (((j : ℕ) : ℝ) / 2) = 2 / 3 ∧ (2 / 3 : ℝ) ≠ 1 / 2 := by
  refine ⟨?_, by norm_num⟩
  simp [Fin.sum_univ_three, wfPsel, Nat.choose]
  norm_num

/-! ## Planted controls (must be caught) -/

end AutoresearchPopGen

#print axioms AutoresearchPopGen.wfP_nonneg
#print axioms AutoresearchPopGen.wfP_sum_one
#print axioms AutoresearchPopGen.binom_mean
#print axioms AutoresearchPopGen.wfP_mean
#print axioms AutoresearchPopGen.wfP_zero
#print axioms AutoresearchPopGen.wfP_last
#print axioms AutoresearchPopGen.wfP_interior_zero_pos
#print axioms AutoresearchPopGen.frac_harmonic
#print axioms AutoresearchPopGen.max_principle
#print axioms AutoresearchPopGen.harmonic_unique
#print axioms AutoresearchPopGen.wright_fisher_fixation_harmonic
#print axioms AutoresearchPopGen.ctl_positive_wfP_value
#print axioms AutoresearchPopGen.ctl_negative_selection_not_harmonic

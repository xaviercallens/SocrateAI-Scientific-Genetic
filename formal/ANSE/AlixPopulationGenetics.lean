import Mathlib.Tactic.Ring
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Linarith
import Mathlib.Data.Nat.Choose.Sum
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Algebra.Order.BigOperators.Group.Finset

/-!
# Two population-genetics claims behind Kernel-Gated Evolution (H11, H12)

Both are stated with the scope corrections an internal review asked for.

* **H12** -- one locus, two alleles, haploid selection with constant positive
  fitnesses: mean fitness never decreases from one generation to the next, with
  the exact increment `p(1-p)(w₁-w₂)²/w̄`. This is NOT Fisher's general
  "fundamental theorem"; it is the one-locus two-allele haploid case, proved
  exactly over `ℚ`.
* **H11** -- for a finite chain on `{0,…,N}` whose interior states each reach
  state `0` in one step with positive probability (true of Wright–Fisher), a
  harmonic function is determined by its boundary values; and for Wright–Fisher,
  `i/N` is harmonic with boundary values `0` and `1`. So `i/N` is THE harmonic
  solution. The further step identifying it with the fixation probability (an
  optional-stopping / absorption argument) is NOT formalised here.

Verification contract: compile, then `#print axioms`; `sorryAx` rejected; axioms
whitelisted to {propext, Classical.choice, Quot.sound}.
-/

open Finset

namespace ANSE.Alix

/-! ## H12: selection never lowers mean fitness (one locus, two alleles, haploid) -/

/-- Mean fitness with allele frequency `p` and allele fitnesses `w₁`, `w₂`. -/
def meanFit (w₁ w₂ p : ℚ) : ℚ := p * w₁ + (1 - p) * w₂

/-- Allele frequency after one round of haploid selection. -/
def selectP (w₁ w₂ p : ℚ) : ℚ := p * w₁ / meanFit w₁ w₂ p

/-- The exact increment: `w̄' - w̄ = p(1-p)(w₁-w₂)² / w̄`. -/
theorem aurelia_mean_fitness_increment (w₁ w₂ p : ℚ) (hm : meanFit w₁ w₂ p ≠ 0) :
    meanFit w₁ w₂ (selectP w₁ w₂ p) - meanFit w₁ w₂ p
      = p * (1 - p) * (w₁ - w₂) ^ 2 / meanFit w₁ w₂ p := by
  unfold selectP
  unfold meanFit at hm ⊢
  field_simp
  ring

/-- H12: with positive fitnesses and `0 ≤ p ≤ 1`, mean fitness does not decrease. -/
theorem aurelia_mean_fitness_nondecreasing (w₁ w₂ p : ℚ)
    (h₁ : 0 < w₁) (h₂ : 0 < w₂) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    meanFit w₁ w₂ p ≤ meanFit w₁ w₂ (selectP w₁ w₂ p) := by
  have hm : 0 < meanFit w₁ w₂ p := by
    unfold meanFit
    rcases eq_or_lt_of_le hp0 with h | h
    · subst h; simp; exact h₂
    · have : 0 ≤ (1 - p) * w₂ := mul_nonneg (by linarith) h₂.le
      nlinarith
  have hinc := aurelia_mean_fitness_increment w₁ w₂ p hm.ne'
  have hnum : 0 ≤ p * (1 - p) * (w₁ - w₂) ^ 2 :=
    mul_nonneg (mul_nonneg hp0 (by linarith)) (sq_nonneg _)
  have : 0 ≤ meanFit w₁ w₂ (selectP w₁ w₂ p) - meanFit w₁ w₂ p := by
    rw [hinc]; exact div_nonneg hnum hm.le
  linarith

/-- The numerator of the increment is the variance of fitness in the population:
`p(1-p)(w₁-w₂)² = p(w₁-w̄)² + (1-p)(w₂-w̄)²`. So H12 reads: one round of selection
raises mean fitness by exactly `Var(w) / w̄` -- progress is proportional to the
fitness variance present, and stops when that variance is gone. -/
theorem kal_increment_is_variance (w₁ w₂ p : ℚ) :
    p * (1 - p) * (w₁ - w₂) ^ 2
      = p * (w₁ - meanFit w₁ w₂ p) ^ 2 + (1 - p) * (w₂ - meanFit w₁ w₂ p) ^ 2 := by
  unfold meanFit; ring

/-! ## H11: a maximum principle, and `i/N` is the harmonic solution for Wright–Fisher -/

/-- One-sided maximum principle. If `d` vanishes at `0` and at `N`, satisfies the
mean-value property at every interior state, and every interior state puts
positive weight on state `0`, then `d ≤ 0` everywhere. -/
theorem ilene_harmonic_le_zero (n : ℕ) (P : Fin (n + 1) → Fin (n + 1) → ℚ) (d : Fin (n + 1) → ℚ)
    (hP : ∀ i j, 0 ≤ P i j) (hrow : ∀ i, ∑ j, P i j = 1)
    (hpos : ∀ i, i ≠ 0 → i ≠ Fin.last n → 0 < P i 0)
    (hharm : ∀ i, i ≠ 0 → i ≠ Fin.last n → d i = ∑ j, P i j * d j)
    (h0 : d 0 = 0) (hN : d (Fin.last n) = 0) :
    ∀ i, d i ≤ 0 := by
  obtain ⟨k, -, hk⟩ := Finset.exists_max_image Finset.univ d Finset.univ_nonempty
  suffices d k ≤ 0 from fun i => le_trans (hk i (Finset.mem_univ i)) this
  by_cases hk0 : k = 0
  · rw [hk0, h0]
  by_cases hkN : k = Fin.last n
  · rw [hkN, hN]
  -- Weighted deviation from the maximum sums to zero ...
  have hsum : ∑ j, P k j * (d k - d j) = 0 := by
    have e1 : ∑ j, P k j * (d k - d j) = (∑ j, P k j) * d k - ∑ j, P k j * d j := by
      rw [Finset.sum_mul, ← Finset.sum_sub_distrib]
      refine Finset.sum_congr rfl fun j _ => by ring
    rw [e1, hrow k, ← hharm k hk0 hkN]; ring
  -- ... and every term is nonnegative, so the term at state 0 vanishes.
  have hterms := (Finset.sum_eq_zero_iff_of_nonneg
    (fun j _ => mul_nonneg (hP k j) (by linarith [hk j (Finset.mem_univ j)]))).1 hsum
  have h0term : P k 0 * (d k - d 0) = 0 := hterms 0 (Finset.mem_univ 0)
  rcases mul_eq_zero.1 h0term with h | h
  · exact absurd h (hpos k hk0 hkN).ne'
  · linarith

/-- Uniqueness: two harmonic functions with the same boundary values coincide. -/
theorem ilene_harmonic_unique (n : ℕ) (P : Fin (n + 1) → Fin (n + 1) → ℚ) (f g : Fin (n + 1) → ℚ)
    (hP : ∀ i j, 0 ≤ P i j) (hrow : ∀ i, ∑ j, P i j = 1)
    (hpos : ∀ i, i ≠ 0 → i ≠ Fin.last n → 0 < P i 0)
    (hf : ∀ i, i ≠ 0 → i ≠ Fin.last n → f i = ∑ j, P i j * f j)
    (hg : ∀ i, i ≠ 0 → i ≠ Fin.last n → g i = ∑ j, P i j * g j)
    (h0 : f 0 = g 0) (hN : f (Fin.last n) = g (Fin.last n)) :
    ∀ i, f i = g i := by
  have harm (s : ℚ) : ∀ i, i ≠ 0 → i ≠ Fin.last n →
      s * (f i - g i) = ∑ j, P i j * (s * (f j - g j)) := by
    intro i hi0 hiN
    rw [hf i hi0 hiN, hg i hi0 hiN, mul_sub, Finset.mul_sum, Finset.mul_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun j _ => by ring
  have up := ilene_harmonic_le_zero n P (fun i => 1 * (f i - g i)) hP hrow hpos (harm 1)
    (by simp [h0]) (by simp [hN])
  have dn := ilene_harmonic_le_zero n P (fun i => (-1) * (f i - g i)) hP hrow hpos (harm (-1))
    (by simp [h0]) (by simp [hN])
  intro i
  have a := up i; have b := dn i
  simp only [one_mul, neg_one_mul] at a b
  linarith

/-- The Wright–Fisher transition from `i` copies to `j` copies in a population of `N`. -/
def wfP (N : ℕ) (i j : Fin (N + 1)) : ℚ :=
  (N.choose j : ℚ) * ((i : ℚ) / N) ^ (j : ℕ) * (1 - (i : ℚ) / N) ^ (N - j)

theorem wfP_nonneg (N : ℕ) (i j : Fin (N + 1)) : 0 ≤ wfP N i j := by
  unfold wfP
  have hi : (i : ℚ) ≤ N := by exact_mod_cast Nat.lt_succ_iff.mp i.isLt
  have h1 : 0 ≤ 1 - (i : ℚ) / N := by
    rcases Nat.eq_zero_or_pos N with h | h
    · simp [h]
    · rw [sub_nonneg, div_le_one (by exact_mod_cast h)]; exact hi
  positivity

theorem wfP_row (N : ℕ) (i : Fin (N + 1)) : ∑ j, wfP N i j = 1 := by
  unfold wfP
  have h := (add_pow ((i : ℚ) / N) (1 - (i : ℚ) / N) N).symm
  rw [show (i : ℚ) / N + (1 - (i : ℚ) / N) = 1 by ring, one_pow] at h
  rw [Fin.sum_univ_eq_sum_range
    (fun j => (N.choose j : ℚ) * ((i : ℚ) / N) ^ j * (1 - (i : ℚ) / N) ^ (N - j)) (N + 1)]
  exact (Finset.sum_congr rfl fun j _ => by ring).trans h

/-- Every interior state reaches state `0` with positive probability. -/
theorem wfP_pos_zero (N : ℕ) (i : Fin (N + 1)) (hiN : i ≠ Fin.last N) : 0 < wfP N i 0 := by
  unfold wfP
  have hlt : (i : ℕ) < N := by
    have := Fin.val_lt_last hiN; simpa using this
  have hNpos : (0 : ℚ) < N := by exact_mod_cast lt_of_le_of_lt (Nat.zero_le _) hlt
  have : (i : ℚ) / N < 1 := by rw [div_lt_one hNpos]; exact_mod_cast hlt
  simp only [Fin.val_zero, Nat.choose_zero_right, Nat.cast_one, pow_zero, one_mul, Nat.sub_zero]
  exact pow_pos (by linarith) _

/-- The binomial mean, `∑ⱼ C(N,j) pʲ qᴺ⁻ʲ · j = N p (p+q)ᴺ⁻¹`, in the form needed. -/
theorem binomial_mean (n : ℕ) (p q : ℚ) :
    ∑ j ∈ Finset.range (n + 2), ((n + 1).choose j : ℚ) * p ^ j * q ^ (n + 1 - j) * j
      = (n + 1) * p * (p + q) ^ n := by
  rw [Finset.sum_range_succ', add_pow]
  simp only [Nat.cast_zero, mul_zero, add_zero]
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  have hc : ((n + 1).choose (j + 1) : ℚ) * (j + 1) = (n + 1) * (n.choose j : ℚ) := by
    have := Nat.add_one_mul_choose_eq n j
    exact_mod_cast this.symm
  rw [Nat.succ_sub_succ, Nat.cast_succ]
  calc ((n + 1).choose (j + 1) : ℚ) * p ^ (j + 1) * q ^ (n - j) * (j + 1)
      = (((n + 1).choose (j + 1) : ℚ) * (j + 1)) * p * (p ^ j * q ^ (n - j)) := by ring
    _ = (n + 1) * p * (p ^ j * q ^ (n - j) * (n.choose j : ℚ)) := by rw [hc]; ring

/-- `j/N` satisfies the mean-value property for Wright–Fisher. -/
theorem ilene_wf_harmonic (n : ℕ) (i : Fin (n + 2)) :
    (i : ℚ) / (n + 1 : ℕ) = ∑ j, wfP (n + 1) i j * ((j : ℚ) / (n + 1 : ℕ)) := by
  unfold wfP
  have hN : ((n + 1 : ℕ) : ℚ) ≠ 0 := by positivity
  rw [Fin.sum_univ_eq_sum_range (fun j : ℕ =>
      ((n + 1).choose j : ℚ) * ((i : ℚ) / (n + 1 : ℕ)) ^ j
        * (1 - (i : ℚ) / (n + 1 : ℕ)) ^ (n + 1 - j) * ((j : ℚ) / (n + 1 : ℕ))) (n + 2)]
  have hm := binomial_mean n ((i : ℚ) / (n + 1 : ℕ)) (1 - (i : ℚ) / (n + 1 : ℕ))
  rw [show (i : ℚ) / (n + 1 : ℕ) + (1 - (i : ℚ) / (n + 1 : ℕ)) = 1 by ring, one_pow] at hm
  have : ∑ j ∈ Finset.range (n + 2), ((n + 1).choose j : ℚ) * ((i : ℚ) / (n + 1 : ℕ)) ^ j
        * (1 - (i : ℚ) / (n + 1 : ℕ)) ^ (n + 1 - j) * ((j : ℚ) / (n + 1 : ℕ))
      = (∑ j ∈ Finset.range (n + 2), ((n + 1).choose j : ℚ) * ((i : ℚ) / (n + 1 : ℕ)) ^ j
        * (1 - (i : ℚ) / (n + 1 : ℕ)) ^ (n + 1 - j) * j) / (n + 1 : ℕ) := by
    simp only [div_eq_mul_inv]
    rw [Finset.sum_mul]; refine Finset.sum_congr rfl fun j _ => by ring
  rw [this, hm]
  push_cast
  field_simp

/-- H11: for Wright–Fisher, any harmonic function equal to `0` at `0` and `1` at
`N` is `i/N`. -/
theorem xav_wf_harmonic_solution_unique (n : ℕ) (h : Fin (n + 2) → ℚ)
    (hharm : ∀ i, i ≠ 0 → i ≠ Fin.last (n + 1) → h i = ∑ j, wfP (n + 1) i j * h j)
    (h0 : h 0 = 0) (hN : h (Fin.last (n + 1)) = 1) :
    ∀ i, h i = (i : ℚ) / (n + 1 : ℕ) := by
  refine ilene_harmonic_unique (n + 1) (wfP (n + 1)) h (fun i => (i : ℚ) / (n + 1 : ℕ))
    (wfP_nonneg (n + 1)) (wfP_row (n + 1)) (fun i _ hiN => wfP_pos_zero (n + 1) i hiN)
    hharm (fun i _ _ => ilene_wf_harmonic n i) ?_ ?_
  · simp [h0]
  · rw [hN]; simp only [Fin.val_last]; push_cast
    exact (div_self (by positivity)).symm

/-! ## Toward fixation: a finite-horizon bound, with no limit taken

The step from "i/N is harmonic" to "i/N is the fixation probability" needs an
absorption limit. Half of it needs none. If `h` has the mean-value property at
EVERY state (boundaries included, as `j/N` does for Wright–Fisher), then the
expected value of `h` is conserved by every step, so after any number of
generations it is still `h i`. Because `h ≥ 0` and `h N = 1`, the mass sitting
on state `N` can never exceed it: the probability of having fixed by time `t` is
at most `i/N`, for every `t`. -/

/-- One step of the chain acting on a distribution over states. -/
def stepDist {n : ℕ} (P : Fin (n + 1) → Fin (n + 1) → ℚ) (μ : Fin (n + 1) → ℚ) :
    Fin (n + 1) → ℚ :=
  fun j => ∑ i, μ i * P i j

/-- A harmonic function's expectation is conserved by one step. -/
theorem alix_expectation_step {n : ℕ} (P : Fin (n + 1) → Fin (n + 1) → ℚ) (h μ : Fin (n + 1) → ℚ)
    (hh : ∀ i, h i = ∑ j, P i j * h j) :
    ∑ j, stepDist P μ j * h j = ∑ i, μ i * h i := by
  unfold stepDist
  calc ∑ j, (∑ i, μ i * P i j) * h j
      = ∑ j, ∑ i, μ i * (P i j * h j) := by
        refine Finset.sum_congr rfl fun j _ => ?_
        rw [Finset.sum_mul]; refine Finset.sum_congr rfl fun i _ => by ring
    _ = ∑ i, ∑ j, μ i * (P i j * h j) := Finset.sum_comm
    _ = ∑ i, μ i * h i := by
        refine Finset.sum_congr rfl fun i _ => ?_
        rw [← Finset.mul_sum, ← hh i]

/-- ... and so by every number of steps. -/
theorem alix_expectation_iterate {n : ℕ} (P : Fin (n + 1) → Fin (n + 1) → ℚ) (h μ : Fin (n + 1) → ℚ)
    (hh : ∀ i, h i = ∑ j, P i j * h j) (t : ℕ) :
    ∑ j, (stepDist P)^[t] μ j * h j = ∑ i, μ i * h i := by
  induction t with
  | zero => rfl
  | succ t ih => rw [Function.iterate_succ_apply', alix_expectation_step P h _ hh, ih]

/-- Steps of a nonnegative matrix keep a distribution nonnegative. -/
theorem alix_iterate_nonneg {n : ℕ} (P : Fin (n + 1) → Fin (n + 1) → ℚ) (μ : Fin (n + 1) → ℚ)
    (hP : ∀ i j, 0 ≤ P i j) (hμ : ∀ i, 0 ≤ μ i) (t : ℕ) :
    ∀ j, 0 ≤ (stepDist P)^[t] μ j := by
  induction t with
  | zero => exact hμ
  | succ t ih =>
    intro j
    rw [Function.iterate_succ_apply']
    exact Finset.sum_nonneg fun i _ => mul_nonneg (ih i) (hP i j)

/-- Wright–Fisher: `j/N` is harmonic at every state, boundaries included. -/
theorem ilene_wf_harmonic_everywhere (n : ℕ) :
    ∀ i : Fin (n + 2), (fun j : Fin (n + 2) => (j : ℚ) / (n + 1 : ℕ)) i
      = ∑ j, wfP (n + 1) i j * (fun j : Fin (n + 2) => (j : ℚ) / (n + 1 : ℕ)) j :=
  fun i => ilene_wf_harmonic n i

/-- H11, toward fixation: starting from `i` copies in a Wright–Fisher population of
`N = n+1`, the probability of having fixed by generation `t` is at most `i/N`,
for every `t`. (The matching lower bound, as `t → ∞`, is the part that needs an
absorption limit and is not formalised.) -/
theorem xav_wf_fixation_by_time_le (n : ℕ) (i : Fin (n + 2)) (t : ℕ) :
    (stepDist (wfP (n + 1)))^[t] (fun k => if k = i then 1 else 0) (Fin.last (n + 1))
      ≤ (i : ℚ) / (n + 1 : ℕ) := by
  set μt := (stepDist (wfP (n + 1)))^[t] (fun k => if k = i then 1 else 0)
  have hcons := alix_expectation_iterate (wfP (n + 1)) (fun j : Fin (n + 2) => (j : ℚ) / (n + 1 : ℕ))
    (fun k => if k = i then 1 else 0) (ilene_wf_harmonic_everywhere n) t
  have hstart : ∑ k : Fin (n + 2), (if k = i then (1 : ℚ) else 0) * ((k : ℚ) / (n + 1 : ℕ))
      = (i : ℚ) / (n + 1 : ℕ) := by
    simp [ite_mul]
  have hnn := alix_iterate_nonneg (wfP (n + 1)) (fun k => if k = i then 1 else 0)
    (wfP_nonneg (n + 1)) (fun k => by split_ifs <;> norm_num) t
  have hle : μt (Fin.last (n + 1)) * (((Fin.last (n + 1) : Fin (n + 2)) : ℚ) / (n + 1 : ℕ))
      ≤ ∑ j, μt j * ((j : ℚ) / (n + 1 : ℕ)) :=
    Finset.single_le_sum (f := fun j => μt j * ((j : ℚ) / (n + 1 : ℕ)))
      (fun j _ => mul_nonneg (hnn j) (by positivity)) (Finset.mem_univ _)
  have hone : (((Fin.last (n + 1) : Fin (n + 2)) : ℚ) / (n + 1 : ℕ)) = 1 := by
    simp only [Fin.val_last]; push_cast; exact div_self (by positivity)
  rw [hone, mul_one] at hle
  have : ∑ j, μt j * ((j : ℚ) / (n + 1 : ℕ)) = (i : ℚ) / (n + 1 : ℕ) := by
    rw [← hstart]; exact hcons
  linarith

end ANSE.Alix

#print axioms ANSE.Alix.alix_expectation_step
#print axioms ANSE.Alix.alix_expectation_iterate
#print axioms ANSE.Alix.alix_iterate_nonneg
#print axioms ANSE.Alix.ilene_wf_harmonic_everywhere
#print axioms ANSE.Alix.xav_wf_fixation_by_time_le
#print axioms ANSE.Alix.aurelia_mean_fitness_increment
#print axioms ANSE.Alix.aurelia_mean_fitness_nondecreasing
#print axioms ANSE.Alix.kal_increment_is_variance
#print axioms ANSE.Alix.ilene_harmonic_le_zero
#print axioms ANSE.Alix.ilene_harmonic_unique
#print axioms ANSE.Alix.wfP_nonneg
#print axioms ANSE.Alix.wfP_row
#print axioms ANSE.Alix.wfP_pos_zero
#print axioms ANSE.Alix.binomial_mean
#print axioms ANSE.Alix.ilene_wf_harmonic
#print axioms ANSE.Alix.xav_wf_harmonic_solution_unique

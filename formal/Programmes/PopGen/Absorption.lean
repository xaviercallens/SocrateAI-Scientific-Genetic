import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Tactic

/-!
# Absorption limits for finite chains on `{0,…,N}` with absorbing ends (P7, programmes 2026-10-25)

General, model-independent part of lane P7. `P : Fin (N+1) → Fin (N+1) → ℝ` is a stochastic kernel
(non-negative entries, rows summing to one) for which `0` and `N` are absorbing. `kPow P n` is its
`n`-step kernel, defined by the left recursion `P^(n+1)(i,j) = ∑ₖ Pⁿ(i,k) P(k,j)` (Chapman–Kolmogorov).

Main result (`fixation_limit`): if `h` is harmonic at the interior states, `h 0 = 0`, `h N = 1`, and the
one-step conditional variance `∑ⱼ P(k,j) (h j - h k)²` is positive at every interior `k`, then

* `Pⁿ(i,N) → h i` and `Pⁿ(i,0) → 1 - h i` as `n → ∞`;
* every other interior-harmonic `g` with `g 0 = 0`, `g N = 1` equals `h` (`harmonic_unique_of_var`).

Proof idea (no measure theory). `h` is preserved by `Pⁿ` (martingale). `h²` gains at least `ε > 0` per
step spent in the interior (second-moment identity), and it is bounded, so the interior masses
`Iₙ(i) = ∑_{0<k<N} Pⁿ(i,k)` have bounded partial sums, hence are summable and tend to `0`. Then
`Pⁿ(i,N) - H·Iₙ ≤ h i ≤ Pⁿ(i,N) + H·Iₙ` squeezes the limit.

This is a known textbook fact about absorbing chains; novelty is not claimed.
-/

namespace ProgPopGen

open Finset Filter Topology

variable {N : ℕ}

/-- `n`-step transition kernel, left recursion `P^(n+1)(i,j) = ∑ₖ Pⁿ(i,k) P(k,j)`. -/
noncomputable def kPow (P : Fin (N + 1) → Fin (N + 1) → ℝ) : ℕ → Fin (N + 1) → Fin (N + 1) → ℝ
  | 0 => fun i j => if i = j then 1 else 0
  | n + 1 => fun i j => ∑ k, kPow P n i k * P k j

/-- Indicator of the interior states `0 < j < N`. -/
noncomputable def intInd (N : ℕ) (j : Fin (N + 1)) : ℝ :=
  if 0 < (j : ℕ) ∧ (j : ℕ) < N then 1 else 0

/-- Mass on the interior after `n` steps from `i`. -/
noncomputable def interiorMass (P : Fin (N + 1) → Fin (N + 1) → ℝ) (n : ℕ) (i : Fin (N + 1)) : ℝ :=
  ∑ k, kPow P n i k * intInd N k

/-- A stochastic kernel on `{0,…,N}` for which `0` and `N` are absorbing. -/
structure AbsorbingChain (P : Fin (N + 1) → Fin (N + 1) → ℝ) : Prop where
  nonneg : ∀ i j, 0 ≤ P i j
  rowSum : ∀ i, ∑ j, P i j = 1
  absorb0 : ∀ j, P 0 j = if j = 0 then 1 else 0
  absorbN : ∀ j, P (Fin.last N) j = if j = Fin.last N then 1 else 0

/-- Harmonic at the interior states (first-step equation for `0 < i < N`). -/
def InteriorHarmonic (P : Fin (N + 1) → Fin (N + 1) → ℝ) (g : Fin (N + 1) → ℝ) : Prop :=
  ∀ i : Fin (N + 1), 0 < (i : ℕ) → (i : ℕ) < N → ∑ j, P i j * g j = g i

theorem kPow_zero_apply (P : Fin (N + 1) → Fin (N + 1) → ℝ) (i j : Fin (N + 1)) :
    kPow P 0 i j = if i = j then 1 else 0 := rfl

theorem kPow_succ_apply (P : Fin (N + 1) → Fin (N + 1) → ℝ) (n : ℕ) (i j : Fin (N + 1)) :
    kPow P (n + 1) i j = ∑ k, kPow P n i k * P k j := rfl

theorem kPow_zero_sum (P : Fin (N + 1) → Fin (N + 1) → ℝ) (i : Fin (N + 1)) (f : Fin (N + 1) → ℝ) :
    ∑ j, kPow P 0 i j * f j = f i := by
  simp [kPow_zero_apply]

theorem kPow_succ_sum (P : Fin (N + 1) → Fin (N + 1) → ℝ) (n : ℕ) (i : Fin (N + 1))
    (f : Fin (N + 1) → ℝ) :
    ∑ j, kPow P (n + 1) i j * f j = ∑ k, kPow P n i k * ∑ j, P k j * f j := by
  simp only [kPow_succ_apply, Finset.sum_mul, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun k _ => Finset.sum_congr rfl (fun j _ => ?_))
  ring

theorem kPow_nonneg {P : Fin (N + 1) → Fin (N + 1) → ℝ} (hP : AbsorbingChain P) :
    ∀ n i j, 0 ≤ kPow P n i j := by
  intro n
  induction n with
  | zero =>
    intro i j
    rw [kPow_zero_apply]
    split_ifs <;> norm_num
  | succ n ih =>
    intro i j
    rw [kPow_succ_apply]
    exact Finset.sum_nonneg (fun k _ => mul_nonneg (ih i k) (hP.nonneg k j))

/-- `Pⁿ` preserves every function harmonic at all states (martingale identity). -/
theorem kPow_harmonic (P : Fin (N + 1) → Fin (N + 1) → ℝ) (g : Fin (N + 1) → ℝ)
    (hg : ∀ k, ∑ j, P k j * g j = g k) : ∀ n i, ∑ j, kPow P n i j * g j = g i := by
  intro n
  induction n with
  | zero => intro i; exact kPow_zero_sum P i g
  | succ n ih =>
    intro i
    rw [kPow_succ_sum]
    simp only [hg]
    exact ih i

theorem kPow_rowSum {P : Fin (N + 1) → Fin (N + 1) → ℝ} (hP : AbsorbingChain P) (n : ℕ)
    (i : Fin (N + 1)) : ∑ j, kPow P n i j = 1 := by
  have := kPow_harmonic P (fun _ => 1) (fun k => by simpa using hP.rowSum k) n i
  simpa using this

theorem fin_trichotomy (j : Fin (N + 1)) :
    j = 0 ∨ j = Fin.last N ∨ (0 < (j : ℕ) ∧ (j : ℕ) < N) := by
  by_cases h0 : (j : ℕ) = 0
  · exact Or.inl (Fin.ext h0)
  by_cases hN : (j : ℕ) = N
  · exact Or.inr (Or.inl (Fin.ext hN))
  · exact Or.inr (Or.inr ⟨Nat.pos_of_ne_zero h0, lt_of_le_of_ne (Nat.lt_succ_iff.mp j.isLt) hN⟩)

theorem zero_ne_last (hN : 0 < N) : (0 : Fin (N + 1)) ≠ Fin.last N := by
  intro h
  have := congrArg Fin.val h
  simp at this
  omega

theorem intInd_nonneg (j : Fin (N + 1)) : 0 ≤ intInd N j := by
  unfold intInd
  split_ifs <;> norm_num

theorem intInd_zero : intInd N 0 = 0 := by
  simp [intInd]

theorem intInd_last : intInd N (Fin.last N) = 0 := by
  simp [intInd]

theorem intInd_of_interior {j : Fin (N + 1)} (h1 : 0 < (j : ℕ)) (h2 : (j : ℕ) < N) :
    intInd N j = 1 := by
  simp [intInd, h1, h2]

theorem interiorMass_nonneg {P : Fin (N + 1) → Fin (N + 1) → ℝ} (hP : AbsorbingChain P) (n : ℕ)
    (i : Fin (N + 1)) : 0 ≤ interiorMass P n i :=
  Finset.sum_nonneg (fun k _ => mul_nonneg (kPow_nonneg hP n i k) (intInd_nonneg k))

/-- Interior harmonicity plus absorbing ends gives harmonicity at every state. -/
theorem harmonic_of_interior {P : Fin (N + 1) → Fin (N + 1) → ℝ} (hP : AbsorbingChain P)
    {g : Fin (N + 1) → ℝ} (hg : InteriorHarmonic P g) : ∀ k, ∑ j, P k j * g j = g k := by
  intro k
  rcases fin_trichotomy k with rfl | rfl | ⟨h1, h2⟩
  · simp [hP.absorb0]
  · simp [hP.absorbN]
  · exact hg k h1 h2

/-- Second-moment identity: `∑ⱼ P(k,j) h(j)² = h(k)² + ∑ⱼ P(k,j) (h j - h k)²` for harmonic `h`. -/
theorem second_moment {P : Fin (N + 1) → Fin (N + 1) → ℝ} (hP : AbsorbingChain P)
    {h : Fin (N + 1) → ℝ} (hall : ∀ k, ∑ j, P k j * h j = h k) (k : Fin (N + 1)) :
    ∑ j, P k j * h j ^ 2 = h k ^ 2 + ∑ j, P k j * (h j - h k) ^ 2 := by
  have e : ∀ j, P k j * (h j - h k) ^ 2
      = P k j * h j ^ 2 - 2 * h k * (P k j * h j) + h k ^ 2 * P k j := fun j => by ring
  rw [Finset.sum_congr rfl (fun j _ => e j), Finset.sum_add_distrib, Finset.sum_sub_distrib,
    ← Finset.mul_sum, ← Finset.mul_sum, hall k, hP.rowSum k]
  ring

/-- Interior mass tends to zero, given one interior-harmonic `h` with positive one-step variance at
every interior state. -/
theorem interiorMass_tendsto_zero {P : Fin (N + 1) → Fin (N + 1) → ℝ} (hP : AbsorbingChain P)
    {h : Fin (N + 1) → ℝ} (hh : InteriorHarmonic P h)
    (hvar : ∀ k : Fin (N + 1), 0 < (k : ℕ) → (k : ℕ) < N → 0 < ∑ j, P k j * (h j - h k) ^ 2)
    (i : Fin (N + 1)) : Tendsto (fun n => interiorMass P n i) atTop (𝓝 0) := by
  obtain ⟨ε, hε, hεle⟩ : ∃ ε > 0, ∀ k : Fin (N + 1), 0 < (k : ℕ) → (k : ℕ) < N →
      ε ≤ ∑ j, P k j * (h j - h k) ^ 2 := by
    by_cases hs : (Finset.univ.filter (fun k : Fin (N + 1) => 0 < (k : ℕ) ∧ (k : ℕ) < N)).Nonempty
    · obtain ⟨m, hm, hmin⟩ := Finset.exists_min_image _ (fun k => ∑ j, P k j * (h j - h k) ^ 2) hs
      rw [Finset.mem_filter] at hm
      exact ⟨_, hvar m hm.2.1 hm.2.2, fun k hk1 hk2 =>
        hmin k (Finset.mem_filter.mpr ⟨Finset.mem_univ _, hk1, hk2⟩)⟩
    · exact ⟨1, one_pos, fun k hk1 hk2 =>
        absurd ⟨k, Finset.mem_filter.mpr ⟨Finset.mem_univ _, hk1, hk2⟩⟩ hs⟩
  have hall := harmonic_of_interior hP hh
  have hsub : ∀ k, h k ^ 2 + ε * intInd N k ≤ ∑ j, P k j * h j ^ 2 := by
    intro k
    rcases fin_trichotomy k with rfl | rfl | ⟨h1, h2⟩
    · simp [hP.absorb0, intInd_zero]
    · simp [hP.absorbN, intInd_last]
    · rw [second_moment hP hall k, intInd_of_interior h1 h2, mul_one]
      linarith [hεle k h1 h2]
  have hacc : ∀ n, h i ^ 2 + ε * ∑ m ∈ range n, interiorMass P m i
      ≤ ∑ j, kPow P n i j * h j ^ 2 := by
    intro n
    induction n with
    | zero =>
      simp only [Finset.range_zero, Finset.sum_empty, mul_zero, add_zero]
      exact le_of_eq (kPow_zero_sum P i (fun j => h j ^ 2)).symm
    | succ n ih =>
      have e1 := kPow_succ_sum P n i (fun j => h j ^ 2)
      rw [e1,Finset.sum_range_succ, mul_add]
      have step : ∑ k, kPow P n i k * (h k ^ 2 + ε * intInd N k)
          ≤ ∑ k, kPow P n i k * ∑ j, P k j * h j ^ 2 :=
        Finset.sum_le_sum (fun k _ => mul_le_mul_of_nonneg_left (hsub k) (kPow_nonneg hP n i k))
      have e2 : ∑ k, kPow P n i k * (h k ^ 2 + ε * intInd N k)
          = ∑ k, kPow P n i k * h k ^ 2 + ε * interiorMass P n i := by
        unfold interiorMass
        rw [Finset.mul_sum, ← Finset.sum_add_distrib]
        exact Finset.sum_congr rfl (fun k _ => by ring)
      linarith
  obtain ⟨G, hG⟩ : ∃ G, ∀ j, h j ^ 2 ≤ G :=
    ⟨∑ j, h j ^ 2, fun j => Finset.single_le_sum (f := fun j => h j ^ 2)
      (fun j _ => sq_nonneg (h j)) (Finset.mem_univ j)⟩
  have hle : ∀ n, ∑ j, kPow P n i j * h j ^ 2 ≤ G := by
    intro n
    calc ∑ j, kPow P n i j * h j ^ 2 ≤ ∑ j, kPow P n i j * G :=
          Finset.sum_le_sum (fun j _ => mul_le_mul_of_nonneg_left (hG j) (kPow_nonneg hP n i j))
      _ = G := by rw [← Finset.sum_mul, kPow_rowSum hP, one_mul]
  have hpart : ∀ n, ∑ m ∈ range n, interiorMass P m i ≤ G / ε := by
    intro n
    rw [le_div_iff₀ hε]
    nlinarith [hacc n, hle n, sq_nonneg (h i)]
  have hsum : Summable (fun n => interiorMass P n i) :=
    summable_of_sum_range_le (fun n => interiorMass_nonneg hP n i) hpart
  exact hsum.tendsto_atTop_zero

/-- Squeeze: if interior mass vanishes, `Pⁿ(i,N)` converges to any interior-harmonic `g` with
boundary values `0` and `1`. -/
theorem absorption_tendsto {P : Fin (N + 1) → Fin (N + 1) → ℝ} (hP : AbsorbingChain P)
    (hN : 0 < N) {g : Fin (N + 1) → ℝ} (hg : InteriorHarmonic P g) (g0 : g 0 = 0)
    (g1 : g (Fin.last N) = 1)
    (hI : ∀ i, Tendsto (fun n => interiorMass P n i) atTop (𝓝 0)) (i : Fin (N + 1)) :
    Tendsto (fun n => kPow P n i (Fin.last N)) atTop (𝓝 (g i)) := by
  have hall := harmonic_of_interior hP hg
  obtain ⟨H, hH⟩ : ∃ H, ∀ j, |g j| ≤ H :=
    ⟨∑ j, |g j|, fun j => Finset.single_le_sum (f := fun j => |g j|)
      (fun j _ => abs_nonneg (g j)) (Finset.mem_univ j)⟩
  have h0l := zero_ne_last hN
  have hpt : ∀ j, (if j = Fin.last N then (1 : ℝ) else 0) - H * intInd N j ≤ g j ∧
      g j ≤ (if j = Fin.last N then (1 : ℝ) else 0) + H * intInd N j := by
    intro j
    rcases fin_trichotomy j with rfl | rfl | ⟨h1, h2⟩
    · simp [h0l, g0, intInd_zero]
    · simp [g1, intInd_last]
    · have hjl : j ≠ Fin.last N := fun e => by
        rw [e, Fin.val_last] at h2; exact lt_irrefl _ h2
      rw [if_neg hjl, intInd_of_interior h1 h2, mul_one]
      have := abs_le.mp (hH j)
      constructor <;> linarith [this.1, this.2]
  have hsplit : ∀ n (s : ℝ), ∑ j, kPow P n i j * ((if j = Fin.last N then (1 : ℝ) else 0)
      + s * intInd N j) = kPow P n i (Fin.last N) + s * interiorMass P n i := by
    intro n s
    have e : ∑ j, kPow P n i j * ((if j = Fin.last N then (1 : ℝ) else 0) + s * intInd N j)
        = ∑ j, kPow P n i j * (if j = Fin.last N then (1 : ℝ) else 0) + s * interiorMass P n i := by
      unfold interiorMass
      rw [Finset.mul_sum, ← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl (fun j _ => by ring)
    rw [e]
    simp [mul_ite]
  have hlow : ∀ n, g i - H * interiorMass P n i ≤ kPow P n i (Fin.last N) := by
    intro n
    have := Finset.sum_le_sum (fun j (_ : j ∈ (Finset.univ : Finset (Fin (N + 1)))) =>
      mul_le_mul_of_nonneg_left (hpt j).2 (kPow_nonneg hP n i j))
    rw [kPow_harmonic P g hall n i, hsplit n H] at this
    linarith
  have hupp : ∀ n, kPow P n i (Fin.last N) ≤ g i + H * interiorMass P n i := by
    intro n
    have := Finset.sum_le_sum (fun j (_ : j ∈ (Finset.univ : Finset (Fin (N + 1)))) =>
      mul_le_mul_of_nonneg_left (hpt j).1 (kPow_nonneg hP n i j))
    have e : ∀ j, (if j = Fin.last N then (1 : ℝ) else 0) - H * intInd N j
        = (if j = Fin.last N then (1 : ℝ) else 0) + (-H) * intInd N j := fun j => by ring
    simp only [e] at this
    rw [kPow_harmonic P g hall n i, hsplit n (-H)] at this
    linarith
  have t1 : Tendsto (fun n => g i - H * interiorMass P n i) atTop (𝓝 (g i)) := by
    have := (tendsto_const_nhds (x := g i)).sub ((hI i).const_mul H)
    simpa using this
  have t2 : Tendsto (fun n => g i + H * interiorMass P n i) atTop (𝓝 (g i)) := by
    have := (tendsto_const_nhds (x := g i)).add ((hI i).const_mul H)
    simpa using this
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le t1 t2 hlow hupp

/-- Conservation of mass: `Pⁿ(i,0) + Pⁿ(i,N) + Iₙ(i) = 1`. -/
theorem mass_split {P : Fin (N + 1) → Fin (N + 1) → ℝ} (hP : AbsorbingChain P) (hN : 0 < N)
    (n : ℕ) (i : Fin (N + 1)) :
    kPow P n i 0 + kPow P n i (Fin.last N) + interiorMass P n i = 1 := by
  have h0l := zero_ne_last hN
  have hpt : ∀ j : Fin (N + 1), (if j = 0 then (1 : ℝ) else 0) + (if j = Fin.last N then 1 else 0)
      + intInd N j = 1 := by
    intro j
    rcases fin_trichotomy j with rfl | rfl | ⟨h1, h2⟩
    · simp [h0l, intInd_zero]
    · simp [h0l.symm, intInd_last]
    · have hj0 : j ≠ 0 := fun e => by rw [e] at h1; simp at h1
      have hjl : j ≠ Fin.last N := fun e => by
        rw [e, Fin.val_last] at h2; exact lt_irrefl _ h2
      simp [hj0, hjl, intInd_of_interior h1 h2]
  calc kPow P n i 0 + kPow P n i (Fin.last N) + interiorMass P n i
      = ∑ j, kPow P n i j * ((if j = 0 then (1 : ℝ) else 0) + (if j = Fin.last N then 1 else 0)
          + intInd N j) := by
        unfold interiorMass
        simp only [mul_add, Finset.sum_add_distrib, mul_ite, mul_one, mul_zero,
          Finset.sum_ite_eq', Finset.mem_univ, if_true]
    _ = ∑ j, kPow P n i j := by simp only [hpt, mul_one]
    _ = 1 := kPow_rowSum hP n i

/-- `Pⁿ(i,N)` is non-decreasing in `n` (state `N` is absorbing), so its limit is the probability of
absorption at `N`. -/
theorem kPow_last_mono {P : Fin (N + 1) → Fin (N + 1) → ℝ} (hP : AbsorbingChain P) (n : ℕ)
    (i : Fin (N + 1)) : kPow P n i (Fin.last N) ≤ kPow P (n + 1) i (Fin.last N) := by
  rw [kPow_succ_apply]
  have := Finset.single_le_sum (f := fun k => kPow P n i k * P k (Fin.last N))
    (fun k _ => mul_nonneg (kPow_nonneg hP n i k) (hP.nonneg _ _)) (Finset.mem_univ (Fin.last N))
  simpa [hP.absorbN] using this

/-- Main general theorem: absorption-limit form of the fixation probability. -/
theorem fixation_limit {P : Fin (N + 1) → Fin (N + 1) → ℝ} (hP : AbsorbingChain P) (hN : 0 < N)
    {h : Fin (N + 1) → ℝ} (hh : InteriorHarmonic P h) (h0 : h 0 = 0) (h1 : h (Fin.last N) = 1)
    (hvar : ∀ k : Fin (N + 1), 0 < (k : ℕ) → (k : ℕ) < N → 0 < ∑ j, P k j * (h j - h k) ^ 2)
    (i : Fin (N + 1)) :
    Tendsto (fun n => kPow P n i (Fin.last N)) atTop (𝓝 (h i)) ∧
    Tendsto (fun n => kPow P n i 0) atTop (𝓝 (1 - h i)) := by
  have hI := interiorMass_tendsto_zero hP hh hvar
  have hlast := absorption_tendsto hP hN hh h0 h1 hI i
  refine ⟨hlast, ?_⟩
  have e : ∀ n, kPow P n i 0 = 1 - kPow P n i (Fin.last N) - interiorMass P n i := by
    intro n
    have := mass_split hP hN n i
    linarith
  have := ((tendsto_const_nhds (x := (1 : ℝ))).sub hlast).sub (hI i)
  rw [sub_zero] at this
  exact this.congr (fun n => (e n).symm)

/-- Uniqueness: any interior-harmonic `g` with boundary values `0`, `1` equals `h`. -/
theorem harmonic_unique_of_var {P : Fin (N + 1) → Fin (N + 1) → ℝ} (hP : AbsorbingChain P)
    (hN : 0 < N) {h : Fin (N + 1) → ℝ} (hh : InteriorHarmonic P h) (h0 : h 0 = 0)
    (h1 : h (Fin.last N) = 1)
    (hvar : ∀ k : Fin (N + 1), 0 < (k : ℕ) → (k : ℕ) < N → 0 < ∑ j, P k j * (h j - h k) ^ 2)
    {g : Fin (N + 1) → ℝ} (hg : InteriorHarmonic P g) (g0 : g 0 = 0) (g1 : g (Fin.last N) = 1) :
    ∀ i, g i = h i := by
  intro i
  have hI := interiorMass_tendsto_zero hP hh hvar
  exact tendsto_nhds_unique (absorption_tendsto hP hN hg g0 g1 hI i)
    (absorption_tendsto hP hN hh h0 h1 hI i)

/-! ## Negative control: the variance hypothesis is necessary

The identity kernel is stochastic with absorbing ends and every function is harmonic for it, yet
`Pⁿ(1,2) = 0` for all `n` (no absorption), so `Pⁿ(1,2) ↛ 1/2 = h(1)` for `h j = j/2`. -/

/-- The identity kernel. -/
noncomputable def idKernel (N : ℕ) (i j : Fin (N + 1)) : ℝ := if i = j then 1 else 0

theorem idKernel_absorbing : AbsorbingChain (idKernel N) where
  nonneg i j := by unfold idKernel; split_ifs <;> norm_num
  rowSum i := by simp [idKernel]
  absorb0 j := by
    by_cases h : j = 0
    · simp [idKernel, h]
    · simp [idKernel, h, Ne.symm h]
  absorbN j := by
    by_cases h : j = Fin.last N
    · simp [idKernel, h]
    · simp [idKernel, h, Ne.symm h]

theorem idKernel_harmonic (g : Fin (N + 1) → ℝ) : InteriorHarmonic (idKernel N) g := by
  intro i _ _
  simp [idKernel]

theorem kPow_idKernel (n : ℕ) (i j : Fin (N + 1)) : kPow (idKernel N) n i j = idKernel N i j := by
  induction n generalizing j with
  | zero => rfl
  | succ n ih =>
    rw [kPow_succ_apply]
    simp only [ih, idKernel, ite_mul, one_mul, zero_mul, Finset.sum_ite_eq, Finset.mem_univ,
      if_true]

theorem ctl_neg_identity_no_absorption :
    InteriorHarmonic (idKernel 2) (fun j => ((j : ℕ) : ℝ) / 2) ∧
    (∀ n, kPow (idKernel 2) n 1 (Fin.last 2) = 0) ∧
    ¬ Tendsto (fun n => kPow (idKernel 2) n 1 (Fin.last 2)) atTop (𝓝 (1 / 2)) := by
  have hz : ∀ n, kPow (idKernel 2) n 1 (Fin.last 2) = 0 := by
    intro n
    rw [kPow_idKernel]
    unfold idKernel
    rw [if_neg (by decide)]
  refine ⟨idKernel_harmonic _, hz, fun ht => ?_⟩
  have h0 : Tendsto (fun n => kPow (idKernel 2) n 1 (Fin.last 2)) atTop (𝓝 0) := by
    simp only [hz]; exact tendsto_const_nhds
  have := tendsto_nhds_unique ht h0
  norm_num at this

end ProgPopGen

#print axioms ProgPopGen.kPow_nonneg
#print axioms ProgPopGen.kPow_harmonic
#print axioms ProgPopGen.kPow_rowSum
#print axioms ProgPopGen.second_moment
#print axioms ProgPopGen.interiorMass_tendsto_zero
#print axioms ProgPopGen.absorption_tendsto
#print axioms ProgPopGen.mass_split
#print axioms ProgPopGen.kPow_last_mono
#print axioms ProgPopGen.fixation_limit
#print axioms ProgPopGen.harmonic_unique_of_var
#print axioms ProgPopGen.ctl_neg_identity_no_absorption

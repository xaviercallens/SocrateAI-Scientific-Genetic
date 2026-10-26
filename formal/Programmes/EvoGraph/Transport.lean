import Programmes.PopGen.Absorption

/-!
# Transport of P7's absorption limit to an arbitrary finite state space (N4, proposal 3, 2026-10-26)

P7 (`Programmes/PopGen/Absorption.lean`, kernel-checked, doi:10.5281/zenodo.23653850) proves the
absorption-limit theorem `ProgPopGen.fixation_limit` for stochastic kernels on `Fin (N+1)` with absorbing
ends `0` and `N`. The configuration chain of the Moran process on a graph lives on `Finset (Fin n)`, with
absorbing states `∅` and `univ`. This file does not re-prove the absorption argument: it builds an
equivalence `e : σ ≃ Fin (M+1)` sending the two absorbing states to `0` and `last M`, conjugates the
kernel, and applies P7's `fixation_limit` verbatim.

* `sPow`: the `t`-step kernel on `σ` (Chapman–Kolmogorov left recursion, as P7's `kPow`).
* `kPow_conj`: P7's `kPow` of the conjugated kernel is `sPow` read through `e`.
* `two_absorbing_limit`: for a stochastic kernel on a finite type with two absorbing states `bot ≠ top`,
  a function harmonic off `{bot, top}` with `h bot = 0`, `h top = 1` and positive one-step variance off
  `{bot, top}` is the limit of `Pᵗ(s, top)` (and `1 - h` the limit of `Pᵗ(s, bot)`).
-/

namespace EvoGraph

open Finset Filter Topology ProgPopGen

variable {σ : Type*} [Fintype σ] [DecidableEq σ]

/-- `t`-step kernel on a finite type, left recursion `P^(t+1)(s,u) = ∑ᵥ Pᵗ(s,v) P(v,u)`. -/
noncomputable def sPow (P : σ → σ → ℝ) : ℕ → σ → σ → ℝ
  | 0 => fun s u => if s = u then 1 else 0
  | t + 1 => fun s u => ∑ v, sPow P t s v * P v u

theorem sPow_succ_apply (P : σ → σ → ℝ) (t : ℕ) (s u : σ) :
    sPow P (t + 1) s u = ∑ v, sPow P t s v * P v u := rfl

/-- A stochastic kernel on `σ` with two distinct absorbing states. -/
structure TwoAbsorbing (P : σ → σ → ℝ) (bot top : σ) : Prop where
  ne : bot ≠ top
  nonneg : ∀ s u, 0 ≤ P s u
  rowSum : ∀ s, ∑ u, P s u = 1
  absorbBot : ∀ u, P bot u = if u = bot then 1 else 0
  absorbTop : ∀ u, P top u = if u = top then 1 else 0

omit [DecidableEq σ] in
/-- An enumeration of `σ` by `Fin (M+1)` with `bot ↦ 0` and `top ↦ last M`. -/
theorem exists_equiv_fin {bot top : σ} (hne : bot ≠ top) :
    ∃ M : ℕ, 0 < M ∧ ∃ e : σ ≃ Fin (M + 1), e bot = 0 ∧ e top = Fin.last M := by
  have hc : 1 < Fintype.card σ := Fintype.one_lt_card_iff.mpr ⟨bot, top, hne⟩
  refine ⟨Fintype.card σ - 1, by omega, ?_⟩
  let e0 : σ ≃ Fin (Fintype.card σ - 1 + 1) :=
    (Fintype.equivFin σ).trans (finCongr (by omega))
  let e1 : σ ≃ Fin (Fintype.card σ - 1 + 1) := e0.trans (Equiv.swap (e0 bot) 0)
  have h1 : e1 bot = 0 := by simp [e1, Equiv.swap_apply_left]
  have h1t : e1 top ≠ 0 := by
    rw [← h1]; exact fun h => hne (e1.injective h).symm
  let e2 : σ ≃ Fin (Fintype.card σ - 1 + 1) := e1.trans (Equiv.swap (e1 top) (Fin.last _))
  refine ⟨e2, ?_, ?_⟩
  · show Equiv.swap (e1 top) (Fin.last _) (e1 bot) = 0
    rw [h1]
    refine Equiv.swap_apply_of_ne_of_ne (Ne.symm h1t) ?_
    intro h
    have := congrArg Fin.val h
    simp at this
    omega
  · show Equiv.swap (e1 top) (Fin.last _) (e1 top) = Fin.last _
    exact Equiv.swap_apply_left _ _

/-- P7's `kPow` of the conjugated kernel equals `sPow` read through the enumeration. -/
theorem kPow_conj {M : ℕ} (e : σ ≃ Fin (M + 1)) (P : σ → σ → ℝ) (t : ℕ) (s u : σ) :
    kPow (fun i j => P (e.symm i) (e.symm j)) t (e s) (e u) = sPow P t s u := by
  induction t generalizing u with
  | zero =>
    show (if e s = e u then (1 : ℝ) else 0) = if s = u then 1 else 0
    simp only [e.apply_eq_iff_eq]
  | succ t ih =>
    rw [kPow_succ_apply, sPow_succ_apply]
    rw [← e.sum_comp (fun k => kPow (fun i j => P (e.symm i) (e.symm j)) t (e s) k *
      P (e.symm k) (e.symm (e u)))]
    simp only [Equiv.symm_apply_apply, ih]

/-- Absorption limit on a finite type with two absorbing states, obtained from P7's
`ProgPopGen.fixation_limit` by conjugation. -/
theorem two_absorbing_limit {P : σ → σ → ℝ} {bot top : σ} (hP : TwoAbsorbing P bot top)
    {h : σ → ℝ} (hh : ∀ s, s ≠ bot → s ≠ top → ∑ u, P s u * h u = h s) (h0 : h bot = 0)
    (h1 : h top = 1) (hvar : ∀ s, s ≠ bot → s ≠ top → 0 < ∑ u, P s u * (h u - h s) ^ 2)
    (s : σ) :
    Tendsto (fun t => sPow P t s top) atTop (𝓝 (h s)) ∧
    Tendsto (fun t => sPow P t s bot) atTop (𝓝 (1 - h s)) := by
  obtain ⟨M, hM, e, he0, heN⟩ := exists_equiv_fin hP.ne
  set Q : Fin (M + 1) → Fin (M + 1) → ℝ := fun i j => P (e.symm i) (e.symm j) with hQ
  have hb : e.symm 0 = bot := by rw [← he0, Equiv.symm_apply_apply]
  have ht : e.symm (Fin.last M) = top := by rw [← heN, Equiv.symm_apply_apply]
  have sumT : ∀ (f : σ → ℝ), ∑ j : Fin (M + 1), f (e.symm j) = ∑ u, f u :=
    fun f => e.symm.sum_comp f
  have hQabs : AbsorbingChain Q := {
    nonneg := fun i j => hP.nonneg _ _
    rowSum := fun i => (sumT (fun u => P (e.symm i) u)).trans (hP.rowSum _)
    absorb0 := fun j => by
      show P (e.symm 0) (e.symm j) = _
      rw [hb, hP.absorbBot]
      congr 1
      exact propext ⟨fun h => by rw [← he0, ← h, Equiv.apply_symm_apply],
        fun h => by rw [h, hb]⟩
    absorbN := fun j => by
      show P (e.symm (Fin.last M)) (e.symm j) = _
      rw [ht, hP.absorbTop]
      congr 1
      exact propext ⟨fun h => by rw [← heN, ← h, Equiv.apply_symm_apply],
        fun h => by rw [h, ht]⟩ }
  have hint : ∀ i : Fin (M + 1), 0 < (i : ℕ) → (i : ℕ) < M → e.symm i ≠ bot ∧ e.symm i ≠ top := by
    intro i h1 h2
    constructor
    · intro h
      have : i = 0 := by rw [← he0, ← h, Equiv.apply_symm_apply]
      rw [this] at h1; simp at h1
    · intro h
      have : i = Fin.last M := by rw [← heN, ← h, Equiv.apply_symm_apply]
      rw [this, Fin.val_last] at h2; exact lt_irrefl _ h2
  have hH : InteriorHarmonic Q (fun j => h (e.symm j)) := by
    intro i h1 h2
    obtain ⟨n0, n1⟩ := hint i h1 h2
    exact (sumT (fun u => P (e.symm i) u * h u)).trans (hh _ n0 n1)
  have hV : ∀ k : Fin (M + 1), 0 < (k : ℕ) → (k : ℕ) < M →
      0 < ∑ j, Q k j * (h (e.symm j) - h (e.symm k)) ^ 2 := by
    intro k h1 h2
    obtain ⟨n0, n1⟩ := hint k h1 h2
    rw [sumT (fun u => P (e.symm k) u * (h u - h (e.symm k)) ^ 2)]
    exact hvar _ n0 n1
  have main := fixation_limit hQabs hM hH (by rw [hb, h0]) (by rw [ht, h1]) hV (e s)
  rw [Equiv.symm_apply_apply] at main
  refine ⟨main.1.congr (fun t => ?_), main.2.congr (fun t => ?_)⟩
  · rw [← heN]; exact kPow_conj e P t s top
  · rw [← he0]; exact kPow_conj e P t s bot

theorem sPow_nonneg {P : σ → σ → ℝ} {bot top : σ} (hP : TwoAbsorbing P bot top) :
    ∀ t s u, 0 ≤ sPow P t s u := by
  intro t
  induction t with
  | zero => intro s u; show 0 ≤ (if s = u then (1 : ℝ) else 0); split_ifs <;> norm_num
  | succ t ih =>
    intro s u
    rw [sPow_succ_apply]
    exact Finset.sum_nonneg (fun v _ => mul_nonneg (ih s v) (hP.nonneg v u))

/-- `Pᵗ(s, top)` is non-decreasing in `t` (`top` is absorbing), so its limit is the probability of
ever being absorbed in `top`. -/
theorem sPow_top_mono {P : σ → σ → ℝ} {bot top : σ} (hP : TwoAbsorbing P bot top) (t : ℕ) (s : σ) :
    sPow P t s top ≤ sPow P (t + 1) s top := by
  rw [sPow_succ_apply]
  have := Finset.single_le_sum (f := fun v => sPow P t s v * P v top)
    (fun v _ => mul_nonneg (sPow_nonneg hP t s v) (hP.nonneg _ _)) (Finset.mem_univ top)
  simpa [hP.absorbTop] using this

end EvoGraph

#print axioms EvoGraph.sPow_top_mono
#print axioms EvoGraph.exists_equiv_fin
#print axioms EvoGraph.kPow_conj
#print axioms EvoGraph.two_absorbing_limit

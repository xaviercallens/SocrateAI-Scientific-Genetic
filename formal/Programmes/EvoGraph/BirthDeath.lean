import Programmes.EvoGraph.Transport
import Programmes.PopGen.MoranModel

/-!
# Birth–death Moran process on a weighted digraph; the isothermal theorem (N4, proposal 3)

Population on the nodes `Fin n` of a weighted directed graph with row-stochastic weight matrix `W`
(`W i j ≥ 0`, `∑ⱼ W i j = 1`). A configuration is the set `S : Finset (Fin n)` of mutant nodes; mutants
have fitness `r > 0`, residents fitness `1`. Selection weights `μ i > 0` (`μ ≡ 1` is the
Lieberman–Hauert–Nowak 2005 birth–death process). One step:

* node `i` is chosen to reproduce with probability `μ i f_S(i) / F(S)`, `F(S) = ∑ₖ μ k f_S(k)`;
* its offspring replaces node `j` with probability `W i j` (`j = i` allowed);
* the new configuration is `upd S i j` (`insert j S` if `i ∈ S`, `S.erase j` otherwise).

The configuration kernel is `bdK μ W r S T = ∑ᵢ ∑ⱼ coef(S,i,j) [upd S i j = T]`.

Results (all for the configuration chain; the mutant count is NOT claimed to be a Markov chain, see
`Examples.lean` for a kernel-checked isothermal counterexample to lumpability):

* `bd_absorbing`: `bdK` is stochastic with absorbing `∅` and `univ` (row-stochastic `W` only).
* `up_eq`, `down_eq`: `P⁺(S) = r·cutOut(S)/F(S)`, `P⁻(S) = cutIn(S)/F(S)`.
* `cut_balance` (stationary `μ`: `∑ᵢ μ i W i j = μ j`): `cutOut(S) = cutIn(S)` for every `S`.
* `forward_bias`: `P⁺(S) = r · P⁻(S)` for every configuration `S`.
* `moranFix_card_harmonic`: `S ↦ moranFix n r |S|` is harmonic for the configuration chain.
* `stationary_fixation`: `Pᵗ(S, univ) → moranFix n r |S|` (Keller–Uğurlu 2024, Thm 3.1, stationary
  selection), under cut positivity (equivalent to strong connectivity; `cutPositive_of_strongly_connected`
  proves the direction used).
* `isothermal_fixation` (μ ≡ 1, `W` doubly stochastic, strongly connected): the isothermal theorem of
  Lieberman–Hauert–Nowak (2005): `Pᵗ(S, univ) → moranFix n r |S|`, for a single mutant at any node
  `(1 - r⁻¹)/(1 - r⁻ⁿ)` (`r ≠ 1`) or `1/n` (`r = 1`).

Known results (LHN 2005; Pattni–Broom–Rychtář–Silvers 2015; Keller–Uğurlu 2024); novelty is not claimed.
The proof is the configuration-level harmonic function / martingale route (no lumpability).
-/

namespace EvoGraph

open Finset Filter Topology ProgPopGen

variable {n : ℕ}

/-- Fitness of node `i` in configuration `S`. -/
noncomputable def fit (r : ℝ) (S : Finset (Fin n)) (i : Fin n) : ℝ := if i ∈ S then r else 1

/-- Indicator of `i ∈ S`. -/
noncomputable def ind (S : Finset (Fin n)) (i : Fin n) : ℝ := if i ∈ S then 1 else 0

/-- Offspring of `i` replaces `j`. -/
def upd (S : Finset (Fin n)) (i j : Fin n) : Finset (Fin n) := if i ∈ S then insert j S else S.erase j

/-- Total selection weight `F(S) = ∑ₖ μ k f_S(k)`. -/
noncomputable def totFit (μ : Fin n → ℝ) (r : ℝ) (S : Finset (Fin n)) : ℝ := ∑ k, μ k * fit r S k

/-- Probability that `i` reproduces and its offspring replaces `j`. -/
noncomputable def coef (μ : Fin n → ℝ) (W : Fin n → Fin n → ℝ) (r : ℝ) (S : Finset (Fin n))
    (i j : Fin n) : ℝ :=
  μ i * fit r S i / totFit μ r S * W i j

/-- Configuration kernel of the birth–death process. -/
noncomputable def bdK (μ : Fin n → ℝ) (W : Fin n → Fin n → ℝ) (r : ℝ) (S T : Finset (Fin n)) : ℝ :=
  ∑ i, ∑ j, coef μ W r S i j * (if upd S i j = T then 1 else 0)

/-- Weight flowing out of `S`. -/
noncomputable def cutOut (μ : Fin n → ℝ) (W : Fin n → Fin n → ℝ) (S : Finset (Fin n)) : ℝ :=
  ∑ i, ∑ j, μ i * W i j * (ind S i * (1 - ind S j))

/-- Weight flowing into `S`. -/
noncomputable def cutIn (μ : Fin n → ℝ) (W : Fin n → Fin n → ℝ) (S : Finset (Fin n)) : ℝ :=
  ∑ i, ∑ j, μ i * W i j * ((1 - ind S i) * ind S j)

/-- Probability that the number of mutants increases by one. -/
noncomputable def upP (μ : Fin n → ℝ) (W : Fin n → Fin n → ℝ) (r : ℝ) (S : Finset (Fin n)) : ℝ :=
  ∑ T, bdK μ W r S T * (if T.card = S.card + 1 then 1 else 0)

/-- Probability that the number of mutants decreases by one. -/
noncomputable def downP (μ : Fin n → ℝ) (W : Fin n → Fin n → ℝ) (r : ℝ) (S : Finset (Fin n)) : ℝ :=
  ∑ T, bdK μ W r S T * (if T.card + 1 = S.card then 1 else 0)

/-- Row-stochastic weights. -/
structure RowStochastic (W : Fin n → Fin n → ℝ) : Prop where
  nonneg : ∀ i j, 0 ≤ W i j
  rowSum : ∀ i, ∑ j, W i j = 1

/-- Cut positivity: every proper non-empty `S` has an edge leaving it. -/
def CutPositive (μ : Fin n → ℝ) (W : Fin n → Fin n → ℝ) : Prop :=
  ∀ S : Finset (Fin n), S ≠ ∅ → S ≠ univ → ∃ i ∈ S, ∃ j ∉ S, 0 < μ i * W i j

/-- Strong connectivity of the digraph `{(a,b) : W a b > 0}`. -/
def StronglyConnected (W : Fin n → Fin n → ℝ) : Prop :=
  ∀ a b, Relation.TransGen (fun x y => 0 < W x y) a b

/-! ## Basic facts -/

theorem fit_pos {r : ℝ} (hr : 0 < r) (S : Finset (Fin n)) (i : Fin n) : 0 < fit r S i := by
  unfold fit; split_ifs <;> linarith

theorem fit_eq (r : ℝ) (S : Finset (Fin n)) (i : Fin n) : fit r S i = 1 + (r - 1) * ind S i := by
  unfold fit ind; split_ifs <;> ring

theorem totFit_pos {μ : Fin n → ℝ} {r : ℝ} (hn : 0 < n) (hr : 0 < r) (hμ : ∀ i, 0 < μ i)
    (S : Finset (Fin n)) : 0 < totFit μ r S := by
  unfold totFit
  have : Nonempty (Fin n) := ⟨⟨0, hn⟩⟩
  exact Finset.sum_pos (fun i _ => mul_pos (hμ i) (fit_pos hr S i)) Finset.univ_nonempty

theorem coef_nonneg {μ : Fin n → ℝ} {W : Fin n → Fin n → ℝ} {r : ℝ} (hn : 0 < n) (hr : 0 < r)
    (hμ : ∀ i, 0 < μ i) (hW : RowStochastic W) (S : Finset (Fin n)) (i j : Fin n) :
    0 ≤ coef μ W r S i j := by
  unfold coef
  exact mul_nonneg (div_nonneg (mul_pos (hμ i) (fit_pos hr S i)).le (totFit_pos hn hr hμ S).le)
    (hW.nonneg i j)

theorem coef_sum {μ : Fin n → ℝ} {W : Fin n → Fin n → ℝ} {r : ℝ} (hn : 0 < n) (hr : 0 < r)
    (hμ : ∀ i, 0 < μ i) (hW : RowStochastic W) (S : Finset (Fin n)) :
    ∑ i, ∑ j, coef μ W r S i j = 1 := by
  have hF := (totFit_pos hn hr hμ S).ne'
  unfold coef
  simp only [← Finset.mul_sum, hW.rowSum, mul_one]
  rw [← Finset.sum_div]
  exact div_self hF

/-- First-step expectation: `∑_T P(S,T) g(T) = ∑ᵢ ∑ⱼ coef(S,i,j) g(upd S i j)`. -/
theorem expect_eq (μ : Fin n → ℝ) (W : Fin n → Fin n → ℝ) (r : ℝ) (S : Finset (Fin n))
    (g : Finset (Fin n) → ℝ) :
    ∑ T, bdK μ W r S T * g T = ∑ i, ∑ j, coef μ W r S i j * g (upd S i j) := by
  unfold bdK
  simp only [Finset.sum_mul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  simp only [mul_assoc, ite_mul, one_mul, zero_mul]
  rw [← Finset.mul_sum, Finset.sum_ite_eq]
  simp

theorem upd_empty (i j : Fin n) : upd ∅ i j = ∅ := by simp [upd]

theorem upd_univ (i j : Fin n) : upd univ i j = univ := by simp [upd]

theorem empty_ne_univ (hn : 0 < n) : (∅ : Finset (Fin n)) ≠ univ := by
  intro h
  have : (⟨0, hn⟩ : Fin n) ∈ (univ : Finset (Fin n)) := Finset.mem_univ _
  rw [← h] at this
  simp at this

/-- The configuration chain is stochastic with absorbing states `∅` and `univ`. -/
theorem bd_absorbing {μ : Fin n → ℝ} {W : Fin n → Fin n → ℝ} {r : ℝ} (hn : 0 < n) (hr : 0 < r)
    (hμ : ∀ i, 0 < μ i) (hW : RowStochastic W) :
    TwoAbsorbing (bdK μ W r) ∅ univ where
  ne := empty_ne_univ hn
  nonneg S T := by
    unfold bdK
    exact Finset.sum_nonneg (fun i _ => Finset.sum_nonneg (fun j _ =>
      mul_nonneg (coef_nonneg hn hr hμ hW S i j) (by split_ifs <;> norm_num)))
  rowSum S := by
    have := expect_eq μ W r S (fun _ => 1)
    simp only [mul_one] at this
    rw [this, coef_sum hn hr hμ hW]
  absorbBot T := by
    unfold bdK
    simp only [upd_empty]
    rw [← Finset.sum_congr rfl (fun i _ => Finset.sum_mul .. ), ← Finset.sum_mul,
      coef_sum hn hr hμ hW, one_mul]
    by_cases h : T = ∅ <;> simp [h, eq_comm]
  absorbTop T := by
    unfold bdK
    simp only [upd_univ]
    rw [← Finset.sum_congr rfl (fun i _ => Finset.sum_mul .. ), ← Finset.sum_mul,
      coef_sum hn hr hμ hW, one_mul]
    by_cases h : T = univ <;> simp [h, eq_comm]

/-! ## Cut balance and the forward bias -/

/-- For stationary `μ` (`∑ᵢ μ i W i j = μ j`) and row-stochastic `W`, the net flow out of any set is `0`. -/
theorem flow_balance {μ : Fin n → ℝ} {W : Fin n → Fin n → ℝ} (hrow : ∀ i, ∑ j, W i j = 1)
    (hstat : ∀ j, ∑ i, μ i * W i j = μ j) (S : Finset (Fin n)) :
    ∑ i, ∑ j, μ i * W i j * (ind S i - ind S j) = 0 := by
  have e : ∀ i j, μ i * W i j * (ind S i - ind S j) = μ i * ind S i * W i j - μ i * W i j * ind S j :=
    fun i j => by ring
  simp only [e, Finset.sum_sub_distrib, ← Finset.mul_sum, hrow, mul_one]
  rw [Finset.sum_comm]
  simp only [← Finset.sum_mul, hstat]
  rw [sub_eq_zero]

theorem cut_balance {μ : Fin n → ℝ} {W : Fin n → Fin n → ℝ} (hrow : ∀ i, ∑ j, W i j = 1)
    (hstat : ∀ j, ∑ i, μ i * W i j = μ j) (S : Finset (Fin n)) : cutOut μ W S = cutIn μ W S := by
  have h := flow_balance hrow hstat S
  unfold cutOut cutIn
  rw [← sub_eq_zero, ← Finset.sum_sub_distrib]
  rw [← h]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl (fun j _ => by ring)

theorem card_upd_in_out {S : Finset (Fin n)} {i j : Fin n} (hi : i ∈ S) (hj : j ∉ S) :
    (upd S i j).card = S.card + 1 := by
  simp [upd, hi, Finset.card_insert_of_notMem hj]

theorem upd_in_in {S : Finset (Fin n)} {i j : Fin n} (hi : i ∈ S) (hj : j ∈ S) : upd S i j = S := by
  simp [upd, hi, Finset.insert_eq_of_mem hj]

theorem card_upd_out_in {S : Finset (Fin n)} {i j : Fin n} (hi : i ∉ S) (hj : j ∈ S) :
    (upd S i j).card + 1 = S.card := by
  simp only [upd, hi, if_false]
  rw [Finset.card_erase_of_mem hj]
  have : 0 < S.card := Finset.card_pos.mpr ⟨j, hj⟩
  omega

theorem upd_out_out {S : Finset (Fin n)} {i j : Fin n} (hi : i ∉ S) (hj : j ∉ S) : upd S i j = S := by
  simp [upd, hi, Finset.erase_eq_of_notMem hj]

/-- `P⁺(S) = r · cutOut(S) / F(S)`. -/
theorem up_eq (μ : Fin n → ℝ) (W : Fin n → Fin n → ℝ) (r : ℝ) (S : Finset (Fin n)) :
    upP μ W r S = r * cutOut μ W S / totFit μ r S := by
  unfold upP
  rw [expect_eq]
  unfold cutOut
  rw [Finset.mul_sum, Finset.sum_div]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [Finset.mul_sum, Finset.sum_div]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  unfold coef fit ind
  by_cases hi : i ∈ S <;> by_cases hj : j ∈ S
  · rw [upd_in_in hi hj]; simp [hi, hj]
  · rw [card_upd_in_out hi hj]; simp [hi, hj]; ring
  · have := card_upd_out_in hi hj
    have h2 : (upd S i j).card ≠ S.card + 1 := by omega
    simp [hi, hj, h2]
  · rw [upd_out_out hi hj]; simp [hi, hj]

/-- `P⁻(S) = cutIn(S) / F(S)`. -/
theorem down_eq (μ : Fin n → ℝ) (W : Fin n → Fin n → ℝ) (r : ℝ) (S : Finset (Fin n)) :
    downP μ W r S = cutIn μ W S / totFit μ r S := by
  unfold downP
  rw [expect_eq]
  unfold cutIn
  rw [Finset.sum_div]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [Finset.sum_div]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  unfold coef fit ind
  by_cases hi : i ∈ S <;> by_cases hj : j ∈ S
  · rw [upd_in_in hi hj]; simp [hi, hj]
  · have h2 : (upd S i j).card + 1 ≠ S.card := by rw [card_upd_in_out hi hj]; omega
    simp [hi, hj, h2]
  · simp [hi, hj, card_upd_out_in hi hj]; ring
  · rw [upd_out_out hi hj]; simp [hi, hj]

/-- Forward bias: for stationary `μ`, `P⁺(S) = r · P⁻(S)` in every configuration. -/
theorem forward_bias {μ : Fin n → ℝ} {W : Fin n → Fin n → ℝ} (r : ℝ) (hrow : ∀ i, ∑ j, W i j = 1)
    (hstat : ∀ j, ∑ i, μ i * W i j = μ j) (S : Finset (Fin n)) :
    upP μ W r S = r * downP μ W r S := by
  rw [up_eq, down_eq, cut_balance hrow hstat S]
  ring

/-! ## The Moran fixation function is harmonic on configurations -/

/-- Increment of the Moran fixation function along one transition. -/
theorem coef_increment {μ : Fin n → ℝ} {W : Fin n → Fin n → ℝ} {r : ℝ} (hr : 0 < r)
    {S : Finset (Fin n)} (hS : 0 < S.card) (i j : Fin n) :
    coef μ W r S i j * (moranFix n r (upd S i j).card - moranFix n r S.card)
      = μ i * W i j / totFit μ r S * (r⁻¹ ^ (S.card - 1) / geomS r n) * (ind S i - ind S j) := by
  obtain ⟨m, hm⟩ : ∃ m, S.card = m + 1 := ⟨S.card - 1, by omega⟩
  have hrr : r * r⁻¹ = 1 := mul_inv_cancel₀ hr.ne'
  unfold coef fit ind
  by_cases hi : i ∈ S <;> by_cases hj : j ∈ S
  · rw [upd_in_in hi hj]; simp [hi, hj]
  · rw [card_upd_in_out hi hj, moranFix_step, hm]
    simp only [hi, hj, if_true, if_false, Nat.add_sub_cancel, pow_succ]
    have : r * (r⁻¹ ^ m * r⁻¹) = r⁻¹ ^ m := by rw [mul_comm (r⁻¹ ^ m), ← mul_assoc, hrr, one_mul]
    linear_combination (μ i * W i j / totFit μ r S / geomS r n * r⁻¹ ^ m) * hrr
  · have hc := card_upd_out_in hi hj
    rw [hm] at hc
    have hc' : (upd S i j).card = m := by omega
    rw [hc', hm]
    have := moranFix_step n r m
    simp only [hi, hj, if_true, if_false, Nat.add_sub_cancel]
    rw [show moranFix n r m - moranFix n r (m + 1) = -(r⁻¹ ^ m / geomS r n) by linarith]
    ring
  · rw [upd_out_out hi hj]; simp [hi, hj]

/-- `S ↦ moranFix n r |S|` satisfies the first-step equation of the configuration chain at every
non-empty `S` (stationary `μ`, row-stochastic `W`). -/
theorem moranFix_card_harmonic {μ : Fin n → ℝ} {W : Fin n → Fin n → ℝ} {r : ℝ} (hn : 0 < n)
    (hr : 0 < r) (hμ : ∀ i, 0 < μ i) (hW : RowStochastic W) (hstat : ∀ j, ∑ i, μ i * W i j = μ j)
    {S : Finset (Fin n)} (hS : 0 < S.card) :
    ∑ T, bdK μ W r S T * moranFix n r T.card = moranFix n r S.card := by
  rw [expect_eq]
  have hsplit : ∀ i j, coef μ W r S i j * moranFix n r (upd S i j).card
      = coef μ W r S i j * (moranFix n r (upd S i j).card - moranFix n r S.card)
        + coef μ W r S i j * moranFix n r S.card := fun i j => by ring
  simp only [hsplit, Finset.sum_add_distrib, ← Finset.sum_mul, coef_sum hn hr hμ hW, one_mul,
    coef_increment hr hS]
  have hfb := flow_balance hW.rowSum hstat S
  have : ∑ i, ∑ j, μ i * W i j / totFit μ r S * (r⁻¹ ^ (S.card - 1) / geomS r n)
      * (ind S i - ind S j)
      = (r⁻¹ ^ (S.card - 1) / geomS r n) / totFit μ r S
        * ∑ i, ∑ j, μ i * W i j * (ind S i - ind S j) := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl (fun j _ => ?_)
    ring
  rw [this, hfb, mul_zero, zero_add]

/-- Positive one-step variance of `moranFix n r |·|` at every proper non-empty `S`, given an edge
leaving `S`. -/
theorem moranFix_card_var_pos {μ : Fin n → ℝ} {W : Fin n → Fin n → ℝ} {r : ℝ} (hn : 0 < n)
    (hr : 0 < r) (hμ : ∀ i, 0 < μ i) (hW : RowStochastic W) {S : Finset (Fin n)}
    {i0 j0 : Fin n} (hi0 : i0 ∈ S) (hj0 : j0 ∉ S) (hw : 0 < μ i0 * W i0 j0) :
    0 < ∑ T, bdK μ W r S T * (moranFix n r T.card - moranFix n r S.card) ^ 2 := by
  rw [expect_eq μ W r S (fun T => (moranFix n r T.card - moranFix n r S.card) ^ 2)]
  have hterm : ∀ i j, 0 ≤ coef μ W r S i j * (moranFix n r (upd S i j).card - moranFix n r S.card) ^ 2 :=
    fun i j => mul_nonneg (coef_nonneg hn hr hμ hW S i j) (sq_nonneg _)
  refine Finset.sum_pos' (fun i _ => Finset.sum_nonneg (fun j _ => hterm i j))
    ⟨i0, Finset.mem_univ _, Finset.sum_pos' (fun j _ => hterm i0 j) ⟨j0, Finset.mem_univ _, ?_⟩⟩
  rw [card_upd_in_out hi0 hj0, moranFix_step]
  have hc : 0 < coef μ W r S i0 j0 := by
    unfold coef
    have hf : fit r S i0 = r := by simp [fit, hi0]
    rw [hf]
    have := totFit_pos hn hr hμ S
    have h' : μ i0 * r / totFit μ r S * W i0 j0 = r / totFit μ r S * (μ i0 * W i0 j0) := by ring
    rw [h']
    exact mul_pos (div_pos hr this) hw
  have hS : 0 < geomS r n := geomS_pos hr hn
  exact mul_pos hc (pow_pos (div_pos (pow_pos (inv_pos.mpr hr) _) hS) 2)

/-! ## Fixation theorems -/

theorem cutPositive_of_strongly_connected {μ : Fin n → ℝ} {W : Fin n → Fin n → ℝ}
    (hμ : ∀ i, 0 < μ i) (hsc : StronglyConnected W) : CutPositive μ W := by
  intro S hS0 hSu
  obtain ⟨a, ha⟩ := Finset.nonempty_iff_ne_empty.mpr hS0
  obtain ⟨b, hb⟩ : ∃ b, b ∉ S := by
    by_contra h
    push Not at h
    exact hSu (Finset.eq_univ_iff_forall.mpr h)
  have key : ∀ x y, Relation.TransGen (fun x y => 0 < W x y) x y → x ∈ S → y ∉ S →
      ∃ i ∈ S, ∃ j ∉ S, 0 < W i j := by
    intro x y hxy
    induction hxy with
    | single h => exact fun hx hy => ⟨_, hx, _, hy, h⟩
    | tail _ hcb ih =>
      rename_i c _ _
      intro hx hy
      by_cases hc : c ∈ S
      · exact ⟨c, hc, _, hy, hcb⟩
      · exact ih hx hc
  obtain ⟨i, hi, j, hj, hw⟩ := key a b (hsc a b) ha hb
  exact ⟨i, hi, j, hj, mul_pos (hμ i) hw⟩

/-- Fixation under stationary selection: the configuration chain started from `S` is absorbed in
`univ` with limiting probability `moranFix n r |S|` (and in `∅` with `1 - moranFix n r |S|`). -/
theorem stationary_fixation {μ : Fin n → ℝ} {W : Fin n → Fin n → ℝ} {r : ℝ} (hn : 0 < n)
    (hr : 0 < r) (hμ : ∀ i, 0 < μ i) (hW : RowStochastic W) (hstat : ∀ j, ∑ i, μ i * W i j = μ j)
    (hcut : CutPositive μ W) (S : Finset (Fin n)) :
    Tendsto (fun t => sPow (bdK μ W r) t S univ) atTop (𝓝 (moranFix n r S.card)) ∧
    Tendsto (fun t => sPow (bdK μ W r) t S ∅) atTop (𝓝 (1 - moranFix n r S.card)) := by
  refine two_absorbing_limit (bd_absorbing hn hr hμ hW) (h := fun T => moranFix n r T.card)
    (fun T h0 _ => moranFix_card_harmonic hn hr hμ hW hstat
      (Finset.card_pos.mpr (Finset.nonempty_iff_ne_empty.mpr h0)))
    (by simp [moranFix_zero]) (by simp [Finset.card_univ, moranFix_last hn hr]) ?_ S
  intro T h0 h1
  obtain ⟨i, hi, j, hj, hw⟩ := hcut T h0 h1
  exact moranFix_card_var_pos hn hr hμ hW hi hj hw

/-- Doubly stochastic (isothermal) weights: rows and columns sum to one. -/
structure Isothermal (W : Fin n → Fin n → ℝ) : Prop extends RowStochastic W where
  colSum : ∀ j, ∑ i, W i j = 1

/-- **Isothermal theorem** (Lieberman–Hauert–Nowak 2005, birth–death updating with uniform
selection `μ ≡ 1`): on a strongly connected isothermal graph the configuration chain started from
mutant set `S` fixes with limiting probability `moranFix n r |S|`. -/
theorem isothermal_fixation {W : Fin n → Fin n → ℝ} {r : ℝ} (hn : 0 < n) (hr : 0 < r)
    (hW : Isothermal W) (hsc : StronglyConnected W) (S : Finset (Fin n)) :
    Tendsto (fun t => sPow (bdK (fun _ => 1) W r) t S univ) atTop (𝓝 (moranFix n r S.card)) := by
  have hμ : ∀ _ : Fin n, (0 : ℝ) < 1 := fun _ => one_pos
  exact (stationary_fixation hn hr hμ hW.toRowStochastic (fun j => by simpa using hW.colSum j)
    (cutPositive_of_strongly_connected hμ hsc) S).1

/-- The fixation sequence `Pᵗ(S, univ)` of the configuration chain is non-decreasing in `t`, so the
limit in `isothermal_fixation` is the probability of ever reaching `univ`. -/
theorem bd_fixation_mono {μ : Fin n → ℝ} {W : Fin n → Fin n → ℝ} {r : ℝ} (hn : 0 < n) (hr : 0 < r)
    (hμ : ∀ i, 0 < μ i) (hW : RowStochastic W) (t : ℕ) (S : Finset (Fin n)) :
    sPow (bdK μ W r) t S univ ≤ sPow (bdK μ W r) (t + 1) S univ :=
  sPow_top_mono (bd_absorbing hn hr hμ hW) t S

/-- Isothermal theorem, single mutant at any node `v`, `r ≠ 1`: limit `(1 - r⁻¹)/(1 - r⁻ⁿ)`. -/
theorem isothermal_fixation_single {W : Fin n → Fin n → ℝ} {r : ℝ} (hn : 0 < n) (hr : 0 < r)
    (hr1 : r ≠ 1) (hW : Isothermal W) (hsc : StronglyConnected W) (v : Fin n) :
    Tendsto (fun t => sPow (bdK (fun _ => 1) W r) t {v} univ) atTop
      (𝓝 ((1 - r⁻¹) / (1 - r⁻¹ ^ n))) := by
  have := isothermal_fixation hn hr hW hsc {v}
  rw [Finset.card_singleton, moranFix_closed hr1, pow_one] at this
  exact this

/-- Isothermal theorem, neutral case `r = 1`: limit `|S|/n`. -/
theorem isothermal_fixation_neutral {W : Fin n → Fin n → ℝ} (hn : 0 < n) (hW : Isothermal W)
    (hsc : StronglyConnected W) (S : Finset (Fin n)) :
    Tendsto (fun t => sPow (bdK (fun _ => 1) W 1) t S univ) atTop (𝓝 ((S.card : ℝ) / n)) := by
  have := isothermal_fixation hn one_pos hW hsc S
  rwa [moranFix_neutral] at this

/-- Isothermal theorem, single mutant placed uniformly at random, `r ≠ 1`. -/
theorem isothermal_fixation_uniform {W : Fin n → Fin n → ℝ} {r : ℝ} (hn : 0 < n) (hr : 0 < r)
    (hr1 : r ≠ 1) (hW : Isothermal W) (hsc : StronglyConnected W) :
    Tendsto (fun t => (∑ v : Fin n, sPow (bdK (fun _ => 1) W r) t {v} univ) / n) atTop
      (𝓝 ((1 - r⁻¹) / (1 - r⁻¹ ^ n))) := by
  have hnr : (n : ℝ) ≠ 0 := by exact_mod_cast hn.ne'
  have hs := tendsto_finsetSum (Finset.univ : Finset (Fin n))
    (fun v _ => isothermal_fixation_single hn hr hr1 hW hsc v)
  have := hs.div_const (n : ℝ)
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at this
  rwa [mul_div_cancel_left₀ _ hnr] at this

end EvoGraph

#print axioms EvoGraph.expect_eq
#print axioms EvoGraph.bd_absorbing
#print axioms EvoGraph.flow_balance
#print axioms EvoGraph.cut_balance
#print axioms EvoGraph.up_eq
#print axioms EvoGraph.down_eq
#print axioms EvoGraph.forward_bias
#print axioms EvoGraph.moranFix_card_harmonic
#print axioms EvoGraph.moranFix_card_var_pos
#print axioms EvoGraph.cutPositive_of_strongly_connected
#print axioms EvoGraph.stationary_fixation
#print axioms EvoGraph.isothermal_fixation
#print axioms EvoGraph.bd_fixation_mono
#print axioms EvoGraph.isothermal_fixation_single
#print axioms EvoGraph.isothermal_fixation_neutral
#print axioms EvoGraph.isothermal_fixation_uniform

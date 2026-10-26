import Programmes.EvoGraph.BirthDeath

/-!
# Positive and negative controls for the isothermal theorem (N4, proposal 3)

Real theorems (planted `sorry`/axiom controls live in `EvoGraphPlantedControls.lean`). All with uniform
selection `μ ≡ 1` (Lieberman–Hauert–Nowak birth–death) unless stated otherwise, fitness `r = 2`.

Positive controls
* `wellMixed_up`, `wellMixed_down`: for `W_H = J/n` (complete graph with loops) the configuration
  chain lumps exactly onto P7's Moran chain: `P⁺(S) = moranUp n r |S|`, `P⁻(S) = moranDown n r |S|`.
* `wellMixed3_fixation`: the isothermal theorem on `W_H`, `n = 3`, gives the P7 value `4/7`.
* `starStationary_fixation`: the star `n = 3` with stationary selection `μ = (1/2, 1/4, 1/4)` fixes
  from the centre with the Moran value `4/7` (stationary-selection theorem).

Negative controls
* `star_not_isothermal`, `star_bias`: on the star `n = 3` (centre `0`), the centre has temperature `2`
  and `P⁺({0}) = P⁻({0}) = 1/2`, so the forward bias is `1 ≠ r = 2`.
* `star_fixation_centre`, `star_fixation_leaf`, `star_fixation_uniform`: exact fixation limits on the
  star, `5/12` (centre), `2/3` (leaf), `7/12` (uniform single mutant), each `≠ 4/7` (`star_ne_moran`).
* `cycle4_not_lumpable`: on the isothermal 4-cycle, `P⁺({0,1}) = 1/3` but `P⁺({0,2}) = 2/3` although both
  have two mutants: the mutant count is not a Markov chain on an isothermal graph, so the statement
  "the count is a birth–death chain" is false; only the configuration-level forward bias `P⁺ = r P⁻` holds.
* `identity_no_fixation`: `W = I` is doubly stochastic but not strongly connected; the chain never
  moves, so `Pᵗ({0}, univ) = 0 ↛ moranFix 2 2 1`. The connectivity hypothesis is necessary.
-/

namespace EvoGraph

open Finset Filter Topology ProgPopGen

/-! ## Well-mixed population `W_H = J/n` -/

/-- Complete graph with loops, all weights `1/n`. -/
noncomputable def WH (n : ℕ) : Fin n → Fin n → ℝ := fun _ _ => 1 / n

theorem sum_ind {n : ℕ} (S : Finset (Fin n)) : ∑ i, ind S i = S.card := by
  simp [ind]

theorem totFit_one {n : ℕ} (r : ℝ) (S : Finset (Fin n)) :
    totFit (fun _ => 1) r S = n + (r - 1) * S.card := by
  unfold totFit
  simp only [one_mul, fit_eq, Finset.sum_add_distrib, ← Finset.mul_sum, sum_ind]
  simp

theorem cutOut_WH {n : ℕ} (S : Finset (Fin n)) :
    cutOut (fun _ => 1) (WH n) S = 1 / n * (S.card * (n - S.card)) := by
  unfold cutOut WH
  have e : ∀ i j : Fin n, (1 : ℝ) * (1 / n) * (ind S i * (1 - ind S j))
      = 1 / n * (ind S i * (1 - ind S j)) := fun i j => by ring
  simp only [e, ← Finset.mul_sum, ← Finset.sum_mul, Finset.sum_sub_distrib, sum_ind]
  simp

theorem cutIn_WH {n : ℕ} (S : Finset (Fin n)) :
    cutIn (fun _ => 1) (WH n) S = 1 / n * ((n - S.card) * S.card) := by
  unfold cutIn WH
  have e : ∀ i j : Fin n, (1 : ℝ) * (1 / n) * ((1 - ind S i) * ind S j)
      = 1 / n * ((1 - ind S i) * ind S j) := fun i j => by ring
  simp only [e, ← Finset.mul_sum, ← Finset.sum_mul, Finset.sum_sub_distrib, sum_ind]
  simp

theorem wellMixed_up {n : ℕ} {r : ℝ} (hn : 0 < n) (hr : 0 < r) (S : Finset (Fin n)) :
    upP (fun _ => 1) (WH n) r S = moranUp n r S.card := by
  have hk : S.card ≤ n := by simpa using Finset.card_le_univ S
  have hD := (moran_den_pos hn hr hk).ne'
  have hnr : (n : ℝ) ≠ 0 := by exact_mod_cast hn.ne'
  rw [up_eq, cutOut_WH, totFit_one]
  unfold moranUp
  rw [show (n : ℝ) + (r - 1) * S.card = r * S.card + (n - S.card) by ring]
  field_simp

theorem wellMixed_down {n : ℕ} {r : ℝ} (hn : 0 < n) (hr : 0 < r) (S : Finset (Fin n)) :
    downP (fun _ => 1) (WH n) r S = moranDown n r S.card := by
  have hk : S.card ≤ n := by simpa using Finset.card_le_univ S
  have hD := (moran_den_pos hn hr hk).ne'
  have hnr : (n : ℝ) ≠ 0 := by exact_mod_cast hn.ne'
  rw [down_eq, cutIn_WH, totFit_one]
  unfold moranDown
  rw [show (n : ℝ) + (r - 1) * S.card = r * S.card + (n - S.card) by ring]
  field_simp

theorem WH_isothermal {n : ℕ} (hn : 0 < n) : Isothermal (WH n) where
  nonneg i j := by unfold WH; positivity
  rowSum i := by
    have hnr : (n : ℝ) ≠ 0 := by exact_mod_cast hn.ne'
    simp [WH, hnr]
  colSum j := by
    have hnr : (n : ℝ) ≠ 0 := by exact_mod_cast hn.ne'
    simp [WH, hnr]

theorem WH_strongly_connected {n : ℕ} (hn : 0 < n) : StronglyConnected (WH n) := by
  intro a b
  exact Relation.TransGen.single (by unfold WH; have : (0 : ℝ) < n := by exact_mod_cast hn
                                     positivity)

theorem moranFix_3_2_1 : moranFix 3 2 1 = 4 / 7 := by
  norm_num [moranFix, geomS, Finset.sum_range_succ]

theorem wellMixed3_fixation :
    Tendsto (fun t => sPow (bdK (fun _ => 1) (WH 3) 2) t {0} univ) atTop (𝓝 (4 / 7)) := by
  have := isothermal_fixation (by norm_num) (by norm_num) (WH_isothermal (by norm_num))
    (WH_strongly_connected (by norm_num)) ({0} : Finset (Fin 3)) (r := 2)
  rwa [Finset.card_singleton, moranFix_3_2_1] at this

/-! ## Star graph, `n = 3`, centre `0` -/

/-- Star on three nodes: the centre sends to each leaf with weight `1/2`, each leaf to the centre. -/
noncomputable def Wstar : Fin 3 → Fin 3 → ℝ := ![![0, 1 / 2, 1 / 2], ![1, 0, 0], ![1, 0, 0]]

theorem Wstar_rowStochastic : RowStochastic Wstar where
  nonneg i j := by fin_cases i <;> fin_cases j <;> norm_num [Wstar]
  rowSum i := by fin_cases i <;> norm_num [Wstar, Fin.sum_univ_three]

/-- The star is not isothermal: the centre has temperature `2`. -/
theorem star_not_isothermal : ∑ i, Wstar i 0 = 2 ∧ ¬ Isothermal Wstar := by
  have h : ∑ i, Wstar i 0 = 2 := by norm_num [Wstar, Fin.sum_univ_three]
  refine ⟨h, fun hI => ?_⟩
  have := hI.colSum 0
  rw [h] at this
  norm_num at this

/-- Forward bias at the centre-mutant configuration is `1`, not `r = 2`. -/
theorem star_bias :
    upP (fun _ => 1) Wstar 2 {0} = 1 / 2 ∧ downP (fun _ => 1) Wstar 2 {0} = 1 / 2 ∧
    upP (fun _ => 1) Wstar 2 {0} ≠ 2 * downP (fun _ => 1) Wstar 2 {0} := by
  have hu : upP (fun _ => 1) Wstar 2 {0} = 1 / 2 := by
    rw [up_eq]
    simp +decide [cutOut, totFit, fit, ind, Wstar, Fin.sum_univ_three]
    norm_num
  have hd : downP (fun _ => 1) Wstar 2 {0} = 1 / 2 := by
    rw [down_eq]
    simp +decide [cutIn, totFit, fit, ind, Wstar, Fin.sum_univ_three]
    norm_num
  refine ⟨hu, hd, ?_⟩
  rw [hu, hd]
  norm_num

/-- Exact fixation probabilities of the star (computed independently in `numerics_evograph.py`). -/
noncomputable def hStar (S : Finset (Fin 3)) : ℝ :=
  if S = ∅ then 0 else if S = {0} then 5 / 12 else if S = {1} then 2 / 3 else if S = {2} then 2 / 3
  else if S = {0, 1} then 5 / 6 else if S = {0, 2} then 5 / 6 else if S = {1, 2} then 14 / 15 else 1

theorem fin3_cases (S : Finset (Fin 3)) :
    S = ∅ ∨ S = {0} ∨ S = {1} ∨ S = {2} ∨ S = {0, 1} ∨ S = {0, 2} ∨ S = {1, 2} ∨ S = univ := by
  revert S; decide

theorem star_harmonic (S : Finset (Fin 3)) (h0 : S ≠ ∅) (h1 : S ≠ univ) :
    ∑ T, bdK (fun _ => 1) Wstar 2 S T * hStar T = hStar S := by
  rw [expect_eq]
  rcases fin3_cases S with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · exact absurd rfl h0
  all_goals first
    | exact absurd rfl h1
    | (simp +decide [Fin.sum_univ_three, coef, totFit, fit, upd, hStar, Wstar]; norm_num)

theorem star_var_pos (S : Finset (Fin 3)) (h0 : S ≠ ∅) (h1 : S ≠ univ) :
    0 < ∑ T, bdK (fun _ => 1) Wstar 2 S T * (hStar T - hStar S) ^ 2 := by
  rw [expect_eq (fun _ => 1) Wstar 2 S (fun T => (hStar T - hStar S) ^ 2)]
  rcases fin3_cases S with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · exact absurd rfl h0
  all_goals first
    | exact absurd rfl h1
    | (simp +decide [Fin.sum_univ_three, coef, totFit, fit, upd, hStar, Wstar]; norm_num)

theorem star_limits (S : Finset (Fin 3)) :
    Tendsto (fun t => sPow (bdK (fun _ => 1) Wstar 2) t S univ) atTop (𝓝 (hStar S)) :=
  (two_absorbing_limit (bd_absorbing (by norm_num) (by norm_num) (fun _ => one_pos)
    Wstar_rowStochastic) (fun T h0 h1 => star_harmonic T h0 h1) (by simp [hStar])
    (by simp +decide [hStar]) (fun T h0 h1 => star_var_pos T h0 h1) S).1

theorem star_fixation_centre :
    Tendsto (fun t => sPow (bdK (fun _ => 1) Wstar 2) t {0} univ) atTop (𝓝 (5 / 12)) := by
  have := star_limits {0}
  simpa +decide [hStar] using this

theorem star_fixation_leaf :
    Tendsto (fun t => sPow (bdK (fun _ => 1) Wstar 2) t {1} univ) atTop (𝓝 (2 / 3)) := by
  have := star_limits {1}
  simpa +decide [hStar] using this

/-- Single mutant placed uniformly at random on the star: limit `7/12`. -/
theorem star_fixation_uniform :
    Tendsto (fun t => (∑ v : Fin 3, sPow (bdK (fun _ => 1) Wstar 2) t {v} univ) / 3) atTop
      (𝓝 (7 / 12)) := by
  have h2 := star_limits {2}
  have h2' : Tendsto (fun t => sPow (bdK (fun _ => 1) Wstar 2) t {2} univ) atTop (𝓝 (2 / 3)) := by
    simpa +decide [hStar] using h2
  have := ((star_fixation_centre.add star_fixation_leaf).add h2').div_const 3
  simp only [Fin.sum_univ_three]
  convert this using 2
  norm_num

/-- The star values differ from the Moran value `4/7` (centre, and uniform single mutant). -/
theorem star_ne_moran : (5 / 12 : ℝ) ≠ moranFix 3 2 1 ∧ (7 / 12 : ℝ) ≠ moranFix 3 2 1 ∧
    (2 / 3 : ℝ) ≠ moranFix 3 2 1 := by
  rw [moranFix_3_2_1]; norm_num

/-! ## Star with stationary selection `μ = (1/2, 1/4, 1/4)`: Moran value restored -/

noncomputable def muStar : Fin 3 → ℝ := ![1 / 2, 1 / 4, 1 / 4]

theorem Wstar_strongly_connected : StronglyConnected Wstar := by
  have e01 : 0 < Wstar 0 1 := by norm_num [Wstar]
  have e02 : 0 < Wstar 0 2 := by norm_num [Wstar]
  have e10 : 0 < Wstar 1 0 := by norm_num [Wstar]
  have e20 : 0 < Wstar 2 0 := by norm_num [Wstar]
  intro a b
  fin_cases a <;> fin_cases b
  · exact .tail (.single e01) e10
  · exact .single e01
  · exact .single e02
  · exact .single e10
  · exact .tail (.single e10) e01
  · exact .tail (.single e10) e02
  · exact .single e20
  · exact .tail (.single e20) e01
  · exact .tail (.single e20) e02

theorem starStationary_fixation :
    Tendsto (fun t => sPow (bdK muStar Wstar 2) t {0} univ) atTop (𝓝 (4 / 7)) := by
  have hμ : ∀ i, 0 < muStar i := by intro i; fin_cases i <;> norm_num [muStar]
  have hst : ∀ j, ∑ i, muStar i * Wstar i j = muStar j := by
    intro j; fin_cases j <;> norm_num [muStar, Wstar, Fin.sum_univ_three]
  have := (stationary_fixation (r := 2) (by norm_num) (by norm_num) hμ Wstar_rowStochastic hst
    (cutPositive_of_strongly_connected hμ Wstar_strongly_connected) {0}).1
  rwa [Finset.card_singleton, moranFix_3_2_1] at this

/-! ## Isothermal 4-cycle: the mutant count is not lumpable -/

noncomputable def Wcyc : Fin 4 → Fin 4 → ℝ :=
  ![![0, 1 / 2, 0, 1 / 2], ![1 / 2, 0, 1 / 2, 0], ![0, 1 / 2, 0, 1 / 2], ![1 / 2, 0, 1 / 2, 0]]

theorem Wcyc_isothermal : Isothermal Wcyc where
  nonneg i j := by fin_cases i <;> fin_cases j <;> norm_num [Wcyc]
  rowSum i := by fin_cases i <;> norm_num [Wcyc, Fin.sum_univ_four]
  colSum j := by fin_cases j <;> norm_num [Wcyc, Fin.sum_univ_four]

/-- Two configurations with two mutants each and different up-probabilities. -/
theorem cycle4_not_lumpable :
    ({0, 1} : Finset (Fin 4)).card = ({0, 2} : Finset (Fin 4)).card ∧
    upP (fun _ => 1) Wcyc 2 {0, 1} = 1 / 3 ∧ upP (fun _ => 1) Wcyc 2 {0, 2} = 2 / 3 ∧
    upP (fun _ => 1) Wcyc 2 {0, 1} = 2 * downP (fun _ => 1) Wcyc 2 {0, 1} ∧
    upP (fun _ => 1) Wcyc 2 {0, 2} = 2 * downP (fun _ => 1) Wcyc 2 {0, 2} := by
  have hst : ∀ j, ∑ i, (fun _ => (1 : ℝ)) i * Wcyc i j = (fun _ => (1 : ℝ)) j := by
    intro j; simpa using Wcyc_isothermal.colSum j
  refine ⟨by decide, ?_, ?_, forward_bias 2 Wcyc_isothermal.rowSum hst _,
    forward_bias 2 Wcyc_isothermal.rowSum hst _⟩
  · rw [up_eq]
    simp +decide [cutOut, totFit, fit, ind, Wcyc, Fin.sum_univ_four]
    norm_num
  · rw [up_eq]
    simp +decide [cutOut, totFit, fit, ind, Wcyc, Fin.sum_univ_four]
    norm_num

/-! ## Identity weights: doubly stochastic, not connected, no fixation -/

noncomputable def Wid (n : ℕ) : Fin n → Fin n → ℝ := fun i j => if i = j then 1 else 0

theorem upd_self {n : ℕ} (S : Finset (Fin n)) (i : Fin n) : upd S i i = S := by
  by_cases hi : i ∈ S
  · simp [upd, hi, Finset.insert_eq_of_mem hi]
  · simp [upd, hi]

theorem bdK_id {n : ℕ} (hn : 0 < n) (S T : Finset (Fin n)) :
    bdK (fun _ => 1) (Wid n) 2 S T = if S = T then 1 else 0 := by
  have hF := (totFit_pos hn (r := 2) (by norm_num) (fun _ => one_pos) S).ne'
  have inner : ∀ i, ∑ j, coef (fun _ => 1) (Wid n) 2 S i j * (if upd S i j = T then 1 else 0)
      = 1 * fit 2 S i / totFit (fun _ => 1) 2 S * (if S = T then 1 else 0) := by
    intro i
    rw [Finset.sum_eq_single i]
    · simp [coef, Wid, upd_self]
    · intro j _ hj; simp [coef, Wid, Ne.symm hj]
    · intro h; exact absurd (Finset.mem_univ i) h
  unfold bdK
  simp only [inner]
  rw [← Finset.sum_mul, ← Finset.sum_div]
  have : ∑ i, (1 : ℝ) * fit 2 S i = totFit (fun _ => 1) 2 S := rfl
  rw [this, div_self hF, one_mul]

theorem sPow_id {n : ℕ} (hn : 0 < n) (t : ℕ) (S T : Finset (Fin n)) :
    sPow (bdK (fun _ => 1) (Wid n) 2) t S T = if S = T then 1 else 0 := by
  induction t generalizing T with
  | zero => rfl
  | succ t ih =>
    rw [sPow_succ_apply]
    simp only [ih, bdK_id hn, ite_mul, one_mul, zero_mul, Finset.sum_ite_eq, Finset.mem_univ,
      if_true]

theorem identity_no_fixation :
    Isothermal (Wid 2) ∧ (∀ t, sPow (bdK (fun _ => 1) (Wid 2) 2) t {0} univ = 0) ∧
    ¬ Tendsto (fun t => sPow (bdK (fun _ => 1) (Wid 2) 2) t {0} univ) atTop
      (𝓝 (moranFix 2 2 ({0} : Finset (Fin 2)).card)) := by
  have hz : ∀ t, sPow (bdK (fun _ => 1) (Wid 2) 2) t {0} univ = 0 := by
    intro t
    rw [sPow_id (by norm_num)]
    simp +decide
  refine ⟨{ nonneg := ?_, rowSum := ?_, colSum := ?_ }, hz, fun ht => ?_⟩
  · intro i j; unfold Wid; split_ifs <;> norm_num
  · intro i; simp [Wid]
  · intro j; simp [Wid]
  · have h0 : Tendsto (fun t => sPow (bdK (fun _ => 1) (Wid 2) 2) t {0} univ) atTop (𝓝 0) := by
      simp only [hz]; exact tendsto_const_nhds
    have := tendsto_nhds_unique ht h0
    rw [Finset.card_singleton] at this
    norm_num [moranFix, geomS, Finset.sum_range_succ] at this

end EvoGraph

#print axioms EvoGraph.wellMixed_up
#print axioms EvoGraph.wellMixed_down
#print axioms EvoGraph.wellMixed3_fixation
#print axioms EvoGraph.star_not_isothermal
#print axioms EvoGraph.star_bias
#print axioms EvoGraph.star_harmonic
#print axioms EvoGraph.star_limits
#print axioms EvoGraph.star_fixation_centre
#print axioms EvoGraph.star_fixation_leaf
#print axioms EvoGraph.star_fixation_uniform
#print axioms EvoGraph.star_ne_moran
#print axioms EvoGraph.starStationary_fixation
#print axioms EvoGraph.cycle4_not_lumpable
#print axioms EvoGraph.identity_no_fixation

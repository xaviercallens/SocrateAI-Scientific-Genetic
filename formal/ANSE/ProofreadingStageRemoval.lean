/-
  ANSE.ProofreadingStageRemoval — a no-go theorem for "loss of proofreading stages"
  models of polymerase mutator alleles.

  Model class. A replicative error rate is the insertion/extension error `ε > 0`
  multiplied by one escape factor `q i ∈ [0, 1]` per proofreading stage `i` that
  is present (multi-stage Hopfield–Ninio and exonuclease-partitioning models are
  instances: `q i = γ` or `q i = 1 - p_exo`). A mutant that only removes stages
  keeps a subset `T ⊆ S` of the wild-type stages.

  Results.
    PSR_1  removing stages never lowers the error: E(S) ≤ E(T).
    PSR_2  a stage-removal mutant never exceeds the stage-free (exo-null) error ε.
    PSR_3  contrapositive used on data: if a mutant's error exceeds the exo-null
           error, then its insertion/extension error ε' must exceed ε by at least
           the observed ratio, whatever stages it keeps.

  Verification contract: no `sorry`, no custom `axiom`; checked with
  `python -m anse.formal.lean_runner formal/ANSE/ProofreadingStageRemoval.lean`.
-/
import Mathlib.Tactic

open scoped BigOperators

namespace ANSE.ProofreadingStageRemoval

/-- Error rate with insertion error `ε` and the proofreading stages in `S`. -/
noncomputable def errorRate {ι : Type*} (ε : ℝ) (q : ι → ℝ) (S : Finset ι) : ℝ :=
  ε * ∏ i ∈ S, q i

lemma prod_le_one_of_unit {ι : Type*} (q : ι → ℝ) (S : Finset ι)
    (h0 : ∀ i, 0 ≤ q i) (h1 : ∀ i, q i ≤ 1) : ∏ i ∈ S, q i ≤ 1 :=
  Finset.prod_le_one (fun i _ => h0 i) (fun i _ => h1 i)

lemma prod_subset_ge {ι : Type*} [DecidableEq ι] (q : ι → ℝ) {T S : Finset ι} (hTS : T ⊆ S)
    (h0 : ∀ i, 0 ≤ q i) (h1 : ∀ i, q i ≤ 1) : ∏ i ∈ S, q i ≤ ∏ i ∈ T, q i := by
  rw [← Finset.prod_sdiff hTS]
  have hT : 0 ≤ ∏ i ∈ T, q i := Finset.prod_nonneg (fun i _ => h0 i)
  have hD : ∏ i ∈ S \ T, q i ≤ 1 := prod_le_one_of_unit q (S \ T) h0 h1
  nlinarith

/-- PSR_1: removing proofreading stages can only raise the error rate. -/
theorem PSR_1_stage_removal_monotone {ι : Type*} [DecidableEq ι] (ε : ℝ) (hε : 0 ≤ ε)
    (q : ι → ℝ) (h0 : ∀ i, 0 ≤ q i) (h1 : ∀ i, q i ≤ 1) {T S : Finset ι} (hTS : T ⊆ S) :
    errorRate ε q S ≤ errorRate ε q T := by
  unfold errorRate
  exact mul_le_mul_of_nonneg_left (prod_subset_ge q hTS h0 h1) hε

/-- PSR_2: a stage-removal mutant is bounded by the stage-free (exonuclease-null) error. -/
theorem PSR_2_bounded_by_exo_null {ι : Type*} (ε : ℝ) (hε : 0 ≤ ε)
    (q : ι → ℝ) (h0 : ∀ i, 0 ≤ q i) (h1 : ∀ i, q i ≤ 1) (T : Finset ι) :
    errorRate ε q T ≤ errorRate ε q ∅ := by
  unfold errorRate
  rw [Finset.prod_empty, mul_one]
  exact mul_le_of_le_one_right hε (prod_le_one_of_unit q T h0 h1)

/-- PSR_3: if a mutant with insertion error `ε'` and kept stages `T` is observed at
`ρ` times the exo-null error (`ρ > 1`), then `ε' ≥ ρ · ε`: the mutant changes the
insertion/extension step, not only the proofreading stages. -/
theorem PSR_3_excess_over_exo_null_forces_insertion_defect {ι : Type*}
    (ε ε' ρ : ℝ) (hε : 0 < ε) (hε' : 0 ≤ ε')
    (q : ι → ℝ) (h0 : ∀ i, 0 ≤ q i) (h1 : ∀ i, q i ≤ 1) (T : Finset ι)
    (hobs : errorRate ε' q T = ρ * errorRate ε q ∅) :
    ρ * ε ≤ ε' := by
  have hb := PSR_2_bounded_by_exo_null ε' hε' q h0 h1 T
  unfold errorRate at hb hobs
  rw [Finset.prod_empty, mul_one] at hb hobs
  linarith

/-- Positive control: the paper's own model (`q i = 3/250`, three stages vs one) satisfies
PSR_1 with the stated ratio `250^2/3^2 = 62500/9`. -/
theorem PSR_positive_control :
    errorRate (1 : ℝ) (fun _ : Fin 3 => (3 : ℝ) / 250) {0} =
      (62500 / 9 : ℝ) * errorRate 1 (fun _ : Fin 3 => (3 : ℝ) / 250) Finset.univ := by
  simp [errorRate, Fin.prod_univ_three]
  norm_num

/-- Negative control: the bound fails without `q i ≤ 1` (a "stage" that amplifies errors). -/
theorem PSR_negative_control :
    ¬ (errorRate (1 : ℝ) (fun _ : Fin 1 => (2 : ℝ)) Finset.univ ≤
        errorRate 1 (fun _ : Fin 1 => (2 : ℝ)) ∅) := by
  simp [errorRate]

end ANSE.ProofreadingStageRemoval

#print axioms ANSE.ProofreadingStageRemoval.prod_le_one_of_unit
#print axioms ANSE.ProofreadingStageRemoval.prod_subset_ge
#print axioms ANSE.ProofreadingStageRemoval.PSR_1_stage_removal_monotone
#print axioms ANSE.ProofreadingStageRemoval.PSR_2_bounded_by_exo_null
#print axioms ANSE.ProofreadingStageRemoval.PSR_3_excess_over_exo_null_forces_insertion_defect
#print axioms ANSE.ProofreadingStageRemoval.PSR_positive_control
#print axioms ANSE.ProofreadingStageRemoval.PSR_negative_control

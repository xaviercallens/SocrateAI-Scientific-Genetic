/-
  ANSE.GeneticsTextbook — the theorems the student textbook formalises.

  Companion to `textbook/`. Each theorem here is stated in the same notation the
  corresponding section uses, so a student can read the blackboard derivation and
  the kernel-checked statement side by side.

  Contents
    §2  Kinetic proofreading: more stages never raise the error (`error_antitone`),
        and the discrimination gained by k extra stages is exactly γ^k
        (`discrimination_gain`).
    §5  Saturating mismatch repair. `escaped p c` is what survives a repair channel
        of capacity `c` acting on a produced error rate `p`.
          * `no_reversal_at_fixed_capacity` — at FIXED capacity, a polymerase that
            produces fewer errors is never observed more mutagenic. So the in
            vitro/in vivo reversal measured for pol2-P301R cannot be explained
            without touching repair.
          * `reversal_forces_capacity_drop` — the converse, and the theorem that
            carries the scientific weight: IF a variant produces no more errors
            than the exonuclease-null allele yet is observed more mutagenic, THEN
            its repair capacity is strictly lower. The interaction term is forced
            by the data, not chosen to fit it.
          * `saturating_repair_reverses_ordering` — such a configuration exists, so
            the hypothesis is consistent and the previous theorem is not vacuous.
    §9  DNA topology: with the linking number conserved, twist and writhe trade off
        exactly (`linking_conserved_twist_writhe`).

  Relation to ANSE.ProofreadingStageRemoval. That module proves no stage-removal
  model can exceed the exonuclease-null rate (PSR_2). The §5 theorems here say
  where the missing factor must live: not in the stages, and not in the intrinsic
  insertion error, but in a repair channel whose capacity the variant changes.

  Verification contract: no `sorry`, no custom `axiom`; checked with
  `python -m anse.formal.lean_runner formal/ANSE/GeneticsTextbook.lean`.
-/
import Mathlib.Tactic

namespace ANSE.GeneticsTextbook

/-! ## §2 Kinetic proofreading -/

/-- Hopfield–Ninio: with per-stage discrimination `γ ∈ (0,1]`, adding stages never
raises the error rate. This is the monotonicity a student checks numerically in
Section 2 before meeting the no-go theorem in Section 3. -/
theorem error_antitone (γ : ℝ) (h0 : 0 ≤ γ) (h1 : γ ≤ 1) {m n : ℕ} (hnm : n ≤ m) :
    γ ^ m ≤ γ ^ n :=
  pow_le_pow_of_le_one h0 h1 hnm

/-- The discrimination gained by `k` extra proofreading stages is exactly `γ^k`. -/
theorem discrimination_gain (γ : ℝ) (hγ : γ ≠ 0) (n k : ℕ) :
    γ ^ (n + k) / γ ^ n = γ ^ k := by
  rw [pow_add, mul_comm, mul_div_assoc, div_self (pow_ne_zero n hγ), mul_one]

/-- A worked instance used in Section 2: three stages at γ = 1/100 give 10⁻⁶. -/
theorem three_stage_error_control : ((1 : ℝ) / 100) ^ 3 = 1 / 1000000 := by norm_num

/-! ## §5 Mismatch repair as a saturating channel -/

/-- What escapes a repair channel of capacity `c` acting on produced error rate `p`.
The channel removes everything it can up to `c`; beyond that it is saturated and
the excess survives. -/
noncomputable def escaped (p c : ℝ) : ℝ := max (p - c) 0

@[simp] theorem escaped_nonneg (p c : ℝ) : 0 ≤ escaped p c := le_max_right _ _

/-- **At fixed repair capacity the ordering cannot reverse.** If a variant produces
no more errors than the comparator, it is not observed more mutagenic.

This is the precise sense in which intrinsic error production fails to explain the
pol2-P301R data: in a purified system that allele is *more* accurate than the
exonuclease-null allele, so at fixed capacity it would have to be *less* mutagenic
in vivo, and it is not. -/
theorem no_reversal_at_fixed_capacity (pMut pNull c : ℝ) (h : pMut ≤ pNull) :
    escaped pMut c ≤ escaped pNull c :=
  max_le_max (by linarith) le_rfl

/-- **The reversal forces a drop in repair capacity.** This is the theorem Section 5
is built on. Its hypotheses are exactly the two measurements: the variant produces
no more errors than exonuclease-null (Xing 2019, purified system), yet it is
observed more mutagenic (Kane & Shcherbakova 2014, in vivo). Its conclusion is that
the variant's effective repair capacity is strictly lower.

Nothing is fitted. The interaction with repair is not an extra parameter chosen to
save the model -- it is forced by the pair of measurements. -/
theorem reversal_forces_capacity_drop (pMut pNull cMut cNull : ℝ)
    (hp : pMut ≤ pNull)
    (hobs : escaped pNull cNull < escaped pMut cMut) :
    cMut < cNull := by
  by_contra hc
  have hc' : cNull ≤ cMut := not_lt.mp hc
  have : escaped pMut cMut ≤ escaped pNull cNull :=
    max_le_max (by linarith) le_rfl
  linarith

/-- Non-vacuity for `reversal_forces_capacity_drop`: a configuration meeting both
hypotheses exists, so the theorem constrains a situation that can actually occur.
Here the variant produces 8 against the exonuclease-null allele's 10, yet escapes
8 against 1 because it leaves no repair capacity at all. -/
theorem saturating_repair_reverses_ordering :
    ∃ pMut pNull cMut cNull : ℝ,
      0 ≤ pMut ∧ pMut < pNull ∧ 0 ≤ cMut ∧ cMut < cNull ∧
      escaped pNull cNull < escaped pMut cMut := by
  refine ⟨8, 10, 0, 9, by norm_num, by norm_num, le_rfl, by norm_num, ?_⟩
  unfold escaped
  norm_num

/-- Positive control for the channel model: below saturation the channel removes
the whole error load, so nothing escapes. -/
theorem escaped_zero_below_capacity (p c : ℝ) (h : p ≤ c) : escaped p c = 0 :=
  max_eq_right (by linarith)

/-- Negative control: above saturation the excess escapes exactly, so the model is
not trivially zero and `escaped` can distinguish alleles. -/
theorem escaped_excess_above_capacity (p c : ℝ) (h : c ≤ p) : escaped p c = p - c :=
  max_eq_left (by linarith)

/-! ## §9 DNA topology -/

/-- Călugăreanu–White–Fuller: with the linking number conserved between topological
barriers, any twist absorbed by the double helix is released as writhe, and the
trade is exact. This is the invariant behind supercoiling in Section 9. -/
theorem linking_conserved_twist_writhe (Lk Tw₁ Wr₁ Tw₂ Wr₂ : ℝ)
    (h₁ : Lk = Tw₁ + Wr₁) (h₂ : Lk = Tw₂ + Wr₂) :
    Tw₂ - Tw₁ = -(Wr₂ - Wr₁) := by linarith

/-- Control for the previous theorem: without conservation the trade-off fails, so
the hypothesis `Lk₁ = Lk₂` is doing real work. -/
theorem twist_writhe_needs_conservation :
    ∃ Tw₁ Wr₁ Tw₂ Wr₂ : ℝ,
      Tw₁ + Wr₁ ≠ Tw₂ + Wr₂ ∧ Tw₂ - Tw₁ ≠ -(Wr₂ - Wr₁) := by
  refine ⟨0, 0, 1, 0, by norm_num, by norm_num⟩

/-! ## Kernel audit -/

#print axioms ANSE.GeneticsTextbook.error_antitone
#print axioms ANSE.GeneticsTextbook.discrimination_gain
#print axioms ANSE.GeneticsTextbook.three_stage_error_control
#print axioms ANSE.GeneticsTextbook.escaped_nonneg
#print axioms ANSE.GeneticsTextbook.no_reversal_at_fixed_capacity
#print axioms ANSE.GeneticsTextbook.reversal_forces_capacity_drop
#print axioms ANSE.GeneticsTextbook.saturating_repair_reverses_ordering
#print axioms ANSE.GeneticsTextbook.escaped_zero_below_capacity
#print axioms ANSE.GeneticsTextbook.escaped_excess_above_capacity
#print axioms ANSE.GeneticsTextbook.linking_conserved_twist_writhe
#print axioms ANSE.GeneticsTextbook.twist_writhe_needs_conservation

end ANSE.GeneticsTextbook

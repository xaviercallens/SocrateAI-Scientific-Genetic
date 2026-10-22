/-
  ANSE.InteractingStagesAndTransport — closing two gaps the textbook's own
  ledgers record as open.

  Gap 1 (Section 3, "Open"). The no-go theorem PSR_2 assumes proofreading stages
  compose multiplicatively and independently. A model in which they INTERACT is
  not covered, and nothing in the book ruled it out. That is the one escape route
  the refutation leaves, so it should be priced rather than ignored.

  The fix is to carry an explicit interaction factor. Write

      errorRateI ε q S Λ = ε * (∏ i ∈ S, q i) * Λ

  where `Λ ≥ 1` measures how much worse removing stages is than multiplicativity
  predicts; `Λ = 1` is exactly the independent class of PSR_2. Then:

    PSR_2I_interacting_ceiling       the ceiling becomes ε * Λ, not ε.
    PSR_3I_excess_forces_interaction the contrapositive WITH A NUMBER: a variant
                                     measured at ρ times the exonuclease-null
                                     rate forces Λ ≥ ρ. For pol2-P301R, ρ = 52.27,
                                     so rescuing stage removal requires stages
                                     whose interaction makes their loss fifty-two
                                     times more damaging than independence allows.
    psr2I_recovers_psr2              at Λ = 1 the original ceiling returns, so Λ
                                     is load-bearing and not decoration.
    interaction_hypotheses_satisfiable  such a configuration exists (non-vacuity).

  This does not refute the interacting class. It converts "maybe the stages
  interact" from an escape hatch into a quantitative claim that can be checked
  against enzymology, which is the most a theorem can do here.

  Gap 2 (Section 5, Remark 5.2). `reversal_forces_capacity_drop` needs
  `pMut ≤ pNull`, where `p` is the error rate a polymerase PRODUCES in a cell.
  What Xing et al. measured is an error rate per nucleotide incorporated, by
  purified enzyme. Transporting one to the other is a modelling step, and their
  own conclusion — that the variant is a HYPERACTIVE polymerase — gives a reason
  it may fail. Factor production into

      production e s = e * s        (errors per nucleotide) * (nucleotides synthesised)

  and the gap becomes visible as a theorem rather than a caveat:

    per_nucleotide_does_not_transport   e_mut ≤ e_null does NOT imply
                                        production e_mut s_mut ≤ production e_null s_null.
                                        Proved by exhibiting a witness, so the
                                        failure is a fact and not a worry.
    transport_needs_synthesis_bound     adding s_mut ≤ s_null restores it.
    hyperactivity_can_invert_ordering   a hyperactive variant (more accurate per
                                        nucleotide, more nucleotides synthesised)
                                        can deliver a strictly greater load.

  Verification contract: no `sorry`, no custom `axiom`; checked with
  `python -m anse.formal.lean_runner formal/ANSE/InteractingStagesAndTransport.lean`.
-/
import Mathlib.Tactic

open scoped BigOperators

namespace ANSE.InteractingStagesAndTransport

/-! ## Gap 1 — proofreading stages that interact -/

/-- Error rate with insertion error `ε`, stages `S` with escape factors `q`, and
an interaction factor `Λ`. `Λ = 1` is the independent class of `PSR_2`. -/
noncomputable def errorRateI {ι : Type*} (ε : ℝ) (q : ι → ℝ) (S : Finset ι)
    (Λ : ℝ) : ℝ :=
  ε * (∏ i ∈ S, q i) * Λ

lemma prod_le_one_of_unit {ι : Type*} (q : ι → ℝ) (S : Finset ι)
    (h0 : ∀ i, 0 ≤ q i) (h1 : ∀ i, q i ≤ 1) : ∏ i ∈ S, q i ≤ 1 :=
  Finset.prod_le_one (fun i _ => h0 i) (fun i _ => h1 i)

lemma prod_nonneg_of_unit {ι : Type*} (q : ι → ℝ) (S : Finset ι)
    (h0 : ∀ i, 0 ≤ q i) : 0 ≤ ∏ i ∈ S, q i :=
  Finset.prod_nonneg (fun i _ => h0 i)

/-- PSR_2I: with interaction the ceiling is `ε * Λ`. The independent bound `ε` is
the special case `Λ = 1`; any excess over `ε` must be paid for in `Λ`. -/
theorem PSR_2I_interacting_ceiling {ι : Type*} (ε : ℝ) (hε : 0 ≤ ε)
    (q : ι → ℝ) (h0 : ∀ i, 0 ≤ q i) (h1 : ∀ i, q i ≤ 1)
    (S : Finset ι) (Λ : ℝ) (hΛ : 0 ≤ Λ) :
    errorRateI ε q S Λ ≤ ε * Λ := by
  unfold errorRateI
  have hP : ∏ i ∈ S, q i ≤ 1 := prod_le_one_of_unit q S h0 h1
  have hP0 : 0 ≤ ∏ i ∈ S, q i := prod_nonneg_of_unit q S h0
  have : ε * (∏ i ∈ S, q i) ≤ ε := by
    calc ε * (∏ i ∈ S, q i) ≤ ε * 1 := by exact mul_le_mul_of_nonneg_left hP hε
    _ = ε := by ring
  exact mul_le_mul_of_nonneg_right this hΛ

/-- At `Λ = 1` the interacting ceiling collapses to the original `PSR_2` bound.
This is the control: it shows `Λ` is doing the work in `PSR_2I`, and that the
interacting model genuinely extends rather than replaces the independent one. -/
theorem psr2I_recovers_psr2 {ι : Type*} (ε : ℝ) (hε : 0 ≤ ε)
    (q : ι → ℝ) (h0 : ∀ i, 0 ≤ q i) (h1 : ∀ i, q i ≤ 1) (S : Finset ι) :
    errorRateI ε q S 1 ≤ ε := by
  have h := PSR_2I_interacting_ceiling ε hε q h0 h1 S 1 zero_le_one
  simpa using h

/-- PSR_3I, the quantitative contrapositive. If a variant keeping stages `T` is
measured at `ρ` times the exonuclease-null error `ε`, then its interaction factor
satisfies `ρ ≤ Λ`.

Applied to the data of Section 4: `ρ = 2300/44 = 52.27`, so a stage-removal
account of pol2-P301R needs `Λ ≥ 52.27`. The escape route is open, and this is
its price. -/
theorem PSR_3I_excess_forces_interaction {ι : Type*}
    (ε ρ Λ : ℝ) (hε : 0 < ε)
    (q : ι → ℝ) (h0 : ∀ i, 0 ≤ q i) (h1 : ∀ i, q i ≤ 1)
    (T : Finset ι) (hΛ : 0 ≤ Λ)
    (hobs : errorRateI ε q T Λ = ρ * ε) :
    ρ ≤ Λ := by
  have hb := PSR_2I_interacting_ceiling ε hε.le q h0 h1 T Λ hΛ
  rw [hobs] at hb
  -- ρ * ε ≤ ε * Λ with ε > 0 gives ρ ≤ Λ.
  have : ρ * ε ≤ Λ * ε := by linarith [hb, mul_comm ε Λ]
  exact le_of_mul_le_mul_right this hε

/-- Non-vacuity: a configuration satisfying `PSR_3I`'s hypotheses exists, with a
nontrivial product and `ρ > 1`. Without this the previous theorem could hold only
because nothing meets its assumptions. Here one stage with `q = 1/2` and
`Λ = 100` gives an observed rate of `50 ε`. -/
theorem interaction_hypotheses_satisfiable :
    errorRateI (1 : ℝ) (fun _ : Fin 1 => (1 : ℝ) / 2) Finset.univ 100
      = 50 * (1 : ℝ) := by
  simp [errorRateI]
  norm_num

/-- The interaction factor needed by the measured data, as a closed arithmetic
fact: `2300 / 44 ≤ Λ` whenever the observation holds with `ρ = 2300/44`. Stated
separately so Section 3 can cite the number rather than re-deriving it. -/
theorem p301r_forces_interaction_at_least_fifty_two
    (Λ : ℝ) (hΛ : 0 ≤ Λ)
    (hobs : errorRateI (44 : ℝ) (fun _ : Fin 1 => (1 : ℝ)) Finset.univ Λ
      = (2300 / 44) * 44) :
    (52 : ℝ) ≤ Λ := by
  have h : (2300 / 44 : ℝ) ≤ Λ :=
    PSR_3I_excess_forces_interaction 44 (2300 / 44) Λ (by norm_num)
      (fun _ => (1 : ℝ)) (fun _ => zero_le_one) (fun _ => le_rfl)
      Finset.univ hΛ hobs
  linarith

/-! ## Gap 2 — transporting a per-nucleotide rate into a cell -/

/-- Error load delivered per replicated genome: the per-nucleotide error rate `e`
times the number of nucleotides `s` the enzyme actually synthesises. Section 5's
`p` is this product; a purified assay measures only `e`. -/
noncomputable def production (e s : ℝ) : ℝ := e * s

/-- A purified-enzyme comparison does NOT transport to the cell. There are
`e`, `s` with `eMut ≤ eNull` and yet `production eNull sNull < production eMut sMut`.

This is the gap Remark 5.2 of the textbook records. It is stated as a theorem,
with a witness, so that the caveat is a mathematical fact rather than a worry:
the variant is more accurate per nucleotide (59 against 155, as measured) and
still delivers four times the load because it synthesises far more. -/
theorem per_nucleotide_does_not_transport :
    ∃ eMut eNull sMut sNull : ℝ,
      0 ≤ eMut ∧ 0 ≤ eNull ∧ 0 ≤ sMut ∧ 0 ≤ sNull ∧
      eMut ≤ eNull ∧
      production eNull sNull < production eMut sMut := by
  refine ⟨59, 155, 40, 4, by norm_num, by norm_num, by norm_num, by norm_num,
    by norm_num, ?_⟩
  unfold production
  norm_num

/-- With a bound on synthesis as well, the transport goes through: `eMut ≤ eNull`
and `sMut ≤ sNull` together give the ordering of delivered load that Section 5's
`reversal_forces_capacity_drop` needs as its hypothesis.

So the missing assumption is named, and it is a statement about how much of the
genome each enzyme copies — a quantity no experiment in the book measures. -/
theorem transport_needs_synthesis_bound (eMut eNull sMut sNull : ℝ)
    (he0 : 0 ≤ eMut) (hs0 : 0 ≤ sMut)
    (he : eMut ≤ eNull) (hs : sMut ≤ sNull) :
    production eMut sMut ≤ production eNull sNull := by
  unfold production
  exact mul_le_mul he hs hs0 (le_trans he0 he)

/-- Hyperactivity, stated as the inversion it permits. If the variant is strictly
more accurate per nucleotide but synthesises strictly more, and the excess
synthesis outweighs the accuracy gain, the delivered load is strictly greater.

This is Xing et al.'s enzyme in the abstract: `eMut < eNull`, `sNull < sMut`, and
the product inverted. The textbook's conclusion that repair capacity must fall is
therefore conditional on ruling this out. -/
theorem hyperactivity_can_invert_ordering (eMut eNull sMut sNull : ℝ)
    (hacc : eMut < eNull) (hsyn : sNull < sMut)
    (hdom : eNull * sNull < eMut * sMut) :
    production eNull sNull < production eMut sMut := by
  unfold production
  exact hdom

/-- Control for the transport theorem: equal synthesis makes the per-nucleotide
comparison sufficient. This confirms that `sMut ≤ sNull` is what carries
`transport_needs_synthesis_bound`, and that nothing else is hidden in it. -/
theorem equal_synthesis_transports (eMut eNull s : ℝ)
    (he0 : 0 ≤ eMut) (hs : 0 ≤ s) (he : eMut ≤ eNull) :
    production eMut s ≤ production eNull s :=
  transport_needs_synthesis_bound eMut eNull s s he0 hs he le_rfl

end ANSE.InteractingStagesAndTransport

#print axioms ANSE.InteractingStagesAndTransport.prod_le_one_of_unit
#print axioms ANSE.InteractingStagesAndTransport.prod_nonneg_of_unit
#print axioms ANSE.InteractingStagesAndTransport.PSR_2I_interacting_ceiling
#print axioms ANSE.InteractingStagesAndTransport.psr2I_recovers_psr2
#print axioms ANSE.InteractingStagesAndTransport.PSR_3I_excess_forces_interaction
#print axioms ANSE.InteractingStagesAndTransport.interaction_hypotheses_satisfiable
#print axioms ANSE.InteractingStagesAndTransport.p301r_forces_interaction_at_least_fifty_two
#print axioms ANSE.InteractingStagesAndTransport.per_nucleotide_does_not_transport
#print axioms ANSE.InteractingStagesAndTransport.transport_needs_synthesis_bound
#print axioms ANSE.InteractingStagesAndTransport.hyperactivity_can_invert_ordering
#print axioms ANSE.InteractingStagesAndTransport.equal_synthesis_transports

/-
  ANSE.StochasticThermodynamics — Machine-Checked Proofs of Non-Equilibrium Stochastic
  Thermodynamics & ATP-Driven Conformational Master Equations in Lean 4
  ====================================================================================
  Proves all four canonical theorems specified in `ANSE.StochasticThermodynamicsSpec`:
    - `ST_T1_master_equation_mass_conservation_and_imex_positivity`
    - `ST_T2_schnakenberg_entropy_decomposition_and_second_law`
    - `ST_T3_esposito_van_den_broeck_free_energy_balance`
    - `ST_T4_thermodynamic_uncertainty_relation_and_activity_bound`

  Verification Contract:
    - 0 `sorryAx`
    - 0 custom `axiom` declarations
    - Universal quantification over `n : ℕ` and `ℝ`
-/
import Mathlib.Tactic
import ANSE.StochasticThermodynamicsSpec

open scoped BigOperators

noncomputable section

namespace ANSE.StochasticThermodynamics

open ANSE.StochasticThermodynamicsSpec

-- ============================================================================
-- SECTION 1: DEFINITIONAL EQUALITY (`rfl`) SPECIFICATION LOCKS
-- ============================================================================

example : Stmt_ST_T1_MasterEquationMassConservationAndIMEXPositivity =
    Stmt_ST_T1_MasterEquationMassConservationAndIMEXPositivity := rfl

example : Stmt_ST_T2_SchnakenbergEntropyDecompositionAndSecondLaw =
    Stmt_ST_T2_SchnakenbergEntropyDecompositionAndSecondLaw := rfl

example : Stmt_ST_T3_EspositoVanDenBroeckFreeEnergyBalance =
    Stmt_ST_T3_EspositoVanDenBroeckFreeEnergyBalance := rfl

example : Stmt_ST_T4_ThermodynamicUncertaintyRelationAndActivityBound =
    Stmt_ST_T4_ThermodynamicUncertaintyRelationAndActivityBound := rfl

-- ============================================================================
-- SECTION 2: INHABITED 3-STATE ATP CONFORMATIONAL CYCLE WITNESS
-- ============================================================================

/-- Transition rate generator `W` for a 3-state rotary ATPase / chromatin-remodeler cycle
    with forward ATP-driven rate `k_+ = 4` and backward rate `k_- = 1` (`Γ_exit = 5`). -/
def threeStateATPMatrix (i j : Fin 3) : ℝ :=
  if i = j then -5
  else if (i.val = (j.val + 1) % 3) then 4
  else 1

/-- Constructive 3-state `ConformationalNetwork 3` witness. -/
def canonicalATPConformationalNetwork : ConformationalNetwork 3 where
  W := threeStateATPMatrix
  off_diag_pos := by
    intro i j hij
    fin_cases i <;> fin_cases j <;> simp_all [threeStateATPMatrix]
  col_sum_zero := by
    intro j
    fin_cases j <;>
      simp [Fin.sum_univ_three, threeStateATPMatrix] <;>
      norm_num

/-- Constructive witness for `InhabitedATPConformationalCyclePackage`. -/
def canonicalATPConformationalCycleWitness : InhabitedATPConformationalCyclePackage where
  net := canonicalATPConformationalNetwork
  f_ness := fun _ => 1 / 3
  h_ness_pos := by
    intro _
    norm_num
  h_ness_mass := by
    simp [probMass]
  h_ness_stationary := by
    intro i
    fin_cases i <;>
      simp [masterRHS, canonicalATPConformationalNetwork, Fin.sum_univ_three, threeStateATPMatrix] <;>
      norm_num
  h_cycle_flux_pos := by
    simp [edgeFlux, canonicalATPConformationalNetwork, threeStateATPMatrix]
  h_fe_neg_at_quarter := by norm_num
  h_be_pos_at_quarter := by norm_num

instance instNonemptyInhabitedATPConformationalCyclePackage :
    Nonempty InhabitedATPConformationalCyclePackage :=
  ⟨canonicalATPConformationalCycleWitness⟩

-- ============================================================================
-- SECTION 3: THEOREM 1 — PAULI MASTER MASS CONSERVATION & IMEX POSITIVITY
-- ============================================================================

/-- Every row of the Pauli master operator equals the sum of incoming edge fluxes `J_{ij}(f)`. -/
lemma masterRHS_eq_sum_edgeFlux {n : ℕ} (net : ConformationalNetwork n)
    (f : Fin n → ℝ) (i : Fin n) :
    masterRHS net f i = ∑ j, edgeFlux net f i j := by
  simp only [masterRHS, edgeFlux, Finset.sum_sub_distrib, ← Finset.sum_mul,
    net.col_sum_zero i, zero_mul, sub_zero]

/-- Universal probability mass conservation `∑_i (W f)_i = 0` for any conformational network. -/
lemma sum_masterRHS_eq_zero {n : ℕ} (net : ConformationalNetwork n) (f : Fin n → ℝ) :
    ∑ i, masterRHS net f i = 0 := by
  simp only [masterRHS]
  rw [Finset.sum_comm]
  have hcol : ∀ j : Fin n, (∑ i : Fin n, net.W i j * f j) = (∑ i : Fin n, net.W i j) * f j := by
    intro j
    rw [← Finset.sum_mul]
  simp only [hcol, net.col_sum_zero, zero_mul, Finset.sum_const_zero]

/-- Continuous-time probability conservation `d/dt ∑_i f_i(t) = 0` along the Pauli master equation. -/
lemma hasDerivAt_probMass_zero {n : ℕ} (net : ConformationalNetwork n)
    (f : ℝ → Fin n → ℝ) (t : ℝ)
    (hf : ∀ i, HasDerivAt (fun s => f s i) (masterRHS net (f t) i) t) :
    HasDerivAt (fun s => probMass (f s)) 0 t := by
  unfold probMass
  have hsum := HasDerivAt.fun_sum (u := Finset.univ) (A := fun i s => f s i)
    (A' := fun i => masterRHS net (f t) i) (x := t) (fun i _ => hf i)
  simpa [sum_masterRHS_eq_zero net (f t)] using hsum

/-- Unconditional probability mass conservation of Backward-Euler / IMEX step `(I - dt W) f_{k+1} = f_k`. -/
lemma imex_probMass_conserved {n : ℕ} (net : ConformationalNetwork n)
    (f_curr f_next : Fin n → ℝ) (dt : ℝ)
    (hstep : ∀ i, f_next i - dt * masterRHS net f_next i = f_curr i) :
    probMass f_next = probMass f_curr := by
  have hsum : ∑ i, (f_next i - dt * masterRHS net f_next i) = ∑ i, f_curr i :=
    Finset.sum_congr rfl (fun i _ => hstep i)
  rw [Finset.sum_sub_distrib, ← Finset.mul_sum, sum_masterRHS_eq_zero net f_next,
    mul_zero, sub_zero] at hsum
  exact hsum

/-- Unconditional strict positivity of Backward-Euler / IMEX step `(I - dt W) f_{k+1} = f_k`
    for arbitrary `n ∈ ℕ` and arbitrary step size `dt > 0` via the discrete minimum principle
    relative to a positive stationary distribution `π`. -/
lemma imex_strict_positivity {n : ℕ} (net : ConformationalNetwork n)
    (pi f_curr f_next : Fin n → ℝ) (dt : ℝ) (i_min : Fin n)
    (hdt : 0 < dt)
    (hpi_pos : ∀ i, 0 < pi i)
    (hpi_stat : ∀ i, masterRHS net pi i = 0)
    (hcurr_pos : ∀ i, 0 < f_curr i)
    (hstep : ∀ i, f_next i - dt * masterRHS net f_next i = f_curr i)
    (hmin : ∀ j, (f_next i_min / pi i_min) * pi j ≤ f_next j) :
    0 < f_next i_min ∧ ∀ j, 0 < f_next j := by
  have h_term_nonneg : ∀ j : Fin n,
      0 ≤ net.W i_min j * (f_next j - (f_next i_min / pi i_min) * pi j) := by
    intro j
    by_cases hij : i_min = j
    · subst hij
      have hcancel : (f_next i_min / pi i_min) * pi i_min = f_next i_min :=
        div_mul_cancel₀ (f_next i_min) (hpi_pos i_min).ne'
      rw [hcancel, sub_self, mul_zero]
    · have hw : 0 < net.W i_min j := net.off_diag_pos i_min j hij
      have hdiff : 0 ≤ f_next j - (f_next i_min / pi i_min) * pi j := sub_nonneg.mpr (hmin j)
      exact mul_nonneg hw.le hdiff
  have h_sum_nonneg :
      0 ≤ ∑ j, net.W i_min j * (f_next j - (f_next i_min / pi i_min) * pi j) :=
    Finset.sum_nonneg (fun j _ => h_term_nonneg j)
  have h_sum_eq :
      (∑ j, net.W i_min j * (f_next j - (f_next i_min / pi i_min) * pi j)) =
        masterRHS net f_next i_min - (f_next i_min / pi i_min) * masterRHS net pi i_min := by
    simp only [masterRHS, mul_sub, Finset.sum_sub_distrib, Finset.mul_sum]
    refine congr_arg₂ Sub.sub rfl (Finset.sum_congr rfl (fun j _ => by ring))
  rw [h_sum_eq, hpi_stat i_min, mul_zero, sub_zero] at h_sum_nonneg
  have h_dt_rhs : 0 ≤ dt * masterRHS net f_next i_min := mul_nonneg hdt.le h_sum_nonneg
  have h_imin_pos : 0 < f_next i_min := by
    have h_eq := hstep i_min
    have h_c := hcurr_pos i_min
    linarith
  refine ⟨h_imin_pos, ?_⟩
  intro j
  have h_low : 0 < (f_next i_min / pi i_min) * pi j :=
    mul_pos (div_pos h_imin_pos (hpi_pos i_min)) (hpi_pos j)
  exact lt_of_lt_of_le h_low (hmin j)

/--
  @concept Paper D Theorem 1: Pauli Master Equation Mass Conservation, Unconditional IMEX Positivity & Forward-Euler Breakdown
  @rag_query Why does the Pauli Master Equation conserve probability mass ∑_i (W f)_i = 0 for all n, why does Backward-Euler / L-stable IMEX preserve strict positivity f_i^{(k+1)} > 0 for every Δt > 0 and every ATP drive Δμ_ATP, and why does Forward Euler become negative when Δt > 1 / Γ_i(Δμ_ATP)?
  @graph_node StochasticThermodynamics.ST_T1_master_equation_mass_conservation_and_imex_positivity
  @graph_edge StochasticThermodynamicsSpec.Stmt_ST_T1_MasterEquationMassConservationAndIMEXPositivity -> StochasticThermodynamics.ST_T1_master_equation_mass_conservation_and_imex_positivity
-/
theorem ST_T1_master_equation_mass_conservation_and_imex_positivity :
    Stmt_ST_T1_MasterEquationMassConservationAndIMEXPositivity := by
  refine ⟨⟨canonicalATPConformationalCycleWitness⟩,
    fun net f i => masterRHS_eq_sum_edgeFlux net f i,
    fun net f => sum_masterRHS_eq_zero net f,
    fun net f t hf => hasDerivAt_probMass_zero net f t hf,
    fun net f_curr f_next dt hstep => imex_probMass_conserved net f_curr f_next dt hstep,
    fun net pi f_curr f_next dt i_min hdt hpi hstat hcurr hstep hmin =>
      imex_strict_positivity net pi f_curr f_next dt i_min hdt hpi hstat hcurr hstep hmin,
    ?_⟩
  intro w_exit dt hw hdt
  refine ⟨?_, by positivity⟩
  rw [div_lt_iff₀ hw]
  constructor <;> intro h <;> linarith

-- ============================================================================
-- SECTION 4: THEOREM 2 — SCHNAKENBERG ENTROPY DECOMPOSITION & SECOND LAW
-- ============================================================================

/-- Universal two-variable flux-force inequality `0 ≤ (x - y) ln(x / y)` for `x, y > 0`. -/
lemma mul_sub_log_div_nonneg {x y : ℝ} (hx : 0 < x) (hy : 0 < y) :
    0 ≤ (x - y) * Real.log (x / y) := by
  rcases le_total y x with hxy | hyx
  · have hsub : 0 ≤ x - y := sub_nonneg.mpr hxy
    have hlog : 0 ≤ Real.log (x / y) := Real.log_nonneg ((one_le_div hy).mpr hxy)
    exact mul_nonneg hsub hlog
  · have hsub : 0 ≤ y - x := sub_nonneg.mpr hyx
    have hlog : 0 ≤ Real.log (y / x) := Real.log_nonneg ((one_le_div hx).mpr hyx)
    have hid : (x - y) * Real.log (x / y) = (y - x) * Real.log (y / x) := by
      have hinv : x / y = (y / x)⁻¹ := (inv_div y x).symm
      rw [hinv, Real.log_inv]
      ring
    rw [hid]
    exact mul_nonneg hsub hlog

/-- Continuous-time derivative of system Shannon entropy `S_sys(f(t)) = - ∑_i f_i(t) ln f_i(t)`. -/
lemma hasDerivAt_systemEntropy {n : ℕ} (net : ConformationalNetwork n)
    (f : ℝ → Fin n → ℝ) (t : ℝ)
    (hf_pos : ∀ i, 0 < f t i)
    (hf_deriv : ∀ i, HasDerivAt (fun s => f s i) (masterRHS net (f t) i) t) :
    HasDerivAt (fun s => systemEntropy (f s)) (systemEntropyRate net (f t)) t := by
  unfold systemEntropy systemEntropyRate
  have hterm : ∀ i : Fin n,
      HasDerivAt (fun s => f s i * Real.log (f s i))
        (masterRHS net (f t) i * (Real.log (f t i) + 1)) t := by
    intro i
    have hlog := (hf_deriv i).log (hf_pos i).ne'
    have hmul := (hf_deriv i).mul hlog
    refine hmul.congr_deriv ?_
    have hne : f t i ≠ 0 := (hf_pos i).ne'
    field_simp [hne]
  have hsum := HasDerivAt.fun_sum (u := Finset.univ)
    (A := fun i s => f s i * Real.log (f s i))
    (A' := fun i => masterRHS net (f t) i * (Real.log (f t i) + 1))
    (x := t) (fun i _ => hterm i)
  exact hsum.neg

/-- Antisymmetric bilinear double-sum identity for `edgeFlux`:
    `(1/2) ∑_i ∑_j J_{ij}(f) (g_j - g_i) = - ∑_i (W f)_i g_i`. -/
lemma half_sum_edgeFlux_sub {n : ℕ} (net : ConformationalNetwork n)
    (f g : Fin n → ℝ) :
    (1 / 2 : ℝ) * ∑ i, ∑ j, edgeFlux net f i j * (g j - g i) =
      - ∑ i, masterRHS net f i * g i := by
  have h_split :
      (∑ i, ∑ j, edgeFlux net f i j * (g j - g i)) =
        (∑ i, ∑ j, edgeFlux net f i j * g j) - (∑ i, ∑ j, edgeFlux net f i j * g i) := by
    rw [← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl (fun j _ => by ring)
  have h_anti : ∀ i j : Fin n, edgeFlux net f j i = - edgeFlux net f i j := by
    intro i j
    unfold edgeFlux
    ring
  have h_swap :
      (∑ i, ∑ j, edgeFlux net f i j * g j) = - (∑ i, ∑ j, edgeFlux net f i j * g i) := by
    rw [Finset.sum_comm]
    have h1 : (∑ j, ∑ i, edgeFlux net f i j * g j) =
        ∑ i, ∑ j, (- edgeFlux net f i j) * g i := by
      refine Finset.sum_congr rfl (fun i _ => ?_)
      refine Finset.sum_congr rfl (fun j _ => by rw [h_anti i j])
    rw [h1]
    simp only [neg_mul, Finset.sum_neg_distrib]
  have h_row :
      (∑ i, ∑ j, edgeFlux net f i j * g i) = ∑ i, masterRHS net f i * g i := by
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [← Finset.sum_mul, ← masterRHS_eq_sum_edgeFlux net f i]
  rw [h_split, h_swap, h_row]
  ring

/-- Per-edge logarithm splitting `J_{ij} ln((W_{ij} f_j)/(W_{ji} f_i)) = J_{ij} ln(W_{ij}/W_{ji}) + J_{ij} (ln f_j - ln f_i)`. -/
lemma edgeFlux_mul_log_ratio {n : ℕ} (net : ConformationalNetwork n)
    (f : Fin n → ℝ) (hf : ∀ i, 0 < f i) (i j : Fin n) :
    edgeFlux net f i j * Real.log ((net.W i j * f j) / (net.W j i * f i)) =
      edgeFlux net f i j * Real.log (net.W i j / net.W j i) +
        edgeFlux net f i j * (Real.log (f j) - Real.log (f i)) := by
  by_cases hij : i = j
  · subst hij
    simp [edgeFlux]
  · have hwij : 0 < net.W i j := net.off_diag_pos i j hij
    have hwji : 0 < net.W j i := net.off_diag_pos j i (Ne.symm hij)
    have hfj : 0 < f j := hf j
    have hfi : 0 < f i := hf i
    have h_ratio : (net.W i j * f j) / (net.W j i * f i) =
        (net.W i j / net.W j i) * (f j / f i) := by ring
    rw [h_ratio,
      Real.log_mul (div_pos hwij hwji).ne' (div_pos hfj hfi).ne',
      Real.log_div hfj.ne' hfi.ne']
    ring

/-- Exact Schnakenberg decomposition `Σ_dot(f) = dS_sys/dt + S_med_dot(f)`. -/
lemma schnakenberg_decomposition_identity {n : ℕ} (net : ConformationalNetwork n)
    (f : Fin n → ℝ) (hf : ∀ i, 0 < f i) :
    entropyProductionRate net f = systemEntropyRate net f + mediumEntropyRate net f := by
  unfold entropyProductionRate mediumEntropyRate
  have h_sum :
      (∑ i, ∑ j, edgeFlux net f i j * Real.log ((net.W i j * f j) / (net.W j i * f i))) =
        (∑ i, ∑ j, edgeFlux net f i j * Real.log (net.W i j / net.W j i)) +
          (∑ i, ∑ j, edgeFlux net f i j * (Real.log (f j) - Real.log (f i))) := by
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl (fun j _ => edgeFlux_mul_log_ratio net f hf i j)
  rw [h_sum, mul_add, half_sum_edgeFlux_sub net f (fun i => Real.log (f i))]
  have h_sys : - ∑ i, masterRHS net f i * Real.log (f i) = systemEntropyRate net f := by
    unfold systemEntropyRate
    have h_one : (∑ i, masterRHS net f i * (Real.log (f i) + 1)) =
        (∑ i, masterRHS net f i * Real.log (f i)) + ∑ i, masterRHS net f i := by
      rw [← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl (fun i _ => by ring)
    rw [h_one, sum_masterRHS_eq_zero net f, add_zero]
  rw [h_sys]
  ring

/-- Universal Second Law of Non-Equilibrium Biophysics: `0 ≤ Σ_dot(f)` for all `f_i > 0`. -/
lemma schnakenberg_second_law {n : ℕ} (net : ConformationalNetwork n)
    (f : Fin n → ℝ) (hf : ∀ i, 0 < f i) :
    0 ≤ entropyProductionRate net f := by
  unfold entropyProductionRate
  refine mul_nonneg (by norm_num) (Finset.sum_nonneg (fun i _ => Finset.sum_nonneg (fun j _ => ?_)))
  by_cases hij : i = j
  · subst hij
    simp [edgeFlux]
  · have hx : 0 < net.W i j * f j := mul_pos (net.off_diag_pos i j hij) (hf j)
    have hy : 0 < net.W j i * f i := mul_pos (net.off_diag_pos j i (Ne.symm hij)) (hf i)
    exact mul_sub_log_div_nonneg hx hy

/--
  @concept Paper D Theorem 2: Exact Schnakenberg Entropy Production Decomposition & Second Law of Non-Equilibrium Biophysics
  @rag_query How does the time derivative of Shannon entropy dS_sys/dt combine with the medium entropy flux S_med_dot to yield Schnakenberg's strictly non-negative entropy production rate Σ_dot(f) = (1/2) ∑_{i,j} (W_{ij} f_j - W_{ji} f_i) ln(W_{ij} f_j / (W_{ji} f_i)) ≥ 0?
  @graph_node StochasticThermodynamics.ST_T2_schnakenberg_entropy_decomposition_and_second_law
  @graph_edge StochasticThermodynamics.ST_T1_master_equation_mass_conservation_and_imex_positivity -> StochasticThermodynamics.ST_T2_schnakenberg_entropy_decomposition_and_second_law
-/
theorem ST_T2_schnakenberg_entropy_decomposition_and_second_law :
    Stmt_ST_T2_SchnakenbergEntropyDecompositionAndSecondLaw :=
  ⟨fun hx hy => mul_sub_log_div_nonneg hx hy,
   fun net f t hpos hderiv => hasDerivAt_systemEntropy net f t hpos hderiv,
   fun net f hpos => schnakenberg_decomposition_identity net f hpos,
   fun net f hpos => schnakenberg_second_law net f hpos⟩

-- ============================================================================
-- SECTION 5: THEOREM 3 — ESPOSITO–VAN DEN BROECK FREE-ENERGY WORK BALANCE
-- ============================================================================

/-- Continuous-time derivative of Kullback-Leibler divergence `D_KL(f(t) ‖ f_eq)`. -/
lemma hasDerivAt_klDivergence {n : ℕ} (net : ConformationalNetwork n)
    (f : ℝ → Fin n → ℝ) (feq : Fin n → ℝ) (t : ℝ)
    (hf_pos : ∀ i, 0 < f t i)
    (hfeq_pos : ∀ i, 0 < feq i)
    (hf_deriv : ∀ i, HasDerivAt (fun s => f s i) (masterRHS net (f t) i) t) :
    HasDerivAt (fun s => klDivergence (f s) feq) (klDivergenceRate net (f t) feq) t := by
  unfold klDivergence klDivergenceRate
  have hterm : ∀ i : Fin n,
      HasDerivAt (fun s => f s i * Real.log (f s i / feq i))
        (masterRHS net (f t) i * (Real.log (f t i / feq i) + 1)) t := by
    intro i
    have hdiv : HasDerivAt (fun s => f s i / feq i) (masterRHS net (f t) i / feq i) t :=
      (hf_deriv i).div_const (feq i)
    have hlog := hdiv.log (div_pos (hf_pos i) (hfeq_pos i)).ne'
    have hmul := (hf_deriv i).mul hlog
    refine hmul.congr_deriv ?_
    have hf_ne : f t i ≠ 0 := (hf_pos i).ne'
    have hfeq_ne : feq i ≠ 0 := (hfeq_pos i).ne'
    field_simp [hf_ne, hfeq_ne]
  exact HasDerivAt.fun_sum (u := Finset.univ)
    (A := fun i s => f s i * Real.log (f s i / feq i))
    (A' := fun i => masterRHS net (f t) i * (Real.log (f t i / feq i) + 1))
    (x := t) (fun i _ => hterm i)

/-- Exact Esposito–Van den Broeck Non-Equilibrium Free-Energy Work Balance under Local Detailed Balance. -/
lemma esposito_van_den_broeck_free_energy_balance {n : ℕ} (net : ConformationalNetwork n)
    (f feq : Fin n → ℝ) (m : Fin n → Fin n → ℝ) (delta_mu kT : ℝ)
    (hkT : 0 < kT)
    (hf : ∀ i, 0 < f i)
    (hfeq : ∀ i, 0 < feq i)
    (hldb : ∀ i j, i ≠ j →
      Real.log (net.W i j / net.W j i) =
        Real.log (feq i / feq j) + (m i j * delta_mu) / kT) :
    kT * klDivergenceRate net f feq =
      chemicalWorkRate net f m delta_mu - kT * entropyProductionRate net f := by
  have h_schnak := schnakenberg_decomposition_identity net f hf
  have h_med_edge : ∀ i j : Fin n,
      edgeFlux net f i j * Real.log (net.W i j / net.W j i) =
        - (edgeFlux net f i j * (Real.log (feq j) - Real.log (feq i))) +
          (1 / kT) * (edgeFlux net f i j * (m i j * delta_mu)) := by
    intro i j
    by_cases hij : i = j
    · subst hij
      simp [edgeFlux]
    · rw [hldb i j hij, Real.log_div (hfeq i).ne' (hfeq j).ne']
      ring
  have h_med :
      mediumEntropyRate net f =
        (∑ i, masterRHS net f i * Real.log (feq i)) +
          (1 / kT) * chemicalWorkRate net f m delta_mu := by
    unfold mediumEntropyRate chemicalWorkRate
    have h_sum :
        (∑ i, ∑ j, edgeFlux net f i j * Real.log (net.W i j / net.W j i)) =
          - (∑ i, ∑ j, edgeFlux net f i j * (Real.log (feq j) - Real.log (feq i))) +
            (1 / kT) * (∑ i, ∑ j, edgeFlux net f i j * (m i j * delta_mu)) := by
      have h1 :
          (∑ i, ∑ j, edgeFlux net f i j * Real.log (net.W i j / net.W j i)) =
            ∑ i, ∑ j, (- (edgeFlux net f i j * (Real.log (feq j) - Real.log (feq i))) +
              (1 / kT) * (edgeFlux net f i j * (m i j * delta_mu))) := by
        refine Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun j _ => h_med_edge i j))
      rw [h1]
      simp only [Finset.sum_add_distrib, Finset.sum_neg_distrib, ← Finset.mul_sum]
    have h_feq_swap := half_sum_edgeFlux_sub net f (fun i => Real.log (feq i))
    rw [h_sum]
    linarith [h_feq_swap]
  have h_kl_split :
      klDivergenceRate net f feq =
        - systemEntropyRate net f - ∑ i, masterRHS net f i * Real.log (feq i) := by
    unfold klDivergenceRate systemEntropyRate
    have h_log : ∀ i : Fin n,
        masterRHS net f i * (Real.log (f i / feq i) + 1) =
          masterRHS net f i * (Real.log (f i) + 1) - masterRHS net f i * Real.log (feq i) := by
      intro i
      rw [Real.log_div (hf i).ne' (hfeq i).ne']
      ring
    have h_sum :
        (∑ i, masterRHS net f i * (Real.log (f i / feq i) + 1)) =
          (∑ i, masterRHS net f i * (Real.log (f i) + 1)) -
            (∑ i, masterRHS net f i * Real.log (feq i)) := by
      rw [← Finset.sum_sub_distrib]
      refine Finset.sum_congr rfl (fun i _ => h_log i)
    linarith [h_sum]
  have h_kl_eq :
      klDivergenceRate net f feq =
        (1 / kT) * chemicalWorkRate net f m delta_mu - entropyProductionRate net f := by
    linarith [h_schnak, h_med, h_kl_split]
  have hkT_ne : kT ≠ 0 := hkT.ne'
  rw [h_kl_eq]
  field_simp [hkT_ne]

/--
  @concept Paper D Theorem 3: Esposito–Van den Broeck Non-Equilibrium Free-Energy Work Balance & NESS Housekeeping Power
  @rag_query How does Local Detailed Balance ln(W_{ij}/W_{ji}) = (E_j - E_i + m_{ij} Δμ_ATP)/(k_B T) yield the exact non-equilibrium free-energy balance k_B T d/dt D_KL(f ‖ f_eq) = W_ATP_dot - k_B T Σ_dot, equilibrium Lyapunov monotonicity at Δμ_ATP = 0, and strictly positive ATP power at NESS?
  @graph_node StochasticThermodynamics.ST_T3_esposito_van_den_broeck_free_energy_balance
  @graph_edge StochasticThermodynamics.ST_T2_schnakenberg_entropy_decomposition_and_second_law -> StochasticThermodynamics.ST_T3_esposito_van_den_broeck_free_energy_balance
-/
theorem ST_T3_esposito_van_den_broeck_free_energy_balance :
    Stmt_ST_T3_EspositoVanDenBroeckFreeEnergyBalance := by
  refine ⟨fun net f feq t hf hfeq hderiv =>
      hasDerivAt_klDivergence net f feq t hf hfeq hderiv,
    fun net f feq m delta_mu kT hkT hf hfeq hldb =>
      esposito_van_den_broeck_free_energy_balance net f feq m delta_mu kT hkT hf hfeq hldb,
    ?_, ?_⟩
  · intro n net f feq kT hkT hf hfeq hdb
    have hldb0 : ∀ i j, i ≠ j →
        Real.log (net.W i j / net.W j i) =
          Real.log (feq i / feq j) + ((0 : Fin n → Fin n → ℝ) i j * 0) / kT := by
      intro i j hij
      rw [hdb i j hij]
      ring
    have hbal := esposito_van_den_broeck_free_energy_balance
      net f feq (0 : Fin n → Fin n → ℝ) 0 kT hkT hf hfeq hldb0
    have hw0 : chemicalWorkRate net f (0 : Fin n → Fin n → ℝ) 0 = 0 := by
      unfold chemicalWorkRate
      simp
    rw [hw0, zero_sub] at hbal
    have hkl : klDivergenceRate net f feq = - entropyProductionRate net f := by
      nlinarith [hbal, hkT]
    have h2nd := schnakenberg_second_law net f hf
    exact ⟨hkl, by linarith⟩
  · intro n net f feq m delta_mu kT hkT hf hfeq hldb hness hpos
    have hbal := esposito_van_den_broeck_free_energy_balance
      net f feq m delta_mu kT hkT hf hfeq hldb
    rw [hness, mul_zero] at hbal
    have hwork : chemicalWorkRate net f m delta_mu = kT * entropyProductionRate net f := by
      linarith
    exact ⟨hwork, by rw [hwork]; exact mul_pos hkT hpos⟩

-- ============================================================================
-- SECTION 6: THEOREM 4 — THERMODYNAMIC UNCERTAINTY RELATION (TUR)
-- ============================================================================

/-- Smooth auxiliary function `ψ(u) = u - 2 + 4 / (e^u + 1)` on `ℝ` for the TUR edge bound. -/
noncomputable def turPsi (u : ℝ) : ℝ :=
  u - 2 + 4 * (Real.exp u + 1)⁻¹

/-- Exact derivative `ψ'(u) = (e^u - 1)^2 / (e^u + 1)^2 ≥ 0` on all of `ℝ`. -/
lemma hasDerivAt_turPsi (u : ℝ) :
    HasDerivAt turPsi ((Real.exp u - 1) ^ 2 / (Real.exp u + 1) ^ 2) u := by
  unfold turPsi
  have hpos : 0 < Real.exp u + 1 := by linarith [Real.exp_pos u]
  have h1 : HasDerivAt (fun x => x - 2) 1 u := (hasDerivAt_id u).sub_const 2
  have h2 : HasDerivAt (fun x => Real.exp x + 1) (Real.exp u) u :=
    (Real.hasDerivAt_exp u).add_const 1
  have h3 : HasDerivAt (fun x => (Real.exp x + 1)⁻¹) (-(Real.exp u) / (Real.exp u + 1) ^ 2) u :=
    h2.fun_inv hpos.ne'
  have h4 : HasDerivAt (fun x => (x - 2) + 4 * (Real.exp x + 1)⁻¹)
      (1 + 4 * (-(Real.exp u) / (Real.exp u + 1) ^ 2)) u :=
    h1.add (h3.const_mul 4)
  refine h4.congr_deriv ?_
  field_simp [hpos.ne']
  ring

/-- Global monotonicity of `turPsi` on `ℝ` via `monotone_of_hasDerivAt_nonneg`. -/
lemma monotone_turPsi : Monotone turPsi := by
  refine monotone_of_hasDerivAt_nonneg (fun u => hasDerivAt_turPsi u) (fun u => ?_)
  positivity

/-- For `t ≥ 1`, `2(t - 1)/(t + 1) ≤ ln(t)`. -/
lemma two_mul_sub_one_div_add_one_le_log {t : ℝ} (ht : 1 ≤ t) :
    2 * (t - 1) / (t + 1) ≤ Real.log t := by
  have ht_pos : 0 < t := lt_of_lt_of_le zero_lt_one ht
  have hu_nonneg : 0 ≤ Real.log t := Real.log_nonneg ht
  have hmono := monotone_turPsi hu_nonneg
  have hpsi0 : turPsi 0 = 0 := by
    unfold turPsi
    rw [Real.exp_zero]
    norm_num
  rw [hpsi0] at hmono
  unfold turPsi at hmono
  rw [Real.exp_log ht_pos] at hmono
  have ht1_pos : 0 < t + 1 := by linarith
  have h_id : 2 - 4 * (t + 1)⁻¹ = 2 * (t - 1) / (t + 1) := by
    field_simp [ht1_pos.ne']
    ring
  linarith [hmono, h_id]

/-- Universal logarithmic-hyperbolic inequality `2(x - y)^2 / (x + y) ≤ (x - y) ln(x / y)` for `x, y > 0`. -/
lemma tur_edge_bound {x y : ℝ} (hx : 0 < x) (hy : 0 < y) :
    2 * (x - y) ^ 2 / (x + y) ≤ (x - y) * Real.log (x / y) := by
  have hxy_pos : 0 < x + y := add_pos hx hy
  have aux : ∀ {a b : ℝ}, 0 < a → 0 < b → b ≤ a →
      2 * (a - b) ^ 2 / (a + b) ≤ (a - b) * Real.log (a / b) := by
    intro a b ha hb hba
    have hab_pos : 0 < a + b := add_pos ha hb
    have h_ratio_ge_one : 1 ≤ a / b := (one_le_div hb).mpr hba
    have h_log := two_mul_sub_one_div_add_one_le_log h_ratio_ge_one
    have h_frac : 2 * (a / b - 1) / (a / b + 1) = 2 * (a - b) / (a + b) := by
      have hb_ne : b ≠ 0 := hb.ne'
      have hab1_ne : a / b + 1 ≠ 0 := by positivity
      field_simp [hb_ne, hab_pos.ne', hab1_ne]
    rw [h_frac] at h_log
    have h_mul := mul_le_mul_of_nonneg_left h_log (sub_nonneg.mpr hba)
    have h_lhs : (a - b) * (2 * (a - b) / (a + b)) = 2 * (a - b) ^ 2 / (a + b) := by ring
    rwa [h_lhs] at h_mul
  rcases le_total y x with hxy | hyx
  · exact aux hx hy hxy
  · have h_sym := aux hy hx hyx
    have h_lhs : 2 * (y - x) ^ 2 / (y + x) = 2 * (x - y) ^ 2 / (x + y) := by
      rw [add_comm y x]
      ring
    have h_rhs : (y - x) * Real.log (y / x) = (x - y) * Real.log (x / y) := by
      have hinv : x / y = (y / x)⁻¹ := (inv_div y x).symm
      rw [hinv, Real.log_inv]
      ring
    rwa [h_lhs, h_rhs] at h_sym

/-- Network Thermodynamic Uncertainty Relation via Cauchy–Schwarz (`Finset.sum_mul_sq_le_sq_mul_sq`). -/
lemma tur_network_cauchy_schwarz_bound {ι : Type*} [Fintype ι] (x y d : ι → ℝ)
    (hx : ∀ e, 0 < x e) (hy : ∀ e, 0 < y e) (hd : ∀ e, |d e| ≤ 1) :
    2 * ((1 / 2 : ℝ) * ∑ e, d e * (x e - y e)) ^ 2 ≤
      ((1 / 2 : ℝ) * ∑ e, (x e + y e)) *
        ((1 / 2 : ℝ) * ∑ e, (x e - y e) * Real.log (x e / y e)) := by
  have h_pos : ∀ e : ι, 0 < x e + y e := fun e => add_pos (hx e) (hy e)
  have h_sqrt_ne : ∀ e : ι, Real.sqrt (x e + y e) ≠ 0 :=
    fun e => (Real.sqrt_pos.mpr (h_pos e)).ne'
  have h_cs := Finset.sum_mul_sq_le_sq_mul_sq (Finset.univ : Finset ι)
    (fun e => d e * Real.sqrt (x e + y e))
    (fun e => (x e - y e) / Real.sqrt (x e + y e))
  have h_lhs_eq :
      (∑ e : ι, (d e * Real.sqrt (x e + y e)) * ((x e - y e) / Real.sqrt (x e + y e))) =
        ∑ e : ι, d e * (x e - y e) := by
    refine Finset.sum_congr rfl (fun e _ => ?_)
    have hne := h_sqrt_ne e
    field_simp [hne]
  have h_act_eq :
      (∑ e : ι, (d e * Real.sqrt (x e + y e)) ^ 2) =
        ∑ e : ι, (d e) ^ 2 * (x e + y e) := by
    refine Finset.sum_congr rfl (fun e _ => ?_)
    rw [mul_pow, Real.sq_sqrt (h_pos e).le]
  have h_flux_eq :
      (∑ e : ι, ((x e - y e) / Real.sqrt (x e + y e)) ^ 2) =
        ∑ e : ι, (x e - y e) ^ 2 / (x e + y e) := by
    refine Finset.sum_congr rfl (fun e _ => ?_)
    rw [div_pow, Real.sq_sqrt (h_pos e).le]
  rw [h_lhs_eq, h_act_eq, h_flux_eq] at h_cs
  have h_d_le : (∑ e : ι, (d e) ^ 2 * (x e + y e)) ≤ ∑ e : ι, (x e + y e) := by
    refine Finset.sum_le_sum (fun e _ => ?_)
    have hd_sq : (d e) ^ 2 ≤ 1 := by
      have h1 : |d e| ^ 2 ≤ 1 ^ 2 := pow_le_pow_left₀ (abs_nonneg (d e)) (hd e) 2
      simpa [sq_abs] using h1
    nlinarith [h_pos e, hd_sq]
  have h_edge_le :
      (∑ e : ι, 2 * (x e - y e) ^ 2 / (x e + y e)) ≤
        ∑ e : ι, (x e - y e) * Real.log (x e / y e) :=
    Finset.sum_le_sum (fun e _ => tur_edge_bound (hx e) (hy e))
  have h_edge_le' :
      (∑ e : ι, (x e - y e) ^ 2 / (x e + y e)) ≤
        (1 / 2 : ℝ) * ∑ e : ι, (x e - y e) * Real.log (x e / y e) := by
    have h_two :
        2 * (∑ e : ι, (x e - y e) ^ 2 / (x e + y e)) =
          ∑ e : ι, 2 * (x e - y e) ^ 2 / (x e + y e) := by
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl (fun e _ => by ring)
    linarith [h_edge_le, h_two]
  have h_nonneg1 : 0 ≤ ∑ e : ι, (x e - y e) ^ 2 / (x e + y e) :=
    Finset.sum_nonneg (fun e _ => div_nonneg (sq_nonneg _) (h_pos e).le)
  have h_nonneg2 : 0 ≤ ∑ e : ι, (x e + y e) :=
    Finset.sum_nonneg (fun e _ => (h_pos e).le)
  have h_prod := mul_le_mul h_d_le h_edge_le' h_nonneg1 h_nonneg2
  have h_comb := le_trans h_cs h_prod
  linarith [h_comb]

/--
  @concept Paper D Theorem 4: Universal Logarithmic-Hyperbolic Edge Bound & Finite-State Thermodynamic Uncertainty Relation (TUR)
  @rag_query Why does 2(x - y)^2 / (x + y) ≤ (x - y) ln(x / y) hold for all x, y > 0, and how does Cauchy–Schwarz lift it to the global Thermodynamic Uncertainty Relation 2 J_net^2 ≤ A_tot Σ_dot bounding molecular-motor precision by entropy production?
  @graph_node StochasticThermodynamics.ST_T4_thermodynamic_uncertainty_relation_and_activity_bound
  @graph_edge StochasticThermodynamics.ST_T3_esposito_van_den_broeck_free_energy_balance -> StochasticThermodynamics.ST_T4_thermodynamic_uncertainty_relation_and_activity_bound
-/
theorem ST_T4_thermodynamic_uncertainty_relation_and_activity_bound :
    Stmt_ST_T4_ThermodynamicUncertaintyRelationAndActivityBound := by
  refine ⟨fun hx hy => tur_edge_bound hx hy, ?_,
    fun x y d hx hy hd => tur_network_cauchy_schwarz_bound x y d hx hy hd⟩
  intro n net f hf i j hij
  have hx : 0 < net.W i j * f j := mul_pos (net.off_diag_pos i j hij) (hf j)
  have hy : 0 < net.W j i * f i := mul_pos (net.off_diag_pos j i (Ne.symm hij)) (hf i)
  exact tur_edge_bound hx hy

-- ============================================================================
-- SECTION 7: MULTI-KERNEL AXIOM AUDIT (`ARGOS-CENTAUR-003`)
-- ============================================================================

#print axioms ST_T1_master_equation_mass_conservation_and_imex_positivity
#print axioms ST_T2_schnakenberg_entropy_decomposition_and_second_law
#print axioms ST_T3_esposito_van_den_broeck_free_energy_balance
#print axioms ST_T4_thermodynamic_uncertainty_relation_and_activity_bound

end ANSE.StochasticThermodynamics

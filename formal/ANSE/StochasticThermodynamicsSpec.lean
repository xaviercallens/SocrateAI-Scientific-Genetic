/-
  ANSE.StochasticThermodynamicsSpec — Formal Specification of Non-Equilibrium Stochastic
  Thermodynamics & ATP-Driven Conformational Master Equations in Lean 4
  ======================================================================================
  Covers the four foundational theorems of non-equilibrium biophysics and molecular-motor
  stochastic thermodynamics (Schnakenberg 1976, Seifert 2005/2012, Esposito & Van den Broeck 2010,
  Barato & Seifert 2015, Gingrich et al. 2016):

    - Theorem 1 (`Stmt_ST_T1_MasterEquationMassConservationAndIMEXPositivity`):
      Continuous-time Pauli Master Equation flux decomposition, universal probability simplex
      mass conservation `∑_i (W f)_i = 0` (`d/dt ∑_i f_i(t) = 0`), unconditional mass conservation
      and strict positivity `f_i^{(k+1)} > 0` of Backward-Euler / L-stable IMEX time integration
      `(I - Δt W) f^{(k+1)} = f^{(k)}` for all `Δt > 0` and all ATP chemical drives `Δμ_ATP ∈ ℝ`,
      contrasted with Forward-Euler positivity breakdown `f_i^{(k+1)} < 0 ↔ Δt > 1 / Γ_i(Δμ_ATP)`.

    - Theorem 2 (`Stmt_ST_T2_SchnakenbergEntropyDecompositionAndSecondLaw`):
      Exact continuous-time Shannon entropy derivative `dS_sys/dt = -∑_i (W f)_i (ln f_i + 1)`
      via `HasDerivAt`, exact Schnakenberg entropy production decomposition
      `Σ_dot(f) = dS_sys/dt + S_med_dot(f)`, and the universal Second Law of Non-Equilibrium
      Biophysics `0 ≤ Σ_dot(f)` for all `n ∈ ℕ` and all positive distributions `f_i > 0`.

    - Theorem 3 (`Stmt_ST_T3_EspositoVanDenBroeckFreeEnergyBalance`):
      Exact continuous-time Kullback-Leibler relative-entropy derivative `d/dt D_KL(f(t) ‖ f_eq)`
      via `HasDerivAt`, and under Local Detailed Balance
      `ln(W_{ij} / W_{ji}) = ln(f_{eq,i} / f_{eq,j}) + m_{ij} Δμ_ATP / (k_B T)`, the exact
      Esposito-Van den Broeck non-equilibrium free-energy work balance
      `k_B T (d/dt D_KL(f ‖ f_eq)) = W_ATP_dot(f) - k_B T Σ_dot(f)`, recovering equilibrium
      Lyapunov decay `d/dt D_KL(f ‖ f_eq) = -Σ_dot(f) ≤ 0` at `Δμ_ATP = 0` and positive
      housekeeping ATP power `W_ATP_dot(f*) = k_B T Σ_dot(f*) > 0` at any non-equilibrium
      steady state (NESS).

    - Theorem 4 (`Stmt_ST_T4_ThermodynamicUncertaintyRelationAndActivityBound`):
      Universal logarithmic-hyperbolic edge inequality `2 (x - y)^2 / (x + y) ≤ (x - y) ln(x / y)`
      for all `x, y > 0` (proved via `monotone_of_hasDerivAt_nonneg` on `ψ(u) = u - 2 + 4/(e^u + 1)`),
      per-edge conformational activity-entropy bound, and Cauchy-Schwarz lifting
      (`Finset.sum_mul_sq_le_sq_mul_sq`) to the finite-state Thermodynamic Uncertainty Relation (TUR)
      `2 J_net^2 ≤ A_tot Σ_dot`.

  Verification Contract:
    - 0 `sorryAx`
    - 0 custom axioms (`axiom`)
    - Universal quantification over `n : ℕ`, `Fin n → ℝ`, and `ℝ`
-/
import Mathlib.Tactic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.Analysis.Calculus.Deriv.Basic
import Mathlib.Analysis.Calculus.Deriv.Inv
import Mathlib.Analysis.Calculus.Deriv.MeanValue

open scoped BigOperators

noncomputable section

namespace ANSE.StochasticThermodynamicsSpec

-- ============================================================================
-- SECTION 1: CORE CONTINUOUS-TIME MARKOV CONFORMATIONAL NETWORK DEFINITIONS
-- ============================================================================

/--
  A continuous-time Markov conformational network on `n` discrete macromolecular microstates.
  `W i j` is the transition rate from microstate `j` to microstate `i`:
    - `off_diag_pos`: strictly positive off-diagonal transition rates `W_{ij} > 0` (`i ≠ j`)
      ensuring microscopic reversibility on connected transitions.
    - `col_sum_zero`: column sums vanish `∑_i W_{ij} = 0`, so `W_{jj} = -∑_{i ≠ j} W_{ij}`.
-/
structure ConformationalNetwork (n : ℕ) where
  W : Fin n → Fin n → ℝ
  off_diag_pos : ∀ i j, i ≠ j → 0 < W i j
  col_sum_zero : ∀ j, ∑ i, W i j = 0

/-- Pauli master equation generator `(W f)_i = ∑_j W_{ij} f_j`. -/
def masterRHS {n : ℕ} (net : ConformationalNetwork n) (f : Fin n → ℝ) (i : Fin n) : ℝ :=
  ∑ j, net.W i j * f j

/-- Total probability mass `∑_i f_i` on the conformational state space `Fin n`. -/
def probMass {n : ℕ} (f : Fin n → ℝ) : ℝ :=
  ∑ i, f i

/-- Probability transition flux `J_{ij}(f) = W_{ij} f_j - W_{ji} f_i` along edge `j → i`. -/
def edgeFlux {n : ℕ} (net : ConformationalNetwork n) (f : Fin n → ℝ) (i j : Fin n) : ℝ :=
  net.W i j * f j - net.W j i * f i

/-- Symmetric dynamical traffic (edge activity) `a_{ij}(f) = W_{ij} f_j + W_{ji} f_i`. -/
def edgeActivity {n : ℕ} (net : ConformationalNetwork n) (f : Fin n → ℝ) (i j : Fin n) : ℝ :=
  net.W i j * f j + net.W j i * f i

/-- Total network dynamical activity `A_tot(f) = (1/2) ∑_i ∑_j (W_{ij} f_j + W_{ji} f_i)`. -/
def totalActivity {n : ℕ} (net : ConformationalNetwork n) (f : Fin n → ℝ) : ℝ :=
  (1 / 2 : ℝ) * ∑ i, ∑ j, edgeActivity net f i j

/-- Generalized antisymmetric conformational current `J_d(f) = (1/2) ∑_i ∑_j d_{ij} J_{ij}(f)`. -/
def generalizedCurrent {n : ℕ} (net : ConformationalNetwork n) (f : Fin n → ℝ)
    (d : Fin n → Fin n → ℝ) : ℝ :=
  (1 / 2 : ℝ) * ∑ i, ∑ j, d i j * edgeFlux net f i j

/-- System Shannon entropy `S_sys(f) = - ∑_i f_i ln(f_i)`. -/
noncomputable def systemEntropy {n : ℕ} (f : Fin n → ℝ) : ℝ :=
  - ∑ i, f i * Real.log (f i)

/-- Instantaneous time derivative of system Shannon entropy along the Pauli master equation:
    `dS_sys/dt = - ∑_i (W f)_i (ln(f_i) + 1)`. -/
noncomputable def systemEntropyRate {n : ℕ} (net : ConformationalNetwork n) (f : Fin n → ℝ) : ℝ :=
  - ∑ i, masterRHS net f i * (Real.log (f i) + 1)

/-- Medium entropy flux rate exported to the thermal/chemical reservoir:
    `S_med_dot(f) = (1/2) ∑_i ∑_j J_{ij}(f) ln(W_{ij} / W_{ji})`. -/
noncomputable def mediumEntropyRate {n : ℕ} (net : ConformationalNetwork n) (f : Fin n → ℝ) : ℝ :=
  (1 / 2 : ℝ) * ∑ i, ∑ j, edgeFlux net f i j * Real.log (net.W i j / net.W j i)

/-- Schnakenberg total entropy production rate:
    `Σ_dot(f) = (1/2) ∑_i ∑_j (W_{ij} f_j - W_{ji} f_i) ln((W_{ij} f_j) / (W_{ji} f_i))`. -/
noncomputable def entropyProductionRate {n : ℕ} (net : ConformationalNetwork n) (f : Fin n → ℝ) : ℝ :=
  (1 / 2 : ℝ) * ∑ i, ∑ j,
    edgeFlux net f i j * Real.log ((net.W i j * f j) / (net.W j i * f i))

/-- Kullback-Leibler divergence (relative entropy) `D_KL(f ‖ f_eq) = ∑_i f_i ln(f_i / f_{eq,i})`. -/
noncomputable def klDivergence {n : ℕ} (f feq : Fin n → ℝ) : ℝ :=
  ∑ i, f i * Real.log (f i / feq i)

/-- Instantaneous time derivative of Kullback-Leibler divergence along the Pauli master equation:
    `d/dt D_KL(f(t) ‖ f_eq) = ∑_i (W f)_i (ln(f_i / f_{eq,i}) + 1)`. -/
noncomputable def klDivergenceRate {n : ℕ} (net : ConformationalNetwork n)
    (f feq : Fin n → ℝ) : ℝ :=
  ∑ i, masterRHS net f i * (Real.log (f i / feq i) + 1)

/-- Chemical work rate (ATP hydrolysis power) injected into the conformational network:
    `W_ATP_dot(f) = (1/2) ∑_i ∑_j J_{ij}(f) (m_{ij} Δμ_ATP)`. -/
def chemicalWorkRate {n : ℕ} (net : ConformationalNetwork n) (f : Fin n → ℝ)
    (m : Fin n → Fin n → ℝ) (delta_mu : ℝ) : ℝ :=
  (1 / 2 : ℝ) * ∑ i, ∑ j, edgeFlux net f i j * (m i j * delta_mu)

-- ============================================================================
-- SECTION 2: MANDATORY INHABITATION PACKAGE FOR 3-STATE ATP MOTOR CYCLE
-- ============================================================================

/--
  Constructively inhabited 3-state ATP-driven molecular motor / chromatin-remodeler
  conformational cycle package (`ARGOS-CENTAUR-002`):
  Certifies a concrete 3-state `ConformationalNetwork 3` driven out of equilibrium by
  ATP hydrolysis (forward rotary rates `k_plus = 4`, backward rates `k_minus = 1`),
  its exact non-equilibrium steady state `f_ness = (1/3, 1/3, 1/3)` with `W f_ness = 0`,
  strictly positive steady-state cyclic current `J_cycle = 1 > 0`, and the exact contrast
  between Forward-Euler negativity (`1 - dt * 5 = -1/4 < 0` at `dt = 1/4`) and
  Backward-Euler positivity (`1 / (1 + dt * 5) = 4/9 > 0`).
-/
structure InhabitedATPConformationalCyclePackage where
  net : ConformationalNetwork 3
  f_ness : Fin 3 → ℝ
  h_ness_pos : ∀ i, 0 < f_ness i
  h_ness_mass : probMass f_ness = 1
  h_ness_stationary : ∀ i, masterRHS net f_ness i = 0
  h_cycle_flux_pos : 0 < edgeFlux net f_ness 1 0
  h_fe_neg_at_quarter : (1 : ℝ) + (1 / 4 : ℝ) * (-5 : ℝ) = -1 / 4 ∧ (-1 / 4 : ℝ) < 0
  h_be_pos_at_quarter : (1 : ℝ) / (1 + (1 / 4 : ℝ) * 5) = 4 / 9 ∧ (0 : ℝ) < 4 / 9

-- ============================================================================
-- SECTION 3: THE FOUR CANONICAL STOCHASTIC THERMODYNAMICS THEOREM SPECS
-- ============================================================================

/--
  @concept Paper D Theorem 1: Pauli Master Equation Mass Conservation, Unconditional IMEX Positivity & Forward-Euler Breakdown
  @rag_query Why does the Pauli Master Equation conserve probability mass ∑_i (W f)_i = 0 for all n, why does Backward-Euler / L-stable IMEX preserve strict positivity f_i^{(k+1)} > 0 for every Δt > 0 and every ATP drive Δμ_ATP, and why does Forward Euler become negative when Δt > 1 / Γ_i(Δμ_ATP)?
  @graph_node StochasticThermodynamicsSpec.Stmt_ST_T1_MasterEquationMassConservationAndIMEXPositivity
-/
def Stmt_ST_T1_MasterEquationMassConservationAndIMEXPositivity : Prop :=
  Nonempty InhabitedATPConformationalCyclePackage ∧
  (∀ {n : ℕ} (net : ConformationalNetwork n) (f : Fin n → ℝ) (i : Fin n),
    masterRHS net f i = ∑ j, edgeFlux net f i j) ∧
  (∀ {n : ℕ} (net : ConformationalNetwork n) (f : Fin n → ℝ),
    ∑ i, masterRHS net f i = 0) ∧
  (∀ {n : ℕ} (net : ConformationalNetwork n) (f : ℝ → Fin n → ℝ) (t : ℝ),
    (∀ i, HasDerivAt (fun s => f s i) (masterRHS net (f t) i) t) →
    HasDerivAt (fun s => probMass (f s)) 0 t) ∧
  (∀ {n : ℕ} (net : ConformationalNetwork n) (f_curr f_next : Fin n → ℝ) (dt : ℝ),
    (∀ i, f_next i - dt * masterRHS net f_next i = f_curr i) →
    probMass f_next = probMass f_curr) ∧
  (∀ {n : ℕ} (net : ConformationalNetwork n) (pi f_curr f_next : Fin n → ℝ)
      (dt : ℝ) (i_min : Fin n),
    0 < dt →
    (∀ i, 0 < pi i) →
    (∀ i, masterRHS net pi i = 0) →
    (∀ i, 0 < f_curr i) →
    (∀ i, f_next i - dt * masterRHS net f_next i = f_curr i) →
    (∀ j, (f_next i_min / pi i_min) * pi j ≤ f_next j) →
    0 < f_next i_min ∧ ∀ j, 0 < f_next j) ∧
  (∀ (w_exit dt : ℝ), 0 < w_exit → 0 < dt →
    ((1 - dt * w_exit < 0 ↔ 1 / w_exit < dt) ∧ 0 < 1 / (1 + dt * w_exit)))

/--
  @concept Paper D Theorem 2: Exact Schnakenberg Entropy Production Decomposition & Second Law of Non-Equilibrium Biophysics
  @rag_query How does the time derivative of Shannon entropy dS_sys/dt combine with the medium entropy flux S_med_dot to yield Schnakenberg's strictly non-negative entropy production rate Σ_dot(f) = (1/2) ∑_{i,j} (W_{ij} f_j - W_{ji} f_i) ln(W_{ij} f_j / (W_{ji} f_i)) ≥ 0?
  @graph_node StochasticThermodynamicsSpec.Stmt_ST_T2_SchnakenbergEntropyDecompositionAndSecondLaw
-/
def Stmt_ST_T2_SchnakenbergEntropyDecompositionAndSecondLaw : Prop :=
  (∀ {x y : ℝ}, 0 < x → 0 < y → 0 ≤ (x - y) * Real.log (x / y)) ∧
  (∀ {n : ℕ} (net : ConformationalNetwork n) (f : ℝ → Fin n → ℝ) (t : ℝ),
    (∀ i, 0 < f t i) →
    (∀ i, HasDerivAt (fun s => f s i) (masterRHS net (f t) i) t) →
    HasDerivAt (fun s => systemEntropy (f s)) (systemEntropyRate net (f t)) t) ∧
  (∀ {n : ℕ} (net : ConformationalNetwork n) (f : Fin n → ℝ),
    (∀ i, 0 < f i) →
    entropyProductionRate net f = systemEntropyRate net f + mediumEntropyRate net f) ∧
  (∀ {n : ℕ} (net : ConformationalNetwork n) (f : Fin n → ℝ),
    (∀ i, 0 < f i) →
    0 ≤ entropyProductionRate net f)

/--
  @concept Paper D Theorem 3: Esposito–Van den Broeck Non-Equilibrium Free-Energy Work Balance & NESS Housekeeping Power
  @rag_query How does Local Detailed Balance ln(W_{ij}/W_{ji}) = (E_j - E_i + m_{ij} Δμ_ATP)/(k_B T) yield the exact non-equilibrium free-energy balance k_B T d/dt D_KL(f ‖ f_eq) = W_ATP_dot - k_B T Σ_dot, equilibrium Lyapunov monotonicity at Δμ_ATP = 0, and strictly positive ATP power at NESS?
  @graph_node StochasticThermodynamicsSpec.Stmt_ST_T3_EspositoVanDenBroeckFreeEnergyBalance
-/
def Stmt_ST_T3_EspositoVanDenBroeckFreeEnergyBalance : Prop :=
  (∀ {n : ℕ} (net : ConformationalNetwork n) (f : ℝ → Fin n → ℝ) (feq : Fin n → ℝ) (t : ℝ),
    (∀ i, 0 < f t i) →
    (∀ i, 0 < feq i) →
    (∀ i, HasDerivAt (fun s => f s i) (masterRHS net (f t) i) t) →
    HasDerivAt (fun s => klDivergence (f s) feq) (klDivergenceRate net (f t) feq) t) ∧
  (∀ {n : ℕ} (net : ConformationalNetwork n) (f feq : Fin n → ℝ)
      (m : Fin n → Fin n → ℝ) (delta_mu kT : ℝ),
    0 < kT →
    (∀ i, 0 < f i) →
    (∀ i, 0 < feq i) →
    (∀ i j, i ≠ j →
      Real.log (net.W i j / net.W j i) =
        Real.log (feq i / feq j) + (m i j * delta_mu) / kT) →
    kT * klDivergenceRate net f feq =
      chemicalWorkRate net f m delta_mu - kT * entropyProductionRate net f) ∧
  (∀ {n : ℕ} (net : ConformationalNetwork n) (f feq : Fin n → ℝ) (kT : ℝ),
    0 < kT →
    (∀ i, 0 < f i) →
    (∀ i, 0 < feq i) →
    (∀ i j, i ≠ j → Real.log (net.W i j / net.W j i) = Real.log (feq i / feq j)) →
    klDivergenceRate net f feq = - entropyProductionRate net f ∧
    klDivergenceRate net f feq ≤ 0) ∧
  (∀ {n : ℕ} (net : ConformationalNetwork n) (f feq : Fin n → ℝ)
      (m : Fin n → Fin n → ℝ) (delta_mu kT : ℝ),
    0 < kT →
    (∀ i, 0 < f i) →
    (∀ i, 0 < feq i) →
    (∀ i j, i ≠ j →
      Real.log (net.W i j / net.W j i) =
        Real.log (feq i / feq j) + (m i j * delta_mu) / kT) →
    klDivergenceRate net f feq = 0 →
    0 < entropyProductionRate net f →
    chemicalWorkRate net f m delta_mu = kT * entropyProductionRate net f ∧
    0 < chemicalWorkRate net f m delta_mu)

/--
  @concept Paper D Theorem 4: Universal Logarithmic-Hyperbolic Edge Bound & Finite-State Thermodynamic Uncertainty Relation (TUR)
  @rag_query Why does 2(x - y)^2 / (x + y) ≤ (x - y) ln(x / y) hold for all x, y > 0, and how does Cauchy–Schwarz lift it to the global Thermodynamic Uncertainty Relation 2 J_net^2 ≤ A_tot Σ_dot bounding molecular-motor precision by entropy production?
  @graph_node StochasticThermodynamicsSpec.Stmt_ST_T4_ThermodynamicUncertaintyRelationAndActivityBound
-/
def Stmt_ST_T4_ThermodynamicUncertaintyRelationAndActivityBound : Prop :=
  (∀ {x y : ℝ}, 0 < x → 0 < y →
    2 * (x - y) ^ 2 / (x + y) ≤ (x - y) * Real.log (x / y)) ∧
  (∀ {n : ℕ} (net : ConformationalNetwork n) (f : Fin n → ℝ),
    (∀ i, 0 < f i) →
    ∀ i j, i ≠ j →
      2 * (edgeFlux net f i j) ^ 2 / edgeActivity net f i j ≤
        edgeFlux net f i j * Real.log ((net.W i j * f j) / (net.W j i * f i))) ∧
  (∀ {m : ℕ} (x y d : Fin m → ℝ),
    (∀ e, 0 < x e) →
    (∀ e, 0 < y e) →
    (∀ e, |d e| ≤ 1) →
    2 * ((1 / 2 : ℝ) * ∑ e, d e * (x e - y e)) ^ 2 ≤
      ((1 / 2 : ℝ) * ∑ e, (x e + y e)) *
        ((1 / 2 : ℝ) * ∑ e, (x e - y e) * Real.log (x e / y e)))

-- ============================================================================
-- SECTION 4: DEEP FLUCTUATION THEOREMS & CHROMATIN POLYMER BIOPHYSICS (`ST.T5` & `ST.T6`)
-- ============================================================================

/-- McGhee–von Hippel 1D hard-rod steric accessibility factor for nucleosomes of footprint
    `l_nuc` base pairs (`l_nuc = 147 bp` for `PDB: 1KX5`) sliding along DNA at density `rho`:
    `g_MvH(l_nuc, rho) = (1 - l_nuc * rho) / (1 - (l_nuc - 1) * rho)`. -/
def mcGheeVonHippelSteric (lNuc rho : ℝ) : ℝ :=
  (1 - lNuc * rho) / (1 - (lNuc - 1) * rho)

/-- 1D Tonks–Girardeau hard-rod dimensionless osmotic pressure `Π / (k_B T) = rho / (1 - l_nuc * rho)`. -/
def tonksGirardeauOsmoticPressure (lNuc rho : ℝ) : ℝ :=
  rho / (1 - lNuc * rho)

/--
  @concept Paper D Theorem 5: Crooks Fluctuation Theorem, Jarzynski Equality, Second-Law Work Dissipation & Schnakenberg Cycle Affinity
  @rag_query How does the Crooks fluctuation relation P_F(γ) = P_R(γ) exp(β(W(γ) - ΔF)) imply the Jarzynski Equality ∑_γ P_F(γ) exp(-β W(γ)) = exp(-β ΔF), why does convexity of exp imply the Second Law ⟨W⟩ ≥ ΔF, and how does Schnakenberg cycle entropy production reduce to J_cyc * Δμ_ATP / (k_B T)?
  @graph_node StochasticThermodynamicsSpec.Stmt_ST_T5_JarzynskiEqualityCrooksFluctuationAndSchnakenbergCycle
-/
def Stmt_ST_T5_JarzynskiEqualityCrooksFluctuationAndSchnakenbergCycle : Prop :=
  (∀ {K : ℕ} (pF pR W : Fin K → ℝ) (beta dF : ℝ),
    (∀ γ, 0 < pR γ) →
    (∀ γ, pF γ = pR γ * Real.exp (beta * (W γ - dF))) →
    (∀ γ, Real.log (pF γ / pR γ) = beta * (W γ - dF)) ∧
    (∑ γ, pR γ = 1 →
      ∑ γ, pF γ * Real.exp (-beta * W γ) = Real.exp (-beta * dF))) ∧
  (∀ {K : ℕ} (pF W : Fin K → ℝ) (beta dF : ℝ),
    0 < beta →
    (∀ γ, 0 ≤ pF γ) →
    ∑ γ, pF γ = 1 →
    ∑ γ, pF γ * Real.exp (-beta * W γ) = Real.exp (-beta * dF) →
    dF ≤ ∑ γ, pF γ * W γ ∧ 0 ≤ (∑ γ, pF γ * W γ) - dF) ∧
  (∀ {m : ℕ} (wPlus wMinus : Fin m → ℝ) (jCyc deltaMu kT : ℝ),
    0 < kT →
    0 ≤ jCyc →
    0 ≤ deltaMu →
    (∑ e, Real.log (wPlus e / wMinus e) = deltaMu / kT) →
    (∑ e, jCyc * Real.log (wPlus e / wMinus e) = jCyc * (deltaMu / kT)) ∧
    0 ≤ ∑ e, jCyc * Real.log (wPlus e / wMinus e))

/--
  @concept Paper D Theorem 6: McGhee–von Hippel / Tonks–Girardeau Nucleosome Lattice Gas, Rouse Chromatin Polymer Dissipation & smFRET Observables
  @rag_query How does the 1D McGhee–von Hippel / Tonks–Girardeau nucleosome lattice gas govern steric slowing g_MvH(147, ρ) and osmotic pressure Π_TG(147, ρ) between euchromatin (200 bp) and heterochromatin (165 bp), and why does a 3D Rouse chromatin bead-spring polymer monotonically dissipate elastic free energy?
  @graph_node StochasticThermodynamicsSpec.Stmt_ST_T6_McGheeVonHippelTonksGirardeauChromatinLatticeAndRousePolymer
-/
def Stmt_ST_T6_McGheeVonHippelTonksGirardeauChromatinLatticeAndRousePolymer : Prop :=
  (∀ (lNuc rho : ℝ), 1 < lNuc → 0 < rho → rho < 1 / lNuc →
    0 < mcGheeVonHippelSteric lNuc rho ∧
    mcGheeVonHippelSteric lNuc rho < 1 ∧
    mcGheeVonHippelSteric lNuc (1 / lNuc) = 0 ∧
    HasDerivAt (mcGheeVonHippelSteric lNuc)
      (-1 / (1 - (lNuc - 1) * rho) ^ 2) rho ∧
    tonksGirardeauOsmoticPressure lNuc rho - rho =
      (lNuc * rho ^ 2) / (1 - lNuc * rho) ∧
    0 < tonksGirardeauOsmoticPressure lNuc rho - rho) ∧
  (∀ {N : ℕ} (U : ℝ → ℝ) (G : Fin N → ℝ) (zeta t : ℝ),
    0 < zeta →
    HasDerivAt U (∑ i, G i * (-(1 / zeta) * G i)) t →
    (∑ i, G i * (-(1 / zeta) * G i)) = -(1 / zeta) * ∑ i, (G i) ^ 2 ∧
    (∑ i, G i * (-(1 / zeta) * G i)) ≤ 0) ∧
  (mcGheeVonHippelSteric 147 (1 / 200) = 53 / 54 ∧
   mcGheeVonHippelSteric 147 (1 / 165) = 18 / 19 ∧
   mcGheeVonHippelSteric 147 (1 / 165) < mcGheeVonHippelSteric 147 (1 / 200) ∧
   tonksGirardeauOsmoticPressure 147 (1 / 165) /
     tonksGirardeauOsmoticPressure 147 (1 / 200) = 53 / 18 ∧
   (2 : ℝ) < 53 / 18)

end ANSE.StochasticThermodynamicsSpec

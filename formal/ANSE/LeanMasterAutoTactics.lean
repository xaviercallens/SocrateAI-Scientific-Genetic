/-
  ANSE.LeanMasterAutoTactics — Domain-Specific Tactic Combinators
  ================================================================
  Integrates SocrateAI-Scientific-Agora-LeanMaster's AutoTactics combinators
  into the ANSE / Elenchus formal verification pipeline.

  Tactics provided:
    - `auto_bound`       : Discharges positivity, boundedness, and inequality goals
    - `auto_eq`          : Discharges algebraic and ring equalities
    - `prove_pos`        : Proves `0 < expr` via `positivity`, `linarith`, `nlinarith`
    - `prove_bounded`    : Proves boundedness obligations via `linarith` / `nlinarith`
    - `conv_rate_lt_one` : Proves Krylov/FGMRES contraction rate `< 1`
    - `geo_decrease`     : Proves geometric residual decrease under contraction
    - `energy_monotone`  : Proves thermodynamic / Hamiltonian energy monotonicity
-/
import Mathlib.Tactic

namespace ANSE.LeanMasterAutoTactics

/--
  @concept LeanMaster Automated Inequality & Positivity Tactic
  @rag_query How does `auto_bound` discharge inequality and positivity subgoals in Lean 4?
  @graph_node LeanMasterAutoTactics.auto_bound
  Attempts in sequence: `positivity`, `norm_num`, `linarith`, `nlinarith`, `omega`, `ring_nf`.
-/
macro "auto_bound" : tactic =>
  `(tactic|
    first
    | positivity
    | norm_num
    | linarith
    | nlinarith
    | omega
    | (ring_nf; try linarith)
    | assumption)

/--
  @concept LeanMaster Automated Algebraic Equality Tactic
  @rag_query How does `auto_eq` discharge algebraic and rational identities in Lean 4?
  @graph_node LeanMasterAutoTactics.auto_eq
  Attempts in sequence: `ring`, `norm_num`, `linarith`, `omega`, `simp; ring`.
-/
macro "auto_eq" : tactic =>
  `(tactic|
    first
    | rfl
    | ring
    | norm_num
    | linarith
    | omega
    | (field_simp; ring)
    | (simp; try ring)
    | assumption)

/--
  @concept LeanMaster Strict Positivity Tactic
  @rag_query How does `prove_pos` certify strict positivity of physical parameters?
  @graph_node LeanMasterAutoTactics.prove_pos
-/
macro "prove_pos" : tactic =>
  `(tactic|
    first
    | positivity
    | linarith
    | nlinarith
    | norm_num
    | assumption)

/--
  @concept LeanMaster Boundedness Tactic
  @rag_query How does `prove_bounded` verify upper and lower bounds in numerical solvers?
  @graph_node LeanMasterAutoTactics.prove_bounded
-/
macro "prove_bounded" : tactic =>
  `(tactic|
    first
    | linarith
    | nlinarith
    | norm_num
    | positivity
    | assumption)

/--
  @concept LeanMaster Krylov Contraction Rate Tactic
  @rag_query How does `conv_rate_lt_one` verify that a contraction factor is strictly less than 1?
  @graph_node LeanMasterAutoTactics.conv_rate_lt_one
-/
macro "conv_rate_lt_one" : tactic =>
  `(tactic|
    first
    | nlinarith
    | linarith
    | norm_num
    | positivity)

/--
  @concept LeanMaster Geometric Residual Decrease Tactic
  @rag_query How does `geo_decrease` prove monotone residual contraction across iterations?
  @graph_node LeanMasterAutoTactics.geo_decrease
-/
macro "geo_decrease" : tactic =>
  `(tactic|
    first
    | nlinarith
    | linarith
    | positivity)

/--
  @concept LeanMaster Thermodynamic Energy Monotonicity Tactic
  @rag_query How does `energy_monotone` prove energy dissipation or conservation?
  @graph_node LeanMasterAutoTactics.energy_monotone
-/
macro "energy_monotone" : tactic =>
  `(tactic|
    first
    | ring
    | linarith
    | nlinarith
    | (field_simp; ring)
    | norm_num)

end ANSE.LeanMasterAutoTactics

import Programmes.PopGen.MoranModel

/-!
# Data certificates for the rusty-SUNDIALS CVODE absorption runs (P7 paper stage, 2026-10-25)

The decimal literals below are the `p_fix` values printed by `fixation_ode`
(rusty-SUNDIALS `crates/qf-popgen`, branch `physics/p7-fixation-2026-10-25`), i.e. floating-point
CVODE output (Tier X). Each theorem lets the Lean kernel check that the literal lies within `1e-11`
(`1e-12` in the neutral case) of the exact value `moranFix N r i` defined in `MoranModel.lean`.
This certifies an arithmetic inequality between a decimal literal and an exact rational; it does
not certify CVODE, and it does not make the literal a theorem about any Markov chain. The
connection between `P^n(i,N)` and the CVODE trajectory (Poisson mixture) is not formalised.

`ctl_neg_wrong_literal` is the negative control: a literal that is not near the exact value is
refused by the same kernel check. Planted `sorry`/`axiom` controls for this module are in
`FixationDataCertControls.lean` (those for the fixation and Price modules are in
`PopGenProgControls.lean` and `PriceControls.lean`).
-/

namespace ProgPopGenCert

open ProgPopGen

theorem cert_N2_r2_i1 : |(0.6666666666666601 : ℝ) - moranFix 2 2 1| < 1e-11 := by
  rw [abs_lt]
  constructor <;> norm_num [moranFix, geomS, Finset.sum_range_succ]

theorem cert_N3_r2_i1 : |(0.5714285714284759 : ℝ) - moranFix 3 2 1| < 1e-11 := by
  rw [abs_lt]
  constructor <;> norm_num [moranFix, geomS, Finset.sum_range_succ]

theorem cert_N3_r2_i2 : |(0.8571428571427792 : ℝ) - moranFix 3 2 2| < 1e-11 := by
  rw [abs_lt]
  constructor <;> norm_num [moranFix, geomS, Finset.sum_range_succ]

theorem cert_N10_r1p5_i1 : |(0.33921585523416337 : ℝ) - moranFix 10 (3 / 2) 1| < 1e-11 := by
  rw [abs_lt]
  constructor <;> norm_num [moranFix, geomS, Finset.sum_range_succ]

theorem cert_N10_r1p5_i5 : |(0.8836363636355173 : ℝ) - moranFix 10 (3 / 2) 5| < 1e-11 := by
  rw [abs_lt]
  constructor <;> norm_num [moranFix, geomS, Finset.sum_range_succ]

/-- Neutral case: `moranFix 10 1 3 = 3/10` (exactly, `moranFix_neutral`). -/
theorem cert_N10_neutral_i3 : |(0.29999999999999977 : ℝ) - moranFix 10 1 3| < 1e-12 := by
  rw [moranFix_neutral]
  rw [abs_lt]
  constructor <;> norm_num

/-- Negative control: the neutral value `0.1` for `i = 1` is not within `1e-3` of the exact
Moran value at `N = 10`, `r = 3/2` (exact value about `0.3392`). -/
theorem ctl_neg_wrong_literal : ¬ (|(0.1 : ℝ) - moranFix 10 (3 / 2) 1| < 1e-3) := by
  rw [not_lt, abs_sub_comm]
  have h : (0.3392 : ℝ) ≤ moranFix 10 (3 / 2) 1 := by
    norm_num [moranFix, geomS, Finset.sum_range_succ]
  have h' : (0.1 : ℝ) ≤ 0.3392 := by norm_num
  rw [abs_of_nonneg (by linarith)]
  norm_num
  linarith

end ProgPopGenCert

#print axioms ProgPopGenCert.cert_N2_r2_i1
#print axioms ProgPopGenCert.cert_N3_r2_i1
#print axioms ProgPopGenCert.cert_N3_r2_i2
#print axioms ProgPopGenCert.cert_N10_r1p5_i1
#print axioms ProgPopGenCert.cert_N10_r1p5_i5
#print axioms ProgPopGenCert.cert_N10_neutral_i3
#print axioms ProgPopGenCert.ctl_neg_wrong_literal

import Programmes.PopGen.MoranModel

/-!
# Planted controls for `FixationDataCert.lean` (P7 paper stage). NEVER imported by a claimed module.

Both theorems claim a FALSE bound (`0.1` is not within `1e-11` of `moranFix 10 (3/2) 1`, about `0.3392`).
`ctl_sorry` is "proved" with `sorry` and `ctl_axiom` from a planted axiom. The file compiles with exit
status 0; the planted controls are caught only by `#print axioms`, which lists `sorryAx` and `ctl_cheat`.
-/

namespace ProgPopGenCertControls

open ProgPopGen

theorem ctl_sorry : |(0.1 : ℝ) - moranFix 10 (3 / 2) 1| < 1e-11 := by
  sorry

axiom ctl_cheat : |(0.1 : ℝ) - moranFix 10 (3 / 2) 1| < 1e-11

theorem ctl_axiom : |(0.1 : ℝ) - moranFix 10 (3 / 2) 1| < 1e-11 := ctl_cheat

end ProgPopGenCertControls

#print axioms ProgPopGenCertControls.ctl_sorry
#print axioms ProgPopGenCertControls.ctl_axiom

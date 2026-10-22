/-
  Certificates emitted by `reticulation` for the kernel to check.

  Source      : deliberately corrupted: far and near distances swapped
  Label       : poc-rejection-control
  Certificates: 1

  Each `CountCertificate` carries integer site counts. Its two `far_*`
  fields are the cross-multiplied inequalities of
  `ANSE.GeneticaReticulata.ratio_lt_ratio_iff`, so the kernel decides
  them exactly. A certificate that fails to typecheck is a defect in
  the search, not a negative result.

  Generated file -- do not edit by hand.
-/
import Mathlib.Tactic
import ANSE.GeneticaReticulata

open ANSE.GeneticaReticulata

namespace Reticulata.Emitted

/-- clonal_002 -- recombinant_of_clonal_003_and_clonal_013 -- clonal_013; exact margin 244000/16000000. -/
def cert0 : ANSE.GeneticaReticulata.CountCertificate (Fin 3) where
  x := 0
  y := 1
  z := 2
  n_xy := 131
  d_xy := 4000
  n_yz := 70
  d_yz := 4000
  n_xz := 69
  d_xz := 4000
  hd_xy := by norm_num
  hd_yz := by norm_num
  hd_xz := by norm_num
  far_xy := by norm_num
  far_yz := by norm_num

theorem cert0_refutes : 0 < ANSE.GeneticaReticulata.margin cert0 := ANSE.GeneticaReticulata.margin_pos cert0

end Reticulata.Emitted

#print axioms Reticulata.Emitted.cert0_refutes

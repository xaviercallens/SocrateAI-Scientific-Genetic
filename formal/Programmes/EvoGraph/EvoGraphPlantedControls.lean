import Programmes.EvoGraph.Examples

/-!
# Planted controls for N4 (must be CAUGHT by the axiom audit)

Not part of any claimed module; nothing imports this file. Each control states a FALSE claim:

* `ctl_sorry`: the star centre fixes with the Moran value `4/7` (kernel-checked value: `5/12`,
  `EvoGraph.star_fixation_centre`). Must report `sorryAx`.
* `ctl_axiom`: the uniform-random single mutant on the star fixes with `4/7` (kernel-checked value
  `7/12`). Must report `EvoGraphControls.ctl_cheat`.
* `ctl_axiom_lumpable`: the mutant count on the isothermal 4-cycle is lumpable, i.e. `P⁺({0,1}) =
  P⁺({0,2})` (kernel-checked: `1/3 ≠ 2/3`, `EvoGraph.cycle4_not_lumpable`). Must report `ctl_cheat`.
-/

namespace EvoGraphControls

open EvoGraph Finset Filter Topology

theorem ctl_sorry :
    Tendsto (fun t => sPow (bdK (fun _ => 1) Wstar 2) t {0} univ) atTop (𝓝 (4 / 7)) := sorry

axiom ctl_cheat : False

theorem ctl_axiom :
    Tendsto (fun t => (∑ v : Fin 3, sPow (bdK (fun _ => 1) Wstar 2) t {v} univ) / 3) atTop
      (𝓝 (4 / 7)) := ctl_cheat.elim

theorem ctl_axiom_lumpable :
    upP (fun _ => 1) Wcyc 2 {0, 1} = upP (fun _ => 1) Wcyc 2 {0, 2} := ctl_cheat.elim

end EvoGraphControls

#print axioms EvoGraphControls.ctl_sorry
#print axioms EvoGraphControls.ctl_axiom
#print axioms EvoGraphControls.ctl_axiom_lumpable

import Programmes.PopGen.WrightFisherFixation
import Programmes.PopGen.MoranFixation

/-!
# Planted controls for P7 (must be CAUGHT by the axiom audit)

Not part of any claimed module; nothing imports this file. `ctl_sorry` must report `sorryAx` and
`ctl_axiom` must report `ProgPopGenControls.ctl_cheat`. Both claim WRONG limits (the exponent-flipped
value `1/7` and the neutral value `1/3` for `N = 3`, `r = 2`, `i = 1`, where the kernel-checked limit is
`4/7`), so an audit that accepted them would accept false statements.
-/

namespace ProgPopGenControls

open ProgPopGen Filter Topology

theorem ctl_sorry :
    Tendsto (fun n => kPow (moranP 3 2) n 1 (Fin.last 3)) atTop (𝓝 (1 / 7)) := sorry

axiom ctl_cheat : False

theorem ctl_axiom :
    Tendsto (fun n => kPow (moranP 3 2) n 1 (Fin.last 3)) atTop (𝓝 (1 / 3)) := ctl_cheat.elim

end ProgPopGenControls

#print axioms ProgPopGenControls.ctl_sorry
#print axioms ProgPopGenControls.ctl_axiom

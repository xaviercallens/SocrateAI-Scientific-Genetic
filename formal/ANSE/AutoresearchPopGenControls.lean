/-! Planted controls for the PopGen audit (kept out of the claimed modules so importers never inherit them). -/

namespace AutoresearchPopGenControls

theorem ctl_sorry : (1 : Nat) = 2 := sorry

axiom ctl_cheat : False

theorem ctl_axiom : (1 : Nat) = 2 := ctl_cheat.elim

end AutoresearchPopGenControls

#print axioms AutoresearchPopGenControls.ctl_sorry
#print axioms AutoresearchPopGenControls.ctl_axiom

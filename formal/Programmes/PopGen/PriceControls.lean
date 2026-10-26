import Programmes.PopGen.Price

/-!
# Planted controls for the P7 Price / Fisher file (must be CAUGHT by the axiom audit)

Not part of any claimed module; nothing imports this file. `ctl_sorry` must report `sorryAx` and
`ctl_axiom` must report `ProgPopGenPriceControls.ctl_cheat`. Both claim WRONG statements:
`ctl_sorry` claims `Δw̄ = Var(w)` (missing the `1/w̄` factor; kernel-checked value `1/6` vs `1/4`,
see `ProgPopGen.Price.ctl_neg_missing_wbar`), and `ctl_axiom` claims the variance formula under
transmission bias (kernel-checked `-1/2` vs `1/6`, see `ProgPopGen.Price.ctl_neg_transmission_bias`).
-/

namespace ProgPopGenPriceControls

open ProgPopGen.Price

theorem ctl_sorry :
    wbarNext (Finset.univ : Finset (Fin 2)) ![1, 1] ![1, 2] - wmean Finset.univ ![1, 1] ![1, 2]
      = wvar (Finset.univ : Finset (Fin 2)) ![1, 1] (![1, 2] : Fin 2 → ℚ) := sorry

axiom ctl_cheat : False

theorem ctl_axiom :
    zbarNext (Finset.univ : Finset (Fin 2)) ![1, 1] ![1, 2] ![1, 1]
        - wmean Finset.univ ![1, 1] ![1, 2]
      = wvar (Finset.univ : Finset (Fin 2)) ![1, 1] (![1, 2] : Fin 2 → ℚ)
          / wmean Finset.univ ![1, 1] ![1, 2] := ctl_cheat.elim

end ProgPopGenPriceControls

#print axioms ProgPopGenPriceControls.ctl_sorry
#print axioms ProgPopGenPriceControls.ctl_axiom

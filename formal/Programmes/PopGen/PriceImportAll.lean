import Autoresearch.PopGenWrightFisher
import Autoresearch.PopGenFisher
import Programmes.PopGen.WrightFisherFixation
import Programmes.PopGen.MoranFixation
import Programmes.PopGen.Price

/-!
# Import-together check for the P7 Price / Fisher module

The Price module imports together with H11, H12 and the P7 fixation modules without name clashes,
and the general haploid monotonicity is instantiated at a concrete 4-type population (non-vacuity).
-/

namespace ProgPopGenPriceCheck

open ProgPopGen.Price

/-- The general haploid result instantiated at a 4-type real population. -/
theorem instance_four_types :
    wmean (Finset.univ : Finset (Fin 4)) ![1, 1, 1, 1] (![1, 2, 3, 4] : Fin 4 → ℝ)
      ≤ wbarNext Finset.univ ![1, 1, 1, 1] ![1, 2, 3, 4] :=
  fisher_haploid_nondecreasing _ _ _ (by intro i _; fin_cases i <;> norm_num)
    (by simp [Fin.sum_univ_four]; norm_num) (by simp [Fin.sum_univ_four]; norm_num)

end ProgPopGenPriceCheck

#print axioms ProgPopGenPriceCheck.instance_four_types
#print axioms ProgPopGen.Price.price_equation
#print axioms ProgPopGen.wright_fisher_fixation_limit
#print axioms AutoresearchPopGen.fisher_mean_fitness_nondecreasing

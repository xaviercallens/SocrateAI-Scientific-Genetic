import Autoresearch.PopGenWrightFisher
import Autoresearch.PopGenFisher
import Programmes.PopGen.WrightFisherFixation
import Programmes.PopGen.MoranFixation

/-!
# Import-together and non-vacuity check for P7

All claimed population-genetics modules (H11, H12 and the P7 modules) import together without name
clashes, and the chain-level theorems are instantiated at concrete parameters, showing that their
hypotheses are satisfiable (non-vacuity).
-/

namespace ProgPopGenCheck

open ProgPopGen AutoresearchPopGen Filter Topology

/-- Hypotheses of the general theorem hold for a concrete Moran chain. -/
theorem nonvacuous_moran_chain : AbsorbingChain (moranP 3 2) :=
  moran_absorbing (by norm_num) (by norm_num)

/-- Hypotheses of the general theorem hold for a concrete Wright–Fisher chain. -/
theorem nonvacuous_wf_chain : AbsorbingChain (wfP 3) :=
  wf_absorbing (by norm_num)

/-- Wright–Fisher, `N = 3`, `i = 1`: `Pⁿ(1,3) → 1/3`. -/
theorem instance_wf_N3 :
    Tendsto (fun n => kPow (wfP 3) n 1 (Fin.last 3)) atTop (𝓝 (1 / 3)) := by
  have h := (wright_fisher_fixation_limit (N := 3) (by norm_num) 1).1
  have hv : (((1 : Fin 4) : ℕ) : ℝ) / (3 : ℕ) = 1 / 3 := by
    rw [show ((1 : Fin 4) : ℕ) = 1 from rfl]
    norm_num
  rwa [hv] at h

/-- Moran, `N = 10`, `r = 3/2`, `i = 1`: `Pⁿ(1,10) → (1 - 2/3)/(1 - (2/3)^10)`. -/
theorem instance_moran_N10 :
    Tendsto (fun n => kPow (moranP 10 (3 / 2)) n 1 (Fin.last 10)) atTop
      (𝓝 ((1 - (2 / 3 : ℝ)) / (1 - (2 / 3) ^ 10))) := by
  have h := moran_fixation_limit_closed (N := 10) (r := 3 / 2) (by norm_num) (by norm_num)
    (by norm_num) 1
  have hv : (1 - (3 / 2 : ℝ)⁻¹ ^ ((1 : Fin 11) : ℕ)) / (1 - (3 / 2 : ℝ)⁻¹ ^ 10)
      = (1 - (2 / 3 : ℝ)) / (1 - (2 / 3) ^ 10) := by
    rw [show ((1 : Fin 11) : ℕ) = 1 from rfl]
    norm_num
  rwa [hv] at h

end ProgPopGenCheck

#print axioms ProgPopGenCheck.nonvacuous_moran_chain
#print axioms ProgPopGenCheck.nonvacuous_wf_chain
#print axioms ProgPopGenCheck.instance_wf_N3
#print axioms ProgPopGenCheck.instance_moran_N10
#print axioms AutoresearchPopGen.harmonic_unique
#print axioms AutoresearchPopGen.fisher_mean_fitness_nondecreasing

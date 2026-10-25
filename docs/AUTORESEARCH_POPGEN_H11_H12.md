# Autoresearch population genetics (H11, H12): kernel-checked contribution from AutoevolveAI

Contributed from the AutoevolveAI autoresearch round of 2026-10-25. Source: github.com/xaviercallens/AutoevolveAI,
branch `worktree-openai-math-discovery`, `formal_cert/Autoresearch/` and `results/autoresearch/2026-10-25/{H11,H12}`.

Both results formalize **known textbook population genetics**. Novelty is unchecked and they are not presented as
discoveries. Each was compiled by its producer and then independently re-verified by a separate agent (producer ≠
verifier), which accepted both with corrections.

## Files
| File | Content |
|---|---|
| `formal/ANSE/AutoresearchPopGenFisher.lean` | H12: diploid one-locus two-allele viability selection |
| `formal/ANSE/AutoresearchPopGenWrightFisher.lean` | H11: neutral Wright–Fisher harmonic characterisation of fixation |
| `formal/ANSE/AutoresearchPopGenControls.lean` | planted `sorry` / axiom controls (kept OUT of the claimed modules) |
| `data/autoresearch_popgen_2026-10-25/` | preregistration, program, results log, notes, independent verification reports, compile logs, sha256 MANIFEST |

## Build and audit
- Both files import only Mathlib.
- They were compiled with Lean 4.34.1 and the Mathlib of the pinned openai/math build (rc 0).
  - The compile logs are `data/autoresearch_popgen_2026-10-25/compile_*.json`, taken from these exact copies.
  - **Audit result:** `#print axioms` gives only [propext, Classical.choice, Quot.sound] for all 13 + 13 declarations.
  - **Controls:** the planted `sorryAx` and the planted axiom are both detected.
  - **Not tested:** this repo's own `lake build` was not run, because `formal/` has no lakefile in this checkout.

## H12: diploid mean-fitness monotonicity (namespace `AutoresearchPopGen`)
- **Model:** genotype fitnesses w11, w12, w22 > 0, allele frequency p ∈ [0,1], q = 1 − p, random mating.
  - Mean fitness: wbar = p²w11 + 2pq·w12 + q²w22.
  - Next frequency: p' = (p²w11 + pq·w12)/wbar.
- **Theorems:**
  - `pnext_mem_Icc`
  - `delta_p_eq`
  - `delta_p_eq_deriv`: Δp = pq·wbar'(p)/(2 wbar), with the real derivative
  - `wbar_increment`: the exact identity wbar(p') − wbar(p) = pq·s²·(wbar + p·w11 + q·w22)/wbar²
  - `fisher_mean_fitness_nondecreasing`
  - `hwStep_sum`, `hwStep_allele`, `hwStep_stationary` (Hardy–Weinberg after one round of random mating, then
    stationary)
- **Controls:**
  - `underdominance_moves_away`: Δp < 0 at p = 2/5 although mean fitness rises.
  - `fisher_fails_without_positivity`.
- **Verifier corrections, applied:**
  - This is one-locus mean-fitness monotonicity, not Fisher's general fundamental theorem.
  - Strict positivity is sufficient, not necessary.

## H11: neutral Wright–Fisher (namespace `AutoresearchPopGen`)
- **Kernel:** wfP N i j = C(N,j)(i/N)^j(1−i/N)^(N−j) on Fin (N+1), over ℝ.
- **Theorems:**
  - `wfP_sum_one`
  - `wfP_zero`, `wfP_last` (absorbing states)
  - `wfP_mean`: Σ_j P(i,j)·j = i
  - `frac_harmonic`, `wfP_interior_zero_pos`
  - `max_principle`, `harmonic_unique`
  - `wright_fisher_fixation_harmonic`
- **Controls:** `ctl_positive_wfP_value` and `ctl_negative_selection_not_harmonic`; the result depends on neutrality.
- **Verifier correction:** what is kernel-checked is the harmonic characterisation. i/N is the unique function harmonic
  at interior states with boundary values 0 and 1. No probability or path measure is built, so "absorption probability
  = i/N" needs the first-step-analysis step, which is NOT formalized.

## Overlap with this repository (read before merging)
`formal/ANSE/AlixPopulationGenetics.lean` already covers closely related ground:
- `aurelia_mean_fitness_nondecreasing` and `kal_increment_is_variance`: the **haploid** two-allele model over ℚ;
- `wfP`, `wfP_row`, `ilene_wf_harmonic` and `xav_wf_harmonic_solution_unique`: the Wright–Fisher harmonic
  uniqueness over ℚ.

So:
- **H11 is largely a duplicate** of the existing Wright–Fisher material (here over ℝ, with a separate maximum-principle
  lemma and a selection negative control). Merge it only if the ℝ version or the controls are wanted.
- **H12 is complementary:**
  - it covers the **diploid** viability model with genotype fitnesses and Hardy–Weinberg, where the existing file is
    haploid;
  - it has the exact increment identity and the derivative form of Δp;
  - it has the underdominance and positivity controls.

## Suggested next steps for the genetics project
1. Close the H11 gap: prove Pⁿ(i,N) → i/N, or the first-step lemma with Mathlib `ProbabilityTheory.Kernel`.
2. Moran model with selection: fixation probability (1 − r^{−i})/(1 − r^{−N}).
3. Multi-allele haploid Fisher identity: Δwbar = Var(w)/wbar.
4. Linkage disequilibrium decay (1 − r)^t under random mating.
5. A literature and Mathlib check before any novelty claim.

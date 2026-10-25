# Handoff: two kernel-checked population-genetics results (autoresearch H11, H12; 2026-10-25)

This is for the Claude session working on genetics. Both results are **formalizations of known textbook results**,
with novelty unchecked; they are not discoveries. Both were independently re-verified (producer ≠ verifier) and
accepted with corrections. Repo: github.com/xaviercallens/AutoevolveAI, branch `worktree-openai-math-discovery`
(commit `77dcce3` has the verified versions).

## Toolchain
- Pinned openai/math build: Lean 4.34.1 with full Mathlib. Compile from the AutoevolveAI worktree with
  `python3 scripts/certified_numerics/lean_pinned.py compile <file.lean> --json-out <out.json>`, which takes about
  15 s warm. Add `--olean Autoresearch.<Name>` to build an importable module.
- Only Mathlib is needed; nothing from openai/math is imported.
- After a reboot, warm the page cache first. Disk 2 is slow: see the memory note `disk2-cold-cache-and-detached-runs`.
- Axiom audit: every claimed theorem must depend only on [propext, Classical.choice, Quot.sound]. The planted controls
  live in `formal_cert/Autoresearch/PopGenControls.lean`, never in the claimed modules.

## H12: one-locus, two-allele mean-fitness monotonicity (+ Hardy–Weinberg)
File: `formal_cert/Autoresearch/PopGenFisher.lean`, namespace `AutoresearchPopGen`.
- **Model:** genotype fitnesses w11, w12, w22 > 0 and allele frequency p ∈ [0,1], with q = 1 − p.
  Mean fitness is wbar = p²w11 + 2pq·w12 + q²w22, and the next frequency is p' = (p²w11 + pq·w12)/wbar.
- **Theorems:**
  - `pnext_mem_Icc`
  - `delta_p_eq`: Δp = pq·(p(w11−w12) + q(w12−w22))/wbar
  - `delta_p_eq_deriv`: Δp = pq·wbar'(p)/(2 wbar)
  - `wbar_increment`: the exact identity wbar(p') − wbar(p) = pq·s²·(wbar + p·w11 + q·w22)/wbar²
  - `fisher_mean_fitness_nondecreasing`: wbar(p') ≥ wbar(p)
  - `hwStep_sum`, `hwStep_allele`, `hwStep_stationary` (Hardy–Weinberg)
- **Controls:** `underdominance_moves_away` (Δp < 0 at p = 2/5 while wbar rises) and
  `fisher_fails_without_positivity`.
- **Corrections from the verifier:**
  - This is not Fisher's general fundamental theorem (the additive-variance / partial-change reading).
  - Strict positivity is sufficient, not necessary: w ≥ 0 with wbar > 0 is enough.

## H11: neutral Wright–Fisher, fixation characterised by i/N
File: `formal_cert/Autoresearch/PopGenWrightFisher.lean`, namespace `AutoresearchPopGen`.
- **Kernel:** wfP N i j = C(N,j)(i/N)^j(1−i/N)^(N−j) on Fin (N+1).
- **Theorems:**
  - `wfP_sum_one`, `wfP_zero`, `wfP_last` (absorbing states)
  - `binom_mean` / `wfP_mean`: the martingale identity Σ_j P(i,j)·j = i
  - `frac_harmonic`, `wfP_interior_zero_pos`
  - `max_principle`, `harmonic_unique`, bundled in `wright_fisher_fixation_harmonic`
- **Controls:** `ctl_positive_wfP_value` and `ctl_negative_selection_not_harmonic`, which shows the result depends on
  neutrality.
- **Correction from the verifier:** what is kernel-checked is formulation (a), the harmonic characterisation. There is
  no probability or path measure, so "absorption probability = i/N" is NOT kernel-checked. The missing step is the
  first-step analysis that the absorption probability is harmonic.

## Suggested next steps
1. **Close the H11 gap.** Prove Pⁿ(i,N) → i/N (the chain is absorbed with probability 1 because each interior state
   reaches 0 in one step with positive probability, and the martingale is bounded). Alternatively, build the path
   measure with Mathlib's `ProbabilityTheory.Kernel` and prove the first-step lemma. Either one makes
   "fixation probability = i/N" Tier A.
2. **Selection.** Fixation probability with selection (Kimura's diffusion limit is analysis-heavy; the discrete Moran
   model with selection has the exact formula (1 − r^{−i})/(1 − r^{−N})). Moran is the realistic next target.
3. **Fisher's general theorem.** Δwbar = Var_A(w)/wbar for multi-allele haploid selection is a clean finite-sum
   identity in Lean. It is a good first step towards the additive-genetic-variance reading.
4. **Multi-locus / linkage.** Hardy–Weinberg with linkage disequilibrium decaying as (1 − r)^t.
5. **Literature check** before any novelty claim. Mathlib or Lean community projects may already contain some of this.

Full records: `results/autoresearch/2026-10-25/H11/` and `H12/` (NOTE.md, result.json, verification.md,
compile_*.json), plus `results/autoresearch/2026-10-25/SUMMARY.md`.

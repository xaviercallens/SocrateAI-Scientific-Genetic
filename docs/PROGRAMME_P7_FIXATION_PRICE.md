# Programme P7: kernel-checked fixation probabilities and the Price equation (2026-10-26)

Paper: doi:10.5281/zenodo.23653850. It is an Elenchus-tiered preprint, independently audited by a separate AI agent
(ACCEPT_WITH_CORRECTIONS, corrections applied). No human has reviewed it, and the 14 Tier A ledger rows stay flagged
LEDGER_UNAUDITED_TIER_A.

## What is new relative to PR #1 (H11/H12)
- **Fixation probability = absorption limit, kernel-checked.** For a general absorbing chain with positive one-step
  variance, Pⁿ(i,N) → h(i). The cases covered:
  - neutral Wright–Fisher: **i/N** (`wright_fisher_fixation_limit`, `fixation_limit`, `absorption_tendsto`,
    `interiorMass_tendsto_zero`). This closes the H11 gap.
  - **Moran with selection**: fixation probability (1 − r^{−i})/(1 − r^{−N}), written S_i/S_N in
    `MoranFixation.lean`. This is the standard model for the chance that a mutation with relative fitness r (for
    example drug resistance) takes over a population of size N.
- **The Price equation** for a finite population, as an exact identity (`Price.lean`). Its haploid corollary is
  Δwbar = Var(w)/wbar with no transmission bias. A negative control shows that transmission bias breaks the variance
  formula. This is not Fisher's general "partial change" theorem.
- `FixationDataCert.lean` holds data certificates for the companion rusty-SUNDIALS crate `qf-popgen`
  (rusty-SUNDIALS PR #80).

## Layout and build
The Lean module root is `formal/`, matching the source layout, so the imports
`Autoresearch.PopGenWrightFisher` / `Programmes.PopGen.*` resolve.
- **Shared files with PR #1.** `formal/Autoresearch/PopGen{Fisher,WrightFisher,Controls}.lean` are the same files
  PR #1 adds under `formal/ANSE/` with different names. If both PRs are merged, keep one copy.
- **How it was compiled.** All files compile with Lean 4.34.1 + Mathlib d13f23b7 (the pinned openai/math build), and
  every claimed theorem is within [propext, Classical.choice, Quot.sound]. The planted controls live only in the
  `*Controls.lean` files.
- **Not tested here.** This repository's own `lake build` was not run; `formal/` has no lakefile.

## Overlap
`formal/ANSE/AlixPopulationGenetics.lean` already proves the harmonic uniqueness for Wright–Fisher over ℚ. The
absorption-limit statements and the Moran-with-selection and Price results here are not in it.

Records: `data/programme_p7_2026-10-26/` (FIXATION.md, PRICE.md, compile logs, MANIFEST with sha256).
Paper: `papers/programme_p7_fixation/`.

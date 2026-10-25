# H12 independent verification (producer ≠ verifier), 2026-10-25

**VERDICT: ACCEPT_WITH_CORRECTIONS.** The theorems are correct and match the prereg. The corrections concern
wording and hygiene.

## 1. Recompile + axioms
- I made a verbatim copy of `PopGenFisher.lean`, appended a verifier block, and compiled it with pinned Lean 4.34.1
  (`lean_pinned.py`, cores 4,5). Result: **rc 0, 12.0 s**.
- My `#print axioms`, verbatim:
  - `fisher_mean_fitness_nondecreasing`, `wbar_increment`, `delta_p_eq_deriv`, `hwStep_stationary` (and every other
    non-control theorem in the file): `[propext, Classical.choice, Quot.sound]`.
  - My planted `VerifierH12.v_planted : wbar 1 1 1 p ≤ 0 := sorry`: `[propext, sorryAx, Classical.choice, Quot.sound]`.
    Caught.
  - The producer's in-file controls: `ctl_sorry` → `[sorryAx]`, `ctl_axiom` → `[AutoresearchPopGen.ctl_cheat]`.

## 2. Statement adequacy
- **The model is the standard one-locus, two-allele viability-selection model with HW zygotes.** Pinned by `rfl`:
  - `wbar = p²w11 + 2p(1−p)w12 + (1−p)²w22`;
  - `pnext = (p²w11 + p(1−p)w12)/wbar`.
  Textbook check, which I proved: genic selection w = ((1+s)², 1+s, 1) gives p' = p(1+s)/(1+ps).
- **The headline is exactly `wbar(p) ≤ wbar(p')`**, assuming `0 < w11, w12, w22` and `0 ≤ p ≤ 1`. This is
  confirmed by `#check`.
- **Non-vacuity.** I proved strict increase at w = (2,1,1), p = 1/2.
- **Division junk.** None: `wbar_pos` gives wbar > 0 on the whole domain, including p = 0 and p = 1, so `pnext`
  never divides by 0.
- **Wright's form.** `delta_p_eq_deriv` uses the real `deriv` of `fun x => wbar … x`, obtained from an explicit
  `HasDerivAt`. It is not a hand-written polynomial, so the prereg's "Δp = pq·wbar'/(2 wbar)" is met.
- **`wbar_increment` is an exact identity.** `wbar(p') − wbar(p) = pq s² (wbar + p w11 + q w22)/wbar²`. Every factor
  is ≥ 0, so the inequality follows honestly.
- **`hwStep_stationary`** holds for every triple and does not need x11 + x12 + x22 = 1, because
  p² + p(1−p) = p identically. The statement is correct. `hwStep` ignores x22, so it presumes normalised input.

## Corrections
1. **Positivity is overstated.** NOTE.md says "the positivity hypotheses are therefore needed". The control
   `fisher_fails_without_positivity` uses w11 = −1, so it only shows that a *negative* fitness breaks the result.
   Strict positivity is not necessary: by `wbar_increment`, w ≥ 0 together with wbar(p) > 0 already gives
   wbar + p w11 + q w22 ≥ wbar > 0. The theorem is correct but narrower than it could be. Reword to "nonnegativity
   is needed; the strict-positivity hypothesis is sufficient but not necessary".
2. **Naming.** The Lean docstring of `fisher_mean_fitness_nondecreasing` says "(iii) Fisher's fundamental theorem",
   and so do the module header and the prereg title. Rename to "one-locus two-allele mean-fitness monotonicity". The
   NOTE scope section already says that this is not Fisher's FTNS (the additive-variance partial change); the Lean
   text should say so too.
3. **Hygiene.** `PopGenFisher.lean` declares `axiom ctl_cheat : False` and a `sorry` theorem inside the claimed module
   and namespace, so any importer inherits a `False` axiom. Its names `AutoresearchPopGen.ctl_sorry`, `ctl_cheat` and
   `ctl_axiom` clash with `PopGenWrightFisher.lean`, so the two modules cannot be imported together. Move the planted
   controls to `PopGenFisherControls.lean`.
4. **"First formalization".** This is unverified: no literature or Mathlib/AFP search was done by producer or
   verifier. Reword to "a formalization (novelty unchecked)".

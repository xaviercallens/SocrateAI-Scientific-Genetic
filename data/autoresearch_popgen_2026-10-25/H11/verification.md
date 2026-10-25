# H11 independent verification (producer ≠ verifier), 2026-10-25

**VERDICT: ACCEPT_WITH_CORRECTIONS.** The Lean theorems are correct and correctly stated. The preregistered claim
("absorption probability at N from i equals i/N") is only **partially** met.

## 1. Recompile + axioms
- I made a verbatim copy of `PopGenWrightFisher.lean`, appended a verifier block, and compiled it with pinned Lean
  4.34.1 (`lean_pinned.py`, cores 4,5). Result: **rc 0, 10.8 s**.
- My `#print axioms`, verbatim:
  - `wright_fisher_fixation_harmonic`, `harmonic_unique`, `frac_harmonic` (and every other non-control theorem in
    the file): `[propext, Classical.choice, Quot.sound]`.
  - My planted `VerifierH11.v_planted : wfP 2 1 1 = 1 := sorry`: `[propext, sorryAx, Classical.choice, Quot.sound]`.
    Caught.
  - The producer's in-file controls: `ctl_sorry` → `[sorryAx]`, `ctl_axiom` → `[AutoresearchPopGen.ctl_cheat]`.

## 2. Statement adequacy
- **The kernel is the neutral Wright–Fisher kernel.** `wfP N i j = C(N,j)(i/N)^j(1−i/N)^(N−j)`, pinned by `rfl`.
  This is binomial resampling of N genes at frequency i/N, with real division. My own values `wfP 3 1 0 = 8/27`,
  `wfP 3 1 2 = 2/9` and `wfP 2 1 2 = 1/4` are proved by `norm_num`. At the boundaries, 0^0 = 1 gives the absorbing
  rows (`wfP_zero`, `wfP_last`).
- **Quantifiers are correct.** The statement requires `0 < N`. Harmonicity `g i = Σ_j P(i,j) g j` is assumed only at
  interior states 0 < i < N, with boundary values `g 0 = 0` and `g (last N) = 1`, and concludes ∀ i, g i = i/N. That
  is the right uniqueness statement. With N = 1 there are no interior states, and the conclusion follows from the
  boundary values alone, which is correct.
- **The hypothesis set is satisfiable.** I instantiated the theorem at N = 3 with g = i/3, using `frac_harmonic`.
- **The maximum principle is genuine.** It uses `wfP_interior_zero_pos` (P(i,0) = (1−i/N)^N > 0) and stochasticity.
  No junk division occurs: N > 0 everywhere it matters.

## 3. Gap: the prereg statement is not fully proved
- No path measure or absorption probability appears anywhere in the file. What is kernel-checked is:
  (i) i/N is harmonic, and (ii) it is the unique interior-harmonic function with boundary values 0 and 1.
- "Fixation probability = i/N" additionally needs first-step analysis: the absorption probability is harmonic at
  interior states. That step is **not formalised**.
- NOTE.md and result.json `formulation_achieved` disclose this correctly.

## Corrections (required)
1. **Mark the result partial.** Label it everywhere as "formulation (a): harmonic characterisation". This applies to
   the result.json `verdict`, any summary, and commit/paper text. Do not write "fixation probability i/N
   kernel-checked". Commit ab72f4d's wording ("fixation i/N as the unique harmonic function (first-step
   identification not formalized)") is acceptable.
2. **Route to close the gap.** Either define a_n(i) = Pⁿ(i, N) and prove `a_n(i) → i/N` (no measure theory needed:
   the martingale identity gives Σ_j Pⁿ(i,j) j = i, and the mass on interior states → 0 because P(i,0) > 0), or
   define the absorption probability through Mathlib's Markov-kernel API and prove first-step harmonicity.
3. **Hygiene.** `PopGenWrightFisher.lean` itself declares `axiom ctl_cheat : False` and `theorem ctl_sorry := sorry`
   inside the claimed module and namespace. Two consequences:
   - any importer inherits a `False` axiom;
   - the file clashes with `PopGenFisher.lean`, which declares the same `AutoresearchPopGen.ctl_*` names, so the two
     modules cannot be imported together.
   Move the planted controls to a separate `*Controls.lean`, as H2 and H6 do.
4. **"First formalization".** The claim is unverified: no literature or Mathlib search was done by producer or
   verifier. Reword to "a formalization (novelty unchecked)".

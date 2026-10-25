# H12: Fisher's theorem, one locus, two alleles (first formalization of a known result)

Lean: `formal_cert/Autoresearch/PopGenFisher.lean`, namespace `AutoresearchPopGen`. Toolchain: pinned openai/math
Lean 4.34.1 + Mathlib. Attempt 1, about 4 min, compile rc 0. Verdict: KEEP, pending producer ≠ verifier (not Tier A).

This is a classical result (Fisher 1930; textbook treatments in Ewens, Hartl & Clark). It is not a discovery.

## Model
Fitnesses `w11, w12, w22 > 0`, `p ∈ [0,1]`, `q = 1 - p`, `wbar = p²w11 + 2pq w12 + q²w22`,
`p' = (p²w11 + pq w12)/wbar`.

## Human proof
1. `wbar > 0`: if `p > 0` the term `p²w11` is positive and the others are ≥ 0. If `p = 0`, `wbar = w22 > 0`.
2. `0 ≤ p' ≤ 1`: the numerator is ≥ 0, and `wbar − numerator = pq w12 + q²w22 ≥ 0`.
3. `Δp = p' − p = pq·s/wbar`, with `s = p(w11 − w12) + q(w12 − w22)`. This follows by clearing denominators.
   Since `wbar(x) = w22 + 2(w12 − w22)x + (w11 − 2w12 + w22)x²`, we have `wbar'(p) = 2s`, which gives
   Wright's form `Δp = pq·wbar'(p)/(2·wbar)`. In Lean, `wbar'` is the real derivative `deriv`, obtained from an
   explicit `HasDerivAt`.
4. Exact increment identity, found with sympy and then checked by the kernel:
   `wbar(p') − wbar(p) = pq·s²·(wbar + p·w11 + q·w22)/wbar²`.
   Every factor is ≥ 0 under the hypotheses, so `wbar(p') ≥ wbar(p)`.
5. Hardy–Weinberg: from `(x11, x12, x22)`, set `p = x11 + x12/2`. The offspring are `(p², 2pq, q²)` and sum to 1.
   The allele frequency is preserved (`p² + pq = p`). Applying the step again gives the same proportions, so HW
   proportions are stationary after one round.

## Controls
- Planted: `ctl_sorry` shows `[sorryAx]` and `ctl_axiom` shows `[AutoresearchPopGen.ctl_cheat]`. Both were caught.
  The grep for sorry/admit/native_decide/axiom matches only these planted lines.
- `underdominance_moves_away`: `w11 = w22 = 1`, `w12 = 1/2` has a fixed interior point at `p = 1/2`. From
  `p = 2/5` the frequency moves away from it (`Δp < 0`), yet `wbar` still rises. So rising mean fitness does not
  mean convergence to a polymorphism.
- `fisher_fails_without_positivity`: `w11 = −1`, `w12 = 1`, `w22 = 0`, `p = 1/2` gives `wbar = 1/4 > 0` and `p' = 0`,
  and `wbar` drops to 0. The positivity hypotheses are therefore needed. By the identity, the sign of the change
  is the sign of `wbar + p·w11 + q·w22`.

## Scope
- This covers one locus, two alleles, discrete non-overlapping generations, viability selection only, and random
  mating with HW zygote proportions. Fitnesses are constant (no frequency dependence).
- This is the mean-fitness monotonicity statement. It is **not** Fisher's general fundamental theorem of natural
  selection (the partial-change / additive-genetic-variance interpretation), and it does not cover multiple loci
  or multiple alleles.
- Producer ≠ verifier is still owed: an independent recompile and an audit of the statement's adequacy.

## Independent verification: ACCEPT_WITH_CORRECTIONS (verification.md)
- Name: one-locus, two-allele **mean-fitness monotonicity**. It is not Fisher's general fundamental theorem (the
  additive-variance / partial-change reading).
- Strict positivity of the fitnesses is sufficient, not necessary (w ≥ 0 with wbar > 0 is enough). The control only
  shows that a negative fitness can break monotonicity.
- Hygiene fixed: the planted controls moved to `PopGenControls.lean` (compile_fixed.json, rc 0, whitelist only).
- Novelty unchecked; no literature search was done.

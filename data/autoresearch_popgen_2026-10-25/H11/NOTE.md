# H11: Neutral Wright–Fisher fixation probability i/N (first formalization of a known result)

Lean: `formal_cert/Autoresearch/PopGenWrightFisher.lean`, namespace `AutoresearchPopGen`. Toolchain: pinned
openai/math Lean 4.34.1 + Mathlib. Attempt 1, about 4 min, compile rc 0. Verdict: KEEP for formulation (a),
pending producer ≠ verifier (not Tier A).

This is a classical result (Fisher, Wright; Ewens, *Mathematical Population Genetics*). It is not a discovery.

## Formulation achieved: (a), the harmonic-function characterisation
The kernel is `wfP N i j = C(N,j)·(i/N)^j·(1 − i/N)^(N−j)` on `Fin (N+1)`, with real division. Proved, for `0 < N`
where needed:
- `wfP_nonneg`, `wfP_sum_one`: the kernel is stochastic (binomial theorem).
- `wfP_zero`, `wfP_last`: states 0 and N are absorbing.
- `binom_mean`, `wfP_mean`: the martingale identity `Σ_j P(i,j)·j = i`.
- `frac_harmonic`: `h(i) = i/N` is harmonic at every state.
- `wfP_interior_zero_pos`: from every interior state, `P(i,0) = (1 − i/N)^N > 0`.
- `max_principle`, `harmonic_unique`, `wright_fisher_fixation_harmonic`: every `g` that is harmonic at interior
  states with `g 0 = 0` and `g N = 1` equals `i/N`.

## Human proof
- Mean: use `k·C(n+1,k) = (n+1)·C(n,k−1)`. Shift the index, then apply the binomial theorem to get
  `Σ_k C(N,k)x^k y^(N−k) k = N·x·(x+y)^(N−1)`. With `x = i/N` and `y = 1 − x` this equals `i`.
- Uniqueness: let `d = g − i/N`. Then `d` is harmonic inside and zero at both ends. Take a maximiser `m` of `d`.
  If `m` is a boundary state, `max d = 0`. If `m` is interior, `Σ_j P(m,j)(d(m) − d(j)) = 0` and every term is ≥ 0,
  so every term vanishes. Because `P(m,0) > 0`, this gives `d(m) = d(0) = 0`. So `d ≤ 0`, and applying the same
  argument to `−d` gives `d = 0`.

## Not formalized (honest gap)
No Markov-chain path measure is constructed. The standard step "the absorption probability at N is harmonic at
interior states, with boundary values 0 and 1" (first-step analysis) is **not** kernel-checked. The statement
"fixation probability = i/N" follows from the kernel-checked uniqueness theorem only together with that
unformalized step. A stronger variant, `lim_n P^n(i,N) = i/N`, would close the gap without measure theory. It was
not attempted.

## Controls
- Planted: `ctl_sorry` shows `[sorryAx]` and `ctl_axiom` shows `[AutoresearchPopGen.ctl_cheat]`. Both were caught.
  The grep matches only these planted lines.
- Positive: `ctl_positive_wfP_value`: `wfP 2 1 1 = 1/2`. This pins definition fidelity: real division, not ℕ
  division.
- Negative: `ctl_negative_selection_not_harmonic`: with selection the kernel samples with probability `ψ ≠ i/N`.
  For `N = 2`, `i = 1`, `ψ = 2/3`, the expected next frequency is `2/3 ≠ 1/2`, so `i/N` is not harmonic and the
  result really depends on neutrality.

## Independent verification: ACCEPT_WITH_CORRECTIONS (verification.md)
- **What is kernel-checked is formulation (a), the harmonic characterisation.** i/N is the unique function harmonic
  at interior states with values 0 and 1 at the boundary. The file has no probability object: the first-step link
  "absorption probability is harmonic" is not formalised. So this file does NOT kernel-check "fixation probability = i/N".
  The gap can be closed by proving P^n(i,N) → i/N, or by a Markov-kernel first-step lemma.
- Hygiene fixed: the planted controls moved to `PopGenControls.lean`. The claimed module no longer declares an axiom,
  and it imports together with PopGenFisher (`PopGenImportBoth.lean`, rc 0; compile_fixed.json, compile_import_both.json).
- Novelty unchecked; no literature search was done ("first formalization" is not claimed).

# P7 WORK1: fixation probabilities (neutral Wright–Fisher, Moran with selection)

Lane P7, stage WORK1, 2026-10-25. Producer: Claude Opus 5.5 subagent. **Producer ≠ verifier is still pending**, so
nothing below is final Tier A until an independent recompile and statement audit is done.

Preregistration: `results/programmes/2026-10-25/preregistration.json`, lane P7. Priors were `H11_gap` 0.5 and
`moran_selection` 0.6. Both statements were attempted, and both compiled on the first budget.

All results are formalizations of **known classical results**. Fixation `i/N` is due to Fisher and Wright; the
Moran formula is due to Moran (1958), as cited in Keller–Uğurlu (arXiv:2403.12598, eq. (4)); see §6. **Novelty is
not claimed.** The search for prior formalizations was limited to two alphaXiv searches (§6). Nothing here depends on
openai/math results; only Mathlib is used.

## 1. What is kernel-checked

**Toolchain.**
- Pinned openai/math Lean 4.34.1 + Mathlib, run through `scripts/certified_numerics/lean_pinned.py` on cores 4-7.
- Every claimed theorem below prints `[propext, Classical.choice, Quot.sound]`, verbatim from the compile JSONs in
  this directory.
- There are no `sorry`, `admit`, `native_decide` or `axiom` tokens in the claimed files. A grep matches only a
  doc-comment line in `MoranModel.lean:237`, plus the planted-controls file.

| File (`formal_cert/Programmes/PopGen/`) | Module | sha256 | Compile record |
|---|---|---|---|
| `Absorption.lean` | `Programmes.PopGen.Absorption` | `86fb5168a6b1efe07b8d050689da88facba0350a798b49df72315d7c283c5f95` | `compile_Absorption.json` (rc 0, 27.5 s) |
| `WrightFisherFixation.lean` | `Programmes.PopGen.WrightFisherFixation` | `149a9328171a6347460e8b2b1294d8073df57e2d6a602cdacf5955a3203262b3` | `compile_WrightFisherFixation.json` (rc 0) |
| `MoranModel.lean` | `Programmes.PopGen.MoranModel` | `b0822dc0ee83505fdffbff00ae2140cfa190dee7c199ecefd156eea251da590f` | `compile_MoranModel.json` (rc 0) |
| `MoranFixation.lean` | `Programmes.PopGen.MoranFixation` | `c49806ac8b14295ddebe032c69b97503e3a033d1cac503530948d72e8bdad8e5` | `compile_MoranFixation.json` (rc 0) |
| `PopGenProgImportAll.lean` | (check only) | `e26cdffa01ab893b92a542540e74f2be07525dd4964f82c3db468a2ec84f23ec` | `compile_import_all.json` (rc 0) |
| `PopGenProgControls.lean` | (planted, never imported) | `685cdd545b9cff73387d958242af254a6d3ffc987868bfa61d5b215b0988c9de` | `compile_controls.json` (rc 0) |

The H11 module it builds on is `formal_cert/Autoresearch/PopGenWrightFisher.lean`, sha256
`9cac6b9d5d0c76be771d8b004777699c539ea384741f7483f362de5884a4e530`. It is unchanged and is imported, not copied.

### 1a. General absorbing-chain theorem (`Absorption.lean`, namespace `ProgPopGen`)

**Objects.**
- `P : Fin (N+1) → Fin (N+1) → ℝ` is the transition kernel.
- `kPow P n` is its `n`-step kernel, defined by the recursion `kPow P 0 i j = [i = j]` and
  `kPow P (n+1) i j = ∑ₖ kPow P n i k · P k j` (Chapman–Kolmogorov).
- `AbsorbingChain P` means: the entries are ≥ 0, every row sums to 1, row 0 is `δ₀`, and row N is `δ_N`.
- `InteriorHarmonic P g` means `∑ⱼ P i j · g j = g i` for every `0 < i < N`.
- `interiorMass P n i = ∑ₖ Pⁿ(i,k)·[0<k<N]`.

**Theorems.**
- `fixation_limit`. Assume `AbsorbingChain P` and `0 < N`, and let `h` satisfy:
  - `h` is interior-harmonic, with `h 0 = 0` and `h N = 1`;
  - the one-step variance `∑ⱼ P(k,j)(h j − h k)²` is `> 0` at every interior `k`.

  Then **`Pⁿ(i,N) → h i`** and **`Pⁿ(i,0) → 1 − h i`**.
- `harmonic_unique_of_var`. Under the same hypotheses, every interior-harmonic `g` with `g 0 = 0` and `g N = 1`
  equals `h`.
- Supporting lemmas:
  - `kPow_nonneg`, `kPow_rowSum`;
  - `kPow_harmonic` (the martingale identity `∑ⱼ Pⁿ(i,j) g j = g i`);
  - `second_moment`;
  - `interiorMass_tendsto_zero`;
  - `absorption_tendsto`;
  - `mass_split` (`Pⁿ(i,0) + Pⁿ(i,N) + Iₙ(i) = 1`);
  - `kPow_last_mono` (`Pⁿ(i,N)` is non-decreasing in n).

### 1b. Neutral Wright–Fisher: the H11 gap is closed (`WrightFisherFixation.lean`)

- `wright_fisher_fixation_limit (hN : 0 < N) (i)`:
  `Tendsto (fun n => kPow (wfP N) n i (Fin.last N)) atTop (𝓝 (i/N))`, and the companion
  `Tendsto (fun n => kPow (wfP N) n i 0) atTop (𝓝 (1 − i/N))`.
- `wfP` is H11's kernel `C(N,j)(i/N)^j(1−i/N)^(N−j)`, unchanged.
- `wf_var_pos` supplies the variance hypothesis. The jump to 0 alone contributes `P(k,0)(k/N)² > 0`, using H11's
  `wfP_interior_zero_pos`.
- `wright_fisher_fixation_mono`: `Pⁿ(i,N) ≤ Pⁿ⁺¹(i,N)`.

### 1c. Moran model with selection (`MoranModel.lean`, `MoranFixation.lean`)

**The model.**
- `moranUp N r i = (r i/(r i + N − i))·((N − i)/N)` and `moranDown N r i = ((N − i)/(r i + N − i))·(i/N)`. These are
  exactly the task's `p_i` and `q_i`, in real arithmetic.
- `moranP N r` is the birth–death kernel built from them (`bdP`): up `p_i`, down `q_i`, stay `1 − p_i − q_i`.
- `moranFix N r i = S_i/S_N`, where `S_n = ∑_{k<n} r^{−k}`.

**ℕ-level lemmas (`MoranModel.lean`).**
- `bdP_sum` gives the first-step expectation of a birth–death kernel.
- `moranUp_eq`: `p_i = r·q_i`.
- `moran_stay_eq`: `1 − p − q = (r i² + (N−i)²)/((r i + N − i)N) ≥ 0`.
- `moranFix_harmonic_nat`: the first-step equation `p_i h(i+1) + q_i h(i−1) + (1−p_i−q_i) h(i) = h(i)` at every
  `0 < i`.
- `moran_var_pos_nat`: the variance is positive at interior i.
- `moranFix_closed` (`r ≠ 1`): `moranFix N r i = (1 − r⁻¹^i)/(1 − r⁻¹^N)`.
- `moranFix_neutral`: `moranFix N 1 i = i/N` exactly.
- `moranFix_tendsto_one`: `moranFix N r i → i/N` as `r → 1`.
- `moran_closed_tendsto_neutral`: `(1 − r⁻¹^i)/(1 − r⁻¹^N) → i/N` as `r → 1` with `r ≠ 1` (`𝓝[≠] 1`). This settles
  the task's "neutral limit r → 1 gives i/N" as a proof, not just a check.

**Chain-level theorems (`MoranFixation.lean`).**
- `moran_absorbing` holds for `0 < N` and `0 < r`.
- `moran_interiorHarmonic`.
- `moran_fixation_unique`: uniqueness among interior-harmonic functions with boundary values 0 and 1.
- `moran_fixation_limit`: `Pⁿ(i,N) → moranFix N r i` and `Pⁿ(i,0) → 1 − moranFix N r i`, for every `r > 0`.
- `moran_fixation_limit_closed` (`r ≠ 1`): `Pⁿ(i,N) → (1 − r^{−i})/(1 − r^{−N})`.
- `moran_fixation_limit_neutral` (`r = 1`): `Pⁿ(i,N) → i/N`.

### 1d. Import-together and non-vacuity (`PopGenProgImportAll.lean`)

- H11 (`Autoresearch.PopGenWrightFisher`), H12 (`Autoresearch.PopGenFisher`) and the P7 modules import together
  without name clashes.
- The hypotheses are satisfiable:
  - `nonvacuous_moran_chain : AbsorbingChain (moranP 3 2)`;
  - `nonvacuous_wf_chain : AbsorbingChain (wfP 3)`.
- Concrete instances:
  - `instance_wf_N3`: `Pⁿ(1,3) → 1/3`;
  - `instance_moran_N10`: for `r = 3/2`, `Pⁿ(1,10) → (1 − 2/3)/(1 − (2/3)^10)`.

### 1e. Attempt log (time-boxed, keep/discard, the kernel as judge)

| File | Attempts | Outcome |
|---|---|---|
| `Absorption.lean` | 2 | Attempt 1: rc 1. `dsimp` made no progress on an already beta-reduced hypothesis; the line was removed. Attempt 2: rc 0. **KEEP** |
| `MoranModel.lean` | 1 | rc 0. **KEEP** |
| `WrightFisherFixation.lean` | 1 | rc 0. **KEEP** |
| `MoranFixation.lean` | 2 | Attempt 1: rc 1. `rw` could not match the higher-order pattern `?f ↑j` in `moran_var_pos`. It was fixed by passing `f` explicitly and closing with `exact`. Attempt 2: rc 0. **KEEP** |
| `PopGenProgImportAll.lean`, `PopGenProgControls.lean` | 1 each | rc 0. Planted controls caught. |

**Axiom-line counts in the final compile records.** Every line is `[propext, Classical.choice, Quot.sound]`:

| Record | Lines |
|---|---|
| `compile_Absorption.json` | 11 |
| `compile_MoranModel.json` | 16 |
| `compile_WrightFisherFixation.json` | 5 |
| `compile_MoranFixation.json` | 8 |
| `compile_import_all.json` | 6 |

The only lines outside the whitelist are the two planted controls in `compile_controls.json`. Total wall time for
the stage was about 75 minutes.

## 2. Bridge to "fixation probability" (stated explicitly, not part of the Tier A claim)

**What is kernel-checked** is the convergence of the n-step kernel `Pⁿ(i,N)`, which the recursion `kPow` defines.

**Identification with probability.**
- `Pⁿ(i,N) = ℙᵢ(Xₙ = N)` is the definition of the n-step law of a time-homogeneous Markov chain with kernel P.
- Because N is absorbing, this equals `ℙᵢ(T_N ≤ n)`. That is consistent with the kernel-checked monotonicity
  `kPow_last_mono`.
- The fixation probability `ℙᵢ(T_N < ∞)` is therefore the limit.

**What is not formalized.** No path space or `ProbabilityTheory.Kernel` trajectory measure is built in Lean. The
identification above is the standard definition (Tier L/C glue). It is not a kernel-checked statement.

**Wording to use:** "the n-step absorption probability at N converges to i/N (resp. (1 − r^{−i})/(1 − r^{−N}));
kernel-checked". Do not write "fixation probability proved in Lean" without this bridge sentence.

## 3. Human proof (step by step, matching the Lean proof)

1. **Martingale.** If `Ph = h` at every state, then `Pⁿh = h`, by induction and one exchange of finite sums.
   Interior harmonicity plus absorbing rows gives `Ph = h` everywhere.
2. **Second moment.** `∑ⱼ P(k,j) h(j)² = h(k)² + ∑ⱼ P(k,j)(h(j) − h(k))²` (expand the square; use `∑P = 1` and
   `Ph = h`).
   - Let `ε > 0` be the minimum of the variance over the finitely many interior states.
   - Then `P(h²) ≥ h² + ε·1_interior` pointwise; at a boundary state the inequality is an equality with
     `1_interior = 0`.
3. **Accumulation.**
   - Applying `Pⁿ`, which has non-negative entries, gives `∑ⱼ Pⁿ(i,j) h(j)² ≥ h(i)² + ε ∑_{m<n} I_m(i)`, where
     `I_m(i)` is the interior mass at time m.
   - The left side is at most `G := ∑ⱼ h(j)²`, since the rows of `Pⁿ` sum to 1.
   - Hence `∑_{m<n} I_m(i) ≤ G/ε` for all n. The non-negative series `∑ I_m(i)` converges, so `I_n(i) → 0`.
4. **Squeeze.**
   - Let `H := ∑ⱼ |h(j)|`. Pointwise, `[j = N] − H·1_interior(j) ≤ h(j) ≤ [j = N] + H·1_interior(j)`.
   - Applying `Pⁿ(i,·)` and using `Pⁿh = h` gives `|Pⁿ(i,N) − h(i)| ≤ H·Iₙ(i) → 0`.
   - Then `Pⁿ(i,0) = 1 − Pⁿ(i,N) − Iₙ(i) → 1 − h(i)`.
5. **Uniqueness.** Any other interior-harmonic g with boundary values 0 and 1 also satisfies `Pⁿ(i,N) → g(i)`, by
   steps 1 and 4 with the same `Iₙ → 0`. Limits are unique, so `g = h`. This replaces the maximum principle, which
   H11 used for Wright–Fisher.
6. **Wright–Fisher.**
   - `h(j) = j/N` is harmonic (H11 `frac_harmonic`, the binomial mean).
   - At interior k, the variance is ≥ `P(k,0)(k/N)² = (1 − k/N)^N (k/N)² > 0`.
7. **Moran.**
   - The ratio is `p_i/q_i = r`.
   - With `S_n = ∑_{k<n} r^{−k}`, `h(i) = S_i/S_N` gives `h(i+1) − h(i) = r^{−i}/S_N`.
   - So `p_i (h(i+1) − h(i)) = r q_i r^{−i}/S_N = q_i r^{−(i−1)}/S_N = q_i (h(i) − h(i−1))`, which is the
     first-step equation.
   - The variance is ≥ `p_i (r^{−i}/S_N)² > 0` at interior i.
   - The geometric sum gives `S_i/S_N = (1 − r^{−i})/(1 − r^{−N})` for `r ≠ 1`. At `r = 1`, `S_i = i`, so the value
     is i/N. `S_i/S_N` is continuous in r at r = 1, which gives the neutral limit.

## 4. Controls

**Planted (in `PopGenProgControls.lean`, caught; verbatim from `compile_controls.json`).**
- `ctl_sorry`, which claims the wrong limit 1/7 with `sorry`, prints
  `[propext, sorryAx, Classical.choice, Quot.sound]`.
- `ctl_axiom`, which claims the wrong limit 1/3 from `axiom ctl_cheat : False`, prints
  `[propext, Classical.choice, ProgPopGenControls.ctl_cheat, Quot.sound]`.

**Positive controls (kernel-checked, axiom-clean).**
- `ctl_pos_wf_two_steps`: `P²(1,2) = 3/8` for Wright–Fisher with N = 2. This is a finite-time value, not yet the
  limit 1/2.
- `ctl_pos_moran_rates`: `p_1 = 1/3` and `q_1 = 1/6` (N = 2, r = 2).
- `ctl_pos_moranFix_N2`: `2/3`, with the closed form agreeing.
- `ctl_pos_moranFix_N3`: `4/7` and `6/7`, and the first-step equation holds with these numbers.

**Negative controls (kernel-checked, axiom-clean).**
- `ctl_neg_identity_no_absorption`:
  - The identity kernel is stochastic and absorbing, and `j/2` is harmonic for it.
  - Yet `Pⁿ(1,2) = 0` for all n, so `Pⁿ(1,2) ↛ 1/2`.
  - So the variance hypothesis of `fixation_limit` cannot be dropped.
- `ctl_neg_flipped_value`: the exponent-flipped formula `(1 − r^i)/(1 − r^N)` gives `1/7 ≠ 4/7` (N = 3, r = 2,
  i = 1).
- `ctl_neg_flipped_not_harmonic`: the flipped formula violates the first-step equation (3/14 ≠ 1/7).
- `ctl_neg_neutral_not_harmonic`: under selection, i/N violates the first-step equation (7/18 ≠ 1/3).
- `ctl_neg_wrong_limits`: at chain level, `Pⁿ(1,3)` for N = 3, r = 2 does **not** tend to 1/3, and does **not** tend
  to 1/7.

## 5. Tier X cross-check (`numerics_fixation.py` → `numerics_fixation.json`)

This check never supports a claim.

**Exact rational arithmetic.**
- `P²(1,2) = 3/8`.
- Solving the Moran first-step system for N = 3, r = 2 exactly gives `[0, 4/7, 6/7, 1]`.
- For N = 10, r = 3/2, the exact solve equals the closed form at every i (`true`) and does not equal the flipped
  formula (`false`).

**mpmath at 60 digits, `P^4096`, N = 10.**
- Wright–Fisher: `|Pⁿ(i,N) − i/N|` ≤ 2e-61, which is at the 60-digit precision floor. Interior mass is
  1.8e-188 to 4.3e-188.
- Moran, r = 3/2: errors are 1.8e-50 to 3.7e-50, and interior mass is 2.0e-50 to 4.0e-50. For i = 1, 3, 7 the flipped
  formula is off by 0.330, 0.674 and 0.674. The neutral i/N is off by 0.239, 0.416 and 0.258.
- Moran, r = 4/5: errors are about 5e-42 to 2e-41.
- Moran, r = 1: errors are about 3e-37 to 7e-37, and the value matches i/N.

## 6. Limitations and honest scope

- **Model scope.** This is the Moran variant in which the individual that dies is drawn uniformly from all N before
  the offspring is added: the task's `p_i`, `q_i`. Other variants have different rates. Hypotheses are `0 < N` and
  `0 < r`.
- **Not formalized.**
  - No path measure: see §2.
  - No time-to-absorption results.
  - Nothing on Kimura's diffusion limit.
  - Nothing on the Price equation or Fisher's general theorem; those are later stages of this lane.
- **Uniqueness route.** For Moran, uniqueness comes from the absorption limit (step 5), not from a maximum
  principle. The kernel accepts either route; the task suggested the maximum principle.
- **Literature (Tier L).** Read through alphaXiv this stage, and pinned:
  - P. Keller and M. Uğurlu, "Fixation probability in Moran-like Processes on graphs", arXiv:2403.12598v2.
    - **Eq. (3), p. 2:** `p⁺(j) = r/(1+(r−1)j/n)·(j/n)(n−j)/n` and `p⁻(j) = 1/(1+(r−1)j/n)·((n−j)/n)(j/n)`. "The
      individual to be replaced is chosen uniformly at random" (p. 2). These are algebraically identical to our
      `moranUp`/`moranDown`, since `r/(1+(r−1)j/n)·j/n = r j/(r j + n − j)`.
    - **Eq. (4), p. 4:** the Moran fixation probabilities `ϱᵢ = i/n` (r = 1) and `(1 − 1/rⁱ)/(1 − 1/rⁿ)` (r ≠ 1).
      §1.1.2 (p. 3) derives them by optional stopping applied to a martingale.
    - Their proof uses optional stopping on a path measure. Ours avoids a path measure; it works on matrix powers.
    - **Source typos** (do not copy them):
      - p. 3 says "For r = 0" where r = 1 is meant.
      - Thm 2.1 on p. 7 prints `(1−1/rⁱ)/(1−1/r)`, missing the ⁿ. Eq. (4) and eq. (23) on p. 10 have it right.
    - They cite Moran (1958), Proc. Camb. Phil. Soc. 54, 60–71, which was **not** read directly. Cite it only "as cited
      in" Keller–Uğurlu, unless the paper stage pins it.
  - E. Marion, "A Formalization of the Ionescu-Tulcea Theorem in Mathlib", arXiv:2506.18616v5.
    - It builds the trajectory measure of a Markov chain in Mathlib: `ProbabilityTheory.Kernel.traj` and
      `partialTraj`, pp. 8–12.
    - The pinned Mathlib contains it: `Mathlib/Probability/Kernel/IonescuTulcea/Traj.lean` exists on disk.
    - The paper does not treat hitting or absorption probabilities.
  - **Novelty.** Not checked beyond these two searches. No Lean or Mathlib formalization of the Wright–Fisher or
    Moran fixation probability was found in them, but that is not a proof of absence.
- **Producer ≠ verifier.**
  - An independent recompile is required before Tier A wording.
  - So is an audit of the statements, especially definition fidelity of `moranUp`/`moranDown`/`bdP` against the
    task formula, and the Fin-indexing of `bdP`.
  - No human audit exists (cf. the `LEDGER_UNAUDITED_TIER_A` flag in earlier papers).
- **Compiler warnings.** Only deprecation warnings (`if_pos`/`if_neg`/`if_true`); no errors.

## 7. Suggested next steps for this lane

- **VERIFY stage.** Independent recompile and statement audit.
- **Optional bridge upgrade (stretch, not started).** For finite discrete state spaces, prove that
  `kPow P n i j = (Kernel.traj κ 0 (δᵢ)) {x | x n = j}`, using the pinned Mathlib's `ProbabilityTheory.Kernel.traj`.
  Then show that `ℙᵢ(∃ n, Xₙ = N)` is the limit. This would move the §2 bridge into Tier A.
- **PAPER stage.** Pin the citations and build the ledger with `papers/ar_*/ledger/build_ledger.py`.
- **Price equation (preregistered prior 0.8).** It is a finite-sum identity and can reuse nothing here except the
  `Fin` sum lemmas.
- **LeanFlow / rusty-SUNDIALS.** A possible contribution is a fail-closed "Lean data-certificate" check that the
  Tier X numerics (`numerics_fixation.json`) agree with the kernel-checked closed form at the sampled
  (N, r, i). It was not attempted in WORK1.

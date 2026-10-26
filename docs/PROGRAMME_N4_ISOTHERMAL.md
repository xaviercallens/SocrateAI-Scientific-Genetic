# Programme N4: the isothermal theorem of evolutionary graph theory, kernel-checked (2026-10-26)

Paper: doi:10.5281/zenodo.23660071. It is an Elenchus-tiered preprint, independently audited by a separate AI agent
(ACCEPT_WITH_CORRECTIONS, corrections applied). It has had no human review.

**Result (Tier A).** Take the birth–death Moran process on a strongly connected weighted graph W, with W row-stochastic
and column sums 1 (isothermal).
- The probability of fixation from any mutant set S, defined as the limit of the t-step kernel, equals the
  well-mixed Moran value moranFix(n, r, |S|).
- For a single mutant at any node this is (1 − r⁻¹)/(1 − r⁻ⁿ); the neutral case gives |S|/n.
- `stationary_fixation` generalizes this to selection weights μ with μW = μ.

**How it is proved.** The per-configuration forward bias satisfies P⁺(S) = r·P⁻(S) (`forward_bias`). From that, the
proof gets harmonicity of S ↦ moranFix(|S|), and then the absorption-limit machinery of programme P7, reused unchanged.

**What the preregistration got wrong, and the known gap.** The mutant count on an isothermal graph is *not* strongly
lumpable: `cycle4_not_lumpable` shows the 4-cycle values 1/3 vs 2/3. This is the gap in the 2005 projection argument
noted by Keller & Uğurlu (arXiv:2403.12598). Pattni et al. (arXiv:1505.06021, Lemma 1) give a valid route.

**Controls:**
- the star graph (not isothermal), with exact fixation values 5/12, 2/3 and 7/12 against Moran's 4/7;
- the disconnected identity matrix;
- planted bad proofs, all caught.

**Gate.** The audit used three checks: the axiom list, the toolchain's `leanchecker` kernel replay (rc 0), and a grep
for `skipKernelTC`/`addDecl`/`native_decide`/`implemented_by`. `#print axioms` alone can be fooled. A forged control was
rejected by the replay.

**Layout.** `formal/Programmes/EvoGraph/*.lean` imports `Programmes.PopGen.*`, which comes from PR #2. **This PR is stacked
on PR #2**, so merge #2 first.

**Scope:** sufficiency only, birth–death updating only, fixation as the limit of the kernel. Novelty is not claimed: this
is a formalization.

Records are in `data/programme_n4_2026-10-26/`, with a MANIFEST of sha256 hashes. The paper is in
`papers/programme_n4_isothermal/`.

# P7 WORK2: the Price equation and the haploid selection increment (Fisher-type)

Lane P7, stage WORK2, 2026-10-25. Producer: Claude Opus 5.5 subagent. **Producer ≠ verifier is still pending**:
nothing below is final Tier A until an independent recompile and statement audit is done.

Preregistration: `results/programmes/2026-10-25/preregistration.json`, lane P7, priors `price` 0.8 and
`fisher_general` 0.35. Outcome against the preregistration:

| Prereg item | Outcome |
|---|---|
| `price` (Price equation) | **Done** on the first budget: exact identity, kernel-checked (Tier A, producer side). |
| `fisher_general` (Fisher's general theorem from Price) | **Not achieved.** Only a special case is proved: the haploid/asexual selection increment with constant type fitnesses and no transmission bias (§3). The general reading (Price 1972 / Ewens "partial change", additive genetic variance) is not formalised. |

Both results formalise **known classical identities** (Price 1970; Price 1972; see §6). **Novelty is not claimed.**
A brief alphaXiv search found no prior proof-assistant formalization of the Price equation; that search was too
small to support a priority claim, so none is made. Nothing here depends on openai/math results; only Mathlib.

## 1. What is kernel-checked

Pinned openai/math Lean 4.34.1 + Mathlib through `scripts/certified_numerics/lean_pinned.py` (cores 4-7). Every
claimed declaration prints `[propext, Classical.choice, Quot.sound]` (verbatim in the compile JSONs). The claimed
files contain no `sorry`, `admit`, `native_decide` or `axiom` token (grep, whole-word: no match).

| File (`formal_cert/Programmes/PopGen/`) | Module | sha256 | Compile record |
|---|---|---|---|
| `Price.lean` | `Programmes.PopGen.Price` | `1309c57a7818da57c9e955764df28fb26fb83c930661bd7c8439994d4d07438f` | `compile_Price.json` (rc 0, 7.1 s, 21 axiom lines, all whitelisted, no warnings) |
| `PriceImportAll.lean` | (check only) | `3e82333747ae4bddb921a4755df83ae27a3dba64b1d555c677ff814fe8f76947` | `compile_price_import_all.json` (rc 0) |
| `PriceControls.lean` | (planted, never imported) | `93bc16d08561f582560849efdd2cc3da493daf20fbf3a3993209c5f15168afc4` | `compile_price_controls.json` (rc 0, controls caught) |

The WORK1 files and their compile JSONs were not modified.

## 2. The Price equation (`ProgPopGen.Price`, over an arbitrary field `K`)

**Model.** A finite index set `s : Finset ι` of types; abundances `n i`, fitness `w i` (offspring per individual),
parental trait `z i`, mean offspring trait `z' i`, all in a field `K`. Type `i` contributes `n i · w i` offspring of
mean trait `z' i`. Definitions:

- `wmean s n f = (∑ n f)/(∑ n)` (weighted population mean),
- `wcov s n f g = (∑ n (f − f̄)(g − ḡ))/(∑ n)` (population covariance, not the sample one), `wvar s n f = wcov s n f f`,
- `zbarNext s n w z' = wmean s (n·w) z' = (∑ n w z')/(∑ n w)` (offspring-generation mean trait).

The offspring-generation definition is the modelling assumption. Given it, Price's equation is an algebraic
identity, not a statement about any particular mechanism.

**`price_equation`.** If `∑ n ≠ 0` and `∑ n w ≠ 0`:

    w̄ · (z̄' − z̄) = Cov_n(w, z) + E_n(w · Δz),    Δz i = z' i − z i.

**Proof route.** First expand both numerators into five raw sums (`cov_num_expand`, `trans_num_expand`, using
pointwise `ring` and `Finset.sum_*_distrib`/`mul_sum`). Then apply the scalar identity `price_scalar`
(`field_simp; ring`).

Also kernel-checked:
- `wcov_eq_moment`: `Cov = E(fg) − E f · E g`.
- `wmean_ne_zero`.
- `price_equation_pos`: the same identity over ℝ under the biological hypotheses `s.Nonempty`, `n i > 0`, `w i > 0`.
  Nonemptiness is necessary for these hypotheses to imply the nonzero totals.

## 3. Fisher-type corollary: which version this is

Set the trait to be fitness, `z = w`, and assume **no transmission bias**, `z' = w`: offspring of type `i` have the
parental type fitness. This covers haploid/asexual (clonal) types with constant type fitnesses and no mutation. The
transmission term is then identically `0` (`w i (w i − w i) = 0`). With `wbarNext s n w = zbarNext s n w w`:

- `fisher_haploid`: `w̄ (w̄' − w̄) = Var_n(w)` (field `K`, hypotheses `∑ n ≠ 0`, `∑ n w ≠ 0`).
- `fisher_haploid_div`: **`Δw̄ = Var_n(w) / w̄`**.
- `wvar_nonneg`, `fisher_haploid_nondecreasing` (ℝ): if `n ≥ 0`, `∑ n > 0` and `∑ n w > 0`, then `w̄ ≤ w̄'`.
- `fisher_haploid_pos` (ℝ): both conclusions under `s.Nonempty`, `n i > 0`, `w i > 0`.

**Scope, stated plainly.** This is the discrete-generation selection increment for haploid/asexual types with
constant fitnesses. It is the discrete analogue of Baez's remark that `d w̄/dt = Var(f)` when the fitness functions
are constant (arXiv:2107.05610v3, p. 4). It is also LaGatta's eq. (2.27) with the environmental term set to zero
(§6). It is **not**:
- Fisher's fundamental theorem in the Price (1972) / Ewens / Edwards "partial change" reading;
- a statement about additive genetic variance (here the variance is the total variance in type fitness, which
  equals the additive variance only trivially, for clonal types);
- the diploid multi-allele case, which is not this equality and is not formalised here. H12's one-locus two-allele
  diploid result, `AutoresearchPopGen.wbar_increment`, has a different exact increment.

Frequency-dependent fitness is also outside scope, because `w` is fixed across the generation.

## 4. Controls

**Positive control, kernel-checked** (each side evaluated by `simp`/`norm_num`, not by the theorem):
- `ctl_pos_three_types`. Three types, `n = (1,2,3)`, `w = (1,2,4)`, `z = (0,1,3)`, `z' = (1,1,2)`. Both sides of the
  Price equation evaluate to `−13/36`.

**Negative controls, kernel-checked:**
- `ctl_neg_transmission_bias`. `n = (1,1)`, `w = z = (1,2)`, but offspring fitness `z' = (1,1)` (transmission bias
  for fitness). Then `Δw̄ = −1/2` while `Var(w)/w̄ = 1/6`, so **the variance formula fails**. Price's equation still
  holds on the same instance: both sides are `−3/4`, with `Cov = 1/4` and `E(wΔz) = −1`. The no-transmission-bias
  hypothesis is therefore what the Fisher-type corollary rests on.
- `ctl_neg_missing_wbar`. `n = (1,1)`, `w = (1,2)`, no bias. `Δw̄ = 1/6 = Var/w̄`, but `Var = 1/4 ≠ Δw̄`. The `1/w̄`
  factor cannot be dropped.

**Non-vacuity, kernel-checked:**
- `nonvacuous_pos`: the positivity hypotheses hold for a concrete population.
- `ProgPopGenPriceCheck.instance_four_types` (in `PriceImportAll.lean`): the monotonicity theorem applied to
  `n = (1,1,1,1)`, `w = (1,2,3,4)`. The same file shows that the Price module imports together with H11, H12 and the
  WORK1 fixation modules without name clashes.

**Planted controls, caught** (`PriceControls.lean`; both claim false statements; verbatim output):

    'ProgPopGenPriceControls.ctl_sorry' depends on axioms: [propext, sorryAx, Classical.choice, Quot.sound]
    'ProgPopGenPriceControls.ctl_axiom' depends on axioms: [propext, Classical.choice, ProgPopGenPriceControls.ctl_cheat, Quot.sound]

- `ctl_sorry` claims `Δw̄ = Var(w)`.
- `ctl_axiom` claims the variance formula under transmission bias.

## 5. Crosswalk to the Genetic repo (`ANSE.Alix`), kernel-checked here

In namespace `ProgPopGen.Price.Alix`, the definitions `meanFit`, `selectP` of
`SocrateAI-Scientific-Genetic/formal/ANSE/AlixPopulationGenetics.lean` are copied verbatim, over ℚ. Two theorems
restate their statements with identical text:
- `aurelia_mean_fitness_increment`,
- `kal_increment_is_variance`.

Both use `s = univ : Finset (Fin 2)`, `n = ![p, 1−p]`, `w = ![w₁, w₂]`, with bridge lemmas `wmean_two`,
`wbarNext_two` and `wvar_two`.
- `aurelia_mean_fitness_increment` is **derived from the general multi-type `fisher_haploid_div`**.
- `kal_increment_is_variance` is only the general `wvar` definition evaluated at `Fin 2`. It is not a consequence of
  Price.
- The third Alix H12 theorem, `aurelia_mean_fitness_nondecreasing`, is not restated. Its multi-type generalisation is
  `fisher_haploid_pos` over ℝ.

The Genetic file itself is not imported. As a reference, it was compiled read-only under the same pin
(`compile_genetic_alix_reference.json`, rc 0, 11 declarations, all `[propext, Classical.choice, Quot.sound]`). The
equality of statements is textual (same definitions, same theorem text, different namespace). No Lean term relates
the two namespaces.

## 6. Literature (Tier L, pinned via alphaXiv)

- T. LaGatta, *The Mathematics of Evolution: The Price Equation, Natural Selection, and Environmental Change*,
  arXiv:2202.10289v4:
  - Discrete Price equation: eq. (1.6), p. 9. It is written for relative fitness `U = W/W̄`; multiplying by `W̄`
    gives our form.
  - Fisher in Price's reading, `∂_NS(U) = var(U)`: eq. (1.8), p. 9.
  - Aggregate form `∂_NS(W) = var(W)/W̄`: eq. (2.24), p. 18.
  - Two-generation identity `W̄′ − W̄ = var(W)/W̄ + E[Δ_w(W,W′)W]/W̄`: eq. (2.27), p. 19. Our `fisher_haploid_div` is
    this with the environmental term zero.
  - p. 4 quotes Price (1972, p. 130): the theorem concerns the change *due to natural selection*, not the total
    change.
- J. C. Baez, *The Fundamental Theorem of Natural Selection*, arXiv:2107.05610v3:
  - Theorem 1, p. 3 (continuous time): `Σ ṗᵢ fᵢ = Var(f)`.
  - p. 4: with constant fitness functions, `d f̄/dt = Var(f)`.
  - pp. 1–3: the Price/Ewens/Edwards partial-change reading, and Ewens' remark that the general theorem "does not
    provide the substantial biological statement that Fisher claimed".
- G. R. Price, *Selection and covariance*, Nature 227 (1970) 520–21 ([Pri70]); G. R. Price, *Fisher's
  'fundamental theorem' made clear*, Ann. Hum. Genet. 36(2) (1972) 129–140 ([Pri72a]). Bibliographic details are as
  listed in LaGatta's reference list (arXiv:2202.10289v4, p. 69); Baez's ref. [14] agrees for Price 1972. The
  originals were not read here, so only these bibliographic entries are verified, not their content.

## 7. Tier X cross-check (never supports a claim)

`price_numerics.py` → `price_numerics.json` (exact `fractions.Fraction`, seeded):
- Lean control instances: `−13/36 = −13/36`; transmission bias `Δw̄ = −1/2`, `Var/w̄ = 1/6`, Price `−3/4 = −3/4`;
  missing factor `Δw̄ = 1/6`, `Var = 1/4`.
- 500 random exact instances (1–8 types): Price equal 500/500, haploid increment equal 500/500.
- float64 on 500 random instances: worst absolute residual `1.42e-14` (rounding only).

## 8. Status and open items

- Tier A on the producer side only. An independent recompile and statement audit is pending.
- `fisher_general` is not achieved; a later budget could formalise the partial-change reading with an explicit
  regression (average-effect) decomposition. Out of scope here.
- Genetic-repo contribution: note at `GENETIC_CONTRIB_NOTE.md`. The PR is the parent session's job.

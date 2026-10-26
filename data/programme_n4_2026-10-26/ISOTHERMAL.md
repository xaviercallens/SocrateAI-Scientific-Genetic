# N4: the isothermal theorem of evolutionary graph theory in Lean 4 (proposal 3, 2026-10-26)

Status: **producer result, not yet audited** (producer ≠ verifier: an independent audit is still owed).
Tier A means the Lean statements below are kernel-checked under pinned Lean 4.34.1 + Mathlib d13f23b7, with
axioms `{propext, Classical.choice, Quot.sound}` only. It certifies the Lean statements, not their biological
reading. Known results (Lieberman, Hauert & Nowak 2005; Pattni, Broom, Rychtář & Silvers 2015; Keller &
Uğurlu 2024). Novelty is not claimed for the mathematics. We are not aware of a prior machine-checked formalization
(no systematic search of Lean/Coq/Isabelle libraries was done). The cycle-4 control formalizes the
non-Markovianity of the count pointed out by Keller–Uğurlu; it is not a new discovery of that gap.

## Literature check (alphaXiv, primary sources read)
- Pattni–Broom–Rychtář–Silvers, arXiv:1505.06021, §2–3, App. A.2/A.4. For BD with fitness-proportional
  birth (their "BDB"), Prop. 2 says that W is ρ-equivalent to Moran for all r iff f_R(W) is a circulation;
  for right-stochastic W that means doubly stochastic, i.e. the isothermal theorem. They assume strong
  connectivity and allow self-loops `w_ii > 0`. Their proof goes through the *jump chain* (holding removed),
  which is lumpable onto ±1 steps with probabilities r/(1+r) and 1/(1+r).
- Keller–Uğurlu, arXiv:2403.12598, §1.2 and §5. LHN's argument projects the configuration chain onto
  the mutant count and treats it as a birth–death chain, but "the projection does not need to be
  Markovian". They prove the theorem with the martingale r^{-M_k} and optional stopping, and generalize it
  to stationary selection (μ W = μ, Thm 3.1).
- The gap is in LHN's projection argument only. Pattni et al.'s jump-chain lumpability argument
  (their Lemma 1) is a valid proof, so the theorem itself was not unproven before Keller–Uğurlu.
- This file follows the configuration-level harmonic-function / martingale route. It does not use
  lumpability of the count.

## Model (`formal_cert/Programmes/EvoGraph/BirthDeath.lean`)
- Nodes `Fin n`. `W : Fin n → Fin n → ℝ` is row-stochastic (`RowStochastic`), self-loops are allowed.
  Selection weights are `μ i > 0`; μ ≡ 1 is LHN's uniform birth–death process. Mutant fitness is r > 0,
  resident fitness is 1.
- The state is `S : Finset (Fin n)` (the set of mutants). Node i reproduces with probability
  μ_i f_S(i)/F(S), where F(S) = Σ_k μ_k f_S(k). Its offspring replaces j with probability W_ij, and the new
  state is `upd S i j` (`insert j S` if i ∈ S, else `S.erase j`).
- `bdK μ W r S T = Σ_i Σ_j coef(S,i,j)·[upd S i j = T]`. `sPow` is its t-step kernel.
- `Isothermal W` means row-stochastic and column sums equal to 1 (doubly stochastic, so every temperature is 1).
- `StronglyConnected W` means that for all a and b, `Relation.TransGen (W · · > 0) a b`.

## Theorem statements (all Tier A)
| Lean name | statement |
|---|---|
| `two_absorbing_limit` (Transport.lean) | stochastic kernel on a finite type with absorbing `bot ≠ top`; h harmonic off {bot,top}, h bot = 0, h top = 1, positive one-step variance off {bot,top} ⇒ Pᵗ(s,top) → h s and Pᵗ(s,bot) → 1 − h s. Obtained by conjugating with an enumeration σ ≃ Fin (M+1) (∅ ↦ 0, univ ↦ last) and applying **P7's `ProgPopGen.fixation_limit` verbatim**. Nothing is re-proved. |
| `bd_absorbing` | `bdK` is stochastic, and ∅ and univ are absorbing (W only row-stochastic) |
| `up_eq`, `down_eq` | P⁺(S) = r·cutOut(S)/F(S), P⁻(S) = cutIn(S)/F(S) |
| `cut_balance` | stationary μ (Σ_i μ_i W_ij = μ_j) ⇒ cutOut(S) = cutIn(S) for every S |
| `forward_bias` | stationary μ ⇒ **P⁺(S) = r·P⁻(S) for every configuration S** |
| `moranFix_card_harmonic` | S ↦ moranFix n r \|S\| (P7's fixation function) is harmonic for the configuration chain at every S ≠ ∅ |
| `moranFix_card_var_pos` | positive one-step variance whenever some edge i∈S → j∉S has μ_i W_ij > 0 |
| `cutPositive_of_strongly_connected` | strong connectivity ⇒ every proper nonempty S has an edge leaving it |
| `stationary_fixation` | n > 0, r > 0, μ > 0, W row-stochastic, μW = μ, cut positive ⇒ Pᵗ(S,univ) → moranFix n r \|S\| and Pᵗ(S,∅) → 1 − moranFix n r \|S\| (Keller–Uğurlu Thm 3.1) |
| `isothermal_fixation` | **isothermal theorem**: μ ≡ 1, W isothermal and strongly connected ⇒ Pᵗ(S,univ) → moranFix n r \|S\| for every initial configuration S |
| `isothermal_fixation_single` | r ≠ 1: a single mutant at **any** node v ⇒ limit (1 − r⁻¹)/(1 − r⁻ⁿ) |
| `isothermal_fixation_neutral` | r = 1 ⇒ limit \|S\|/n |
| `isothermal_fixation_uniform` | a single mutant at a uniformly random node, r ≠ 1 ⇒ (1/n)Σ_v Pᵗ({v},univ) → (1 − r⁻¹)/(1 − r⁻ⁿ) |

The sequence Pᵗ(S, univ) is non-decreasing in t because univ is absorbing (`sPow_top_mono` in Transport.lean,
`bd_fixation_mono` in BirthDeath.lean, both Tier A), so its limit is the probability of ever reaching univ.

## Statement precision: configuration process versus count process (deviation from the preregistered wording)
The preregistration says that the number of mutants "is a birth–death chain" with step ratio r. **That is
false as a Markov-chain statement**, and the false part is not formalized. Kernel-checked
`cycle4_not_lumpable` shows it: on the isothermal 4-cycle (W_{i,i±1} = 1/2) with r = 2,
P⁺({0,1}) = 1/3 but P⁺({0,2}) = 2/3, although both configurations have two mutants. So the transition
probabilities of the count are not functions of the count. What *is* true, and what is proved, is the
per-configuration forward bias P⁺(S) = r·P⁻(S) (`forward_bias`), which holds in both of those
configurations. That is enough for `moranFix n r |S|` to be harmonic on configurations, and the absorption limit
follows from it. This is the gap in the LHN proof identified by Keller–Uğurlu (2024). The independent
numeric twin `numerics_evograph.json → cycle4_up_down_two_mutants` gives (1/3, 1/6) and (2/3, 1/3).

## Controls
Real controls are theorems in `Examples.lean`. All use r = 2 and μ ≡ 1 unless noted.
- **Positive: complete graph equals well-mixed Moran.** `wellMixed_up` and `wellMixed_down` hold for every n and
  every S: for W_H = J/n, P⁺(S) = moranUp n r |S| and P⁻(S) = moranDown n r |S|. Here the configuration chain
  lumps exactly onto P7's Moran chain. `wellMixed3_fixation`: Pᵗ({0},univ) → 4/7, which is P7's value.
- **Positive: stationary selection.** `starStationary_fixation`: on the star with μ = (1/2, 1/4, 1/4), the limit
  from the centre is 4/7, the Moran value.
- **Negative: star n = 3** (centre 0, W_0j = 1/2, W_j0 = 1), which is not isothermal:
  - `star_not_isothermal`: the centre temperature is 2.
  - `star_bias`: P⁺({0}) = P⁻({0}) = 1/2, so the ratio is 1 ≠ r.
  - Exact fixation limits, through the same transported P7 machinery with the explicit harmonic table
    (0, 5/12, 2/3, 2/3, 5/6, 5/6, 14/15, 1): `star_fixation_centre` gives **5/12**, `star_fixation_leaf`
    gives **2/3**, and `star_fixation_uniform` gives **7/12**.
  - `star_ne_moran`: 5/12, 7/12 and 2/3 are each ≠ 4/7 = moranFix 3 2 1. The uniform single mutant on the
    star (7/12 > 4/7) is amplified, as expected for a star under BD.
- **Negative: the connectivity hypothesis is necessary.** `identity_no_fixation`: W = I is isothermal, but
  Pᵗ({0},univ) = 0 for all t, so the limit is not moranFix 2 2 1 = 2/3.
- **Negative: no lumpability.** `cycle4_not_lumpable`, described above.
- **Planted** (`EvoGraphPlantedControls.lean`, imported by nothing). Each states a false claim, and all three
  are CAUGHT:
  - `ctl_sorry` (star centre → 4/7) reports `sorryAx`.
  - `ctl_axiom` (star uniform → 4/7) reports `EvoGraphControls.ctl_cheat`.
  - `ctl_axiom_lumpable` (P⁺({0,1}) = P⁺({0,2}) on the cycle) reports `ctl_cheat`.
- **Independent numerics** (`numerics_evograph.py`, exact rationals, Gaussian elimination on all 2ⁿ states):
  - W_H with n = 3 and n = 4, the 4-cycle, and a random rational doubly stochastic 4×4 matrix: every
    configuration equals the Moran value (4/7, 8/15).
  - Negative cases differ: star n = 3 (5/12, 2/3, 2/3, average 7/12), star n = 4 (392/1265, 168/253, average
    728/1265 ≠ 8/15), and a random row-stochastic, non-doubly-stochastic matrix.
  - The star with stationary μ gives 4/7. The Lean star table matches the numerics exactly.

## `#print axioms` (verbatim, `compile_ImportAll.json`, rc 0)
```
'EvoGraph.two_absorbing_limit' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.sPow_top_mono' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.bd_fixation_mono' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.forward_bias' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.cut_balance' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.moranFix_card_harmonic' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.stationary_fixation' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.isothermal_fixation' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.isothermal_fixation_single' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.isothermal_fixation_neutral' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.isothermal_fixation_uniform' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.wellMixed_up' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.wellMixed_down' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.wellMixed3_fixation' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.star_not_isothermal' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.star_bias' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.star_fixation_centre' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.star_fixation_leaf' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.star_fixation_uniform' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.star_ne_moran' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.starStationary_fixation' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.cycle4_not_lumpable' depends on axioms: [propext, Classical.choice, Quot.sound]
'EvoGraph.identity_no_fixation' depends on axioms: [propext, Classical.choice, Quot.sound]
```
Planted (`compile_PlantedControls.json`):
```
'EvoGraphControls.ctl_sorry' depends on axioms: [propext, sorryAx, Classical.choice, Quot.sound]
'EvoGraphControls.ctl_axiom' depends on axioms: [propext, Classical.choice, EvoGraphControls.ctl_cheat, Quot.sound]
'EvoGraphControls.ctl_axiom_lumpable' depends on axioms: [propext, Classical.choice, EvoGraphControls.ctl_cheat, Quot.sound]
```

## Proved versus stated (scope)
- **Proved (Tier A).** Everything in the theorem table and the controls.
  - The absorption-limit step is P7's kernel-checked `fixation_limit`, reused through conjugation. P7's
    `moranFix`, `moranFix_step`, `moranFix_closed`, `moranFix_neutral`, `geomS_pos`, `moranUp` and
    `moranDown` are reused unchanged.
  - P7 files are unmodified: `Absorption.lean` sha256 86fb5168…, `MoranModel.lean` b0822dc0….
- **Limits of the statement:**
  - "Fixation probability" means lim_t Pᵗ(S,univ) in the n-step kernel sense (as in P7). There is no
    measure-theoretic path space.
  - Strong connectivity ⇒ cut positivity is proved. The converse (equivalence) is only stated in prose and
    is not needed.
  - Only the sufficiency direction is proved (isothermal ⇒ Moran). The converse "Moran for all r ⇒
    isothermal" (Pattni et al. Prop. 2, (b) ⇒ (a)) is not formalized.
  - BD with fitness at birth only. Death–birth and other dynamics are not formalized.
  - Exact star values are formalized for n = 3 only. The n = 4 star values are numerics only (Tier C).
- **Not claimed:** that the count process is Markov. It is refuted on an isothermal graph.

## Files and builds
- Lean files (pinned build, `taskset -c 4,5`, every rc 0):
  - `formal_cert/Programmes/EvoGraph/Transport.lean`, 5.7 s, olean `Programmes.EvoGraph.Transport`
  - `formal_cert/Programmes/EvoGraph/BirthDeath.lean`, 8.1 s
  - `formal_cert/Programmes/EvoGraph/Examples.lean`, 12.6 s
  - `formal_cert/Programmes/EvoGraph/EvoGraphImportAll.lean`
  - `formal_cert/Programmes/EvoGraph/EvoGraphPlantedControls.lean`
- Compile records: `compile_Transport.json`, `compile_BirthDeath.json`, `compile_Examples.json`,
  `compile_ImportAll.json`, `compile_PlantedControls.json`.
- A sorry/axiom/native_decide grep over the three claimed files finds only a docstring mention.

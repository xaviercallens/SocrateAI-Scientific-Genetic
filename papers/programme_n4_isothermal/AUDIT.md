# N4 audit: the isothermal theorem (birth–death) kernel-checked in Lean 4

Auditor: an independent adversarial AI auditor (Claude agent). It did not produce the paper, the Lean files or the ledger (producer != verifier). This is **not a human audit**. Date: 2026-10-26.
Scope: `papers/p3_isothermal_theorem/` (main.tex/pdf, ledger, publication_meta.json, gen/), `results/proposal3/2026-10-26/N4/`, `formal_cert/Programmes/EvoGraph/*.lean`, and the P7 modules it reuses (`Programmes/PopGen/Absorption.lean`, `MoranModel.lean`; `MoranFixation.lean` is **not** imported by EvoGraph). Preregistration: `results/proposal3/2026-10-26/preregistration.json`, item N4.
Auditor scratch (not durable): `/tmp/n4audit/` (own olean dir, compile JSONs, audit Lean files, numerics script).

## VERDICT: ACCEPT_WITH_CORRECTIONS

The Lean results are what the paper says they are. I recompiled them from source into my own olean directory with no producer oleans, ran my own `#print axioms`, ran a kernel replay with `leanchecker`, planted `sorry` and a `debug.skipKernelTC` forgery (both caught), and recomputed every quoted number with my own exact-rational script. The statements model the Lieberman–Hauert–Nowak birth–death process correctly, and the deviation from the preregistration is reported honestly. The two blocking items are both about stale review-status text and how the soundness gate is described. Neither is a mathematical defect.

## 1. Kernel checks (all run by the auditor)

Toolchain: Lean 4.34.1 (`.../elan/toolchains/leanprover--lean4---v4.34.1`), Mathlib `d13f23b723b8a846827a245b89c10fc7d3f11612` (pinned `lake-manifest.json`). The olean directory in `scripts/certified_numerics/lean_pinned.py` is hard-coded (`OUT = formal_cert/.olean_out`); there is no CLI option for it. I ran the unmodified script through a wrapper (`/tmp/n4audit/mycompile.py`) that sets `lean_pinned.OUT = /tmp/n4audit/olean`. The wrapper printed the `LEAN_PATH` entries that contain `Programmes/`, and the only one was `/tmp/n4audit/olean`, so no producer olean was visible.

| step | module | rc |
|---|---|---|
| P7 | `Programmes.PopGen.Absorption` | 0 |
| P7 | `Programmes.PopGen.MoranModel` | 0 |
| N4 | `Programmes.EvoGraph.Transport` | 0 |
| N4 | `Programmes.EvoGraph.BirthDeath` | 0 |
| N4 | `Programmes.EvoGraph.Examples` | 0 |
| N4 | `EvoGraphImportAll.lean` (no olean) | 0; 23/23 lines `[propext, Classical.choice, Quot.sound]` |
| N4 | `EvoGraphPlantedControls.lean` | 0; `ctl_sorry` → `sorryAx`; `ctl_axiom`, `ctl_axiom_lumpable` → `EvoGraphControls.ctl_cheat` (caught) |

I ran `#print axioms` myself (`/tmp/n4audit/src/AuditN4.lean`, rc 0) on `isothermal_fixation`, `isothermal_fixation_single`, `stationary_fixation`, `forward_bias`, `cycle4_not_lumpable`, `star_ne_moran`, `identity_no_fixation`, `star_fixation_{centre,leaf,uniform}`, `isothermal_fixation_{uniform,neutral}` and `two_absorbing_limit`. Every one reports `[propext, Classical.choice, Quot.sound]`.

The same file also contains three checks of my own:
- **Restatement.** `isothermal_fixation_single`, restated in my own typing with the closed form `(1 - 1/r)/(1 - 1/r^n)` and the hypotheses spelled out (non-negativity, row sums, column sums, `TransGen` connectivity), is closed by the producer theorem.
- **Unfolding.** `bdK μ W r S T` is checked by `rfl` against the explicit double sum `∑ i, ∑ j, μ i f_S(i)/(∑ k μ k f_S(k)) · W i j · [upd = T]`.
- **Non-vacuity.** `Isothermal (WH 3) ∧ StronglyConnected (WH 3)` holds, so the hypotheses can be satisfied.

**Planted sorry** (`/tmp/n4audit/src/AuditPlanted.lean`, rc 0). This file states two false claims: that the star centre fixes with probability `moranFix 3 2 1`, and that the star has forward bias `r`. Both report `sorryAx`, so they were caught, and the exit code was 0. As the paper says, rc is not the gate.

**Forbidden-construct grep** over EvoGraph/*.lean, Absorption.lean and MoranModel.lean for `skipKernelTC|addDecl|debug|native_decide|implemented_by|#exit|unsafe|extern|set_option|macro|elab|axiom|opaque|sorry|admit|run_cmd|ofReduceBool|@[`. The only hits are the documented planted controls in `EvoGraphPlantedControls.lean`, which nothing imports, and comment text. `simp +decide` and `decide` use kernel reduction, not `native_decide`.

**Kernel replay.** I ran `leanchecker Programmes.PopGen.Absorption Programmes.PopGen.MoranModel Programmes.EvoGraph.Transport Programmes.EvoGraph.BirthDeath Programmes.EvoGraph.Examples` on my oleans, and it returned **rc 0**.

**Negative control for the replay** (`/tmp/n4audit/src/AuditForged.lean`). I added `AuditForged.forged_false : False := True.intro` under `set_option debug.skipKernelTC true` and used it to "prove" that the star centre fixes with 4/7. `#print axioms` reports **only the whitelist** `[propext, Classical.choice, Quot.sound]`, so the bypass is real. `leanchecker AuditForged` exits **1** with "(kernel) declaration type mismatch, 'AuditForged.forged_false' has type True but it is expected to have type False". So the replay check is not vacuous, and the replayed N4 modules contain no such forgery.

Limitation: the replay was not `--fresh`. Mathlib and core were loaded from the pinned build and not re-checked.

**Hashes.** These match `publication_meta.json`, `result.json` and the paper's Table 3. A diff of commit `8ff5912` against HEAD, restricted to EvoGraph/ and PopGen/, is empty, and the working tree has no changes there.

**Listings.** All 45 statement and definition blocks displayed in `gen/*.lean` (proofs elided) occur verbatim, after whitespace normalisation, in the hashed sources. Checked with `/tmp/n4audit/check_listings.py`: 0 mismatches. The statement lock is in `results/proposal3/2026-10-26/N4/STATEMENT_LOCK.json`: file sha256s, 17 headline statements and 18 definitions copied verbatim.

## 2. Statement adequacy

- **`bdK` is the LHN birth–death process on a weighted digraph.** Node `i` reproduces with probability `μ_i f_S(i)/F(S)`, where `f = r` for mutants and `1` for residents. The offspring replaces `j` with probability `W i j`, and `j = i` is allowed. The new state is `insert j S` if `i ∈ S` and `S.erase j` otherwise. With `μ ≡ 1` this is BDB in Pattni et al. Table 2, `b_i = F_i/ΣF`, `d_ij = w_ij` for right-stochastic W. It also matches Keller–Uğurlu Def. 2.1, with `π` replaced by an unnormalised positive `μ`. Checked by `rfl` (§1), and the star and cycle transition numbers were reproduced independently (§3). `upP` and `downP` are the total probabilities of moving to a configuration with `|T| = |S| ± 1`. The only states that can be reached are `S`, `S ∪ {j}` and `S \ {j}`, so these are exactly `P⁺_S` and `P⁻_S` in Pattni's sense.
- **`Isothermal` means row-stochastic, non-negative, and every column sum equal to 1.** LHN and Keller–Uğurlu (p. 7) define isothermal as equal temperatures. For row-stochastic `W` this is equivalent to being doubly stochastic: Keller–Uğurlu p. 7 says "W is isothermal if and only if W is bistochastic". The Lean condition is exactly that. General non-negative weights are excluded and the paper says so (§Scope).
- **`StronglyConnected` means `∀ a b, TransGen (W>0) a b`.** For `n ≥ 2` this is equivalent to the usual definition (every node lies on a cycle). For `n = 1` it needs `W 0 0 > 0`, which row-stochasticity forces. So it is not stronger than the standard notion.
- **Fixation** means `Tendsto (t ↦ sPow (bdK …) t S univ) atTop (𝓝 ρ)`. `sPow` is the Chapman–Kolmogorov t-step kernel. `bd_fixation_mono` and `sPow_top_mono` show the sequence is non-decreasing, so the limit is P(ever reaching `univ`). This identification is Tier C, and the paper labels it that way. The claim is not vacuous: under the same `sPow` limit, the star gives 5/12 ≠ 4/7, and `W = I` gives 0 ≠ 2/3.
- **Every hypothesis is used.** `hn` is used for `totFit_pos` and `∅ ≠ univ`. `hr` is used for positivity and the step lemma. `hW.colSum` provides stationarity of `μ ≡ 1`. `hsc` is used through `cutPositive_of_strongly_connected`, which gives the variance positivity needed by P7's `fixation_limit`. `identity_no_fixation` shows that dropping `hsc` makes the conclusion false.
- **`isothermal_fixation_single`** gives exactly `(1 - r⁻¹)/(1 - r⁻¹ ^ n)` for `r ≠ 1` and every node `v`. My restatement with `1/r`, `1/r^n` (§1) is accepted. For `n = 1` the value is 1, which is correct.
- **`stationary_fixation`** has the same content as Keller–Uğurlu Thm 3.1. They assume strong connectivity and normalised `π`. Lean assumes `CutPositive` and any positive stationary `μ`; `CutPositive` is equivalent to strong connectivity and only the direction actually used is proved, as the paper says. Lean's conclusion is stated per configuration. A KU initial law on `X_j` is a convex combination of configurations, so it follows by averaging. The paper says "same hypotheses and conclusion (up to normalisation)", which is fair.

## 3. Independent numerics (auditor's own exact-rational script, `/tmp/n4audit/audit_numerics.py`)

The script builds the 2ⁿ-state chain from bitmasks and solves `h = Kh` exactly with `Fraction` Gauss–Jordan. It shares no code with `numerics_evograph.py`.

- **Star n=3, r=2, μ≡1:** table `{0}:5/12, {1}:2/3, {2}:2/3, {0,1}:5/6, {0,2}:5/6, {1,2}:14/15`. This is identical to `hStar`. Centre **5/12**, leaf **2/3**, uniform **7/12**. Moran `ρ₃(1) = 4/7`. `P⁺({0}) = P⁻({0}) = 1/2`.
- **Star with stationary μ = (1/2, 1/4, 1/4):** centre 4/7.
- **4-cycle:** `(P⁺, P⁻)` is (1/3, 1/6) at {0,1} and (2/3, 1/3) at {0,2}, so **1/3 versus 2/3**. Every one of the 16 configurations fixes with the Moran value (single mutant 8/15).
- **`W_H` (n = 3):** 4/7. **`W = I` (n = 2):** every row is a delta, so the chain never moves; `ρ₂(1) = 2/3`.
- **Star n=4 (Tier X in the paper):** 392/1265 (centre), 168/253 (leaf), mean 728/1265. These match the paper.

## 4. Human proof, deviation, citations

- **Human proof.** Steps 1–8 are correct. I checked the cut formulas, cut balance, the harmonicity computation `E[Δh] = (r^{-k}/S_n)(P⁺ − rP⁻)`, and the 4-cycle arithmetic (out = in = 1 and 2, F = 6). Two minor wording issues are listed below as NB4 and NB5.
- **Deviation section.** It is honest and correct. The preregistered "count is a birth–death chain" clause is reported as false for the chain with holding steps. Only strong lumpability is refuted, and the paper says so explicitly. The jump-chain remark is correct: `P⁺/(P⁺+P⁻) = 2/3` at both configurations, with holding probabilities 1/2 and 0. The prereg's own negative control (a star that differs on an explicit instance) is delivered.
- **Keller–Uğurlu arXiv:2403.12598v2 (read via alphaXiv):**
  - The abstract does say "Nature Communications" and the reference list does say *Nature* 433:312–316. Verified.
  - "remains unproven till today" (p. 4): verified.
  - "it does not follow that (M_k) is a Markov Chain" (p. 7): verified.
  - Thm 3.1 (p. 9), stationary selection: verified.
  - Thm 3.2 (p. 10), `p⁻/p⁺` constant iff `μ = π`: verified.
  - The source's Thm 2.1 typo (denominator `1 − 1/r`) was correctly not copied.
  - The quoted phrase "the projection does not need to be Markovian" appears in the **abstract only**. P. 4 has different words. See NB1.
- **Pattni et al. arXiv:1505.06021v1 (read via alphaXiv):**
  - Lemma 1 (stated p. 11, proof pp. 11–12) lumps the holding-removed matrix S′ onto the count via Kemeny–Snell Thm 6.3.2. Verified.
  - Prop. 2 (p. 6) states the equivalence with `f_R(W)` being a circulation. Verified.
  - The remark that this is the isothermal theorem of LHN is in **§3.1 on p. 7**, not p. 6. See NB2.
  - Their "isothermal" (`W_I`: T⁺_i = T⁻_j) and "circulation" both reduce to doubly stochastic for right-stochastic W, consistent with the paper.
- **LHN 2005 not read.** This is disclosed in the abstract, the bibliography, §Literature and publication_meta. Moran 1958 and Kemeny–Snell are also disclosed as not read.
- **"Our reading" (Tier C).** The paper does not endorse KU's "unproven" claim, on the strength of Pattni's Lemma 1. That reading is defensible: Pattni's lumpability condition holds because every non-absorbing S has `P⁺ + P⁻ > 0` under strong connectivity. It is correctly labelled Tier C.

## 5. Overclaim check

- No novelty is claimed. "Not aware of an earlier machine-checked formalization" is qualified with "no systematic search", which is acceptable.
- The biology section is labelled Tier C, explicitly limits itself to the model, and says that 7/12 > 4/7 is one data point. No tissue claim is made. Acceptable.
- Every Tier A tag is attached to Lean statements, and the converse direction, death–birth, and general weights are listed as not proved.
- `result.json` (producer file, outside my write scope) glosses `cycle4_not_lumpable` as "(count not Markov)", which is stronger than the paper's careful "strong lumpability only". See NB7.

## 6. publication_meta.json

The title, sha256s, Mathlib commit, ledger counts (11 A / 4 L / 9 C / 2 X = 26, verified in `ledger.json`), the gate result (0 blocking findings, 11 `LEDGER_UNAUDITED_TIER_A` flags, exit 1, verified in `gate_output.txt`), the related identifiers, and the not-read disclosures are all accurate. The **status** sentences ("no independent recompile, no independent AI audit", "Producer result, NOT independently audited") become stale once this audit exists. See B1.

## BLOCKING

**B1. The review-status text is stale.**
- **main.tex, the "AI-assistance disclosure and review status" paragraph (lines 41–46).** Replace "There has been no independent recompile, no independent AI audit and no human audit … a producer≠verifier audit is still owed" with: "An independent adversarial AI auditor (producer≠verifier, not a human) recompiled the EvoGraph and reused P7 modules from source into its own olean directory, re-ran `#print axioms`, ran a `leanchecker` kernel replay (rc 0, with a `debug.skipKernelTC` forgery control rejected), recomputed the star and 4-cycle values with an independent exact-rational script, and checked both cited sources; verdict ACCEPT_WITH_CORRECTIONS (`AUDIT.md`). There has been no human audit."
- **Scope bullet "No independent verification" (line 203).** Change it the same way.
- **The final line "Read every Tier A tag …".** Keep it, since adequacy has still not been audited by a human.
- **publication_meta.json, `description` "Status" paragraph.** Replace "There has been no independent recompile, no independent AI audit and no human audit …" with the same summary.
- **publication_meta.json, `notes`.** Replace "Producer result, NOT independently audited" with "Independently AI-audited (ACCEPT_WITH_CORRECTIONS, papers/p3_isothermal_theorem/AUDIT.md; statement lock results/proposal3/2026-10-26/N4/STATEMENT_LOCK.json); no human audit."
- **Ledger.** Leave the audit fields null, because they mean a human audit. The 11 flags therefore remain and should still be reported.

**B2. The soundness gate is mis-stated.** main.tex line 179 says "exit status is not the gate; the axiom list is." The axiom list is not sufficient on its own. A theorem added under `debug.skipKernelTC` shows only whitelist axioms; I demonstrated this on an N4 statement in §1, and proposal-2 A2 found the same. Change the sentence to: "Exit status is not the gate. The gate is the axiom list **together with** a kernel replay (`leanchecker`) of the modules and a grep for `debug.skipKernelTC`/`addDecl`/`native_decide`/`implemented_by`; the independent audit ran all three (`AUDIT.md`)." Add the same gate definition to the Tier A definition in the "Tiers" paragraph (line 66) and to the abstract's Tier A sentence. In publication_meta's Tier A sentence, add "kernel replay (leanchecker) by an independent auditor".

## NONBLOCKING

- **NB1. Quote attribution.** main.tex line 80 puts "the projection does not need to be Markovian" at "abstract and p. 4". That exact phrase is in the abstract. P. 4 reads "the projection of a Markov Chain is not necessarily Markov, unless under lumpability". Fix: "(arXiv:2403.12598v2, abstract; similarly p. 4)".
- **NB2. Page number.** In line 185 and bib entry `PBRS`, Pattni et al. note that Prop. 2 is the isothermal theorem of LHN in §3.1, **p. 7**, not p. 6. Fix: "Proposition 2 (p. 6; identified with the isothermal theorem of [LHN] in §3.1, p. 7)".
- **NB3. Reproducibility recipe is incomplete** (line 211). Without `--olean`, later files cannot import earlier ones. The olean directory is hard-coded to `formal_cert/.olean_out`. The P7 modules must be compiled first. Fix: "order: `Programmes/PopGen/Absorption.lean`, `Programmes/PopGen/MoranModel.lean`, then Transport, BirthDeath, Examples, each with `--olean Programmes.<Dir>.<Name>`; then the import-all and planted files (no `--olean`); oleans go to `formal_cert/.olean_out` (hard-coded); kernel replay: `leanchecker Programmes.EvoGraph.Examples …` with the same `LEAN_PATH`." Also note that the paper-stage recompile (`paper_stage/recompile_*.json`, `module: null`) reused the producer's oleans for the imports.
- **NB4. Step 5 sketch** (line 127). "∑_{t<N} ∑_u P^t(s,u) v(u) ≤ 1" assumes `0 ≤ h ≤ 1`. That holds for `ρ_n(|S|)` but is not a hypothesis of `two_absorbing_limit`. Fix: "≤ max_u h(u)² (= 1 here, since 0 ≤ ρ_n ≤ 1)".
- **NB5. Step 6 wording** (line 132). "the last vertex of the path in S before the path first leaves it" is awkward. Fix: "let ℓ be the least index with x_{ℓ+1} ∉ S (it exists since x_0 ∈ S, x_m ∉ S); then x_ℓ ∈ S".
- **NB6. Kernel replay not mentioned.** Add a Reproducibility bullet naming the `leanchecker` replay and its forged-declaration control (rejected), and state that the replay was not `--fresh`: Mathlib and core were loaded from the pinned build, not replayed.
- **NB7. `result.json` wording** (producer artifact; I did not edit it). The `cycle4_not_lumpable` gloss "(count not Markov)" should read "(count not strongly lumpable; holding chain)", in line with paper §1 and EG-C-0001.
- **NB8. Optional.** State explicitly in §Model that `Isothermal` ⇔ equal temperatures holds *because* W is row-stochastic (Keller–Uğurlu p. 7). The paper only cites EG-L-0002 for "equivalently doubly stochastic".

## Evidence summary

- Auditor compile JSONs: `/tmp/n4audit/out/{Absorption,MoranModel,Transport,BirthDeath,Examples,ImportAll,Planted,AuditN4,AuditPlanted}.json`
- Forged control: `/tmp/n4audit/src/AuditForged.lean`, `leanchecker` rc 1
- Numerics: `/tmp/n4audit/out/numerics.json`
- Statement lock: `results/proposal3/2026-10-26/N4/STATEMENT_LOCK.json` (sha256 475c3c97fa18ef0b379ca1fadcdcd3c3063be05e06e3186722e52fc2226b9fa9)

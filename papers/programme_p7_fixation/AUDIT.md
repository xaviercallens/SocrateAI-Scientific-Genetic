# Independent audit: P7 "Fixation probabilities as absorption limits of matrix powers, kernel-checked"

- **Auditor:** an AI agent (Claude, Opus-family). I am not the producer and not the paper-stage agent. This is the "producer ≠ verifier" stage of the Elenchus protocol. It is **not** a human audit.
- **Date:** 2026-10-25.
- **Object audited:**
  - `papers/prog_population_genetics_fixation/`: `main.tex`, `main.pdf` (sha256 `ce4e89b3…`, 14 pp.), `ledger/`, `publication_meta.json`, `data/`, `gen/`;
  - `results/programmes/2026-10-25/P7/`;
  - `formal_cert/Programmes/PopGen/*.lean` (11 files);
  - rusty-SUNDIALS worktree `/home/callensxavier_gmail_com/rusty-SUNDIALS-wt-p7`, branch `physics/p7-fixation-2026-10-25`, commit `2b716c58fc7175798ad54dd9363938a3df37d9b0`, compared against `origin/main` = `7630f7e6436f6791a8abc9fea05b76868305bc35`.
- **Write scope:** I edited nothing outside this file. My scratch work is under `/home/callensxavier_gmail_com/programmes-2026-10-25/P7/audit/`.

## Verdict: ACCEPT_WITH_CORRECTIONS

- Every mathematical claim I checked holds at the tier it is given:
  - all 14 Tier A rows recompile from source under my own build, with the allowlisted axioms only;
  - the planted controls are caught;
  - every number traces to a result file;
  - the CVODE run reproduces byte for byte;
  - all four arXiv citations say what the paper says they say.
- The corrections are wording and status updates. None changes a result.
- The most important correction is that the "no independent verifier / no AI audit" language is now out of date (see C1).

## 1. Lean: independent recompile, own axiom audit, own planted controls

**Method.**
- My driver is `audit/audit_lean.py` (sha256 `3d609e69…`). It rebuilt every module **from source** into my own olean directory, `audit/olean/`.
- The producer's `formal_cert/.olean_out` was **not** on `LEAN_PATH`. The script asserts this before running.
- Rebuilt modules: H11 `Autoresearch.PopGenWrightFisher`, H12 `Autoresearch.PopGenFisher`, and the five P7 modules (`Absorption`, `MoranModel`, `WrightFisherFixation`, `MoranFixation`, `Price`).
- It then compiled the two import-all files, `FixationDataCert`, the three planted-controls files, and three files of my own.
- Toolchain: pinned Lean 4.34.1 and Mathlib from `/mnt/disks/disk-socrateai-local-1/callensxavier_home_data/openai_math_pinned/lean`, run on cores 4–7.

**Results** (from `audit/audit_lean.log` and `audit/out/*.json`):

| file | sha256 (first 16) | rc | `#print axioms` lines | outside allowlist |
|---|---|---|---|---|
| H11 PopGenWrightFisher | 9cac6b9d5d0c76be | 0 | 13 | 0 |
| H12 PopGenFisher | 1625d54caeee4472 | 0 | 13 | 0 |
| Absorption | 86fb5168a6b1efe0 | 0 | 11 | 0 |
| MoranModel | b0822dc0ee83505f | 0 | 16 | 0 |
| WrightFisherFixation | 149a9328171a6347 | 0 | 5 | 0 |
| MoranFixation | c49806ac8b14295d | 0 | 8 | 0 |
| Price | 1309c57a7818da57 | 0 | 21 | 0 |
| PopGenProgImportAll | e26cdffa01ab893b | 0 | 6 | 0 |
| PriceImportAll | 3e82333747ae4bdd | 0 | 4 | 0 |
| FixationDataCert | 79ced614149be765 | 0 | 7 | 0 |
| PopGenProgControls (planted) | 685cdd545b9cff73 | 0 | 2 | 2 (`sorryAx`; `ProgPopGenControls.ctl_cheat`) |
| PriceControls (planted) | 93bc16d08561f582 | 0 | 2 | 2 (`sorryAx`; `ProgPopGenPriceControls.ctl_cheat`) |
| FixationDataCertControls (planted) | 098c2ed751a65886 | 0 | 2 | 2 (`sorryAx`; `ctl_cheat`) |
| **audit** AuditAxioms (25a677b1…) | | 0 | 38 | 0 |
| **audit** AuditPlanted (158e49ec…) | | 0 | 3 | 3 (caught, see below) |
| **audit** AuditMutation (2c482463…) | | **1** | 0 | kernel rejects (`unsolved goals ⊢ False`) |

**What the table shows.**
- The sha256 of every audited file equals the value in `publication_meta.json` and in the paper-stage `summary.json`. The H11 hash `9cac6b9d…` also matches.
- `audit/check_evidence.py` compares every `{path, sha256}` pair recorded inside the 35 ledger evidence blobs with the file on disk: 103 references match, 0 mismatch, 0 missing.
- The pinned `lake-manifest.json` gives Mathlib rev `d13f23b723b8a846827a245b89c10fc7d3f11612`, as the paper and meta state.
- The axiom line counts equal the paper's Table 1 (11, 5, 16, 8, 21, 7, 6, 4, 2, 2, 2).
- Every claimed declaration prints exactly `[propext, Classical.choice, Quot.sound]`.
- The only warnings are of two kinds:
  - deprecation warnings (`if_pos`/`if_neg`/`if_true`) in `Absorption` and `MoranModel`, as the paper says;
  - the expected `declaration uses 'sorry'` in the planted files.

**`AuditAxioms`.** It imports the claimed modules and prints axioms for 36 key declarations. These include:
- `fixation_limit`, `harmonic_unique_of_var`, `wright_fisher_fixation_limit`;
- all the `moran_fixation_*` theorems and `moran_closed_tendsto_neutral`;
- `price_equation`, the `fisher_haploid*` theorems, both `Alix.*` restatements, and every `ctl_pos`/`ctl_neg` theorem.

It also `#check`s the main statements; the statements printed match the paper's listings. Two non-vacuity theorems of my own also check:
- `audit_moran_N2`: `kPow (moranP 2 2) n 1 2 → 2/3`, derived from the claimed theorem;
- `audit_moran_row`: the Moran row at N=2, r=2, i=1 is (1/6, 1/2, 1/3).

**`AuditPlanted`.** My own planted controls are all caught by `#print axioms`, while the file still exits rc 0:
- a `sorry` on a false neutral-limit claim → `sorryAx`;
- a planted axiom on a false Price claim → `P7AuditPlanted.audit_cheat`;
- a `sorry` hidden inside a lemma used by an honest-looking theorem → `sorryAx` is inherited.

This confirms that exit status is not the gate and that the axiom audit is.

**`AuditMutation`.** It changes the positive control `moranFix 3 2 1 = 4/7` to `= 1/3`. The kernel rejects it (rc 1). The checked values are therefore sensitive to the arithmetic.

**Forbidden tokens.** A whole-word grep for `sorry|admit|native_decide|axiom` over the 8 claimed or check-only files finds no use, only two doc-comment mentions: `FixationDataCert.lean:15` and `MoranModel.lean:237`. The paper says "a doc-comment line"; there are two (C6).

**Statement adequacy** (AI judgement only, not a human audit):
- `AbsorbingChain` encodes a stochastic kernel with point-mass rows at 0 and N.
- `kPow` is the left Chapman–Kolmogorov recursion.
- The variance hypothesis is the one-step conditional variance.
- `moranUp`/`moranDown` match Keller–Uğurlu eq. (3) exactly. I verified this symbolically; see §3.
- `bdP` uses ℕ-valued `Fin` coercions correctly. `p N = 0` and `q 0 = 0` make the out-of-range jumps vanish, and `bdP_sum` handles the edge cases.
- `wfP` (H11) is the binomial kernel `C(N,j)(i/N)^j(1-i/N)^(N-j)`.
- `zbarNext` is the offspring mean `Σ n w z' / Σ n w`. As the paper says, this definition is the modelling assumption.
- I found no vacuity: `AbsorbingChain (moranP 3 2)`, `AbsorbingChain (wfP 3)` and my N=2 instance hold.

## 2. Numbers trace to result files; spot re-runs

**CVODE re-run.**
- `git status` in the rusty worktree is clean.
- I rebuilt `qf-popgen` at `2b716c5`: `cargo build --release -j 2 --offline`, finished.
- I re-ran `fixation_ode` into `audit/cvode_rerun/`. Exit 0; summary `{"rows":115,"accepted":106,"refused":9,"controls_or_runs_not_as_expected":0,…}`.
- `results.jsonl` sha256 `9860803d…8c4a` is **byte-identical** to `papers/…/data/results.jsonl` and to `crates/qf-popgen/results/results.jsonl`.

**Macros and tables.** `audit/audit_numbers.py` recomputes every `gen_numbers.tex` macro from `data/results.jsonl`. All match: nRows 115, nAcc 106, nRef 9, maxErr 3.1e-10, maxInterior 6.3e-10, maxDrift 9.8e-14, nMoranRuns 78, nWFRuns 28. Every `table_cvode.tex` row also matches (starts and max err). The exception is the `wf_neutral` row, which prints N = 10 although it pools N = 10 and N = 20 (C4).

**Exact check.** I computed exact `Fraction` values of S_i/S_N, or i/N, for all 106 accepted rows. The largest |p_fix − exact| is 3.13e-10, consistent with maxErr.

**Data-certificate literals.**
- All six literals in `FixationDataCert.lean` are the `p_fix` doubles of the matching (N, r, i) rows.
- My errors against the exact values are 6.57e-15, 9.55e-14, 7.79e-14, 6.49e-13, 8.46e-13 and 2.30e-16. Ledger PG-X-0004 gives 6.6e-15, 9.5e-14, 7.8e-14, 6.5e-13, 8.5e-13 and **2.2e-16**.
- The last value is the float64 difference. The decimal literal as written in Lean differs from 3/10 by 2.30e-16 (C3). The Lean bound (<1e-12) is unaffected.

**Independent matrix powers.**
- Moran N=3, r=2, i=1, with exact rational P^n computed by my own kernel constructor: P^n(1,3) is 0, 0.0889, 0.5098, 0.571427, 0.5714285714260 and 0.5714285714285714 at n = 1, 2, 10, 50, 100, 400. This is monotone, and the gap to 4/7 at n=400 is 6.1e-47.
- Moran N=10, r=3/2, i=1, with mpmath at 50 digits: P^2048(1,10) differs from the exact closed form by 1.6e-25.
- Negative controls on the same run: the gap to the neutral value 1/10 is 0.23922, and to the flipped formula 0.33039.

**Other traced numbers.**
- `numerics_fixation.json`:
  - Wright–Fisher errors ≤1.94e-61 (paper: ≤2e-61);
  - Moran r=3/2 errors ≤3.64e-50 (paper: ≤3.7e-50);
  - r=0.8 errors ≤2.17e-41 (ledger: ≤2.2e-41);
  - r=1 errors 6.6062e-37, which the ledger gives as "≤6.6e-37"; rounding goes the wrong way (C7);
  - flipped and neutral gaps 0.33–0.67 and 0.24–0.42.
- `price_numerics.json`, recomputed exactly by me: −13/36 = −13/36; −1/2 versus 1/6 with Price −3/4 = −3/4; 1/6 versus 1/4.
- `independent_check.json`: expm differences ≤1.02e-13, matching the paper. I did not re-run `gen_all.py`, because it writes into the paper directory.
- The negative-control values in §4 (0.339/0.716/0.958 against 0.0088/0.042/0.284 and 0.1/0.3/0.7) match `results.jsonl`.

**Rust build hygiene, re-run by me.**
- `cargo test -p qf-popgen -- --test-threads=2`: 6 passed, 0 failed.
- `cargo clippy --workspace --features experimental -- -D warnings`: finished with no warnings.
- `cargo fmt --all --check`: rc 0.
- The diff against `origin/main` is a one-line workspace-member addition in `Cargo.toml` plus the new crate `crates/qf-popgen`. No CVODE core file changed, as the paper states.
- `git branch -a --contains 2b716c5` lists only the local branch, so it is not pushed.

## 3. Human proofs, step by step

`audit/audit_sympy.txt` records sympy residuals of 0 for:
- p − r·q;
- the stay-probability identity `1 − p − q = (r i² + (N−i)²)/((r i + N − i) N)`;
- Keller–Uğurlu p± versus `moranUp`/`moranDown` (PG-C-0004);
- the first-step equation for h = S_i/S_N;
- the Price five-sum scalar identity and the raw-moment form of Cov.

It also confirms the limits (1−x^{-i})/(1−x^{-N}) → i/N as x → 1 for three (i, N).

| step | verdict |
|---|---|
| 1 martingale `P^n h = h`, row sums, nonnegativity | correct (point-mass rows give `Ph = h` at the ends; induction) |
| 2 second-moment identity | correct (expand the square; uses row sum 1 and harmonicity) |
| 3 summability of the interior mass | correct: `P(h²) ≥ h² + ε·1_int` pointwise, apply `P^m` (nonnegative), telescope, bound `P^n(h²) ≤ Σ h(j)²` by row sums; ε is a min over a finite set (the empty-interior case is handled in Lean) |
| 4 squeeze | correct: `|h(j) − 1_{j=N}| ≤ H·1_int(j)`; mass split gives the limit at 0 |
| 5 uniqueness by uniqueness of limits; monotonicity | correct; I_n → 0 depends only on P (given one h with v > 0), so it applies to g |
| 6 Wright–Fisher variance | correct: `P(k,0)(k/N)² = (1−k/N)^N (k/N)² > 0` |
| 7 Moran stochasticity | correct (sympy) |
| 8 first-step equation ⇔ p_i = r q_i; variance ≥ p_i (r^{-i}/S_N)² > 0 | correct (sympy) |
| 9 closed form; continuity at r = 1 | correct |
| 10 Price | correct (sympy) |
| 11 haploid increment (z = z' = w) | correct; the scope remark is accurate |

The two Tier C bridges are correctly labelled as not formalised:
- PG-C-0001: `kPow` versus `P_i(T_N ≤ n)`.
- PG-C-0002: the Poisson mixture for the uniformised ODE. Monotone bounded convergence of P^n(i,N) gives convergence of the Poisson mixture.

The human proofs contain no gap.

## 4. Tiers and overclaim

- Tier assignments are honest:
  - Tier A attaches to Lean statements about matrix powers;
  - the path-measure reading is Tier C;
  - CVODE results are Tier X and are explicitly said "not to support any claim";
  - the data certificates are described as arithmetic between a literal and an exact rational, not as a certificate of CVODE.
- "No novelty claimed" appears in the abstract, PG-C-0006 and the Lean headers. There is no "breakthrough", "first" or "proved" language beyond Tier A.
- The prereg item `fisher_general` (prior 0.35) is reported as **not achieved** in the abstract, Remark (Scope of Step 11), the Limitations and PG-C-0005.
- The prereg file is unchanged since commit `a549c2f`. Its P7 entries are H11_gap, moran_selection, price and fisher_general, and the paper reports all four honestly.
- The "H11 gap" is closed at the matrix-power level only. The paper says so: "closes a gap" sits next to "identification … is Tier C".
- Nothing depends on openai/math (PG-C-0008). The modules import only Mathlib, H11 and H12.

## 5. Citations

**alphaXiv page text.** I read the page text and matched every Tier L excerpt by hand:
- Keller–Uğurlu arXiv:2403.12598v2:
  - p. 2 has eq. (3) and "The individual to be replaced is chosen uniformly at random", plus the "r = 0" typo;
  - p. 3 has the optional-stopping argument on r^{-M_k};
  - p. 4 has eq. (4);
  - p. 7 has Theorem 2.1 with the denominator typo `1−1/r`;
  - p. 18 has [Mor58].
- LaGatta arXiv:2202.10289v4: eqs. (1.6) and (1.8) on p. 9, (2.24) on p. 18, (2.27) on p. 19, the Price (1972) quote on p. 4, and [Pri70] and [Pri72a] on p. 69.
- Baez arXiv:2107.05610v3: p. 1 ("only true under very restrictive conditions"), p. 2 (the Ewens quote), and p. 4 (constant fitness gives dF̄/dt = Var(f)).
- Marion arXiv:2506.18616v5: p. 1 ("fully integrated to Mathlib"), p. 5 (the Markov chain from a constant kernel), p. 9 (`partialTraj`) and p. 12 (`def traj`).

**`citation_check.py`.** Its scope is theorem-kind only, so I ran it on KU Theorem 3.1 (p. 9).
- **Primary source.** I fetched `https://arxiv.org/pdf/2403.12598v2` (18 pp., sha256 `1ef6cce2…`) into a fresh directory and ran `pdftotext` on it, giving `audit/cite/primary/ku.txt` (sha256 `c2753bbb…`).
  - positive claim (the theorem's hypothesis sentence): **ACCEPTED**, exit 0;
  - negative control, the same quote with "Definition 2.1" changed to "Definition 2.2": **REFUSED** with `CITE_QUOTE_ABSENT`, exit 1.
- **Earlier check on weaker evidence.** Before that, I ran the tool on my transcription of the alphaXiv extraction (`audit/cite/ku_p9.txt`). The positive claim was ACCEPTED, and a negative control with the exponent flipped was REFUSED.
- **Scope.** The displayed formula itself is scattered by PDF extraction, so the tool checks the hypothesis sentence. The formula was compared by eye on the alphaXiv and primary page text. The other excerpts are prose or displayed equations, outside the tool's scope (as the ledger note states), and I matched them by hand as above.

**Limits of this check.** Except for KU, the sources were read through alphaXiv's PDF extraction. Moran (1958) and Price (1970, 1972) were not read by the producer or by me, and the paper says so.

**Minor points.**
- KU p. 3 also has a sign slip, `r^i = E[r^{M_0}] … r^n`, where `r^{-i}` and `r^{-n}` are meant. The paper's "source typos on pp. 2 and 7" list does not include it. This is optional to add (C8).
- PG-L-0003 and §1 say `Kernel.traj` "is present in the pinned build [L, PG-L-0003]". Presence in the pinned build is a repository fact: `ProbabilityTheory.Kernel.traj` is at `Mathlib/Probability/Kernel/IonescuTulcea/Traj.lean:518` in the pinned Mathlib, which I checked by grep. Marion supports only "in Mathlib". The wording should separate the two (C5).

**Genetic crosswalk.**
- `meanFit`, `selectP`, `aurelia_mean_fitness_increment` and `kal_increment_is_variance` are textually identical to `/home/callensxavier_gmail_com/SocrateAI-Scientific-Genetic-pr/formal/ANSE/AlixPopulationGenetics.lean` (lines 37, 40, 43–45, 72–74).
- That repository has no Price result. Grep finds "price" only as English in `InteractingStagesAndTransport.lean`.

## 6. Data provenance

- `data/results.jsonl` and `summary.json` come from `fixation_ode` at commit `2b716c5`, and I reproduced them byte for byte.
- `data/numerics_fixation.json` and `data/price_numerics.json` are byte-identical to the copies in `results/programmes/2026-10-25/P7/`.
- The ledger evidence blobs carry the sha256 of their files. The Elenchus gate, which I re-ran, gives output identical to `ledger/gate_output.txt`: 0 blocking findings, 14 `LEDGER_UNAUDITED_TIER_A` flags, exit 1.
- H11 and H12 have separate-AI-verifier verdicts of ACCEPT_WITH_CORRECTIONS, recorded in `results/autoresearch/2026-10-25/SUMMARY.md` and `H11/verification.md`, as the paper says.

## 7. publication_meta.json

- The title, description, keywords and notes agree with the paper, and the hashes are correct.
- The key set equals that of all five `ar_*` metas and contains every field `scripts/publish_ar_paper.py` reads (title, description, keywords, version, related_identifiers, notes).
- The five related DOIs, 10.5281/zenodo.23646103/05/07/11/17, are the companion papers' **reserved** DOIs (title pages of `ar_*/main.tex`). If those depositions are not published first, these `references` relations will not resolve. The parent should publish them first or accept this.
- The description and notes repeat the stale status statement in C1.

## Corrections (all nonblocking)

**C1. Status language is now out of date.** This audit is a separate-AI recompile from source.
- Places to update:
  - abstract disclosure: "no independent verifier … has yet recompiled them, and no AI audit of this preprint has been done";
  - Tiers paragraph: "producer-checked and re-compiled by the same AI pipeline only";
  - Limitations bullet: "Producer ≠ verifier is pending";
  - PG-X-0008 context;
  - `publication_meta.json` description ("no independent verifier has recompiled them") and notes ("no independent verifier yet").
- Suggested replacement: "recompiled from source in a separate build directory by an independent AI verifier, with its own axiom audit and planted controls (AUDIT.md, verdict ACCEPT_WITH_CORRECTIONS); no human audit."
- The 14 `LEDGER_UNAUDITED_TIER_A` flags must stay, and the ledger `audit` fields must stay null. An AI audit is not the human statement audit that flag asks for, so the gate should still exit 1.
- Suggestion: add a new ledger row, e.g. PG-X-0009, "independent AI recompile from source", with `audit/out/*.json` as evidence. This makes the recompile traceable without touching the `audit` fields.

**C2. Listing 11 caption.** It says "two of the seven statements", but three are shown (`cert_N3_r2_i1`, `cert_N10_r1p5_i1`, `ctl_neg_wrong_literal`). Separately, the "verbatim" listings substitute glyphs (`𝓝` → `\mathcal{N}`, `ℝ` → `\mathbb{R}`, `¬`). Say "verbatim up to glyph rendering".

**C3. Literal provenance (§4, PG-X-0004).**
- The Lean literals are the shortest round-trip decimal forms of the `p_fix` doubles. They are not the digit strings `fixation_ode` writes; for example, `0.6666666666666601` is written as `6.66666666666660079e-1` in the JSONL.
- The (10,1,3) error is 2.30e-16 for the decimal literal and 2.2e-16 as a float64 difference. Say which one is meant.

**C4. `table_cvode.tex` `wf_neutral` row.** The N column shows 10 but the row pools N = 10 and N = 20 (28 starts). Print "10, 20" or split the row.

**C5. PG-L-0003 and §1.** Split "present in the pinned build" (a repository fact, checked by grep at the Mathlib path above) from the Tier L Marion citation, which supports only "in Mathlib". Optionally write the full name `ProbabilityTheory.Kernel.traj`.

**C6. §5 forbidden-token sentence.** It says "a doc-comment line mentions the word". There are two lines: `FixationDataCert.lean:15` and `MoranModel.lean:237`.

**C7. PG-X-0005.** It says "r = 1 errors ≤ 6.6e-37", but the file has 6.6062e-37. Write ≤ 6.7e-37 or 6.6e-37 (rounded).

**C8 (optional). Cosmetic items.**
- KU p. 3 sign slip: optionally add it to the typo list.
- `FixationDataCert.lean` doc comment, lines 15–16, points its planted controls to `PopGenProgControls.lean`/`PriceControls.lean`; they are in `FixationDataCertControls.lean`. Fixing this changes the file hash, so re-hash and recompile if you edit it.
- PG-C-0008's wording "import only Mathlib (plus H11, H12 and Mathlib)" is garbled.
- The identity-kernel negative-control row in `results.jsonl` carries a meaningless `r: 1.5`, because r is hard-coded in the `negative` closure.

**C9 (Rust, nonblocking, for a follow-up).** `certify_absorption` checks the interior mass as a signed sum and the total mass drift. It never checks that each component satisfies π_k ≥ −tol. A run whose interior entries cancel could pass the interior test. This does not affect any reported run, since mass drift ≤ 9.8e-14 and the expm cross-check agrees to 1e-13, but a fail-closed certificate should add a componentwise check.

## Limits of this audit

- The verifier is an AI from the same model family as the producer. Statement adequacy, definition fidelity and the review of the hand proofs are AI judgements; no human has audited them.
- I did not re-run `gen_all.py`, `gen_listings.py` or `numerics_fixation.py`, which write into the paper or results directories. I recomputed their key outputs independently instead (§2).
- Citations were checked against alphaXiv's PDF text extraction.

## Auditor's scratch artefacts

Everything is under `/home/callensxavier_gmail_com/programmes-2026-10-25/P7/audit/`:
- `audit_lean.py` (3d609e69…), `audit_lean.log`, `out/*.json`, `olean/`;
- `src/AuditAxioms.lean` (25a677b1…), `src/AuditPlanted.lean` (158e49ec…), `src/AuditMutation.lean` (2c482463…);
- `audit_numbers.py` (3dd11fcf…) and `audit_numbers.json` (dc6e24b8…);
- `audit_sympy.py` (e8b174d9…) and `audit_sympy.txt`;
- `cvode_rerun/` and `gate_rerun.txt`;
- `check_evidence.py`;
- `cite/primary/ku.pdf` and `cite/primary/ku.txt`, `cite/claim_primary_pos.json`, `cite/claim_primary_neg.json`;
- `cite/ku_p9.txt`, `cite/claim_pos.json`, `cite/claim_neg_flipped.json`.

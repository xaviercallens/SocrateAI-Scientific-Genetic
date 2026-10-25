# Autoresearch program, 2026-10-25 (Karpathy-style loop over Lean 4 hypotheses)

This adapts Karpathy's *autoresearch* idea: an agent repeatedly changes one artifact, runs a fixed-budget experiment,
reads one metric, and keeps or discards the change. Here the artifact is a Lean 4 file and the judge is the Lean kernel.

## Loop (per hypothesis)
1. **Preregistration** (`preregistration.json`): statement, prior, kill criterion, toolchain and budget. It is
   committed before any attempt and never edited afterwards. Outcomes go in `results.tsv` and in the per-hypothesis
   folder.
2. **Falsification gate first.** Where numbers can refute the statement (complete-monotonicity sign checks, a
   finite-range search), run them before any Lean work. A refutation is a result, and it is logged.
3. **Attempt** within a fixed budget: **45 min wall clock on 2 cores** per attempt, at most 2 attempts per
   hypothesis in this round. A stuck attempt is discarded, not extended.
4. **Metric (binary), KEEP only if all hold:**
   - the file compiles with the toolchain named in the preregistration (rc 0);
   - `#print axioms` of every claimed theorem is within {propext, Classical.choice, Quot.sound};
   - the planted controls are still caught (`ctl_sorry` shows sorryAx; `ctl_axiom` shows its axiom);
   - a grep for sorry, admit, native_decide or a user `axiom` finds nothing.
   Otherwise DISCARD and log the error class.
5. **Producer ≠ verifier.** A KEEP becomes Tier A only after an independent agent recompiles it and audits statement
   adequacy.

## Rules
- Corollaries of `universal_energy_minimum` use the pinned openai/math build (Lean 4.34.1,
  `scripts/certified_numerics/lean_pinned.py`). QuantumFluids and LeanMaster work uses their own v4.34.0-rc2
  environments. Toolchains are never mixed.
- At most 3 concurrent attempts, cores 4–7 (cores 0–3 run the jellium study). One prover model on the T4 at a time.
- Every Tier A claim is conditional on upstream's theorem (AI-generated, not peer reviewed) and on definition
  fidelity.
- "Breakthrough" wording is banned. The genetics items are first formalizations of known results, not discoveries.
  Novelty claims need a literature search first.

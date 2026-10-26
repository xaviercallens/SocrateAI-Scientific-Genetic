# N4 paper revision (2026-10-26), after AUDIT.md (ACCEPT_WITH_CORRECTIONS)

Lean files untouched (locked; `git diff` on `formal_cert/Programmes/EvoGraph` is empty).

## Blocking
- B1 (stale review status): replaced in main.tex (AI-assistance disclosure, Scope bullet now "Independent AI audit, no human audit"), publication_meta.json (`description` Status paragraph, `notes`), and the ledger source (`ledger/build_ledger.py`: STATUS string and row EG-C-0009; ledger.json regenerated, 26 rows). Summary cites AUDIT.md and STATEMENT_LOCK.json. Ledger audit fields remain null (human audit), so the 11 `LEDGER_UNAUDITED_TIER_A` flags remain.
- B2 (gate definition): the gate is now axiom list AND `leanchecker` kernel replay AND a grep for `skipKernelTC`/`addDecl`/`native_decide`/`implemented_by`, with the reason (`#print axioms` can be fooled by `debug.skipKernelTC`). Applied in the abstract Tier A sentence, the Tiers paragraph, the "Verbatim axiom transcripts" paragraph, and publication_meta's Tier A sentence.

## Non-blocking
- NB1 quote attribution: "abstract; similarly p. 4".
- NB2: Prop. 2 p. 6, identified with LHN isothermal theorem in section 3.1 p. 7 (literature item and bib entry).
- NB3: reproducibility recipe corrected (`--olean`, P7 modules first, hard-coded `formal_cert/.olean_out`, paper-stage recompile reused producer oleans).
- NB4: Step 5 bound is `max_u h(u)^2` (= 1 since 0 <= rho_n <= 1).
- NB5: Step 6 uses the least index ell with x_{ell+1} not in S.
- NB6: new Reproducibility bullet "Kernel replay" with the forged-control rejection (rc 1) and the not-`--fresh` limitation.
- NB7: `results/proposal3/2026-10-26/N4/result.json` gloss changed to "(count not strongly lumpable; holding chain)".
- NB8: Model section states that equal temperatures is equivalent to column sums 1 because W is row-stochastic (Keller-Ugurlu p. 7).

## Results
- PDF rebuilt with `build_pdf.py`: 12 pages, 0 undefined references, 0 overfull boxes.
- Ledger gate (`ledger/run_gate.py`): 0 blocking findings, 11 `LEDGER_UNAUDITED_TIER_A` flags, exit code 1 (expected while audit fields are null).
- No human audit exists; LHN 2005, Moran 1958 and Kemeny-Snell 1960 remain unread (disclosed).

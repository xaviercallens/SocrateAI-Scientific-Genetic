# Autoresearch round 2026-10-25: 13 preregistered hypotheses, results

Preregistration `924c79d` came before any attempt. The program is `program.md`, a Karpathy-style keep/discard loop with
the Lean kernel as judge, and every attempt is logged in `results.tsv`. Priors were for "a quick time-boxed attempt
yields the stated result", not for "breakthrough". No result below is a major breakthrough. The most useful ones are
new kernel-checked corollaries for physically important interactions, and kernel checks of the owner's own published or
roadmapped claims.

| id | result | prior | outcome | verification |
|---|---|---|---|---|
| H1 | Bilayer dipolar-exciton potential 2(1/r − 1/√(r²+d²)) admissible ⇒ triangular lattice optimal (classical) | 0.70 | KEEP, 1 attempt, 4 min | ACCEPT_WITH_CORRECTIONS |
| H1b | Its lattice energy is finite, so H1 is informative (unlike Coulomb) | 0.55 | KEEP | ACCEPT_WITH_CORRECTIONS |
| H2 | Planar Yukawa e^{−κr}/r admissible ⇒ triangular optimal; finite energy; reusable product closure `admissible_mul` | 0.35 | KEEP, 11 min | ACCEPT |
| H3 | Rytova–Keldysh (TMD) potential completely monotone in r² (explicit Gaussian-mixture weight) | 0.55 | KEEP, Tier C (proof + numerics 2e-21) | producer only |
| H4 | Montgomery lattice theorem | 0.08 | not attempted (scoped) | — |
| H5 | Riesz s = 2 endpoint | 0.02 | not attempted (open) | — |
| H6 | Dyadic triangular Hilbert C* ≥ 3: product rule and tensor powers, kernel-checked | 0.30 | KEEP, 7 min | ACCEPT |
| H7 | σ = 3 cutoff-uniform enstrophy bound 2DΩ ≤ H₃ + (2E)^{3/2} (QuantumFluids) | 0.30 | KEEP, 17 min | producer only |
| H8 | Bogoliubov S(k) ≤ √(k/2k*), equality iff k = k* (QuantumFluids) | 0.45 | KEEP on an out-of-protocol 3rd run (2 infra discards) | ACCEPT_WITH_CORRECTIONS |
| H9 | NavierStokesAndEuler breakdown theorem vs Fefferman's official Clay statement | 0.35 | FAITHFUL_WITH_CAVEATS (statement audit only) | — |
| H10 | Darmon–Merel on the FLT corpus | 0.01 | untested (corpus unbuilt) | — |
| H11 | Neutral Wright–Fisher: i/N is the unique harmonic function (fixation characterisation) | 0.35 | KEEP (formulation (a) only) | ACCEPT_WITH_CORRECTIONS |
| H12 | One-locus two-allele mean-fitness monotonicity + Hardy–Weinberg | 0.60 | KEEP | ACCEPT_WITH_CORRECTIONS |
| H13 | Tadpole uniqueness on LeanMaster constants (integer arithmetic) | 0.35 | KEEP, 2nd attempt (1st infra) | ACCEPT_WITH_CORRECTIONS |

**Calibration.** Of the 10 Lean/proof items attempted, 9 were kept within protocol (H1, H1b, H2, H3, H6, H7, H11,
H12, H13). H8 was kept only on an out-of-protocol infrastructure retry. 8 of them are independently verified: H1,
H1b, H2, H6, H8, H11, H12 and H13; H3 and H7 are producer-only. (Commit 6633bb8's message says "7"; 8 is correct.) The
mean prior was 0.45, so the priors were **underconfident**. The attempted items were corollaries or specializations
of existing machinery (upstream's universal optimality, or the owner's own lemmas), which made them easier than
priced. The round says nothing about breakthrough odds, because nothing hard was attempted: H4, H5 and H10 were not
attempted at this budget.

**What every Tier A item rests on.**
- The H1/H1b/H2 corollaries are conditional on `OAI.AtomicTriangular.universal_energy_minimum` (openai/math,
  AI-generated, not peer reviewed) and on how faithfully its definitions match the physics. They are classical, at
  zero temperature and density one (other densities via D2).
- Novelty: the complete-monotonicity steps are elementary or classical. What no paper found states is the explicit
  kernel-checked corollaries for the bilayer exciton and Keldysh potentials. H11/H12 are formalizations of textbook
  results, with novelty unchecked.
- H3 (Keldysh) and Coulomb-like potentials have an infinite raw lattice energy. The theorem then reads ∞ ≤ ∞, and
  the physically relevant jellium version is not covered.

**Infrastructure.** Disk 2 was I/O-starved by another session's `gcloud storage rsync` (about 18 workers) on top of
a cold page cache after a reboot. That cost H8 two attempts and H13 one. The fix is to warm the `.olean` closure with
sequential reads (`tmp/warm_siblings.py`) and run `lean` directly with LEAN_PATH.

**H9 caveat, which is the highest-stakes line here.** "FAITHFUL_WITH_CAVEATS" is a statement-only comparison:
- by one auditor, with no second fidelity audit;
- with Mathlib definitions read in a v4.32.2 copy, not the pinned rc2;
- on a corpus that was **not** rebuilt here and whose proofs and axioms were not checked by us.

Nothing here checks the corpus's proof, and the Clay Institute's award has its own process. It must not be read as
"Clay Navier–Stokes settled".

**Open next steps.**
1. H9: a second independent fidelity audit, rebuilding the NavierStokesAndEuler corpus (rc2 Mathlib in qf-lake,
   already warm), and `#print axioms` on `navier_stokes_breakdown_R3` and `_periodic`.
2. Independent verification of H3 and H7.
3. H3L: a Lean lemma "positive Gaussian mixture ⇒ admissible" (prior 0.30).
4. H11: close the first-step gap with P^n(i,N) → i/N.
5. H7: H₃ conservation along the σ = 3 flow, in Lean.
6. A literature pass before any public novelty wording (H1, H3).

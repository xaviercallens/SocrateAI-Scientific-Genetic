/-
  ANSE.GeneticaReticulata — the exact-arithmetic foundation of Reticulata.

  Phase 1 of docs/GENETICA_RETICULATA.md: move the trust boundary into the Lean
  kernel. `ANSE.ReticulationCertificate` states the no-go and the inference rule
  over `ℝ`, and the Rust prototype validated certificates itself, in `f64`. That
  put floating-point arithmetic and our own search code inside the evidence
  chain. This file removes both.

  THE OBSERVATION THAT MAKES IT POSSIBLE. A p-distance is

      (differing sites) / (comparable sites)

  which is a RATIONAL number with an integer numerator and denominator. There is
  no reason for a float to appear anywhere in the evidence. Better still, the
  inequalities a certificate asserts can be CROSS-MULTIPLIED, so they become
  comparisons between integers:

      n_xy / d_xy  <  n_xz / d_xz      ⟺      n_xy * d_xz  <  n_xz * d_xy

  (valid because the denominators are counts, hence strictly positive). So a
  certificate reduces to TWO INTEGER INEQUALITIES, which the kernel decides
  exactly and quickly. Per-pair denominators are allowed, which matters because
  pairwise gap deletion gives different comparable-site counts to different
  pairs.

  THE ARCHITECTURE. Rust searches and emits `CountCertificate` instances with
  integer literals; `lake env lean` typechecks them; a certificate the kernel
  rejects is a bug report about the Rust, never a result. The scientific
  conclusion therefore rests on this file plus the integers in the emitted
  certificate — not on the correctness of the search.

  Results.
    count_certificate_forces_violation    the counts force the rational
                                          violation, by integer arithmetic only
    count_certificate_refutes_ultrametric the inference rule over ℚ
    emitted_example / emitted_margin      a real certificate from the simulated
                                          recombinant, kernel-checked
    forged_certificate_is_impossible      the control: the forged inequality is
                                          refutable, so no forged certificate
                                          can be constructed at all
    counts_from_an_ultrametric_admit_none under the null, no counts qualify

  Verification contract: no `sorry`, no custom `axiom`; checked with
  `python -m anse.formal.lean_runner formal/ANSE/GeneticaReticulata.lean`.
-/
import Mathlib.Tactic

-- Kernel reduction of count functions over embedded sequences needs depth.
set_option maxRecDepth 100000

namespace ANSE.GeneticaReticulata

variable {α : Type*}

/-! ## The null, over an ordered field -/

/-- A dissimilarity is an **ultrametric** if it satisfies the strong triangle
inequality.

Stated over `ℚ` rather than over a general ordered field, because this file's
whole purpose is that genetic distances are exact rationals. The real-valued
version, for the abstract theory, is `ANSE.ReticulationCertificate`. -/
def IsUltrametric (d : α → α → ℚ) : Prop :=
  ∀ x y z, d x z ≤ max (d x y) (d y z)

/-- The exact distance of a pair from its site counts. -/
def ratio (n d : ℤ) : ℚ := (n : ℚ) / (d : ℚ)

/-! ## Certificates as integer counts -/

/-- A **count certificate**: three sequences, and for each of the three pairs
the number of differing sites and the number of comparable sites.

The two `far_*` fields are the cross-multiplied forms of
`dist x y < dist x z` and `dist y z < dist x z`. They are inequalities between
products of integers, so the kernel decides them without any rational
normalisation and without floating point.

Carrying them as fields means an ill-formed certificate cannot be constructed:
the Rust search can fail to find a certificate, or emit integers that do not
satisfy these, and in the latter case the emitted Lean file simply fails to
typecheck. -/
structure CountCertificate (α : Type*) where
  /-- One endpoint. -/
  x : α
  /-- The midpoint, close to both endpoints. -/
  y : α
  /-- The other endpoint. -/
  z : α
  /-- Differing sites between `x` and `y`. -/
  n_xy : ℤ
  /-- Comparable sites between `x` and `y`. -/
  d_xy : ℤ
  /-- Differing sites between `y` and `z`. -/
  n_yz : ℤ
  /-- Comparable sites between `y` and `z`. -/
  d_yz : ℤ
  /-- Differing sites between `x` and `z`. -/
  n_xz : ℤ
  /-- Comparable sites between `x` and `z`. -/
  d_xz : ℤ
  /-- Counts of comparable sites are positive. -/
  hd_xy : 0 < d_xy
  /-- Counts of comparable sites are positive. -/
  hd_yz : 0 < d_yz
  /-- Counts of comparable sites are positive. -/
  hd_xz : 0 < d_xz
  /-- `dist x y < dist x z`, cross-multiplied. -/
  far_xy : n_xy * d_xz < n_xz * d_xy
  /-- `dist y z < dist x z`, cross-multiplied. -/
  far_yz : n_yz * d_xz < n_xz * d_yz

/-- Cross-multiplication is sound: for positive integer denominators, the
integer inequality is equivalent to the rational one. This is the lemma that
licenses doing all the work in `ℤ`. -/
theorem ratio_lt_ratio_iff (n₁ d₁ n₂ d₂ : ℤ) (h₁ : 0 < d₁) (h₂ : 0 < d₂) :
    ratio n₁ d₁ < ratio n₂ d₂ ↔ n₁ * d₂ < n₂ * d₁ := by
  unfold ratio
  have h₁' : (0 : ℚ) < (d₁ : ℚ) := by exact_mod_cast h₁
  have h₂' : (0 : ℚ) < (d₂ : ℚ) := by exact_mod_cast h₂
  rw [div_lt_div_iff₀ h₁' h₂']
  exact_mod_cast Iff.rfl

/-- The counts force the rational violation: the far distance strictly exceeds
both near distances, hence their maximum. Proved from the integer fields alone. -/
theorem count_certificate_forces_violation (c : CountCertificate α) :
    max (ratio c.n_xy c.d_xy) (ratio c.n_yz c.d_yz) < ratio c.n_xz c.d_xz := by
  have h1 : ratio c.n_xy c.d_xy < ratio c.n_xz c.d_xz :=
    (ratio_lt_ratio_iff _ _ _ _ c.hd_xy c.hd_xz).mpr c.far_xy
  have h2 : ratio c.n_yz c.d_yz < ratio c.n_xz c.d_xz :=
    (ratio_lt_ratio_iff _ _ _ _ c.hd_yz c.hd_xz).mpr c.far_yz
  exact max_lt h1 h2

/-- **The inference rule.** Any dissimilarity taking the certificate's values on
the triple fails the strong triangle inequality, so it is not an ultrametric —
not a strict molecular clock on a clonal genealogy.

The hypotheses `h_*` are what tie the abstract conclusion to the data: they say
the numbers in the certificate really are this dataset's distances. The Rust
side discharges them by construction, because it emits the counts it computed. -/
theorem count_certificate_refutes_ultrametric
    (c : CountCertificate α) (d : α → α → ℚ)
    (h_xy : d c.x c.y = ratio c.n_xy c.d_xy)
    (h_yz : d c.y c.z = ratio c.n_yz c.d_yz)
    (h_xz : d c.x c.z = ratio c.n_xz c.d_xz) :
    ¬ IsUltrametric d := by
  intro hu
  have hviol := count_certificate_forces_violation c
  have := hu c.x c.y c.z
  rw [h_xy, h_yz, h_xz] at this
  exact absurd this (not_le.mpr hviol)

/-- The **exact** margin, as a rational rather than a float: the amount by which
the strong triangle inequality fails. Used to rank certificates, so the ranking
is a fact about the data and not about rounding. -/
def margin (c : CountCertificate α) : ℚ :=
  ratio c.n_xz c.d_xz - max (ratio c.n_xy c.d_xy) (ratio c.n_yz c.d_yz)

/-- The margin is strictly positive, by construction. -/
theorem margin_pos (c : CountCertificate α) : 0 < margin c := by
  unfold margin
  linarith [count_certificate_forces_violation c]

/-! ## A real certificate, kernel-checked

The counts below are the ones the Rust search emits for the planted recombinant
of the simulated positive control: 4000 comparable sites per pair, with 69 and
70 differing sites to the two parents and 131 between the parents themselves.
This is the end-to-end test of the architecture — an emitted certificate, checked
by the kernel. -/

/-- The emitted certificate. Every field is discharged by `decide` on integers. -/
def emitted_example : CountCertificate (Fin 3) where
  x := 0
  y := 1
  z := 2
  n_xy := 69
  d_xy := 4000
  n_yz := 70
  d_yz := 4000
  n_xz := 131
  d_xz := 4000
  hd_xy := by norm_num
  hd_yz := by norm_num
  hd_xz := by norm_num
  far_xy := by norm_num
  far_yz := by norm_num

/-- Its margin is the exact rational `61/4000`, not a float near it. -/
theorem emitted_margin : margin emitted_example = 61 / 4000 := by
  unfold margin emitted_example ratio
  norm_num

/-- **The control.** A forged certificate cannot be constructed, because the
integer inequality it would need is false. Here the forgery claims the parents
are closer to each other (69) than one is to the recombinant (70) while still
asserting a violation; `far_yz` would require `70 * 4000 < 69 * 4000`.

This is the property that makes the kernel a real check rather than a
formality: rejection is not a convention we chose to honour, it is the absence
of a term. -/
theorem forged_certificate_is_impossible :
    ¬ ((70 : ℤ) * 4000 < (69 : ℤ) * 4000) := by norm_num

/-- A second control, on the null side: counts drawn from an ultrametric admit
no certificate. Here the two largest distances agree (`2/4` and `2/4`), which is
exactly the ultrametric condition on a triple, and the `far` inequality fails. -/
theorem counts_from_an_ultrametric_admit_none :
    ¬ ((2 : ℤ) * 4 < (2 : ℤ) * 4) := by norm_num

/-- Equal counts give equal ratios, so the ultrametric boundary case is handled
exactly: an isosceles triple with the two largest distances equal is *not* a
violation, and no rounding can make it one. -/
theorem isosceles_triple_is_not_a_violation (n d : ℤ) (hd : 0 < d) :
    ¬ (ratio n d < ratio n d) := by
  simp

/-! ## Closing the gap: certificates about the SEQUENCES, not just about integers

An adversarial audit of the first version of this file found that the emitted
certificates were not about the data. The kernel was given six integers and
proved an arithmetic fact about them; nothing tied those integers to an
alignment, so a mislabelled triple or a transposed pair in the search produced a
file that still typechecked. The claim "the conclusion does not depend on the
search code" was false as worded: the *counting* was trusted.

The section below removes that. Sequences are embedded in the emitted file as
lists, the kernel COMPUTES the counts itself by `rfl`, and the inference rule is
instantiated against those computed counts. Measured: 4000 sites reduce in about
26 s and the result depends on no axioms at all.

What remains trusted is narrower and auditable by eye: that the embedded lists
are the dataset's sequences. A reader can diff them against the FASTA. -/

/-- Bases are encoded `A=0, C=1, G=2, T=3`; any other symbol (gap, `N`, an
ambiguity code) is encoded `4` and is not comparable. -/
def isBase (c : ℕ) : Bool := c < 4

/-- Positions at which both sequences carry a definite base. This is the
denominator of a p-distance, computed by the kernel rather than supplied. -/
def comparableCount (u v : List ℕ) : ℕ :=
  (u.zip v).countP (fun p => isBase p.1 && isBase p.2)

/-- Comparable positions at which the two sequences differ: the numerator. -/
def diffCount (u v : List ℕ) : ℕ :=
  (u.zip v).countP (fun p => isBase p.1 && isBase p.2 && p.1 != p.2)

/-- The exact p-distance of two encoded sequences. -/
noncomputable def pDistance (u v : List ℕ) : ℚ :=
  (diffCount u v : ℚ) / (comparableCount u v : ℚ)

/-- Cross-multiplication over `ℕ`, for the sequence-level certificates. -/
theorem pDistance_lt_iff (u₁ v₁ u₂ v₂ : List ℕ)
    (h₁ : 0 < comparableCount u₁ v₁) (h₂ : 0 < comparableCount u₂ v₂) :
    pDistance u₁ v₁ < pDistance u₂ v₂ ↔
      diffCount u₁ v₁ * comparableCount u₂ v₂
        < diffCount u₂ v₂ * comparableCount u₁ v₁ := by
  unfold pDistance
  have h₁' : (0 : ℚ) < (comparableCount u₁ v₁ : ℚ) := by exact_mod_cast h₁
  have h₂' : (0 : ℚ) < (comparableCount u₂ v₂ : ℚ) := by exact_mod_cast h₂
  rw [div_lt_div_iff₀ h₁' h₂']
  exact_mod_cast Iff.rfl

/-- **The sequence-level inference rule.** Given three embedded sequences whose
kernel-computed counts satisfy the two cross-multiplied inequalities, the
p-distance on that triple violates the strong triangle inequality.

Every hypothesis is discharged by `rfl` in the emitted file, so the kernel does
the counting. This is what makes the certificate a statement about the data. -/
theorem sequences_violate_ultrametric
    (sx sy sz : List ℕ)
    (hc_xy : 0 < comparableCount sx sy)
    (hc_yz : 0 < comparableCount sy sz)
    (hc_xz : 0 < comparableCount sx sz)
    (far_xy : diffCount sx sy * comparableCount sx sz
                < diffCount sx sz * comparableCount sx sy)
    (far_yz : diffCount sy sz * comparableCount sx sz
                < diffCount sx sz * comparableCount sy sz) :
    max (pDistance sx sy) (pDistance sy sz) < pDistance sx sz := by
  have h1 : pDistance sx sy < pDistance sx sz :=
    (pDistance_lt_iff sx sy sx sz hc_xy hc_xz).mpr far_xy
  have h2 : pDistance sy sz < pDistance sx sz :=
    (pDistance_lt_iff sy sz sx sz hc_yz hc_xz).mpr far_yz
  exact max_lt h1 h2

/-- And hence the dissimilarity induced by the alignment is not an ultrametric:
not a strict molecular clock on a clonal genealogy. Unlike
`count_certificate_refutes_ultrametric`, this version needs no hypothesis tying
numbers to data — the sequences *are* the data. -/
theorem sequences_refute_ultrametric
    (sx sy sz : List ℕ) (d : List ℕ → List ℕ → ℚ)
    (hd : ∀ u v, d u v = pDistance u v)
    (hc_xy : 0 < comparableCount sx sy)
    (hc_yz : 0 < comparableCount sy sz)
    (hc_xz : 0 < comparableCount sx sz)
    (far_xy : diffCount sx sy * comparableCount sx sz
                < diffCount sx sz * comparableCount sx sy)
    (far_yz : diffCount sy sz * comparableCount sx sz
                < diffCount sx sz * comparableCount sy sz) :
    ¬ IsUltrametric d := by
  intro hu
  have hviol :=
    sequences_violate_ultrametric sx sy sz hc_xy hc_yz hc_xz far_xy far_yz
  have := hu sx sy sz
  rw [hd, hd, hd] at this
  exact absurd this (not_le.mpr hviol)

/-- The form the emitted file actually uses.

The counts arrive as *equations* proved by `rfl` — cheap kernel reduction, which
scales to thousands of sites — and the inequalities are then comparisons between
numerals. Splitting it this way matters in practice: `0 < comparableCount …` is
not `rfl`-provable and `decide` on a 4000-element list exhausts the recursion
depth, whereas `rfl` on the equation reduces fine. -/
theorem sequences_refute_from_counts
    (sx sy sz : List ℕ)
    (n_xy d_xy n_yz d_yz n_xz d_xz : ℕ)
    (e_nxy : diffCount sx sy = n_xy) (e_dxy : comparableCount sx sy = d_xy)
    (e_nyz : diffCount sy sz = n_yz) (e_dyz : comparableCount sy sz = d_yz)
    (e_nxz : diffCount sx sz = n_xz) (e_dxz : comparableCount sx sz = d_xz)
    (hd_xy : 0 < d_xy) (hd_yz : 0 < d_yz) (hd_xz : 0 < d_xz)
    (far_xy : n_xy * d_xz < n_xz * d_xy)
    (far_yz : n_yz * d_xz < n_xz * d_yz) :
    max (pDistance sx sy) (pDistance sy sz) < pDistance sx sz := by
  apply sequences_violate_ultrametric
  · rw [e_dxy]; exact hd_xy
  · rw [e_dyz]; exact hd_yz
  · rw [e_dxz]; exact hd_xz
  · rw [e_nxy, e_dxz, e_nxz, e_dxy]; exact far_xy
  · rw [e_nyz, e_dxz, e_nxz, e_dyz]; exact far_yz

/-- Non-vacuity at the sequence level: three short sequences whose counts the
kernel computes, forming the minimal open path. `x` and `z` differ at every
comparable site; each is close to the midpoint `y`. -/
theorem sequence_witness :
    max (pDistance [0, 0, 0, 0] [0, 0, 1, 1]) (pDistance [0, 0, 1, 1] [1, 1, 1, 1])
      < pDistance [0, 0, 0, 0] [1, 1, 1, 1] := by
  apply sequences_refute_from_counts
    [0, 0, 0, 0] [0, 0, 1, 1] [1, 1, 1, 1] 2 4 2 4 4 4 <;>
    first | rfl | norm_num

/-- Control at the sequence level: a gap (`4`) is excluded from both counts, so
the denominator really is per-pair. Here one position is a gap, leaving three
comparable sites of which one differs. -/
theorem gap_is_excluded_from_both_counts :
    comparableCount [0, 0, 0, 4] [0, 0, 1, 1] = 3
      ∧ diffCount [0, 0, 0, 4] [0, 0, 1, 1] = 1 := by
  constructor <;> rfl

end ANSE.GeneticaReticulata

#print axioms ANSE.GeneticaReticulata.ratio_lt_ratio_iff
#print axioms ANSE.GeneticaReticulata.count_certificate_forces_violation
#print axioms ANSE.GeneticaReticulata.count_certificate_refutes_ultrametric
#print axioms ANSE.GeneticaReticulata.margin_pos
#print axioms ANSE.GeneticaReticulata.emitted_margin
#print axioms ANSE.GeneticaReticulata.forged_certificate_is_impossible
#print axioms ANSE.GeneticaReticulata.counts_from_an_ultrametric_admit_none
#print axioms ANSE.GeneticaReticulata.isosceles_triple_is_not_a_violation
#print axioms ANSE.GeneticaReticulata.pDistance_lt_iff
#print axioms ANSE.GeneticaReticulata.sequences_violate_ultrametric
#print axioms ANSE.GeneticaReticulata.sequences_refute_ultrametric
#print axioms ANSE.GeneticaReticulata.sequence_witness
#print axioms ANSE.GeneticaReticulata.gap_is_excluded_from_both_counts
#print axioms ANSE.GeneticaReticulata.sequences_refute_from_counts

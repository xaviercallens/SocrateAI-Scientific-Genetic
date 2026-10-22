/-
  ANSE.GeneticaReticulataTree — the tree-level null, and complete deletion.

  Phase 2/4 of docs/GENETICA_RETICULATA.md. Two problems with the ultrametric
  certificate of `ANSE.GeneticaReticulata` are addressed here, and both were
  named by that file's own adversarial audit rather than discovered later.

  PROBLEM 1 — THE NULL IS TOO STRONG TO BE INTERESTING ON ITS OWN.
  `IsUltrametric` encodes "clonal genealogy AND strict molecular clock". Real
  clocks are not strict, so a violation is ambiguous between reticulation and
  rate variation. The classical weaker null is the FOUR-POINT CONDITION, which
  encodes "tree" without "clock": branch lengths may vary freely. Buneman's
  theorem (1974) says a metric is realised by a weighted tree exactly when every
  quadruple satisfies it, and that step is CITED here, not proved — the same
  epistemic status the ultrametric file gives to "strict clock ⇒ ultrametric".

  What IS proved here is the nesting, `fourPoint_of_ultrametric`, and therefore
  its contrapositive `quartet_violation_refutes_ultrametric`: a quartet
  certificate refutes strictly more than a triple certificate does. That is the
  formal content of the experimental design — a signal that survives the
  relaxation of the clock is not explained by the clock.

  PROBLEM 2 — PAIRWISE DELETION OF MISSING DATA CAN MANUFACTURE A CERTIFICATE.
  Finding 6 of the audit: two pairs compared over different sets of columns have
  incommensurable distances, and that alone can break the ultrametric inequality
  on an exactly clock-like genealogy. The fix is complete deletion WITHIN the
  triple (or quadruple): count only columns where every member carries a base.
  Then all distances in the certificate share one denominator, the certificate
  becomes a comparison between integers with no cross-multiplication at all, and
  the confound cannot arise.

  `pairwise_deletion_manufactures_what_complete_deletion_denies` is the kernel's
  own witness that the two policies really do disagree: three explicit sequences
  on which pairwise deletion yields a valid violation and complete deletion
  yields none. The deletion policy is therefore part of the claim, not an
  implementation detail, and this file makes that checkable rather than asserted.

  Results.
    fourPoint_of_ultrametric          every ultrametric satisfies the four-point
                                      condition: the nulls are nested
    quartet_violation_refutes_ultrametric   so a quartet certificate is strictly
                                      stronger evidence than a triple one
    quartet_counts_force_violation    six integers force the rational violation
    quartet_refutes_four_point        the inference rule
    quartet_refute_from_counts        the emitted form; kernel counts the data
    complete_triple_refute_from_counts  the triple test under complete deletion
    pairwise_deletion_manufactures_…  finding 6, witnessed in the kernel
    tree_quartet_admits_no_certificate  the control: a tree quartet yields none
    reciprocal_recombinants_witness   non-vacuity: the mosaic signature

  Verification contract: no `sorry`, no custom `axiom`; checked with
  `python -m anse.formal.lean_runner formal/ANSE/GeneticaReticulataTree.lean`.
-/
import ANSE.GeneticaReticulata

-- Kernel reduction of count functions over embedded sequences needs depth.
set_option maxRecDepth 100000
-- The four-point nesting is proved by case analysis over 256 branches.
set_option maxHeartbeats 1000000

namespace ANSE.GeneticaReticulata

variable {α : Type*}

/-! ## The tree-level null, and its nesting inside the clock-level null -/

/-- A dissimilarity satisfies the **four-point condition** if, for every
quadruple, no pairing's sum of distances strictly exceeds both of the other two
pairings' sums.

This is the tree-level null. Unlike `IsUltrametric` it says nothing about rates:
a weighted tree with wildly unequal branch lengths satisfies it. By Buneman's
theorem a metric satisfying it on every quadruple is realised by a weighted tree
— **cited, not proved here**, exactly as "a strict clock on a clonal genealogy
gives an ultrametric" is cited in `ANSE.GeneticaReticulata`. -/
def IsFourPoint (d : α → α → ℚ) : Prop :=
  ∀ a b c e, d a b + d c e ≤ max (d a c + d b e) (d a e + d b c)

/-- **The nulls are nested.** Every symmetric ultrametric satisfies the
four-point condition.

Proved rather than cited, because the whole experimental design rests on it: it
is what makes "the ultrametric test fires but the four-point test does not" a
coherent reading (rate variation on a tree) instead of a contradiction. -/
theorem fourPoint_of_ultrametric {d : α → α → ℚ}
    (hsym : ∀ x y, d x y = d y x) (hu : IsUltrametric d) : IsFourPoint d := by
  intro a b c e
  have h1 : d a b ≤ max (d a c) (d b c) := by
    have := hu a c b; rwa [hsym c b] at this
  have h2 : d a b ≤ max (d a e) (d b e) := by
    have := hu a e b; rwa [hsym e b] at this
  have h3 : d c e ≤ max (d a c) (d a e) := by
    have := hu c a e; rwa [hsym c a] at this
  have h4 : d c e ≤ max (d b c) (d b e) := by
    have := hu c b e; rwa [hsym c b] at this
  rw [le_max_iff] at h1 h2 h3 h4
  rw [le_max_iff]
  rcases h1 with h1 | h1 <;> rcases h2 with h2 | h2 <;>
    rcases h3 with h3 | h3 <;> rcases h4 with h4 | h4 <;>
    rcases le_total (d a c) (d a e) with k2 | k2 <;>
    rcases le_total (d a c) (d b c) with k3 | k3 <;>
    rcases le_total (d b e) (d a e) with k4 | k4 <;>
    rcases le_total (d b e) (d b c) with k5 | k5 <;>
    first
      | (left; linarith)
      | (right; linarith)

/-- **The contrapositive, and the point of the whole file.** A quadruple that
refutes the four-point condition also refutes ultrametricity, so a quartet
certificate is strictly stronger evidence than a triple certificate: it cannot
be explained by relaxing the molecular clock, because the four-point condition
never assumed one.

The converse fails, and that asymmetry is the experiment: a dataset whose triple
certificates are large while its quartet certificates stay at the noise floor is
evidence of rate variation on a tree, not of reticulation. -/
theorem quartet_violation_refutes_ultrametric {d : α → α → ℚ}
    (hsym : ∀ x y, d x y = d y x) (h : ¬ IsFourPoint d) : ¬ IsUltrametric d :=
  fun hu => h (fourPoint_of_ultrametric hsym hu)

/-! ## Complete deletion: one denominator, hence integer certificates

Under complete deletion within the subset, every distance in a certificate is
`(differing sites) / D` for the *same* `D`. The certificate's inequalities then
need no cross-multiplication: they are comparisons between the numerators. This
is both a smaller trust surface and the fix for the audit's finding 6. -/

/-- A distance under complete deletion: differing sites over the subset's own
common count of comparable sites. -/
noncomputable def ratioN (n D : ℕ) : ℚ := (n : ℚ) / (D : ℚ)

/-- With a common positive denominator, the rational comparison *is* the integer
comparison. No cross-product appears. -/
theorem ratioN_lt_ratioN_iff (a b D : ℕ) (hD : 0 < D) :
    ratioN a D < ratioN b D ↔ a < b := by
  unfold ratioN
  have hD' : (0 : ℚ) < (D : ℚ) := by exact_mod_cast hD
  rw [div_lt_div_iff₀ hD' hD']
  constructor
  · intro h
    have : (a : ℚ) < (b : ℚ) := lt_of_mul_lt_mul_right h (le_of_lt hD')
    exact_mod_cast this
  · intro h
    have h' : (a : ℚ) < (b : ℚ) := by exact_mod_cast h
    exact mul_lt_mul_of_pos_right h' hD'

/-- Sums share the denominator too, which is what makes the four-point test an
integer test. -/
theorem ratioN_add (a b D : ℕ) : ratioN a D + ratioN b D = ratioN (a + b) D := by
  unfold ratioN
  push_cast
  ring

/-! ### Counting over a triple, complete-deletion style -/

/-- Positions where all three sequences carry a definite base. The single
denominator of a complete-deletion triple. -/
def comparable3 (u v w : List ℕ) : ℕ :=
  (u.zip (v.zip w)).countP (fun p => isBase p.1 && isBase p.2.1 && isBase p.2.2)

/-- Differing sites between the first and second sequence, counted only where all
three are comparable. -/
def diff3xy (u v w : List ℕ) : ℕ :=
  (u.zip (v.zip w)).countP
    (fun p => isBase p.1 && isBase p.2.1 && isBase p.2.2 && p.1 != p.2.1)

/-- Differing sites between the second and third, over the same columns. -/
def diff3yz (u v w : List ℕ) : ℕ :=
  (u.zip (v.zip w)).countP
    (fun p => isBase p.1 && isBase p.2.1 && isBase p.2.2 && p.2.1 != p.2.2)

/-- Differing sites between the first and third, over the same columns. -/
def diff3xz (u v w : List ℕ) : ℕ :=
  (u.zip (v.zip w)).countP
    (fun p => isBase p.1 && isBase p.2.1 && isBase p.2.2 && p.1 != p.2.2)

/-- **The complete-deletion triple certificate.** Two integer comparisons on one
denominator force the ultrametric violation. Compare
`count_certificate_forces_violation`, which needs cross-products because
pairwise deletion gives the three pairs three different denominators. -/
theorem complete_triple_forces_violation (n_xy n_yz n_xz D : ℕ) (hD : 0 < D)
    (far_xy : n_xy < n_xz) (far_yz : n_yz < n_xz) :
    max (ratioN n_xy D) (ratioN n_yz D) < ratioN n_xz D :=
  max_lt ((ratioN_lt_ratioN_iff _ _ _ hD).mpr far_xy)
    ((ratioN_lt_ratioN_iff _ _ _ hD).mpr far_yz)

/-- The emitted form for a complete-deletion triple: the counts arrive as
equations the kernel discharges by `rfl`, so the kernel does the counting, and
the inequalities are comparisons between numerals. -/
theorem complete_triple_refute_from_counts
    (sx sy sz : List ℕ) (n_xy n_yz n_xz D : ℕ)
    (e_xy : diff3xy sx sy sz = n_xy)
    (e_yz : diff3yz sx sy sz = n_yz)
    (e_xz : diff3xz sx sy sz = n_xz)
    (e_D : comparable3 sx sy sz = D)
    (hD : 0 < D) (far_xy : n_xy < n_xz) (far_yz : n_yz < n_xz) :
    max (ratioN (diff3xy sx sy sz) (comparable3 sx sy sz))
        (ratioN (diff3yz sx sy sz) (comparable3 sx sy sz))
      < ratioN (diff3xz sx sy sz) (comparable3 sx sy sz) := by
  rw [e_xy, e_yz, e_xz, e_D]
  exact complete_triple_forces_violation n_xy n_yz n_xz D hD far_xy far_yz

/-- And hence no ultrametric takes these values on the triple. -/
theorem complete_triple_refutes_ultrametric
    (d : α → α → ℚ) (x y z : α) (n_xy n_yz n_xz D : ℕ) (hD : 0 < D)
    (h_xy : d x y = ratioN n_xy D) (h_yz : d y z = ratioN n_yz D)
    (h_xz : d x z = ratioN n_xz D)
    (far_xy : n_xy < n_xz) (far_yz : n_yz < n_xz) :
    ¬ IsUltrametric d := by
  intro hu
  have hviol := complete_triple_forces_violation n_xy n_yz n_xz D hD far_xy far_yz
  have := hu x y z
  rw [h_xy, h_yz, h_xz] at this
  exact absurd this (not_le.mpr hviol)

/-! ### Finding 6, witnessed in the kernel

Three sequences on which the two deletion policies disagree. Columns 7 and 8 are
gapped in `y` only, so they enter the `x`–`z` comparison and neither of the
others: pairwise deletion then reports `d x z = 4/8` against `d x y = d y z =
2/6`, a valid violation. Restricted to the six columns where all three carry a
base, all three distances are `2/6` and there is no violation at all.

Both halves are computed by the kernel from the sequences, so this is not an
illustration of the audit's finding — it is the finding, checked. -/

/-- First sequence of the confounded triple. -/
def confoundX : List ℕ := [0, 0, 0, 0, 0, 0, 0, 0]

/-- Second: gapped (`4`) at the last two positions. -/
def confoundY : List ℕ := [1, 1, 0, 0, 0, 0, 4, 4]

/-- Third: differs from `confoundX` at exactly those two gapped positions. -/
def confoundZ : List ℕ := [1, 0, 1, 0, 0, 0, 1, 1]

/-- **Pairwise deletion manufactures a violation that complete deletion denies.**

The first conjunct is a genuine pairwise-deletion certificate — the kernel
accepts it, and would accept it on real data. The remaining conjuncts show the
complete-deletion counts are `2, 2, 2` out of `6`, so the far distance does not
strictly exceed the near ones and no certificate exists.

Consequence for the method, stated plainly: a kernel-accepted pairwise-deletion
certificate on gappy data may be a true statement about the distance matrix and
a false statement about the genealogy. The kernel cannot detect that, because the
inequality it is given is true. Complete deletion within the triple is therefore
the default for gappy input, not an option. -/
theorem pairwise_deletion_manufactures_what_complete_deletion_denies :
    max (pDistance confoundX confoundY) (pDistance confoundY confoundZ)
        < pDistance confoundX confoundZ
    ∧ comparable3 confoundX confoundY confoundZ = 6
    ∧ diff3xy confoundX confoundY confoundZ = 2
    ∧ diff3yz confoundX confoundY confoundZ = 2
    ∧ diff3xz confoundX confoundY confoundZ = 2
    ∧ ¬ (max (ratioN 2 6) (ratioN 2 6) < ratioN 2 6) := by
  refine ⟨?_, by rfl, by rfl, by rfl, by rfl, ?_⟩
  · apply sequences_refute_from_counts confoundX confoundY confoundZ 2 6 2 6 4 8 <;>
      first | rfl | norm_num
  · simp

/-! ## Quartet certificates: the four-point test as integer comparisons -/

/-- The shape a four-point violation must have: the proposed far pairing's sum
strictly exceeds both others. Named so the control below can state its absence
without repeating the inequalities. -/
def fourPointViolates (s_far s_b s_c : ℕ) : Prop := s_b < s_far ∧ s_c < s_far

/-- **Six integers force the four-point violation.** Complete deletion within the
quadruple gives all six pairs one denominator, so the sums compare as integers. -/
theorem quartet_counts_force_violation
    (n_ab n_ce n_ac n_be n_ae n_bc D : ℕ) (hD : 0 < D)
    (h : fourPointViolates (n_ab + n_ce) (n_ac + n_be) (n_ae + n_bc)) :
    max (ratioN n_ac D + ratioN n_be D) (ratioN n_ae D + ratioN n_bc D)
      < ratioN n_ab D + ratioN n_ce D := by
  rw [ratioN_add, ratioN_add, ratioN_add]
  exact max_lt ((ratioN_lt_ratioN_iff _ _ _ hD).mpr h.1)
    ((ratioN_lt_ratioN_iff _ _ _ hD).mpr h.2)

/-- **The inference rule at the tree level.** A dissimilarity taking the
certificate's values on the quadruple fails the four-point condition, so it is
not realised by any weighted tree (Buneman) — regardless of how the rates vary
along that tree. -/
theorem quartet_refutes_four_point
    (d : α → α → ℚ) (a b c e : α) (n_ab n_ce n_ac n_be n_ae n_bc D : ℕ) (hD : 0 < D)
    (h_ab : d a b = ratioN n_ab D) (h_ce : d c e = ratioN n_ce D)
    (h_ac : d a c = ratioN n_ac D) (h_be : d b e = ratioN n_be D)
    (h_ae : d a e = ratioN n_ae D) (h_bc : d b c = ratioN n_bc D)
    (h : fourPointViolates (n_ab + n_ce) (n_ac + n_be) (n_ae + n_bc)) :
    ¬ IsFourPoint d := by
  intro hfp
  have hviol := quartet_counts_force_violation n_ab n_ce n_ac n_be n_ae n_bc D hD h
  have := hfp a b c e
  rw [h_ab, h_ce, h_ac, h_be, h_ae, h_bc] at this
  exact absurd this (not_le.mpr hviol)

/-- And therefore, by the nesting, it is not an ultrametric either: one quartet
certificate refutes both nulls at once. -/
theorem quartet_refutes_both_nulls
    (d : α → α → ℚ) (hsym : ∀ x y, d x y = d y x)
    (a b c e : α) (n_ab n_ce n_ac n_be n_ae n_bc D : ℕ) (hD : 0 < D)
    (h_ab : d a b = ratioN n_ab D) (h_ce : d c e = ratioN n_ce D)
    (h_ac : d a c = ratioN n_ac D) (h_be : d b e = ratioN n_be D)
    (h_ae : d a e = ratioN n_ae D) (h_bc : d b c = ratioN n_bc D)
    (h : fourPointViolates (n_ab + n_ce) (n_ac + n_be) (n_ae + n_bc)) :
    ¬ IsFourPoint d ∧ ¬ IsUltrametric d :=
  let hfp := quartet_refutes_four_point d a b c e n_ab n_ce n_ac n_be n_ae n_bc D hD
    h_ab h_ce h_ac h_be h_ae h_bc h
  ⟨hfp, quartet_violation_refutes_ultrametric hsym hfp⟩

/-! ### Counting over a quadruple -/

/-- Positions where all four sequences carry a definite base: the quadruple's
single denominator. -/
def comparable4 (u v w t : List ℕ) : ℕ :=
  (u.zip (v.zip (w.zip t))).countP
    (fun p => isBase p.1 && isBase p.2.1 && isBase p.2.2.1 && isBase p.2.2.2)

/-- Differing sites for the pair `(1,2)` over quadruple-complete columns. -/
def diff4ab (u v w t : List ℕ) : ℕ :=
  (u.zip (v.zip (w.zip t))).countP
    (fun p => isBase p.1 && isBase p.2.1 && isBase p.2.2.1 && isBase p.2.2.2
              && p.1 != p.2.1)

/-- Differing sites for the pair `(3,4)`. -/
def diff4ce (u v w t : List ℕ) : ℕ :=
  (u.zip (v.zip (w.zip t))).countP
    (fun p => isBase p.1 && isBase p.2.1 && isBase p.2.2.1 && isBase p.2.2.2
              && p.2.2.1 != p.2.2.2)

/-- Differing sites for the pair `(1,3)`. -/
def diff4ac (u v w t : List ℕ) : ℕ :=
  (u.zip (v.zip (w.zip t))).countP
    (fun p => isBase p.1 && isBase p.2.1 && isBase p.2.2.1 && isBase p.2.2.2
              && p.1 != p.2.2.1)

/-- Differing sites for the pair `(2,4)`. -/
def diff4be (u v w t : List ℕ) : ℕ :=
  (u.zip (v.zip (w.zip t))).countP
    (fun p => isBase p.1 && isBase p.2.1 && isBase p.2.2.1 && isBase p.2.2.2
              && p.2.1 != p.2.2.2)

/-- Differing sites for the pair `(1,4)`. -/
def diff4ae (u v w t : List ℕ) : ℕ :=
  (u.zip (v.zip (w.zip t))).countP
    (fun p => isBase p.1 && isBase p.2.1 && isBase p.2.2.1 && isBase p.2.2.2
              && p.1 != p.2.2.2)

/-- Differing sites for the pair `(2,3)`. -/
def diff4bc (u v w t : List ℕ) : ℕ :=
  (u.zip (v.zip (w.zip t))).countP
    (fun p => isBase p.1 && isBase p.2.1 && isBase p.2.2.1 && isBase p.2.2.2
              && p.2.1 != p.2.2.1)

/-- **The emitted quartet form.** Seven counts arrive as `rfl`-provable equations
— so the kernel reads the alignment — and the verdict is one comparison between
sums of numerals. -/
theorem quartet_refute_from_counts
    (sa sb sc se : List ℕ) (n_ab n_ce n_ac n_be n_ae n_bc D : ℕ)
    (e_ab : diff4ab sa sb sc se = n_ab) (e_ce : diff4ce sa sb sc se = n_ce)
    (e_ac : diff4ac sa sb sc se = n_ac) (e_be : diff4be sa sb sc se = n_be)
    (e_ae : diff4ae sa sb sc se = n_ae) (e_bc : diff4bc sa sb sc se = n_bc)
    (e_D : comparable4 sa sb sc se = D) (hD : 0 < D)
    (h : fourPointViolates (n_ab + n_ce) (n_ac + n_be) (n_ae + n_bc)) :
    max (ratioN (diff4ac sa sb sc se) (comparable4 sa sb sc se)
          + ratioN (diff4be sa sb sc se) (comparable4 sa sb sc se))
        (ratioN (diff4ae sa sb sc se) (comparable4 sa sb sc se)
          + ratioN (diff4bc sa sb sc se) (comparable4 sa sb sc se))
      < ratioN (diff4ab sa sb sc se) (comparable4 sa sb sc se)
          + ratioN (diff4ce sa sb sc se) (comparable4 sa sb sc se) := by
  rw [e_ab, e_ce, e_ac, e_be, e_ae, e_bc, e_D]
  exact quartet_counts_force_violation n_ab n_ce n_ac n_be n_ae n_bc D hD h

/-! ### Non-vacuity and controls -/

/-- Parent `a`. -/
def mosaicA : List ℕ := [0, 0, 0, 0]
/-- Parent `b`, differing from `a` at every site. -/
def mosaicB : List ℕ := [1, 1, 1, 1]
/-- A mosaic of the two parents: first half from `b`, second from `a`. -/
def mosaicC : List ℕ := [1, 1, 0, 0]
/-- The reciprocal mosaic: first half from `a`, second from `b`. -/
def mosaicE : List ℕ := [0, 0, 1, 1]

/-- **Non-vacuity, and the biological signature.** Two parents and their two
reciprocal mosaics. Every count is computed by the kernel. The parent pairing's
sum is `8/4` while both mosaic pairings give `4/4`, so the four-point condition
fails by the exact rational `1`.

This is the shape recombination makes and a tree cannot: each mosaic is close to
both parents, so the two "crossing" pairings are cheap while the parents stay
far apart. -/
theorem reciprocal_recombinants_witness :
    comparable4 mosaicA mosaicB mosaicC mosaicE = 4
    ∧ diff4ab mosaicA mosaicB mosaicC mosaicE = 4
    ∧ diff4ce mosaicA mosaicB mosaicC mosaicE = 4
    ∧ diff4ac mosaicA mosaicB mosaicC mosaicE = 2
    ∧ diff4be mosaicA mosaicB mosaicC mosaicE = 2
    ∧ diff4ae mosaicA mosaicB mosaicC mosaicE = 2
    ∧ diff4bc mosaicA mosaicB mosaicC mosaicE = 2
    ∧ max (ratioN 2 4 + ratioN 2 4) (ratioN 2 4 + ratioN 2 4)
        < ratioN 4 4 + ratioN 4 4 := by
  refine ⟨by rfl, by rfl, by rfl, by rfl, by rfl, by rfl, by rfl, ?_⟩
  apply quartet_counts_force_violation 4 4 2 2 2 2 4 (by norm_num)
  exact ⟨by norm_num, by norm_num⟩

/-- Parent-side leaf of a tree quartet. -/
def treeA : List ℕ := [0, 0, 0, 0, 0, 0]
/-- Its sister, one substitution away. -/
def treeB : List ℕ := [1, 0, 0, 0, 0, 0]
/-- A leaf of the far clade. -/
def treeC : List ℕ := [0, 0, 1, 1, 1, 1]
/-- Its sister. -/
def treeE : List ℕ := [0, 1, 1, 1, 1, 1]

/-- **The control.** A quartet drawn from a tree admits no certificate, whichever
pairing is proposed as the far one. The kernel computes the three pairing sums as
`2, 10, 10`: the two largest are equal, which is exactly the four-point
condition, and `fourPointViolates` fails for every choice.

Without this the witness above would be uninformative — a test that fires on tree
data as well would have no content. -/
theorem tree_quartet_admits_no_certificate :
    diff4ab treeA treeB treeC treeE + diff4ce treeA treeB treeC treeE = 2
    ∧ diff4ac treeA treeB treeC treeE + diff4be treeA treeB treeC treeE = 10
    ∧ diff4ae treeA treeB treeC treeE + diff4bc treeA treeB treeC treeE = 10
    ∧ ¬ fourPointViolates 2 10 10
    ∧ ¬ fourPointViolates 10 2 10
    ∧ ¬ fourPointViolates 10 10 2 := by
  refine ⟨by rfl, by rfl, by rfl, ?_, ?_, ?_⟩ <;>
    rintro ⟨h1, h2⟩ <;> omega

/-- A forged quartet certificate cannot be constructed: the inequality it would
need is refutable. Rejection is the absence of a term, not a convention. -/
theorem forged_quartet_is_impossible : ¬ fourPointViolates 4 4 2 := by
  rintro ⟨h1, _⟩
  omega

/-- A gap in any one of the four sequences removes the column from all six
counts, so the quadruple really does share one denominator. Here position 4 is
gapped in the third sequence, leaving three comparable columns. -/
theorem a_gap_anywhere_shrinks_the_common_denominator :
    comparable4 [0, 0, 0, 0] [1, 1, 1, 1] [1, 1, 0, 4] [0, 0, 1, 1] = 3
    ∧ diff4ab [0, 0, 0, 0] [1, 1, 1, 1] [1, 1, 0, 4] [0, 0, 1, 1] = 3 := by
  constructor <;> rfl

/-!
### The Gromov product: the hinge of Buneman's converse

Buneman's theorem is cited in this file, not proved. The constructive proof in
the literature does not attack the tree directly; it factors through the
**Gromov product** at a chosen root, reducing the four-point condition on `d` to
the ultrametric condition on a derived function, and only then builds an
equidistant tree by cherry-picking. The reduction is the hinge, and unlike the
tree construction it is pure algebra, so it can be proved here with no appeal to
anything.

Define, at a root `r`,
`gromov d r x y = (d x y - d r x - d r y) / 2`.
This is the negated Gromov product, a dissimilarity rather than a similarity,
and it is *not* a metric: it takes negative values. What matters is only whether
it satisfies the strong (max) triangle inequality, so that condition is named
separately rather than reusing `IsUltrametric`, which would wrongly suggest
non-negativity.

`fourPoint_iff_gromov_strongTriangle` is the equivalence. Both directions are
relabelling plus symmetry — set `(a,b,c,e) := (x,z,r,y)` and read it off — and
the kernel is a better judge of the sign conventions than prose is, which is the
reason to state it here rather than cite it.

What this does and does not buy. It does *not* prove Buneman's theorem: the
equidistant-tree construction, and the re-rooting that turns it back into a tree
metric, are still ahead. It does replace the first step of that proof with a
machine-checked one, and it connects this development to the ultrametric notion
the rest of the file is built on — the two nulls of §10.2 turn out to be the
same condition seen from two places, which is a sharper statement of the nesting
than `fourPoint_of_ultrametric` alone.
-/

/-- The negated Gromov product at a root. Not a metric: it may be negative. -/
noncomputable def gromov (d : α → α → ℚ) (r x y : α) : ℚ :=
  (d x y - d r x - d r y) / 2

/-- The strong (max) triangle inequality, with no positivity assumed.

Stated separately from `IsUltrametric` on purpose: `gromov` takes negative
values, so calling it an ultrametric would import a hypothesis it does not
satisfy. The literature's "the Gromov product is an ultrametric" means exactly
this inequality. -/
def StrongTriangle (f : α → α → ℚ) : Prop :=
  ∀ x y z, f x z ≤ max (f x y) (f y z)

/-- **The hinge of Buneman's converse.** A symmetric `d` satisfies the four-point
condition exactly when its Gromov product at every root satisfies the strong
triangle inequality.

The tree construction that follows from the right-hand side is still cited; this
equivalence is not. -/
theorem fourPoint_iff_gromov_strongTriangle {d : α → α → ℚ}
    (hsym : ∀ x y, d x y = d y x) :
    IsFourPoint d ↔ ∀ r, StrongTriangle (gromov d r) := by
  constructor
  · intro h r x y z
    -- Instantiate the four-point condition at (x, z, r, y).
    have key := h x z r y
    unfold gromov
    rw [le_max_iff]
    rw [le_max_iff] at key
    rcases key with key | key
    · right
      rw [hsym x r] at key
      rw [hsym z y] at key
      linarith
    · left
      rw [hsym z r] at key
      linarith
  · intro h a b c e
    -- Read the four-point condition back off the root c, with the remaining
    -- three in the order (a, e, b): the Gromov product at c of the pair (a, b)
    -- is what carries `d a b + d c e`.
    have key := h c a e b
    unfold gromov at key
    rw [le_max_iff] at key
    rw [le_max_iff]
    rcases key with key | key
    · right
      rw [hsym c b] at key
      linarith
    · left
      rw [hsym e b] at key
      rw [hsym c a] at key
      linarith

/-- The direction the certificate uses: a quartet violation gives a root at which
the Gromov product fails the strong triangle inequality. -/
theorem gromov_strongTriangle_fails_of_quartet_violation {d : α → α → ℚ}
    (hsym : ∀ x y, d x y = d y x) (h : ¬ IsFourPoint d) :
    ¬ (∀ r, StrongTriangle (gromov d r)) :=
  fun hg => h ((fourPoint_iff_gromov_strongTriangle hsym).mpr hg)

end ANSE.GeneticaReticulata

#print axioms ANSE.GeneticaReticulata.fourPoint_of_ultrametric
#print axioms ANSE.GeneticaReticulata.quartet_violation_refutes_ultrametric
#print axioms ANSE.GeneticaReticulata.ratioN_lt_ratioN_iff
#print axioms ANSE.GeneticaReticulata.ratioN_add
#print axioms ANSE.GeneticaReticulata.complete_triple_forces_violation
#print axioms ANSE.GeneticaReticulata.complete_triple_refute_from_counts
#print axioms ANSE.GeneticaReticulata.complete_triple_refutes_ultrametric
#print axioms ANSE.GeneticaReticulata.pairwise_deletion_manufactures_what_complete_deletion_denies
#print axioms ANSE.GeneticaReticulata.quartet_counts_force_violation
#print axioms ANSE.GeneticaReticulata.quartet_refutes_four_point
#print axioms ANSE.GeneticaReticulata.quartet_refutes_both_nulls
#print axioms ANSE.GeneticaReticulata.quartet_refute_from_counts
#print axioms ANSE.GeneticaReticulata.reciprocal_recombinants_witness
#print axioms ANSE.GeneticaReticulata.tree_quartet_admits_no_certificate
#print axioms ANSE.GeneticaReticulata.forged_quartet_is_impossible
#print axioms ANSE.GeneticaReticulata.a_gap_anywhere_shrinks_the_common_denominator
#print axioms ANSE.GeneticaReticulata.fourPoint_iff_gromov_strongTriangle
#print axioms ANSE.GeneticaReticulata.gromov_strongTriangle_fails_of_quartet_violation

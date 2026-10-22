/-
  ANSE.ReticulationCertificate — a replacement for the failed topological
  approach of Section 8, with a null that is proved rather than shuffled.

  WHY THE OLD APPROACH FAILED. Section 8 embedded ONE sequence by the chaos game
  and measured the topology of the resulting cloud. The chaos game contracts any
  sequence onto a self-similar attractor, so the holes measured were the
  attractor's, not the sequence's. Having no mathematical null, we had to build a
  computational one out of shuffles -- which cost us a p-value floor, a
  frame-sampling artefact (§11.3), and a statistic with no proved meaning. The
  honest result was negative.

  THE REPLACEMENT. Measure the topology of a POPULATION of sequences under a
  genetic-distance metric, in the manner of Chan, Carlsson and Rabadan (PNAS
  2013). Here H_1 has a meaning: loops indicate non-tree-like (reticulate)
  evolution -- recombination, reassortment, horizontal transfer. And crucially,
  the null is a THEOREM rather than an ensemble of shuffles:

      an ultrametric has no 1-cycles in its Vietoris-Rips filtration,
      at ANY scale.

  So a cycle is not "surprising relative to shuffles"; it is a CERTIFICATE that
  the data is not ultrametric. No p-value, no floor, no shuffle artefacts.

  Biologically, an ultrametric is exactly a strict molecular clock on a clonal
  (tree) genealogy: all pairwise distances are determined by divergence times
  from a common ancestor, so the two largest distances in any triple are equal.

  THE ARCHITECTURE this file is built for. The Rust crate searches for
  certificates; this file proves the inference rule that turns a certificate into
  a conclusion; and a certificate is three sequence labels and three numbers, so
  a reader can recheck it by hand. The scientific claim therefore does NOT
  depend on trusting the search code -- only on the arithmetic of the triple it
  returns. A bug in the Rust can produce a non-certificate, which fails to
  validate; it cannot produce a false conclusion that validates.

  Results.
    ultrametric_iff_no_positive_margin   defect = 0 iff ultrametric. The
                                         statistic is parameter-free and exact.
    ultrametric_no_induced_path          at every scale, two close pairs force
                                         the third close: the Rips graph is a
                                         disjoint union of cliques, so no cycle
                                         can be unfilled and beta_1 = 0.
    Certificate / certificate_refutes_ultrametric
                                         the inference rule, with a strictly
                                         positive margin naming the violation.
    certificate_exists_nontrivially      non-vacuity: a certificate exists, so
                                         the rule is not about an empty class.
    ultrametric_admits_no_certificate    the control: under the null nothing is
                                         findable, at any scale. A test that
                                         fired on ultrametric data would be
                                         measuring its own pipeline.
    abs_not_ultrametric                  the condition is strictly stronger than
                                         being a metric, so the null has content.

  Verification contract: no `sorry`, no custom `axiom`; checked with
  `python -m anse.formal.lean_runner formal/ANSE/ReticulationCertificate.lean`.
-/
import Mathlib.Tactic

namespace ANSE.ReticulationCertificate

variable {α : Type*}

/-! ## The null: ultrametricity -/

/-- A dissimilarity is an **ultrametric** if it satisfies the strong triangle
inequality. On genetic distances this is a strict molecular clock on a clonal
genealogy: in every triple, the two largest distances agree. -/
def IsUltrametric (d : α → α → ℝ) : Prop :=
  ∀ x y z, d x z ≤ max (d x y) (d y z)

/-- The **reticulation margin** of an ordered triple: by how much the strong
triangle inequality fails. Positive exactly when the triple is a violation. -/
noncomputable def margin (d : α → α → ℝ) (x y z : α) : ℝ :=
  d x z - max (d x y) (d y z)

/-- The statistic is exact and parameter-free: the margin is nowhere positive
precisely when the data satisfies the null. Nothing is fitted, and there is no
scale, threshold or ensemble in the statement. -/
theorem ultrametric_iff_no_positive_margin (d : α → α → ℝ) :
    IsUltrametric d ↔ ∀ x y z, margin d x y z ≤ 0 := by
  unfold IsUltrametric margin
  constructor
  · intro h x y z
    have := h x y z
    linarith
  · intro h x y z
    have := h x y z
    linarith

/-! ## The no-go: an ultrametric has no unfilled cycles, at any scale -/

/-- At every scale `ε`, an ultrametric admits no *induced path on three
vertices*: if `x` is within `ε` of `y` and `y` within `ε` of `z`, then `x` is
within `ε` of `z`.

This is the whole no-go. The relation `d · · ≤ ε` is therefore transitive, and
being reflexive and symmetric for a genetic distance it is an equivalence
relation: the Vietoris-Rips graph at scale `ε` is a disjoint union of complete
graphs. Every cycle has all of its chords, so no cycle is unfilled and
`β₁ = 0` at every scale simultaneously. -/
theorem ultrametric_no_induced_path (d : α → α → ℝ) (hu : IsUltrametric d)
    (ε : ℝ) (x y z : α) (hxy : d x y ≤ ε) (hyz : d y z ≤ ε) :
    d x z ≤ ε :=
  le_trans (hu x y z) (max_le hxy hyz)

/-- Transitivity of the Rips relation, stated as such: this is the property that
makes the components cliques, and hence the reason no cycle can be unfilled. -/
theorem rips_relation_transitive (d : α → α → ℝ) (hu : IsUltrametric d) (ε : ℝ) :
    ∀ x y z : α, d x y ≤ ε → d y z ≤ ε → d x z ≤ ε :=
  fun _ _ _ hxy hyz => ultrametric_no_induced_path d hu ε _ _ _ hxy hyz

/-! ## The certificate -/

/-- A **reticulation certificate**: three points and a scale at which two pairs
are close but the third is not. This is precisely what the Rust search emits,
and it is small enough to recheck by hand.

Carrying the three inequalities as fields means an ill-formed certificate cannot
be constructed, so a bug in the search can fail to produce a certificate but
cannot produce one that is wrong. -/
structure Certificate (d : α → α → ℝ) where
  /-- One endpoint of the open path. -/
  x : α
  /-- The midpoint: close to both endpoints. -/
  y : α
  /-- The other endpoint. -/
  z : α
  /-- The scale at which the path is open. -/
  scale : ℝ
  /-- `x` and `y` are within `scale`. -/
  close_xy : d x y ≤ scale
  /-- `y` and `z` are within `scale`. -/
  close_yz : d y z ≤ scale
  /-- …but `x` and `z` are not: the path is not filled. -/
  far_xz : scale < d x z

/-- **The inference rule.** A certificate refutes the null, and exhibits a
strictly positive margin quantifying the failure.

This is the theorem the whole design rests on. The conclusion about the data
depends on this proof and on the three numbers in the certificate — not on the
correctness of the program that found it. -/
theorem certificate_refutes_ultrametric (d : α → α → ℝ) (c : Certificate d) :
    ¬ IsUltrametric d ∧ 0 < margin d c.x c.y c.z := by
  have hmax : max (d c.x c.y) (d c.y c.z) ≤ c.scale :=
    max_le c.close_xy c.close_yz
  refine ⟨fun hu => ?_, ?_⟩
  · exact absurd
      (ultrametric_no_induced_path d hu c.scale c.x c.y c.z c.close_xy c.close_yz)
      (not_le.mpr c.far_xz)
  · unfold margin
    linarith [c.far_xz]

/-- The margin is a lower bound on how far the data is from the null: it is at
least the gap between the scale and the far distance. Used by the Rust side to
rank certificates, so the ranking is meaningful rather than conventional. -/
theorem margin_ge_gap (d : α → α → ℝ) (c : Certificate d) :
    c.scale - max (d c.x c.y) (d c.y c.z) + (d c.x c.z - c.scale)
      ≤ margin d c.x c.y c.z := by
  unfold margin
  linarith

/-! ## Non-vacuity and controls -/

/-- The minimal open path, realised concretely on the line: two adjacent pairs
within scale `1`, endpoints at distance `2`. This is the shape a recombinant and
its two parents make — each parent close to the recombinant, far from each
other — which is why such a configuration is the elementary signal. -/
noncomputable def openPath : ℝ → ℝ → ℝ := fun x y => |x - y|

/-- **Non-vacuity.** A certificate exists, so
`certificate_refutes_ultrametric` is not a statement about an empty class. A
conditional whose antecedent nothing satisfies proves nothing, and this project
has shipped that mistake before. -/
noncomputable def certificate_exists_nontrivially : Certificate openPath where
  x := 0
  y := 1
  z := 2
  scale := 1
  close_xy := by norm_num [openPath]
  close_yz := by norm_num [openPath]
  far_xz := by norm_num [openPath]

/-- The witness really does refute the null, with margin `2 - 1 = 1`. -/
theorem witness_refutes_with_margin :
    ¬ IsUltrametric openPath ∧ 0 < margin openPath 0 1 2 :=
  certificate_refutes_ultrametric openPath certificate_exists_nontrivially

/-- **The control.** Under the null, no certificate exists at any scale. A
search that returned something on ultrametric data would be reporting a property
of itself rather than of the data, which is the defect this project keeps
finding; this theorem is what rules it out. -/
theorem ultrametric_admits_no_certificate (d : α → α → ℝ) (hu : IsUltrametric d) :
    IsEmpty (Certificate d) := by
  constructor
  intro c
  exact absurd
    (ultrametric_no_induced_path d hu c.scale c.x c.y c.z c.close_xy c.close_yz)
    (not_le.mpr c.far_xz)

/-- A constant dissimilarity is an ultrametric, so the null class is itself
non-empty: the control above is not vacuous either. -/
theorem constant_is_ultrametric (r : ℝ) :
    IsUltrametric (fun _ _ : α => r) := by
  intro x y z
  simp

/-- The null has content: being a metric is not enough. The absolute-value
metric on the line satisfies the triangle inequality but not the strong one, so
`IsUltrametric` is strictly stronger and the no-go is not vacuously about all
distances. -/
theorem abs_not_ultrametric : ¬ IsUltrametric (fun x y : ℝ => |x - y|) := by
  intro h
  have := h 0 1 2
  norm_num at this

/-! ## The second, weaker null: tree metrics

An ultrametric is a clonal genealogy under a strict clock. Dropping the clock
but keeping the tree gives the four-point condition. Stated here because the
Rust side computes both, which separates *rate variation* (ultrametric fails,
four-point holds) from *reticulation* (four-point fails). Only the ultrametric
no-go is proved in this file; the corresponding statement for tree metrics is
cited in the textbook, not proved here. -/

/-- The **four-point condition**: `d` is realisable as path lengths in a tree.

Defined here because the Rust side computes its defect alongside the ultrametric
one, and the two together separate *rate variation* from *reticulation*: a
clonal genealogy with unequal branch rates fails `IsUltrametric` but satisfies
`IsFourPoint`, while recombination fails both.

Only the ultrametric no-go is proved in this file. That tree metrics also have
acyclic Vietoris-Rips filtrations is a result we cite rather than prove --- it
needs Gromov hyperbolicity machinery well beyond this file --- so the textbook
records it as cited, and the Rust side labels four-point defects as *suggestive*
where it labels ultrametric certificates as *proved*. The asymmetry is
deliberate and is the point of keeping the two statistics apart. -/
def IsFourPoint (d : α → α → ℝ) : Prop :=
  ∀ x y z w, d x y + d z w ≤ max (d x z + d y w) (d x w + d y z)

/-- The four-point defect of a quadruple: the amount by which the condition
fails. The Rust side maximises this, and reports it as a separate, weaker
statistic from the certified ultrametric margin. -/
noncomputable def fourPointDefect (d : α → α → ℝ) (x y z w : α) : ℝ :=
  d x y + d z w - max (d x z + d y w) (d x w + d y z)

/-- As with the ultrametric statistic, the four-point defect is exact and
parameter-free: it is nowhere positive precisely when the condition holds. This
much *is* proved, and it is all the Rust side needs in order to compute the
statistic honestly; what is not proved here is the topological consequence. -/
theorem four_point_iff_no_positive_defect (d : α → α → ℝ) :
    IsFourPoint d ↔ ∀ x y z w, fourPointDefect d x y z w ≤ 0 := by
  unfold IsFourPoint fourPointDefect
  constructor
  · intro h x y z w
    have := h x y z w
    linarith
  · intro h x y z w
    have := h x y z w
    linarith

end ANSE.ReticulationCertificate

#print axioms ANSE.ReticulationCertificate.ultrametric_iff_no_positive_margin
#print axioms ANSE.ReticulationCertificate.ultrametric_no_induced_path
#print axioms ANSE.ReticulationCertificate.rips_relation_transitive
#print axioms ANSE.ReticulationCertificate.certificate_refutes_ultrametric
#print axioms ANSE.ReticulationCertificate.margin_ge_gap
#print axioms ANSE.ReticulationCertificate.certificate_exists_nontrivially
#print axioms ANSE.ReticulationCertificate.witness_refutes_with_margin
#print axioms ANSE.ReticulationCertificate.ultrametric_admits_no_certificate
#print axioms ANSE.ReticulationCertificate.constant_is_ultrametric
#print axioms ANSE.ReticulationCertificate.abs_not_ultrametric
#print axioms ANSE.ReticulationCertificate.four_point_iff_no_positive_defect

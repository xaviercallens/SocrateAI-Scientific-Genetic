import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Analysis.Calculus.Deriv.Add
import Mathlib.Tactic

/-!
# One-locus, two-allele mean-fitness monotonicity (classical; often called Fisher's theorem for this model) (H12, autoresearch 2026-10-25)

First formalization attempt of a KNOWN classical result (Fisher 1930; e.g. Hartl & Clark, Ewens §2):
discrete generations, viability selection, random mating (Hardy–Weinberg proportions at zygote stage).

Genotype fitnesses `w11 w12 w22 > 0`, allele frequency `p ∈ [0,1]`, `q = 1 - p`,
`wbar p = p² w11 + 2pq w12 + q² w22`, `p' = (p² w11 + p q w12) / wbar p`.
-/

namespace AutoresearchPopGen

/-- Mean fitness under Hardy–Weinberg proportions. -/
noncomputable def wbar (w11 w12 w22 p : ℝ) : ℝ :=
  p ^ 2 * w11 + 2 * p * (1 - p) * w12 + (1 - p) ^ 2 * w22

/-- Allele frequency after one generation of viability selection. -/
noncomputable def pnext (w11 w12 w22 p : ℝ) : ℝ :=
  (p ^ 2 * w11 + p * (1 - p) * w12) / wbar w11 w12 w22 p

theorem wbar_pos {w11 w12 w22 p : ℝ} (h11 : 0 < w11) (h12 : 0 < w12) (h22 : 0 < w22)
    (hp0 : 0 ≤ p) (hp1 : p ≤ 1) : 0 < wbar w11 w12 w22 p := by
  unfold wbar
  rcases eq_or_lt_of_le hp0 with h | h
  · subst h; norm_num; exact h22
  · have hq : 0 ≤ 1 - p := by linarith
    have t1 : 0 < p ^ 2 * w11 := by positivity
    have t2 : 0 ≤ 2 * p * (1 - p) * w12 := by positivity
    have t3 : 0 ≤ (1 - p) ^ 2 * w22 := by positivity
    linarith

/-- (i) The next-generation allele frequency stays in `[0,1]`. -/
theorem pnext_mem_Icc {w11 w12 w22 p : ℝ} (h11 : 0 < w11) (h12 : 0 < w12) (h22 : 0 < w22)
    (hp0 : 0 ≤ p) (hp1 : p ≤ 1) : pnext w11 w12 w22 p ∈ Set.Icc (0 : ℝ) 1 := by
  have hw := wbar_pos h11 h12 h22 hp0 hp1
  have hq : 0 ≤ 1 - p := by linarith
  unfold pnext
  constructor
  · apply div_nonneg _ hw.le
    positivity
  · rw [div_le_one hw]
    unfold wbar
    have : 0 ≤ p * (1 - p) * w12 := by positivity
    have : 0 ≤ (1 - p) ^ 2 * w22 := by positivity
    nlinarith

/-- (ii, first form) `Δp = p q (p (w11 - w12) + q (w12 - w22)) / wbar`. -/
theorem delta_p_eq {w11 w12 w22 p : ℝ} (h11 : 0 < w11) (h12 : 0 < w12) (h22 : 0 < w22)
    (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    pnext w11 w12 w22 p - p
      = p * (1 - p) * (p * (w11 - w12) + (1 - p) * (w12 - w22)) / wbar w11 w12 w22 p := by
  have hw := (wbar_pos h11 h12 h22 hp0 hp1).ne'
  unfold pnext
  field_simp
  unfold wbar
  ring

/-- The derivative of mean fitness in `p` (as a real function). -/
theorem hasDerivAt_wbar (w11 w12 w22 p : ℝ) :
    HasDerivAt (fun x => wbar w11 w12 w22 x)
      (2 * (p * (w11 - w12) + (1 - p) * (w12 - w22))) p := by
  have h : HasDerivAt (fun x : ℝ => w22 + 2 * (w12 - w22) * x + (w11 - 2 * w12 + w22) * x ^ 2)
      (2 * (w12 - w22) * 1 + (w11 - 2 * w12 + w22) * (↑(2:ℕ) * p ^ (2 - 1))) p := by
    exact (((hasDerivAt_id p).const_mul (2 * (w12 - w22))).const_add w22).add
      ((hasDerivAt_pow 2 p).const_mul (w11 - 2 * w12 + w22))
  have e : (fun x : ℝ => wbar w11 w12 w22 x)
      = (fun x : ℝ => w22 + 2 * (w12 - w22) * x + (w11 - 2 * w12 + w22) * x ^ 2) := by
    funext x; unfold wbar; ring
  rw [e]
  convert h using 1
  push_cast
  ring

theorem deriv_wbar (w11 w12 w22 p : ℝ) :
    deriv (fun x => wbar w11 w12 w22 x) p
      = 2 * (p * (w11 - w12) + (1 - p) * (w12 - w22)) :=
  (hasDerivAt_wbar w11 w12 w22 p).deriv

/-- (ii) Wright's form: `Δp = p q wbar'(p) / (2 wbar)`. -/
theorem delta_p_eq_deriv {w11 w12 w22 p : ℝ} (h11 : 0 < w11) (h12 : 0 < w12) (h22 : 0 < w22)
    (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    pnext w11 w12 w22 p - p
      = p * (1 - p) * deriv (fun x => wbar w11 w12 w22 x) p / (2 * wbar w11 w12 w22 p) := by
  have hw := (wbar_pos h11 h12 h22 hp0 hp1).ne'
  rw [delta_p_eq h11 h12 h22 hp0 hp1, deriv_wbar]
  field_simp

/-- Exact increment identity (derived with sympy, then kernel-checked):
`wbar p' - wbar p = p q s² (wbar p + p w11 + q w22) / wbar p²`, `s = p(w11-w12) + q(w12-w22)`. -/
theorem wbar_increment {w11 w12 w22 p : ℝ} (h11 : 0 < w11) (h12 : 0 < w12) (h22 : 0 < w22)
    (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    wbar w11 w12 w22 (pnext w11 w12 w22 p) - wbar w11 w12 w22 p
      = p * (1 - p) * (p * (w11 - w12) + (1 - p) * (w12 - w22)) ^ 2
          * (wbar w11 w12 w22 p + p * w11 + (1 - p) * w22) / (wbar w11 w12 w22 p) ^ 2 := by
  have hw := (wbar_pos h11 h12 h22 hp0 hp1).ne'
  unfold pnext
  rw [eq_div_iff (pow_ne_zero 2 hw)]
  set W := wbar w11 w12 w22 p with hW
  unfold wbar
  field_simp
  rw [hW]; unfold wbar
  ring

/-- (iii) One-locus, two-allele mean-fitness monotonicity: mean fitness does not decrease. -/
theorem fisher_mean_fitness_nondecreasing {w11 w12 w22 p : ℝ}
    (h11 : 0 < w11) (h12 : 0 < w12) (h22 : 0 < w22) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    wbar w11 w12 w22 p ≤ wbar w11 w12 w22 (pnext w11 w12 w22 p) := by
  have hw := wbar_pos h11 h12 h22 hp0 hp1
  have hinc := wbar_increment h11 h12 h22 hp0 hp1
  have hq : 0 ≤ 1 - p := by linarith
  have : 0 ≤ p * (1 - p) * (p * (w11 - w12) + (1 - p) * (w12 - w22)) ^ 2
          * (wbar w11 w12 w22 p + p * w11 + (1 - p) * w22) / (wbar w11 w12 w22 p) ^ 2 := by
    positivity
  linarith

/-! ## Hardy–Weinberg -/

/-- One round of random mating from genotype frequencies `(x11, x12, x22)`:
allele frequency `p = x11 + x12/2`, offspring in proportions `(p², 2pq, q²)`. -/
noncomputable def hwStep (x : ℝ × ℝ × ℝ) : ℝ × ℝ × ℝ :=
  let p := x.1 + x.2.1 / 2
  (p ^ 2, 2 * p * (1 - p), (1 - p) ^ 2)

theorem hwStep_sum (x : ℝ × ℝ × ℝ) :
    (hwStep x).1 + (hwStep x).2.1 + (hwStep x).2.2 = 1 := by
  simp only [hwStep]; ring

/-- Allele frequencies are preserved by random mating: the offspring A1 frequency equals the
parental `x11 + x12/2`, and (when parental genotype frequencies sum to 1) the offspring A2 frequency
equals the parental `x22 + x12/2`. -/
theorem hwStep_allele {x11 x12 x22 : ℝ} (hs : x11 + x12 + x22 = 1) :
    (hwStep (x11, x12, x22)).1 + (hwStep (x11, x12, x22)).2.1 / 2 = x11 + x12 / 2 ∧
    (hwStep (x11, x12, x22)).2.2 + (hwStep (x11, x12, x22)).2.1 / 2 = x22 + x12 / 2 := by
  simp only [hwStep]
  constructor
  · ring
  · linear_combination (-1 : ℝ) * hs

/-- Hardy–Weinberg proportions are reached after one round and are then stationary. -/
theorem hwStep_stationary (x : ℝ × ℝ × ℝ) : hwStep (hwStep x) = hwStep x := by
  obtain ⟨x11, x12, x22⟩ := x
  simp only [hwStep, Prod.mk.injEq]
  refine ⟨by ring, by ring, by ring⟩

/-! ## Negative controls (informative) -/

/-- Underdominance (w12 < w11, w22): w11 = w22 = 1, w12 = 1/2 has an interior equilibrium at p = 1/2,
yet from p = 2/5 the allele frequency moves AWAY from it (Δp < 0) while mean fitness still rises.
Refutes the naive reading "selection pushes p towards the polymorphic equilibrium". -/
theorem underdominance_moves_away :
    pnext 1 (1/2) 1 (2/5) - 2/5 < 0 ∧ pnext 1 (1/2) 1 (1/2) = 1/2 ∧
      wbar 1 (1/2) 1 (2/5) < wbar 1 (1/2) 1 (pnext 1 (1/2) 1 (2/5)) := by
  refine ⟨?_, ?_, ?_⟩ <;> norm_num [pnext, wbar]

/-- Positivity hypotheses are needed: with a negative homozygote fitness, w11 = -1, w12 = 1,
w22 = 0, p = 1/2 (still wbar = 1/4 > 0), one generation gives p' = 0 and mean fitness DROPS to 0.
By `wbar_increment` the sign is that of `wbar + p w11 + q w22`, which positivity of the fitnesses
guarantees; this instance shows the conclusion is false once that hypothesis is dropped. -/
theorem fisher_fails_without_positivity :
    0 < wbar (-1) 1 0 (1/2) ∧
      wbar (-1) 1 0 (pnext (-1) 1 0 (1/2)) < wbar (-1) 1 0 (1/2) := by
  refine ⟨?_, ?_⟩ <;> norm_num [pnext, wbar]

/-! ## Planted controls (must be caught) -/

end AutoresearchPopGen

#print axioms AutoresearchPopGen.wbar_pos
#print axioms AutoresearchPopGen.pnext_mem_Icc
#print axioms AutoresearchPopGen.delta_p_eq
#print axioms AutoresearchPopGen.hasDerivAt_wbar
#print axioms AutoresearchPopGen.deriv_wbar
#print axioms AutoresearchPopGen.delta_p_eq_deriv
#print axioms AutoresearchPopGen.wbar_increment
#print axioms AutoresearchPopGen.fisher_mean_fitness_nondecreasing
#print axioms AutoresearchPopGen.hwStep_sum
#print axioms AutoresearchPopGen.hwStep_allele
#print axioms AutoresearchPopGen.hwStep_stationary
#print axioms AutoresearchPopGen.underdominance_moves_away
#print axioms AutoresearchPopGen.fisher_fails_without_positivity

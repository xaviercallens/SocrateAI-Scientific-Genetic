import Mathlib.Tactic

/-  Keystone feasibility check for the proposed architecture.

    A p-distance is (differing sites) / (comparable sites) -- a RATIONAL number
    with an integer numerator and denominator. So a certificate can be stated
    over ℚ with integer literals and decided by the kernel exactly, with no
    floating point anywhere in the trust path.

    Below: can Lean check a certificate emitted as integer counts? -/

structure QCertificate where
  d_xy_num : ℤ
  d_yz_num : ℤ
  d_xz_num : ℤ
  sites    : ℤ
  sites_pos : 0 < sites
  -- scale = max of the two near distances, as a rational
  close_xy : (d_xy_num : ℚ) / sites ≤ (max d_xy_num d_yz_num : ℚ) / sites
  close_yz : (d_yz_num : ℚ) / sites ≤ (max d_xy_num d_yz_num : ℚ) / sites
  far_xz   : (max d_xy_num d_yz_num : ℚ) / sites < (d_xz_num : ℚ) / sites

/-- A certificate the Rust side would emit for the simulated recombinant:
    4000 comparable sites, d_xy = 69/4000, d_yz = 70/4000, d_xz = 131/4000. -/
def emitted : QCertificate where
  d_xy_num := 69
  d_yz_num := 70
  d_xz_num := 131
  sites := 4000
  sites_pos := by norm_num
  close_xy := by norm_num
  close_yz := by norm_num
  far_xz := by norm_num

/-- And the margin is exactly rational, not approximately: 61/4000. -/
theorem emitted_margin :
    ((131 : ℚ) / 4000) - ((70 : ℚ) / 4000) = 61 / 4000 := by norm_num

/-- A FORGED certificate must fail to typecheck. Demonstrated by showing the
    inequality it would need is false, so no such term can be built. -/
theorem forgery_is_impossible :
    ¬ (((70 : ℚ) / 4000) < ((69 : ℚ) / 4000)) := by norm_num

#print axioms emitted
#print axioms emitted_margin
#print axioms forgery_is_impossible

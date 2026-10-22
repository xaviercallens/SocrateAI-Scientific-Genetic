/-
  A certificate ABOUT THE SEQUENCES, for the kernel to check.

  Source : simulated recombinant, seed 2026
  Triple : clonal_002 -- recombinant_of_clonal_003_and_clonal_013 -- clonal_013

  The three sequences are embedded below. Every count is proved by
  `rfl`, so the KERNEL computes the differing and comparable sites --
  they are not supplied by the search. The only remaining trusted
  step is that these lists are the dataset's sequences, which a
  reader can diff against the alignment.

  Encoding: A=0 C=1 G=2 T=3, anything else 4 (not comparable).
  Generated file -- do not edit by hand.
-/
import Mathlib.Tactic
import ANSE.GeneticaReticulata

set_option maxRecDepth 1000000

namespace Reticulata.EmittedSequences

def sx : List Nat := [3,1,2,2,1,3,2,1,0,2,1,0,0,1,1,1,1,3,3,1,0,0,0,2,1,1,3,2,0,0,2,2,2,3,1,3,2,0,2,2,2,3,1,0,1,0,1,2,2,0,0,0,3,2,1,1,2,1,1,2,3,3,3,1,1,0,1,1,1,2,0,2,2,2,1,3,2,0,3,2,1,3,2,3,1,1,3,0,2,2,2,3,1,1,0,0,0,0,0,2,2,1,3,0,1,3,0,1,0,1,2,2,0,3,0,3,0,2,0,1,2,0,3,0,3,2,1,0,2,3,2,2,3,2,3,2,1,2,0,1,0,1,1,0,3,0,3,2,0,1,0,2,1,0,2,3,3,1,1,0,1,3,1,0,1,3,2,1,2,1,3,0,3,1,1,3,3,1,2,1,2,2,0,1,1,3,2,0,0,1,2,1,3,2,1,3,0,2,0,1,3,1,3,1,2,0,1,3,2,2,1,1,0,3,2,2,2,1,0,3,2,2,2,2,0,0,2,2,0,2,1,0,0,2,0,2,1,3,2,2,1,3,2,2,2,3,0,2,2,2,3,1,1,3,0,0,2,0,2,0,2,2,3,1,1,2,0,1,2,2,1,1,2,2,3,1,1,3,0,2,3,2,0,2,0,3,2,0,0,0,2,1,0,2,1,1,0,1,3,0,0,3,2,2,2,1,2,1,0,2,1,0,2,0,0,1,2,2,3,2,1,2,0,

-- [ ... 25255 bytes total; the three embedded sequences are elided here. The full file is produced by
--   reticulation poc --formal-dir <formal> --out <dir>
-- and is what the kernel actually checks. ... ]

,2,0,1,1,1,1,2,2,3,0,1,3,1,1,2,0,3,1,3,3,3,2,0,0,3,0,1,1,0,0,2,2,1,1,2,0,2,0,1,2,0,1,0,0,2,0,1,1,3,0,2,0,2,2,1,0,3,2,2,1,0,3,3,2,2,1,2,1,3,2,0,1,2,1,0,1,1,2,0,1,2,3,2,2,1,2,2,1,0,0,2,2,0,2,0,0,0,3,3,3,1,3,3,2,3,2,1,1,2,2,3,2,0,3,0,0,2,1,3,3,0,0]

/-- The kernel computes the counts and derives the violation. -/
theorem violates :
    max (ANSE.GeneticaReticulata.pDistance sx sy)
        (ANSE.GeneticaReticulata.pDistance sy sz)
      < ANSE.GeneticaReticulata.pDistance sx sz := by
  apply ANSE.GeneticaReticulata.sequences_refute_from_counts
    sx sy sz 69 4000 70 4000 131 4000 <;>
    first | rfl | norm_num

end Reticulata.EmittedSequences

#print axioms Reticulata.EmittedSequences.violates

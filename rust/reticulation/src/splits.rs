//! Certifying a **split of the whole population**, not just a quadruple.
//!
//! The quadruple-level result (`subset::buneman_indices`) established that a
//! quartet certificate *is* the split-decomposition condition for a split system
//! to fail to be a tree. That is a statement about four taxa. What a biologist
//! reads is a *split*: this group of sequences against that group. This module
//! promotes the certificate to that level.
//!
//! # The one design decision that matters
//!
//! Bandelt and Dress define a split's isolation index as a **minimum over the
//! quadruples straddling it**. Taken literally that is unusable here, because
//! under complete deletion *within the subset* every quadruple has its own
//! denominator — the columns at which those particular four sequences all carry
//! a definite base. A minimum over quantities with different denominators is an
//! incommensurable comparison, the defect this project has caught in its own
//! apparatus more than once.
//!
//! So the split certificate is **a logical conjunction of per-quadruple integer
//! comparisons**, never a numeric minimum:
//!
//! > The split `A|B` is certified iff for **every** quadruple `{a,a',b,b'}` with
//! > `a,a' ∈ A` and `b,b' ∈ B`, the split's own pairing sum is strictly less
//! > than the larger of the two crossing sums — each comparison exact, over that
//! > quadruple's own denominator.
//!
//! Nothing needs a common metric, nothing is averaged, and every conjunct is an
//! integer comparison the Lean kernel can check. The split-level statement
//! inherits its exactness from the quadruple level rather than paying for it.
//!
//! # Why this is the right generalisation
//!
//! On an exact tree metric the two crossing sums of any quadruple are equal and
//! larger than the supported one, so the tree's own splits are certified
//! (`own = s₁ < s₂ = s₃`) and every other split is not (`own = s₂ = s₃ = max`,
//! so `own < max` is false). Off a tree metric a quadruple's three sums are
//! generically distinct and **two** pairings can be non-maximal at once — which
//! is the same two-positive-indices condition the quadruple experiment measured.
//!
//! # The theorem this module exists to exhibit
//!
//! If two *incompatible* splits are both certified, then the four-point
//! condition **must** fail, and a quadruple witnessing it can be constructed
//! rather than searched for. Incompatibility means all four intersections
//! `A∩C, A∩D, B∩C, B∩D` are non-empty; take one taxon from each. On that
//! quadruple the two splits induce two *different* pairings, both certified
//! non-maximal, so the third pairing is strictly the largest — which is exactly
//! a four-point violation. Proved in `GeneticaReticulataTree.lean`
//! (`third_pairing_strictly_largest_of_two_certified`), and
//! `witness_quadruple` below constructs the witness the proof describes.
//!
//! That is the lift: a population-level claim about two groups of sequences,
//! backed by a quadruple-level certificate that a kernel can check.

use crate::subset::Alignment;

/// A bipartition of the taxon set, held as a bitmask of one side.
///
/// Canonical form always contains taxon 0, so a split has exactly one
/// representation and two splits are equal iff their masks are.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Split {
    mask: u64,
    taxa: usize,
}

impl Split {
    /// Canonicalise a mask into the side containing taxon 0.
    pub fn new(mask: u64, taxa: usize) -> Self {
        let full = if taxa >= 64 { u64::MAX } else { (1u64 << taxa) - 1 };
        let mask = mask & full;
        let mask = if mask & 1 == 1 { mask } else { full & !mask };
        Split { mask, taxa }
    }

    /// The side containing taxon 0.
    pub fn left(&self) -> u64 {
        self.mask
    }

    /// The other side.
    pub fn right(&self) -> u64 {
        let full = if self.taxa >= 64 {
            u64::MAX
        } else {
            (1u64 << self.taxa) - 1
        };
        full & !self.mask
    }

    /// Members of each side, as indices.
    pub fn sides(&self) -> (Vec<usize>, Vec<usize>) {
        let (mut a, mut b) = (Vec::new(), Vec::new());
        for i in 0..self.taxa {
            if self.mask >> i & 1 == 1 {
                a.push(i);
            } else {
                b.push(i);
            }
        }
        (a, b)
    }

    /// The smaller side's size, which is the split's "depth".
    pub fn smaller_side(&self) -> usize {
        let l = self.left().count_ones() as usize;
        l.min(self.taxa - l)
    }

    /// Both sides have at least two taxa, so straddling quadruples exist.
    ///
    /// A split isolating a single taxon has no straddling quadruple and so is
    /// trivially "certified" by an empty conjunction. Those carry no evidence
    /// and are excluded rather than reported as free certificates.
    pub fn is_informative(&self) -> bool {
        self.smaller_side() >= 2
    }

    /// Compatibility in the split-system sense: at least one of the four
    /// intersections is empty.
    ///
    /// By the splits-equivalence theorem a set of splits is displayed by a tree
    /// exactly when it is pairwise compatible, so an incompatible pair is a
    /// population-level statement that no tree displays both.
    pub fn compatible_with(&self, other: &Split) -> bool {
        let (a, b) = (self.left(), self.right());
        let (c, d) = (other.left(), other.right());
        a & c == 0 || a & d == 0 || b & c == 0 || b & d == 0
    }

    /// One taxon from each of the four intersections, when incompatible.
    ///
    /// This is the witness the Lean theorem describes: on this quadruple the two
    /// splits induce two different pairings, so if both splits are certified the
    /// four-point condition fails here. Returned as
    /// `(a ∈ A∩C, b ∈ A∩D, c ∈ B∩C, e ∈ B∩D)`.
    pub fn witness_quadruple(&self, other: &Split) -> Option<(usize, usize, usize, usize)> {
        if self.compatible_with(other) {
            return None;
        }
        let (a, b) = (self.left(), self.right());
        let (c, d) = (other.left(), other.right());
        let first = |m: u64| -> Option<usize> {
            if m == 0 {
                None
            } else {
                Some(m.trailing_zeros() as usize)
            }
        };
        Some((
            first(a & c)?,
            first(a & d)?,
            first(b & c)?,
            first(b & d)?,
        ))
    }
}

/// Every informative split of `taxa` taxa, canonicalised and deduplicated.
pub fn enumerate_splits(taxa: usize) -> Vec<Split> {
    assert!(taxa >= 4, "a split certificate needs four taxa");
    assert!(taxa < 64, "bitmask splits are limited to 63 taxa");
    let mut seen = std::collections::BTreeSet::new();
    for mask in 1..(1u64 << taxa) {
        let split = Split::new(mask, taxa);
        if split.is_informative() {
            seen.insert(split);
        }
    }
    seen.into_iter().collect()
}

/// The verdict on one split, with the evidence that produced it.
#[derive(Debug, Clone, serde::Serialize)]
pub struct SplitVerdict {
    pub left: Vec<String>,
    pub right: Vec<String>,
    pub smaller_side: usize,
    /// Straddling quadruples examined: the size of the conjunction.
    pub conjuncts: usize,
    /// Conjuncts on which the split's own pairing was strictly non-maximal.
    pub conjuncts_satisfied: usize,
    /// Straddling quadruples with no column where all four are definite. The
    /// comparison is undefined there, so the split is not certified -- reported
    /// rather than silently skipped.
    pub conjuncts_undecidable: usize,
    /// Certified iff every conjunct is satisfied and none is undecidable.
    pub certified: bool,
}

/// A noise floor for the conjunction, as an exact rational rate.
///
/// The unthresholded criterion `own < max(cross)` has no floor, and the first
/// run of this experiment showed what that costs: on a clonal population with
/// no reticulation the two crossing sums of a quadruple differed by 2 counts in
/// 2000 -- pure sampling noise -- and the split was certified on that margin. On
/// an exact tree metric those two sums are EQUAL, so their observed difference
/// is noise of order sqrt(L) and any criterion that reads its sign will fire on
/// a tree.
///
/// So a conjunct must clear a *rate*: the margin, as a fraction of that
/// quadruple's own denominator, must exceed `num/den`. The comparison stays an
/// exact integer one -- `margin * den > num * common` -- so nothing is
/// introduced that a kernel cannot check.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub struct Floor {
    pub num: i64,
    pub den: i64,
}

impl Floor {
    /// No floor: the original criterion, kept so the failure stays reproducible.
    pub const ZERO: Floor = Floor { num: 0, den: 1 };

    pub fn new(num: i64, den: i64) -> Result<Self, String> {
        if den <= 0 {
            return Err(format!("floor denominator must be positive, got {den}"));
        }
        if num < 0 {
            return Err(format!("floor numerator must be non-negative, got {num}"));
        }
        Ok(Floor { num, den })
    }

    /// `margin / common > num / den`, as integers.
    pub fn cleared_by(&self, margin: i64, common: i64) -> bool {
        margin * self.den > self.num * common
    }
}

/// Certify one split: a conjunction of exact integer comparisons.
pub fn certify_split(aln: &Alignment, split: &Split) -> SplitVerdict {
    certify_split_at(aln, split, Floor::ZERO)
}

/// Certify one split with a noise floor on every conjunct.
pub fn certify_split_at(aln: &Alignment, split: &Split, floor: Floor) -> SplitVerdict {
    let (left, right) = split.sides();
    let mut conjuncts = 0usize;
    let mut satisfied = 0usize;
    let mut undecidable = 0usize;

    for ia in 0..left.len() {
        for ib in (ia + 1)..left.len() {
            for ja in 0..right.len() {
                for jb in (ja + 1)..right.len() {
                    conjuncts += 1;
                    let counts = aln.quadruple_counts(left[ia], left[ib], right[ja], right[jb]);
                    if counts.common <= 0 {
                        undecidable += 1;
                        continue;
                    }
                    // `pairing_sums()[0]` is (1,2)+(3,4), i.e. the split's own
                    // pairing, because the quadruple was passed in that order.
                    let sums = counts.pairing_sums();
                    let own = sums[0];
                    let crossing_max = sums[1].max(sums[2]);
                    if floor.cleared_by(crossing_max - own, counts.common) {
                        satisfied += 1;
                    }
                }
            }
        }
    }

    let label = |idx: &[usize]| -> Vec<String> {
        idx.iter().map(|&i| aln.labels[i].clone()).collect()
    };
    SplitVerdict {
        left: label(&left),
        right: label(&right),
        smaller_side: split.smaller_side(),
        conjuncts,
        conjuncts_satisfied: satisfied,
        conjuncts_undecidable: undecidable,
        certified: conjuncts > 0 && satisfied == conjuncts && undecidable == 0,
    }
}

/// An incompatible pair of certified splits, with its witness quadruple.
#[derive(Debug, Clone, serde::Serialize)]
pub struct IncompatiblePair {
    pub first: (Vec<String>, Vec<String>),
    pub second: (Vec<String>, Vec<String>),
    /// The quadruple the Lean theorem constructs, in witness order.
    pub witness: [String; 4],
    /// The witness quadruple's three pairing sums, over its own denominator.
    pub witness_common: i64,
    pub witness_pairing_sums: [i64; 3],
    /// Checked, not assumed: the four-point condition really does fail here.
    pub four_point_fails_on_witness: bool,
}

/// The whole split system of a population, and whether it is a tree.
#[derive(Debug, Clone, serde::Serialize)]
pub struct SplitSystem {
    /// The noise floor every conjunct had to clear.
    pub floor: Floor,
    pub taxa: usize,
    pub informative_splits: usize,
    pub certified: Vec<SplitVerdict>,
    pub certified_count: usize,
    /// Pairwise compatible, i.e. displayed by some tree.
    pub pairwise_compatible: bool,
    pub incompatible_pairs: Vec<IncompatiblePair>,
    /// Every incompatible pair's witness really violates four-point. This is a
    /// check on the theorem, not a restatement of it.
    pub every_witness_violates_four_point: bool,
}

/// Enumerate, certify, and test the certified system for tree-likeness.
pub fn split_system(aln: &Alignment) -> Result<SplitSystem, String> {
    split_system_at(aln, Floor::ZERO)
}

/// The same, with a noise floor on every conjunct.
pub fn split_system_at(aln: &Alignment, floor: Floor) -> Result<SplitSystem, String> {
    let taxa = aln.len();
    if taxa < 4 {
        return Err(format!("need four taxa, got {taxa}"));
    }
    if taxa >= 64 {
        return Err(format!("bitmask splits are limited to 63 taxa, got {taxa}"));
    }
    let splits = enumerate_splits(taxa);
    let mut certified_splits = Vec::new();
    let mut verdicts = Vec::new();
    for split in &splits {
        let verdict = certify_split_at(aln, split, floor);
        if verdict.certified {
            certified_splits.push(*split);
            verdicts.push(verdict);
        }
    }

    let mut pairs = Vec::new();
    let mut all_witnesses_violate = true;
    for i in 0..certified_splits.len() {
        for j in (i + 1)..certified_splits.len() {
            let (first, second) = (certified_splits[i], certified_splits[j]);
            if first.compatible_with(&second) {
                continue;
            }
            let (a, b, c, e) = first
                .witness_quadruple(&second)
                .ok_or("incompatible split pair without a witness")?;
            let counts = aln.quadruple_counts(a, b, c, e);
            let sums = counts.pairing_sums();
            // Four-point fails iff some pairing strictly exceeds the other two.
            let largest = sums[0].max(sums[1]).max(sums[2]);
            let strictly_largest = sums.iter().filter(|&&s| s == largest).count() == 1;
            let fails = counts.common > 0 && strictly_largest;
            if !fails {
                all_witnesses_violate = false;
            }
            let name = |i: usize| aln.labels[i].clone();
            let (fl, fr) = first.sides();
            let (sl, sr) = second.sides();
            let label = |idx: &[usize]| -> Vec<String> {
                idx.iter().map(|&i| aln.labels[i].clone()).collect()
            };
            pairs.push(IncompatiblePair {
                first: (label(&fl), label(&fr)),
                second: (label(&sl), label(&sr)),
                witness: [name(a), name(b), name(c), name(e)],
                witness_common: counts.common,
                witness_pairing_sums: sums,
                four_point_fails_on_witness: fails,
            });
        }
    }

    Ok(SplitSystem {
        floor,
        taxa,
        informative_splits: splits.len(),
        certified_count: verdicts.len(),
        certified: verdicts,
        pairwise_compatible: pairs.is_empty(),
        every_witness_violates_four_point: all_witnesses_violate,
        incompatible_pairs: pairs,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln(spec: &[(&str, &str)]) -> Alignment {
        Alignment::new(
            spec.iter().map(|(n, _)| n.to_string()).collect(),
            spec.iter().map(|(_, s)| s.as_bytes().to_vec()).collect(),
        )
        .expect("well-formed alignment")
    }

    /// Six sites realising the tree `((a,b),(c,d))` with every pendant edge of
    /// length one and the internal edge of length two. One site per unit of
    /// branch length, so the pairing sums are `[4, 8, 8]` by construction: the
    /// two crossing sums are EQUAL, which is the signature of a tree metric.
    fn tree_quadruple() -> Alignment {
        aln(&[
            ("a", "CAAAAA"),
            ("b", "ACAAAA"),
            ("c", "AACACC"),
            ("d", "AAACCC"),
        ])
    }

    /// Four sites of `AACC` plus two of `ACAC`, giving pairing sums
    /// `[4, 8, 12]` -- all distinct, with a unique maximum, so the four-point
    /// condition fails.
    fn non_tree_quadruple() -> Alignment {
        aln(&[
            ("a", "AAAAAA"),
            ("b", "AAAACC"),
            ("c", "CCCCAA"),
            ("d", "CCCCCC"),
        ])
    }

    #[test]
    fn a_four_set_has_exactly_three_informative_splits() {
        let splits = enumerate_splits(4);
        assert_eq!(splits.len(), 3, "the three balanced bipartitions");
        assert!(
            splits.iter().all(|s| s.smaller_side() == 2),
            "each has two taxa a side"
        );
    }

    #[test]
    fn the_hand_built_sums_are_what_the_construction_intended() {
        let sums = tree_quadruple().quadruple_counts(0, 1, 2, 3).pairing_sums();
        assert_eq!(sums, [4, 8, 8], "tree metric: crossing sums equal");
        let sums = non_tree_quadruple()
            .quadruple_counts(0, 1, 2, 3)
            .pairing_sums();
        assert_eq!(sums, [4, 8, 12], "non-tree: a unique maximum");
    }

    #[test]
    fn nested_splits_are_compatible_and_crossing_splits_are_not() {
        // {0,1} | {2,3,4,5} against {0,1,2} | {3,4,5}: nested, so compatible.
        let nested_a = Split::new(0b000011, 6);
        let nested_b = Split::new(0b000111, 6);
        assert!(
            nested_a.compatible_with(&nested_b),
            "one side of each is contained in a side of the other"
        );
        // {0,1} | {2,3,4,5} against {0,2} | {1,3,4,5}: crossing.
        let crossing = Split::new(0b000101, 6);
        assert!(
            !nested_a.compatible_with(&crossing),
            "all four intersections are non-empty"
        );
        assert!(
            nested_a.compatible_with(&nested_a),
            "a split never crosses itself"
        );
    }

    #[test]
    fn a_tree_metric_certifies_exactly_the_tree_split() {
        let system = split_system(&tree_quadruple()).expect("four taxa");
        assert_eq!(
            system.certified_count, 1,
            "only the supported split clears the comparison"
        );
        let certified = &system.certified[0];
        let mut left = certified.left.clone();
        left.sort();
        assert_eq!(left, vec!["a".to_string(), "b".to_string()]);
        assert!(system.pairwise_compatible, "one split is trivially a tree");
        assert!(system.incompatible_pairs.is_empty());
    }

    #[test]
    fn a_four_point_violation_certifies_two_incompatible_splits() {
        let system = split_system(&non_tree_quadruple()).expect("four taxa");
        assert_eq!(
            system.certified_count, 2,
            "the smallest and middle sums both clear the comparison"
        );
        assert!(
            !system.pairwise_compatible,
            "no tree displays two splits of a four-set"
        );
        assert_eq!(system.incompatible_pairs.len(), 1);
        let pair = &system.incompatible_pairs[0];
        assert!(
            pair.four_point_fails_on_witness,
            "the constructed witness really violates four-point"
        );
        assert_eq!(pair.witness_pairing_sums, [4, 8, 12]);
    }

    #[test]
    fn the_witness_takes_one_taxon_from_each_intersection() {
        let first = Split::new(0b000011, 6); // {0,1} | {2,3,4,5}
        let second = Split::new(0b000101, 6); // {0,2} | {1,3,4,5}
        let (a, b, c, e) = first
            .witness_quadruple(&second)
            .expect("the pair is incompatible");
        assert_eq!([a, b, c, e], [0, 1, 2, 3], "one from each intersection");
        let mut distinct = vec![a, b, c, e];
        distinct.sort_unstable();
        distinct.dedup();
        assert_eq!(distinct.len(), 4, "four distinct taxa");
    }

    #[test]
    fn a_singleton_split_carries_no_evidence_and_is_excluded() {
        // {0} | {1,2,3} has no straddling quadruple, so an unguarded
        // conjunction over its quadruples would be vacuously true.
        let singleton = Split::new(0b0001, 4);
        assert!(
            !singleton.is_informative(),
            "a split isolating one taxon is excluded"
        );
        let verdict = certify_split(&non_tree_quadruple(), &singleton);
        assert_eq!(verdict.conjuncts, 0, "no straddling quadruples exist");
        assert!(
            !verdict.certified,
            "an empty conjunction must not count as a certificate"
        );
    }
}

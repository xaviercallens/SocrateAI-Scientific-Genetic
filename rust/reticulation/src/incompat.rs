//! The incumbent this method has to be compared against: pairwise homoplasy.
//!
//! # Why this module exists
//!
//! Paper F lists, among its own limitations, that no head-to-head against the
//! established tests had been run, and that until one had, the four-point
//! certificate's standing relative to the incumbents was unestablished. This
//! module removes that limitation for the one incumbent that matters most to
//! the comparison.
//!
//! The relevant incumbent is **not** one of the many recombination tests that
//! assume a clock. It is the pairwise homoplasy index of Bruen, Philippe and
//! Bryant (*Genetics* 172:2665–2681, 2006), whose null is **site-order
//! exchangeability** — no clock, no tree, no rate model. That makes it the
//! natural competitor: like the four-point condition it is immune, by
//! construction, to the relaxed-clock confound that destroyed the triple
//! statistic, so a comparison between them is about power rather than about who
//! survives the obvious control.
//!
//! # The statistic
//!
//! Two biallelic sites are **incompatible** when all four combinations of their
//! two states occur among the sequences — the four-gamete test. Under a single
//! tree, compatibility between two sites does not depend on how far apart they
//! sit in the alignment. Under recombination it does: neighbouring sites tend
//! to share a local history and so to be compatible, while distant sites may
//! come from different histories. The statistic is therefore the mean
//! compatibility of *neighbouring* informative sites, and the null is obtained
//! by permuting the site order.
//!
//! Two deliberate restrictions, both documented rather than hidden:
//!
//! * Only **biallelic parsimony-informative** columns are used — exactly two
//!   distinct bases, each carried by at least two sequences. The four-gamete
//!   test is then exact and needs no refinement for multistate characters.
//!   Bruen et al. define a refined incompatibility score for the general case;
//!   this is the restricted version, and it is named as such rather than
//!   presented as their statistic in full.
//! * Columns with a gap in any member of the subset are skipped, which is
//!   complete deletion within the subset — the same policy the rest of this
//!   crate uses, so the two statistics see the same data.
//!
//! # Exactness
//!
//! Φ_w is a count of compatible pairs over a count of pairs, so it is an exact
//! rational; and the pair count is invariant under permutation of the site
//! order, so the permutation test compares *numerators* and is exact integer
//! arithmetic throughout. No floating point decides anything here. That is not
//! a property of the original test, which is usually reported as a normal
//! approximation; it is a property of this implementation, and it means the
//! head-to-head is decided the same way the certificate is.

use serde::{Deserialize, Serialize};

use crate::exact::Ratio;
use crate::simulate::SplitMix64;
use crate::subset::Alignment;

fn is_base(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

/// One biallelic informative column, reduced to a bit per sequence.
///
/// `mask` has bit `i` set when sequence `i` carries the second of the column's
/// two states. With at most 64 sequences this is a single word, and the
/// four-gamete test becomes four bitwise tests.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BiallelicSite {
    /// Column index in the alignment.
    pub column: usize,
    /// One bit per sequence: set for the second state.
    pub mask: u64,
}

/// The outcome of a pairwise-homoplasy test.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhiTest {
    /// Biallelic parsimony-informative columns found.
    pub informative_sites: usize,
    /// Neighbouring pairs scored, i.e. pairs within the window.
    pub pairs_scored: usize,
    /// Compatible pairs among those: the statistic's numerator.
    pub compatible_pairs: usize,
    /// Φ_w as an exact rational: compatible pairs over pairs scored.
    pub phi_w_num: i64,
    /// Its denominator.
    pub phi_w_den: i64,
    /// Φ_w as a decimal, for reading only.
    pub phi_w: f64,
    /// Permutations of the site order drawn for the null.
    pub permutations: usize,
    /// Permutations whose statistic was at least the observed one.
    pub permutations_at_least_observed: usize,
    /// `(permutations_at_least_observed + 1) / (permutations + 1)`, the standard
    /// conservative permutation p-value. Exact: both parts are counts.
    pub p_value: f64,
    /// The smallest p-value this many permutations can produce. A p-value equal
    /// to the floor means "as extreme as this test can report", not "zero".
    pub p_value_floor: f64,
    /// Mean Φ_w over the permutations, for reading: the level the statistic sits
    /// at when site order carries no information.
    pub permuted_mean_phi_w: f64,
}

impl PhiTest {
    /// Φ_w as an exact ratio.
    pub fn phi_w_exact(&self) -> Ratio {
        Ratio {
            num: self.phi_w_num,
            den: self.phi_w_den,
        }
    }
}

/// Extract the biallelic parsimony-informative columns of an alignment.
///
/// Complete deletion within the subset: a column with a gap or an ambiguity in
/// any sequence is skipped entirely.
pub fn informative_sites(aln: &Alignment) -> Result<Vec<BiallelicSite>, String> {
    let n = aln.len();
    if n > 64 {
        return Err(format!(
            "this implementation packs one bit per sequence into a u64, so it \
             handles at most 64 sequences; got {n}. Widen the mask deliberately \
             rather than silently truncating the population."
        ));
    }
    let width = aln.width();
    let mut sites = Vec::new();
    for column in 0..width {
        let mut states: Vec<(u8, usize)> = Vec::new();
        let mut usable = true;
        for row in &aln.sequences {
            let base = row[column];
            if !is_base(base) {
                usable = false;
                break;
            }
            let base = base.to_ascii_uppercase();
            match states.iter_mut().find(|(b, _)| *b == base) {
                Some((_, count)) => *count += 1,
                None => states.push((base, 1)),
            }
        }
        if !usable || states.len() != 2 {
            continue;
        }
        // Parsimony-informative: both states carried at least twice. A singleton
        // state is uninformative about any split and would inflate the
        // compatible-pair count.
        if states.iter().any(|(_, count)| *count < 2) {
            continue;
        }
        let second = states[1].0;
        let mut mask = 0u64;
        for (index, row) in aln.sequences.iter().enumerate() {
            if row[column].to_ascii_uppercase() == second {
                mask |= 1u64 << index;
            }
        }
        sites.push(BiallelicSite { column, mask });
    }
    Ok(sites)
}

/// The four-gamete test: are these two biallelic sites compatible?
///
/// Incompatible exactly when all four state combinations occur, which on the bit
/// representation is four non-empty intersections.
pub fn compatible(left: u64, right: u64, present: u64) -> bool {
    let both = left & right;
    let only_left = left & !right & present;
    let only_right = right & !left & present;
    let neither = !left & !right & present;
    !(both != 0 && only_left != 0 && only_right != 0 && neither != 0)
}

fn phi_numerator(sites: &[BiallelicSite], order: &[usize], window: usize, present: u64) -> (usize, usize) {
    let mut compatible_pairs = 0usize;
    let mut pairs = 0usize;
    let count = order.len();
    for i in 0..count {
        let upper = (i + window).min(count - 1);
        for j in (i + 1)..=upper {
            pairs += 1;
            if compatible(sites[order[i]].mask, sites[order[j]].mask, present) {
                compatible_pairs += 1;
            }
        }
    }
    (compatible_pairs, pairs)
}

/// Run the restricted pairwise-homoplasy test.
///
/// `window` is the neighbourhood in informative-site index space (Bruen et al.
/// use 100 by default). `permutations` draws of the site order give the null.
///
/// # Errors
///
/// More than 64 sequences, or fewer than `window + 2` informative sites, in
/// which case the window covers the whole alignment and the statistic cannot
/// distinguish neighbouring pairs from distant ones — reported rather than
/// returned as a number that means nothing.
pub fn phi_test(
    aln: &Alignment,
    window: usize,
    permutations: usize,
    seed: u64,
) -> Result<PhiTest, String> {
    let sites = informative_sites(aln)?;
    if sites.len() < window + 2 {
        return Err(format!(
            "{} biallelic informative sites for a window of {window}: the window \
             would span the whole alignment, so Phi_w would equal the mean \
             compatibility of every pair and the permutation null would be \
             identical to the observed value by construction. Lower the window \
             deliberately.",
            sites.len()
        ));
    }
    let present = if aln.len() >= 64 {
        u64::MAX
    } else {
        (1u64 << aln.len()) - 1
    };

    let identity: Vec<usize> = (0..sites.len()).collect();
    let (compatible_pairs, pairs) = phi_numerator(&sites, &identity, window, present);
    if pairs == 0 {
        return Err("no neighbouring pairs to score".to_string());
    }

    // The permutation null. The pair count is invariant under permutation, so
    // only the numerators are compared and the comparison is exact.
    let mut rng = SplitMix64::new(seed);
    let mut order = identity.clone();
    let mut at_least = 0usize;
    let mut numerator_total = 0u64;
    for _ in 0..permutations {
        // Fisher-Yates, drawing from the generator the rest of the crate uses.
        for i in (1..order.len()).rev() {
            let j = rng.below(i + 1);
            order.swap(i, j);
        }
        let (permuted, permuted_pairs) = phi_numerator(&sites, &order, window, present);
        debug_assert_eq!(permuted_pairs, pairs, "the pair count must be invariant");
        numerator_total += permuted as u64;
        if permuted >= compatible_pairs {
            at_least += 1;
        }
    }

    let p_value = (at_least + 1) as f64 / (permutations + 1) as f64;
    Ok(PhiTest {
        informative_sites: sites.len(),
        pairs_scored: pairs,
        compatible_pairs,
        phi_w_num: compatible_pairs as i64,
        phi_w_den: pairs as i64,
        phi_w: compatible_pairs as f64 / pairs as f64,
        permutations,
        permutations_at_least_observed: at_least,
        p_value,
        p_value_floor: 1.0 / (permutations + 1) as f64,
        permuted_mean_phi_w: if permutations == 0 {
            0.0
        } else {
            numerator_total as f64 / (permutations as f64 * pairs as f64)
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln(rows: &[&str]) -> Alignment {
        let labels: Vec<String> = (0..rows.len()).map(|i| format!("s{i}")).collect();
        Alignment::new(labels, rows.iter().map(|r| r.bytes().collect()).collect()).unwrap()
    }

    #[test]
    fn the_four_gamete_test_separates_compatible_from_incompatible() {
        let present = 0b1111u64;
        // Sites splitting {0,1} | {2,3} and {0,1} | {2,3}: nested, compatible.
        assert!(compatible(0b1100, 0b1100, present));
        // Nested: {2,3} inside {1,2,3}.
        assert!(compatible(0b1110, 0b1100, present));
        // All four gametes present: 00, 01, 10, 11 across the four sequences.
        assert!(
            !compatible(0b1010, 0b1100, present),
            "a site pair showing all four combinations must be incompatible"
        );
        // Disjoint non-empty splits with a remainder are still incompatible only
        // if all four cells are occupied; with three sequences they cannot be.
        let three = 0b111u64;
        assert!(compatible(0b110, 0b011, three));
    }

    #[test]
    fn informative_sites_skips_gaps_singletons_and_constant_columns() {
        // Column 0 constant; column 1 singleton; column 2 has a gap;
        // column 3 is biallelic 2/2.
        let a = aln(&["AAAA", "AC-A", "AAGC", "AAAC"]);
        let sites = informative_sites(&a).unwrap();
        assert_eq!(
            sites.len(),
            1,
            "only the last column is biallelic, gap-free and informative, got {sites:?}"
        );
        assert_eq!(sites[0].column, 3);
        // Sequences 2 and 3 carry the second state seen at that column.
        assert_eq!(sites[0].mask.count_ones(), 2);
    }

    #[test]
    fn the_window_must_be_smaller_than_the_site_count() {
        let a = aln(&["AACC", "AACC", "GGTT", "GGTT"]);
        let refused = phi_test(&a, 100, 10, 1);
        assert!(refused.is_err(), "a window spanning everything is refused");
        assert!(refused.unwrap_err().contains("window"));
    }

    #[test]
    fn phi_w_is_an_exact_ratio_of_counts() {
        use crate::simulate::clonal_tree;
        let pop = clonal_tree(2000, 3, 0.01, 11);
        let a = Alignment::new(pop.labels.clone(), pop.sequences.clone()).unwrap();
        let test = phi_test(&a, 20, 50, 7).expect("enough informative sites");
        assert_eq!(test.phi_w_den as usize, test.pairs_scored);
        assert_eq!(test.phi_w_num as usize, test.compatible_pairs);
        assert!(
            (0.0..=1.0).contains(&test.phi_w),
            "a mean of indicators is in [0,1], got {}",
            test.phi_w
        );
        assert!(test.p_value >= test.p_value_floor);
    }

    #[test]
    fn reciprocal_mosaics_raise_phi_w_above_its_permuted_level() {
        // The statistic's whole premise: recombination makes NEIGHBOURING sites
        // more compatible than the site order alone would explain. On a clonal
        // tree there is no such neighbourhood structure, so the observed value
        // should sit near the permuted mean; with mosaics it should sit above.
        use crate::simulate::{clonal_tree, with_reciprocal_recombinants};
        let clonal = clonal_tree(4000, 4, 0.01, 5);
        let (mosaic, _) =
            with_reciprocal_recombinants(clonal_tree(4000, 4, 0.01, 5), 0.5, 77);

        let null = Alignment::new(clonal.labels.clone(), clonal.sequences.clone()).unwrap();
        let positive = Alignment::new(mosaic.labels.clone(), mosaic.sequences.clone()).unwrap();
        let a = phi_test(&null, 50, 199, 3).expect("searchable");
        let b = phi_test(&positive, 50, 199, 3).expect("searchable");

        let lift_null = a.phi_w - a.permuted_mean_phi_w;
        let lift_positive = b.phi_w - b.permuted_mean_phi_w;
        assert!(
            lift_positive > lift_null,
            "mosaics must lift Phi_w above its permuted level more than a clonal \
             tree does: {lift_positive} against {lift_null}"
        );
        assert!(
            b.p_value <= a.p_value,
            "the positive control must not have a larger p-value than the null: \
             {} against {}",
            b.p_value,
            a.p_value
        );
    }
}

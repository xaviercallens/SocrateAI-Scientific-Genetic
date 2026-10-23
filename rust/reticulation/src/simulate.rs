//! Mechanistic nulls: generate clonal and recombinant populations.
//!
//! The ultrametric no-go gives the *exact* boundary of the null, but finite
//! sequences make distances stochastic, so real clonal data is never exactly
//! ultrametric and small margins occur under the null. Calibration is therefore
//! still needed. What changed relative to the textbook's Section 8 is the kind
//! of null available: a **generative model of clonal evolution** at the same
//! sequence length and divergence, rather than a permutation of the observed
//! data whose meaning was unclear.
//!
//! This is an improvement in interpretability, not an escape from calibration,
//! and the crate documentation says so in those terms.
//!
//! The controls this module exists to supply:
//!
//! * **negative** — [`clonal_tree`]: strictly clonal descent. Margins should be
//!   small, set by sequence length. A detector that fires here is measuring its
//!   own sampling noise.
//! * **positive** — [`with_recombination`]: one recombinant formed from two
//!   diverged parents. Must produce a large-margin certificate, or the detector
//!   has no power.

use crate::distance::{DistanceMatrix, DistanceModel};

/// A deterministic, well-mixed 64-bit generator (SplitMix64).
///
/// Written out rather than pulled in as a dependency so that a reported seed
/// reproduces a run exactly, with no dependence on a crate's version.
pub struct SplitMix64(u64);

impl SplitMix64 {
    /// Seed the generator.
    pub fn new(seed: u64) -> Self {
        SplitMix64(seed)
    }

    /// Next raw 64 bits.
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }

    /// Uniform in `[0, 1)`.
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform integer in `[0, n)`.
    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Mutate a sequence, substituting each site independently with probability `mu`.
fn mutate(seq: &[u8], mu: f64, rng: &mut SplitMix64) -> Vec<u8> {
    seq.iter()
        .map(|&b| {
            if rng.next_f64() < mu {
                // Substitute to a different base, so the rate is the rate.
                let mut nb = BASES[rng.below(4)];
                while nb == b {
                    nb = BASES[rng.below(4)];
                }
                nb
            } else {
                b
            }
        })
        .collect()
}

/// A simulated population: labelled, aligned sequences.
#[derive(Debug, Clone)]
pub struct Population {
    /// Sequence labels.
    pub labels: Vec<String>,
    /// Aligned sequences, all of the same length.
    pub sequences: Vec<Vec<u8>>,
}

impl Population {
    /// Pairwise distances under the given model.
    pub fn distances(&self, model: DistanceModel) -> Result<DistanceMatrix, String> {
        DistanceMatrix::from_alignment(self.labels.clone(), &self.sequences, model)
    }
}

/// **Negative control.** Strictly clonal descent on a balanced binary tree.
///
/// Starting from one random founder, each generation doubles the population and
/// every lineage accumulates substitutions at the same rate `mu_per_generation`.
/// Equal rates on a balanced tree is exactly the strict-clock, clonal case whose
/// expected distances are ultrametric; the realised distances deviate only by
/// multinomial sampling at `length` sites.
///
/// Returns the `2^generations` leaves.
pub fn clonal_tree(
    length: usize,
    generations: u32,
    mu_per_generation: f64,
    seed: u64,
) -> Population {
    let mut rng = SplitMix64::new(seed);
    let founder: Vec<u8> = (0..length).map(|_| BASES[rng.below(4)]).collect();

    let mut current = vec![founder];
    for _ in 0..generations {
        let mut next = Vec::with_capacity(current.len() * 2);
        for parent in &current {
            next.push(mutate(parent, mu_per_generation, &mut rng));
            next.push(mutate(parent, mu_per_generation, &mut rng));
        }
        current = next;
    }
    let labels = (0..current.len()).map(|i| format!("clonal_{i:03}")).collect();
    Population {
        labels,
        sequences: current,
    }
}

/// **Positive control.** A clonal population plus one recombinant.
///
/// Two leaves from opposite halves of the tree — hence maximally diverged — are
/// chosen as parents, and a recombinant is formed by taking the first
/// `breakpoint` fraction of its sequence from one and the remainder from the
/// other. The recombinant is close to each parent over the segment it inherited,
/// while the parents stay far from each other: the open path the certificate
/// detects.
///
/// Returns the population with the recombinant appended, and its label.
pub fn with_recombination(
    mut population: Population,
    breakpoint: f64,
    seed: u64,
) -> (Population, String) {
    let n = population.sequences.len();
    assert!(n >= 2, "need two parents");
    let mut rng = SplitMix64::new(seed);
    // Opposite halves of a balanced tree are the most diverged pairs.
    let left = rng.below(n / 2);
    let right = n / 2 + rng.below(n - n / 2);

    let length = population.sequences[left].len();
    let cut = ((breakpoint.clamp(0.05, 0.95)) * length as f64) as usize;
    let mut recombinant = Vec::with_capacity(length);
    recombinant.extend_from_slice(&population.sequences[left][..cut]);
    recombinant.extend_from_slice(&population.sequences[right][cut..]);

    let label = format!(
        "recombinant_of_{}_and_{}",
        population.labels[left], population.labels[right]
    );
    population.labels.push(label.clone());
    population.sequences.push(recombinant);
    (population, label)
}

/// **The adversarial control.** Clonal descent under a *relaxed* clock.
///
/// Every lineage independently rescales its substitution rate at each
/// generation, by `exp(spread * u)` with `u` uniform on `[-1, 1]`. The genealogy
/// is still exactly a tree and there is no reticulation whatsoever — only the
/// clock is broken.
///
/// This is the dataset that decides whether the ultrametric certificate measures
/// anything biological. `IsUltrametric` encodes "tree AND strict clock", so a
/// relaxed clock must violate it; if the resulting margins are comparable to a
/// recombinant's, then a triple certificate does not distinguish reticulation
/// from rate variation and the four-point test has to carry the method alone.
/// Risk 2 of `docs/GENETICA_RETICULATA.md` §7, made into an experiment.
///
/// Additivity on the tree is preserved, so the four-point condition should
/// survive. That contrast is the prediction being tested.
///
/// # The multipliers are normalised, deliberately
///
/// `exp(spread * u)` has mean `sinh(spread) / spread > 1`, so an unnormalised
/// relaxed clock is also a *faster* clock: the population diverges further and
/// p-distances begin to saturate. Saturation breaks additivity on its own, so an
/// unnormalised control would confound rate variation with saturation and a
/// four-point violation could not be attributed to either. The multipliers are
/// therefore divided by their mean, which leaves the expected branch rate — and
/// hence the expected divergence — equal to the strict-clock case while the
/// realised rates vary widely. [`clonal_tree`] at a high `mu` is the separate
/// probe for saturation.
pub fn clonal_tree_relaxed(
    length: usize,
    generations: u32,
    mu_per_generation: f64,
    spread: f64,
    seed: u64,
) -> Population {
    let mut rng = SplitMix64::new(seed);
    let founder: Vec<u8> = (0..length).map(|_| BASES[rng.below(4)]).collect();

    // E[exp(spread * u)] for u uniform on [-1, 1]; 1 in the limit spread -> 0.
    let mean_multiplier = if spread.abs() < 1e-12 {
        1.0
    } else {
        spread.sinh() / spread
    };

    // Each lineage carries its own rate multiplier, inherited and rescaled.
    let mut current = vec![(founder, 1.0f64)];
    for _ in 0..generations {
        let mut next = Vec::with_capacity(current.len() * 2);
        for (parent, rate) in &current {
            for _ in 0..2 {
                let u = 2.0 * rng.next_f64() - 1.0;
                let child_rate = rate * (spread * u).exp() / mean_multiplier;
                // Clamp to a probability: a rate above 1 is not a rate.
                let mu = (mu_per_generation * child_rate).clamp(0.0, 1.0);
                next.push((mutate(parent, mu, &mut rng), child_rate));
            }
        }
        current = next;
    }
    let labels = (0..current.len())
        .map(|i| format!("relaxed_{i:03}"))
        .collect();
    Population {
        labels,
        sequences: current.into_iter().map(|(s, _)| s).collect(),
    }
}

/// **The quartet-shaped positive control.** Two reciprocal mosaics of one pair.
///
/// A single recombinant breaks the ultrametric inequality on triples. The
/// four-point condition is about quadruples, and the shape that violates it most
/// cleanly is a pair of *reciprocal* mosaics: each is close to both parents, so
/// the two crossing pairings are cheap while the parents stay far apart. This is
/// the configuration `reciprocal_recombinants_witness` proves non-vacuous in
/// Lean.
///
/// Returns the population with both mosaics appended, and their labels.
pub fn with_reciprocal_recombinants(
    mut population: Population,
    breakpoint: f64,
    seed: u64,
) -> (Population, [String; 2]) {
    let n = population.sequences.len();
    assert!(n >= 2, "need two parents");
    let mut rng = SplitMix64::new(seed);
    let left = rng.below(n / 2);
    let right = n / 2 + rng.below(n - n / 2);

    let length = population.sequences[left].len();
    let cut = (breakpoint.clamp(0.05, 0.95) * length as f64) as usize;

    let mut first = Vec::with_capacity(length);
    first.extend_from_slice(&population.sequences[left][..cut]);
    first.extend_from_slice(&population.sequences[right][cut..]);
    let mut second = Vec::with_capacity(length);
    second.extend_from_slice(&population.sequences[right][..cut]);
    second.extend_from_slice(&population.sequences[left][cut..]);

    let labels = [
        format!(
            "mosaic_{}_then_{}",
            population.labels[left], population.labels[right]
        ),
        format!(
            "mosaic_{}_then_{}",
            population.labels[right], population.labels[left]
        ),
    ];
    population.labels.push(labels[0].clone());
    population.sequences.push(first);
    population.labels.push(labels[1].clone());
    population.sequences.push(second);
    (population, labels)
}

/// A population simulated with deletion events, carrying its own true alignment.
///
/// The point of this structure is that it holds *both* what an aligner would be
/// handed ([`IndelPopulation::raw`], gap-free sequences of unequal length) and
/// what the right answer is ([`IndelPopulation::truth`], rows in the founder's
/// coordinate frame). Phase 3 needs the pair, because the alignment step is the
/// one part of the pipeline whose errors the Lean certificate cannot catch: the
/// kernel verifies counts against the sequences it is given, and an alignment
/// error changes which sequences those are.
pub struct IndelPopulation {
    /// Labels, shared by `truth` and `raw` in the same order.
    pub labels: Vec<String>,
    /// The true alignment: rows of width `length`, `-` where a lineage has lost
    /// the site. This is ground truth, not an inference.
    pub truth: Population,
    /// The same sequences with gaps removed: unequal lengths, which is what an
    /// aligner actually sees.
    pub raw: Vec<Vec<u8>>,
    /// Number of deletion events that were realised along the tree.
    pub deletion_events: usize,
}

/// Substitute, leaving gaps gapped.
///
/// A site already lost to a deletion cannot mutate, and letting it would make
/// the truth frame incoherent.
fn mutate_preserving_gaps(seq: &[u8], mu: f64, rng: &mut SplitMix64) -> Vec<u8> {
    seq.iter()
        .map(|&b| {
            if b == b'-' {
                return b'-';
            }
            if rng.next_f64() < mu {
                let mut nb = BASES[rng.below(4)];
                while nb == b {
                    nb = BASES[rng.below(4)];
                }
                nb
            } else {
                b
            }
        })
        .collect()
}

/// **The aligner's control.** Clonal, relaxed-clock descent with deletions, and
/// the true alignment that goes with it.
///
/// Deletions are frame coordinates set to `-`, inherited by every descendant:
/// once a lineage has lost a site it stays lost, which is what makes the founder
/// frame a genuine alignment rather than a bookkeeping convenience. Insertions
/// are not simulated directly, and do not need to be: a deletion in the lineage
/// that is later chosen as alignment reference presents to the aligner exactly
/// as an insertion in every other sequence, so both directions are exercised.
///
/// There is **no reticulation here at all**. Every quartet violation found in
/// the inferred alignment of this population, and absent from its true
/// alignment, was manufactured by the aligner.
///
/// The deletion half of the generator's parameters.
pub struct DeletionModel {
    /// Per-lineage, per-generation probability of one deletion event.
    pub rate: f64,
    /// Expected deletion length in sites.
    pub mean_length: usize,
    /// A lineage refuses an event that would take its gap load past this
    /// fraction, so the truth frame stays comparable across the population.
    pub max_gap_fraction: f64,
}

pub fn clonal_tree_relaxed_with_deletions(
    length: usize,
    generations: u32,
    mu_per_generation: f64,
    spread: f64,
    deletions: &DeletionModel,
    seed: u64,
) -> IndelPopulation {
    let DeletionModel {
        rate: deletion_rate,
        mean_length: mean_deletion_length,
        max_gap_fraction,
    } = *deletions;
    let mut rng = SplitMix64::new(seed);
    let founder: Vec<u8> = (0..length).map(|_| BASES[rng.below(4)]).collect();
    let mean_multiplier = if spread.abs() < 1e-12 {
        1.0
    } else {
        spread.sinh() / spread
    };
    let cap = ((max_gap_fraction.clamp(0.0, 0.45)) * length as f64) as usize;

    let mut events = 0usize;
    let mut current = vec![(founder, 1.0f64)];
    for _ in 0..generations {
        let mut next = Vec::with_capacity(current.len() * 2);
        for (parent, rate) in &current {
            for _ in 0..2 {
                let u = 2.0 * rng.next_f64() - 1.0;
                let child_rate = rate * (spread * u).exp() / mean_multiplier;
                let mu = (mu_per_generation * child_rate).clamp(0.0, 1.0);
                let mut child = mutate_preserving_gaps(parent, mu, &mut rng);
                if deletion_rate > 0.0 && rng.next_f64() < deletion_rate {
                    // Geometric-ish length: mean `mean_deletion_length`, at
                    // least one site.
                    let span = 1 + rng.below(2 * mean_deletion_length.max(1));
                    let span = span.min(length);
                    let start = rng.below(length - span + 1);
                    let already = child.iter().filter(|&&b| b == b'-').count();
                    if already + span <= cap {
                        for site in child[start..start + span].iter_mut() {
                            *site = b'-';
                        }
                        events += 1;
                    }
                }
                next.push((child, child_rate));
            }
        }
        current = next;
    }

    let labels: Vec<String> = (0..current.len())
        .map(|i| format!("indel_{i:03}"))
        .collect();
    let rows: Vec<Vec<u8>> = current.into_iter().map(|(s, _)| s).collect();
    let raw: Vec<Vec<u8>> = rows
        .iter()
        .map(|row| row.iter().copied().filter(|&b| b != b'-').collect())
        .collect();
    IndelPopulation {
        labels: labels.clone(),
        truth: Population {
            labels,
            sequences: rows,
        },
        raw,
        deletion_events: events,
    }
}

/// Replace a contiguous block of each sequence with gaps.
///
/// Contiguous missing blocks are the realistic pattern — partial genes,
/// truncated reads, unsequenced termini — and they are exactly the pattern that
/// gives different pairs different comparable-site counts, which is the confound
/// of the audit's finding 6. One sequence is left intact so that every triple
/// retains some common coverage.
///
/// `fraction` is the share of each sequence to blank; the total gap load must
/// stay below half or [`crate::exact::validate_alphabet`] refuses the input.
pub fn inject_gap_blocks(mut population: Population, fraction: f64, seed: u64) -> Population {
    let mut rng = SplitMix64::new(seed);
    let length = population.sequences[0].len();
    let block = ((fraction.clamp(0.0, 0.9) * length as f64) as usize).min(length);
    if block == 0 {
        return population;
    }
    for seq in population.sequences.iter_mut().skip(1) {
        let start = rng.below(length - block + 1);
        for b in seq[start..start + block].iter_mut() {
            *b = b'-';
        }
    }
    population
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::certificate::find_certificates;

    #[test]
    fn generator_is_deterministic_and_well_spread() {
        let mut a = SplitMix64::new(42);
        let mut b = SplitMix64::new(42);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
        // Crude uniformity: the mean of many draws should be near 1/2.
        let mut r = SplitMix64::new(7);
        let mean: f64 = (0..20_000).map(|_| r.next_f64()).sum::<f64>() / 20_000.0;
        assert!((mean - 0.5).abs() < 0.02, "mean was {mean}");
    }

    #[test]
    fn mutation_rate_is_approximately_realised() {
        let mut rng = SplitMix64::new(3);
        let seq = vec![b'A'; 20_000];
        let mutated = mutate(&seq, 0.05, &mut rng);
        let changed = seq.iter().zip(&mutated).filter(|(a, b)| a != b).count();
        let realised = changed as f64 / seq.len() as f64;
        assert!(
            (realised - 0.05).abs() < 0.005,
            "realised rate {realised} far from 0.05"
        );
    }

    #[test]
    fn clonal_tree_has_the_expected_size_and_alignment() {
        let pop = clonal_tree(200, 4, 0.01, 11);
        assert_eq!(pop.sequences.len(), 16);
        assert!(pop.sequences.iter().all(|s| s.len() == 200));
        assert_eq!(pop.labels.len(), 16);
    }

    #[test]
    fn negative_control_margins_stay_small() {
        // Strictly clonal, long sequences: the realised distances are close to
        // ultrametric, so the largest margin should be a small fraction of the
        // distances themselves. This is the control that makes a positive
        // finding meaningful; if it fails, the detector is reading noise.
        let pop = clonal_tree(4000, 4, 0.004, 2026);
        let dm = pop.distances(DistanceModel::PDistance).unwrap();
        let report = find_certificates(&dm, 3);
        let typical = dm.get(0, dm.len() - 1).max(1e-9);
        assert!(
            report.ultrametric_defect < 0.5 * typical,
            "clonal defect {} is not small against typical distance {typical}",
            report.ultrametric_defect
        );
    }

    #[test]
    fn positive_control_is_detected_with_a_large_margin() {
        // A recombinant between diverged parents must produce a certificate
        // whose margin dwarfs the clonal baseline. Without this the negative
        // result above would be uninformative.
        let pop = clonal_tree(4000, 4, 0.004, 2026);
        let dm_clonal = pop.distances(DistanceModel::PDistance).unwrap();
        let clonal_defect = find_certificates(&dm_clonal, 1).ultrametric_defect;

        let pop = clonal_tree(4000, 4, 0.004, 2026);
        let (recombined, label) = with_recombination(pop, 0.5, 99);
        let dm = recombined.distances(DistanceModel::PDistance).unwrap();
        let report = find_certificates(&dm, 5);

        assert!(
            report.ultrametric_defect > 2.0 * clonal_defect,
            "recombinant defect {} should dwarf clonal {clonal_defect}",
            report.ultrametric_defect
        );
        assert!(report.all_certificates_validated);
        // The recombinant should appear as the midpoint of a top certificate:
        // it is the point close to both diverged parents.
        let names: Vec<&str> = report.certificates.iter().map(|c| c.y.as_str()).collect();
        assert!(
            names.iter().any(|y| *y == label.as_str()),
            "the recombinant {label} should be a midpoint; midpoints were {names:?}"
        );
    }

    /// Largest-to-smallest ratio among pairs separated by the root.
    ///
    /// Not the ratio over *all* pairs: on a balanced tree even a perfect clock
    /// makes sisters closer than cousins, so that ratio is large by design and
    /// says nothing about the clock. Pairs in opposite halves all cross the same
    /// number of branches, so under a strict clock their distances agree up to
    /// multinomial noise — and under a relaxed clock they do not. That is what
    /// "the clock is relaxed" means operationally.
    fn cross_root_spread(pop: &Population) -> f64 {
        let dm = pop.distances(DistanceModel::PDistance).unwrap();
        let half = dm.len() / 2;
        let mut lo = f64::INFINITY;
        let mut hi = 0.0f64;
        for i in 0..half {
            for j in half..dm.len() {
                lo = lo.min(dm.get(i, j));
                hi = hi.max(dm.get(i, j));
            }
        }
        hi / lo
    }

    #[test]
    fn the_relaxed_clock_really_relaxes_the_clock() {
        // The adversarial control has to be adversarial: if the rate spread did
        // not change the realised distances, the experiment that uses it would
        // be measuring nothing.
        let strict = cross_root_spread(&clonal_tree(4000, 4, 0.004, 2026));
        let relaxed = cross_root_spread(&clonal_tree_relaxed(4000, 4, 0.004, 1.0, 2026));
        assert!(
            relaxed > 2.0 * strict,
            "relaxed cross-root spread {relaxed} should dwarf strict {strict}"
        );
        // A spread of zero must put the clock back, so the parameter is the only
        // thing doing the work.
        let flat = cross_root_spread(&clonal_tree_relaxed(4000, 4, 0.004, 0.0, 2026));
        assert!(
            flat < 1.5,
            "with zero spread the cross-root distances should agree, got {flat}"
        );
        assert!(
            strict < 1.5,
            "the strict-clock baseline should itself be tight, got {strict}"
        );
    }

    #[test]
    fn the_relaxed_clock_is_divergence_matched_to_the_strict_one() {
        // Without normalising the multipliers, a relaxed clock is also a faster
        // clock, and a four-point violation could be attributed to saturation
        // rather than to rate variation. The control has to vary the rates
        // WITHOUT varying the expected divergence, so this is load-bearing.
        let mean = |pop: &Population| {
            let dm = pop.distances(DistanceModel::PDistance).unwrap();
            let mut s = 0.0;
            let mut n = 0.0;
            for i in 0..dm.len() {
                for j in (i + 1)..dm.len() {
                    s += dm.get(i, j);
                    n += 1.0;
                }
            }
            s / n
        };
        // Average over seeds: a single pair of populations is one draw, and the
        // claim is about the expectation.
        let (mut strict, mut relaxed) = (0.0, 0.0);
        for seed in 0..8u64 {
            strict += mean(&clonal_tree(4000, 4, 0.004, 1000 + seed));
            relaxed += mean(&clonal_tree_relaxed(4000, 4, 0.004, 1.0, 1000 + seed));
        }
        let ratio = relaxed / strict;
        assert!(
            (ratio - 1.0).abs() < 0.10,
            "mean divergence should match within 10%, got ratio {ratio}"
        );
    }

    #[test]
    fn reciprocal_mosaics_are_each_close_to_both_parents() {
        let pop = clonal_tree(4000, 4, 0.004, 2026);
        let n_before = pop.sequences.len();
        let (pop, mosaics) = with_reciprocal_recombinants(pop, 0.5, 99);
        assert_eq!(pop.sequences.len(), n_before + 2);
        let dm = pop.distances(DistanceModel::PDistance).unwrap();
        let idx = |name: &str| pop.labels.iter().position(|l| l == name).unwrap();
        // Recover the parents from the generated labels.
        let parts: Vec<&str> = mosaics[0]
            .strip_prefix("mosaic_")
            .unwrap()
            .split("_then_")
            .collect();
        let (p, q) = (idx(parts[0]), idx(parts[1]));
        let parental = dm.get(p, q);
        for m in &mosaics {
            let mi = idx(m);
            assert!(
                dm.get(mi, p) < parental && dm.get(mi, q) < parental,
                "{m} must sit inside the parental distance {parental}"
            );
        }
        // The two mosaics are reciprocal, so they are about as far apart as the
        // parents: that is what makes the parent pairing dominate.
        assert!(dm.get(idx(&mosaics[0]), idx(&mosaics[1])) > 0.5 * parental);
    }

    #[test]
    fn gap_blocks_give_pairs_different_comparable_counts() {
        let pop = clonal_tree(400, 3, 0.01, 5);
        let intact = pop.sequences[0].clone();
        let gapped = inject_gap_blocks(pop, 0.25, 11);
        assert_eq!(
            gapped.sequences[0], intact,
            "one sequence is left intact so triples keep common coverage"
        );
        let dash = |s: &Vec<u8>| s.iter().filter(|b| **b == b'-').count();
        assert_eq!(dash(&gapped.sequences[1]), 100, "25% of 400 sites");
        // The confound the gaps exist to create: per-pair denominators differ.
        let d = crate::exact::ExactDistances::from_alignment(
            gapped.labels.clone(),
            &gapped.sequences,
        )
        .unwrap();
        let dens: std::collections::BTreeSet<i64> = (0..d.len())
            .flat_map(|i| ((i + 1)..d.len()).map(move |j| (i, j)))
            .map(|(i, j)| d.ratio(i, j).den)
            .collect();
        assert!(
            dens.len() > 1,
            "gap blocks must produce differing denominators, got {dens:?}"
        );
    }

    #[test]
    fn with_no_deletion_rate_the_truth_frame_is_gapless_and_raw_equals_truth() {
        let pop = clonal_tree_relaxed_with_deletions(
            500,
            3,
            0.004,
            1.0,
            &DeletionModel {
                rate: 0.0,
                mean_length: 10,
                max_gap_fraction: 0.3,
            },
            11,
        );
        assert_eq!(pop.deletion_events, 0, "no events at rate 0");
        assert_eq!(pop.truth.sequences.len(), 8);
        assert!(
            pop.truth.sequences.iter().all(|s| !s.contains(&b'-')),
            "the truth frame must be gapless when nothing is deleted"
        );
        assert_eq!(
            pop.raw, pop.truth.sequences,
            "stripping gaps from a gapless frame is the identity"
        );
    }

    #[test]
    fn deletions_appear_in_the_truth_frame_and_shorten_the_raw_sequences() {
        let pop = clonal_tree_relaxed_with_deletions(
            1000,
            4,
            0.004,
            1.0,
            &DeletionModel {
                rate: 0.6,
                mean_length: 20,
                max_gap_fraction: 0.3,
            },
            5,
        );
        assert!(
            pop.deletion_events > 0,
            "a rate of 0.6 over 4 generations must realise some events"
        );
        assert!(
            pop.truth.sequences.iter().all(|s| s.len() == 1000),
            "every truth row keeps the founder's width"
        );
        let shortened = pop.raw.iter().filter(|s| s.len() < 1000).count();
        assert!(
            shortened > 0,
            "deletions must shorten at least one raw sequence"
        );
        for (row, raw) in pop.truth.sequences.iter().zip(&pop.raw) {
            let gaps = row.iter().filter(|&&b| b == b'-').count();
            assert_eq!(
                row.len() - gaps,
                raw.len(),
                "the raw sequence is exactly the truth row minus its gaps"
            );
            assert!(
                gaps <= 300,
                "the 30% cap must be respected, got {gaps} gaps"
            );
        }
    }
}

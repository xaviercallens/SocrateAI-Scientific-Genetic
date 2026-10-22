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
}

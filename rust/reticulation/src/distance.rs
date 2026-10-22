//! Pairwise genetic distances over an alignment.
//!
//! Only models whose assumptions we can state are included. The p-distance is
//! the raw proportion of differing sites; Jukes–Cantor corrects it for multiple
//! substitutions at the same site under the assumption that all substitutions
//! are equally likely, which is wrong in detail but is the standard first
//! correction and is monotone in p, so it cannot create or destroy an
//! ultrametric violation on its own (see `jukes_cantor_is_monotone`).

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Which substitution model converts site differences into a distance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DistanceModel {
    /// Raw proportion of differing comparable sites.
    PDistance,
    /// Jukes–Cantor correction `-3/4 ln(1 - 4p/3)`, saturating below `p = 3/4`.
    JukesCantor,
}

impl DistanceModel {
    /// Apply the model to a raw proportion of differing sites.
    ///
    /// Jukes–Cantor diverges as `p -> 3/4`; beyond that the correction is
    /// undefined. We saturate rather than return infinity, and report the
    /// saturation count so a caller cannot silently analyse saturated data.
    pub fn apply(self, p: f64) -> (f64, bool) {
        match self {
            DistanceModel::PDistance => (p, false),
            DistanceModel::JukesCantor => {
                let arg = 1.0 - (4.0 / 3.0) * p;
                if arg <= 1e-12 {
                    // Saturated: the model cannot distinguish these sequences.
                    (3.0, true)
                } else {
                    (-0.75 * arg.ln(), false)
                }
            }
        }
    }
}

/// A symmetric distance matrix with zero diagonal, plus its provenance.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DistanceMatrix {
    /// Sequence labels, in the order the matrix indexes them.
    pub labels: Vec<String>,
    /// Row-major `n * n` distances.
    pub values: Vec<f64>,
    /// The model used.
    pub model: DistanceModel,
    /// Alignment length in columns.
    pub alignment_length: usize,
    /// How many pairs saturated the model, if any. Non-zero invalidates
    /// quantitative use of those pairs and is reported rather than hidden.
    pub saturated_pairs: usize,
}

impl DistanceMatrix {
    /// Number of sequences.
    pub fn len(&self) -> usize {
        self.labels.len()
    }

    /// Whether the matrix is empty.
    pub fn is_empty(&self) -> bool {
        self.labels.is_empty()
    }

    /// Distance between two indices.
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.values[i * self.labels.len() + j]
    }

    /// Build from an alignment: equal-length sequences over `ACGT-N`.
    ///
    /// Sites where either sequence carries a gap or an ambiguity are skipped
    /// pairwise, which is the usual convention. A pair with no comparable sites
    /// is given distance zero and counted as saturated, because "no information"
    /// must not masquerade as "identical".
    ///
    /// # Errors
    ///
    /// Returns an error if fewer than three sequences are supplied (a triple is
    /// the smallest object that can carry a certificate) or if the sequences are
    /// not all the same length.
    pub fn from_alignment(
        labels: Vec<String>,
        sequences: &[Vec<u8>],
        model: DistanceModel,
    ) -> Result<Self, String> {
        if sequences.len() < 3 {
            return Err(format!(
                "need at least 3 sequences to form a triple, got {}",
                sequences.len()
            ));
        }
        if labels.len() != sequences.len() {
            return Err(format!(
                "{} labels for {} sequences",
                labels.len(),
                sequences.len()
            ));
        }
        let length = sequences[0].len();
        if let Some(bad) = sequences.iter().position(|s| s.len() != length) {
            return Err(format!(
                "sequence {bad} has length {} but the first has {length}; \
                 this expects an alignment, not raw reads",
                sequences[bad].len()
            ));
        }

        let n = sequences.len();
        // Upper triangle in parallel; each entry is independent.
        let pairs: Vec<(usize, usize, f64, bool)> = (0..n)
            .into_par_iter()
            .flat_map_iter(|i| (i + 1..n).map(move |j| (i, j)))
            .map(|(i, j)| {
                let (diff, comparable) = count_differences(&sequences[i], &sequences[j]);
                if comparable == 0 {
                    return (i, j, 0.0, true);
                }
                let p = diff as f64 / comparable as f64;
                let (d, saturated) = model.apply(p);
                (i, j, d, saturated)
            })
            .collect();

        let mut values = vec![0.0; n * n];
        let mut saturated_pairs = 0usize;
        for (i, j, d, saturated) in pairs {
            values[i * n + j] = d;
            values[j * n + i] = d;
            if saturated {
                saturated_pairs += 1;
            }
        }

        Ok(DistanceMatrix {
            labels,
            values,
            model,
            alignment_length: length,
            saturated_pairs,
        })
    }

    /// Build directly from a symmetric matrix, for testing and for data that
    /// arrives as distances rather than as an alignment.
    ///
    /// # Errors
    ///
    /// Returns an error if the matrix is not square, not symmetric, or has a
    /// non-zero diagonal. These are checked rather than assumed because every
    /// theorem downstream depends on them.
    pub fn from_values(labels: Vec<String>, values: Vec<f64>) -> Result<Self, String> {
        let n = labels.len();
        if values.len() != n * n {
            return Err(format!("{} values for {n} labels", values.len()));
        }
        for i in 0..n {
            if values[i * n + i] != 0.0 {
                return Err(format!("diagonal entry {i} is {}, not 0", values[i * n + i]));
            }
            for j in 0..n {
                let (a, b) = (values[i * n + j], values[j * n + i]);
                if (a - b).abs() > 1e-12 {
                    return Err(format!("not symmetric at ({i},{j}): {a} vs {b}"));
                }
            }
        }
        Ok(DistanceMatrix {
            labels,
            values,
            model: DistanceModel::PDistance,
            alignment_length: 0,
            saturated_pairs: 0,
        })
    }
}

/// Count differing and comparable sites between two equal-length sequences.
///
/// Comparable means both sides carry one of `ACGT` (case-insensitive). Gaps
/// (`-`), `N` and any other symbol are skipped.
fn count_differences(a: &[u8], b: &[u8]) -> (usize, usize) {
    let mut diff = 0usize;
    let mut comparable = 0usize;
    for (&x, &y) in a.iter().zip(b.iter()) {
        let (x, y) = (x.to_ascii_uppercase(), y.to_ascii_uppercase());
        if !matches!(x, b'A' | b'C' | b'G' | b'T') || !matches!(y, b'A' | b'C' | b'G' | b'T') {
            continue;
        }
        comparable += 1;
        if x != y {
            diff += 1;
        }
    }
    (diff, comparable)
}

/// Jukes–Cantor is strictly increasing in `p` below saturation.
///
/// This matters for the certificate logic: a monotone reparameterisation of the
/// distances cannot turn an ultrametric into a non-ultrametric or the reverse,
/// because the ultrametric condition depends only on the *ordering* of the three
/// distances in a triple. So the choice between the two models cannot manufacture
/// a certificate, and this function exists so the claim is testable.
pub fn jukes_cantor_is_monotone(p1: f64, p2: f64) -> bool {
    let (d1, s1) = DistanceModel::JukesCantor.apply(p1);
    let (d2, s2) = DistanceModel::JukesCantor.apply(p2);
    if s1 || s2 {
        return true;
    }
    (p1 < p2) == (d1 < d2) || (p1 - p2).abs() < 1e-15
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_sequences_have_zero_distance() {
        let seqs = vec![b"ACGTACGT".to_vec(), b"ACGTACGT".to_vec(), b"ACGTACGT".to_vec()];
        let labels = vec!["a".into(), "b".into(), "c".into()];
        let dm = DistanceMatrix::from_alignment(labels, &seqs, DistanceModel::PDistance).unwrap();
        assert_eq!(dm.get(0, 1), 0.0);
        assert_eq!(dm.get(0, 0), 0.0);
        assert_eq!(dm.saturated_pairs, 0);
    }

    #[test]
    fn p_distance_counts_differing_sites() {
        // 2 of 8 sites differ.
        let seqs = vec![
            b"ACGTACGT".to_vec(),
            b"ACGTACGA".to_vec(),
            b"TCGTACGT".to_vec(),
        ];
        let labels = vec!["a".into(), "b".into(), "c".into()];
        let dm = DistanceMatrix::from_alignment(labels, &seqs, DistanceModel::PDistance).unwrap();
        assert!((dm.get(0, 1) - 1.0 / 8.0).abs() < 1e-12);
        assert!((dm.get(0, 2) - 1.0 / 8.0).abs() < 1e-12);
        assert!((dm.get(1, 2) - 2.0 / 8.0).abs() < 1e-12);
    }

    #[test]
    fn gaps_and_ambiguities_are_skipped_pairwise() {
        let seqs = vec![
            b"ACGT".to_vec(),
            b"AN-T".to_vec(),
            b"ACGT".to_vec(),
        ];
        let labels = vec!["a".into(), "b".into(), "c".into()];
        let dm = DistanceMatrix::from_alignment(labels, &seqs, DistanceModel::PDistance).unwrap();
        // Only sites 0 and 3 are comparable for the pair (0,1), and both agree.
        assert_eq!(dm.get(0, 1), 0.0);
        assert_eq!(dm.saturated_pairs, 0);
    }

    #[test]
    fn no_comparable_sites_is_flagged_not_called_identical() {
        let seqs = vec![b"NN".to_vec(), b"NN".to_vec(), b"AC".to_vec()];
        let labels = vec!["a".into(), "b".into(), "c".into()];
        let dm = DistanceMatrix::from_alignment(labels, &seqs, DistanceModel::PDistance).unwrap();
        assert!(dm.saturated_pairs >= 1, "must flag uninformative pairs");
    }

    #[test]
    fn unequal_lengths_are_rejected() {
        let seqs = vec![b"ACGT".to_vec(), b"ACG".to_vec(), b"ACGT".to_vec()];
        let labels = vec!["a".into(), "b".into(), "c".into()];
        assert!(DistanceMatrix::from_alignment(labels, &seqs, DistanceModel::PDistance).is_err());
    }

    #[test]
    fn fewer_than_three_sequences_is_rejected() {
        let seqs = vec![b"ACGT".to_vec(), b"ACGA".to_vec()];
        let labels = vec!["a".into(), "b".into()];
        assert!(DistanceMatrix::from_alignment(labels, &seqs, DistanceModel::PDistance).is_err());
    }

    #[test]
    fn asymmetric_input_is_rejected() {
        let labels = vec!["a".into(), "b".into(), "c".into()];
        let values = vec![0.0, 1.0, 2.0, 9.0, 0.0, 1.0, 2.0, 1.0, 0.0];
        assert!(DistanceMatrix::from_values(labels, values).is_err());
    }

    #[test]
    fn jukes_cantor_monotone_and_saturating() {
        assert!(jukes_cantor_is_monotone(0.01, 0.10));
        assert!(jukes_cantor_is_monotone(0.10, 0.30));
        let (_, saturated) = DistanceModel::JukesCantor.apply(0.80);
        assert!(saturated, "p > 3/4 must saturate rather than return NaN");
        let (d, _) = DistanceModel::JukesCantor.apply(0.0);
        assert_eq!(d, 0.0);
    }
}

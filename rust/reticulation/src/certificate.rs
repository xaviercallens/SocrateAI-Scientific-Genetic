//! Reticulation certificates: searching for proved violations of the null.
//!
//! # The statistic
//!
//! For an unordered triple, write the three pairwise distances sorted as
//! `d_min <= d_mid <= d_max`. The ultrametric condition says the largest is at
//! most the larger of the other two, i.e. `d_max <= d_mid`. So the triple
//! violates the null exactly when
//!
//! ```text
//! margin = d_max - d_mid > 0
//! ```
//!
//! and the *ultrametric defect* of the whole matrix is the maximum margin over
//! all triples. It is zero precisely when the matrix is an ultrametric — proved
//! as `ultrametric_iff_no_positive_margin` in
//! `formal/ANSE/ReticulationCertificate.lean`. There is no scale parameter, no
//! threshold and no ensemble in that statement.
//!
//! A violating triple maps exactly onto the fields of the Lean `Certificate`
//! structure: take `x` and `z` to be the pair realising `d_max`, `y` the
//! remaining point, and `scale = d_mid`. Then `d(x,y) <= scale`,
//! `d(y,z) <= scale` and `scale < d(x,z)`, which are its three hypotheses.
//!
//! # Honest limits
//!
//! Finite sequences make distances stochastic, so margins are positive under the
//! null too. The defect is an exact statistic about the matrix, not a test
//! against sampling noise. See [`crate::simulate`] for the calibration, and the
//! crate-level docs for why we consider a mechanistic null an improvement over a
//! shuffle but not a substitute for calibration.

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::distance::DistanceMatrix;

/// A proved violation of ultrametricity: an open path at a scale.
///
/// The field names and the three inequalities match
/// `ANSE.ReticulationCertificate.Certificate` so that a reader can transcribe
/// one into the theorem and recheck it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Certificate {
    /// Label of one endpoint.
    pub x: String,
    /// Label of the midpoint, close to both endpoints.
    pub y: String,
    /// Label of the other endpoint.
    pub z: String,
    /// Indices into the distance matrix, for machine rechecking.
    pub indices: [usize; 3],
    /// The scale at which the path is open: `d_mid`.
    pub scale: f64,
    /// `d(x,y) <= scale`.
    pub d_xy: f64,
    /// `d(y,z) <= scale`.
    pub d_yz: f64,
    /// `d(x,z) > scale`.
    pub d_xz: f64,
    /// `d_xz - scale`, strictly positive.
    pub margin: f64,
}

impl Certificate {
    /// Recheck the three inequalities. The search is not trusted; this is.
    ///
    /// Returns `Ok(())` only if the object really is a certificate in the sense
    /// of the Lean theorem. Called on every certificate before it is reported,
    /// so a bug in the search cannot emit a non-certificate.
    pub fn validate(&self) -> Result<(), String> {
        if !(self.d_xy <= self.scale) {
            return Err(format!("d_xy = {} exceeds scale {}", self.d_xy, self.scale));
        }
        if !(self.d_yz <= self.scale) {
            return Err(format!("d_yz = {} exceeds scale {}", self.d_yz, self.scale));
        }
        if !(self.scale < self.d_xz) {
            return Err(format!(
                "d_xz = {} does not exceed scale {}",
                self.d_xz, self.scale
            ));
        }
        let expected = self.d_xz - self.scale;
        if (expected - self.margin).abs() > 1e-12 {
            return Err(format!(
                "margin {} does not equal d_xz - scale = {expected}",
                self.margin
            ));
        }
        Ok(())
    }
}

/// A four-point (tree-metric) violation. Weaker evidence than a [`Certificate`]:
/// the topological consequence for tree metrics is cited in the textbook, not
/// proved in Lean, so these are reported as *suggestive* and never as proved.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FourPointDefect {
    /// The four labels involved.
    pub labels: [String; 4],
    /// Indices into the distance matrix.
    pub indices: [usize; 4],
    /// `S_max - S_mid` over the three pairings, strictly positive.
    pub defect: f64,
    /// Whether every quadruple was examined, or a sample.
    pub exhaustive: bool,
}

/// The result of a search, with its controls and its own limits recorded.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CertificateReport {
    /// Number of sequences.
    pub n_sequences: usize,
    /// Maximum margin over all triples. Zero iff the matrix is an ultrametric.
    pub ultrametric_defect: f64,
    /// How many triples violate the null at all.
    pub violating_triples: usize,
    /// Total number of triples examined (always exhaustive).
    pub triples_examined: usize,
    /// The strongest certificates, most severe first.
    pub certificates: Vec<Certificate>,
    /// Maximum four-point defect, or `None` if not computed.
    pub four_point_defect: Option<f64>,
    /// The strongest four-point violations.
    pub four_point: Vec<FourPointDefect>,
    /// Whether the four-point search was exhaustive.
    pub four_point_exhaustive: bool,
    /// Every reported certificate passed [`Certificate::validate`].
    pub all_certificates_validated: bool,
    /// Pairs that saturated the distance model; non-zero makes margins suspect.
    pub saturated_pairs: usize,
}

/// Search every triple for ultrametric violations.
///
/// Exhaustive and exact: `O(n^3)` comparisons, parallel over the first index.
/// `top_k` bounds how many certificates are *returned*, not how many are
/// counted, and the bound is recorded in the report so a truncated list cannot
/// read as a complete one.
pub fn find_certificates(dm: &DistanceMatrix, top_k: usize) -> CertificateReport {
    let n = dm.len();

    // Each triple independently; collect (margin, certificate) for violations.
    let per_i: Vec<(f64, usize, Vec<Certificate>)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut best = 0.0f64;
            let mut count = 0usize;
            let mut found: Vec<Certificate> = Vec::new();
            for j in (i + 1)..n {
                for k in (j + 1)..n {
                    let (d_ij, d_ik, d_jk) = (dm.get(i, j), dm.get(i, k), dm.get(j, k));
                    // Identify the largest distance; its endpoints are x and z,
                    // and the remaining point is the midpoint y. The scale is
                    // the second largest, which is then an upper bound for both
                    // of the midpoint's distances.
                    let (x, y, z, d_xz) = if d_ij >= d_ik && d_ij >= d_jk {
                        (i, k, j, d_ij)
                    } else if d_ik >= d_ij && d_ik >= d_jk {
                        (i, j, k, d_ik)
                    } else {
                        (j, i, k, d_jk)
                    };
                    let d_xy = dm.get(x, y);
                    let d_yz = dm.get(y, z);
                    let scale = d_xy.max(d_yz);
                    let margin = d_xz - scale;
                    if margin > 0.0 {
                        count += 1;
                        if margin > best {
                            best = margin;
                        }
                        found.push(Certificate {
                            x: dm.labels[x].clone(),
                            y: dm.labels[y].clone(),
                            z: dm.labels[z].clone(),
                            indices: [x, y, z],
                            scale,
                            d_xy,
                            d_yz,
                            d_xz,
                            margin,
                        });
                        // Keep the local list bounded so memory stays linear in
                        // top_k rather than in the number of violations.
                        if found.len() > 8 * top_k.max(1) {
                            found.sort_by(|a, b| b.margin.total_cmp(&a.margin));
                            found.truncate(top_k.max(1));
                        }
                    }
                }
            }
            (best, count, found)
        })
        .collect();

    let ultrametric_defect = per_i.iter().map(|(b, _, _)| *b).fold(0.0f64, f64::max);
    let violating_triples = per_i.iter().map(|(_, c, _)| *c).sum();
    let mut certificates: Vec<Certificate> =
        per_i.into_iter().flat_map(|(_, _, f)| f).collect();
    certificates.sort_by(|a, b| b.margin.total_cmp(&a.margin));
    certificates.truncate(top_k);

    // Nothing is reported that has not been rechecked against the Lean
    // theorem's hypotheses.
    let all_validated = certificates.iter().all(|c| c.validate().is_ok());
    certificates.retain(|c| c.validate().is_ok());

    let triples_examined = if n >= 3 { n * (n - 1) * (n - 2) / 6 } else { 0 };

    CertificateReport {
        n_sequences: n,
        ultrametric_defect,
        violating_triples,
        triples_examined,
        certificates,
        four_point_defect: None,
        four_point: Vec::new(),
        four_point_exhaustive: false,
        all_certificates_validated: all_validated,
        saturated_pairs: dm.saturated_pairs,
    }
}

/// Add the four-point (tree-metric) statistic to a report.
///
/// Exhaustive when `n <= max_exhaustive_n`, otherwise a deterministic stride
/// sample. The `exhaustive` flag is recorded; a sampled defect is a lower bound
/// on the true one and is labelled as such rather than quoted as the maximum.
pub fn add_four_point(
    report: &mut CertificateReport,
    dm: &DistanceMatrix,
    top_k: usize,
    max_exhaustive_n: usize,
) {
    let n = dm.len();
    if n < 4 {
        report.four_point_defect = Some(0.0);
        report.four_point_exhaustive = true;
        return;
    }
    let exhaustive = n <= max_exhaustive_n;
    // A deterministic stride keeps the sample reproducible without an RNG. The
    // stride is coprime with small periods so it does not lock onto a pattern
    // in the input ordering -- the same mistake as a codon-aligned stride.
    let stride = if exhaustive { 1 } else { 7 };

    let per_i: Vec<(f64, Vec<FourPointDefect>)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut best = 0.0f64;
            let mut found: Vec<FourPointDefect> = Vec::new();
            let mut counter = 0usize;
            for j in (i + 1)..n {
                for k in (j + 1)..n {
                    for l in (k + 1)..n {
                        counter += 1;
                        if counter % stride != 0 {
                            continue;
                        }
                        let s1 = dm.get(i, j) + dm.get(k, l);
                        let s2 = dm.get(i, k) + dm.get(j, l);
                        let s3 = dm.get(i, l) + dm.get(j, k);
                        let mut sums = [s1, s2, s3];
                        sums.sort_by(f64::total_cmp);
                        let defect = sums[2] - sums[1];
                        if defect > 0.0 {
                            if defect > best {
                                best = defect;
                            }
                            found.push(FourPointDefect {
                                labels: [
                                    dm.labels[i].clone(),
                                    dm.labels[j].clone(),
                                    dm.labels[k].clone(),
                                    dm.labels[l].clone(),
                                ],
                                indices: [i, j, k, l],
                                defect,
                                exhaustive,
                            });
                            if found.len() > 8 * top_k.max(1) {
                                found.sort_by(|a, b| b.defect.total_cmp(&a.defect));
                                found.truncate(top_k.max(1));
                            }
                        }
                    }
                }
            }
            (best, found)
        })
        .collect();

    let defect = per_i.iter().map(|(b, _)| *b).fold(0.0f64, f64::max);
    let mut all: Vec<FourPointDefect> = per_i.into_iter().flat_map(|(_, f)| f).collect();
    all.sort_by(|a, b| b.defect.total_cmp(&a.defect));
    all.truncate(top_k);

    report.four_point_defect = Some(defect);
    report.four_point = all;
    report.four_point_exhaustive = exhaustive;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::distance::DistanceMatrix;

    fn matrix(labels: &[&str], values: Vec<f64>) -> DistanceMatrix {
        DistanceMatrix::from_values(labels.iter().map(|s| s.to_string()).collect(), values).unwrap()
    }

    /// A genuine ultrametric: the two largest distances in every triple agree.
    fn ultrametric_3() -> DistanceMatrix {
        matrix(&["a", "b", "c"], vec![0.0, 1.0, 2.0, 1.0, 0.0, 2.0, 2.0, 2.0, 0.0])
    }

    /// The minimal open path: endpoints far, both close to the midpoint.
    fn open_path_3() -> DistanceMatrix {
        matrix(&["a", "m", "b"], vec![0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0])
    }

    #[test]
    fn ultrametric_yields_no_certificate() {
        // The control. This is `ultrametric_admits_no_certificate` in Lean.
        let report = find_certificates(&ultrametric_3(), 10);
        assert_eq!(report.ultrametric_defect, 0.0);
        assert_eq!(report.violating_triples, 0);
        assert!(report.certificates.is_empty());
        assert_eq!(report.triples_examined, 1);
    }

    #[test]
    fn open_path_yields_a_validated_certificate() {
        // The positive control, and the Lean non-vacuity witness.
        let report = find_certificates(&open_path_3(), 10);
        assert_eq!(report.violating_triples, 1);
        assert!((report.ultrametric_defect - 2.0).abs() < 1e-12);
        let c = &report.certificates[0];
        assert_eq!(c.y, "m", "the midpoint must be the recombinant-like point");
        assert!((c.scale - 1.0).abs() < 1e-12);
        assert!((c.d_xz - 3.0).abs() < 1e-12);
        assert!((c.margin - 2.0).abs() < 1e-12);
        assert!(c.validate().is_ok());
        assert!(report.all_certificates_validated);
    }

    #[test]
    fn certificate_fields_satisfy_the_lean_hypotheses() {
        // The three inequalities are exactly Certificate's fields in Lean.
        let report = find_certificates(&open_path_3(), 10);
        for c in &report.certificates {
            assert!(c.d_xy <= c.scale);
            assert!(c.d_yz <= c.scale);
            assert!(c.scale < c.d_xz);
            assert!(c.margin > 0.0);
        }
    }

    #[test]
    fn validate_rejects_a_forged_certificate() {
        // A search bug must not be able to emit a passing non-certificate.
        let mut forged = Certificate {
            x: "a".into(),
            y: "m".into(),
            z: "b".into(),
            indices: [0, 1, 2],
            scale: 5.0,
            d_xy: 1.0,
            d_yz: 1.0,
            d_xz: 3.0,
            margin: 2.0,
        };
        assert!(forged.validate().is_err(), "scale >= d_xz must be rejected");
        forged.scale = 1.0;
        forged.margin = 99.0;
        assert!(forged.validate().is_err(), "wrong margin must be rejected");
        forged.margin = 2.0;
        assert!(forged.validate().is_ok());
    }

    #[test]
    fn defect_is_zero_exactly_on_ultrametrics() {
        // `ultrametric_iff_no_positive_margin`, checked numerically.
        assert_eq!(find_certificates(&ultrametric_3(), 1).ultrametric_defect, 0.0);
        assert!(find_certificates(&open_path_3(), 1).ultrametric_defect > 0.0);
    }

    #[test]
    fn four_point_vanishes_on_an_additive_metric() {
        // Four points on a line: 0, 1, 2, 3. A line is a tree, so the
        // four-point condition holds and the defect must be zero.
        let labels = vec!["p0".into(), "p1".into(), "p2".into(), "p3".into()];
        let pts: [f64; 4] = [0.0, 1.0, 2.0, 3.0];
        let mut values = vec![0.0; 16];
        for i in 0..4 {
            for j in 0..4 {
                values[i * 4 + j] = (pts[i] - pts[j]).abs();
            }
        }
        let dm = DistanceMatrix::from_values(labels, values).unwrap();
        let mut report = find_certificates(&dm, 5);
        add_four_point(&mut report, &dm, 5, 64);
        assert!(report.four_point_exhaustive);
        assert!(
            report.four_point_defect.unwrap() < 1e-12,
            "a line is a tree metric; defect was {:?}",
            report.four_point_defect
        );
        // But the line is NOT an ultrametric, so certificates exist. This is
        // the pair of statistics doing their job: rate variation without
        // reticulation.
        assert!(report.ultrametric_defect > 0.0);
    }

    #[test]
    fn four_point_detects_a_non_tree_metric() {
        // The 4-cycle metric: opposite corners at distance 2, adjacent at 1.
        // This is the standard non-additive example.
        let labels = vec!["a".into(), "b".into(), "c".into(), "d".into()];
        let values = vec![
            0.0, 1.0, 2.0, 1.0, //
            1.0, 0.0, 1.0, 2.0, //
            2.0, 1.0, 0.0, 1.0, //
            1.0, 2.0, 1.0, 0.0,
        ];
        let dm = DistanceMatrix::from_values(labels, values).unwrap();
        let mut report = find_certificates(&dm, 5);
        add_four_point(&mut report, &dm, 5, 64);
        assert!(
            report.four_point_defect.unwrap() > 0.0,
            "the 4-cycle is not a tree metric"
        );
    }
}

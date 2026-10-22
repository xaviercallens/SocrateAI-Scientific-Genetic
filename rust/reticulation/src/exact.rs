//! Exact arithmetic on the evidence path. No floating point reaches a
//! reported inequality.
//!
//! A p-distance is `(differing sites) / (comparable sites)` — a rational with
//! an integer numerator and denominator. Pairwise gap deletion gives different
//! pairs different comparable-site counts, so denominators differ per pair; the
//! comparisons are therefore **cross-multiplied** into integer comparisons,
//! exactly as `ANSE.GeneticaReticulata.ratio_lt_ratio_iff` licenses.
//!
//! All cross-products are computed in `i128`. With counts bounded by alignment
//! length, a product of two counts cannot approach `i128::MAX`, and
//! [`Ratio::checked_cmp`] returns `None` rather than wrapping if it ever could.

use serde::{Deserialize, Serialize};

/// An exact non-negative rational: `num / den` with `den > 0`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Ratio {
    /// Numerator — differing sites.
    pub num: i64,
    /// Denominator — comparable sites, strictly positive.
    pub den: i64,
}

impl Ratio {
    /// Build a ratio, rejecting a non-positive denominator.
    pub fn new(num: i64, den: i64) -> Result<Self, String> {
        if den <= 0 {
            return Err(format!("denominator {den} must be positive"));
        }
        if num < 0 {
            return Err(format!("numerator {num} must be non-negative"));
        }
        Ok(Ratio { num, den })
    }

    /// Exact comparison by cross-multiplication, or `None` on overflow.
    ///
    /// Overflow is reported rather than wrapped: a silently wrapped comparison
    /// would put a wrong inequality into the evidence chain, which is the one
    /// failure this module exists to prevent.
    pub fn checked_cmp(self, other: Ratio) -> Option<std::cmp::Ordering> {
        let lhs = (self.num as i128).checked_mul(other.den as i128)?;
        let rhs = (other.num as i128).checked_mul(self.den as i128)?;
        Some(lhs.cmp(&rhs))
    }

    /// Exact `<`. Panics only on overflow, which `checked_cmp` makes detectable
    /// and which cannot occur for counts bounded by an alignment length.
    pub fn lt_exact(self, other: Ratio) -> bool {
        self.checked_cmp(other)
            .expect("count product overflowed i128")
            .is_lt()
    }

    /// Exact difference `self - other`, as a rational.
    pub fn sub_exact(self, other: Ratio) -> Option<Ratio> {
        let num = (self.num as i128)
            .checked_mul(other.den as i128)?
            .checked_sub((other.num as i128).checked_mul(self.den as i128)?)?;
        let den = (self.den as i128).checked_mul(other.den as i128)?;
        if num < 0 {
            return None;
        }
        Some(Ratio {
            num: i64::try_from(num).ok()?,
            den: i64::try_from(den).ok()?,
        })
    }

    /// For display and ranking only — never for a reported inequality.
    pub fn to_f64(self) -> f64 {
        self.num as f64 / self.den as f64
    }
}

/// Pairwise site counts for an alignment: the exact distance matrix.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExactDistances {
    /// Sequence labels.
    pub labels: Vec<String>,
    /// Row-major differing-site counts.
    pub diff: Vec<i64>,
    /// Row-major comparable-site counts.
    pub comparable: Vec<i64>,
    /// Alignment width.
    pub alignment_length: usize,
}

impl ExactDistances {
    /// Number of sequences.
    pub fn len(&self) -> usize {
        self.labels.len()
    }

    /// Whether there are no sequences.
    pub fn is_empty(&self) -> bool {
        self.labels.is_empty()
    }

    /// The exact distance of a pair.
    pub fn ratio(&self, i: usize, j: usize) -> Ratio {
        let n = self.labels.len();
        Ratio {
            num: self.diff[i * n + j],
            den: self.comparable[i * n + j],
        }
    }

    /// Whether a pair has any comparable sites at all.
    pub fn is_informative(&self, i: usize, j: usize) -> bool {
        let n = self.labels.len();
        self.comparable[i * n + j] > 0
    }

    /// Build from an alignment. Sites where either sequence is not `ACGT` are
    /// skipped for that pair, so denominators are per-pair.
    ///
    /// # Errors
    ///
    /// Fewer than three sequences, mismatched labels, or unequal lengths.
    pub fn from_alignment(labels: Vec<String>, sequences: &[Vec<u8>]) -> Result<Self, String> {
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
        let width = sequences[0].len();
        if let Some(bad) = sequences.iter().position(|s| s.len() != width) {
            return Err(format!(
                "sequence {bad} has length {} but the first has {width}; an \
                 alignment is required",
                sequences[bad].len()
            ));
        }

        validate_alphabet(sequences)?;

        let n = sequences.len();
        let mut diff = vec![0i64; n * n];
        let mut comparable = vec![0i64; n * n];
        for i in 0..n {
            for j in (i + 1)..n {
                let (d, c) = count_pair(&sequences[i], &sequences[j]);
                diff[i * n + j] = d;
                diff[j * n + i] = d;
                comparable[i * n + j] = c;
                comparable[j * n + i] = c;
            }
            comparable[i * n + i] = width as i64;
        }
        Ok(ExactDistances {
            labels,
            diff,
            comparable,
            alignment_length: width,
        })
    }

    /// Build directly from counts, for tests.
    ///
    /// # Errors
    ///
    /// Non-square input or a non-positive off-diagonal comparable count.
    pub fn from_counts(
        labels: Vec<String>,
        diff: Vec<i64>,
        comparable: Vec<i64>,
    ) -> Result<Self, String> {
        let n = labels.len();
        if diff.len() != n * n || comparable.len() != n * n {
            return Err("count vectors must be n*n".into());
        }
        for i in 0..n {
            for j in 0..n {
                if i != j && comparable[i * n + j] <= 0 {
                    return Err(format!("pair ({i},{j}) has no comparable sites"));
                }
                if diff[i * n + j] != diff[j * n + i]
                    || comparable[i * n + j] != comparable[j * n + i]
                {
                    return Err(format!("counts not symmetric at ({i},{j})"));
                }
            }
        }
        Ok(ExactDistances {
            labels,
            diff,
            comparable,
            alignment_length: 0,
        })
    }
}

/// Reject input that is not DNA, instead of silently analysing its A/C/G/T columns.
///
/// Two defects an audit found: an RNA alignment had every `U` dropped as
/// non-`ACGT`, and a protein alignment was analysed over whatever A/C/G/T
/// columns it happened to contain. Both produced confident output from input the
/// model does not apply to. Silence is the bug; this function makes it loud.
pub fn validate_alphabet(sequences: &[Vec<u8>]) -> Result<(), String> {
    let mut acgt = 0usize;
    let mut u = 0usize;
    let mut t = 0usize;
    let mut other = 0usize;
    let mut total = 0usize;
    for seq in sequences {
        for &b in seq {
            total += 1;
            match b.to_ascii_uppercase() {
                b'A' | b'C' | b'G' => acgt += 1,
                b'T' => {
                    acgt += 1;
                    t += 1;
                }
                b'U' => u += 1,
                b'N' | b'-' | b'.' | b'?' => {}
                _ => other += 1,
            }
        }
    }
    if total == 0 {
        return Err("empty input".into());
    }
    if u > 0 && t == 0 {
        return Err(format!(
            "input looks like RNA: {u} `U` and no `T`. This computes DNA \
             p-distances and would silently drop every `U` as non-comparable. \
             Convert U to T first, deliberately."
        ));
    }
    // Amino-acid alphabets are dominated by symbols outside ACGTN-.
    if other * 10 > total {
        return Err(format!(
            "{other} of {total} characters are outside ACGTUN-. This looks like \
             protein or an unsupported alphabet; analysing only its A/C/G/T \
             columns would be meaningless. Refusing rather than guessing."
        ));
    }
    if acgt * 2 < total {
        return Err(format!(
            "only {acgt} of {total} characters are definite DNA bases; the \
             alignment is more than half gaps or ambiguities, so pairwise \
             deletion would leave distances incomparable between pairs."
        ));
    }
    Ok(())
}

/// Differing and comparable sites for one pair.
fn count_pair(a: &[u8], b: &[u8]) -> (i64, i64) {
    let mut diff = 0i64;
    let mut comparable = 0i64;
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

/// A certificate carrying only integers: exactly the fields of
/// `ANSE.GeneticaReticulata.CountCertificate`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CountCertificate {
    /// Label of one endpoint.
    pub x: String,
    /// Label of the midpoint.
    pub y: String,
    /// Label of the other endpoint.
    pub z: String,
    /// Indices into the matrix.
    pub indices: [usize; 3],
    /// Differing sites, `x`–`y`.
    pub n_xy: i64,
    /// Comparable sites, `x`–`y`.
    pub d_xy: i64,
    /// Differing sites, `y`–`z`.
    pub n_yz: i64,
    /// Comparable sites, `y`–`z`.
    pub d_yz: i64,
    /// Differing sites, `x`–`z`.
    pub n_xz: i64,
    /// Comparable sites, `x`–`z`.
    pub d_xz: i64,
    /// Exact margin numerator.
    pub margin_num: i64,
    /// Exact margin denominator.
    pub margin_den: i64,
}

impl CountCertificate {
    /// Re-derive the two cross-multiplied inequalities the Lean structure
    /// requires. The search is not trusted; this is, and the kernel is more so.
    pub fn validate(&self) -> Result<(), String> {
        if self.d_xy <= 0 || self.d_yz <= 0 || self.d_xz <= 0 {
            return Err("a comparable-site count is not positive".into());
        }
        let far_xy =
            (self.n_xy as i128) * (self.d_xz as i128) < (self.n_xz as i128) * (self.d_xy as i128);
        let far_yz =
            (self.n_yz as i128) * (self.d_xz as i128) < (self.n_xz as i128) * (self.d_yz as i128);
        if !far_xy {
            return Err(format!(
                "far_xy fails: {} * {} is not < {} * {}",
                self.n_xy, self.d_xz, self.n_xz, self.d_xy
            ));
        }
        if !far_yz {
            return Err(format!(
                "far_yz fails: {} * {} is not < {} * {}",
                self.n_yz, self.d_xz, self.n_xz, self.d_yz
            ));
        }
        Ok(())
    }

    /// The exact margin as an `f64`, for ranking and display only.
    pub fn margin_f64(&self) -> f64 {
        self.margin_num as f64 / self.margin_den as f64
    }
}

/// Exhaustively search every triple for count certificates.
///
/// Exact throughout: every comparison is an integer cross-product. Triples with
/// an uninformative pair are skipped and counted, because "no comparable sites"
/// must not masquerade as "distance zero".
pub fn find_count_certificates(
    dist: &ExactDistances,
    top_k: usize,
) -> (Vec<CountCertificate>, usize, usize) {
    let n = dist.len();
    let mut found: Vec<CountCertificate> = Vec::new();
    let mut violating = 0usize;
    let mut skipped = 0usize;

    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                if !dist.is_informative(i, j)
                    || !dist.is_informative(i, k)
                    || !dist.is_informative(j, k)
                {
                    skipped += 1;
                    continue;
                }
                let (r_ij, r_ik, r_jk) = (dist.ratio(i, j), dist.ratio(i, k), dist.ratio(j, k));
                // The largest distance's endpoints are x and z; the third point
                // is the midpoint y.
                let (x, y, z) = if !r_ij.lt_exact(r_ik) && !r_ij.lt_exact(r_jk) {
                    (i, k, j)
                } else if !r_ik.lt_exact(r_ij) && !r_ik.lt_exact(r_jk) {
                    (i, j, k)
                } else {
                    (j, i, k)
                };
                let (r_xy, r_yz, r_xz) = (dist.ratio(x, y), dist.ratio(y, z), dist.ratio(x, z));
                // A violation requires the far distance to strictly exceed BOTH
                // near distances -- exactly the two Lean fields.
                if !(r_xy.lt_exact(r_xz) && r_yz.lt_exact(r_xz)) {
                    continue;
                }
                violating += 1;
                let near = if r_xy.lt_exact(r_yz) { r_yz } else { r_xy };
                let Some(margin) = r_xz.sub_exact(near) else {
                    continue; // overflow: refuse to report rather than guess
                };
                let cert = CountCertificate {
                    x: dist.labels[x].clone(),
                    y: dist.labels[y].clone(),
                    z: dist.labels[z].clone(),
                    indices: [x, y, z],
                    n_xy: r_xy.num,
                    d_xy: r_xy.den,
                    n_yz: r_yz.num,
                    d_yz: r_yz.den,
                    n_xz: r_xz.num,
                    d_xz: r_xz.den,
                    margin_num: margin.num,
                    margin_den: margin.den,
                };
                if cert.validate().is_ok() {
                    found.push(cert);
                }
            }
        }
    }

    // Rank by exact margin, descending, comparing rationals by cross-product.
    found.sort_by(|a, b| {
        let ra = Ratio {
            num: a.margin_num,
            den: a.margin_den,
        };
        let rb = Ratio {
            num: b.margin_num,
            den: b.margin_den,
        };
        rb.checked_cmp(ra).unwrap_or(std::cmp::Ordering::Equal)
    });
    found.truncate(top_k);
    (found, violating, skipped)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ratio_rejects_bad_denominators() {
        assert!(Ratio::new(1, 0).is_err());
        assert!(Ratio::new(1, -3).is_err());
        assert!(Ratio::new(-1, 3).is_err());
        assert!(Ratio::new(0, 3).is_ok());
    }

    #[test]
    fn cross_multiplied_comparison_is_exact_where_floats_are_not() {
        // 1/3 vs 1/3 written with different denominators: equal exactly.
        let a = Ratio { num: 1, den: 3 };
        let b = Ratio { num: 333_333_333, den: 999_999_999 };
        assert_eq!(a.checked_cmp(b), Some(std::cmp::Ordering::Equal));
        assert!(!a.lt_exact(b) && !b.lt_exact(a));
        // A float comparison of these is not reliably equal; the point of the
        // module is that we never ask it to be.
        let c = Ratio { num: 333_333_334, den: 999_999_999 };
        assert!(a.lt_exact(c), "1/3 < 333333334/999999999 exactly");
    }

    #[test]
    fn comparison_cannot_overflow_by_construction() {
        // An audit flagged the previous version of this test as unfalsifiable,
        // and it was right: the product of two i64 always fits in i128, so
        // `checked_cmp` can never return None. Stating that as the invariant is
        // honest; asserting `is_some()` and calling it an overflow test was not.
        let big = Ratio { num: i64::MAX, den: i64::MAX };
        assert_eq!(big.checked_cmp(big), Some(std::cmp::Ordering::Equal));
    }

    #[test]
    fn subtraction_reports_overflow_instead_of_wrapping() {
        // `sub_exact` CAN overflow, because it narrows back to i64. This is the
        // path where reporting rather than wrapping actually matters.
        let a = Ratio { num: i64::MAX, den: i64::MAX };
        let b = Ratio { num: 1, den: 2 };
        assert!(
            a.sub_exact(b).is_none(),
            "a result too large for i64 must be refused, not wrapped"
        );
        // And a negative difference is refused rather than returned unsigned.
        let small = Ratio { num: 1, den: 4 };
        let large = Ratio { num: 3, den: 4 };
        assert!(small.sub_exact(large).is_none());
    }

    #[test]
    fn subtraction_is_exact() {
        let a = Ratio { num: 131, den: 4000 };
        let b = Ratio { num: 70, den: 4000 };
        let d = a.sub_exact(b).unwrap();
        // 61/4000 up to the unreduced representation.
        assert_eq!(d.num as i128 * 4000, 61 * d.den as i128);
    }

    fn counts(labels: &[&str], diff: Vec<i64>, comp: Vec<i64>) -> ExactDistances {
        ExactDistances::from_counts(
            labels.iter().map(|s| s.to_string()).collect(),
            diff,
            comp,
        )
        .unwrap()
    }

    #[test]
    fn ultrametric_counts_yield_no_certificate() {
        // Two largest equal: the ultrametric condition on a triple.
        let d = counts(
            &["a", "b", "c"],
            vec![0, 1, 2, 1, 0, 2, 2, 2, 0],
            vec![4, 4, 4, 4, 4, 4, 4, 4, 4],
        );
        let (certs, violating, skipped) = find_count_certificates(&d, 10);
        assert!(certs.is_empty());
        assert_eq!(violating, 0);
        assert_eq!(skipped, 0);
    }

    #[test]
    fn open_path_counts_yield_a_validated_certificate() {
        let d = counts(
            &["a", "m", "b"],
            vec![0, 1, 3, 1, 0, 1, 3, 1, 0],
            vec![4, 4, 4, 4, 4, 4, 4, 4, 4],
        );
        let (certs, violating, _) = find_count_certificates(&d, 10);
        assert_eq!(violating, 1);
        assert_eq!(certs.len(), 1);
        let c = &certs[0];
        assert_eq!(c.y, "m");
        assert_eq!((c.n_xz, c.d_xz), (3, 4));
        assert!(c.validate().is_ok());
        // Exact margin 3/4 - 1/4 = 1/2.
        assert_eq!(c.margin_num as i128 * 2, c.margin_den as i128);
    }

    #[test]
    fn differing_denominators_are_handled_exactly() {
        // The gap case: pairs have different comparable-site counts. 2/10 and
        // 2/10 against 3/9 = 1/3. Violation, computed without a common
        // denominator anywhere.
        let d = counts(
            &["a", "m", "b"],
            vec![0, 2, 3, 2, 0, 2, 3, 2, 0],
            vec![9, 10, 9, 10, 9, 10, 9, 10, 9],
        );
        let (certs, violating, _) = find_count_certificates(&d, 10);
        assert_eq!(violating, 1, "3/9 exceeds both 2/10 values");
        assert!(certs[0].validate().is_ok());
        assert_ne!(certs[0].d_xy, certs[0].d_xz, "denominators really do differ");
    }

    #[test]
    fn uninformative_pairs_are_skipped_not_treated_as_zero() {
        // `a` and `b` occupy disjoint halves of the alignment, so that PAIR has
        // no comparable column, while the input as a whole is ordinary DNA and
        // passes alphabet validation. The triple must be skipped, not scored as
        // distance zero.
        let seqs = vec![
            b"ACGTA-----".to_vec(),
            b"-----ACGTA".to_vec(),
            b"ACGTAACGTA".to_vec(),
        ];
        let labels = vec!["a".into(), "b".into(), "c".into()];
        let d = ExactDistances::from_alignment(labels, &seqs).unwrap();
        assert_eq!(d.ratio(0, 1).den, 0, "a and b share no comparable column");
        assert!(!d.is_informative(0, 1));
        let (certs, _, skipped) = find_count_certificates(&d, 10);
        assert_eq!(skipped, 1, "the disjoint pair makes the triple unusable");
        assert!(certs.is_empty());
    }

    #[test]
    fn validate_rejects_forged_integers() {
        let mut c = CountCertificate {
            x: "a".into(),
            y: "m".into(),
            z: "b".into(),
            indices: [0, 1, 2],
            n_xy: 1,
            d_xy: 4,
            n_yz: 1,
            d_yz: 4,
            n_xz: 3,
            d_xz: 4,
            margin_num: 2,
            margin_den: 4,
        };
        assert!(c.validate().is_ok());
        // Forge the far distance downward: no longer a violation.
        c.n_xz = 1;
        assert!(c.validate().is_err());
        // Forge a denominator to zero.
        c.n_xz = 3;
        c.d_xy = 0;
        assert!(c.validate().is_err());
    }

    #[test]
    fn alignment_counts_match_a_hand_count() {
        let seqs = vec![
            b"ACGTACGT".to_vec(),
            b"ACGTACGA".to_vec(),
            b"TCGTACGT".to_vec(),
        ];
        let labels = vec!["a".into(), "b".into(), "c".into()];
        let d = ExactDistances::from_alignment(labels, &seqs).unwrap();
        assert_eq!(d.ratio(0, 1), Ratio { num: 1, den: 8 });
        assert_eq!(d.ratio(0, 2), Ratio { num: 1, den: 8 });
        assert_eq!(d.ratio(1, 2), Ratio { num: 2, den: 8 });
    }
}

#[cfg(test)]
mod alphabet_tests {
    use super::*;

    #[test]
    fn rna_is_refused_rather_than_silently_stripped() {
        let seqs = vec![
            b"ACGUACGU".to_vec(),
            b"ACGUACGA".to_vec(),
            b"UCGUACGU".to_vec(),
        ];
        let err = validate_alphabet(&seqs).expect_err("RNA must be refused");
        assert!(err.contains("RNA"), "{err}");
    }

    #[test]
    fn protein_is_refused() {
        let seqs = vec![
            b"MKLPQRSTVWY".to_vec(),
            b"MKLPQRSTVWF".to_vec(),
            b"MKLPQRSTVWH".to_vec(),
        ];
        let err = validate_alphabet(&seqs).expect_err("protein must be refused");
        assert!(err.contains("protein") || err.contains("alphabet"), "{err}");
    }

    #[test]
    fn mostly_gapped_input_is_refused() {
        let seqs = vec![
            b"A-------".to_vec(),
            b"A-------".to_vec(),
            b"C-------".to_vec(),
        ];
        assert!(validate_alphabet(&seqs).is_err());
    }

    #[test]
    fn ordinary_dna_with_some_gaps_is_accepted() {
        let seqs = vec![
            b"ACGTACGTN-".to_vec(),
            b"ACGTACGAN-".to_vec(),
            b"TCGTACGTN-".to_vec(),
        ];
        assert!(validate_alphabet(&seqs).is_ok());
    }

    #[test]
    fn dna_containing_both_t_and_u_is_allowed_through() {
        // Mixed T and U is more likely a sequencing artefact than RNA; we do
        // not refuse it, but the U positions are incomparable. Documented so the
        // behaviour is a choice rather than an accident.
        let seqs = vec![
            b"ACGTACGU".to_vec(),
            b"ACGTACGA".to_vec(),
            b"TCGTACGT".to_vec(),
        ];
        assert!(validate_alphabet(&seqs).is_ok());
    }

    #[test]
    fn validation_is_wired_into_the_constructor() {
        // The check must actually run, not merely exist.
        let seqs = vec![
            b"ACGUACGU".to_vec(),
            b"ACGUACGA".to_vec(),
            b"UCGUACGU".to_vec(),
        ];
        let labels = vec!["a".into(), "b".into(), "c".into()];
        assert!(ExactDistances::from_alignment(labels, &seqs).is_err());
    }
}

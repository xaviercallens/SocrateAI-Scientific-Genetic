//! Complete deletion within the subset, and the four-point (tree-level) test.
//!
//! Two defects of the pairwise-deletion ultrametric certificate motivate this
//! module, and an adversarial audit of the Phase 1 prototype named both.
//!
//! # Complete deletion within the subset
//!
//! Under pairwise deletion, two pairs of sequences are compared over *different
//! sets of columns*. Distances computed over different columns are not
//! commensurable, and their ratio can violate the ultrametric inequality even
//! when the genealogy is exactly a clock-like tree. The kernel cannot detect
//! this: the inequality it is handed is true. Finding 6 of the audit.
//!
//! The fix is to restrict a certificate's counting to the columns at which
//! *every member of the subset* carries a definite base. All distances in the
//! certificate then share one denominator, so
//!
//! * the certificate's inequalities become comparisons between integers, with
//!   no cross-multiplication anywhere, and
//! * the confound cannot arise, because there is only one set of columns.
//!
//! `ANSE.GeneticaReticulata.pairwise_deletion_manufactures_what_complete_deletion_denies`
//! is the kernel's own witness that the two policies really disagree. The price
//! is power: the common denominator is smaller than any pairwise one, so margins
//! are coarser. That is the right trade for data with gaps.
//!
//! # The four-point condition
//!
//! `IsUltrametric` encodes "clonal genealogy AND strict molecular clock", so a
//! violation is ambiguous between reticulation and rate variation. The four-point
//! condition encodes "tree" *without* "clock": branch lengths may vary freely.
//! For a quadruple, the three ways of pairing it give three sums of distances,
//! and a tree metric always leaves the two largest equal. A quadruple whose
//! largest sum strictly exceeds both others therefore refutes treeness however
//! the rates vary along the tree.
//!
//! `fourPoint_of_ultrametric` proves the two nulls are nested, so a quartet
//! certificate refutes strictly more than a triple certificate does. That
//! theorem is what makes the experimental design coherent: a dataset with large
//! triple margins and floor-level quartet margins is evidence of rate variation
//! on a tree, not of reticulation.
//!
//! That a metric satisfying the four-point condition on every quadruple is
//! realised by a weighted tree is **Buneman's theorem, cited not proved** — the
//! same status "a strict clock gives an ultrametric" has in [`crate::exact`].

use serde::{Deserialize, Serialize};

use crate::exact::{validate_alphabet, Ratio};

/// An alignment: labelled sequences of equal length, validated as DNA.
///
/// Carries the sequences rather than a distance matrix, because complete
/// deletion within a subset cannot be read off a precomputed matrix — the
/// denominator depends on which subset is being examined.
#[derive(Debug, Clone)]
pub struct Alignment {
    /// Sequence labels.
    pub labels: Vec<String>,
    /// Aligned sequences, all of the same length.
    pub sequences: Vec<Vec<u8>>,
}

/// Counts for one triple under complete deletion within that triple.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct TripleCounts {
    /// Columns at which all three carry a definite base: the one denominator.
    pub common: i64,
    /// Differing sites, first–second.
    pub n_xy: i64,
    /// Differing sites, second–third.
    pub n_yz: i64,
    /// Differing sites, first–third.
    pub n_xz: i64,
}

/// Counts for one quadruple under complete deletion within that quadruple.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct QuadrupleCounts {
    /// Columns at which all four carry a definite base: the one denominator.
    pub common: i64,
    /// Differing sites for the pair `(1,2)`.
    pub n_ab: i64,
    /// Differing sites for the pair `(3,4)`.
    pub n_ce: i64,
    /// Differing sites for the pair `(1,3)`.
    pub n_ac: i64,
    /// Differing sites for the pair `(2,4)`.
    pub n_be: i64,
    /// Differing sites for the pair `(1,4)`.
    pub n_ae: i64,
    /// Differing sites for the pair `(2,3)`.
    pub n_bc: i64,
}

impl QuadrupleCounts {
    /// The three pairing sums, in the order `(ab|ce)`, `(ac|be)`, `(ae|bc)`.
    ///
    /// Sums, not averages: with one denominator the four-point comparison is an
    /// integer comparison between these three numbers.
    pub fn pairing_sums(&self) -> [i64; 3] {
        [
            self.n_ab + self.n_ce,
            self.n_ac + self.n_be,
            self.n_ae + self.n_bc,
        ]
    }
}

fn is_base(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

impl Alignment {
    /// Build an alignment, validating length agreement and the alphabet.
    ///
    /// # Errors
    ///
    /// Fewer than three sequences, a label/sequence count mismatch, unequal
    /// lengths, or input that is not DNA (RNA and protein are refused rather
    /// than silently analysed; see [`validate_alphabet`]).
    pub fn new(labels: Vec<String>, sequences: Vec<Vec<u8>>) -> Result<Self, String> {
        if sequences.len() < 3 {
            return Err(format!(
                "need at least 3 sequences, got {}",
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
        validate_alphabet(&sequences)?;
        Ok(Alignment { labels, sequences })
    }

    /// Number of sequences.
    pub fn len(&self) -> usize {
        self.labels.len()
    }

    /// Whether the alignment is empty.
    pub fn is_empty(&self) -> bool {
        self.labels.is_empty()
    }

    /// Alignment width.
    pub fn width(&self) -> usize {
        self.sequences.first().map_or(0, Vec::len)
    }

    /// Counts for a triple, over the columns where all three carry a base.
    pub fn triple_counts(&self, i: usize, j: usize, k: usize) -> TripleCounts {
        let (u, v, w) = (&self.sequences[i], &self.sequences[j], &self.sequences[k]);
        let mut c = TripleCounts {
            common: 0,
            n_xy: 0,
            n_yz: 0,
            n_xz: 0,
        };
        for p in 0..u.len() {
            let (a, b, d) = (u[p], v[p], w[p]);
            if !is_base(a) || !is_base(b) || !is_base(d) {
                continue;
            }
            let (a, b, d) = (
                a.to_ascii_uppercase(),
                b.to_ascii_uppercase(),
                d.to_ascii_uppercase(),
            );
            c.common += 1;
            if a != b {
                c.n_xy += 1;
            }
            if b != d {
                c.n_yz += 1;
            }
            if a != d {
                c.n_xz += 1;
            }
        }
        c
    }

    /// Counts for a quadruple, over the columns where all four carry a base.
    pub fn quadruple_counts(&self, i: usize, j: usize, k: usize, l: usize) -> QuadrupleCounts {
        let (u, v, w, t) = (
            &self.sequences[i],
            &self.sequences[j],
            &self.sequences[k],
            &self.sequences[l],
        );
        let mut c = QuadrupleCounts {
            common: 0,
            n_ab: 0,
            n_ce: 0,
            n_ac: 0,
            n_be: 0,
            n_ae: 0,
            n_bc: 0,
        };
        for p in 0..u.len() {
            if !is_base(u[p]) || !is_base(v[p]) || !is_base(w[p]) || !is_base(t[p]) {
                continue;
            }
            let a = u[p].to_ascii_uppercase();
            let b = v[p].to_ascii_uppercase();
            let cc = w[p].to_ascii_uppercase();
            let e = t[p].to_ascii_uppercase();
            c.common += 1;
            if a != b {
                c.n_ab += 1;
            }
            if cc != e {
                c.n_ce += 1;
            }
            if a != cc {
                c.n_ac += 1;
            }
            if b != e {
                c.n_be += 1;
            }
            if a != e {
                c.n_ae += 1;
            }
            if b != cc {
                c.n_bc += 1;
            }
        }
        c
    }
}

/// A triple certificate under complete deletion: one denominator, so the
/// obligations are comparisons between integers.
///
/// Mirrors the arguments of
/// `ANSE.GeneticaReticulata.complete_triple_refute_from_counts`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CompleteTripleCertificate {
    /// One endpoint.
    pub x: String,
    /// The midpoint, close to both endpoints.
    pub y: String,
    /// The other endpoint.
    pub z: String,
    /// Indices into the alignment, in the order `x, y, z`.
    pub indices: [usize; 3],
    /// The triple's common comparable-site count.
    pub common: i64,
    /// Differing sites, `x`–`y`.
    pub n_xy: i64,
    /// Differing sites, `y`–`z`.
    pub n_yz: i64,
    /// Differing sites, `x`–`z`.
    pub n_xz: i64,
    /// Exact margin numerator: `n_xz - max(n_xy, n_yz)`.
    pub margin_num: i64,
    /// Exact margin denominator: the common count.
    pub margin_den: i64,
}

impl CompleteTripleCertificate {
    /// Re-derive the obligations. No cross-products appear: one denominator
    /// means the rational comparison *is* the integer comparison, which is the
    /// content of `ratioN_lt_ratioN_iff`.
    pub fn validate(&self) -> Result<(), String> {
        if self.common <= 0 {
            return Err(format!(
                "the triple has {} comparable columns in common; with none, no \
                 complete-deletion distance is defined",
                self.common
            ));
        }
        if self.n_xy >= self.n_xz {
            return Err(format!(
                "far_xy fails: {} is not < {}",
                self.n_xy, self.n_xz
            ));
        }
        if self.n_yz >= self.n_xz {
            return Err(format!(
                "far_yz fails: {} is not < {}",
                self.n_yz, self.n_xz
            ));
        }
        let expected = self.n_xz - self.n_xy.max(self.n_yz);
        if self.margin_num != expected || self.margin_den != self.common {
            return Err(format!(
                "margin {}/{} is not the exact {expected}/{}",
                self.margin_num, self.margin_den, self.common
            ));
        }
        Ok(())
    }

    /// The exact margin as a ratio.
    pub fn margin(&self) -> Ratio {
        Ratio {
            num: self.margin_num,
            den: self.margin_den,
        }
    }
}

/// A quartet certificate: the four-point condition refuted on one quadruple.
///
/// The labels are ordered so that `(a,b)` against `(c,e)` is the **dominant**
/// pairing — the one whose distance sum strictly exceeds both others. Mirrors
/// `ANSE.GeneticaReticulata.quartet_refute_from_counts`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct QuartetCertificate {
    /// First member of the dominant pairing.
    pub a: String,
    /// Its partner: `a` and `b` are the far pair.
    pub b: String,
    /// First member of the other side of the dominant pairing.
    pub c: String,
    /// Its partner.
    pub e: String,
    /// Indices into the alignment, in the order `a, b, c, e`.
    pub indices: [usize; 4],
    /// The quadruple's common comparable-site count: the one denominator.
    pub common: i64,
    /// Differing sites for each of the six pairs, in the Lean field order.
    pub n_ab: i64,
    /// Differing sites, `c`–`e`.
    pub n_ce: i64,
    /// Differing sites, `a`–`c`.
    pub n_ac: i64,
    /// Differing sites, `b`–`e`.
    pub n_be: i64,
    /// Differing sites, `a`–`e`.
    pub n_ae: i64,
    /// Differing sites, `b`–`c`.
    pub n_bc: i64,
    /// Exact margin numerator: dominant sum minus the second largest.
    pub margin_num: i64,
    /// Exact margin denominator: the common count.
    pub margin_den: i64,
}

impl QuartetCertificate {
    /// Re-derive `fourPointViolates` and the margin, in integers.
    pub fn validate(&self) -> Result<(), String> {
        if self.common <= 0 {
            return Err(format!(
                "the quadruple has {} columns comparable in all four sequences",
                self.common
            ));
        }
        let s_far = self.n_ab + self.n_ce;
        let s_b = self.n_ac + self.n_be;
        let s_c = self.n_ae + self.n_bc;
        if s_b >= s_far {
            return Err(format!(
                "fourPointViolates fails: pairing (ac|be) sums to {s_b}, which is \
                 not < the claimed dominant {s_far}"
            ));
        }
        if s_c >= s_far {
            return Err(format!(
                "fourPointViolates fails: pairing (ae|bc) sums to {s_c}, which is \
                 not < the claimed dominant {s_far}"
            ));
        }
        let expected = s_far - s_b.max(s_c);
        if self.margin_num != expected || self.margin_den != self.common {
            return Err(format!(
                "margin {}/{} is not the exact {expected}/{}",
                self.margin_num, self.margin_den, self.common
            ));
        }
        Ok(())
    }

    /// The exact margin as a ratio.
    pub fn margin(&self) -> Ratio {
        Ratio {
            num: self.margin_num,
            den: self.margin_den,
        }
    }
}

/// Summaries of the *whole* quartet distribution, not only its maximum.
///
/// # Why this exists
///
/// The certificate this module is built around is a maximum over quadruples:
/// the single most-violating quadruple, exhibited as a witness the Lean kernel
/// can check. That is what makes a certificate possible, and it is also an
/// extreme-value statistic, whose null level grows with the number of subsets
/// examined. The Phase 2 grid shows it: the quartet floor moves from 5–7 to
/// 9–10 between 16 and 32 taxa. So margins are not comparable across
/// populations of different size, which is the incommensurability trap this
/// project keeps rediscovering — here latent in the primary statistic itself.
///
/// The treelikeness literature took the other branch. Its statistics are
/// *averages* over quadruples, and in the one published comparison of
/// treelikeness statistics used as hypothesis tests, the averaged variants
/// (mean δq, tree proportion) held a correct false-positive rate while the
/// mode-based and proportion-of-resolved-quartets variants failed as tests.
/// This structure computes the averaged summaries so the two can be compared on
/// the same data, which is a comparison nobody appears to have run.
///
/// # What these numbers are not
///
/// They are **floating point, and they are not certificates**. Every verdict in
/// [`CompleteTripleCertificate`] and [`QuartetCertificate`] is decided by
/// integer comparison and checkable by the kernel; nothing here is. A mean over
/// quadruples has no single witness to exhibit, which is exactly the trade: the
/// averaged statistic is better behaved as a test and cannot be certified, and
/// the maximum can be certified and is worse behaved. Reporting both is the
/// only honest option until the trade has been measured.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QuartetDistribution {
    /// Quadruples examined.
    pub examined: usize,
    /// Quadruples skipped for having no column comparable in all four.
    pub skipped_no_common_columns: usize,
    /// Quadruples whose dominant pairing strictly exceeds both others.
    pub violating: usize,
    /// `violating / examined`. The statistic reported to *fail* as a test in the
    /// published comparison, included because it is cheap and because its
    /// failure is itself informative.
    pub proportion_violating: f64,
    /// Mean of `δq = (m1 − m2) / (m1 − m3)` over the examined quadruples, where
    /// `m1 ≥ m2 ≥ m3` are the three pairing sums. Scale-invariant, zero exactly
    /// when the quadruple is additive. The averaged statistic.
    pub mean_delta_q: f64,
    /// The largest δq seen. Comparable to the certificate's own statistic in
    /// being a maximum, but normalised differently.
    pub max_delta_q: f64,
    /// Mean of `((m1 − m2) / D)²` over quadruples, divided by the squared mean
    /// pairwise p-distance: the Q-residual. Scale-*dependent* before that
    /// division, which is why the division is part of the definition.
    pub q_residual: f64,
    /// Mean pairwise p-distance under pairwise deletion, the Q-residual's
    /// normaliser, reported so the division can be checked.
    pub mean_p_distance: f64,
}

/// Compute the averaged quartet statistics alongside the certificate's maximum.
///
/// Exhaustive over quadruples, with the same refusal-to-sample ceiling as
/// [`find_quartet_certificates`]: a subsampled combination set would change what
/// was examined without saying so.
///
/// # Errors
///
/// Fewer than four sequences, a quadruple count above `max_quadruples`, or a
/// malformed alignment.
pub fn quartet_distribution(
    aln: &Alignment,
    max_quadruples: usize,
) -> Result<QuartetDistribution, String> {
    let n = aln.len();
    if n < 4 {
        return Err(format!(
            "need at least 4 sequences for a four-point test, got {n}"
        ));
    }
    let total = n * (n - 1) * (n - 2) * (n - 3) / 24;
    if total > max_quadruples {
        return Err(format!(
            "{n} sequences give {total} quadruples, above the ceiling of              {max_quadruples}. Raise it deliberately: this search is exhaustive              and will not sample."
        ));
    }

    // The Q-residual's normaliser, under pairwise deletion over the whole
    // alignment rather than per quadruple, following the published definition.
    let distances =
        crate::exact::ExactDistances::from_alignment(aln.labels.clone(), &aln.sequences)?;
    let mut distance_sum = 0.0f64;
    let mut pairs = 0.0f64;
    for i in 0..n {
        for j in (i + 1)..n {
            distance_sum += distances.ratio(i, j).to_f64();
            pairs += 1.0;
        }
    }
    let mean_p_distance = if pairs > 0.0 {
        distance_sum / pairs
    } else {
        0.0
    };

    let mut examined = 0usize;
    let mut skipped = 0usize;
    let mut violating = 0usize;
    let mut delta_sum = 0.0f64;
    let mut delta_max = 0.0f64;
    let mut squared_sum = 0.0f64;

    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                for l in (k + 1)..n {
                    let counts = aln.quadruple_counts(i, j, k, l);
                    if counts.common == 0 {
                        skipped += 1;
                        continue;
                    }
                    examined += 1;
                    let mut sums = counts.pairing_sums();
                    sums.sort_unstable_by(|a, b| b.cmp(a));
                    let (m1, m2, m3) = (sums[0], sums[1], sums[2]);
                    if m1 > m2 {
                        violating += 1;
                    }
                    // δq is 0 on an additive quadruple, where m1 == m2; the
                    // denominator vanishes only when all three tie, and then
                    // the quadruple is additive too.
                    let delta = if m1 > m3 {
                        (m1 - m2) as f64 / (m1 - m3) as f64
                    } else {
                        0.0
                    };
                    delta_sum += delta;
                    if delta > delta_max {
                        delta_max = delta;
                    }
                    let gap = (m1 - m2) as f64 / counts.common as f64;
                    squared_sum += gap * gap;
                }
            }
        }
    }

    let denominator = examined.max(1) as f64;
    let q_residual = if mean_p_distance > 0.0 {
        (squared_sum / denominator) / (mean_p_distance * mean_p_distance)
    } else {
        0.0
    };
    Ok(QuartetDistribution {
        examined,
        skipped_no_common_columns: skipped,
        violating,
        proportion_violating: violating as f64 / denominator,
        mean_delta_q: delta_sum / denominator,
        max_delta_q: delta_max,
        q_residual,
        mean_p_distance,
    })
}

/// Rank certificates by exact margin, descending, and keep the top `k`.
fn rank_by_margin<T, F: Fn(&T) -> Ratio>(mut v: Vec<T>, margin: F, top_k: usize) -> Vec<T> {
    v.sort_by(|p, q| {
        margin(q)
            .checked_cmp(margin(p))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    v.truncate(top_k);
    v
}

/// Outcome of a complete-deletion triple search.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TripleSearch {
    /// Top certificates by exact margin.
    pub certificates: Vec<CompleteTripleCertificate>,
    /// Triples examined.
    pub examined: usize,
    /// Triples that violated the strong triangle inequality.
    pub violating: usize,
    /// Triples skipped because the three sequences share no comparable column.
    pub skipped_no_common_columns: usize,
    /// Every emitted certificate re-validated.
    pub all_validated: bool,
}

/// Outcome of a quartet search.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QuartetSearch {
    /// Top certificates by exact margin.
    pub certificates: Vec<QuartetCertificate>,
    /// Quadruples examined.
    pub examined: usize,
    /// Quadruples whose dominant pairing strictly exceeded both others.
    pub violating: usize,
    /// Quadruples skipped because the four share no comparable column.
    pub skipped_no_common_columns: usize,
    /// Every emitted certificate re-validated.
    pub all_validated: bool,
}

/// Search every triple under complete deletion within the triple.
///
/// Exhaustive and exact. The midpoint `y` is the sequence not in the most
/// distant pair; a violation requires the far distance to strictly exceed both
/// near ones, which is exactly the pair of Lean obligations.
pub fn find_complete_triple_certificates(
    aln: &Alignment,
    top_k: usize,
) -> Result<TripleSearch, String> {
    let n = aln.len();
    let mut found = Vec::new();
    let mut examined = 0usize;
    let mut violating = 0usize;
    let mut skipped = 0usize;

    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                examined += 1;
                let c = aln.triple_counts(i, j, k);
                if c.common == 0 {
                    skipped += 1;
                    continue;
                }
                // Order so that (x, z) is the most distant pair and y the
                // midpoint. Counts share a denominator, so this is a comparison
                // between integers.
                let (xi, yi, zi, n_xy, n_yz, n_xz) = if c.n_xy >= c.n_yz && c.n_xy >= c.n_xz {
                    // i--j is the far pair; k is the midpoint.
                    (i, k, j, c.n_xz, c.n_yz, c.n_xy)
                } else if c.n_xz >= c.n_xy && c.n_xz >= c.n_yz {
                    // i--k is the far pair; j is the midpoint.
                    (i, j, k, c.n_xy, c.n_yz, c.n_xz)
                } else {
                    // j--k is the far pair; i is the midpoint.
                    (j, i, k, c.n_xy, c.n_xz, c.n_yz)
                };
                if n_xy >= n_xz || n_yz >= n_xz {
                    continue;
                }
                violating += 1;
                let cert = CompleteTripleCertificate {
                    x: aln.labels[xi].clone(),
                    y: aln.labels[yi].clone(),
                    z: aln.labels[zi].clone(),
                    indices: [xi, yi, zi],
                    common: c.common,
                    n_xy,
                    n_yz,
                    n_xz,
                    margin_num: n_xz - n_xy.max(n_yz),
                    margin_den: c.common,
                };
                cert.validate()
                    .map_err(|e| format!("search produced an invalid certificate: {e}"))?;
                found.push(cert);
            }
        }
    }

    let certificates = rank_by_margin(found, CompleteTripleCertificate::margin, top_k);
    let all_validated = certificates.iter().all(|c| c.validate().is_ok());
    Ok(TripleSearch {
        certificates,
        examined,
        violating,
        skipped_no_common_columns: skipped,
        all_validated,
    })
}

/// Search every quadruple for a four-point violation.
///
/// Exhaustive, with an explicit ceiling: above `max_quadruples` the search
/// **refuses** rather than sampling. Stride sampling over combinations is how
/// the textbook's Section 11.3 produced a false positive, so a cap that silently
/// changes what was examined is not an option here.
///
/// # Errors
///
/// If the number of quadruples exceeds `max_quadruples`.
pub fn find_quartet_certificates(
    aln: &Alignment,
    top_k: usize,
    max_quadruples: usize,
) -> Result<QuartetSearch, String> {
    let n = aln.len();
    if n < 4 {
        return Err(format!(
            "need at least 4 sequences for a four-point test, got {n}"
        ));
    }
    let total = n * (n - 1) * (n - 2) * (n - 3) / 24;
    if total > max_quadruples {
        return Err(format!(
            "{n} sequences give {total} quadruples, above the ceiling of \
             {max_quadruples}. Raise it deliberately: this search is exhaustive \
             and will not sample, because a subsampled combination set would \
             change what was examined without saying so."
        ));
    }

    let mut found = Vec::new();
    let mut examined = 0usize;
    let mut violating = 0usize;
    let mut skipped = 0usize;

    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                for l in (k + 1)..n {
                    examined += 1;
                    let c = aln.quadruple_counts(i, j, k, l);
                    if c.common == 0 {
                        skipped += 1;
                        continue;
                    }
                    let sums = c.pairing_sums();
                    // Which pairing dominates? Ties mean no violation.
                    let (best, _) = sums
                        .iter()
                        .enumerate()
                        .max_by_key(|(_, s)| **s)
                        .expect("three sums");
                    let others: Vec<i64> = (0..3).filter(|t| *t != best).map(|t| sums[t]).collect();
                    if others.iter().any(|s| *s >= sums[best]) {
                        continue;
                    }
                    violating += 1;
                    // Relabel so the dominant pairing is (a,b)|(c,e), then
                    // recount in that order: the certificate's field names must
                    // mean what Lean's field names mean.
                    let order = match best {
                        0 => [i, j, k, l],
                        1 => [i, k, j, l],
                        _ => [i, l, j, k],
                    };
                    let q = aln.quadruple_counts(order[0], order[1], order[2], order[3]);
                    let s = q.pairing_sums();
                    let cert = QuartetCertificate {
                        a: aln.labels[order[0]].clone(),
                        b: aln.labels[order[1]].clone(),
                        c: aln.labels[order[2]].clone(),
                        e: aln.labels[order[3]].clone(),
                        indices: order,
                        common: q.common,
                        n_ab: q.n_ab,
                        n_ce: q.n_ce,
                        n_ac: q.n_ac,
                        n_be: q.n_be,
                        n_ae: q.n_ae,
                        n_bc: q.n_bc,
                        margin_num: s[0] - s[1].max(s[2]),
                        margin_den: q.common,
                    };
                    cert.validate()
                        .map_err(|e| format!("search produced an invalid quartet: {e}"))?;
                    found.push(cert);
                }
            }
        }
    }

    let certificates = rank_by_margin(found, QuartetCertificate::margin, top_k);
    let all_validated = certificates.iter().all(|c| c.validate().is_ok());
    Ok(QuartetSearch {
        certificates,
        examined,
        violating,
        skipped_no_common_columns: skipped,
        all_validated,
    })
}

/// Per-triple comparison of the two deletion policies.
///
/// The statistic that measures the audit's finding 6. Comparing *margin
/// magnitudes* between the policies does not, and a first version of the
/// experiment made that mistake: complete deletion shrinks the denominator, so a
/// one-site difference is worth more and its margins are *larger* even where it
/// is strictly more conservative about what counts as a violation. The two
/// numbers are not commensurable, for the same reason two pairwise distances over
/// different columns are not.
///
/// What is commensurable is the per-triple verdict. `pairwise_only` counts
/// triples that pairwise deletion calls a violation and complete deletion does
/// not: certificates manufactured by the incomparability of denominators.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolicyComparison {
    /// Triples examined.
    pub triples: usize,
    /// Triples both policies call a violation.
    pub both: usize,
    /// **Manufactured:** a violation under pairwise deletion only.
    pub pairwise_only: usize,
    /// A violation under complete deletion only. Possible in principle — the
    /// common columns are a different sample — and reported rather than assumed
    /// away.
    pub complete_only: usize,
    /// Triples neither policy calls a violation.
    pub neither: usize,
    /// Triples with no column comparable in all three.
    pub skipped_no_common_columns: usize,
    /// One manufactured triple, named, so a reader can inspect it.
    pub example_pairwise_only: Option<[String; 3]>,
}

/// Compare the deletion policies triple by triple.
///
/// On an alignment with no missing data the two policies coincide exactly, so
/// `pairwise_only` and `complete_only` must both be zero — which is the control
/// that gives a non-zero count on gappy data its meaning.
///
/// # Errors
///
/// If the pairwise distance matrix cannot be built from the alignment.
pub fn compare_deletion_policies(aln: &Alignment) -> Result<PolicyComparison, String> {
    let pw = crate::exact::ExactDistances::from_alignment(aln.labels.clone(), &aln.sequences)?;
    let n = aln.len();
    let mut out = PolicyComparison {
        triples: 0,
        both: 0,
        pairwise_only: 0,
        complete_only: 0,
        neither: 0,
        skipped_no_common_columns: 0,
        example_pairwise_only: None,
    };

    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                out.triples += 1;
                let c = aln.triple_counts(i, j, k);
                if c.common == 0 {
                    out.skipped_no_common_columns += 1;
                    continue;
                }
                // Complete deletion: one denominator, so integers decide.
                let m = c.n_xy.max(c.n_yz).max(c.n_xz);
                let complete_violates = [c.n_xy, c.n_yz, c.n_xz]
                    .iter()
                    .filter(|v| **v == m)
                    .count()
                    == 1;

                // Pairwise deletion: three denominators, so cross-products do.
                let pairwise_violates = if pw.is_informative(i, j)
                    && pw.is_informative(i, k)
                    && pw.is_informative(j, k)
                {
                    let r = [pw.ratio(i, j), pw.ratio(j, k), pw.ratio(i, k)];
                    (0..3).any(|t| {
                        (0..3)
                            .filter(|s| *s != t)
                            .all(|s| r[s].lt_exact(r[t]))
                    })
                } else {
                    false
                };

                match (pairwise_violates, complete_violates) {
                    (true, true) => out.both += 1,
                    (true, false) => {
                        out.pairwise_only += 1;
                        if out.example_pairwise_only.is_none() {
                            out.example_pairwise_only = Some([
                                aln.labels[i].clone(),
                                aln.labels[j].clone(),
                                aln.labels[k].clone(),
                            ]);
                        }
                    }
                    (false, true) => out.complete_only += 1,
                    (false, false) => out.neither += 1,
                }
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln(seqs: &[&[u8]]) -> Alignment {
        let labels: Vec<String> = (0..seqs.len()).map(|i| format!("s{i}")).collect();
        Alignment::new(labels, seqs.iter().map(|s| s.to_vec()).collect()).unwrap()
    }

    #[test]
    fn complete_deletion_drops_a_column_gapped_in_any_member() {
        // Position 4 is gapped in the second sequence, so it leaves the triple's
        // common set entirely -- including the first--third comparison, which
        // pairwise deletion would have kept.
        let a = aln(&[b"AAAA", b"AAC-", b"AACC"]);
        let c = a.triple_counts(0, 1, 2);
        assert_eq!(c.common, 3, "the gapped column must leave all three counts");
        assert_eq!(c.n_xz, 1, "A vs C at position 3 only");
        // Pairwise, the first--third pair would have had 4 columns and 2 diffs.
        let pairwise_diff = b"AAAA"
            .iter()
            .zip(b"AACC")
            .filter(|(x, y)| x != y)
            .count();
        assert_eq!(pairwise_diff, 2, "pairwise keeps the column, complete does not");
    }

    #[test]
    fn complete_deletion_denies_a_certificate_pairwise_deletion_grants() {
        // The same three sequences as the Lean witness
        // `pairwise_deletion_manufactures_what_complete_deletion_denies`:
        // gapping the midpoint at two columns inflates only the far distance.
        let seqs: Vec<Vec<u8>> = vec![
            b"AAAAAAAA".to_vec(),
            b"CCAAAA--".to_vec(),
            b"CACAAACC".to_vec(),
        ];
        let labels = vec!["x".into(), "y".into(), "z".into()];

        // Pairwise deletion: a valid certificate.
        let pw = crate::exact::ExactDistances::from_alignment(labels.clone(), &seqs).unwrap();
        let (pw_certs, pw_violating, _) = crate::exact::find_count_certificates(&pw, 5);
        assert_eq!(pw_violating, 1, "2/6, 2/6 against 4/8 is a violation");
        assert!(pw_certs[0].validate().is_ok());

        // Complete deletion within the triple: no certificate at all.
        let a = Alignment::new(labels, seqs).unwrap();
        let counts = a.triple_counts(0, 1, 2);
        assert_eq!(
            (counts.common, counts.n_xy, counts.n_yz, counts.n_xz),
            (6, 2, 2, 2),
            "all three distances are 2/6 once the gapped columns are dropped"
        );
        let search = find_complete_triple_certificates(&a, 5).unwrap();
        assert_eq!(
            search.violating, 0,
            "the violation was an artefact of the differing denominators"
        );
        assert!(search.certificates.is_empty());
    }

    #[test]
    fn the_policies_coincide_exactly_when_nothing_is_missing() {
        // The control that gives a non-zero manufactured count its meaning: with
        // no gaps, complete deletion within the triple IS pairwise deletion, so
        // the two policies must agree on every triple. If this can fail, the
        // comparison statistic is measuring something else.
        let pop = crate::simulate::clonal_tree(600, 3, 0.02, 4242);
        let a = Alignment::new(pop.labels, pop.sequences).unwrap();
        let cmp = compare_deletion_policies(&a).unwrap();
        assert_eq!(cmp.triples, 56);
        assert_eq!(cmp.pairwise_only, 0, "no gaps, so nothing can be manufactured");
        assert_eq!(cmp.complete_only, 0);
        assert!(cmp.both > 0, "clonal data still violates by sampling noise");
        assert_eq!(cmp.both + cmp.neither, cmp.triples);
    }

    #[test]
    fn missing_data_manufactures_a_violation_the_comparison_detects() {
        let seqs: Vec<Vec<u8>> = vec![
            b"AAAAAAAA".to_vec(),
            b"CCAAAA--".to_vec(),
            b"CACAAACC".to_vec(),
        ];
        let a = Alignment::new(vec!["x".into(), "y".into(), "z".into()], seqs).unwrap();
        let cmp = compare_deletion_policies(&a).unwrap();
        assert_eq!(cmp.triples, 1);
        assert_eq!(cmp.pairwise_only, 1, "the certificate is an artefact of the gaps");
        assert_eq!(cmp.both, 0);
        assert_eq!(cmp.complete_only, 0);
        assert_eq!(
            cmp.example_pairwise_only,
            Some(["x".into(), "y".into(), "z".into()])
        );
    }

    #[test]
    fn an_ultrametric_triple_yields_no_complete_certificate() {
        // Three equidistant sequences: the two largest distances are equal, so
        // the strong triangle inequality holds with equality.
        let a = aln(&[b"AAAA", b"CCAA", b"CAAC"]);
        let c = a.triple_counts(0, 1, 2);
        assert_eq!((c.n_xy, c.n_yz, c.n_xz), (2, 2, 2));
        let s = find_complete_triple_certificates(&a, 5).unwrap();
        assert_eq!(s.violating, 0);
        assert_eq!(s.examined, 1);
    }

    #[test]
    fn an_open_path_yields_a_complete_certificate_with_the_right_midpoint() {
        let a = aln(&[b"AAAA", b"AACC", b"CCCC"]);
        let s = find_complete_triple_certificates(&a, 5).unwrap();
        assert_eq!(s.violating, 1);
        let cert = &s.certificates[0];
        assert_eq!(cert.y, "s1", "the mosaic is the midpoint");
        assert_eq!((cert.n_xy, cert.n_yz, cert.n_xz), (2, 2, 4));
        // Exact margin 4/4 - 2/4 = 2/4.
        assert_eq!((cert.margin_num, cert.margin_den), (2, 4));
        assert!(cert.validate().is_ok());
    }

    #[test]
    fn reciprocal_mosaics_violate_the_four_point_condition() {
        // Two parents and their two reciprocal mosaics: the parent pairing sums
        // to 8/4 while both crossing pairings sum to 4/4. This is the shape
        // recombination makes and a tree cannot.
        let a = aln(&[b"AAAA", b"CCCC", b"CCAA", b"AACC"]);
        let s = find_quartet_certificates(&a, 5, 1000).unwrap();
        assert_eq!(s.examined, 1);
        assert_eq!(s.violating, 1);
        let cert = &s.certificates[0];
        assert_eq!(cert.common, 4);
        assert_eq!(cert.n_ab + cert.n_ce, 8, "the dominant pairing");
        assert_eq!(cert.n_ac + cert.n_be, 4);
        assert_eq!(cert.n_ae + cert.n_bc, 4);
        assert_eq!((cert.margin_num, cert.margin_den), (4, 4));
        // The far pair must be the two parents, not a parent and a mosaic.
        let far = [cert.a.as_str(), cert.b.as_str()];
        assert!(
            far.contains(&"s0") && far.contains(&"s1"),
            "the dominant pairing should name the parents, got {far:?}"
        );
        assert!(cert.validate().is_ok());
    }

    #[test]
    fn a_tree_quartet_yields_no_quartet_certificate() {
        // Two sister pairs separated by a long internal branch: the two largest
        // pairing sums are equal, which is the four-point condition. This is the
        // control that gives the test above its content.
        let a = aln(&[b"AAAAAA", b"CAAAAA", b"AACCCC", b"ACCCCC"]);
        let q = a.quadruple_counts(0, 1, 2, 3);
        assert_eq!(q.pairing_sums(), [2, 10, 10], "two largest equal");
        let s = find_quartet_certificates(&a, 5, 1000).unwrap();
        assert_eq!(s.violating, 0);
        assert!(s.certificates.is_empty());
    }

    #[test]
    fn an_additive_population_has_zero_mean_delta_q() {
        // Four identical sequences: every pairing sum is 0, so every quadruple
        // ties and is additive. The averaged statistic must be exactly zero,
        // and no quadruple may be called violating.
        let s: Vec<u8> = "ACGTTGCAACGTTGCA".bytes().collect();
        let labels: Vec<String> = (0..4).map(|i| format!("s{i}")).collect();
        let aln = Alignment::new(labels, vec![s.clone(), s.clone(), s.clone(), s]).unwrap();
        let d = quartet_distribution(&aln, 100).expect("four sequences give one quadruple");
        assert_eq!(d.examined, 1);
        assert_eq!(d.violating, 0, "an additive quadruple is not a violation");
        assert_eq!(d.mean_delta_q, 0.0);
        assert_eq!(d.max_delta_q, 0.0);
        assert_eq!(d.proportion_violating, 0.0);
    }

    #[test]
    fn the_averaged_statistics_separate_mosaics_from_a_clonal_tree() {
        // The head-to-head the certificate's maximum is usually reported on,
        // run on the averaged statistics instead: both must rank the
        // reciprocal-mosaic positive control above the clonal null.
        use crate::simulate::{clonal_tree, with_reciprocal_recombinants};
        let clonal = clonal_tree(2000, 3, 0.004, 7);
        let (mosaic, _) = with_reciprocal_recombinants(clonal_tree(2000, 3, 0.004, 7), 0.5, 99);

        let null = Alignment::new(clonal.labels.clone(), clonal.sequences.clone()).unwrap();
        let positive = Alignment::new(mosaic.labels.clone(), mosaic.sequences.clone()).unwrap();
        let a = quartet_distribution(&null, 20_000).expect("searchable");
        let b = quartet_distribution(&positive, 20_000).expect("searchable");

        assert!(
            b.mean_delta_q > a.mean_delta_q,
            "mean delta_q must be larger with reticulation: {} against {}",
            b.mean_delta_q,
            a.mean_delta_q
        );
        assert!(
            b.max_delta_q > a.max_delta_q,
            "max delta_q must be larger with reticulation: {} against {}",
            b.max_delta_q,
            a.max_delta_q
        );
        // Sanity on the ranges, which are definitional rather than empirical.
        for d in [&a, &b] {
            assert!(
                (0.0..=1.0).contains(&d.mean_delta_q),
                "delta_q is a ratio in [0,1], got {}",
                d.mean_delta_q
            );
            assert!(d.max_delta_q >= d.mean_delta_q, "a max is at least a mean");
            assert!(d.q_residual >= 0.0, "a mean of squares is non-negative");
            assert!(d.mean_p_distance > 0.0, "these populations are not identical");
        }
    }

    #[test]
    fn quartet_search_refuses_to_sample_silently() {
        let labels: Vec<String> = (0..40).map(|i| format!("s{i}")).collect();
        let seqs: Vec<Vec<u8>> = (0..40).map(|_| b"ACGTACGT".to_vec()).collect();
        let a = Alignment::new(labels, seqs).unwrap();
        let err = find_quartet_certificates(&a, 5, 100).expect_err("must refuse");
        assert!(err.contains("91390"), "{err}");
        assert!(err.contains("will not sample"), "{err}");
        // And it proceeds when the ceiling is raised deliberately.
        assert!(find_quartet_certificates(&a, 5, 100_000).is_ok());
    }

    #[test]
    fn validation_rejects_a_forged_quartet() {
        let a = aln(&[b"AAAA", b"CCCC", b"CCAA", b"AACC"]);
        let s = find_quartet_certificates(&a, 5, 1000).unwrap();
        let mut bad = s.certificates[0].clone();
        assert!(bad.validate().is_ok());
        // Swap the dominant pairing with one of the others: the claimed
        // domination becomes false.
        std::mem::swap(&mut bad.n_ab, &mut bad.n_ac);
        std::mem::swap(&mut bad.n_ce, &mut bad.n_be);
        assert!(bad.validate().is_err());
        // Forge the margin alone, leaving the inequalities intact.
        let mut wrong_margin = s.certificates[0].clone();
        wrong_margin.margin_num += 1;
        assert!(wrong_margin.validate().is_err());
    }

    #[test]
    fn validation_rejects_a_forged_complete_triple() {
        let a = aln(&[b"AAAA", b"AACC", b"CCCC"]);
        let s = find_complete_triple_certificates(&a, 5).unwrap();
        let mut bad = s.certificates[0].clone();
        assert!(bad.validate().is_ok());
        bad.n_xz = bad.n_xy;
        assert!(bad.validate().is_err());
        let mut no_columns = s.certificates[0].clone();
        no_columns.common = 0;
        assert!(no_columns.validate().is_err());
    }

    #[test]
    fn disjoint_coverage_is_skipped_not_scored() {
        // No column has a base in all three, so the triple has no
        // complete-deletion distance and must be skipped rather than scored.
        let seqs: Vec<Vec<u8>> = vec![
            b"ACGTACGT----".to_vec(),
            b"----ACGTACGT".to_vec(),
            b"ACGT----ACGT".to_vec(),
        ];
        let a = Alignment::new(vec!["p".into(), "q".into(), "r".into()], seqs).unwrap();
        assert_eq!(a.triple_counts(0, 1, 2).common, 0);
        let s = find_complete_triple_certificates(&a, 5).unwrap();
        assert_eq!(s.skipped_no_common_columns, 1);
        assert_eq!(s.violating, 0);
    }
}

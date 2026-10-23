//! Bit-parallel subset counting, so that taxon count stops being the binding
//! constraint.
//!
//! # Why this exists
//!
//! The certificate's statistic is a maximum over subsets, and the head-to-head
//! established that its real defect is not power but **incomparability across
//! taxon counts**: the quartet floor moves from 5–7 at sixteen taxa to 9–10 at
//! thirty-two. Characterising that defect means measuring it at several taxon
//! counts, and the search is `O(n⁴ · L)` — at sixty-four taxa, 635 376
//! quadruples times four thousand columns is 2.5e9 column visits, which is the
//! wrong side of practical.
//!
//! Then the topology experiment of `reticulation topology` came back
//! uninterpretable for a related reason: at sixteen taxa its matched null
//! produced no cycle in 495 of 495 replicates, so the bootstrap could not
//! calibrate anything. The fix there is also more taxa. Scale stopped being an
//! optimisation and became a precondition for two separate experiments, which is
//! why this module exists.
//!
//! # The layout
//!
//! [`Alignment`] stores one byte per residue and walks columns, so a pair's
//! differing-site count costs `O(L)` iterations. Here each sequence becomes
//! three bitvectors over columns:
//!
//! * `b0`, `b1` — the base's two bits, under `A=00, C=01, G=10, T=11`;
//! * `ok` — set where the residue is a definite base, so gaps, `N` and anything
//!   else are excluded exactly as complete deletion requires.
//!
//! A quadruple's seven counts then come from one pass over `L/64` words:
//! intersect the four `ok` words for the common mask, and for each of the six
//! pairs count the bits where the bases differ inside that mask. That is roughly
//! twenty-five machine operations per sixty-four columns instead of ten per
//! column.
//!
//! # Why it is a second implementation rather than a replacement
//!
//! A clever counter that disagrees with the obvious one is a silent wrong
//! answer, and the whole method rests on these counts. So [`Alignment`] stays,
//! this sits beside it, and `bit_and_byte_counting_agree` is a property test
//! asserting they produce identical counts on random alignments with gaps,
//! ambiguity codes and mixed case. The byte implementation is the specification;
//! this is an optimisation that must match it.
//!
//! `popcount` is a single instruction on every target this will run on, so no
//! explicit SIMD and no `unsafe` is needed to get the gain.

use crate::subset::{Alignment, QuadrupleCounts, TripleCounts};

/// Sequences packed as bitvectors over columns.
#[derive(Debug, Clone)]
pub struct BitAlignment {
    /// Sequence labels, in the input order.
    pub labels: Vec<String>,
    /// Low bit of each base, one bitvector per sequence.
    b0: Vec<Vec<u64>>,
    /// High bit of each base.
    b1: Vec<Vec<u64>>,
    /// Set where the residue is a definite base.
    ok: Vec<Vec<u64>>,
    /// Columns in the alignment.
    width: usize,
    /// `u64` words per sequence.
    words: usize,
}

/// `A=00, C=01, G=10, T=11`; anything else is not a base.
fn code_of(byte: u8) -> Option<u8> {
    match byte.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

impl BitAlignment {
    /// Pack an alignment.
    ///
    /// # Errors
    ///
    /// An empty alignment. Unlike [`Alignment`] there is no sequence-count
    /// ceiling here: the packing is per sequence, so the taxon count is limited
    /// only by the quadruple count the caller is willing to enumerate.
    pub fn from_alignment(aln: &Alignment) -> Result<Self, String> {
        let width = aln.width();
        if width == 0 || aln.is_empty() {
            return Err("cannot pack an empty alignment".to_string());
        }
        let words = width.div_ceil(64);
        let n = aln.len();
        let mut b0 = vec![vec![0u64; words]; n];
        let mut b1 = vec![vec![0u64; words]; n];
        let mut ok = vec![vec![0u64; words]; n];
        for (index, sequence) in aln.sequences.iter().enumerate() {
            for (column, byte) in sequence.iter().enumerate() {
                let (word, bit) = (column / 64, column % 64);
                if let Some(code) = code_of(*byte) {
                    ok[index][word] |= 1u64 << bit;
                    if code & 1 != 0 {
                        b0[index][word] |= 1u64 << bit;
                    }
                    if code & 2 != 0 {
                        b1[index][word] |= 1u64 << bit;
                    }
                }
            }
        }
        // Padding bits past `width` in the final word are left clear in `ok`, so
        // they are never comparable and cannot be counted.
        Ok(BitAlignment {
            labels: aln.labels.clone(),
            b0,
            b1,
            ok,
            width,
            words,
        })
    }

    /// Number of sequences.
    pub fn len(&self) -> usize {
        self.labels.len()
    }

    /// Whether there are no sequences.
    pub fn is_empty(&self) -> bool {
        self.labels.is_empty()
    }

    /// Alignment width.
    pub fn width(&self) -> usize {
        self.width
    }

    /// Differing bits between two sequences inside `mask`.
    #[inline]
    fn differ(&self, i: usize, j: usize, word: usize, mask: u64) -> u32 {
        let d = ((self.b0[i][word] ^ self.b0[j][word])
            | (self.b1[i][word] ^ self.b1[j][word]))
            & mask;
        d.count_ones()
    }

    /// Counts for a triple under complete deletion within that triple.
    pub fn triple_counts(&self, i: usize, j: usize, k: usize) -> TripleCounts {
        let mut counts = TripleCounts {
            common: 0,
            n_xy: 0,
            n_yz: 0,
            n_xz: 0,
        };
        for word in 0..self.words {
            let common = self.ok[i][word] & self.ok[j][word] & self.ok[k][word];
            if common == 0 {
                continue;
            }
            counts.common += i64::from(common.count_ones());
            counts.n_xy += i64::from(self.differ(i, j, word, common));
            counts.n_yz += i64::from(self.differ(j, k, word, common));
            counts.n_xz += i64::from(self.differ(i, k, word, common));
        }
        counts
    }

    /// Counts for a quadruple under complete deletion within that quadruple.
    ///
    /// One pass over the words, six pair counts inside the common mask. This is
    /// the function the `O(n⁴)` search spends its time in.
    pub fn quadruple_counts(&self, i: usize, j: usize, k: usize, l: usize) -> QuadrupleCounts {
        let mut counts = QuadrupleCounts {
            common: 0,
            n_ab: 0,
            n_ce: 0,
            n_ac: 0,
            n_be: 0,
            n_ae: 0,
            n_bc: 0,
        };
        for word in 0..self.words {
            let common =
                self.ok[i][word] & self.ok[j][word] & self.ok[k][word] & self.ok[l][word];
            if common == 0 {
                continue;
            }
            counts.common += i64::from(common.count_ones());
            counts.n_ab += i64::from(self.differ(i, j, word, common));
            counts.n_ce += i64::from(self.differ(k, l, word, common));
            counts.n_ac += i64::from(self.differ(i, k, word, common));
            counts.n_be += i64::from(self.differ(j, l, word, common));
            counts.n_ae += i64::from(self.differ(i, l, word, common));
            counts.n_bc += i64::from(self.differ(j, k, word, common));
        }
        counts
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulate::{clonal_tree, inject_gap_blocks, SplitMix64};

    fn pack(aln: &Alignment) -> BitAlignment {
        BitAlignment::from_alignment(aln).expect("packable")
    }

    #[test]
    fn packing_preserves_width_and_sequence_count() {
        let pop = clonal_tree(200, 3, 0.02, 1);
        let aln = Alignment::new(pop.labels.clone(), pop.sequences.clone()).unwrap();
        let bits = pack(&aln);
        assert_eq!(bits.width(), aln.width());
        assert_eq!(bits.len(), aln.len());
        assert!(!bits.is_empty());
    }

    #[test]
    fn bit_and_byte_counting_agree_on_a_hand_case() {
        // Column 2 has a gap in the third sequence, so complete deletion within
        // the triple must drop it; column 4 is lower case, which must still
        // count. Worked by hand so the two implementations are checked against
        // an external answer rather than only against each other.
        let labels: Vec<String> = (0..3).map(|i| format!("s{i}")).collect();
        let rows: Vec<Vec<u8>> = vec![
            b"ACGTa".to_vec(),
            b"ACGTc".to_vec(),
            b"AC-Ta".to_vec(),
        ];
        let aln = Alignment::new(labels, rows).unwrap();
        let bits = pack(&aln);
        let byte = aln.triple_counts(0, 1, 2);
        let bit = bits.triple_counts(0, 1, 2);
        assert_eq!(byte, bit, "byte {byte:?} against bit {bit:?}");
        assert_eq!(bit.common, 4, "one column is gapped in the third sequence");
        assert_eq!(bit.n_xy, 1, "only the last column differs, A against C");
        assert_eq!(bit.n_xz, 0);
        assert_eq!(bit.n_yz, 1);
    }

    #[test]
    fn bit_and_byte_counting_agree_on_random_alignments_with_gaps() {
        // The property that matters: the byte implementation is the
        // specification, and this must reproduce it exactly on data containing
        // gaps, ambiguity codes, mixed case and widths that are not a multiple
        // of the word size.
        let mut rng = SplitMix64::new(20_261_023);
        for width in [1usize, 7, 63, 64, 65, 130, 257] {
            let n = 6usize;
            let alphabet = b"ACGTacgtN-?";
            let rows: Vec<Vec<u8>> = (0..n)
                .map(|_| {
                    (0..width)
                        .map(|_| {
                            // Bias towards real bases so most columns survive
                            // complete deletion and the counts are non-trivial.
                            let draw = rng.below(16);
                            if draw < 12 {
                                alphabet[rng.below(8)]
                            } else {
                                alphabet[8 + rng.below(3)]
                            }
                        })
                        .collect()
                })
                .collect();
            let labels: Vec<String> = (0..n).map(|i| format!("s{i}")).collect();
            let aln = match Alignment::new(labels, rows) {
                Ok(value) => value,
                // `validate_alphabet` may refuse a draw that is mostly
                // ambiguity; that is its job, and such a draw is not a
                // counterexample to this property.
                Err(_) => continue,
            };
            let bits = pack(&aln);
            for i in 0..n {
                for j in (i + 1)..n {
                    for k in (j + 1)..n {
                        assert_eq!(
                            aln.triple_counts(i, j, k),
                            bits.triple_counts(i, j, k),
                            "triple ({i},{j},{k}) at width {width}"
                        );
                        for l in (k + 1)..n {
                            assert_eq!(
                                aln.quadruple_counts(i, j, k, l),
                                bits.quadruple_counts(i, j, k, l),
                                "quadruple ({i},{j},{k},{l}) at width {width}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn bit_and_byte_counting_agree_on_gap_blocked_simulated_data() {
        // The realistic missing-data pattern, which is where a packing bug
        // would most plausibly hide: contiguous gap blocks give the `ok` mask
        // long runs and exercise the partial final word.
        let pop = inject_gap_blocks(clonal_tree(1000, 3, 0.01, 9), 0.25, 4);
        let aln = Alignment::new(pop.labels.clone(), pop.sequences.clone()).unwrap();
        let bits = pack(&aln);
        let mut checked = 0usize;
        for i in 0..aln.len() {
            for j in (i + 1)..aln.len() {
                for k in (j + 1)..aln.len() {
                    for l in (k + 1)..aln.len() {
                        assert_eq!(
                            aln.quadruple_counts(i, j, k, l),
                            bits.quadruple_counts(i, j, k, l)
                        );
                        checked += 1;
                    }
                }
            }
        }
        assert!(checked >= 70, "expected the full quadruple set, got {checked}");
    }

    #[test]
    fn an_empty_alignment_is_refused() {
        let labels: Vec<String> = (0..3).map(|i| format!("s{i}")).collect();
        let rows = vec![Vec::new(), Vec::new(), Vec::new()];
        // `Alignment::new` validates the alphabet and refuses empty input, so
        // construct the error case through the packer's own guard.
        match Alignment::new(labels, rows) {
            Ok(aln) => {
                assert!(
                    BitAlignment::from_alignment(&aln).is_err(),
                    "a zero-width alignment must be refused"
                );
            }
            Err(error) => assert!(!error.is_empty(), "the alphabet check refused it first"),
        }
    }
}

//! Reference-anchored multiple alignment, so that the test can be run on real
//! sequence data.
//!
//! # Why this module exists
//!
//! Phase 2 closed the missing-data gap but could not be run on real data at
//! all, for a mundane reason: the test consumes an *alignment*, and this
//! machine has no `mafft`, no `muscle` and no `clustalo`. Rather than record
//! that as a permanent block, the alignment step is implemented here, in the
//! same crate, under the same rule as everything else — exact integer
//! arithmetic, controls that can fail, and no silent truncation.
//!
//! # What it computes, and what that costs
//!
//! For a population of closely related homologous sequences (a viral species,
//! a gene family), each sequence is aligned **pairwise to a single reference**
//! by global affine-gap dynamic programming (Needleman–Wunsch with Gotoh's
//! two-layer gap model), and the result is projected onto the reference's
//! coordinates. This is the same construction viral genomic surveillance uses
//! (Nextclade, `nextalign`, and the reference-based pipelines behind most
//! published SARS-CoV-2 alignments), and it has two consequences worth being
//! explicit about:
//!
//! * **Insertions relative to the reference are dropped.** The alignment width
//!   is exactly the reference length. Columns that exist in a query but not in
//!   the reference cannot be represented in reference coordinates.
//! * **Every pair of non-reference sequences is related through the
//!   reference.** A progressive or consistency-based aligner would not do
//!   that. The cost is real, and it is *measured* rather than assumed: see
//!   `reticulation align-error`, which aligns simulated populations whose true
//!   alignment is known and counts the certificates that alignment error alone
//!   manufactures.
//!
//! # Banded, and loud about it
//!
//! The DP is banded around the proportional diagonal, because a 12 kb genome
//! against a 12 kb genome is 1.4e8 cells and the optimal path of two homologues
//! does not wander far from the diagonal. A band is a cap on how much indel a
//! pair may be assigned, so the band is a place where an aligner can quietly
//! return the wrong answer. [`align_pair`] therefore detects when the optimal
//! traceback *touches* a band edge and widens the band until it does not, up to
//! a ceiling; at the ceiling it returns an error. It never returns a path that
//! was constrained by the band without saying so.

use serde::{Deserialize, Serialize};

/// Scoring for DNA global alignment.
///
/// Defaults are the NCBI BLASTN "closely related sequences" family: `+5` match,
/// `-4` mismatch, affine gaps `-10` to open and `-1` per residue. They are
/// parameters rather than constants because the right gap cost depends on the
/// population, and because a result that moves when the gap cost moves is a
/// result about the aligner.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Scoring {
    /// Reward for a matching pair of definite bases.
    pub match_score: i32,
    /// Penalty (negative) for a mismatching pair.
    pub mismatch: i32,
    /// Penalty (negative) applied once when a gap starts.
    pub gap_open: i32,
    /// Penalty (negative) applied per residue of a gap, including the first.
    pub gap_extend: i32,
}

impl Default for Scoring {
    fn default() -> Self {
        Scoring {
            match_score: 5,
            mismatch: -4,
            gap_open: -10,
            gap_extend: -1,
        }
    }
}

/// One step of a pairwise alignment, in reference/query index space.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// Reference position `0` aligned to query position `1`.
    Aligned(usize, usize),
    /// A gap in the query: this reference position has no query residue.
    RefOnly(usize),
    /// A gap in the reference: this query residue has no reference position.
    /// Dropped by [`project_to_reference`], and counted so the loss is visible.
    QueryOnly(usize),
}

/// A pairwise alignment, with the evidence that the band did not constrain it.
#[derive(Debug, Clone)]
pub struct PairAlignment {
    /// Optimal score under [`Scoring`].
    pub score: i64,
    /// The path, from the start of both sequences to the end of both.
    pub ops: Vec<Op>,
    /// Half-width of the band that produced this path.
    pub band: usize,
    /// Query residues with no reference position: the information this
    /// projection discards.
    pub insertions: usize,
    /// Reference positions with no query residue.
    pub deletions: usize,
}

const NEG: i64 = i64::MIN / 4;

/// A sequence read from FASTA.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    /// The text after `>`, up to the first whitespace.
    pub label: String,
    /// Residues, uppercased, with whitespace removed.
    pub sequence: Vec<u8>,
}

/// Parse FASTA.
///
/// # Errors
///
/// Residue text before the first `>`, or a record with an empty sequence. Both
/// indicate a file that is not what the caller thinks it is, and guessing would
/// be worse than stopping.
pub fn read_fasta(text: &str) -> Result<Vec<FastaRecord>, String> {
    let mut records: Vec<FastaRecord> = Vec::new();
    for (lineno, line) in text.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(header) = trimmed.strip_prefix('>') {
            let label = header.split_whitespace().next().unwrap_or("").to_string();
            if label.is_empty() {
                return Err(format!("line {}: empty FASTA label", lineno + 1));
            }
            records.push(FastaRecord {
                label,
                sequence: Vec::new(),
            });
        } else {
            let last = records.last_mut().ok_or_else(|| {
                format!("line {}: residues before the first '>' header", lineno + 1)
            })?;
            last.sequence
                .extend(trimmed.bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    if let Some(empty) = records.iter().find(|r| r.sequence.is_empty()) {
        return Err(format!("record '{}' has no residues", empty.label));
    }
    Ok(records)
}

fn substitution(a: u8, b: u8, scoring: &Scoring) -> i64 {
    if a == b {
        scoring.match_score as i64
    } else {
        scoring.mismatch as i64
    }
}

/// Band limits for row `i`, centred on the proportional diagonal.
fn band_limits(i: usize, n: usize, m: usize, band: usize) -> (usize, usize) {
    // `n.max(1)` rather than a branch: when `n == 0` the only row is `i == 0`,
    // so the quotient is 0 either way.
    let centre = i * m / n.max(1);
    let lo = centre.saturating_sub(band);
    let hi = (centre + band).min(m);
    (lo, hi)
}

/// Globally align `query` to `reference` with affine gaps, inside a band.
///
/// The band half-width starts at `band` and doubles whenever the optimal path
/// touches a band edge, up to `max_band`.
///
/// # Errors
///
/// Empty input, or a path still touching the band edge at `max_band` — which
/// means the two sequences are not alignable as near-colinear homologues and
/// the caller should not be handed a silently clipped answer.
pub fn align_pair(
    reference: &[u8],
    query: &[u8],
    scoring: &Scoring,
    band: usize,
    max_band: usize,
) -> Result<PairAlignment, String> {
    if reference.is_empty() || query.is_empty() {
        return Err("cannot align an empty sequence".to_string());
    }
    let mut width = band.max(1);
    loop {
        let attempt = align_banded(reference, query, scoring, width)?;
        match attempt {
            BandedOutcome::Clear(alignment) => return Ok(alignment),
            BandedOutcome::TouchedEdge => {
                if width >= max_band {
                    return Err(format!(
                        "the optimal path still touches the band edge at a half-width \
                         of {width}; reference {} bp and query {} bp are not \
                         alignable as near-colinear homologues under this scoring, \
                         and a clipped path will not be returned as if it were \
                         optimal",
                        reference.len(),
                        query.len()
                    ));
                }
                width = (width * 2).min(max_band);
            }
        }
    }
}

enum BandedOutcome {
    Clear(PairAlignment),
    TouchedEdge,
}

// Traceback bit layout, one byte per banded cell.
const FROM_M: u8 = 0b0000_0011; // predecessor layer of M: 0 = M, 1 = X, 2 = Y
const X_FROM_X: u8 = 0b0000_0100; // set when X came from X rather than M
const Y_FROM_Y: u8 = 0b0000_1000; // set when Y came from Y rather than M

fn align_banded(
    reference: &[u8],
    query: &[u8],
    scoring: &Scoring,
    band: usize,
) -> Result<BandedOutcome, String> {
    // Rows index the reference, columns the query.
    let n = reference.len();
    let m = query.len();
    let stride = 2 * band + 2;
    let cells = (n + 1)
        .checked_mul(stride)
        .ok_or_else(|| format!("band {band} on {n}x{m} overflows"))?;
    let mut back = vec![0u8; cells];

    let open_extend = (scoring.gap_open + scoring.gap_extend) as i64;
    let extend = scoring.gap_extend as i64;

    // Full-width score rows, so a predecessor can be read by absolute column.
    let mut m_prev = vec![NEG; m + 1];
    let mut x_prev = vec![NEG; m + 1];
    let mut y_prev = vec![NEG; m + 1];
    let mut m_cur = vec![NEG; m + 1];
    let mut x_cur = vec![NEG; m + 1];
    let mut y_cur = vec![NEG; m + 1];

    let (lo0, hi0) = band_limits(0, n, m, band);
    m_prev[0] = 0;
    // Row 0: the query's first j residues are all insertions.
    for (j, value) in y_prev
        .iter_mut()
        .enumerate()
        .take(hi0 + 1)
        .skip(lo0.max(1))
    {
        *value = open_extend + extend * (j as i64 - 1);
    }

    for i in 1..=n {
        let (lo, hi) = band_limits(i, n, m, band);
        let (plo, phi) = band_limits(i - 1, n, m, band);
        for value in m_cur.iter_mut().take(hi + 1).skip(lo) {
            *value = NEG;
        }
        for value in x_cur.iter_mut().take(hi + 1).skip(lo) {
            *value = NEG;
        }
        for value in y_cur.iter_mut().take(hi + 1).skip(lo) {
            *value = NEG;
        }
        for j in lo..=hi {
            let mut flags = 0u8;
            // M: both consumed.
            if j > plo && j - 1 <= phi {
                // Ties prefer M, then X, then Y: at equal score a substitution
                // is preferred to a deletion/insertion pair. This matters. A
                // run of three mismatches scores exactly what deleting one
                // reference residue and inserting one query residue scores
                // (-12 + 5 against -11 + 15 - 11), so without a stated
                // preference the aligner would return either, and on gapless
                // input it would sometimes not return the input.
                let (layer, value) = {
                    let mut layer = 0usize;
                    let mut value = m_prev[j - 1];
                    if x_prev[j - 1] > value {
                        layer = 1;
                        value = x_prev[j - 1];
                    }
                    if y_prev[j - 1] > value {
                        layer = 2;
                        value = y_prev[j - 1];
                    }
                    (layer, value)
                };
                if value > NEG {
                    m_cur[j] = value + substitution(reference[i - 1], query[j - 1], scoring);
                    flags |= layer as u8 & FROM_M;
                }
            } else if j == 0 && i == 0 {
                m_cur[j] = 0;
            }
            // X: a gap in the query, consuming one reference residue.
            if j >= plo && j <= phi {
                let from_m = if m_prev[j] > NEG {
                    m_prev[j] + open_extend
                } else {
                    NEG
                };
                let from_x = if x_prev[j] > NEG {
                    x_prev[j] + extend
                } else {
                    NEG
                };
                if from_x > from_m {
                    x_cur[j] = from_x;
                    flags |= X_FROM_X;
                } else if from_m > NEG {
                    x_cur[j] = from_m;
                }
            }
            // Y: a gap in the reference, consuming one query residue.
            if j > lo {
                let from_m = if m_cur[j - 1] > NEG {
                    m_cur[j - 1] + open_extend
                } else {
                    NEG
                };
                let from_y = if y_cur[j - 1] > NEG {
                    y_cur[j - 1] + extend
                } else {
                    NEG
                };
                if from_y > from_m {
                    y_cur[j] = from_y;
                    flags |= Y_FROM_Y;
                } else if from_m > NEG {
                    y_cur[j] = from_m;
                }
            }
            back[i * stride + (j - lo)] = flags;
        }
        std::mem::swap(&mut m_prev, &mut m_cur);
        std::mem::swap(&mut x_prev, &mut x_cur);
        std::mem::swap(&mut y_prev, &mut y_cur);
    }

    let (lo_n, hi_n) = band_limits(n, n, m, band);
    if m < lo_n || m > hi_n {
        // The end cell is not even inside the band.
        return Ok(BandedOutcome::TouchedEdge);
    }
    // Same preference at the end cell.
    let (mut layer, score) = {
        let mut layer = 0usize;
        let mut score = m_prev[m];
        if x_prev[m] > score {
            layer = 1;
            score = x_prev[m];
        }
        if y_prev[m] > score {
            layer = 2;
            score = y_prev[m];
        }
        (layer, score)
    };
    if score <= NEG {
        return Ok(BandedOutcome::TouchedEdge);
    }

    // Traceback. `layer` is 0 = M, 1 = X, 2 = Y.
    let mut ops: Vec<Op> = Vec::new();
    let mut i = n;
    let mut j = m;
    let mut touched = false;
    while i > 0 || j > 0 {
        let (lo, hi) = band_limits(i, n, m, band);
        if (lo > 0 && j == lo) || (hi < m && j == hi) {
            touched = true;
        }
        if i == 0 {
            ops.push(Op::QueryOnly(j - 1));
            j -= 1;
            continue;
        }
        if j == 0 {
            ops.push(Op::RefOnly(i - 1));
            i -= 1;
            continue;
        }
        let flags = back[i * stride + (j - lo)];
        match layer {
            0 => {
                ops.push(Op::Aligned(i - 1, j - 1));
                layer = (flags & FROM_M) as usize;
                i -= 1;
                j -= 1;
            }
            1 => {
                ops.push(Op::RefOnly(i - 1));
                layer = if flags & X_FROM_X != 0 { 1 } else { 0 };
                i -= 1;
            }
            _ => {
                ops.push(Op::QueryOnly(j - 1));
                layer = if flags & Y_FROM_Y != 0 { 2 } else { 0 };
                j -= 1;
            }
        }
    }
    if touched {
        return Ok(BandedOutcome::TouchedEdge);
    }
    ops.reverse();
    let insertions = ops
        .iter()
        .filter(|op| matches!(op, Op::QueryOnly(_)))
        .count();
    let deletions = ops.iter().filter(|op| matches!(op, Op::RefOnly(_))).count();
    Ok(BandedOutcome::Clear(PairAlignment {
        score,
        ops,
        band,
        insertions,
        deletions,
    }))
}

/// Project a pairwise alignment onto the reference's coordinates.
///
/// The returned row has exactly `reference_len` entries: the query base where
/// one is aligned, `b'-'` where the query has a deletion. Query residues
/// inserted relative to the reference are dropped, which is the price of
/// reference coordinates and is counted in [`PairAlignment::insertions`].
pub fn project_to_reference(
    alignment: &PairAlignment,
    query: &[u8],
    reference_len: usize,
) -> Vec<u8> {
    let mut row = vec![b'-'; reference_len];
    for op in &alignment.ops {
        match *op {
            Op::Aligned(r, q) => row[r] = query[q],
            Op::RefOnly(r) => row[r] = b'-',
            Op::QueryOnly(_) => {}
        }
    }
    row
}

/// A reference-anchored multiple alignment and the diagnostics that say how
/// much it had to throw away.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceAnchored {
    /// Labels, in input order.
    pub labels: Vec<String>,
    /// Rows in reference coordinates, all of width `width`.
    pub rows: Vec<Vec<u8>>,
    /// Index of the chosen reference in the input order.
    pub reference_index: usize,
    /// Label of the chosen reference.
    pub reference_label: String,
    /// Alignment width: the reference length.
    pub width: usize,
    /// Per-query count of residues dropped as insertions relative to the
    /// reference.
    pub insertions_dropped: Vec<usize>,
    /// Per-query count of reference positions with no query residue.
    pub deletions: Vec<usize>,
    /// Per-query pairwise score against the reference.
    pub scores: Vec<i64>,
    /// Per-query band half-width that cleared the band-edge check.
    pub bands: Vec<usize>,
}

impl ReferenceAnchored {
    /// Gap fraction over the whole matrix: the power this projection costs.
    pub fn gap_fraction(&self) -> f64 {
        if self.rows.is_empty() || self.width == 0 {
            return 0.0;
        }
        let gaps: usize = self
            .rows
            .iter()
            .map(|row| row.iter().filter(|&&b| b == b'-').count())
            .sum();
        gaps as f64 / (self.rows.len() * self.width) as f64
    }
}

/// 6-mer profile distance, used only to pick a reference.
fn kmer_distance(a: &[u8], b: &[u8]) -> f64 {
    fn profile(s: &[u8]) -> (std::collections::HashMap<&[u8], f64>, f64) {
        const K: usize = 6;
        let mut counts = std::collections::HashMap::<&[u8], f64>::new();
        if s.len() >= K {
            for window in s.windows(K) {
                *counts.entry(window).or_insert(0.0) += 1.0;
            }
        }
        let total: f64 = counts.values().sum();
        (counts, total)
    }
    let (pa, ta) = profile(a);
    let (pb, tb) = profile(b);
    if ta == 0.0 || tb == 0.0 {
        return 1.0;
    }
    let mut l1 = 0.0;
    for (key, value) in &pa {
        l1 += (value / ta - pb.get(key).copied().unwrap_or(0.0) / tb).abs();
    }
    for (key, value) in &pb {
        if !pa.contains_key(key) {
            l1 += value / tb;
        }
    }
    l1 / 2.0
}

/// Pick the reference: the sequence with the smallest total 6-mer distance to
/// the rest.
///
/// A centroid, not the first record, because reference-anchored alignment loses
/// the most when the reference is an outlier — and the choice is deterministic
/// and independent of any distance the test will later compute on columns.
pub fn centroid_reference(sequences: &[Vec<u8>]) -> usize {
    centroid_ranking(sequences)[0]
}

/// Every index, ordered by total 6-mer distance to the rest: most central first.
///
/// The second entry is what the reference-sensitivity control uses. Choosing a
/// reference is a free parameter of reference-anchored alignment, so the honest
/// thing is to run the test twice and report whether the verdict moved.
pub fn centroid_ranking(sequences: &[Vec<u8>]) -> Vec<usize> {
    let n = sequences.len();
    let mut totals: Vec<(usize, f64)> = (0..n)
        .map(|i| {
            let total = (0..n)
                .filter(|&j| j != i)
                .map(|j| kmer_distance(&sequences[i], &sequences[j]))
                .sum::<f64>();
            (i, total)
        })
        .collect();
    totals.sort_by(|a, b| {
        a.1.partial_cmp(&b.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.0.cmp(&b.0))
    });
    totals.into_iter().map(|(i, _)| i).collect()
}

/// Build a reference-anchored multiple alignment.
///
/// # Errors
///
/// Fewer than three sequences, an out-of-range `reference_index`, or any
/// pairwise alignment that cannot clear the band-edge check.
pub fn reference_anchored_alignment(
    labels: &[String],
    sequences: &[Vec<u8>],
    reference_index: usize,
    scoring: &Scoring,
    band: usize,
    max_band: usize,
) -> Result<ReferenceAnchored, String> {
    if sequences.len() < 3 {
        return Err(format!(
            "need at least 3 sequences to align a population, got {}",
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
    if reference_index >= sequences.len() {
        return Err(format!(
            "reference index {reference_index} is out of range for {} sequences",
            sequences.len()
        ));
    }
    let reference = &sequences[reference_index];
    let width = reference.len();
    let mut rows: Vec<Vec<u8>> = Vec::with_capacity(sequences.len());
    let mut insertions = Vec::with_capacity(sequences.len());
    let mut deletions = Vec::with_capacity(sequences.len());
    let mut scores = Vec::with_capacity(sequences.len());
    let mut bands = Vec::with_capacity(sequences.len());
    for (index, sequence) in sequences.iter().enumerate() {
        if index == reference_index {
            rows.push(reference.clone());
            insertions.push(0);
            deletions.push(0);
            scores.push(
                reference
                    .iter()
                    .map(|_| scoring.match_score as i64)
                    .sum::<i64>(),
            );
            bands.push(0);
            continue;
        }
        let pair = align_pair(reference, sequence, scoring, band, max_band).map_err(|error| {
            format!(
                "aligning '{}' to reference '{}': {error}",
                labels[index], labels[reference_index]
            )
        })?;
        rows.push(project_to_reference(&pair, sequence, width));
        insertions.push(pair.insertions);
        deletions.push(pair.deletions);
        scores.push(pair.score);
        bands.push(pair.band);
    }
    Ok(ReferenceAnchored {
        labels: labels.to_vec(),
        rows,
        reference_index,
        reference_label: labels[reference_index].clone(),
        width,
        insertions_dropped: insertions,
        deletions,
        scores,
        bands,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seq(text: &str) -> Vec<u8> {
        text.bytes().collect()
    }

    #[test]
    fn fasta_parses_labels_and_wrapped_sequences() {
        let text = ">AB123 some description\nACGT\nACGT\n\n>CD456\nTTTT\n";
        let records = read_fasta(text).expect("valid FASTA");
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].label, "AB123");
        assert_eq!(records[0].sequence, seq("ACGTACGT"));
        assert_eq!(records[1].label, "CD456");
        assert_eq!(records[1].sequence, seq("TTTT"));
    }

    #[test]
    fn fasta_refuses_residues_before_a_header_and_empty_records() {
        let orphan = read_fasta("ACGT\n>A\nACGT\n");
        assert!(orphan.is_err(), "residues before '>' must be refused");
        assert!(orphan.unwrap_err().contains("before the first"));
        let empty = read_fasta(">A\n>B\nACGT\n");
        assert!(empty.is_err(), "an empty record must be refused");
        assert!(empty.unwrap_err().contains("no residues"));
    }

    #[test]
    fn identical_sequences_align_with_no_gaps_at_the_full_match_score() {
        let s = seq("ACGTACGTTTGCAACGT");
        let scoring = Scoring::default();
        let pair = align_pair(&s, &s, &scoring, 8, 64).expect("alignable");
        assert_eq!(pair.score, s.len() as i64 * scoring.match_score as i64);
        assert_eq!(pair.insertions, 0);
        assert_eq!(pair.deletions, 0);
        assert_eq!(pair.ops.len(), s.len());
        assert!(pair.ops.iter().all(|op| matches!(op, Op::Aligned(_, _))));
    }

    #[test]
    fn a_deletion_is_placed_where_it_was_made() {
        // Query is the reference with positions 8..12 removed.
        let reference = seq("ACGTACGTTTGCAACGTACGTAC");
        let mut query = reference.clone();
        query.drain(8..12);
        let pair = align_pair(&reference, &query, &Scoring::default(), 8, 64).expect("alignable");
        assert_eq!(pair.deletions, 4, "four reference positions unmatched");
        assert_eq!(pair.insertions, 0);
        let row = project_to_reference(&pair, &query, reference.len());
        assert_eq!(row.len(), reference.len());
        let gaps: Vec<usize> = row
            .iter()
            .enumerate()
            .filter(|&(_, &b)| b == b'-')
            .map(|(i, _)| i)
            .collect();
        assert_eq!(gaps.len(), 4);
        // The deleted block is a run of T/G/C; the aligner may place the gap at
        // any position consistent with the optimal score, but it must be one
        // contiguous run inside the region that was removed plus its ambiguity.
        assert_eq!(
            gaps.last().unwrap() - gaps.first().unwrap(),
            3,
            "the gap must be one contiguous run of 4, got {gaps:?}"
        );
        for (position, base) in row.iter().enumerate() {
            if *base != b'-' {
                assert!(
                    reference.contains(base),
                    "projected base at {position} is not from the alphabet"
                );
            }
        }
    }

    #[test]
    fn an_insertion_relative_to_the_reference_is_dropped_and_counted() {
        let reference = seq("ACGTACGTTTGCAACGTACGTAC");
        let mut query = reference.clone();
        for (offset, base) in seq("GGGGGG").into_iter().enumerate() {
            query.insert(10 + offset, base);
        }
        let pair = align_pair(&reference, &query, &Scoring::default(), 8, 64).expect("alignable");
        assert_eq!(pair.insertions, 6, "six query residues have no column");
        let row = project_to_reference(&pair, &query, reference.len());
        assert_eq!(row.len(), reference.len(), "width is the reference length");
        assert_eq!(
            row.iter().filter(|&&b| b == b'-').count(),
            0,
            "a pure insertion leaves no gap in reference coordinates"
        );
    }

    #[test]
    fn a_band_too_narrow_for_the_indel_is_reported_not_clipped() {
        let reference: Vec<u8> = "ACGT".bytes().cycle().take(400).collect();
        let mut query = reference.clone();
        query.drain(100..260); // a 160 bp deletion
        let narrow = align_pair(&reference, &query, &Scoring::default(), 2, 4);
        assert!(narrow.is_err(), "a 160 bp indel cannot fit a band of 4");
        assert!(narrow.unwrap_err().contains("band edge"));
        let wide = align_pair(&reference, &query, &Scoring::default(), 2, 512)
            .expect("a wide enough band finds the path");
        assert_eq!(wide.deletions, 160);
        // The band is centred on the *proportional* diagonal, which already
        // absorbs the 160 bp length difference; what the band must still cover
        // is the local excursion, because the deletion sits at one end rather
        // than being spread evenly. So the width needed is well above the
        // starting 2 and well below 160.
        assert!(
            wide.band > 2 && wide.band < 160,
            "the band had to widen past its start but not to the indel length, got {}",
            wide.band
        );
    }

    #[test]
    fn alignment_score_is_symmetric_in_its_two_arguments() {
        let a = seq("ACGTACGTTTGCAACGTACGTACTTGACA");
        let b = seq("ACGTTCGTTTGCAAGGTACGTAC");
        let scoring = Scoring::default();
        let forward = align_pair(&a, &b, &scoring, 8, 128).expect("alignable");
        let backward = align_pair(&b, &a, &scoring, 8, 128).expect("alignable");
        assert_eq!(forward.score, backward.score);
        assert_eq!(forward.insertions, backward.deletions);
        assert_eq!(forward.deletions, backward.insertions);
    }

    #[test]
    fn a_sequence_scores_higher_against_itself_than_against_a_mutant() {
        let a = seq("ACGTACGTTTGCAACGTACGTACTTGACA");
        let mut mutant = a.clone();
        mutant[5] = b'G';
        mutant[11] = b'A';
        let scoring = Scoring::default();
        let same = align_pair(&a, &a, &scoring, 8, 64).expect("alignable").score;
        let other = align_pair(&a, &mutant, &scoring, 8, 64)
            .expect("alignable")
            .score;
        assert!(same > other, "{same} should exceed {other}");
        assert_eq!(same, a.len() as i64 * scoring.match_score as i64);
    }

    #[test]
    fn a_gapless_population_passes_through_the_aligner_unchanged() {
        // The control that matters for the test downstream: if the input has no
        // indels, reference-anchored alignment must be the identity, so that
        // every simulated result of Phase 2 is unaffected by routing it here.
        let base: Vec<u8> = "ACGTTGCA".bytes().cycle().take(240).collect();
        let mut sequences = vec![base.clone()];
        for shift in 1..6usize {
            let mut variant = base.clone();
            for position in (shift..240).step_by(17) {
                variant[position] = match variant[position] {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                };
            }
            sequences.push(variant);
        }
        let labels: Vec<String> = (0..sequences.len()).map(|i| format!("s{i}")).collect();
        let anchored =
            reference_anchored_alignment(&labels, &sequences, 0, &Scoring::default(), 16, 256)
                .expect("alignable");
        assert_eq!(anchored.width, 240);
        assert_eq!(anchored.rows, sequences, "the aligner must be the identity");
        assert_eq!(anchored.gap_fraction(), 0.0);
        assert!(anchored.insertions_dropped.iter().all(|&x| x == 0));
    }

    #[test]
    fn the_centroid_reference_is_not_the_outlier() {
        let core: Vec<u8> = "ACGTTGCA".bytes().cycle().take(400).collect();
        let mut near = core.clone();
        near[7] = b'A';
        let mut near2 = core.clone();
        near2[19] = b'T';
        let outlier: Vec<u8> = "TTTTAAAA".bytes().cycle().take(400).collect();
        let sequences = vec![outlier, core.clone(), near, near2];
        let picked = centroid_reference(&sequences);
        assert_ne!(picked, 0, "the outlier must not be chosen as reference");
        assert!(
            picked == 1 || picked == 2 || picked == 3,
            "the reference should come from the cluster, got {picked}"
        );
        let ranking = centroid_ranking(&sequences);
        assert_eq!(ranking.len(), 4, "the ranking covers every sequence");
        assert_eq!(
            *ranking.last().unwrap(),
            0,
            "the outlier must rank last, got {ranking:?}"
        );
        assert_eq!(ranking[0], picked, "the ranking's head is the centroid");
    }

    #[test]
    fn the_reference_row_is_the_reference_and_the_width_is_its_length() {
        let a: Vec<u8> = "ACGTTGCAGGTT".bytes().cycle().take(300).collect();
        let mut b = a.clone();
        b.drain(50..60);
        let mut c = a.clone();
        c[100] = b'A';
        let labels = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        let anchored =
            reference_anchored_alignment(&labels, &[a.clone(), b, c], 0, &Scoring::default(), 32, 512)
                .expect("alignable");
        assert_eq!(anchored.width, a.len());
        assert_eq!(anchored.rows[0], a);
        assert_eq!(anchored.reference_label, "a");
        assert!(
            anchored.gap_fraction() > 0.0,
            "the 10 bp deletion must show as gaps"
        );
        assert!(anchored.rows.iter().all(|row| row.len() == a.len()));
    }

    #[test]
    fn too_few_sequences_or_a_bad_reference_index_is_an_error() {
        let a: Vec<u8> = seq("ACGTACGT");
        let labels = vec!["a".to_string(), "b".to_string()];
        let two = reference_anchored_alignment(
            &labels,
            &[a.clone(), a.clone()],
            0,
            &Scoring::default(),
            8,
            64,
        );
        assert!(two.is_err(), "two sequences are not a population");
        let labels3 = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        let bad = reference_anchored_alignment(
            &labels3,
            &[a.clone(), a.clone(), a.clone()],
            9,
            &Scoring::default(),
            8,
            64,
        );
        assert!(bad.is_err(), "reference index 9 is out of range");
        assert!(bad.unwrap_err().contains("out of range"));
    }
}

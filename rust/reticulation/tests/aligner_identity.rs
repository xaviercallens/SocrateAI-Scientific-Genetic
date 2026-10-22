//! The aligner's identity control, as an integration test.
//!
//! Prediction A1 of `reticulation align-error` is that reference-anchored
//! alignment is the identity on gapless input, so that no Phase 2 result
//! changes by routing it through the aligner. That claim is cheap to state and
//! was false on the first implementation: with a `-10 / -1` affine gap, a run
//! of three mismatches scores exactly what deleting one reference residue and
//! inserting one query residue scores, and the aligner returned the tied
//! alternative at 1 of 5 seeds. The fix was a stated tie preference (M over X
//! over Y: a substitution over a gap pair), and this test is the regression.

use reticulation::align::{
    align_pair, centroid_reference, reference_anchored_alignment, Op, Scoring,
};
use reticulation::simulate::{clonal_tree_relaxed_with_deletions, DeletionModel};

#[test]
fn the_aligner_is_the_identity_on_gapless_simulated_populations() {
    let scoring = Scoring::default();
    for seed in 1..=5u64 {
        let pop = clonal_tree_relaxed_with_deletions(
            2000,
            4,
            0.004,
            1.0,
            &DeletionModel {
                rate: 0.0,
                mean_length: 30,
                max_gap_fraction: 0.2,
            },
            seed,
        );
        let reference = centroid_reference(&pop.raw);
        let anchored =
            reference_anchored_alignment(&pop.labels, &pop.raw, reference, &scoring, 64, 2048, true)
                .expect("gapless equal-length input is alignable");
        assert_eq!(
            anchored.rows, pop.truth.sequences,
            "seed {seed}: the aligner must return gapless input unchanged"
        );
        assert_eq!(
            anchored.insertions_dropped.iter().sum::<usize>(),
            0,
            "seed {seed}: nothing may be dropped as an insertion"
        );
        assert_eq!(anchored.gap_fraction(), 0.0, "seed {seed}: no gaps");
    }
}

#[test]
fn a_three_mismatch_cluster_is_not_turned_into_a_tied_gap_pair() {
    // The exact configuration that broke A1. Reference and query are the same
    // length and differ over three consecutive sites, so the diagonal is one of
    // two optimal alignments and the aligner must return it.
    let reference: Vec<u8> = b"CCCTTGGGGAACGCTAGTCATAGA".to_vec();
    let query: Vec<u8> = b"CCCTTGGGGAACCTGAGTCATAGA".to_vec();
    let scoring = Scoring::default();
    let pair = align_pair(&reference, &query, &scoring, 16, 256).expect("alignable");

    let mismatches = reference
        .iter()
        .zip(&query)
        .filter(|(left, right)| left != right)
        .count();
    assert_eq!(mismatches, 3, "the fixture must carry three mismatches");
    let diagonal = (reference.len() - mismatches) as i64 * scoring.match_score as i64
        + mismatches as i64 * scoring.mismatch as i64;
    assert_eq!(
        pair.score, diagonal,
        "the diagonal must be optimal, so the optimum must equal its score"
    );
    assert_eq!(pair.insertions, 0, "no insertion may be introduced");
    assert_eq!(pair.deletions, 0, "no deletion may be introduced");
    assert!(
        pair.ops.iter().all(|op| matches!(op, Op::Aligned(_, _))),
        "every step must be an aligned pair"
    );
}

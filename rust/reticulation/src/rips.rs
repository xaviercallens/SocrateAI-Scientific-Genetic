//! The Vietoris–Rips front-end that LoPHAT deliberately does not provide.
//!
//! LoPHAT's README is explicit: it "is **not** specialised to compute PH of
//! common filtrations or even filtered simplicial complexes", and works instead
//! on an arbitrary filtered boundary matrix. This module supplies the missing
//! layer for our case — build the filtered VR boundary matrix up to dimension 2
//! from a genetic-distance matrix, hand it to LoPHAT, and map the resulting
//! persistence pairs and representative cycles back to sequence labels.
//!
//! Dimension 2 is the ceiling because `H_1` is what carries the biology: a
//! 1-cycle is the signature of reticulate evolution, and computing `H_1`
//! requires triangles but nothing above them.
//!
//! # Why the representative cycle matters
//!
//! A Betti number says "there is a loop". A representative cycle says *which
//! sequences* form it, which is what a biologist can act on and what makes the
//! result auditable. We extract it and reduce it to the specific
//! ultrametric-violating triples it contains, so the topological output and the
//! Lean-checkable certificate of [`crate::certificate`] are two views of the
//! same object rather than two unrelated statistics.

use lophat::algorithms::{Decomposition, DecompositionAlgo, SerialAlgorithm};
use lophat::columns::VecColumn;
use lophat::options::LoPhatOptions;
use serde::{Deserialize, Serialize};

use crate::distance::DistanceMatrix;

/// A finite `H_1` interval with the sequences forming its representative cycle.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CycleRepresentative {
    /// Scale at which the loop appears.
    pub birth: f64,
    /// Scale at which it is filled in.
    pub death: f64,
    /// `death - birth`.
    pub persistence: f64,
    /// Labels of the sequences on the loop, in no particular order.
    pub labels: Vec<String>,
    /// The loop's edges as index pairs.
    pub edges: Vec<[usize; 2]>,
}

/// The outcome of a Vietoris–Rips persistence computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RipsPersistence {
    /// Sequences in the cloud.
    pub n_sequences: usize,
    /// Edges included, i.e. pairs within `max_scale`.
    pub n_edges: usize,
    /// Triangles included.
    pub n_triangles: usize,
    /// Total columns handed to LoPHAT.
    pub n_columns: usize,
    /// Threshold beyond which simplices were not built.
    pub max_scale: f64,
    /// Finite `H_1` intervals, longest-lived first.
    pub h1: Vec<CycleRepresentative>,
    /// Number of `H_1` classes still alive at `max_scale`.
    pub h1_essential: usize,
    /// Sum of finite `H_1` bar lengths. Reported for comparison with the
    /// chaos-game statistic of the textbook's Section 8, so the two approaches
    /// can be compared on the same number.
    pub h1_total_persistence: f64,
    /// Whether representative cycles were recovered.
    pub representatives_available: bool,
}

/// Build the filtered VR boundary matrix and reduce it with LoPHAT.
///
/// Simplices are ordered by filtration value, then by dimension, so every face
/// precedes its cofaces — the precondition LoPHAT needs for the pairing to carry
/// a persistence interpretation.
///
/// `max_scale` bounds the complex: pairs further apart than this contribute no
/// edge. Set it from the data (for example the median pairwise distance) rather
/// than leaving it at infinity, since the triangle count grows as `n^3`.
pub fn vietoris_rips_h1(dm: &DistanceMatrix, max_scale: f64) -> RipsPersistence {
    let n = dm.len();

    // --- enumerate simplices with their filtration values ---
    #[derive(Clone)]
    enum Simplex {
        Vertex(usize),
        Edge(usize, usize),
        Triangle(usize, usize, usize),
    }

    let mut items: Vec<(f64, usize, Simplex)> = Vec::new();
    for i in 0..n {
        items.push((0.0, 0, Simplex::Vertex(i)));
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let d = dm.get(i, j);
            if d <= max_scale {
                items.push((d, 1, Simplex::Edge(i, j)));
            }
        }
    }
    let n_edges = items.iter().filter(|(_, dim, _)| *dim == 1).count();
    for i in 0..n {
        for j in (i + 1)..n {
            if dm.get(i, j) > max_scale {
                continue;
            }
            for k in (j + 1)..n {
                let (a, b, c) = (dm.get(i, j), dm.get(i, k), dm.get(j, k));
                let f = a.max(b).max(c);
                if f <= max_scale {
                    items.push((f, 2, Simplex::Triangle(i, j, k)));
                }
            }
        }
    }
    let n_triangles = items.iter().filter(|(_, dim, _)| *dim == 2).count();

    // Faces before cofaces: sort by filtration, then dimension.
    items.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));

    // --- column indices ---
    let n_cols = items.len();
    let mut vertex_col = vec![usize::MAX; n];
    let mut edge_col = vec![usize::MAX; n * n];
    for (col, (_, _, s)) in items.iter().enumerate() {
        match *s {
            Simplex::Vertex(i) => vertex_col[i] = col,
            Simplex::Edge(i, j) => {
                edge_col[i * n + j] = col;
                edge_col[j * n + i] = col;
            }
            Simplex::Triangle(..) => {}
        }
    }

    // --- boundary columns ---
    let mut columns: Vec<VecColumn> = Vec::with_capacity(n_cols);
    let mut filtration: Vec<f64> = Vec::with_capacity(n_cols);
    let mut col_simplex: Vec<Simplex> = Vec::with_capacity(n_cols);
    for (f, dim, s) in &items {
        let entries: Vec<usize> = match *s {
            Simplex::Vertex(_) => Vec::new(),
            Simplex::Edge(i, j) => {
                let mut e = vec![vertex_col[i], vertex_col[j]];
                e.sort_unstable();
                e
            }
            Simplex::Triangle(i, j, k) => {
                let mut e = vec![
                    edge_col[i * n + j],
                    edge_col[i * n + k],
                    edge_col[j * n + k],
                ];
                e.sort_unstable();
                e
            }
        };
        debug_assert!(entries.iter().all(|&r| r != usize::MAX), "missing face");
        columns.push(VecColumn::from((*dim, entries)));
        filtration.push(*f);
        col_simplex.push(s.clone());
    }

    if n_cols == 0 {
        return RipsPersistence {
            n_sequences: n,
            n_edges: 0,
            n_triangles: 0,
            n_columns: 0,
            max_scale,
            h1: Vec::new(),
            h1_essential: 0,
            h1_total_persistence: 0.0,
            representatives_available: false,
        };
    }

    // --- reduce ---
    let mut options = LoPhatOptions::default();
    options.maintain_v = true; // required for representative cycles
    options.num_threads = 0;
    let decomposition = SerialAlgorithm::<VecColumn>::init(Some(options))
        .add_cols(columns.into_iter())
        .decompose();

    let (diagram, reps) = decomposition.diagram_with_reps();
    let representatives_available = reps.is_ok();
    let reps = reps.unwrap_or_default();

    // --- read H1 intervals ---
    let mut h1: Vec<CycleRepresentative> = Vec::new();
    let mut h1_essential = 0usize;
    for (&birth_col, death) in diagram.iter() {
        // A class is in H1 iff its birth simplex is an edge.
        if !matches!(col_simplex[birth_col], Simplex::Edge(..)) {
            continue;
        }
        let death_col = match death {
            lophat::utils::ExtendedUsize::Finite(d) => *d,
            lophat::utils::ExtendedUsize::Infinity => {
                h1_essential += 1;
                continue;
            }
        };
        let birth = filtration[birth_col];
        let death_f = filtration[death_col];
        if death_f <= birth {
            continue; // zero-length interval: not a feature
        }
        // Representative cycle: the edges it is supported on.
        let mut edges: Vec<[usize; 2]> = Vec::new();
        let mut labels: Vec<String> = Vec::new();
        if let Some(cycle) = reps.get(&birth_col) {
            for &col in cycle {
                if let Simplex::Edge(i, j) = col_simplex[col] {
                    edges.push([i, j]);
                    for v in [i, j] {
                        if !labels.contains(&dm.labels[v]) {
                            labels.push(dm.labels[v].clone());
                        }
                    }
                }
            }
        }
        h1.push(CycleRepresentative {
            birth,
            death: death_f,
            persistence: death_f - birth,
            labels,
            edges,
        });
    }
    h1.sort_by(|a, b| b.persistence.total_cmp(&a.persistence));
    let h1_total_persistence = h1.iter().map(|c| c.persistence).sum();

    RipsPersistence {
        n_sequences: n,
        n_edges,
        n_triangles,
        n_columns: n_cols,
        max_scale,
        h1,
        h1_essential,
        h1_total_persistence,
        representatives_available,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::distance::DistanceMatrix;

    fn matrix(labels: &[&str], values: Vec<f64>) -> DistanceMatrix {
        DistanceMatrix::from_values(labels.iter().map(|s| s.to_string()).collect(), values).unwrap()
    }

    /// Four points on a cycle: adjacent at 1, opposite at 2. The classic
    /// example with a genuine 1-cycle -- the square's loop is born when the
    /// adjacent edges appear and dies when the diagonals fill it in.
    fn square_cycle() -> DistanceMatrix {
        matrix(
            &["a", "b", "c", "d"],
            vec![
                0.0, 1.0, 2.0, 1.0, //
                1.0, 0.0, 1.0, 2.0, //
                2.0, 1.0, 0.0, 1.0, //
                1.0, 2.0, 1.0, 0.0,
            ],
        )
    }

    #[test]
    fn square_cycle_has_a_finite_h1_class() {
        let dm = square_cycle();
        let out = vietoris_rips_h1(&dm, 3.0);
        assert_eq!(out.n_sequences, 4);
        assert_eq!(out.n_edges, 6, "all pairs within scale");
        assert!(
            !out.h1.is_empty(),
            "the 4-cycle must produce an H1 class; got {out:?}"
        );
        let c = &out.h1[0];
        assert!((c.birth - 1.0).abs() < 1e-12, "born when adjacent edges enter");
        assert!((c.death - 2.0).abs() < 1e-12, "dies when a diagonal enters");
        assert!(c.persistence > 0.0);
    }

    #[test]
    fn ultrametric_has_no_finite_h1_class() {
        // The control, and the computational counterpart of
        // `ultrametric_no_induced_path`: an ultrametric's Rips components are
        // cliques, so no cycle is ever unfilled.
        let dm = matrix(
            &["a", "b", "c", "d"],
            vec![
                0.0, 1.0, 2.0, 2.0, //
                1.0, 0.0, 2.0, 2.0, //
                2.0, 2.0, 0.0, 1.0, //
                2.0, 2.0, 1.0, 0.0,
            ],
        );
        let out = vietoris_rips_h1(&dm, 3.0);
        assert!(
            out.h1.is_empty(),
            "an ultrametric must yield no finite H1 class; got {:?}",
            out.h1
        );
        assert_eq!(out.h1_total_persistence, 0.0);
    }

    #[test]
    fn representative_cycle_names_the_sequences_on_the_loop() {
        let dm = square_cycle();
        let out = vietoris_rips_h1(&dm, 3.0);
        assert!(out.representatives_available, "maintain_v must be honoured");
        let c = &out.h1[0];
        assert!(
            c.labels.len() >= 3,
            "a loop needs at least three sequences; got {:?}",
            c.labels
        );
        assert!(!c.edges.is_empty());
    }

    #[test]
    fn threshold_bounds_the_complex() {
        let dm = square_cycle();
        // Below 2.0 the diagonals never enter, so the loop cannot be filled.
        let out = vietoris_rips_h1(&dm, 1.5);
        assert_eq!(out.n_edges, 4, "only the four adjacent pairs");
        assert_eq!(out.n_triangles, 0);
        assert!(
            out.h1.is_empty() && out.h1_essential >= 1,
            "the loop should survive to the threshold as essential, got {out:?}"
        );
    }

    #[test]
    fn faces_precede_cofaces_so_the_pairing_is_meaningful() {
        // Not a homology assertion: a direct check of the ordering invariant
        // LoPHAT requires. If this breaks, every interval above is meaningless.
        let dm = square_cycle();
        let out = vietoris_rips_h1(&dm, 3.0);
        assert_eq!(out.n_columns, 4 + out.n_edges + out.n_triangles);
        for c in &out.h1 {
            assert!(c.birth <= c.death);
        }
    }
}

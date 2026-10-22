//! Proof-carrying detection of reticulate evolution.
//!
//! # Why this crate exists
//!
//! Section 8 of the accompanying textbook embedded a *single* DNA sequence by the
//! chaos game and measured the topology of the resulting point cloud. That was a
//! mistake of a specific kind: the chaos game contracts any sequence onto a
//! self-similar attractor, so the holes being measured were the attractor's and
//! not the sequence's. With no mathematical null available, the test had to be
//! calibrated against shuffles — which bought a p-value floor, a frame-sampling
//! artefact, and a statistic with no proved meaning. The honest result was
//! negative.
//!
//! This crate implements the alternative. Following Chan, Carlsson and Rabadan
//! (PNAS 2013), topology is computed on a *population* of sequences under a
//! genetic-distance metric, where a 1-cycle has an interpretation: it indicates
//! non-tree-like, reticulate evolution — recombination, reassortment, horizontal
//! transfer.
//!
//! # The part that is new here
//!
//! The null is a **theorem**, not an ensemble. An ultrametric — a strict
//! molecular clock on a clonal genealogy — admits no unfilled cycle at any
//! scale, because `d(·,·) <= eps` is then transitive and every Vietoris-Rips
//! component is a clique. So an unfilled triangle is not "surprising relative to
//! shuffles"; it is a **certificate** that the data is not ultrametric.
//!
//! Each certificate is three labels and three distances. The inference rule from
//! certificate to conclusion is proved in
//! `formal/ANSE/ReticulationCertificate.lean`
//! (`certificate_refutes_ultrametric`), its non-vacuity is witnessed there, and
//! the control — that no certificate exists under the null — is
//! `ultrametric_admits_no_certificate`.
//!
//! # What the kernel does and does not vouch for
//!
//! An adversarial audit of the first version of this crate falsified the claim
//! we had written here, which was that "the scientific conclusion does not
//! depend on the correctness of this crate". It does not hold in that form, and
//! the corrected statement is worth being precise about.
//!
//! The integer-only certificates of [`certificate`] and [`exact`] are checked
//! *given the counts*: the kernel proves an arithmetic fact about six integers.
//! Nothing in such a certificate ties those integers to an alignment, so a
//! mislabelled triple or a transposed pair in the search would produce a file
//! that still typechecks. The counting was trusted, and we had said otherwise.
//!
//! [`emit_lean::render_with_sequences`] closes that gap. It embeds the three
//! sequences and the kernel computes the differing and comparable sites itself
//! by `rfl`, then instantiates
//! `ANSE.GeneticaReticulata.sequences_refute_from_counts` against those computed
//! counts. Measured at 4000 sites: about 50 s, depending on no axioms beyond the
//! whitelist.
//!
//! So the accurate claim is:
//!
//! * **Excluded by the kernel:** a false inequality, a wrong count, a
//!   transposed pair, a miscomputed margin — anything downstream of the
//!   sequences, when the sequence-level form is used.
//! * **Still trusted:** that the embedded lists are the dataset's sequences.
//!   That is checkable by diffing them against the alignment, and is far
//!   narrower than trusting a search, but it is not nothing.
//! * **Not a kernel matter at all:** whether a violation of ultrametricity
//!   means recombination. See the limitation below.
//!
//! # Which statistic to report, and why it is not the obvious one
//!
//! [`exact`] implements the ultrametric (clock-level) test; [`subset`]
//! implements the four-point (tree-level) test. `reticulation experiment`
//! compares them on seven simulated populations at five seeds each, and the
//! result is one-sided enough to change the recommendation.
//!
//! On a **clonal tree with a relaxed clock and no reticulation at all** —
//! divergence-matched to the strict-clock control, so only the clock differs —
//! the ultrametric margin reaches `106/4000`, larger than every
//! single-recombinant margin in the positive control, and exceeds twice the
//! clonal floor at 3 of 5 seeds. The four-point margin on the same data does not
//! move: 0 of 5 seeds exceed twice the floor, and the worst-seed ratio is 0.14.
//! Meanwhile the four-point test separates its own positive control by 16.6× at
//! its weakest seed, against 1.72× for the ultrametric test.
//!
//! So: **report the quartet statistic as primary.** Keep the triple statistic —
//! it is cheaper, and a triple violation with *no* quartet violation is
//! informative in its own right, being the signature of rate variation on a tree
//! — but a triple certificate alone is not evidence of reticulation.
//!
//! The licence for reading it that way is a theorem, not an assumption:
//! `ANSE.GeneticaReticulata.fourPoint_of_ultrametric` proves every symmetric
//! ultrametric satisfies the four-point condition, so the two nulls are nested
//! and a quartet certificate refutes strictly more than a triple one.
//!
//! # Missing data: complete deletion within the subset
//!
//! Pairwise deletion gives different pairs different comparable-site counts, and
//! that alone can break the ultrametric inequality on an exactly clock-like
//! genealogy — the audit's finding 6. Measured on gap-blocked but strictly clonal
//! data: pairwise deletion calls **33–52 of 560 triples** violations that
//! complete deletion within the triple denies, at every seed. With no missing
//! data the two policies agree on all 560, as they must.
//!
//! [`subset::Alignment`] counts only the columns where every member of the
//! subset carries a base. All distances in a certificate then share one
//! denominator, so the obligations become comparisons between integers rather
//! than cross-products — a smaller trust surface *and* the fix, together. The
//! cost is power: at a 25% gap load roughly 1000–1700 of 4000 columns survive
//! within a triple.
//!
//! # What this does NOT solve
//!
//! Real sequence data is never *exactly* ultrametric. Finite sequence length
//! makes pairwise distances stochastic, so certificates exist under the null
//! too, with small margins. The theorem gives the exact boundary of the null; it
//! does not give the sampling distribution of the margin. Calibration is
//! therefore still required, and [`simulate`] provides it — but by a
//! **mechanistic generative model** of clonal evolution at the same sequence
//! length and divergence, rather than by permuting the observed data. That is a
//! real improvement in interpretability, not an escape from calibration, and it
//! should not be described as one.

pub mod certificate;
pub mod distance;
pub mod emit_lean;
pub mod exact;
pub mod rips;
pub mod simulate;
pub mod subset;

pub use certificate::{Certificate, CertificateReport, FourPointDefect};
pub use distance::{DistanceMatrix, DistanceModel};
pub use exact::{CountCertificate, ExactDistances, Ratio};
pub use rips::{CycleRepresentative, RipsPersistence};
pub use subset::{
    Alignment, CompleteTripleCertificate, QuadrupleCounts, QuartetCertificate, QuartetSearch,
    TripleCounts, TripleSearch,
};

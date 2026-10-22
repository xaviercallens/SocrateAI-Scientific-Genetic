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
//! The consequence for trust is the point of the design: **the scientific
//! conclusion does not depend on the correctness of this crate.** A bug here can
//! fail to find a certificate, or emit a triple that is not one; it cannot
//! produce a false conclusion that survives checking, because the certificate is
//! small enough to recheck by hand and its three inequalities are exactly the
//! hypotheses of the Lean theorem.
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
pub mod rips;
pub mod simulate;

pub use certificate::{Certificate, CertificateReport, FourPointDefect};
pub use distance::{DistanceMatrix, DistanceModel};
pub use rips::{CycleRepresentative, RipsPersistence};

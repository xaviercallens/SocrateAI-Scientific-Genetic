//! Emit certificates as a Lean 4 file for the kernel to typecheck.
//!
//! This is the step that moves the trust boundary. The Rust search produces
//! integers; this module writes them as `CountCertificate` instances; and
//! `lake env lean` decides whether each one is admissible. A certificate the
//! kernel rejects is a bug report about this crate, never a result.
//!
//! Nothing here needs to be trusted either. If the emitter writes the wrong
//! integers, the kernel either rejects the file (the inequalities fail) or
//! accepts a certificate about numbers that are not this dataset's — which the
//! accompanying provenance block makes checkable by a reader. The emitter
//! cannot make a false statement true.

use crate::exact::CountCertificate;

/// Render certificates as a self-contained Lean file.
///
/// `module_label` appears in the provenance comment, and `source` records where
/// the counts came from, so the emitted file says what it is about rather than
/// being an anonymous list of integers.
pub fn render(certs: &[CountCertificate], module_label: &str, source: &str) -> String {
    let mut out = String::new();
    out.push_str("/-\n");
    out.push_str("  Certificates emitted by `reticulation` for the kernel to check.\n\n");
    out.push_str(&format!("  Source      : {}\n", sanitise_label(source)));
    out.push_str(&format!("  Label       : {}\n", sanitise_label(module_label)));
    out.push_str(&format!("  Certificates: {}\n\n", certs.len()));
    out.push_str("  Each `CountCertificate` carries integer site counts. Its two `far_*`\n");
    out.push_str("  fields are the cross-multiplied inequalities of\n");
    out.push_str("  `ANSE.GeneticaReticulata.ratio_lt_ratio_iff`, so the kernel decides\n");
    out.push_str("  them exactly. A certificate that fails to typecheck is a defect in\n");
    out.push_str("  the search, not a negative result.\n\n");
    out.push_str("  Generated file -- do not edit by hand.\n-/\n");
    out.push_str("import Mathlib.Tactic\n");
    out.push_str("import ANSE.GeneticaReticulata\n\n");
    out.push_str("open ANSE.GeneticaReticulata\n\n");
    out.push_str("namespace Reticulata.Emitted\n\n");

    for (i, c) in certs.iter().enumerate() {
        // Labels come from FASTA headers, which are untrusted: `-/` would
        // close this comment and have the rest parsed as Lean.
        out.push_str(&format!(
            "/-- {} -- {} -- {}; exact margin {}/{}. -/\n",
            sanitise_label(&c.x), sanitise_label(&c.y), sanitise_label(&c.z),
            c.margin_num, c.margin_den
        ));
            // Fully qualified: a local definition cannot shadow the imported one.
        out.push_str(&format!(
            "def cert{i} : ANSE.GeneticaReticulata.CountCertificate (Fin 3) where\n"
        ));
        out.push_str("  x := 0\n  y := 1\n  z := 2\n");
        out.push_str(&format!("  n_xy := {}\n  d_xy := {}\n", c.n_xy, c.d_xy));
        out.push_str(&format!("  n_yz := {}\n  d_yz := {}\n", c.n_yz, c.d_yz));
        out.push_str(&format!("  n_xz := {}\n  d_xz := {}\n", c.n_xz, c.d_xz));
        out.push_str("  hd_xy := by norm_num\n");
        out.push_str("  hd_yz := by norm_num\n");
        out.push_str("  hd_xz := by norm_num\n");
        out.push_str("  far_xy := by norm_num\n");
        out.push_str("  far_yz := by norm_num\n\n");
        // Force the kernel to evaluate the consequence, not merely the fields.
        out.push_str(&format!(
            "theorem cert{i}_refutes : 0 < ANSE.GeneticaReticulata.margin cert{i} := "
        ));
        out.push_str(&format!(
            "ANSE.GeneticaReticulata.margin_pos cert{i}\n\n"
        ));
    }

    out.push_str("end Reticulata.Emitted\n\n");
    for i in 0..certs.len() {
        out.push_str(&format!(
            "#print axioms Reticulata.Emitted.cert{i}_refutes\n"
        ));
    }
    out
}

/// Encode a sequence for Lean: `A=0 C=1 G=2 T=3`, everything else `4`.
///
/// `4` is "not a definite base", which `ANSE.GeneticaReticulata.isBase`
/// excludes from both counts. Gaps, `N`, and IUPAC ambiguity codes all land
/// here, which is the pairwise-deletion convention.
pub fn encode(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        })
        .collect()
}

/// Replace characters that could escape a Lean doc comment.
///
/// A sequence label comes from a FASTA header, which is untrusted input. `-/`
/// closes a block comment, so an unsanitised label could end the comment early
/// and have the remainder parsed as Lean. Found by an adversarial audit of the
/// first version of this module.
pub fn sanitise_label(label: &str) -> String {
    label
        .chars()
        .map(|c| match c {
            // Both `-/` and `/-` need a slash, so removing slashes closes the
            // comment-escape without mangling ordinary hyphenated labels.
            '/' | '\\' | '\n' | '\r' | '`' => '_',
            c if c.is_control() => '_',
            c => c,
        })
        .collect()
}

/// Render a certificate **about the sequences**, so the kernel does the counting.
///
/// This is the form that closes the gap an audit found in [`render`]: there, the
/// kernel was handed six integers and proved an arithmetic fact about them, with
/// nothing tying those integers to an alignment. Here the three sequences are
/// embedded as lists, every count is proved by `rfl` — kernel reduction — and
/// the inference rule `sequences_refute_from_counts` is instantiated against
/// the computed counts.
///
/// What remains trusted is narrower and checkable by eye: that the embedded
/// lists are the dataset's sequences. A reader can diff them against the FASTA.
pub fn render_with_sequences(
    cert: &CountCertificate,
    sx: &[u8],
    sy: &[u8],
    sz: &[u8],
    source: &str,
) -> String {
    let list = |s: &[u8]| {
        let mut out = String::from("[");
        for (i, b) in s.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push_str(&b.to_string());
        }
        out.push(']');
        out
    };

    let mut out = String::new();
    out.push_str("/-\n");
    out.push_str("  A certificate ABOUT THE SEQUENCES, for the kernel to check.\n\n");
    out.push_str(&format!("  Source : {}\n", sanitise_label(source)));
    out.push_str(&format!(
        "  Triple : {} -- {} -- {}\n",
        sanitise_label(&cert.x),
        sanitise_label(&cert.y),
        sanitise_label(&cert.z)
    ));
    out.push_str("\n  The three sequences are embedded below. Every count is proved by\n");
    out.push_str("  `rfl`, so the KERNEL computes the differing and comparable sites --\n");
    out.push_str("  they are not supplied by the search. The only remaining trusted\n");
    out.push_str("  step is that these lists are the dataset's sequences, which a\n");
    out.push_str("  reader can diff against the alignment.\n\n");
    out.push_str("  Encoding: A=0 C=1 G=2 T=3, anything else 4 (not comparable).\n");
    out.push_str("  Generated file -- do not edit by hand.\n-/\n");
    out.push_str("import Mathlib.Tactic\n");
    out.push_str("import ANSE.GeneticaReticulata\n\n");
    out.push_str("set_option maxRecDepth 1000000\n\n");
    out.push_str("namespace Reticulata.EmittedSequences\n\n");
    out.push_str(&format!("def sx : List Nat := {}\n", list(sx)));
    out.push_str(&format!("def sy : List Nat := {}\n", list(sy)));
    out.push_str(&format!("def sz : List Nat := {}\n\n", list(sz)));

    // Fully qualified, so a local definition cannot shadow the imported one.
    out.push_str("/-- The kernel computes the counts and derives the violation. -/\n");
    out.push_str("theorem violates :\n");
    out.push_str("    max (ANSE.GeneticaReticulata.pDistance sx sy)\n");
    out.push_str("        (ANSE.GeneticaReticulata.pDistance sy sz)\n");
    out.push_str("      < ANSE.GeneticaReticulata.pDistance sx sz := by\n");
    out.push_str("  apply ANSE.GeneticaReticulata.sequences_refute_from_counts\n");
    out.push_str(&format!(
        "    sx sy sz {} {} {} {} {} {} <;>\n",
        cert.n_xy, cert.d_xy, cert.n_yz, cert.d_yz, cert.n_xz, cert.d_xz
    ));
    out.push_str("    first | rfl | norm_num\n\n");
    out.push_str("end Reticulata.EmittedSequences\n\n");
    out.push_str("#print axioms Reticulata.EmittedSequences.violates\n");
    out
}

/// Encode a sequence as a Lean `List Nat` literal.
fn list_literal(s: &[u8]) -> String {
    let mut out = String::from("[");
    for (i, b) in s.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(&b.to_string());
    }
    out.push(']');
    out
}

/// Render a **complete-deletion** triple certificate about the sequences.
///
/// The difference from [`render_with_sequences`] is the deletion policy, and it
/// is the whole point: here the kernel counts only the columns at which all
/// three sequences carry a base, so the three distances share one denominator
/// and the obligations are comparisons between integers rather than
/// cross-products. That removes the confound of the audit's finding 6, which
/// `ANSE.GeneticaReticulata.pairwise_deletion_manufactures_what_complete_deletion_denies`
/// shows is real.
pub fn render_complete_triple(
    cert: &crate::subset::CompleteTripleCertificate,
    sx: &[u8],
    sy: &[u8],
    sz: &[u8],
    source: &str,
) -> String {
    let q = "ANSE.GeneticaReticulata";
    let mut out = String::new();
    out.push_str("/-\n");
    out.push_str("  A COMPLETE-DELETION triple certificate, for the kernel to check.\n\n");
    out.push_str(&format!("  Source : {}\n", sanitise_label(source)));
    out.push_str(&format!(
        "  Triple : {} -- {} -- {}\n",
        sanitise_label(&cert.x),
        sanitise_label(&cert.y),
        sanitise_label(&cert.z)
    ));
    out.push_str(&format!(
        "  Counts : {} / {} / {} differing over {} columns common to all three\n",
        cert.n_xy, cert.n_yz, cert.n_xz, cert.common
    ));
    out.push_str("\n  Only the columns at which ALL THREE sequences carry a definite base\n");
    out.push_str("  are counted, so the three distances share one denominator and the\n");
    out.push_str("  certificate's obligations are integer comparisons -- no\n");
    out.push_str("  cross-multiplication. Under pairwise deletion, two pairs are\n");
    out.push_str("  compared over different columns, which can break the ultrametric\n");
    out.push_str("  inequality on an exactly clock-like genealogy.\n\n");
    out.push_str("  Encoding: A=0 C=1 G=2 T=3, anything else 4 (not comparable).\n");
    out.push_str("  Generated file -- do not edit by hand.\n-/\n");
    out.push_str("import Mathlib.Tactic\n");
    out.push_str("import ANSE.GeneticaReticulataTree\n\n");
    out.push_str("set_option maxRecDepth 1000000\n\n");
    out.push_str("namespace Reticulata.EmittedComplete\n\n");
    out.push_str(&format!("def sx : List Nat := {}\n", list_literal(sx)));
    out.push_str(&format!("def sy : List Nat := {}\n", list_literal(sy)));
    out.push_str(&format!("def sz : List Nat := {}\n\n", list_literal(sz)));
    out.push_str("/-- The kernel computes the triple-complete counts and derives the\n");
    out.push_str("violation of the strong triangle inequality. -/\n");
    out.push_str("theorem violates :\n");
    out.push_str(&format!(
        "    max ({q}.ratioN ({q}.diff3xy sx sy sz) ({q}.comparable3 sx sy sz))\n"
    ));
    out.push_str(&format!(
        "        ({q}.ratioN ({q}.diff3yz sx sy sz) ({q}.comparable3 sx sy sz))\n"
    ));
    out.push_str(&format!(
        "      < {q}.ratioN ({q}.diff3xz sx sy sz) ({q}.comparable3 sx sy sz) := by\n"
    ));
    out.push_str(&format!(
        "  apply {q}.complete_triple_refute_from_counts\n"
    ));
    out.push_str(&format!(
        "    sx sy sz {} {} {} {} <;>\n",
        cert.n_xy, cert.n_yz, cert.n_xz, cert.common
    ));
    out.push_str("    first | rfl | norm_num\n\n");
    out.push_str("end Reticulata.EmittedComplete\n\n");
    out.push_str("#print axioms Reticulata.EmittedComplete.violates\n");
    out
}

/// Render a **quartet** certificate about the sequences: the four-point test.
///
/// This is the certificate that survives relaxation of the molecular clock. The
/// kernel counts the columns common to all four sequences, forms the three
/// pairing sums, and checks that the one named here strictly exceeds both
/// others. By `fourPoint_of_ultrametric` that refutes ultrametricity too, so one
/// quartet certificate refutes both nulls.
pub fn render_quartet(
    cert: &crate::subset::QuartetCertificate,
    sa: &[u8],
    sb: &[u8],
    sc: &[u8],
    se: &[u8],
    source: &str,
) -> String {
    let q = "ANSE.GeneticaReticulata";
    let mut out = String::new();
    out.push_str("/-\n");
    out.push_str("  A QUARTET certificate: the four-point condition refuted.\n\n");
    out.push_str(&format!("  Source   : {}\n", sanitise_label(source)));
    out.push_str(&format!(
        "  Quadruple: ({} , {}) | ({} , {})\n",
        sanitise_label(&cert.a),
        sanitise_label(&cert.b),
        sanitise_label(&cert.c),
        sanitise_label(&cert.e)
    ));
    out.push_str(&format!(
        "  Sums     : dominant {} vs {} and {}, over {} common columns\n",
        cert.n_ab + cert.n_ce,
        cert.n_ac + cert.n_be,
        cert.n_ae + cert.n_bc,
        cert.common
    ));
    out.push_str("\n  The four-point condition says a weighted tree leaves the two largest\n");
    out.push_str("  pairing sums EQUAL, whatever the branch lengths. It therefore does\n");
    out.push_str("  not assume a molecular clock, and a violation cannot be explained\n");
    out.push_str("  by rate variation. That a metric satisfying it on every quadruple\n");
    out.push_str("  is realised by a tree is Buneman's theorem: CITED, not proved.\n\n");
    out.push_str("  Encoding: A=0 C=1 G=2 T=3, anything else 4 (not comparable).\n");
    out.push_str("  Generated file -- do not edit by hand.\n-/\n");
    out.push_str("import Mathlib.Tactic\n");
    out.push_str("import ANSE.GeneticaReticulataTree\n\n");
    out.push_str("set_option maxRecDepth 1000000\n\n");
    out.push_str("namespace Reticulata.EmittedQuartet\n\n");
    out.push_str(&format!("def sa : List Nat := {}\n", list_literal(sa)));
    out.push_str(&format!("def sb : List Nat := {}\n", list_literal(sb)));
    out.push_str(&format!("def sc : List Nat := {}\n", list_literal(sc)));
    out.push_str(&format!("def se : List Nat := {}\n\n", list_literal(se)));
    out.push_str("/-- The kernel computes all six counts and the common denominator,\n");
    out.push_str("then checks that the named pairing strictly dominates. -/\n");
    out.push_str("theorem violates :\n");
    out.push_str(&format!(
        "    max ({q}.ratioN ({q}.diff4ac sa sb sc se) ({q}.comparable4 sa sb sc se)\n"
    ));
    out.push_str(&format!(
        "          + {q}.ratioN ({q}.diff4be sa sb sc se) ({q}.comparable4 sa sb sc se))\n"
    ));
    out.push_str(&format!(
        "        ({q}.ratioN ({q}.diff4ae sa sb sc se) ({q}.comparable4 sa sb sc se)\n"
    ));
    out.push_str(&format!(
        "          + {q}.ratioN ({q}.diff4bc sa sb sc se) ({q}.comparable4 sa sb sc se))\n"
    ));
    out.push_str(&format!(
        "      < {q}.ratioN ({q}.diff4ab sa sb sc se) ({q}.comparable4 sa sb sc se)\n"
    ));
    out.push_str(&format!(
        "          + {q}.ratioN ({q}.diff4ce sa sb sc se) ({q}.comparable4 sa sb sc se) := by\n"
    ));
    out.push_str(&format!("  apply {q}.quartet_refute_from_counts\n"));
    out.push_str(&format!(
        "    sa sb sc se {} {} {} {} {} {} {}\n",
        cert.n_ab, cert.n_ce, cert.n_ac, cert.n_be, cert.n_ae, cert.n_bc, cert.common
    ));
    // The anonymous constructor is tried first: it fails on the count equations
    // and on `0 < D`, and is the only thing that closes `fourPointViolates`.
    out.push_str("  all_goals first\n");
    out.push_str("    | exact ⟨by norm_num, by norm_num⟩\n");
    out.push_str("    | rfl\n");
    out.push_str("    | norm_num\n\n");
    out.push_str("end Reticulata.EmittedQuartet\n\n");
    out.push_str("#print axioms Reticulata.EmittedQuartet.violates\n");
    out
}

/// Corrupt a complete-deletion triple certificate, for the rejection control.
///
/// Swapping the far and near numerators inverts `n_xy < n_xz`, so the
/// obligation is false by construction and no amount of kernel work can accept
/// it. The swap cannot accidentally yield a still-valid certificate.
pub fn corrupt_complete_triple(
    cert: &crate::subset::CompleteTripleCertificate,
) -> crate::subset::CompleteTripleCertificate {
    let mut bad = cert.clone();
    std::mem::swap(&mut bad.n_xz, &mut bad.n_xy);
    bad
}

/// Corrupt a quartet certificate, for the rejection control.
///
/// Exchanging the dominant pairing with one of the dominated ones makes the
/// claimed domination false: the new `n_ab + n_ce` is the old second-largest sum
/// and the new `n_ac + n_be` is the old largest, so `fourPointViolates` fails on
/// its first conjunct. Invalid by construction, never by luck.
pub fn corrupt_quartet(
    cert: &crate::subset::QuartetCertificate,
) -> crate::subset::QuartetCertificate {
    let mut bad = cert.clone();
    std::mem::swap(&mut bad.n_ab, &mut bad.n_ac);
    std::mem::swap(&mut bad.n_ce, &mut bad.n_be);
    bad
}

/// Deliberately corrupt one certificate, for the rejection control.
///
/// A verifier that accepts everything is worthless, so the proof of concept
/// must show the kernel *rejecting* something. Swapping the far and near
/// numerators inverts the inequality the certificate asserts, so `far_xy`
/// becomes false and the file cannot typecheck.
pub fn corrupt(cert: &CountCertificate) -> CountCertificate {
    let mut bad = cert.clone();
    std::mem::swap(&mut bad.n_xz, &mut bad.n_xy);
    std::mem::swap(&mut bad.d_xz, &mut bad.d_xy);
    bad
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> CountCertificate {
        CountCertificate {
            x: "a".into(),
            y: "m".into(),
            z: "b".into(),
            indices: [0, 1, 2],
            n_xy: 69,
            d_xy: 4000,
            n_yz: 70,
            d_yz: 4000,
            n_xz: 131,
            d_xz: 4000,
            margin_num: 61,
            margin_den: 4000,
        }
    }

    #[test]
    fn rendered_file_has_the_structure_lean_expects() {
        let text = render(&[sample()], "test", "unit test");
        assert!(text.contains("import ANSE.GeneticaReticulata"));
        assert!(text.contains(
            "def cert0 : ANSE.GeneticaReticulata.CountCertificate (Fin 3) where"
        ));
        assert!(text.contains("n_xz := 131"));
        assert!(text.contains("far_xy := by norm_num"));
        // The consequence must be forced, not just the fields.
        assert!(text.contains("ANSE.GeneticaReticulata.margin_pos cert0"));
        assert!(text.contains("#print axioms Reticulata.Emitted.cert0_refutes"));
    }

    #[test]
    fn rendered_file_records_its_provenance() {
        let text = render(&[sample()], "poc-label", "simulated recombinant");
        assert!(text.contains("poc-label"));
        assert!(text.contains("simulated recombinant"));
        assert!(text.contains("do not edit by hand"));
    }

    #[test]
    fn corruption_really_breaks_the_inequality() {
        // The rejection control must actually be invalid, or the PoC proves
        // nothing about the kernel.
        let good = sample();
        assert!(good.validate().is_ok());
        let bad = corrupt(&good);
        assert!(
            bad.validate().is_err(),
            "the corrupted certificate must be invalid, else the kernel \
             rejection test is vacuous"
        );
    }

    #[test]
    fn empty_input_renders_a_valid_but_empty_module() {
        let text = render(&[], "empty", "none");
        assert!(text.contains("namespace Reticulata.Emitted"));
        assert!(!text.contains("def cert0"));
    }
}

#[cfg(test)]
mod audit_fix_tests {
    use super::*;

    /// A FASTA header containing `-/` must not be able to close the Lean doc
    /// comment. Found by an adversarial audit; the label reached the emitter
    /// unescaped.
    #[test]
    fn label_cannot_escape_a_lean_comment() {
        for hostile in ["a-/axiom", "x/-nested", "back`tick", "line\nbreak", "sl/ash"] {
            let clean = sanitise_label(hostile);
            assert!(!clean.contains("-/"), "{hostile} -> {clean}");
            assert!(!clean.contains("/-"), "{hostile} -> {clean}");
            assert!(!clean.contains('/'), "{hostile} -> {clean}");
            assert!(!clean.contains('`'), "{hostile} -> {clean}");
            assert!(!clean.contains('\n'), "{hostile} -> {clean}");
        }
        // Ordinary labels must survive legibly.
        assert_eq!(sanitise_label("clonal_003"), "clonal_003");
        assert_eq!(sanitise_label("poc-label"), "poc-label");
    }

    /// The hostile label must be neutralised in the rendered file too, not just
    /// by the helper.
    #[test]
    fn rendered_file_neutralises_hostile_labels() {
        let c = CountCertificate {
            x: "a-/axiom".into(),
            y: "m".into(),
            z: "b".into(),
            indices: [0, 1, 2],
            n_xy: 1, d_xy: 4, n_yz: 1, d_yz: 4, n_xz: 3, d_xz: 4,
            margin_num: 2, margin_den: 4,
        };
        let text = render(&[c], "lbl", "src");
        // Exactly one `-/` for the header plus one per certificate doc comment.
        // A hostile label adding a third would mean it escaped.
        assert_eq!(
            text.matches("-/").count(),
            2,
            "one header close plus one doc-comment close, and no more"
        );
        // Every comment opener is matched, so nothing after a label is code.
        assert_eq!(text.matches("/-").count(), text.matches("-/").count());
    }

    /// Emitted names must be fully qualified, so a local definition in a
    /// hand-edited or third-party file cannot shadow the imported structure.
    #[test]
    fn emitted_names_are_fully_qualified() {
        let c = CountCertificate {
            x: "a".into(), y: "m".into(), z: "b".into(),
            indices: [0, 1, 2],
            n_xy: 1, d_xy: 4, n_yz: 1, d_yz: 4, n_xz: 3, d_xz: 4,
            margin_num: 2, margin_den: 4,
        };
        let text = render(&[c], "lbl", "src");
        assert!(text.contains("ANSE.GeneticaReticulata.CountCertificate"));
        assert!(text.contains("ANSE.GeneticaReticulata.margin_pos"));
        assert!(text.contains("0 < ANSE.GeneticaReticulata.margin cert0"));
    }

    /// The sequence-level renderer embeds the data, so the kernel does the
    /// counting. This is the fix for the audit's headline finding.
    #[test]
    fn sequence_renderer_embeds_the_data_and_qualifies_names() {
        let c = CountCertificate {
            x: "a".into(), y: "m".into(), z: "b".into(),
            indices: [0, 1, 2],
            n_xy: 2, d_xy: 4, n_yz: 2, d_yz: 4, n_xz: 4, d_xz: 4,
            margin_num: 2, margin_den: 4,
        };
        let sx = encode(b"AAAA");
        let sy = encode(b"AACC");
        let sz = encode(b"CCCC");
        let text = render_with_sequences(&c, &sx, &sy, &sz, "unit test");
        assert!(text.contains("def sx : List Nat := [0,0,0,0]"));
        assert!(text.contains("def sz : List Nat := [1,1,1,1]"));
        assert!(text.contains("ANSE.GeneticaReticulata.sequences_refute_from_counts"));
        assert!(text.contains("ANSE.GeneticaReticulata.pDistance"));
        // The counts must appear as arguments, for the kernel to check by rfl.
        assert!(text.contains("sx sy sz 2 4 2 4 4 4"));
    }

    /// The complete-deletion renderer must name the tree module's theorem and
    /// pass the single common denominator, not three per-pair ones.
    #[test]
    fn complete_triple_renderer_passes_one_denominator() {
        let aln = crate::subset::Alignment::new(
            vec!["x".into(), "y".into(), "z".into()],
            vec![b"AAAA".to_vec(), b"AACC".to_vec(), b"CCCC".to_vec()],
        )
        .unwrap();
        let search = crate::subset::find_complete_triple_certificates(&aln, 1).unwrap();
        let c = &search.certificates[0];
        let text = render_complete_triple(
            c,
            &encode(b"AAAA"),
            &encode(b"AACC"),
            &encode(b"CCCC"),
            "unit test",
        );
        assert!(text.contains("import ANSE.GeneticaReticulataTree"));
        assert!(text.contains("ANSE.GeneticaReticulata.complete_triple_refute_from_counts"));
        assert!(text.contains("ANSE.GeneticaReticulata.diff3xy sx sy sz"));
        // Three numerators and ONE denominator: 2 2 4 4.
        assert!(text.contains("sx sy sz 2 2 4 4"), "{text}");
        assert!(text.contains("def sy : List Nat := [0,0,1,1]"));
    }

    /// The quartet renderer must emit the six counts in the Lean field order and
    /// discharge `fourPointViolates` with the anonymous constructor.
    #[test]
    fn quartet_renderer_emits_six_counts_and_the_violation() {
        let aln = crate::subset::Alignment::new(
            vec!["p".into(), "q".into(), "m1".into(), "m2".into()],
            vec![
                b"AAAA".to_vec(),
                b"CCCC".to_vec(),
                b"CCAA".to_vec(),
                b"AACC".to_vec(),
            ],
        )
        .unwrap();
        let search = crate::subset::find_quartet_certificates(&aln, 1, 100).unwrap();
        let c = &search.certificates[0];
        let text = render_quartet(
            c,
            &encode(b"AAAA"),
            &encode(b"CCCC"),
            &encode(b"CCAA"),
            &encode(b"AACC"),
            "unit test",
        );
        assert!(text.contains("ANSE.GeneticaReticulata.quartet_refute_from_counts"));
        assert!(text.contains("ANSE.GeneticaReticulata.comparable4 sa sb sc se"));
        // n_ab n_ce n_ac n_be n_ae n_bc D = 4 4 2 2 2 2 4.
        assert!(text.contains("sa sb sc se 4 4 2 2 2 2 4"), "{text}");
        assert!(text.contains("exact ⟨by norm_num, by norm_num⟩"));
        assert!(text.contains("Buneman's theorem: CITED, not proved"));
    }

    /// Both corruptions must be invalid by construction, or the rejection
    /// controls that use them prove nothing about the kernel.
    #[test]
    fn the_new_corruptions_really_break_their_obligations() {
        let aln = crate::subset::Alignment::new(
            vec!["p".into(), "q".into(), "m1".into(), "m2".into()],
            vec![
                b"AAAA".to_vec(),
                b"CCCC".to_vec(),
                b"CCAA".to_vec(),
                b"AACC".to_vec(),
            ],
        )
        .unwrap();
        let quartet = &crate::subset::find_quartet_certificates(&aln, 1, 100)
            .unwrap()
            .certificates[0];
        assert!(quartet.validate().is_ok());
        assert!(corrupt_quartet(quartet).validate().is_err());

        let triple_aln = crate::subset::Alignment::new(
            vec!["x".into(), "y".into(), "z".into()],
            vec![b"AAAA".to_vec(), b"AACC".to_vec(), b"CCCC".to_vec()],
        )
        .unwrap();
        let triple = &crate::subset::find_complete_triple_certificates(&triple_aln, 1)
            .unwrap()
            .certificates[0];
        assert!(triple.validate().is_ok());
        assert!(corrupt_complete_triple(triple).validate().is_err());
    }

    /// Gaps and ambiguity codes encode to 4, which Lean's `isBase` excludes.
    #[test]
    fn encoding_maps_non_bases_to_the_incomparable_code() {
        assert_eq!(encode(b"ACGT"), vec![0, 1, 2, 3]);
        assert_eq!(encode(b"acgt"), vec![0, 1, 2, 3]);
        for other in [b"N".as_slice(), b"-", b"U", b"R", b"?"] {
            assert_eq!(encode(other), vec![4], "{:?} must be incomparable", other);
        }
    }
}

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

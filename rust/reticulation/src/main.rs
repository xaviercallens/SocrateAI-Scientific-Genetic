//! `reticulation` — search a sequence population for proved reticulation certificates.
//!
//! Subcommands:
//!
//! ```text
//! reticulation controls [--length N] [--generations G] [--mu M] [--seed S]
//!     Run the negative (clonal) and positive (recombinant) controls and emit
//!     JSON. Exits non-zero if either control fails, so a broken detector
//!     cannot report a result.
//!
//! reticulation fasta <FILE> [--model p|jc] [--top K]
//!     Read an aligned FASTA, search every triple, emit certificates and the
//!     Vietoris-Rips H1 intervals with their representative cycles.
//! ```
//!
//! Every certificate is rechecked against the hypotheses of
//! `ANSE.ReticulationCertificate.certificate_refutes_ultrametric` before it is
//! printed.

use std::process::ExitCode;

use reticulation::certificate::{add_four_point, find_certificates};
use reticulation::distance::{DistanceMatrix, DistanceModel};
use reticulation::rips::vietoris_rips_h1;
use reticulation::simulate::{clonal_tree, with_recombination};

fn median(mut v: Vec<f64>) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(f64::total_cmp);
    v[v.len() / 2]
}

/// Median off-diagonal distance: a defensible default Rips threshold.
fn default_scale(dm: &DistanceMatrix) -> f64 {
    let n = dm.len();
    let mut all = Vec::with_capacity(n * (n - 1) / 2);
    for i in 0..n {
        for j in (i + 1)..n {
            all.push(dm.get(i, j));
        }
    }
    median(all)
}

fn parse_flag<T: std::str::FromStr>(args: &[String], name: &str, default: T) -> T {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn read_fasta(path: &str) -> Result<(Vec<String>, Vec<Vec<u8>>), String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let mut labels = Vec::new();
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(name) = line.strip_prefix('>') {
            labels.push(name.split_whitespace().next().unwrap_or(name).to_string());
            seqs.push(Vec::new());
        } else {
            match seqs.last_mut() {
                Some(s) => s.extend_from_slice(line.as_bytes()),
                None => return Err("sequence data before the first '>' header".into()),
            }
        }
    }
    if labels.is_empty() {
        return Err(format!("{path}: no FASTA records found"));
    }
    Ok((labels, seqs))
}

fn run_controls(args: &[String]) -> Result<bool, String> {
    let length: usize = parse_flag(args, "--length", 4000);
    let generations: u32 = parse_flag(args, "--generations", 4);
    let mu: f64 = parse_flag(args, "--mu", 0.004);
    let seed: u64 = parse_flag(args, "--seed", 2026);

    // Negative control: strictly clonal.
    let clonal = clonal_tree(length, generations, mu, seed);
    let dm_clonal = clonal.distances(DistanceModel::PDistance)?;
    let clonal_report = find_certificates(&dm_clonal, 3);

    // Positive control: the same population plus one recombinant.
    let base = clonal_tree(length, generations, mu, seed);
    let (recombined, recombinant_label) = with_recombination(base, 0.5, seed ^ 0xABCD);
    let dm_rec = recombined.distances(DistanceModel::PDistance)?;
    let rec_report = find_certificates(&dm_rec, 5);

    let typical = default_scale(&dm_clonal);
    let negative_passes = clonal_report.ultrametric_defect < 0.5 * typical.max(1e-12);
    let positive_passes = rec_report.ultrametric_defect
        > 2.0 * clonal_report.ultrametric_defect.max(1e-12)
        && rec_report
            .certificates
            .iter()
            .any(|c| c.y == recombinant_label);

    let out = serde_json::json!({
        "description": "Controls for the reticulation certificate search. The \
                        negative control is strictly clonal descent under a \
                        strict clock; the positive control adds one recombinant \
                        between maximally diverged parents.",
        "parameters": {
            "alignment_length": length, "generations": generations,
            "mu_per_generation": mu, "seed": seed,
            "model": "p-distance",
        },
        "negative_control_clonal": {
            "n_sequences": clonal_report.n_sequences,
            "ultrametric_defect": clonal_report.ultrametric_defect,
            "violating_triples": clonal_report.violating_triples,
            "triples_examined": clonal_report.triples_examined,
            "typical_distance": typical,
            "passes": negative_passes,
            "criterion": "defect < 0.5 * median pairwise distance",
            "note": "finite sequences make clonal data only approximately \
                     ultrametric, so a small defect here is expected and is the \
                     noise floor the positive control must clear",
        },
        "positive_control_recombinant": {
            "n_sequences": rec_report.n_sequences,
            "ultrametric_defect": rec_report.ultrametric_defect,
            "violating_triples": rec_report.violating_triples,
            "recombinant_label": recombinant_label,
            "recombinant_is_a_midpoint": rec_report
                .certificates.iter().any(|c| c.y == recombinant_label),
            "all_certificates_validated": rec_report.all_certificates_validated,
            "top_certificates": rec_report.certificates,
            "passes": positive_passes,
            "criterion": "defect > 2x clonal defect AND the recombinant is the \
                          midpoint of a top certificate",
        },
        "verdict": {
            "negative_passes": negative_passes,
            "positive_passes": positive_passes,
            "lean_module": "ANSE.ReticulationCertificate",
            "inference_rule": "certificate_refutes_ultrametric",
            "control_theorem": "ultrametric_admits_no_certificate",
        }
    });
    println!("{}", serde_json::to_string_pretty(&out).unwrap());

    if !negative_passes {
        eprintln!(
            "FAIL negative control: clonal defect {} is not small against the \
             typical distance {typical}. The detector is reading sampling noise \
             and no positive result would be usable.",
            clonal_report.ultrametric_defect
        );
        return Ok(false);
    }
    if !positive_passes {
        eprintln!(
            "FAIL positive control: a recombinant between diverged parents was \
             not detected (defect {} against clonal {}). The detector has no \
             power and a null result would be uninformative.",
            rec_report.ultrametric_defect, clonal_report.ultrametric_defect
        );
        return Ok(false);
    }
    Ok(true)
}

fn run_fasta(args: &[String]) -> Result<bool, String> {
    let path = args
        .iter()
        .skip(1)
        .find(|a| !a.starts_with("--"))
        .ok_or("usage: reticulation fasta <FILE> [--model p|jc] [--top K]")?;
    let model = match parse_flag::<String>(args, "--model", "p".into()).as_str() {
        "jc" => DistanceModel::JukesCantor,
        _ => DistanceModel::PDistance,
    };
    let top: usize = parse_flag(args, "--top", 10);

    let (labels, seqs) = read_fasta(path)?;
    let dm = DistanceMatrix::from_alignment(labels, &seqs, model)?;
    let mut report = find_certificates(&dm, top);
    add_four_point(&mut report, &dm, top, 48);

    let scale = parse_flag(args, "--scale", default_scale(&dm));
    let persistence = vietoris_rips_h1(&dm, scale);

    let out = serde_json::json!({
        "input": path,
        "model": model,
        "alignment_length": dm.alignment_length,
        "certificates": report,
        "rips_persistence": persistence,
        "interpretation": {
            "proved": "A validated certificate proves the distance matrix is not \
                       an ultrametric, i.e. not a strict clock on a clonal \
                       genealogy (ANSE.ReticulationCertificate.\
                       certificate_refutes_ultrametric).",
            "cited_not_proved": "That tree metrics likewise have acyclic \
                                 Vietoris-Rips filtrations is cited, not proved \
                                 here, so four-point defects are suggestive only.",
            "not_established": "Finite sequences make clonal data only \
                                approximately ultrametric. A certificate alone \
                                is not evidence of recombination unless its \
                                margin exceeds the clonal noise floor; run \
                                `reticulation controls` at the same length and \
                                divergence to obtain it.",
        }
    });
    println!("{}", serde_json::to_string_pretty(&out).unwrap());
    Ok(true)
}

/// Typecheck an emitted Lean file with the project toolchain.
///
/// Returns `(accepted, first_error_line)`. We invoke `lake env lean` from the
/// repository's `formal` directory, which is the same path the project's Lean
/// gate uses, so the kernel and the axiom whitelist are the ones already in use.
fn kernel_check(formal_dir: &std::path::Path, file: &std::path::Path) -> (bool, String) {
    let out = std::process::Command::new("lake")
        .args(["env", "lean", &file.to_string_lossy()])
        .current_dir(formal_dir)
        .output();
    match out {
        Ok(o) => {
            let text = format!(
                "{}{}",
                String::from_utf8_lossy(&o.stdout),
                String::from_utf8_lossy(&o.stderr)
            );
            let first_error = text
                .lines()
                .find(|l| l.contains("error"))
                .unwrap_or("")
                .trim()
                .to_string();
            // Accepted only if the toolchain exited 0, reported no error, and
            // every axiom the file depends on is whitelisted. Checking only for
            // `sorryAx` was a gap an audit found: a file could declare its own
            // axiom and still pass. Compilation alone is never sufficient.
            let whitelist = ["propext", "Classical.choice", "Quot.sound"];
            let mut axioms_ok = true;
            for line in text.lines() {
                if let Some(rest) = line.split("depends on axioms:").nth(1) {
                    for ax in rest
                        .trim()
                        .trim_start_matches('[')
                        .trim_end_matches(']')
                        .split(',')
                    {
                        let ax = ax.trim();
                        if !ax.is_empty() && !whitelist.contains(&ax) {
                            axioms_ok = false;
                        }
                    }
                }
            }
            let accepted = o.status.success()
                && !text.contains("error")
                && !text.contains("sorryAx")
                && axioms_ok;
            (accepted, first_error)
        }
        Err(e) => (false, format!("could not run lake: {e}")),
    }
}

/// The proof of concept: Rust searches, the Lean kernel decides.
///
/// Runs the full chain on simulated data with known ground truth, then runs the
/// rejection control. Both directions are required: a verifier that accepts
/// everything would prove nothing.
fn run_poc(args: &[String]) -> Result<bool, String> {
    use reticulation::emit_lean;
    use reticulation::exact::{find_count_certificates, ExactDistances};

    let length: usize = parse_flag(args, "--length", 4000);
    let generations: u32 = parse_flag(args, "--generations", 4);
    let mu: f64 = parse_flag(args, "--mu", 0.004);
    let seed: u64 = parse_flag(args, "--seed", 2026);
    let top: usize = parse_flag(args, "--top", 3);
    let formal_dir: String = parse_flag(
        args,
        "--formal-dir",
        "../../../../../formal".to_string(),
    );
    let out_dir: String = parse_flag(args, "--out", "/tmp".to_string());

    let formal = std::path::Path::new(&formal_dir)
        .canonicalize()
        .map_err(|e| format!("--formal-dir {formal_dir}: {e}"))?;

    // --- clonal baseline, for the noise floor ---
    let clonal = clonal_tree(length, generations, mu, seed);
    let dm_clonal = ExactDistances::from_alignment(clonal.labels.clone(), &clonal.sequences)?;
    let (clonal_certs, clonal_violating, _) = find_count_certificates(&dm_clonal, 1);
    let clonal_floor = clonal_certs.first().map(|c| c.margin_f64()).unwrap_or(0.0);
    // Keep the clonal margin as an exact rational: the PASS/FAIL gate below
    // must not be decided in floating point. An audit found the original gate
    // compared f64 ratios, which is a verdict depending on rounding.
    let clonal_exact = clonal_certs
        .first()
        .map(|c| (c.margin_num as i128, c.margin_den as i128));

    // --- the planted recombinant ---
    let base = clonal_tree(length, generations, mu, seed);
    let (recombined, recombinant_label) = with_recombination(base, 0.5, seed ^ 0xABCD);
    let dist = ExactDistances::from_alignment(recombined.labels.clone(), &recombined.sequences)?;
    let (certs, violating, skipped) = find_count_certificates(&dist, top);
    if certs.is_empty() {
        return Err("no certificate found on the positive control".into());
    }

    // --- emit and kernel-check the genuine certificates ---
    let good_path = std::path::Path::new(&out_dir).join("reticulata_emitted.lean");
    std::fs::write(
        &good_path,
        emit_lean::render(
            &certs,
            "poc-positive-control",
            &format!(
                "simulated clonal tree ({length} nt, {generations} generations, \
                 mu={mu}, seed={seed}) plus one recombinant"
            ),
        ),
    )
    .map_err(|e| format!("{}: {e}", good_path.display()))?;
    let (good_accepted, good_err) = kernel_check(&formal, &good_path);

    // --- the audit fix: a certificate ABOUT THE SEQUENCES ---
    // An adversarial audit found that the integer-only certificates above are
    // not about the data: the kernel proves an arithmetic fact about six
    // numbers, with nothing tying them to an alignment. Here the three
    // sequences are embedded and the kernel computes the counts itself by
    // `rfl`, which is what makes the certificate a statement about the data.
    let top = &certs[0];
    let [ix, iy, iz] = top.indices;
    let seq_path = std::path::Path::new(&out_dir).join("reticulata_sequences.lean");
    std::fs::write(
        &seq_path,
        emit_lean::render_with_sequences(
            top,
            &emit_lean::encode(&recombined.sequences[ix]),
            &emit_lean::encode(&recombined.sequences[iy]),
            &emit_lean::encode(&recombined.sequences[iz]),
            &format!("simulated recombinant, seed {seed}"),
        ),
    )
    .map_err(|e| format!("{}: {e}", seq_path.display()))?;
    let seq_start = std::time::Instant::now();
    let (seq_accepted, seq_err) = kernel_check(&formal, &seq_path);
    let seq_secs = seq_start.elapsed().as_secs_f64();

    // --- the rejection control: a deliberately corrupted certificate ---
    let bad = emit_lean::corrupt(&certs[0]);
    let bad_is_invalid = bad.validate().is_err();
    let bad_path = std::path::Path::new(&out_dir).join("reticulata_corrupted.lean");
    std::fs::write(
        &bad_path,
        emit_lean::render(
            std::slice::from_ref(&bad),
            "poc-rejection-control",
            "deliberately corrupted: far and near distances swapped",
        ),
    )
    .map_err(|e| format!("{}: {e}", bad_path.display()))?;
    let (bad_accepted, bad_err) = kernel_check(&formal, &bad_path);

    let recombinant_is_midpoint = certs.iter().any(|c| c.y == recombinant_label);
    let ratio = if clonal_floor > 0.0 {
        certs[0].margin_f64() / clonal_floor
    } else {
        f64::INFINITY
    };

    let accept_passes = good_accepted && seq_accepted;
    let reject_passes = bad_is_invalid && !bad_accepted;
    // Exact: margin_pos / margin_clonal > 2  <=>  num_p * den_c > 2 * num_c * den_p
    let exceeds_twice_floor = match clonal_exact {
        Some((cn, cd)) if cn > 0 => {
            (certs[0].margin_num as i128) * cd > 2 * cn * (certs[0].margin_den as i128)
        }
        // No clonal violation at all: any positive margin exceeds the floor.
        _ => certs[0].margin_num > 0,
    };
    let discriminates = exceeds_twice_floor && recombinant_is_midpoint;

    let report = serde_json::json!({
        "proof_of_concept": "Genetica Reticulata Phase 1 -- the trust boundary \
                             moved into the Lean kernel. Rust emits integer \
                             count certificates; `lake env lean` decides them.",
        "parameters": {
            "alignment_length": length, "generations": generations,
            "mu_per_generation": mu, "seed": seed, "top_k": top,
        },
        "exactness": {
            "arithmetic": "integer cross-products in i128; no floating point on \
                           the evidence path",
            "denominators_per_pair": true,
            "note": "f64 appears only in the ranking display below",
        },
        "clonal_baseline": {
            "violating_triples": clonal_violating,
            "largest_margin": clonal_floor,
            "note": "finite sequences are never exactly ultrametric, so clonal \
                     data yields valid certificates with small margins; this is \
                     the noise floor, not a detection",
        },
        "positive_control": {
            "n_sequences": dist.len(),
            "violating_triples": violating,
            "skipped_uninformative_triples": skipped,
            "recombinant_label": recombinant_label,
            "recombinant_is_a_certificate_midpoint": recombinant_is_midpoint,
            "largest_margin_exact": format!("{}/{}", certs[0].margin_num, certs[0].margin_den),
            "largest_margin": certs[0].margin_f64(),
            "margin_ratio_to_clonal_floor": ratio,
            "exceeds_twice_clonal_floor_exact": exceeds_twice_floor,
            "gate_arithmetic": "exact integer cross-product; the displayed ratio \
                                is f64 for reading only",
            "certificates": certs,
        },
        "kernel_counts_the_data": {
            "emitted_file": seq_path.display().to_string(),
            "accepted": seq_accepted,
            "first_error": seq_err,
            "seconds": seq_secs,
            "lean_theorem": "sequences_refute_from_counts",
            "what_changed": "the three sequences are embedded and every count is \
                             proved by `rfl`, so the kernel computes the differing \
                             and comparable sites rather than being handed them. \
                             This closes the gap an adversarial audit found in the \
                             integer-only form.",
            "what_is_still_trusted": "that the embedded lists are the dataset's \
                                      sequences -- checkable by diffing against \
                                      the alignment, and far narrower than \
                                      trusting the counting code.",
        },
        "kernel_verdict": {
            "emitted_file": good_path.display().to_string(),
            "genuine_certificates_accepted": good_accepted,
            "genuine_first_error": good_err,
            "corrupted_file": bad_path.display().to_string(),
            "corrupted_rejected_by_rust_validate": bad_is_invalid,
            "corrupted_rejected_by_kernel": !bad_accepted,
            "corrupted_first_error": bad_err,
            "lean_module": "ANSE.GeneticaReticulata",
            "inference_rule": "count_certificate_refutes_ultrametric",
        },
        "verdict": {
            "acceptance_control_passes": accept_passes,
            "kernel_counted_the_data": seq_accepted,
            "rejection_control_passes": reject_passes,
            "discriminates_recombinant": discriminates,
            "criterion": "the kernel must ACCEPT genuine certificates AND REJECT \
                          a corrupted one; and the planted recombinant must be a \
                          certificate midpoint with margin > 2x the clonal floor",
        },
        "what_this_does_not_show": "That a certificate implies recombination. It \
            implies only that the data is not an ultrametric, which finite \
            sampling alone achieves. Calibration against the clonal baseline is \
            what carries the biological reading, and it remains simulation-based.",
    });
    println!("{}", serde_json::to_string_pretty(&report).unwrap());

    if !accept_passes {
        eprintln!(
            "FAIL acceptance control: the kernel rejected a genuine certificate. \
             First error: {good_err}"
        );
        return Ok(false);
    }
    if !reject_passes {
        eprintln!(
            "FAIL rejection control: the kernel ACCEPTED a corrupted certificate. \
             The verifier is not checking anything and no result above is usable."
        );
        return Ok(false);
    }
    if !discriminates {
        eprintln!(
            "FAIL discrimination: margin ratio {ratio:.2} or the recombinant was \
             not a midpoint. The architecture works but has no power here."
        );
        return Ok(false);
    }
    Ok(true)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.first().map(String::as_str) {
        Some("controls") => run_controls(&args),
        Some("poc") => run_poc(&args),
        Some("fasta") => run_fasta(&args),
        _ => {
            eprintln!(
                "usage:\n  reticulation controls [--length N] [--generations G] \
                 [--mu M] [--seed S]\n  reticulation fasta <FILE> [--model p|jc] \
                 [--top K] [--scale X]"
            );
            return ExitCode::from(2);
        }
    };
    match result {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::FAILURE,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

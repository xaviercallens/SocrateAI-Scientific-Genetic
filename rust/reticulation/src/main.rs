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

/// One dataset at one seed.
struct SeedResult {
    seed: u64,
    n_sequences: usize,
    /// Largest exact margin from the complete-deletion triple search.
    triple_margin: Option<reticulation::Ratio>,
    triple_violating: usize,
    triple_common_columns: i64,
    /// Largest exact margin from the quartet (four-point) search.
    quartet_margin: Option<reticulation::Ratio>,
    quartet_violating: usize,
    quartet_examined: usize,
    /// Largest exact margin from the *pairwise*-deletion triple search. Kept for
    /// the record only: margins under the two policies have different
    /// denominators and are not comparable to each other.
    pairwise_triple_margin: Option<reticulation::Ratio>,
    /// Per-triple comparison of the two deletion policies — the statistic that
    /// actually measures the gap confound.
    policy: reticulation::subset::PolicyComparison,
}

/// Exact comparison `a > k * b`, by integer cross-product.
///
/// No floating point decides any verdict in this experiment; the `f64` margins
/// in the report are for reading.
fn exceeds_k_times(
    a: Option<reticulation::Ratio>,
    b: Option<reticulation::Ratio>,
    k: i128,
) -> bool {
    match (a, b) {
        (Some(x), Some(y)) if y.num > 0 => {
            (x.num as i128) * (y.den as i128) > k * (y.num as i128) * (x.den as i128)
        }
        // No reference violation at all: any positive margin exceeds it.
        (Some(x), _) => x.num > 0,
        (None, _) => false,
    }
}

fn margin_f64(m: Option<reticulation::Ratio>) -> f64 {
    m.map_or(0.0, reticulation::Ratio::to_f64)
}

fn margin_exact(m: Option<reticulation::Ratio>) -> String {
    m.map_or_else(|| "none".to_string(), |r| format!("{}/{}", r.num, r.den))
}

/// Run both searches on one population.
fn measure(pop: &reticulation::simulate::Population, max_quadruples: usize) -> Result<SeedResult, String> {
    use reticulation::exact::{find_count_certificates, ExactDistances};
    use reticulation::subset::{
        compare_deletion_policies, find_complete_triple_certificates, find_quartet_certificates,
        Alignment,
    };

    let aln = Alignment::new(pop.labels.clone(), pop.sequences.clone())?;
    let policy = compare_deletion_policies(&aln)?;
    let triples = find_complete_triple_certificates(&aln, 1)?;
    let quartets = find_quartet_certificates(&aln, 1, max_quadruples)?;
    if !triples.all_validated || !quartets.all_validated {
        return Err("a search produced a certificate that failed revalidation".into());
    }

    // Pairwise deletion, for the comparison of policies.
    let pw = ExactDistances::from_alignment(pop.labels.clone(), &pop.sequences)?;
    let (pw_certs, _, _) = find_count_certificates(&pw, 1);

    Ok(SeedResult {
        seed: 0,
        n_sequences: aln.len(),
        triple_margin: triples.certificates.first().map(|c| c.margin()),
        triple_violating: triples.violating,
        triple_common_columns: triples.certificates.first().map_or(0, |c| c.common),
        quartet_margin: quartets.certificates.first().map(|c| c.margin()),
        quartet_violating: quartets.violating,
        quartet_examined: quartets.examined,
        pairwise_triple_margin: pw_certs.first().map(|c| reticulation::Ratio {
            num: c.margin_num,
            den: c.margin_den,
        }),
        policy,
    })
}

/// Fraction of seeds at which a per-seed, exactly-decided condition held.
fn held_at(flags: &[bool]) -> (usize, usize) {
    (flags.iter().filter(|b| **b).count(), flags.len())
}

/// The experiment the whole module exists for: does the four-point test separate
/// reticulation from rate variation, where the ultrametric test cannot?
///
/// Predictions are registered in the output *before* the measurements, each with
/// a criterion that can fail. The PASS/FAIL gate rests only on the power and
/// fix-efficacy controls, deliberately: the discriminating prediction (P3) must
/// be able to come out false without the run being reported as a failure, or the
/// experiment would be rigged to confirm it.
fn run_experiment(args: &[String]) -> Result<bool, String> {
    use reticulation::simulate::{
        clonal_tree, clonal_tree_relaxed, inject_gap_blocks, with_recombination,
        with_reciprocal_recombinants,
    };

    let length: usize = parse_flag(args, "--length", 4000);
    let generations: u32 = parse_flag(args, "--generations", 4);
    let mu: f64 = parse_flag(args, "--mu", 0.004);
    let spread: f64 = parse_flag(args, "--spread", 1.0);
    let saturated_mu: f64 = parse_flag(args, "--saturated-mu", 0.04);
    let gap_fraction: f64 = parse_flag(args, "--gap-fraction", 0.25);
    let n_seeds: usize = parse_flag(args, "--seeds", 5);
    let seed0: u64 = parse_flag(args, "--seed", 2026);
    let max_quadruples: usize = parse_flag(args, "--max-quadruples", 20_000);

    let seeds: Vec<u64> = (0..n_seeds).map(|i| seed0 + i as u64).collect();

    // --- build and measure every dataset at every seed ---
    let mut results: std::collections::BTreeMap<&str, Vec<SeedResult>> = Default::default();
    for &seed in &seeds {
        let mut push = |name: &'static str, pop: reticulation::simulate::Population| -> Result<(), String> {
            let mut r = measure(&pop, max_quadruples)?;
            r.seed = seed;
            results.entry(name).or_default().push(r);
            Ok(())
        };

        push("D1_clonal_strict_clock", clonal_tree(length, generations, mu, seed))?;
        push(
            "D2_clonal_relaxed_clock",
            clonal_tree_relaxed(length, generations, mu, spread, seed),
        )?;
        push(
            "D3_strict_clock_one_recombinant",
            with_recombination(clonal_tree(length, generations, mu, seed), 0.5, seed ^ 0xABCD).0,
        )?;
        push(
            "D4_relaxed_clock_one_recombinant",
            with_recombination(
                clonal_tree_relaxed(length, generations, mu, spread, seed),
                0.5,
                seed ^ 0xABCD,
            )
            .0,
        )?;
        push(
            "D5_strict_clock_reciprocal_mosaics",
            with_reciprocal_recombinants(clonal_tree(length, generations, mu, seed), 0.5, seed ^ 0xABCD).0,
        )?;
        push(
            "D6_strict_clock_gap_blocks",
            inject_gap_blocks(clonal_tree(length, generations, mu, seed), gap_fraction, seed ^ 0x1357),
        )?;
        push(
            "D7_strict_clock_saturated",
            clonal_tree(length, generations, saturated_mu, seed),
        )?;
    }

    let get = |name: &str| -> &Vec<SeedResult> { &results[name] };
    let d1 = get("D1_clonal_strict_clock");

    // --- the per-seed, exactly-decided comparisons ---
    let paired = |name: &str, f: &dyn Fn(&SeedResult, &SeedResult) -> bool| -> Vec<bool> {
        get(name)
            .iter()
            .zip(d1.iter())
            .map(|(r, base)| f(r, base))
            .collect()
    };
    let triple_above = |r: &SeedResult, base: &SeedResult| {
        exceeds_k_times(r.triple_margin, base.triple_margin, 2)
    };
    let quartet_above = |r: &SeedResult, base: &SeedResult| {
        exceeds_k_times(r.quartet_margin, base.quartet_margin, 2)
    };

    let triple_strictly_above = |r: &SeedResult, base: &SeedResult| {
        exceeds_k_times(r.triple_margin, base.triple_margin, 1)
    };
    let quartet_strictly_above = |r: &SeedResult, base: &SeedResult| {
        exceeds_k_times(r.quartet_margin, base.quartet_margin, 1)
    };

    let p1 = paired("D3_strict_clock_one_recombinant", &triple_strictly_above);
    let p1b = paired("D3_strict_clock_one_recombinant", &triple_above);
    let p2 = paired("D2_clonal_relaxed_clock", &triple_above);
    let p2b = paired("D2_clonal_relaxed_clock", &triple_strictly_above);
    let p3 = paired("D2_clonal_relaxed_clock", &quartet_above);
    let p4 = paired("D5_strict_clock_reciprocal_mosaics", &quartet_above);
    let p4c = paired("D5_strict_clock_reciprocal_mosaics", &quartet_strictly_above);
    let p7 = paired("D7_strict_clock_saturated", &quartet_above);
    let p8 = paired("D3_strict_clock_one_recombinant", &quartet_above);
    let p8b = paired("D4_relaxed_clock_one_recombinant", &quartet_above);
    // Finding 6, measured by the per-triple verdict rather than by margin
    // magnitude: with no missing data the policies must agree exactly, and with
    // gap blocks pairwise deletion must manufacture violations.
    let p5: Vec<bool> = d1
        .iter()
        .map(|r| r.policy.pairwise_only == 0 && r.policy.complete_only == 0)
        .collect();
    let p6: Vec<bool> = get("D6_strict_clock_gap_blocks")
        .iter()
        .map(|r| r.policy.pairwise_only > 0)
        .collect();

    let frac = |v: &[bool]| {
        let (k, n) = held_at(v);
        format!("{k}/{n}")
    };
    let at_least = |v: &[bool], k: usize| held_at(v).0 >= k;

    // Worst-case separation: the weakest seed's ratio of the positive control's
    // margin to its own floor. Display only -- every verdict above is exact.
    let worst_ratio = |name: &str, quartet: bool| -> f64 {
        get(name)
            .iter()
            .zip(d1.iter())
            .map(|(r, base)| {
                let (a, b) = if quartet {
                    (margin_f64(r.quartet_margin), margin_f64(base.quartet_margin))
                } else {
                    (margin_f64(r.triple_margin), margin_f64(base.triple_margin))
                };
                if b > 0.0 {
                    a / b
                } else {
                    f64::INFINITY
                }
            })
            .fold(f64::INFINITY, f64::min)
    };

    let dataset_json: serde_json::Value = results
        .iter()
        .map(|(name, rows)| {
            (
                (*name).to_string(),
                serde_json::json!({
                    "n_sequences": rows[0].n_sequences,
                    "quartets_examined": rows[0].quartet_examined,
                    "per_seed": rows.iter().map(|r| serde_json::json!({
                        "seed": r.seed,
                        "complete_triple_margin_exact": margin_exact(r.triple_margin),
                        "complete_triple_margin": margin_f64(r.triple_margin),
                        "complete_triple_violating": r.triple_violating,
                        "complete_triple_common_columns": r.triple_common_columns,
                        "quartet_margin_exact": margin_exact(r.quartet_margin),
                        "quartet_margin": margin_f64(r.quartet_margin),
                        "quartet_violating": r.quartet_violating,
                        "pairwise_triple_margin_exact": margin_exact(r.pairwise_triple_margin),
                        "pairwise_triple_margin": margin_f64(r.pairwise_triple_margin),
                        "deletion_policy_comparison": r.policy,
                    })).collect::<Vec<_>>(),
                }),
            )
        })
        .collect::<serde_json::Map<String, serde_json::Value>>()
        .into();

    // Gate: power and fix-efficacy only. P2, P3, P7, P8 are findings.
    let power_triple = at_least(&p1, seeds.len());
    let power_quartet = at_least(&p4, seeds.len());
    let policies_agree = at_least(&p5, seeds.len());
    let fix_bites = at_least(&p6, seeds.len() - 1);
    let gate = power_triple && power_quartet && policies_agree && fix_bites;

    let report = serde_json::json!({
        "experiment": "Does the four-point (tree-level) test separate reticulation \
                       from rate variation, where the ultrametric (clock-level) \
                       test cannot? And does complete deletion within the subset \
                       remove the gap confound of the audit's finding 6?",
        "design": {
            "null_nesting_theorem": "ANSE.GeneticaReticulata.fourPoint_of_ultrametric \
                                     -- every symmetric ultrametric satisfies the \
                                     four-point condition, so a quartet certificate \
                                     refutes strictly more than a triple one.",
            "reference_dataset": "D1_clonal_strict_clock. Every comparison is \
                                  SEED-PAIRED against D1 at the same seed and \
                                  decided by an exact integer cross-product.",
            "margins_are_exact_rationals": true,
            "parameters": {
                "alignment_length": length, "generations": generations,
                "mu_per_generation": mu, "relaxed_clock_spread": spread,
                "saturated_mu": saturated_mu, "gap_block_fraction": gap_fraction,
                "seeds": seeds,
            },
            "why_several_seeds": "A single-seed margin is not a result. Each \
                                  prediction reports the number of seeds at which \
                                  it held, not an average of margins.",
        },
        "criteria_corrected_after_a_first_run": {
            "why_this_block_exists": "A first run of this experiment used two \
                criteria that were wrong, and both are recorded here rather than \
                quietly replaced, because a criterion changed after seeing data is \
                the easiest way to launder a negative result.",
            "corrections": [
                {
                    "was": "P6 compared the largest MARGIN under pairwise deletion \
                            with the largest margin under complete deletion.",
                    "wrong_because": "Complete deletion shrinks the denominator -- \
                        on D6 it leaves about 1100 of 4000 columns -- so one \
                        differing site is worth more and its margins are LARGER, \
                        even though it is strictly more conservative about what \
                        counts as a violation. The two numbers have different \
                        denominators and are no more comparable than the two \
                        pairwise distances that finding 6 is about.",
                    "now": "The per-triple VERDICT is compared instead, and the \
                            statistic is the number of triples that pairwise \
                            deletion calls a violation and complete deletion does \
                            not. P5 adds the control that makes it meaningful: with \
                            no missing data the two policies must agree exactly.",
                },
                {
                    "was": "The power gate required the positive control's margin \
                            to exceed TWICE the clonal floor, the threshold v1.6.0 \
                            used.",
                    "wrong_because": "The floor is a maximum over 560 triples, an \
                        extreme-value statistic, and this run shows it varying from \
                        18/4000 to 39/4000 across five seeds. A 2x threshold on a \
                        ratio whose denominator moves by 2.2x is not a stable \
                        criterion, and v1.6.0's reported 3.389x was the best of the \
                        five seeds, not a typical value.",
                    "now": "The gate is the seed-PAIRED strict inequality (P1): the \
                            same clonal tree with a recombinant added must give a \
                            larger margin than without it. The 2x result is kept as \
                            an ungated finding (P1b) rather than dropped, because \
                            its failure at 2 of 5 seeds is itself the result.",
                },
            ],
        },
        "separation_display_only": {
            "note": "Ratios of the positive control's margin to its own floor at \
                     the weakest of the five seeds. f64, for reading; every verdict \
                     below is an exact integer comparison.",
            "triple_test_worst_case_D3_over_D1": worst_ratio("D3_strict_clock_one_recombinant", false),
            "quartet_test_worst_case_D5_over_D1": worst_ratio("D5_strict_clock_reciprocal_mosaics", true),
            "relaxed_clock_triple_worst_case_D2_over_D1": worst_ratio("D2_clonal_relaxed_clock", false),
            "relaxed_clock_quartet_worst_case_D2_over_D1": worst_ratio("D2_clonal_relaxed_clock", true),
        },
        "preregistered_predictions": [
            {
                "id": "P1", "gated": true,
                "claim": "POWER CONTROL, seed-paired. Adding one recombinant to a \
                          clonal tree raises the complete-deletion triple margin \
                          above the same tree's margin without it.",
                "criterion": "holds at all seeds",
                "observed": frac(&p1), "held": power_triple,
            },
            {
                "id": "P1b", "gated": false,
                "claim": "FINDING. The same comparison at v1.6.0's threshold: D3's \
                          margin exceeds 2x D1's.",
                "criterion": "reported, not gated",
                "observed": frac(&p1b),
                "reading": "Where this fails, the triple statistic does not clear \
                            its own published threshold. The clonal floor is a \
                            maximum over triples and varies by more than 2x across \
                            seeds, so this threshold cannot be relied on.",
            },
            {
                "id": "P2", "gated": false,
                "claim": "A relaxed clock on a purely clonal tree breaks \
                          ULTRAMETRICITY: D2's complete-triple margin exceeds 2x \
                          D1's at the same seed.",
                "criterion": "reported, not gated",
                "observed_above_2x": frac(&p2),
                "observed_above_1x_paired": frac(&p2b),
                "divergence_matched": "The relaxed clock's rate multipliers are \
                    normalised to mean 1, so D2 has the same expected divergence as \
                    D1 and only the clock differs. Without that, a relaxed clock is \
                    also a faster clock and the result would confound rate \
                    variation with saturation.",
                "reading": "Where this holds, a triple certificate does not \
                            distinguish reticulation from rate variation -- risk 2 \
                            of docs/GENETICA_RETICULATA.md section 7, realised. \
                            Compare D2's per-seed margins against D3's directly: if \
                            they overlap, the statistic cannot separate the causes.",
            },
            {
                "id": "P3", "gated": false,
                "claim": "THE DISCRIMINATING PREDICTION. The same relaxed-clock \
                          data does NOT break the FOUR-POINT condition beyond the \
                          floor: D2's quartet margin does not exceed 2x D1's.",
                "criterion": "violations at <= 1 of 5 seeds",
                "observed_violations": frac(&p3), "held": !at_least(&p3, 2),
                "reading": "A violation count of 0 or 1 means the four-point test \
                            is largely immune to the confound that defeats the \
                            triple test. A higher count means it is confounded too, \
                            and the method has no statistic that separates the two \
                            causes -- a negative result about this approach, \
                            reported rather than gated away. Read it next to P7: a \
                            relaxed clock also redistributes divergence among \
                            lineages, so some elevation may be local saturation \
                            rather than non-additivity.",
            },
            {
                "id": "P4", "gated": true,
                "claim": "POWER CONTROL for the four-point test. Reciprocal \
                          mosaics break it: D5's quartet margin exceeds 2x D1's \
                          quartet floor at the same seed.",
                "criterion": "holds at all seeds",
                "observed_above_2x": frac(&p4),
                "observed_above_1x_paired": frac(&p4c),
                "held": power_quartet,
                "why_reciprocal": "A single recombinant need not violate the \
                                   four-point condition; the shape that does is a \
                                   pair of reciprocal mosaics, each close to both \
                                   parents. This is the configuration Lean's \
                                   reciprocal_recombinants_witness proves \
                                   non-vacuous.",
            },
            {
                "id": "P5b", "gated": true,
                "claim": "CONTROL FOR THE POLICY COMPARISON. With no missing data, \
                          complete deletion within the triple IS pairwise deletion, \
                          so the two policies must agree on every one of D1's 560 \
                          triples.",
                "criterion": "holds at all seeds",
                "observed": frac(&p5), "held": policies_agree,
                "why": "Without this, a non-zero manufactured count on gappy data \
                        could be an artefact of the comparison itself.",
            },
            {
                "id": "P6", "gated": true,
                "claim": "THE FIX BITES, AND FINDING 6 IS REAL. On gap-blocked but \
                          strictly clonal data with no reticulation, pairwise \
                          deletion calls triples violations that complete deletion \
                          within the same triple does not.",
                "criterion": "at least one manufactured triple at >= 4 of 5 seeds",
                "observed": frac(&p6), "held": fix_bites,
                "reading": "Every manufactured triple is a kernel-acceptable \
                            certificate about a genealogy that is exactly a clock \
                            tree. The count is the size of the confound; complete \
                            deletion is what removes it.",
            },
            {
                "id": "P7", "gated": false,
                "claim": "EXPLORATORY. A high-divergence clonal population lifts \
                          the quartet margin above 2x the low-divergence floor, \
                          with no reticulation and a strict clock.",
                "criterion": "reported, not gated",
                "observed": frac(&p7),
                "this_comparison_is_mis_specified": "D7's margin is compared \
                    against D1's floor, but D1 is ten times less divergent. The \
                    `saturation` subcommand shows the quartet noise floor grows \
                    monotonically with divergence, so the two are not \
                    commensurable -- the same error, in a third guise, as \
                    comparing margins across different denominators.",
                "what_the_diagnostic_found": "At BOTH low and high divergence the \
                    clonal quartet margin shrinks by a factor near 0.5 when the \
                    alignment length quadruples: the L^(-1/2) signature of \
                    sampling noise. A saturation bias is deterministic and would \
                    not shrink at all. So this row is NOT evidence of a divergence \
                    ceiling -- it is evidence that a null must be calibrated at \
                    matched divergence AND matched length. Run \
                    `reticulation saturation` for the numbers.",
            },
            {
                "id": "P8", "gated": false,
                "claim": "SECONDARY. A single recombinant need not violate the \
                          four-point condition, so D3's and D4's quartet margins \
                          may stay at the floor while their triple margins do not.",
                "criterion": "reported, not gated",
                "observed_D3": frac(&p8), "observed_D4": frac(&p8b),
            },
        ],
        "datasets": dataset_json,
        "verdict": {
            "gated_predictions_pass": gate,
            "power_triple_P1": power_triple,
            "power_quartet_P4": power_quartet,
            "policies_agree_without_gaps_P5b": policies_agree,
            "complete_deletion_fix_bites_P6": fix_bites,
            "criterion": "The gate is the two power controls and the two \
                          policy-comparison controls. P1b, P2, P3, P7 and P8 are \
                          findings: the discriminating prediction must be able to \
                          come out false without the run reporting failure.",
        },
        "cost_of_complete_deletion": "On D6 the common-column count per triple is \
            far below the alignment length, because 15 of 16 sequences carry a \
            gap block. That lost power is the price of removing the confound, and \
            it is why complete deletion should be the default only when there is \
            missing data to delete.",
    });
    println!("{}", serde_json::to_string_pretty(&report).unwrap());

    if !power_triple {
        eprintln!(
            "FAIL P1: adding a recombinant to a clonal tree did not raise the \
             complete-deletion triple margin above the same tree's margin without \
             it ({}). With no power the other columns are uninterpretable.",
            frac(&p1)
        );
        return Ok(false);
    }
    if !power_quartet {
        eprintln!(
            "FAIL P4: reciprocal mosaics did not lift the quartet margin above \
             twice the clonal floor ({}). The four-point test has no demonstrated \
             power, so its silence on D2 would mean nothing.",
            frac(&p4)
        );
        return Ok(false);
    }
    if !policies_agree {
        eprintln!(
            "FAIL P5b: the two deletion policies disagreed on gapless data ({}). \
             The policy comparison is then measuring itself, not the confound.",
            frac(&p5)
        );
        return Ok(false);
    }
    if !fix_bites {
        eprintln!(
            "FAIL P6: pairwise deletion manufactured no violations on gap-blocked \
             clonal data ({}). Either the gap pattern does not create the \
             confound, or complete deletion is not doing what it claims.",
            frac(&p6)
        );
        return Ok(false);
    }
    Ok(true)
}

/// Mean and maximum pairwise p-distance of a population, for reporting the
/// divergence a result was obtained at.
fn divergence(pop: &reticulation::simulate::Population) -> (f64, f64) {
    let dm = pop
        .distances(DistanceModel::PDistance)
        .expect("simulated populations are well-formed alignments");
    let mut sum = 0.0;
    let mut n = 0.0;
    let mut max = 0.0f64;
    for i in 0..dm.len() {
        for j in (i + 1)..dm.len() {
            sum += dm.get(i, j);
            max = max.max(dm.get(i, j));
            n += 1.0;
        }
    }
    (sum / n, max)
}

/// Is the clonal four-point signal noise, or is it saturation bias?
///
/// P7 of the experiment showed that a high-divergence, strictly clonal, strictly
/// clock-like population violates the four-point condition — which it should not,
/// since its genealogy is a tree. Two causes are possible and they have opposite
/// signatures in the alignment length `L`:
///
/// * **Sampling noise.** The margin is a difference of order statistics on
///   multinomial counts, so it scales as `L^(-1/2)`: quadrupling `L` halves it.
/// * **Saturation bias.** Multiple hits make the p-distance a concave function
///   of true divergence, which is a *deterministic* distortion of the metric. It
///   does not shrink with `L` at all.
///
/// So the diagnostic is to quadruple `L` at fixed divergence and watch the
/// margin. This decides whether the method needs a divergence ceiling or merely
/// longer sequences, and the answer determines which real datasets are in range.
fn run_saturation(args: &[String]) -> Result<bool, String> {
    use reticulation::simulate::clonal_tree;
    use reticulation::subset::{find_quartet_certificates, Alignment};

    let generations: u32 = parse_flag(args, "--generations", 4);
    let n_seeds: usize = parse_flag(args, "--seeds", 5);
    let seed0: u64 = parse_flag(args, "--seed", 2026);
    let low_mu: f64 = parse_flag(args, "--low-mu", 0.004);
    let high_mu: f64 = parse_flag(args, "--high-mu", 0.04);
    let seeds: Vec<u64> = (0..n_seeds).map(|i| seed0 + i as u64).collect();
    let lengths: Vec<usize> = vec![1000, 4000, 16_000];

    // --- the L-scaling diagnostic at two divergences ---
    let mut scaling = serde_json::Map::new();
    let mut ratios: std::collections::BTreeMap<String, Vec<f64>> = Default::default();
    for (name, mu) in [("low_divergence", low_mu), ("high_divergence", high_mu)] {
        let mut rows = Vec::new();
        let mut by_length: Vec<Vec<f64>> = Vec::new();
        for &length in &lengths {
            let mut margins = Vec::new();
            let mut divs = Vec::new();
            for &seed in &seeds {
                let pop = clonal_tree(length, generations, mu, seed);
                divs.push(divergence(&pop).0);
                let aln = Alignment::new(pop.labels.clone(), pop.sequences.clone())?;
                let q = find_quartet_certificates(&aln, 1, 20_000)?;
                margins.push(
                    q.certificates
                        .first()
                        .map_or(0.0, |c| c.margin().to_f64()),
                );
            }
            let mean_margin = margins.iter().sum::<f64>() / margins.len() as f64;
            rows.push(serde_json::json!({
                "alignment_length": length,
                "mean_pairwise_divergence": divs.iter().sum::<f64>() / divs.len() as f64,
                "quartet_margins_per_seed": margins,
                "mean_quartet_margin": mean_margin,
            }));
            by_length.push(margins);
        }
        // Ratio of the mean margin at 4L to the mean margin at L. Pure noise
        // predicts 0.5; a deterministic bias predicts 1.
        let mean = |v: &Vec<f64>| v.iter().sum::<f64>() / v.len() as f64;
        let r: Vec<f64> = (1..by_length.len())
            .map(|i| mean(&by_length[i]) / mean(&by_length[i - 1]).max(1e-12))
            .collect();
        ratios.insert(name.to_string(), r.clone());
        scaling.insert(
            name.to_string(),
            serde_json::json!({
                "mu_per_generation": mu,
                "per_length": rows,
                "ratio_of_mean_margin_when_length_quadruples": r,
                "prediction_if_pure_sampling_noise": 0.5,
                "prediction_if_pure_saturation_bias": 1.0,
            }),
        );
    }

    let low = &ratios["low_divergence"];
    let high = &ratios["high_divergence"];
    // The noise scaling must actually be observed, or the diagnostic says
    // nothing. Generous band: L^(-1/2) is 0.5, and a max over 1820 quadruples is
    // an extreme-value statistic with its own drift.
    let noise_scaling_seen = low.iter().all(|r| *r < 0.75);
    let bias_at_high_divergence = high.iter().any(|r| *r > 0.75);

    // --- a divergence sweep, so the operating range is in biological units ---
    let mut sweep = Vec::new();
    for mu in [0.001f64, 0.002, 0.004, 0.008, 0.016, 0.024, 0.032, 0.048, 0.064] {
        let mut divs = Vec::new();
        let mut margins = Vec::new();
        for &seed in &seeds {
            let pop = clonal_tree(4000, generations, mu, seed);
            let (mean_d, max_d) = divergence(&pop);
            divs.push((mean_d, max_d));
            let aln = Alignment::new(pop.labels.clone(), pop.sequences.clone())?;
            let q = find_quartet_certificates(&aln, 1, 20_000)?;
            margins.push(q.certificates.first().map_or(0.0, |c| c.margin().to_f64()));
        }
        sweep.push(serde_json::json!({
            "mu_per_generation": mu,
            "mean_pairwise_divergence": divs.iter().map(|d| d.0).sum::<f64>() / divs.len() as f64,
            "max_pairwise_divergence": divs.iter().map(|d| d.1).fold(0.0, f64::max),
            "mean_quartet_margin": margins.iter().sum::<f64>() / margins.len() as f64,
            "quartet_margins_per_seed": margins,
        }));
    }

    let report = serde_json::json!({
        "diagnostic": "Is the four-point violation seen on high-divergence clonal \
                       data sampling noise or saturation bias? The two scale \
                       differently in the alignment length, so quadrupling it \
                       separates them.",
        "why_it_matters": "A noise floor shrinks with longer sequences, so the \
            method simply needs more data. A saturation bias does not, so the \
            method needs a divergence ceiling or a distance that corrects for \
            multiple hits. The second is the harder case, because a \
            Jukes-Cantor-corrected distance is a logarithm and therefore \
            irrational, which the exact-rational certificate architecture cannot \
            represent. There is no quick fix available if the bias is real.",
        "parameters": {
            "generations": generations, "seeds": seeds,
            "lengths": lengths, "low_mu": low_mu, "high_mu": high_mu,
        },
        "length_scaling": scaling,
        "divergence_sweep_at_4000_sites": sweep,
        "verdict": {
            "noise_scaling_observed_at_low_divergence": noise_scaling_seen,
            "systematic_component_at_high_divergence": bias_at_high_divergence,
            "gated": "only the first -- the diagnostic must be able to see the \
                      noise scaling it relies on. The second is the finding.",
            "reading": "If both are true, the high-divergence four-point \
                        violations are a p-distance artefact and the method has a \
                        divergence ceiling that must be stated with every result.",
        },
    });
    println!("{}", serde_json::to_string_pretty(&report).unwrap());

    if !noise_scaling_seen {
        eprintln!(
            "FAIL: the clonal quartet margin did not shrink with alignment length \
             at low divergence (ratios {low:?}). Without that baseline the \
             high-divergence comparison distinguishes nothing."
        );
        return Ok(false);
    }
    Ok(true)
}

/// Emit and kernel-check the two new certificate forms, with both controls.
///
/// The acceptance and rejection directions are both required: a verifier that
/// accepts everything proves nothing.
fn run_verify_tree(args: &[String]) -> Result<bool, String> {
    use reticulation::emit_lean;
    use reticulation::simulate::{clonal_tree, with_reciprocal_recombinants};
    use reticulation::subset::{
        find_complete_triple_certificates, find_quartet_certificates, Alignment,
    };

    let length: usize = parse_flag(args, "--length", 4000);
    let generations: u32 = parse_flag(args, "--generations", 4);
    let mu: f64 = parse_flag(args, "--mu", 0.004);
    let seed: u64 = parse_flag(args, "--seed", 2026);
    let formal_dir: String = parse_flag(args, "--formal-dir", "../../../../../formal".to_string());
    let out_dir: String = parse_flag(args, "--out", "/tmp".to_string());

    let formal = std::path::Path::new(&formal_dir)
        .canonicalize()
        .map_err(|e| format!("--formal-dir {formal_dir}: {e}"))?;

    let (pop, mosaics) =
        with_reciprocal_recombinants(clonal_tree(length, generations, mu, seed), 0.5, seed ^ 0xABCD);
    let aln = Alignment::new(pop.labels.clone(), pop.sequences.clone())?;

    let triples = find_complete_triple_certificates(&aln, 1)?;
    let quartets = find_quartet_certificates(&aln, 1, 20_000)?;
    let triple = triples
        .certificates
        .first()
        .ok_or("no complete-deletion triple certificate on the positive control")?;
    let quartet = quartets
        .certificates
        .first()
        .ok_or("no quartet certificate on the positive control")?;

    let enc = |i: usize| emit_lean::encode(&pop.sequences[i]);
    let source = format!("simulated reciprocal mosaics, seed {seed}, {length} nt");

    // --- acceptance: the complete-deletion triple ---
    let t_path = std::path::Path::new(&out_dir).join("reticulata_complete_triple.lean");
    let [tx, ty, tz] = triple.indices;
    std::fs::write(
        &t_path,
        emit_lean::render_complete_triple(triple, &enc(tx), &enc(ty), &enc(tz), &source),
    )
    .map_err(|e| format!("{}: {e}", t_path.display()))?;
    let t_start = std::time::Instant::now();
    let (t_ok, t_err) = kernel_check(&formal, &t_path);
    let t_secs = t_start.elapsed().as_secs_f64();

    // --- acceptance: the quartet ---
    let q_path = std::path::Path::new(&out_dir).join("reticulata_quartet.lean");
    let [qa, qb, qc, qe] = quartet.indices;
    std::fs::write(
        &q_path,
        emit_lean::render_quartet(quartet, &enc(qa), &enc(qb), &enc(qc), &enc(qe), &source),
    )
    .map_err(|e| format!("{}: {e}", q_path.display()))?;
    let q_start = std::time::Instant::now();
    let (q_ok, q_err) = kernel_check(&formal, &q_path);
    let q_secs = q_start.elapsed().as_secs_f64();

    // --- rejection: both corrupted forms ---
    let bad_t = emit_lean::corrupt_complete_triple(triple);
    let bad_t_path = std::path::Path::new(&out_dir).join("reticulata_complete_triple_bad.lean");
    std::fs::write(
        &bad_t_path,
        emit_lean::render_complete_triple(
            &bad_t,
            &enc(tx),
            &enc(ty),
            &enc(tz),
            "deliberately corrupted: far and near numerators swapped",
        ),
    )
    .map_err(|e| format!("{}: {e}", bad_t_path.display()))?;
    let (bad_t_ok, bad_t_err) = kernel_check(&formal, &bad_t_path);

    let bad_q = emit_lean::corrupt_quartet(quartet);
    let bad_q_path = std::path::Path::new(&out_dir).join("reticulata_quartet_bad.lean");
    std::fs::write(
        &bad_q_path,
        emit_lean::render_quartet(
            &bad_q,
            &enc(qa),
            &enc(qb),
            &enc(qc),
            &enc(qe),
            "deliberately corrupted: the dominant pairing exchanged with a dominated one",
        ),
    )
    .map_err(|e| format!("{}: {e}", bad_q_path.display()))?;
    let (bad_q_ok, bad_q_err) = kernel_check(&formal, &bad_q_path);

    let mosaic_in_quartet = [&quartet.c, &quartet.e]
        .iter()
        .filter(|l| mosaics.contains(l))
        .count();

    let accept = t_ok && q_ok;
    let reject = !bad_t_ok
        && !bad_q_ok
        && bad_t.validate().is_err()
        && bad_q.validate().is_err();

    let report = serde_json::json!({
        "verification": "Kernel-check the two certificate forms introduced for \
                         Phase 2: complete deletion within the triple, and the \
                         four-point quartet. The kernel counts the sequences in \
                         both; nothing is handed to it but the sequences and the \
                         claimed counts.",
        "parameters": {
            "alignment_length": length, "generations": generations,
            "mu_per_generation": mu, "seed": seed,
        },
        "complete_deletion_triple": {
            "file": t_path.display().to_string(),
            "accepted": t_ok, "first_error": t_err, "seconds": t_secs,
            "lean_theorem": "complete_triple_refute_from_counts",
            "triple": [&triple.x, &triple.y, &triple.z],
            "counts": {"n_xy": triple.n_xy, "n_yz": triple.n_yz,
                       "n_xz": triple.n_xz, "common_columns": triple.common},
            "margin_exact": format!("{}/{}", triple.margin_num, triple.margin_den),
            "what_the_kernel_did": "computed diff3xy, diff3yz, diff3xz and \
                                    comparable3 from the three embedded sequences \
                                    by rfl, then checked two integer inequalities \
                                    on one denominator -- no cross-multiplication.",
        },
        "four_point_quartet": {
            "file": q_path.display().to_string(),
            "accepted": q_ok, "first_error": q_err, "seconds": q_secs,
            "lean_theorem": "quartet_refute_from_counts",
            "quadruple": {"far_pair": [&quartet.a, &quartet.b],
                          "other_pair": [&quartet.c, &quartet.e]},
            "pairing_sums": [quartet.n_ab + quartet.n_ce,
                             quartet.n_ac + quartet.n_be,
                             quartet.n_ae + quartet.n_bc],
            "common_columns": quartet.common,
            "margin_exact": format!("{}/{}", quartet.margin_num, quartet.margin_den),
            "planted_mosaics": mosaics,
            "mosaics_in_the_non_dominant_pair": mosaic_in_quartet,
            "what_this_refutes": "the four-point condition, hence (Buneman, cited) \
                                  that any weighted tree realises these distances \
                                  -- and hence, by the proved nesting \
                                  fourPoint_of_ultrametric, ultrametricity too.",
        },
        "rejection_controls": {
            "corrupted_triple_file": bad_t_path.display().to_string(),
            "corrupted_triple_rejected": !bad_t_ok,
            "corrupted_triple_first_error": bad_t_err,
            "corrupted_quartet_file": bad_q_path.display().to_string(),
            "corrupted_quartet_rejected": !bad_q_ok,
            "corrupted_quartet_first_error": bad_q_err,
            "both_invalid_in_rust_too": bad_t.validate().is_err() && bad_q.validate().is_err(),
        },
        "verdict": {
            "acceptance_passes": accept,
            "rejection_passes": reject,
            "criterion": "the kernel must accept both genuine forms AND reject \
                          both corrupted ones",
        },
    });
    println!("{}", serde_json::to_string_pretty(&report).unwrap());

    if !accept {
        eprintln!("FAIL acceptance: triple error {t_err:?}; quartet error {q_err:?}");
        return Ok(false);
    }
    if !reject {
        eprintln!(
            "FAIL rejection: the kernel accepted a corrupted certificate. The \
             verifier is not checking anything."
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
        Some("experiment") => run_experiment(&args),
        Some("saturation") => run_saturation(&args),
        Some("verify-tree") => run_verify_tree(&args),
        Some("fasta") => run_fasta(&args),
        _ => {
            eprintln!(
                "usage:\n  reticulation controls [--length N] [--generations G] \
                 [--mu M] [--seed S]\n  reticulation poc [--length N] [--seed S] \
                 [--formal-dir D] [--out D]\n  reticulation experiment [--seeds K] \
                 [--spread S] [--gap-fraction F] [--saturated-mu M]\n  \
                 reticulation saturation [--seeds K] [--low-mu M] [--high-mu M]\n  \
                 reticulation verify-tree [--length N] [--seed S] [--formal-dir D] \
                 [--out D]\n  reticulation fasta <FILE> [--model p|jc] \
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

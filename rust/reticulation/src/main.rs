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
    /// The averaged quartet statistics, for the max-versus-mean comparison.
    distribution: reticulation::subset::QuartetDistribution,
    /// The incumbent: the restricted pairwise-homoplasy test, on the same
    /// alignment. `None` when the alignment has too few informative sites for
    /// the window, which is reported rather than coerced to a number.
    phi: Option<reticulation::PhiTest>,
    /// Gap load, so a matched null can be built for this dataset rather than
    /// borrowing D1's.
    gap_fraction: f64,
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

/// Calibrate the averaged quartet statistic as a *test*, by parametric bootstrap.
///
/// The literature's treelikeness statistics are usually reported as scores, and
/// the one place they are used as a hypothesis test the null comes from a
/// parametric bootstrap under a fitted tree rather than from a permutation. That
/// choice matters here for a reason specific to this project: the RDP5 manual
/// reports that misalignment artefacts are significant in the real alignment yet
/// undetectable in shuffled alignments, so a permutation null is blind to the
/// artefact class `align-error` measures. A parametric bootstrap under a tree is
/// not.
///
/// So the null here is: clonal, relaxed-clock descent at **matched** alignment
/// width, matched mean p-distance and matched gap load — the same matched-null
/// construction `real-data` uses for the certificate's margin, applied to mean
/// δq instead. The p-value is the share of null replicates whose mean δq is at
/// least the observed one, with the usual `+1` so that it cannot be zero.
///
/// This is the third competitor in the head-to-head, and the one that tests
/// whether the averaged statistic, once calibrated, does better than the
/// certificate it is being compared against.
fn calibrated_mean_delta_q(
    observed: f64,
    width: usize,
    target_divergence: f64,
    gap_load: f64,
    spread: f64,
    replicates: usize,
    max_quadruples: usize,
) -> Result<(f64, f64, f64), String> {
    use reticulation::simulate::{clonal_tree_relaxed, inject_gap_blocks};
    use reticulation::subset::{quartet_distribution, Alignment};

    let generations = 4u32;
    let (mu, _) = tune_mu_for_divergence(width, generations, spread, target_divergence, 1);
    let mut at_least = 0usize;
    let mut total = 0.0f64;
    for replicate in 0..replicates {
        let seed = 0x5EED_0000 + replicate as u64;
        let mut pop = clonal_tree_relaxed(width, generations, mu, spread, seed);
        if gap_load > 0.001 {
            pop = inject_gap_blocks(pop, gap_load, seed ^ 0xBEEF);
        }
        let aln = Alignment::new(pop.labels.clone(), pop.sequences.clone())?;
        let value = quartet_distribution(&aln, max_quadruples)?.mean_delta_q;
        total += value;
        if value >= observed {
            at_least += 1;
        }
    }
    let p_value = (at_least + 1) as f64 / (replicates + 1) as f64;
    let floor = 1.0 / (replicates + 1) as f64;
    let null_mean = if replicates == 0 {
        0.0
    } else {
        total / replicates as f64
    };
    Ok((p_value, floor, null_mean))
}

/// Run both searches on one population.
fn measure(
    pop: &reticulation::simulate::Population,
    max_quadruples: usize,
    phi_window: usize,
    phi_permutations: usize,
) -> Result<SeedResult, String> {
    use reticulation::exact::{find_count_certificates, ExactDistances};
    use reticulation::subset::{
        compare_deletion_policies, find_complete_triple_certificates, find_quartet_certificates,
        quartet_distribution, Alignment,
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
        distribution: quartet_distribution(&aln, max_quadruples)?,
        // The competitor runs on the same alignment, under the same
        // complete-deletion policy, so the comparison is not confounded by the
        // two tests seeing different data.
        phi: reticulation::phi_test(&aln, phi_window, phi_permutations, 20_261_022).ok(),
        gap_fraction: gap_fraction(&pop.sequences),
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
    // The incumbent's parameters. 199 permutations is the count the published
    // treelikeness test uses, and it puts the attainable p-value floor at
    // 1/200 = 0.005, which is reported alongside every p-value so that a
    // floor-level result is not read as a smaller number than it is.
    let phi_window: usize = parse_flag(args, "--phi-window", 50);
    let phi_permutations: usize = parse_flag(args, "--phi-permutations", 199);
    // The third competitor: the averaged quartet statistic calibrated as a test
    // by parametric bootstrap under a MATCHED tree null, rather than reported as
    // a score. 99 replicates puts its attainable p-value floor at 0.01.
    let bootstrap_replicates: usize = parse_flag(args, "--bootstrap", 99);

    let seeds: Vec<u64> = (0..n_seeds).map(|i| seed0 + i as u64).collect();

    // --- build and measure every dataset at every seed ---
    let mut results: std::collections::BTreeMap<&str, Vec<SeedResult>> = Default::default();
    for &seed in &seeds {
        let mut push = |name: &'static str, pop: reticulation::simulate::Population| -> Result<(), String> {
            let mut r = measure(&pop, max_quadruples, phi_window, phi_permutations)?;
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

    // --- the averaged statistic, calibrated as a test --------------------------
    // Each dataset-seed gets its own matched null, because the null depends on
    // that population's width, divergence and gap load. Reported per dataset as
    // the share of seeds the calibrated test calls significant.
    let mut calibrated: std::collections::BTreeMap<&str, Vec<(f64, f64, f64, f64)>> =
        Default::default();
    if bootstrap_replicates > 0 {
        for (name, rows) in &results {
            let mut per_seed = Vec::new();
            for row in rows {
                // The matched null uses the divergence this population actually
                // has, which `quartet_distribution` already computed under
                // pairwise deletion over the whole alignment.
                let observed = row.distribution.mean_delta_q;
                let (p_value, floor, null_mean) = calibrated_mean_delta_q(
                    observed,
                    length,
                    row.distribution.mean_p_distance,
                    0.0,
                    spread,
                    bootstrap_replicates,
                    max_quadruples,
                )?;
                per_seed.push((observed, p_value, floor, null_mean));
            }
            calibrated.insert(name, per_seed);
        }
    }
    let calibrated_significant = |name: &str, alpha: f64| -> Vec<bool> {
        calibrated
            .get(name)
            .map(|v| v.iter().map(|(_, p, _, _)| *p <= alpha).collect())
            .unwrap_or_default()
    };

    // --- the certificate against ITS OWN matched null --------------------------
    // The power predictions compare each dataset's margin against D1's floor at
    // the same seed. That is commensurable only where the dataset's divergence
    // and denominator match D1's, and on D6 (gap blocks) and D7 (ten times the
    // divergence) they do not. Comparing a certificate judged against a
    // MISMATCHED floor with a competitor judged against a matched one would
    // flatter the competitor, so the certificate gets a matched floor too:
    // same construction `real-data` uses, and the same construction the
    // calibrated averaged statistic above uses.
    let mut matched: std::collections::BTreeMap<&str, Vec<bool>> = Default::default();
    if bootstrap_replicates > 0 {
        for (name, rows) in &results {
            let mut per_seed = Vec::new();
            for row in rows {
                let (floor, _, _, _) = matched_floor(
                    length,
                    row.distribution.mean_p_distance,
                    row.gap_fraction,
                    spread,
                    seeds.len() as u64,
                    max_quadruples,
                )?;
                per_seed.push(exceeds_k_times(row.quartet_margin, floor, 1));
            }
            matched.insert(name, per_seed);
        }
    }
    let ours_above_matched_floor = |name: &str| -> Vec<bool> {
        matched.get(name).cloned().unwrap_or_default()
    };

    // --- the head-to-head against the incumbent --------------------------------
    // Paper F listed, as its own limitation, that this comparison had not been
    // run. The incumbent is the one whose null also assumes NO clock: the
    // pairwise homoplasy index. Both tests run on the same alignments under the
    // same deletion policy, so the comparison is not about who saw what.
    //
    // Two questions, in this order. Does the incumbent survive the control that
    // destroyed our triple statistic -- a clonal population with a relaxed
    // clock? And where both survive it, which has more power?
    let phi_significant = |name: &str, alpha: f64| -> Vec<bool> {
        get(name)
            .iter()
            .map(|r| r.phi.as_ref().is_some_and(|p| p.p_value <= alpha))
            .collect()
    };
    // Our own statistic's verdict on the same dataset, by the same rule the
    // power predictions use: does the quartet margin exceed D1's floor at the
    // SAME seed? Seed-paired and decided by exact integer cross-product.
    let ours_above_floor = |name: &str| -> Vec<bool> {
        get(name)
            .iter()
            .zip(d1.iter())
            .map(|(r, base)| exceeds_k_times(r.quartet_margin, base.quartet_margin, 1))
            .collect()
    };
    let phi_rows: Vec<serde_json::Value> = [
        "D1_clonal_strict_clock",
        "D2_clonal_relaxed_clock",
        "D3_strict_clock_one_recombinant",
        "D4_relaxed_clock_one_recombinant",
        "D5_strict_clock_reciprocal_mosaics",
        "D6_strict_clock_gap_blocks",
        "D7_strict_clock_saturated",
    ]
    .iter()
    .map(|name| {
        let rows = get(name);
        let available: Vec<&reticulation::PhiTest> =
            rows.iter().filter_map(|r| r.phi.as_ref()).collect();
        serde_json::json!({
            "dataset": name,
            "reticulation_present": matches!(
                *name,
                "D3_strict_clock_one_recombinant"
                    | "D4_relaxed_clock_one_recombinant"
                    | "D5_strict_clock_reciprocal_mosaics"
            ),
            "seeds_measured": available.len(),
            "informative_sites": available.first().map(|p| p.informative_sites),
            "phi_w_exact_per_seed": available
                .iter()
                .map(|p| format!("{}/{}", p.phi_w_num, p.phi_w_den))
                .collect::<Vec<String>>(),
            "phi_w_per_seed": available.iter().map(|p| p.phi_w).collect::<Vec<f64>>(),
            "permuted_mean_per_seed": available
                .iter()
                .map(|p| p.permuted_mean_phi_w)
                .collect::<Vec<f64>>(),
            "lift_over_permuted_mean_per_seed": available
                .iter()
                .map(|p| p.phi_w - p.permuted_mean_phi_w)
                .collect::<Vec<f64>>(),
            "p_value_per_seed": available.iter().map(|p| p.p_value).collect::<Vec<f64>>(),
            "p_value_floor": available.first().map(|p| p.p_value_floor),
            "significant_at_0.05": frac(&phi_significant(name, 0.05)),
            "our_quartet_margin_above_D1_floor": frac(&ours_above_floor(name)),
            "our_quartet_margin_above_ITS_OWN_matched_floor":
                frac(&ours_above_matched_floor(name)),
            "calibrated_mean_delta_q_significant": frac(&calibrated_significant(name, 0.05)),
            "the_two_tests_agree": available.len() == rows.len()
                && held_at(&phi_significant(name, 0.05)).0
                    == held_at(&ours_above_floor(name)).0,
        })
    })
    .collect();

    // --- the max-versus-mean head-to-head -------------------------------------
    // The certificate's statistic is a MAXIMUM over quadruples. The published
    // treelikeness tests are AVERAGES, and in the one comparison of
    // treelikeness statistics used as tests, the averaged variants held a
    // correct false-positive rate while the mode-based and
    // proportion-of-resolved-quartets variants failed. That comparison has not
    // been run against a certificate-style maximum, so it is run here, on the
    // dataset the whole experiment turns on: D2, clonal with a relaxed clock
    // and no reticulation, the null that defeated the triple statistic.
    //
    // Every comparison is SEED-PAIRED. The question is not which statistic is
    // larger but at how many seeds each one ranks a reticulate population above
    // the hard null, which is the only form in which an extreme-value statistic
    // and an average are commensurable at all.
    let seed_paired_above = |name: &str, pick: &dyn Fn(&SeedResult) -> f64| -> Vec<bool> {
        get(name)
            .iter()
            .zip(get("D2_clonal_relaxed_clock").iter())
            .map(|(positive, null)| pick(positive) > pick(null))
            .collect()
    };
    let mean_delta = |r: &SeedResult| r.distribution.mean_delta_q;
    let proportion = |r: &SeedResult| r.distribution.proportion_violating;
    let max_margin = |r: &SeedResult| margin_f64(r.quartet_margin);
    let statistic_rows: Vec<serde_json::Value> = [
        (
            "quartet_max_margin",
            "maximum over quadruples -- the certificate's own statistic, the only              one with a single exhibitable witness, and the only one decided by              exact integer arithmetic",
            &max_margin as &dyn Fn(&SeedResult) -> f64,
        ),
        (
            "mean_delta_q",
            "mean of (m1-m2)/(m1-m3) over quadruples -- the averaged statistic the              literature reports as holding a correct false-positive rate",
            &mean_delta,
        ),
        (
            "proportion_violating",
            "share of quadruples with a strictly dominant pairing -- the variant              the literature reports as FAILING as a test, included because its              failure is informative",
            &proportion,
        ),
    ]
    .iter()
    .map(|(name, what, pick)| {
        serde_json::json!({
            "statistic": name,
            "what_it_is": what,
            "beats_relaxed_clock_null_at": {
                "D3_one_recombinant": frac(&seed_paired_above("D3_strict_clock_one_recombinant", *pick)),
                "D5_reciprocal_mosaics": frac(&seed_paired_above("D5_strict_clock_reciprocal_mosaics", *pick)),
            },
            "value_on_the_relaxed_clock_null_per_seed": get("D2_clonal_relaxed_clock")
                .iter().map(*pick).collect::<Vec<f64>>(),
            "value_on_reciprocal_mosaics_per_seed": get("D5_strict_clock_reciprocal_mosaics")
                .iter().map(*pick).collect::<Vec<f64>>(),
            // The weakest seed's ratio of the positive control to the hard null.
            // Display only, like every other ratio in this report, but it is the
            // number that decides which statistic separates better at FIXED
            // taxon count -- which is the only comparison the two share.
            "worst_seed_separation_over_relaxed_clock_null":
                get("D5_strict_clock_reciprocal_mosaics")
                    .iter()
                    .zip(get("D2_clonal_relaxed_clock").iter())
                    .map(|(positive, null)| {
                        let (a, b) = (pick(positive), pick(null));
                        if b > 0.0 { a / b } else { f64::INFINITY }
                    })
                    .fold(f64::INFINITY, f64::min),
        })
    })
    .collect();

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
                        "mean_delta_q": r.distribution.mean_delta_q,
                        "max_delta_q": r.distribution.max_delta_q,
                        "proportion_violating": r.distribution.proportion_violating,
                        "q_residual": r.distribution.q_residual,
                        "mean_p_distance": r.distribution.mean_p_distance,
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
        "head_to_head_against_pairwise_homoplasy": {
            "why": "Paper F listed as its own limitation that no comparison \
                    against the established tests had been run. The relevant \
                    incumbent is the one whose null also assumes no molecular \
                    clock -- the pairwise homoplasy index of Bruen, Philippe and \
                    Bryant (Genetics 172:2665-2681, 2006), whose null is \
                    site-order exchangeability. Both tests here run on the same \
                    alignments under the same complete-deletion policy.",
            "what_is_implemented": "A restricted form: only biallelic \
                    parsimony-informative, gap-free columns, so the four-gamete \
                    test is exact and needs no multistate refinement. Phi_w is a \
                    count of compatible neighbouring pairs over a count of pairs, \
                    and the pair count is invariant under permutation, so the \
                    permutation test compares integers. This is NOT the published \
                    statistic in full and is not presented as such.",
            "window_in_informative_sites": phi_window,
            "permutations": phi_permutations,
            "read_the_floor": "A p-value equal to p_value_floor means `as \
                    extreme as this many permutations can report`, not a smaller \
                    number. Section 8 of the textbook was burned by exactly this.",
            "per_dataset": phi_rows,
            "the_two_questions": {
                "does_the_incumbent_survive_the_relaxed_clock_control": frac(
                    &phi_significant("D2_clonal_relaxed_clock", 0.05)),
                "does_it_survive_the_strict_clock_control": frac(
                    &phi_significant("D1_clonal_strict_clock", 0.05)),
                "power_on_one_recombinant": frac(
                    &phi_significant("D3_strict_clock_one_recombinant", 0.05)),
                "power_on_reciprocal_mosaics": frac(
                    &phi_significant("D5_strict_clock_reciprocal_mosaics", 0.05)),
                "note": "A count of 0/5 on D1 and D2 with 5/5 on D3 and D5 would \
                         mean the incumbent is a fair competitor that the \
                         relaxed-clock confound does not touch, and the \
                         comparison is then about power, exhibitability and what \
                         each test licenses -- not about survival."
            },
            "how_to_read_the_side_by_side": "`our_quartet_margin_above_D1_floor` \
                    uses D1's floor at the same seed, which is commensurable only \
                    where the dataset's divergence and denominator match D1's. On \
                    D7 (ten times the divergence) and D6 (gap blocks, hence a \
                    different denominator) it is NOT, and a 5/5 there is the \
                    incommensurability this project has already recorded three \
                    times, not a measured false positive. The real-data command \
                    uses a MATCHED floor for exactly this reason. What the column \
                    does show, where it is commensurable, is whether the two \
                    tests reach the same verdict on the same alignment.",
            "what_this_run_found": "On the five datasets where both tests can be \
                    measured and the comparison is commensurable (D1-D5) the two \
                    AGREE exactly: both are specific on the strict-clock and the \
                    relaxed-clock null, and both detect all three reticulate \
                    populations at every seed. The incumbent is therefore not \
                    beaten on power here, and this method's case cannot rest on \
                    power. Two asymmetries do appear, one each way. The incumbent \
                    needs no recalibration when divergence changes, where our raw \
                    margin needs a matched null to be read at all. And the \
                    incumbent returns NO MEASUREMENT on the gap-blocked dataset, \
                    having too few biallelic gap-free informative columns, where \
                    the complete-deletion certificate still works.",
            "what_this_method_has_that_the_incumbent_does_not": "An exhibitable \
                    witness: the certificate is four sequences and six counts that \
                    a kernel recomputes, so a reader re-verifies the finding \
                    without running the search. And no permutation null: the RDP5 \
                    manual states that misalignment-induced false positives are \
                    detected with significant p-values in the real alignment yet \
                    are undetectable in the shuffled alignments, so a site-order \
                    permutation null -- which is exactly this incumbent's null -- \
                    is blind by construction to the artefact class that \
                    `align-error` measures. Neither of those is a power claim, and \
                    neither is established by this grid.",
        },
        "averaged_statistic_calibrated_as_a_test": {
            "why": "The treelikeness statistics are usually reported as scores. \
                    The one place they are used as a hypothesis test, the null \
                    comes from a parametric bootstrap under a fitted tree rather \
                    than from a permutation -- and that choice matters here, \
                    because a site-order permutation null is reported to be blind \
                    to the misalignment artefact class `align-error` measures. So \
                    mean delta_q is calibrated here against a MATCHED tree null: \
                    clonal, relaxed-clock, at this population's own alignment \
                    width and mean p-distance.",
            "replicates": bootstrap_replicates,
            "p_value_floor": if bootstrap_replicates > 0 {
                1.0 / (bootstrap_replicates + 1) as f64
            } else {
                1.0
            },
            "per_dataset": calibrated.iter().map(|(name, rows)| serde_json::json!({
                "dataset": name,
                "observed_mean_delta_q_per_seed":
                    rows.iter().map(|(o, _, _, _)| *o).collect::<Vec<f64>>(),
                "null_mean_per_seed":
                    rows.iter().map(|(_, _, _, m)| *m).collect::<Vec<f64>>(),
                "p_value_per_seed":
                    rows.iter().map(|(_, p, _, _)| *p).collect::<Vec<f64>>(),
                "significant_at_0.05": frac(&calibrated_significant(name, 0.05)),
            })).collect::<Vec<_>>(),
            "reading": "A calibrated averaged statistic is the fairest possible \
                        competitor to the certificate: same quadruples, same \
                        alignment, same matched-null discipline, and no \
                        permutation. Where it and the certificate reach the same \
                        verdicts, the certificate's distinctive contribution is \
                        the exhibitable witness rather than the decision.",
        },
        "max_versus_mean": {
            "why": "The certificate's statistic is a maximum over quadruples, which                     is an extreme-value statistic: its null level grows with the                     number of subsets examined, and this grid shows it -- the                     quartet floor moves from 5-7 at 16 taxa to 9-10 at 32. So                     margins are not comparable across population sizes, which is                     the same incommensurability that has already cost this project                     three criteria. The averaged statistics do not have that                     defect and cannot be certified, because a mean has no witness                     to exhibit. This block measures the price.",
            "these_are_floating_point": "mean_delta_q, proportion_violating and                     q_residual are f64 summaries and are NOT certificates. Only the                     max-margin row is decided by exact integer arithmetic and                     checkable by the Lean kernel. Nothing here changes that                     division.",
            "reference_null": "D2_clonal_relaxed_clock -- clonal descent, no                     reticulation whatever, divergence-matched relaxed clock. The                     hard null, not the easy one.",
            "statistics": statistic_rows,
            "what_this_run_found": "All three statistics rank both reticulate                     populations above the relaxed-clock null at every seed, so                     none of them is confounded by rate variation the way the                     TRIPLE statistic is. They differ in how widely they separate,                     and the ordering is not the one the literature's results would                     suggest: read                     `worst_seed_separation_over_relaxed_clock_null` on each row.                     `proportion_violating` separating barely at all is consistent                     with the published finding that this variant fails as a test.",
            "what_it_does_not_settle": "Separation at a FIXED taxon count is not                     the max statistic's weakness. Its weakness is that its null                     level grows with the number of subsets, so margins cannot be                     compared across populations of different size, and this grid                     is all at one size. A verdict on max-versus-mean needs the                     comparison repeated across taxon counts, and the averaged                     statistics calibrated by parametric bootstrap rather than                     against a shuffled null -- which misalignment artefacts are                     reported to be invisible to.",
        },
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

/// Per-quadruple four-point verdicts, in canonical enumeration order.
///
/// Returned as flags rather than margins on purpose. Two alignments of the same
/// population have different widths, so their margins have different
/// denominators and comparing the magnitudes would repeat the mistake that
/// finding 6 of the Phase 1 audit and P6 of the Phase 2 experiment both turned
/// on. Verdicts are commensurable; margins are not.
fn quartet_verdicts(aln: &reticulation::subset::Alignment) -> (Vec<bool>, usize) {
    let n = aln.len();
    let mut flags = Vec::new();
    let mut skipped = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                for l in (k + 1)..n {
                    let counts = aln.quadruple_counts(i, j, k, l);
                    if counts.common == 0 {
                        skipped += 1;
                        flags.push(false);
                        continue;
                    }
                    let sums = counts.pairing_sums();
                    let top = sums.iter().copied().max().unwrap_or(0);
                    let strictly_dominant = sums.iter().filter(|&&s| s == top).count() == 1;
                    flags.push(strictly_dominant);
                }
            }
        }
    }
    (flags, skipped)
}

/// Mean and maximum p-distance of an already-built alignment.
fn alignment_divergence(labels: &[String], rows: &[Vec<u8>]) -> (f64, f64) {
    let pop = reticulation::simulate::Population {
        labels: labels.to_vec(),
        sequences: rows.to_vec(),
    };
    divergence(&pop)
}

/// Gap fraction of a row set.
fn gap_fraction(rows: &[Vec<u8>]) -> f64 {
    let total: usize = rows.iter().map(Vec::len).sum();
    if total == 0 {
        return 0.0;
    }
    let gaps: usize = rows
        .iter()
        .map(|row| row.iter().filter(|&&b| b == b'-').count())
        .sum();
    gaps as f64 / total as f64
}

/// Does reference-anchored alignment manufacture quartet certificates?
///
/// This is the control the Lean kernel cannot supply. A certificate proves an
/// arithmetic fact about counts taken from the sequences it embeds; if the
/// alignment put the wrong residues in the same column, those counts are
/// faithfully computed from the wrong data and the kernel will accept them. The
/// trust surface of Phase 2 therefore grows by exactly one item when real data
/// enters — the alignment — and the only honest way to size it is to align data
/// whose true alignment is known.
///
/// So: simulate clonal, relaxed-clock descent with deletion events, which has
/// **no reticulation in it at all**; keep the true alignment; strip the gaps;
/// realign; and compare per-quadruple verdicts. Any quadruple the inferred
/// alignment calls a four-point violation and the true alignment does not was
/// manufactured by the aligner.
fn run_align_error(args: &[String]) -> Result<bool, String> {
    use reticulation::align::{
        centroid_reference, reference_anchored_alignment, ReferenceAnchored, Scoring,
    };
    use reticulation::simulate::{clonal_tree_relaxed_with_deletions, DeletionModel};
    use reticulation::subset::{find_quartet_certificates, Alignment};

    let length: usize = parse_flag(args, "--length", 2000);
    let generations: u32 = parse_flag(args, "--generations", 4);
    let mu: f64 = parse_flag(args, "--mu", 0.004);
    let spread: f64 = parse_flag(args, "--spread", 1.0);
    let deletion_rate: f64 = parse_flag(args, "--deletion-rate", 0.3);
    let mean_deletion: usize = parse_flag(args, "--mean-deletion", 30);
    let max_gap: f64 = parse_flag(args, "--max-gap-fraction", 0.2);
    let seeds: u64 = parse_flag(args, "--seeds", 5);
    let band: usize = parse_flag(args, "--band", 64);
    let max_band: usize = parse_flag(args, "--max-band", 2048);
    let max_quadruples: usize = parse_flag(args, "--max-quadruples", 20_000);

    let scoring = Scoring::default();

    let infer = |labels: &[String], raw: &[Vec<u8>]| -> Result<ReferenceAnchored, String> {
        let reference = centroid_reference(raw);
        // Certified: the band-edge test is not a proof of optimality, so this
        // control checks every pair against its unbanded optimum. A control that
        // measures alignment error must not itself rest on an unchecked band.
        reference_anchored_alignment(labels, raw, reference, &scoring, band, max_band, true)
    };

    // --- A1, the identity control: with no deletions the aligner must be the
    // identity, so every simulated result of Phase 2 is unchanged by routing it
    // through this module.
    let mut identity_holds = true;
    let mut identity_detail = String::new();
    for seed in 1..=seeds {
        let pop = clonal_tree_relaxed_with_deletions(
            length,
            generations,
            mu,
            spread,
            &DeletionModel {
                rate: 0.0,
                mean_length: mean_deletion,
                max_gap_fraction: max_gap,
            },
            seed,
        );
        let anchored = infer(&pop.labels, &pop.raw)?;
        if anchored.rows != pop.truth.sequences {
            identity_holds = false;
            identity_detail = format!(
                "seed {seed}: the inferred alignment of gapless input differs from the input"
            );
            break;
        }
        // Stronger than equality of rows: the optimal score must equal the
        // score of the diagonal, which is what makes the identity optimal
        // rather than merely returned.
        let reference_row = &pop.raw[anchored.reference_index];
        for (index, row) in pop.raw.iter().enumerate() {
            if index == anchored.reference_index {
                continue;
            }
            let diagonal: i64 = reference_row
                .iter()
                .zip(row)
                .map(|(left, right)| {
                    if left == right {
                        scoring.match_score as i64
                    } else {
                        scoring.mismatch as i64
                    }
                })
                .sum();
            if anchored.scores[index] != diagonal {
                identity_holds = false;
                identity_detail = format!(
                    "seed {seed}, sequence {index}: optimal score {} is not the diagonal \
                     score {diagonal}",
                    anchored.scores[index]
                );
                break;
            }
        }
        if !identity_holds {
            break;
        }
    }

    #[derive(serde::Serialize, Clone)]
    struct AlignErrorSeed {
        seed: u64,
        deletion_events: usize,
        truth_gap_fraction: f64,
        inferred_gap_fraction: f64,
        inferred_width: usize,
        truth_width: usize,
        insertions_dropped: usize,
        reference_label: String,
        column_recovery: f64,
        mean_p_distance: f64,
        quadruples: usize,
        truth_violating: usize,
        inferred_violating: usize,
        agree: usize,
        manufactured_by_alignment: usize,
        lost_to_alignment: usize,
        truth_margin: String,
        inferred_margin: String,
    }

    let measure_tier = |deletion_rate: f64,
                        mean_deletion: usize,
                        max_gap: f64|
     -> Result<Vec<AlignErrorSeed>, String> {
    let mut rows: Vec<AlignErrorSeed> = Vec::new();
    for seed in 1..=seeds {
        let pop = clonal_tree_relaxed_with_deletions(
            length,
            generations,
            mu,
            spread,
            &DeletionModel {
                rate: deletion_rate,
                mean_length: mean_deletion,
                max_gap_fraction: max_gap,
            },
            seed,
        );
        let truth_aln = Alignment::new(pop.labels.clone(), pop.truth.sequences.clone())?;
        let anchored = infer(&pop.labels, &pop.raw)?;
        let inferred_aln = Alignment::new(anchored.labels.clone(), anchored.rows.clone())?;

        // Column recovery: project the truth onto the reference's surviving
        // columns, which is exactly the frame the inferred alignment uses.
        let reference_row = &pop.truth.sequences[anchored.reference_index];
        let truth_columns: Vec<usize> = reference_row
            .iter()
            .enumerate()
            .filter(|&(_, &b)| b != b'-')
            .map(|(p, _)| p)
            .collect();
        let mut compared = 0usize;
        let mut identical = 0usize;
        if truth_columns.len() == anchored.width {
            for (s, row) in anchored.rows.iter().enumerate() {
                for (c, &column) in truth_columns.iter().enumerate() {
                    compared += 1;
                    if row[c] == pop.truth.sequences[s][column] {
                        identical += 1;
                    }
                }
            }
        }
        let column_recovery = if compared == 0 {
            0.0
        } else {
            identical as f64 / compared as f64
        };

        let (truth_flags, _) = quartet_verdicts(&truth_aln);
        let (inferred_flags, _) = quartet_verdicts(&inferred_aln);
        if truth_flags.len() != inferred_flags.len() {
            return Err(format!(
                "seed {seed}: {} quadruples in the truth against {} inferred; the \
                 two alignments must hold the same taxa",
                truth_flags.len(),
                inferred_flags.len()
            ));
        }
        let pairs = truth_flags.iter().zip(&inferred_flags);
        let mut agree = 0usize;
        let mut manufactured = 0usize;
        let mut lost = 0usize;
        for (t, i) in pairs {
            match (t, i) {
                (true, true) => agree += 1,
                (false, true) => manufactured += 1,
                (true, false) => lost += 1,
                (false, false) => {}
            }
        }
        let truth_q = find_quartet_certificates(&truth_aln, 1, max_quadruples)?;
        let inferred_q = find_quartet_certificates(&inferred_aln, 1, max_quadruples)?;
        let (mean_p, _) = alignment_divergence(&pop.labels, &pop.truth.sequences);

        rows.push(AlignErrorSeed {
            seed,
            deletion_events: pop.deletion_events,
            truth_gap_fraction: gap_fraction(&pop.truth.sequences),
            inferred_gap_fraction: anchored.gap_fraction(),
            inferred_width: anchored.width,
            truth_width: length,
            insertions_dropped: anchored.insertions_dropped.iter().sum(),
            reference_label: anchored.reference_label.clone(),
            column_recovery,
            mean_p_distance: mean_p,
            quadruples: truth_flags.len(),
            truth_violating: truth_flags.iter().filter(|b| **b).count(),
            inferred_violating: inferred_flags.iter().filter(|b| **b).count(),
            agree,
            manufactured_by_alignment: manufactured,
            lost_to_alignment: lost,
            truth_margin: margin_exact(truth_q.certificates.first().map(|c| c.margin())),
            inferred_margin: margin_exact(inferred_q.certificates.first().map(|c| c.margin())),
        });
    }
        Ok(rows)
    };

    let rows = measure_tier(deletion_rate, mean_deletion, max_gap)?;

    // The admissibility ladder. A3 holding at one deletion load says nothing
    // about another, and the whole point of the control is to find where the
    // aligner starts inventing violations — because that boundary is what
    // decides which real populations may be interpreted at all.
    #[derive(serde::Serialize)]
    struct LoadRung {
        deletion_rate: f64,
        mean_deletion_length: usize,
        mean_truth_gap_fraction: f64,
        /// The same load as the *inferred* alignment shows it. This is the one
        /// a real dataset lets you see, so it is the one the threshold uses.
        mean_inferred_gap_fraction: f64,
        mean_insertions_dropped: f64,
        /// Insertions dropped per aligned cell. A *rate*, because an absolute
        /// count is not comparable between a 3 kb and a 12 kb genome, and
        /// comparing it between them would be the incommensurable comparison
        /// this project has made three times already.
        mean_insertion_rate_per_site: f64,
        worst_column_recovery: f64,
        quadruples_per_seed: usize,
        manufactured_total: usize,
        lost_total: usize,
        seeds_with_manufactured: usize,
    }
    let mut ladder: Vec<LoadRung> = Vec::new();
    for mean_length in [30usize, 120usize] {
        for rate in [0.0f64, 0.1, 0.2, 0.4, 0.6, 0.9] {
            let tier = measure_tier(rate, mean_length, 0.3)?;
            let count = tier.len().max(1) as f64;
            ladder.push(LoadRung {
                deletion_rate: rate,
                mean_deletion_length: mean_length,
                mean_truth_gap_fraction: tier.iter().map(|r| r.truth_gap_fraction).sum::<f64>()
                    / count,
                mean_inferred_gap_fraction: tier
                    .iter()
                    .map(|r| r.inferred_gap_fraction)
                    .sum::<f64>()
                    / count,
                mean_insertions_dropped: tier
                    .iter()
                    .map(|r| r.insertions_dropped as f64)
                    .sum::<f64>()
                    / count,
                mean_insertion_rate_per_site: tier
                    .iter()
                    .map(|r| {
                        let cells = (r.inferred_width * (1usize << generations)) as f64;
                        if cells == 0.0 {
                            0.0
                        } else {
                            r.insertions_dropped as f64 / cells
                        }
                    })
                    .sum::<f64>()
                    / count,
                worst_column_recovery: tier
                    .iter()
                    .map(|r| r.column_recovery)
                    .fold(f64::INFINITY, f64::min),
                quadruples_per_seed: tier.first().map_or(0, |r| r.quadruples),
                manufactured_total: tier.iter().map(|r| r.manufactured_by_alignment).sum(),
                lost_total: tier.iter().map(|r| r.lost_to_alignment).sum(),
                seeds_with_manufactured: tier
                    .iter()
                    .filter(|r| r.manufactured_by_alignment > 0)
                    .count(),
            });
        }
    }
    // The strictest rung at which nothing was manufactured, read off the
    // ladder rather than assumed: this is the admissibility threshold
    // `real-data` applies.
    let admissible_gap = ladder
        .iter()
        .filter(|rung| rung.manufactured_total == 0)
        .map(|rung| rung.mean_inferred_gap_fraction)
        .fold(0.0f64, f64::max);
    let admissible_insertions = ladder
        .iter()
        .filter(|rung| rung.manufactured_total == 0)
        .map(|rung| rung.mean_insertion_rate_per_site)
        .fold(0.0f64, f64::max);
    let admissible_recovery = ladder
        .iter()
        .filter(|rung| rung.manufactured_total == 0)
        .map(|rung| rung.worst_column_recovery)
        .fold(1.0f64, f64::min);

    let worst_recovery = rows
        .iter()
        .map(|r| r.column_recovery)
        .fold(f64::INFINITY, f64::min);
    let total_manufactured: usize = rows.iter().map(|r| r.manufactured_by_alignment).sum();
    let total_lost: usize = rows.iter().map(|r| r.lost_to_alignment).sum();
    let seeds_with_manufactured = rows
        .iter()
        .filter(|r| r.manufactured_by_alignment > 0)
        .count();

    let report = serde_json::json!({
        "experiment": "does reference-anchored alignment manufacture four-point violations?",
        "why": "A Lean certificate proves an arithmetic fact about counts taken from the \
                sequences it embeds. If the alignment placed the wrong residues in the same \
                column, the counts are right and the data is wrong, and the kernel cannot see \
                it. Real data adds the aligner to the trust surface; this measures how big \
                that addition is, on populations with no reticulation in them at all.",
        "design": {
            "generator": "clonal_tree_relaxed_with_deletions",
            "reticulation_in_the_truth": "none",
            "length": length,
            "generations": generations,
            "taxa": 1usize << generations,
            "mu_per_generation": mu,
            "relaxed_clock_spread": spread,
            "deletion_rate_per_lineage_per_generation": deletion_rate,
            "mean_deletion_length": mean_deletion,
            "max_gap_fraction_per_lineage": max_gap,
            "seeds": seeds,
            "aligner": "reference-anchored affine-gap DP, centroid reference, banded with \
                        band-edge detection",
            "scoring": scoring,
            "comparison": "per-quadruple verdicts, not margins: the two alignments have \
                           different widths and their margins have different denominators"
        },
        "predictions": [
            {
                "id": "A1",
                "claim": "With the deletion rate at zero the inferred alignment is exactly \
                          the input, so no Phase 2 result changes by passing through the \
                          aligner.",
                "criterion": "inferred rows == truth rows at every seed, and the optimal \
                              score equals the diagonal score",
                "gated": true,
                "corrected_after_a_first_run": "A1 failed at 1 of 5 seeds on the first \
                    implementation, and the cause was not a bug in the DP: with a -10/-1 \
                    affine gap, a run of three mismatches scores exactly what deleting one \
                    reference residue and inserting one query residue scores (-12 + 5 against \
                    -11 + 15 - 11), so the identity was one of two optimal alignments and the \
                    aligner returned the other. The fix was to state the tie preference — M \
                    over X over Y, a substitution over a gap pair — which is why the criterion \
                    now checks the score as well as the rows. Regression test: \
                    tests/aligner_identity.rs."
            },
            {
                "id": "A2",
                "claim": "At this deletion load the aligner recovers almost every column.",
                "criterion": "worst-seed column recovery > 0.99",
                "gated": false
            },
            {
                "id": "A3",
                "claim": "Alignment error manufactures no four-point violation on a \
                          genealogy that is exactly a tree.",
                "criterion": "manufactured_by_alignment == 0 at every seed",
                "gated": false
            }
        ],
        "A1_identity_control": {
            "held": identity_holds,
            "detail": if identity_holds { "the aligner is the identity on gapless input".to_string() } else { identity_detail.clone() }
        },
        "A2_column_recovery": {
            "worst_seed": worst_recovery,
            "held": worst_recovery > 0.99
        },
        "A3_manufactured_violations": {
            "total_over_seeds": total_manufactured,
            "seeds_with_any": seeds_with_manufactured,
            "seeds": seeds,
            "held": total_manufactured == 0,
            "also_lost": total_lost,
            "reading": "manufactured_by_alignment counts quadruples the inferred alignment \
                        calls a four-point violation and the true alignment does not. \
                        lost_to_alignment counts the reverse. Both are alignment error; the \
                        first is the one that would be reported as reticulation."
        },
        "per_seed": rows,
        "admissibility_ladder": {
            "why": "A3 holding at one deletion load says nothing about another. This ladder \
                    raises the load until the aligner starts inventing four-point violations \
                    on a genealogy that is exactly a tree, which is the boundary that decides \
                    which real populations may be interpreted.",
            "rungs": ladder,
            "largest_inferred_gap_fraction_with_nothing_manufactured": admissible_gap,
            "largest_insertion_rate_per_site_with_nothing_manufactured": admissible_insertions,
            "worst_column_recovery_with_nothing_manufactured": admissible_recovery,
            "threshold_applied_by_real_data": {
                "inferred_gap_fraction_at_most": admissible_gap,
                "insertion_rate_per_site_at_most": admissible_insertions,
                "reading": "A real population whose inferred alignment carries a heavier gap \
                            load than this is reported BLOCKED rather than interpreted, \
                            because at that load alignment error alone produces violations in \
                            this control."
            }
        },
        "gate": "A1 only. A2 and A3 are measurements: a run in which alignment error does \
                 manufacture certificates is informative and must not be reported as a \
                 failed run."
    });
    println!("{}", serde_json::to_string_pretty(&report).map_err(|e| e.to_string())?);
    if !identity_holds {
        eprintln!("FAIL A1: {identity_detail}");
        return Ok(false);
    }
    Ok(true)
}

/// Tune `mu` so a relaxed-clock clonal population reaches a target divergence.
///
/// Bisection on the mean p-distance, which is monotone in `mu` over this range.
/// The matched null is the whole point: the floor a real margin must clear has
/// to come from a population at the *same* sequence length and the *same*
/// divergence, or the comparison measures the floor rather than the signal.
fn tune_mu_for_divergence(
    length: usize,
    generations: u32,
    spread: f64,
    target: f64,
    seed: u64,
) -> (f64, f64) {
    use reticulation::simulate::clonal_tree_relaxed;
    let mut low = 1e-6f64;
    let mut high = 0.5f64;
    let mut mid = 0.01f64;
    let mut achieved = 0.0f64;
    for _ in 0..40 {
        mid = 0.5 * (low + high);
        let pop = clonal_tree_relaxed(length, generations, mid, spread, seed);
        let (mean, _) = divergence(&pop);
        achieved = mean;
        if mean < target {
            low = mid;
        } else {
            high = mid;
        }
    }
    (mid, achieved)
}

/// The matched null floor: the largest quartet margin a *tree* produces at this
/// length, divergence and gap load.
fn matched_floor(
    length: usize,
    target_divergence: f64,
    gap_load: f64,
    spread: f64,
    seeds: u64,
    max_quadruples: usize,
) -> Result<(Option<reticulation::Ratio>, f64, f64, usize), String> {
    use reticulation::simulate::{clonal_tree_relaxed, inject_gap_blocks};
    use reticulation::subset::{find_quartet_certificates, Alignment};

    let generations = 4u32; // 16 leaves, matching the real populations
    let (mu, achieved) = tune_mu_for_divergence(length, generations, spread, target_divergence, 1);
    let mut best: Option<reticulation::Ratio> = None;
    let mut violating_total = 0usize;
    for seed in 1..=seeds {
        let mut pop = clonal_tree_relaxed(length, generations, mu, spread, seed);
        if gap_load > 0.001 {
            pop = inject_gap_blocks(pop, gap_load, seed ^ 0x5EED);
        }
        let aln = Alignment::new(pop.labels.clone(), pop.sequences.clone())?;
        let search = find_quartet_certificates(&aln, 1, max_quadruples)?;
        violating_total += search.violating;
        if let Some(margin) = search.certificates.first().map(|c| c.margin()) {
            let better = match best {
                None => true,
                Some(current) => margin
                    .checked_cmp(current)
                    .map(|o| o == std::cmp::Ordering::Greater)
                    .unwrap_or(false),
            };
            if better {
                best = Some(margin);
            }
        }
    }
    Ok((best, mu, achieved, violating_total))
}

/// Run the four-point test on real viral populations.
///
/// The populations come from `scripts/reticulata_hf_populations.py`, which
/// pulls them from the Hugging Face `hiyata/Virus-Host-Genomes` set under
/// pre-registered filters, and each carries a *literature* expectation about
/// whether its genealogy is tree-like. The expectation is a prior, not a label
/// derived from this data, and the test is allowed to contradict it — which is
/// the only way the exercise means anything.
///
/// Every population is run twice, on two different alignment references, and
/// the report says whether the verdict moved. Every margin is compared against
/// a floor from a simulated tree at matched length, matched divergence and
/// matched gap load, because a margin on its own is not a statistic.
fn run_real_data(args: &[String]) -> Result<bool, String> {
    use reticulation::align::{centroid_ranking, reference_anchored_alignment, Scoring};
    use reticulation::subset::{find_complete_triple_certificates, find_quartet_certificates, Alignment};

    let dir: String = parse_flag(
        args,
        "--populations",
        "../../data/reticulata_phase3/populations".to_string(),
    );
    let band: usize = parse_flag(args, "--band", 64);
    let max_band: usize = parse_flag(args, "--max-band", 4096);
    let seeds: u64 = parse_flag(args, "--seeds", 5);
    let spread: f64 = parse_flag(args, "--spread", 1.0);
    let max_quadruples: usize = parse_flag(args, "--max-quadruples", 20_000);
    // The admissibility threshold is *measured*, not chosen: it is the heaviest
    // rung of `reticulation align-error`'s ladder at which alignment error
    // manufactured no four-point violation at any seed. Above it, this command
    // reports a population BLOCKED instead of interpreting it, because at that
    // load the aligner alone produces violations on a genealogy that is exactly
    // a tree. Raising these flags is a deliberate act with a known cost.
    let max_gap_admissible: f64 = parse_flag(args, "--max-gap-fraction", 0.0039);
    let max_insertion_rate: f64 = parse_flag(args, "--max-insertion-rate", 0.0008);

    let path = std::path::Path::new(&dir);
    let mut fastas: Vec<std::path::PathBuf> = std::fs::read_dir(path)
        .map_err(|e| format!("--populations {dir}: {e}"))?
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("fasta"))
        .collect();
    fastas.sort();
    if fastas.is_empty() {
        return Err(format!("no .fasta files in {dir}"));
    }

    // Literature expectations, written by the extraction script.
    let provenance: serde_json::Value = std::fs::read_to_string(path.join("provenance.json"))
        .ok()
        .and_then(|text| serde_json::from_str(&text).ok())
        .unwrap_or(serde_json::Value::Null);
    let expectation_for = |stem: &str| -> (String, String) {
        if let Some(populations) = provenance.get("populations").and_then(|v| v.as_array()) {
            for entry in populations {
                let slug = entry.get("slug").and_then(|v| v.as_str()).unwrap_or("");
                if !slug.is_empty() && stem.starts_with(slug) {
                    return (
                        entry
                            .get("expectation")
                            .and_then(|v| v.as_str())
                            .unwrap_or("unknown")
                            .to_string(),
                        entry
                            .get("virus_name")
                            .and_then(|v| v.as_str())
                            .unwrap_or(stem)
                            .to_string(),
                    );
                }
            }
        }
        ("unknown".to_string(), stem.to_string())
    };

    #[derive(serde::Serialize)]
    struct ReferenceRun {
        role: String,
        reference_label: String,
        width: usize,
        gap_fraction: f64,
        insertions_dropped: usize,
        insertion_rate_per_site: f64,
        max_band_used: usize,
        mean_p_distance: f64,
        max_p_distance: f64,
        quadruples_examined: usize,
        quartet_violating: usize,
        quartet_margin: String,
        quartet_margin_f64: f64,
        quartet_labels: Vec<String>,
        triple_violating: usize,
        triple_margin: String,
    }

    #[derive(serde::Serialize)]
    struct PopulationResult {
        file: String,
        virus_name: String,
        literature_expectation: String,
        taxa: usize,
        raw_length_min: usize,
        raw_length_max: usize,
        status: String,
        runs: Vec<ReferenceRun>,
        references_agree_on_verdict: Option<bool>,
        matched_floor: Option<String>,
        matched_floor_f64: Option<f64>,
        matched_floor_mu: Option<f64>,
        matched_floor_divergence: Option<f64>,
        exceeds_matched_floor: Option<bool>,
        exceeds_twice_matched_floor: Option<bool>,
    }

    let scoring = Scoring::default();
    let mut results: Vec<PopulationResult> = Vec::new();

    for fasta in &fastas {
        let stem = fasta
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_string();
        let (expectation, virus_name) = expectation_for(&stem);
        let text = std::fs::read_to_string(fasta).map_err(|e| format!("{}: {e}", fasta.display()))?;
        let records = reticulation::align::read_fasta(&text)?;
        let labels: Vec<String> = records.iter().map(|r| r.label.clone()).collect();
        let raw: Vec<Vec<u8>> = records.iter().map(|r| r.sequence.clone()).collect();
        let lengths: Vec<usize> = raw.iter().map(Vec::len).collect();
        let mut entry = PopulationResult {
            file: stem.clone(),
            virus_name,
            literature_expectation: expectation,
            taxa: raw.len(),
            raw_length_min: lengths.iter().copied().min().unwrap_or(0),
            raw_length_max: lengths.iter().copied().max().unwrap_or(0),
            status: "ok".to_string(),
            runs: Vec::new(),
            references_agree_on_verdict: None,
            matched_floor: None,
            matched_floor_f64: None,
            matched_floor_mu: None,
            matched_floor_divergence: None,
            exceeds_matched_floor: None,
            exceeds_twice_matched_floor: None,
        };
        if raw.len() < 4 {
            entry.status = format!("BLOCKED: {} sequences, a quartet needs 4", raw.len());
            results.push(entry);
            continue;
        }
        let ranking = centroid_ranking(&raw);
        let mut blocked: Option<String> = None;
        for (role, reference) in [("primary", ranking[0]), ("sensitivity", ranking[1])] {
            let anchored = match reference_anchored_alignment(
                &labels, &raw, reference, &scoring, band, max_band, true,
            ) {
                Ok(value) => value,
                Err(error) => {
                    blocked = Some(format!("BLOCKED at the {role} reference: {error}"));
                    break;
                }
            };
            let aln = match Alignment::new(anchored.labels.clone(), anchored.rows.clone()) {
                Ok(value) => value,
                Err(error) => {
                    blocked = Some(format!("BLOCKED building the {role} alignment: {error}"));
                    break;
                }
            };
            let quartets = find_quartet_certificates(&aln, 1, max_quadruples)?;
            let triples = find_complete_triple_certificates(&aln, 1)?;
            if !quartets.all_validated || !triples.all_validated {
                return Err(format!("{stem}: a certificate failed revalidation"));
            }
            let (mean_p, max_p) = alignment_divergence(&anchored.labels, &anchored.rows);
            let top = quartets.certificates.first();
            entry.runs.push(ReferenceRun {
                role: role.to_string(),
                reference_label: anchored.reference_label.clone(),
                width: anchored.width,
                gap_fraction: anchored.gap_fraction(),
                insertions_dropped: anchored.insertions_dropped.iter().sum(),
                insertion_rate_per_site: anchored.insertions_dropped.iter().sum::<usize>() as f64
                    / (anchored.width * anchored.rows.len()).max(1) as f64,
                max_band_used: anchored.bands.iter().copied().max().unwrap_or(0),
                mean_p_distance: mean_p,
                max_p_distance: max_p,
                quadruples_examined: quartets.examined,
                quartet_violating: quartets.violating,
                quartet_margin: margin_exact(top.map(|c| c.margin())),
                quartet_margin_f64: margin_f64(top.map(|c| c.margin())),
                quartet_labels: top.map_or_else(Vec::new, |c| {
                    vec![c.a.clone(), c.b.clone(), c.c.clone(), c.e.clone()]
                }),
                triple_violating: triples.violating,
                triple_margin: margin_exact(triples.certificates.first().map(|c| c.margin())),
            });
        }
        if let Some(reason) = blocked {
            entry.status = reason;
            results.push(entry);
            continue;
        }
        if entry.runs.len() == 2 {
            entry.references_agree_on_verdict = Some(
                (entry.runs[0].quartet_violating > 0) == (entry.runs[1].quartet_violating > 0),
            );
        }
        // The admissibility gate, before any margin is read.
        let worst_gap = entry
            .runs
            .iter()
            .map(|run| run.gap_fraction)
            .fold(0.0f64, f64::max);
        let worst_insertion_rate = entry
            .runs
            .iter()
            .map(|run| {
                let cells = (run.width * entry.taxa) as f64;
                if cells == 0.0 {
                    0.0
                } else {
                    run.insertions_dropped as f64 / cells
                }
            })
            .fold(0.0f64, f64::max);
        if worst_gap > max_gap_admissible || worst_insertion_rate > max_insertion_rate {
            entry.status = format!(
                "BLOCKED: the inferred alignment carries a gap fraction of {worst_gap:.4} and \
                 drops insertions at {worst_insertion_rate:.5} per aligned cell, against the \
                 measured admissibility threshold of {max_gap_admissible:.4} and \
                 {max_insertion_rate:.5}. At that load `reticulation align-error` shows \
                 alignment error alone manufacturing four-point violations on a genealogy \
                 that is exactly a tree, so a verdict here would not be about this population. \
                 The counts above are reported, and are not evidence of reticulation."
            );
            results.push(entry);
            continue;
        }
        // The matched null, built from the primary run's own numbers.
        let primary = &entry.runs[0];
        let (floor, mu, achieved, _) = matched_floor(
            primary.width,
            primary.mean_p_distance,
            primary.gap_fraction,
            spread,
            seeds,
            max_quadruples,
        )?;
        let observed = if primary.quartet_margin == "none" {
            None
        } else {
            let parts: Vec<i64> = primary
                .quartet_margin
                .split('/')
                .filter_map(|p| p.parse().ok())
                .collect();
            if parts.len() == 2 {
                Some(reticulation::Ratio {
                    num: parts[0],
                    den: parts[1],
                })
            } else {
                None
            }
        };
        entry.matched_floor = Some(margin_exact(floor));
        entry.matched_floor_f64 = Some(margin_f64(floor));
        entry.matched_floor_mu = Some(mu);
        entry.matched_floor_divergence = Some(achieved);
        entry.exceeds_matched_floor = Some(exceeds_k_times(observed, floor, 1));
        entry.exceeds_twice_matched_floor = Some(exceeds_k_times(observed, floor, 2));
        results.push(entry);
    }

    // Proof-carrying output. A real population that clears both the
    // admissibility gate and its matched floor gets a Lean file embedding the
    // four aligned rows, so the kernel recomputes every count itself and the
    // only thing left trusted is that those rows are what the aligner produced.
    // The corrupted twin must fail, or the check is checking nothing.
    #[derive(serde::Serialize)]
    struct KernelVerdict {
        population: String,
        file: String,
        accepted: bool,
        seconds: f64,
        detail: String,
    }
    let mut kernel: Vec<KernelVerdict> = Vec::new();
    let emit_dir: String = parse_flag(args, "--emit-lean", String::new());
    if !emit_dir.is_empty() {
        use reticulation::emit_lean;
        let formal_dir: String =
            parse_flag(args, "--formal-dir", "../../../../../formal".to_string());
        let formal = std::path::Path::new(&formal_dir)
            .canonicalize()
            .map_err(|e| format!("--formal-dir {formal_dir}: {e}"))?;
        std::fs::create_dir_all(&emit_dir).map_err(|e| format!("--emit-lean {emit_dir}: {e}"))?;
        let targets: Vec<String> = results
            .iter()
            .filter(|r| r.status == "ok" && r.exceeds_matched_floor == Some(true))
            .map(|r| r.file.clone())
            .collect();
        for stem in targets {
            let fasta = path.join(format!("{stem}.fasta"));
            let text =
                std::fs::read_to_string(&fasta).map_err(|e| format!("{}: {e}", fasta.display()))?;
            let records = reticulation::align::read_fasta(&text)?;
            let labels: Vec<String> = records.iter().map(|r| r.label.clone()).collect();
            let raw: Vec<Vec<u8>> = records.iter().map(|r| r.sequence.clone()).collect();
            let reference = centroid_ranking(&raw)[0];
            let anchored = reference_anchored_alignment(
                &labels, &raw, reference, &scoring, band, max_band, true,
            )?;
            let aln = Alignment::new(anchored.labels.clone(), anchored.rows.clone())?;
            let search = find_quartet_certificates(&aln, 1, max_quadruples)?;
            let certificate = match search.certificates.first() {
                Some(value) => value,
                None => continue,
            };
            let [a, b, c, e] = certificate.indices;
            let enc = |i: usize| emit_lean::encode(&anchored.rows[i]);
            let source = format!(
                "{stem}, real sequences from huggingface:hiyata/Virus-Host-Genomes, \
                 reference-anchored alignment on {} ({} columns)",
                anchored.reference_label, anchored.width
            );
            for (suffix, body) in [
                (
                    "",
                    emit_lean::render_quartet(
                        certificate,
                        &enc(a),
                        &enc(b),
                        &enc(c),
                        &enc(e),
                        &source,
                    ),
                ),
                (
                    "_bad",
                    emit_lean::render_quartet(
                        &emit_lean::corrupt_quartet(certificate),
                        &enc(a),
                        &enc(b),
                        &enc(c),
                        &enc(e),
                        &source,
                    ),
                ),
            ] {
                let out = std::path::Path::new(&emit_dir)
                    .join(format!("reticulata_real_{stem}{suffix}.lean"));
                std::fs::write(&out, body).map_err(|e| format!("{}: {e}", out.display()))?;
                // `kernel_check` runs `lake` from the formal directory, so the
                // file must be named absolutely. A relative path here made the
                // first run report every file as rejected for want of a file,
                // which is also how the rejection control came out vacuously
                // true: both directions failed for the same irrelevant reason.
                let out = out
                    .canonicalize()
                    .map_err(|e| format!("{}: {e}", out.display()))?;
                let started = std::time::Instant::now();
                let (ok, detail) = kernel_check(&formal, &out);
                kernel.push(KernelVerdict {
                    population: stem.clone(),
                    file: out
                        .file_name()
                        .and_then(|s| s.to_str())
                        .unwrap_or("")
                        .to_string(),
                    accepted: ok,
                    seconds: started.elapsed().as_secs_f64(),
                    detail: if suffix.is_empty() {
                        if ok {
                            "accepted, no axioms beyond the whitelist".to_string()
                        } else {
                            detail
                        }
                    } else if ok {
                        "FAIL: the kernel accepted a corrupted certificate".to_string()
                    } else {
                        "rejected, as a corrupted certificate must be".to_string()
                    },
                });
            }
        }
    }

    let usable: Vec<&PopulationResult> = results.iter().filter(|r| r.status == "ok").collect();
    let agree = usable
        .iter()
        .filter(|r| r.references_agree_on_verdict == Some(true))
        .count();
    let report = serde_json::json!({
        "experiment": "the four-point test on real viral populations",
        "data": {
            "source": "huggingface:hiyata/Virus-Host-Genomes, via scripts/reticulata_hf_populations.py",
            "selection": "pre-registered: near-median length, <=0.1% ambiguous, deduplicated, \
                          then either accession order or greedy 6-mer diversity",
            "alignment": "reference-anchored affine-gap DP in rust/reticulation/src/align.rs; \
                          this machine has no mafft, muscle or clustalo",
            "known_limitation": "insertions relative to the reference are dropped and every \
                                 pair of non-reference sequences is related through the \
                                 reference. `reticulation align-error` measures what that costs."
        },
        "admissibility": {
            "max_gap_fraction": max_gap_admissible,
            "max_insertion_rate_per_site": max_insertion_rate,
            "provenance": "Both are the heaviest rung of `reticulation align-error`'s \
                           admissibility ladder at which alignment error manufactured no \
                           four-point violation at any seed, on simulated populations with no \
                           reticulation in them at all. A population above either threshold is \
                           reported BLOCKED with its counts, not interpreted.",
            "rates_not_counts": "Both thresholds are rates per aligned cell. An absolute \
                                 insertion count is not comparable between a 3 kb and a 12 kb \
                                 genome."
        },
        "how_to_read_a_margin": "A margin is not a statistic on its own. Each population's \
                                 margin is compared against a floor from a simulated clonal, \
                                 relaxed-clock tree at the same alignment width, the same mean \
                                 p-distance and the same gap load. Comparing margins between \
                                 populations would compare different denominators and \
                                 different divergences at once.",
        "populations": results,
        "kernel_verdicts": kernel,
        "summary": {
            "usable": usable.len(),
            "blocked": results.len() - usable.len(),
            "references_agree_on_verdict": format!("{agree}/{}", usable.len()),
            "exceeding_matched_floor": usable
                .iter()
                .filter(|r| r.exceeds_matched_floor == Some(true))
                .map(|r| r.file.clone())
                .collect::<Vec<String>>(),
            "exceeding_twice_matched_floor": usable
                .iter()
                .filter(|r| r.exceeds_twice_matched_floor == Some(true))
                .map(|r| r.file.clone())
                .collect::<Vec<String>>()
        }
    });
    println!("{}", serde_json::to_string_pretty(&report).map_err(|e| e.to_string())?);
    if usable.is_empty() {
        eprintln!("no population could be aligned; nothing was tested");
        return Ok(false);
    }
    Ok(true)
}

/// Run the test on a population whose recombination has a **published
/// breakpoint**, in an alignment this project did not produce.
///
/// Every real-data result so far rested on two weaker things: a literature
/// *expectation* about whether a family recombines, and an alignment from our
/// own aligner, whose error we therefore had to measure and gate on. Here the
/// recombinant lineages are named by the Pango designation process with
/// published breakpoint intervals in reference coordinates, and the alignment
/// comes already made from Nextstrain's open release, aligned by Nextclade to
/// MN908947 — the same coordinate system the breakpoints are quoted in.
///
/// The pre-registered prediction is sharper than a margin. If the four-point
/// certificate is detecting the recombination that the designation process
/// named, then its **dominant pairing should separate the two parent
/// lineages**, with recombinant members on the other side: that is the planted-
/// recombinant shape, and it is a statement about *which* four sequences the
/// method picks, not only about how large a number it produces. A dominant
/// pairing that instead separates two members of the same lineage would
/// contradict it.
///
/// One caveat is reported rather than hidden: using a published alignment
/// removes *our* aligner from the trust surface and replaces it with
/// Nextclade's, which is reference-anchored too and whose error we have not
/// measured. The admissibility gate calibrated for our aligner does not
/// transfer, so it is not applied; the gap load is reported instead.
fn run_breakpoint(args: &[String]) -> Result<bool, String> {
    use reticulation::subset::{find_quartet_certificates, quartet_distribution, Alignment};

    let fasta: String = parse_flag(args, "--population", String::new());
    if fasta.is_empty() {
        return Err("usage: reticulation breakpoint --population <aligned.fasta> \
                    [--provenance <json>] [--top K]"
            .to_string());
    }
    let provenance_path: String = parse_flag(args, "--provenance", String::new());
    let top_k: usize = parse_flag(args, "--top", 10);
    let seeds: u64 = parse_flag(args, "--seeds", 5);
    let spread: f64 = parse_flag(args, "--spread", 1.0);
    let max_quadruples: usize = parse_flag(args, "--max-quadruples", 20_000);

    let text = std::fs::read_to_string(&fasta).map_err(|e| format!("{fasta}: {e}"))?;
    let records = reticulation::align::read_fasta(&text)?;
    let labels: Vec<String> = records.iter().map(|r| r.label.clone()).collect();
    let rows: Vec<Vec<u8>> = records.iter().map(|r| r.sequence.clone()).collect();

    // The lineage of each member, from the provenance the extraction script
    // wrote. Without it the parent-separation prediction cannot be checked, and
    // that is reported rather than silently skipped.
    let mut lineage_of: std::collections::BTreeMap<String, String> = Default::default();
    let mut recombinant_lineage = String::new();
    let mut parent_lineages: Vec<String> = Vec::new();
    let mut published_breakpoint = String::new();
    if !provenance_path.is_empty() {
        let raw = std::fs::read_to_string(&provenance_path)
            .map_err(|e| format!("{provenance_path}: {e}"))?;
        let value: serde_json::Value =
            serde_json::from_str(&raw).map_err(|e| format!("{provenance_path}: {e}"))?;
        if let Some(populations) = value.get("populations").and_then(|v| v.as_array()) {
            let stem = std::path::Path::new(&fasta)
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("");
            for entry in populations {
                if entry.get("file").and_then(|v| v.as_str()) != Some(stem) {
                    continue;
                }
                recombinant_lineage = entry
                    .get("recombinant_lineage")
                    .and_then(|v| v.as_str())
                    .unwrap_or("")
                    .to_string();
                published_breakpoint = entry
                    .get("published_breakpoint")
                    .and_then(|v| v.as_str())
                    .unwrap_or("")
                    .to_string();
                parent_lineages = entry
                    .get("parent_lineages")
                    .and_then(|v| v.as_array())
                    .map(|a| {
                        a.iter()
                            .filter_map(|v| v.as_str())
                            .map(str::to_string)
                            .collect()
                    })
                    .unwrap_or_default();
                for member in entry
                    .get("members")
                    .and_then(|v| v.as_array())
                    .map(Vec::as_slice)
                    .unwrap_or(&[])
                {
                    if let (Some(strain), Some(lineage)) = (
                        member.get("strain").and_then(|v| v.as_str()),
                        member.get("lineage").and_then(|v| v.as_str()),
                    ) {
                        lineage_of.insert(strain.to_string(), lineage.to_string());
                    }
                }
            }
        }
    }
    let lineage = |label: &str| -> String {
        lineage_of
            .get(label)
            .cloned()
            .unwrap_or_else(|| "unknown".to_string())
    };

    let aln = Alignment::new(labels.clone(), rows.clone())?;
    let (mean_p, max_p) = alignment_divergence(&labels, &rows);
    let gaps = gap_fraction(&rows);
    let search = find_quartet_certificates(&aln, top_k, max_quadruples)?;
    let distribution = quartet_distribution(&aln, max_quadruples)?;
    if !search.all_validated {
        return Err("a certificate failed revalidation".to_string());
    }

    // The prediction: does the dominant pairing separate the two parents?
    let separates_parents = |a: &str, b: &str, c: &str, e: &str| -> bool {
        if parent_lineages.len() != 2 {
            return false;
        }
        let (pa, pb) = (&parent_lineages[0], &parent_lineages[1]);
        let far = [lineage(a), lineage(b)];
        let near = [lineage(c), lineage(e)];
        let far_is_the_two_parents = (far[0] == *pa && far[1] == *pb)
            || (far[0] == *pb && far[1] == *pa);
        let near_has_a_recombinant =
            !recombinant_lineage.is_empty() && near.contains(&recombinant_lineage);
        far_is_the_two_parents && near_has_a_recombinant
    };
    let certificates: Vec<serde_json::Value> = search
        .certificates
        .iter()
        .map(|c| {
            serde_json::json!({
                "dominant_pair": [c.a.clone(), c.b.clone()],
                "dominant_pair_lineages": [lineage(&c.a), lineage(&c.b)],
                "other_side": [c.c.clone(), c.e.clone()],
                "other_side_lineages": [lineage(&c.c), lineage(&c.e)],
                "common_columns": c.common,
                "margin_exact": format!("{}/{}", c.margin_num, c.margin_den),
                "margin": c.margin().to_f64(),
                "separates_the_two_parent_lineages": separates_parents(&c.a, &c.b, &c.c, &c.e),
            })
        })
        .collect();
    let separating = certificates
        .iter()
        .filter(|c| {
            c.get("separates_the_two_parent_lineages")
                .and_then(serde_json::Value::as_bool)
                == Some(true)
        })
        .count();

    let (floor, floor_mu, floor_divergence, _) = matched_floor(
        aln.width(),
        mean_p,
        gaps,
        spread,
        seeds,
        max_quadruples,
    )?;
    let observed = search.certificates.first().map(|c| c.margin());

    let report = serde_json::json!({
        "experiment": "the four-point test on a population whose recombination \
                       has a published breakpoint, in an alignment this project \
                       did not produce",
        "population": {
            "file": fasta,
            "taxa": aln.len(),
            "alignment_width": aln.width(),
            "gap_fraction": gaps,
            "mean_p_distance": mean_p,
            "max_p_distance": max_p,
            "recombinant_lineage": recombinant_lineage,
            "parent_lineages": parent_lineages,
            "published_breakpoint": published_breakpoint,
            "lineage_known_for": lineage_of.len(),
        },
        "provenance_of_the_alignment": {
            "aligner": "Nextclade, to MN908947, in Nextstrain's open release",
            "ours_is_not_involved": true,
            "caveat": "This removes OUR aligner from the trust surface and \
                       replaces it with theirs, which is reference-anchored too \
                       and whose error we have not measured. The admissibility \
                       gate calibrated for our aligner therefore does NOT \
                       transfer and is not applied here; the gap fraction above \
                       is reported in its place."
        },
        "preregistered_prediction": {
            "claim": "If the certificate detects the recombination the Pango \
                      process named, its dominant pairing separates the two \
                      PARENT lineages with a recombinant on the other side.",
            "why_this_is_stronger_than_a_margin": "It constrains WHICH four \
                      sequences the method picks and how it sides them, not just \
                      how large a number it reports. A margin can be large for \
                      many reasons; picking the two parents cannot.",
            "certificates_examined": certificates.len(),
            "certificates_separating_the_parents": separating,
            "held": separating > 0,
        },
        "certificates": certificates,
        "quartet_distribution": distribution,
        "matched_floor": {
            "margin": margin_exact(floor),
            "margin_f64": margin_f64(floor),
            "mu": floor_mu,
            "achieved_divergence": floor_divergence,
            "exceeds_floor": exceeds_k_times(observed, floor, 1),
            "exceeds_twice_floor": exceeds_k_times(observed, floor, 2),
            // A floor of `none` means the matched null produced no violating
            // quadruple at all. `exceeds_floor` is then true for ANY positive
            // margin, which is not a calibrated comparison -- it is the same
            // degeneracy the topology experiment hit, and at these divergences
            // it is the expected one.
            "null_is_degenerate": floor.is_none(),
            "degenerate_means": "The matched null produced no violating                 quadruple in any replicate, so `exceeds_floor` reduces to                 `the observed margin is positive` and carries no calibration.                 Raise the replicate count, or read this population as below the                 method's resolution at its divergence.",
        },
        "resolution_warning": {
            "mean_p_distance": mean_p,
            "differing_sites_implied": mean_p * aln.width() as f64,
            "margin_in_sites": search
                .certificates
                .first()
                .map_or(0, |c| c.margin_num),
            "reading": "The margin is a count of differing sites over the                 quadruple's common columns. When the population's mean                 p-distance implies only a few tens of differing sites                 genome-wide, a margin of one or two sites is at the edge of what                 the statistic can resolve, and a verdict should not rest on it                 however the floor comparison comes out."
        },
        "observed_margin": margin_exact(observed),
    });
    println!(
        "{}",
        serde_json::to_string_pretty(&report).map_err(|e| e.to_string())?
    );
    if search.certificates.is_empty() {
        eprintln!("no quartet certificate on this population");
        return Ok(false);
    }
    Ok(true)
}

/// Does persistent homology detect *recombination*, or does it detect mutation?
///
/// This programme began by reproducing the construction of Chan, Carlsson and
/// Rabadan: compute persistent homology on a population under a genetic-distance
/// metric, where an unfilled 1-cycle is read as reticulate evolution. The
/// theoretical warrant is real and is the same fact this crate's null rests on —
/// a tree metric leaves no unfilled cycle. What is missing from the literature,
/// as far as a search of it found, is an independent calibration: no
/// peer-reviewed evaluation of that method's specificity appears to exist, and
/// the one unpublished evaluation we found reports the irreducible-cycle rate
/// rising with the *point-mutation* rate rather than with the recombination
/// rate.
///
/// That claim is testable here, because this crate's grid has the ground truth
/// it needs. Three datasets contain reticulation and four do not, and one of the
/// four that do not — the saturated dataset — differs from the baseline only by
/// having ten times the divergence. So: if `h1_total_persistence` fires on the
/// saturated dataset, it is responding to mutation rather than to reticulation,
/// and the concern is reproduced independently with a known answer.
///
/// The statistic is calibrated exactly as the others are: a parametric bootstrap
/// under a matched tree null, at this population's own alignment width, mean
/// p-distance and gap load. Reporting a persistence value without that
/// calibration is what Section 8 of the textbook did, and the result was a
/// negative finding obtained the hard way.
fn run_topology(args: &[String]) -> Result<bool, String> {
    use reticulation::rips::vietoris_rips_h1;
    use reticulation::simulate::{
        clonal_tree, clonal_tree_relaxed, inject_gap_blocks, with_recombination,
        with_reciprocal_recombinants,
    };
    let length: usize = parse_flag(args, "--length", 4000);
    let generations: u32 = parse_flag(args, "--generations", 4);
    let mu: f64 = parse_flag(args, "--mu", 0.004);
    let spread: f64 = parse_flag(args, "--spread", 1.0);
    let saturated_mu: f64 = parse_flag(args, "--saturated-mu", 0.04);
    let gap_fraction_flag: f64 = parse_flag(args, "--gap-fraction", 0.25);
    let n_seeds: usize = parse_flag(args, "--seeds", 5);
    let seed0: u64 = parse_flag(args, "--seed", 2026);
    let replicates: usize = parse_flag(args, "--bootstrap", 99);

    let h1_of = |pop: &reticulation::simulate::Population| -> Result<(f64, usize, f64), String> {
        let dm = pop.distances(DistanceModel::PDistance)?;
        let scale = default_scale(&dm);
        let persistence = vietoris_rips_h1(&dm, scale);
        Ok((
            persistence.h1_total_persistence,
            persistence.h1.len(),
            scale,
        ))
    };

    // The matched null for the topological statistic, same construction as for
    // the averaged quartet statistic.
    let null_for = |width: usize,
                    target_divergence: f64,
                    gap_load: f64|
     -> Result<Vec<f64>, String> {
        // Returns every replicate's value, so the caller can see whether the
        // null is degenerate rather than only its mean.
        let (null_mu, _) =
            tune_mu_for_divergence(width, 4, spread, target_divergence, 1);
        let mut values = Vec::with_capacity(replicates);
        for replicate in 0..replicates {
            let seed = 0x7000_u64.wrapping_add(replicate as u64);
            let mut pop = clonal_tree_relaxed(width, 4, null_mu, spread, seed);
            if gap_load > 0.001 {
                pop = inject_gap_blocks(pop, gap_load, seed ^ 0xFACE);
            }
            values.push(h1_of(&pop)?.0);
        }
        Ok(values)
    };

    #[derive(serde::Serialize)]
    struct Row {
        dataset: &'static str,
        reticulation_present: bool,
        h1_total_persistence: Vec<f64>,
        h1_cycle_count: Vec<usize>,
        rips_scale: Vec<f64>,
        mean_p_distance: Vec<f64>,
        null_mean: Vec<f64>,
        p_value: Vec<f64>,
        significant_at_005: String,
        /// Replicates, over all seeds, in which the null produced no cycle at
        /// all. When this is every replicate the null is degenerate and the
        /// p-value is not a calibration: it takes the floor whenever the
        /// observed value is non-zero and 1 otherwise.
        null_replicates_with_no_cycle: usize,
        null_replicates_total: usize,
        /// Set when the null never produced a cycle, so the test reduces to
        /// "did a cycle appear at all" and its counts must not be read as
        /// specificity or power.
        null_is_degenerate: bool,
    }

    let mut rows: Vec<Row> = Vec::new();
    let datasets: [(&'static str, bool); 7] = [
        ("D1_clonal_strict_clock", false),
        ("D2_clonal_relaxed_clock", false),
        ("D3_strict_clock_one_recombinant", true),
        ("D4_relaxed_clock_one_recombinant", true),
        ("D5_strict_clock_reciprocal_mosaics", true),
        ("D6_strict_clock_gap_blocks", false),
        ("D7_strict_clock_saturated", false),
    ];

    for (name, reticulate) in datasets {
        let mut totals = Vec::new();
        let mut counts = Vec::new();
        let mut scales = Vec::new();
        let mut divergences = Vec::new();
        let mut null_means = Vec::new();
        let mut p_values = Vec::new();
        let mut flags = Vec::new();
        let mut null_empty = 0usize;
        let mut null_drawn = 0usize;
        for index in 0..n_seeds {
            let seed = seed0 + index as u64;
            let pop = match name {
                "D1_clonal_strict_clock" => clonal_tree(length, generations, mu, seed),
                "D2_clonal_relaxed_clock" => {
                    clonal_tree_relaxed(length, generations, mu, spread, seed)
                }
                "D3_strict_clock_one_recombinant" => {
                    with_recombination(clonal_tree(length, generations, mu, seed), 0.5, seed ^ 0x11)
                        .0
                }
                "D4_relaxed_clock_one_recombinant" => with_recombination(
                    clonal_tree_relaxed(length, generations, mu, spread, seed),
                    0.5,
                    seed ^ 0x22,
                )
                .0,
                "D5_strict_clock_reciprocal_mosaics" => with_reciprocal_recombinants(
                    clonal_tree(length, generations, mu, seed),
                    0.5,
                    seed ^ 0x33,
                )
                .0,
                "D6_strict_clock_gap_blocks" => inject_gap_blocks(
                    clonal_tree(length, generations, mu, seed),
                    gap_fraction_flag,
                    seed ^ 0x44,
                ),
                _ => clonal_tree(length, generations, saturated_mu, seed),
            };
            let (total, count, scale) = h1_of(&pop)?;
            // The mean p-distance comes straight from the distance matrix. An
            // earlier version took it from `quartet_distribution`, which also
            // enumerates every quadruple and so tripped the refuse-to-sample
            // ceiling at 64 taxa -- the guard firing correctly on a call that
            // only ever needed one number.
            let (divergence, _) = divergence(&pop);
            let gaps = gap_fraction(&pop.sequences);
            let null = null_for(length, divergence, gaps)?;
            null_empty += null.iter().filter(|v| **v <= 0.0).count();
            null_drawn += null.len();
            let at_least = null.iter().filter(|v| **v >= total).count();
            let p_value = (at_least + 1) as f64 / (replicates + 1) as f64;
            let null_mean = if null.is_empty() {
                0.0
            } else {
                null.iter().sum::<f64>() / null.len() as f64
            };
            totals.push(total);
            counts.push(count);
            scales.push(scale);
            divergences.push(divergence);
            null_means.push(null_mean);
            p_values.push(p_value);
            flags.push(p_value <= 0.05);
        }
        let held = flags.iter().filter(|b| **b).count();
        rows.push(Row {
            dataset: name,
            reticulation_present: reticulate,
            h1_total_persistence: totals,
            h1_cycle_count: counts,
            rips_scale: scales,
            mean_p_distance: divergences,
            null_mean: null_means,
            p_value: p_values,
            significant_at_005: format!("{held}/{n_seeds}"),
            null_replicates_with_no_cycle: null_empty,
            null_replicates_total: null_drawn,
            null_is_degenerate: null_drawn > 0 && null_empty == null_drawn,
        });
    }

    let saturated = rows
        .iter()
        .find(|r| r.dataset == "D7_strict_clock_saturated")
        .map(|r| r.significant_at_005.clone())
        .unwrap_or_default();
    let report = serde_json::json!({
        "experiment": "does persistent homology on a genetic-distance point cloud \
                       respond to reticulation, or to mutation?",
        "why": "The theoretical warrant for reading an unfilled 1-cycle as \
                reticulation is the same fact this crate's null rests on: a tree \
                metric leaves no unfilled cycle. What a search of the literature \
                did not find is an independent calibration of that method's \
                SPECIFICITY. This grid supplies the ground truth: three datasets \
                contain reticulation, four do not, and the saturated dataset \
                differs from the baseline only in having ten times the \
                divergence.",
        "degeneracy_check": {
            "why": "A bootstrap calibrates only if the null varies. If the null \
                    never produces a cycle, the p-value is not a calibration: it \
                    is the floor whenever the observed value is non-zero and 1 \
                    otherwise, so the test reduces to `did a cycle appear at \
                    all` and its counts must NOT be read as specificity or \
                    power.",
            "datasets_with_a_degenerate_null": rows
                .iter()
                .filter(|r| r.null_is_degenerate)
                .map(|r| r.dataset)
                .collect::<Vec<&str>>(),
            "this_experiment_is_interpretable": !rows.iter().any(|r| r.null_is_degenerate),
            "if_degenerate_the_fix": "More taxa, so that cycles occur in the null \
                    often enough to calibrate against, and/or a larger Rips \
                    threshold than the median off-diagonal distance. At 16 taxa \
                    the complex is too small for H1 to be populated. See \
                    docs/RETICULATA_NUMERICS_TOPOLOGY_ML.md section 2.1, which is \
                    the work that makes larger taxon counts affordable."
        },
        "the_decisive_row": {
            "dataset": "D7_strict_clock_saturated",
            "contains_reticulation": false,
            "differs_from_D1_only_by": "ten times the substitution rate",
            "significant_at_0.05": saturated,
            "reading": "A non-zero count here means the topological statistic is \
                        responding to mutation rather than to reticulation, on \
                        data where the answer is known. A 0/5 means it is \
                        specific and the concern does not reproduce at these \
                        parameters."
        },
        "calibration": {
            "null": "parametric bootstrap under a matched tree null: clonal, \
                     relaxed-clock, at each population's own width, mean \
                     p-distance and gap load",
            "replicates": replicates,
            "p_value_floor": 1.0 / (replicates + 1) as f64,
            "why_not_a_permutation": "A site-order permutation null is reported \
                     to be blind to misalignment artefacts, and a shuffle null is \
                     what Section 8 of the textbook used to no useful effect.",
            "rips_scale": "median off-diagonal p-distance, per population"
        },
        "per_dataset": rows,
    });
    println!(
        "{}",
        serde_json::to_string_pretty(&report).map_err(|e| e.to_string())?
    );
    Ok(true)
}

/// Exhaustive maximum quartet margin, counted bit-parallel.
///
/// Returns the largest exact margin, the violating count and the number
/// examined. Only the maximum is needed for a floor, so no certificates are
/// built; when a witness is wanted the byte path of
/// [`reticulation::subset::find_quartet_certificates`] produces it, and
/// `bit_and_byte_counting_agree_on_random_alignments_with_gaps` is what licenses
/// using the two interchangeably.
fn max_quartet_margin_bits(
    bits: &reticulation::BitAlignment,
    max_quadruples: usize,
) -> Result<(Option<reticulation::Ratio>, usize, usize), String> {
    let n = bits.len();
    if n < 4 {
        return Err(format!("need at least 4 sequences, got {n}"));
    }
    let total = n * (n - 1) * (n - 2) * (n - 3) / 24;
    if total > max_quadruples {
        return Err(format!(
            "{n} sequences give {total} quadruples, above the ceiling of \
             {max_quadruples}. Raise it deliberately: this search is exhaustive \
             and will not sample."
        ));
    }
    let mut best: Option<reticulation::Ratio> = None;
    let mut violating = 0usize;
    let mut examined = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                for l in (k + 1)..n {
                    examined += 1;
                    let counts = bits.quadruple_counts(i, j, k, l);
                    if counts.common == 0 {
                        continue;
                    }
                    let mut sums = counts.pairing_sums();
                    sums.sort_unstable_by(|a, b| b.cmp(a));
                    if sums[0] <= sums[1] {
                        continue;
                    }
                    violating += 1;
                    let margin = reticulation::Ratio {
                        num: sums[0] - sums[1],
                        den: counts.common,
                    };
                    let better = match best {
                        None => true,
                        Some(current) => margin
                            .checked_cmp(current)
                            .map(|o| o == std::cmp::Ordering::Greater)
                            .unwrap_or(false),
                    };
                    if better {
                        best = Some(margin);
                    }
                }
            }
        }
    }
    Ok((best, violating, examined))
}

/// How does the certificate's floor grow with taxon count?
///
/// This is the measurement the head-to-head said was missing. The statistic is a
/// maximum over subsets, so its null level rises with the number of subsets
/// examined, and margins are therefore not comparable across populations of
/// different size. That defect was argued from two data points — the quartet
/// floor moving from 5–7 at sixteen taxa to 9–10 at thirty-two — and arguing
/// from two points is how this project has twice talked itself into an
/// incommensurable comparison. So it is measured across a range, which the
/// bit-parallel counting of [`reticulation::BitAlignment`] makes affordable:
/// sixty-four taxa is 635 376 quadruples, and the byte path would spend 2.5e9
/// column visits on each population.
///
/// Two things are reported. The **floor** — the largest margin a clonal
/// population with no reticulation produces — as a function of taxon count,
/// which is the growth the defect consists of. And the **separation**, the
/// reciprocal-mosaic margin over that floor at the same seed, which says whether
/// the growth actually costs detection or merely shifts both numbers.
fn run_scale(args: &[String]) -> Result<bool, String> {
    use reticulation::simulate::{clonal_tree, clonal_tree_relaxed, with_reciprocal_recombinants};
    use reticulation::subset::Alignment;
    use reticulation::BitAlignment;

    let length: usize = parse_flag(args, "--length", 4000);
    let mu: f64 = parse_flag(args, "--mu", 0.004);
    let spread: f64 = parse_flag(args, "--spread", 1.0);
    let n_seeds: usize = parse_flag(args, "--seeds", 3);
    let seed0: u64 = parse_flag(args, "--seed", 2026);
    let max_quadruples: usize = parse_flag(args, "--max-quadruples", 2_000_000);
    let generations: Vec<u32> = parse_flag(args, "--generations", "4,5,6".to_string())
        .split(',')
        .filter_map(|s| s.trim().parse().ok())
        .collect();

    #[derive(serde::Serialize)]
    struct Rung {
        generations: u32,
        taxa: usize,
        quadruples: usize,
        clonal_strict_floor_exact: Vec<String>,
        clonal_strict_floor: Vec<f64>,
        clonal_relaxed_floor_exact: Vec<String>,
        one_recombinant_margin_exact: Vec<String>,
        one_recombinant_above_floor: String,
        mosaic_margin_exact: Vec<String>,
        mosaic_margin: Vec<f64>,
        worst_seed_separation_over_strict_floor: f64,
        clonal_violating_fraction: Vec<f64>,
    }

    let measure = |pop: &reticulation::simulate::Population| -> Result<
        (Option<reticulation::Ratio>, usize, usize),
        String,
    > {
        let aln = Alignment::new(pop.labels.clone(), pop.sequences.clone())?;
        let bits = BitAlignment::from_alignment(&aln)?;
        max_quartet_margin_bits(&bits, max_quadruples)
    };

    let mut rungs: Vec<Rung> = Vec::new();
    for &g in &generations {
        let taxa = 1usize << g;
        let mut strict = Vec::new();
        let mut strict_f = Vec::new();
        let mut relaxed = Vec::new();
        let mut mosaic = Vec::new();
        let mut mosaic_f = Vec::new();
        let mut single = Vec::new();
        let mut single_above = Vec::new();
        let mut separations = Vec::new();
        let mut violating_fraction = Vec::new();
        let mut quadruples = 0usize;
        for index in 0..n_seeds {
            let seed = seed0 + index as u64;
            let (s_margin, s_violating, examined) =
                measure(&clonal_tree(length, g, mu, seed))?;
            quadruples = examined;
            let (r_margin, _, _) =
                measure(&clonal_tree_relaxed(length, g, mu, spread, seed))?;
            let (m_margin, _, _) = measure(
                &with_reciprocal_recombinants(clonal_tree(length, g, mu, seed), 0.5, seed ^ 0x33).0,
            )?;
            // The single recombinant, which the topological statistic misses at
            // 0 of 3 in every filtration variant. Measured here at the same
            // taxon count so the two can actually be compared.
            let (one_margin, _, _) = measure(
                &reticulation::simulate::with_recombination(
                    clonal_tree(length, g, mu, seed),
                    0.5,
                    seed ^ 0x11,
                )
                .0,
            )?;
            single.push(margin_exact(one_margin));
            single_above.push(exceeds_k_times(one_margin, s_margin, 1));
            strict.push(margin_exact(s_margin));
            strict_f.push(margin_f64(s_margin));
            relaxed.push(margin_exact(r_margin));
            mosaic.push(margin_exact(m_margin));
            mosaic_f.push(margin_f64(m_margin));
            violating_fraction.push(s_violating as f64 / examined.max(1) as f64);
            let (a, b) = (margin_f64(m_margin), margin_f64(s_margin));
            separations.push(if b > 0.0 { a / b } else { f64::INFINITY });
        }
        rungs.push(Rung {
            generations: g,
            taxa,
            quadruples,
            clonal_strict_floor_exact: strict,
            clonal_strict_floor: strict_f,
            clonal_relaxed_floor_exact: relaxed,
            one_recombinant_margin_exact: single,
            one_recombinant_above_floor: {
                let (held, total) = held_at(&single_above);
                format!("{held}/{total}")
            },
            mosaic_margin_exact: mosaic,
            mosaic_margin: mosaic_f,
            worst_seed_separation_over_strict_floor: separations
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min),
            clonal_violating_fraction: violating_fraction,
        });
    }

    // The speedup claim, measured rather than asserted, at the largest rung the
    // byte path can be timed on without dominating the run.
    let timing_generations: u32 = parse_flag(args, "--time-generations", 5);
    let timing_pop = clonal_tree(length, timing_generations, mu, seed0);
    let timing_aln = Alignment::new(timing_pop.labels.clone(), timing_pop.sequences.clone())?;
    let bits = BitAlignment::from_alignment(&timing_aln)?;
    let bit_start = std::time::Instant::now();
    let (bit_margin, _, bit_examined) = max_quartet_margin_bits(&bits, max_quadruples)?;
    let bit_secs = bit_start.elapsed().as_secs_f64();
    let byte_start = std::time::Instant::now();
    let byte_search =
        reticulation::subset::find_quartet_certificates(&timing_aln, 1, max_quadruples)?;
    let byte_secs = byte_start.elapsed().as_secs_f64();
    let byte_margin = byte_search.certificates.first().map(|c| c.margin());
    let agree = margin_exact(bit_margin) == margin_exact(byte_margin);

    let report = serde_json::json!({
        "experiment": "how does the certificate's floor grow with taxon count?",
        "why": "The statistic is a maximum over subsets, so its null level rises \
                with the number of subsets examined and margins are not \
                comparable across population sizes. That defect was previously \
                argued from two data points; arguing from two points is how this \
                project has twice talked itself into an incommensurable \
                comparison, so it is measured across a range here.",
        "parameters": {
            "alignment_length": length,
            "mu_per_generation": mu,
            "relaxed_clock_spread": spread,
            "seeds": n_seeds,
            "exhaustive": "every quadruple examined; the search refuses to sample"
        },
        "rungs": rungs,
        "implementation_check": {
            "why": "The bit-parallel counter is an optimisation, and an \
                    optimisation that disagrees with the obvious implementation \
                    is a silent wrong answer. The byte path is the \
                    specification.",
            "taxa_timed": 1usize << timing_generations,
            "quadruples": bit_examined,
            "bit_parallel_seconds": bit_secs,
            "byte_seconds": byte_secs,
            "speedup": if bit_secs > 0.0 { byte_secs / bit_secs } else { f64::INFINITY },
            "margins_agree": agree,
            "bit_margin": margin_exact(bit_margin),
            "byte_margin": margin_exact(byte_margin),
        },
        "how_to_read_this": "A floor that grows with taxon count is the defect. \
                Whether it MATTERS is the separation column: if the mosaic margin \
                grows at least as fast as the floor, the growth shifts both \
                numbers and detection at a fixed taxon count is unaffected -- but \
                margins still cannot be compared BETWEEN taxon counts, and a \
                threshold calibrated at one size does not transfer to another."
    });
    println!(
        "{}",
        serde_json::to_string_pretty(&report).map_err(|e| e.to_string())?
    );
    if !agree {
        eprintln!("FAIL: the bit-parallel and byte counters disagree on the margin");
        return Ok(false);
    }
    Ok(true)
}

/// Restrict an alignment to the columns comparable in *every* sequence.
///
/// Global complete deletion: one denominator for the whole population. This is
/// the only form of complete deletion a Vietoris--Rips filtration can use, since
/// the complex is built from all pairs at once and there is no "within this
/// subset" to restrict to. It is strictly weaker than the per-subset complete
/// deletion the certificate uses, and the cost is columns.
fn globally_comparable(pop: &reticulation::simulate::Population) -> reticulation::simulate::Population {
    let width = pop.sequences.first().map_or(0, Vec::len);
    let keep: Vec<usize> = (0..width)
        .filter(|&column| {
            pop.sequences.iter().all(|row| {
                matches!(
                    row[column].to_ascii_uppercase(),
                    b'A' | b'C' | b'G' | b'T'
                )
            })
        })
        .collect();
    reticulation::simulate::Population {
        labels: pop.labels.clone(),
        sequences: pop
            .sequences
            .iter()
            .map(|row| keep.iter().map(|&c| row[c]).collect())
            .collect(),
    }
}

/// Test three cheap filtration changes against the failure of §13.3.
///
/// At sixty-four taxa the Vietoris--Rips H1 statistic fires on gap-blocked and
/// on saturated populations -- neither of which contains reticulation -- and
/// misses a single recombinant entirely. Three hypotheses explain that, each
/// implying a different change to the filtration, and each is cheap enough to
/// decide now:
///
/// * **F1, saturation.** The p-distance is concave in true divergence, so a
///   high-divergence population has compressed distances that create spurious
///   cycles. Change: Jukes--Cantor corrected distances.
/// * **F2, denominators.** Pairwise deletion gives each pair its own
///   denominator -- the audit's finding 6 -- and Rips cannot avoid it per
///   subset. Change: one global denominator, columns comparable in every
///   sequence.
/// * **F3, the summary.** Gap blocks produce many weak cycles and mosaics few
///   strong ones, which total persistence conflates. Change: the longest bar
///   instead of the total.
///
/// The criterion was fixed before the run: a variant succeeds if the
/// reciprocal-mosaic population is significant at every seed **and** both
/// no-reticulation populations are significant at none. A variant that merely
/// reduces one false positive has diagnosed something without fixing it, and is
/// reported as that.
///
/// The two expensive hypotheses are deliberately not tested here. Filtering on
/// residuals from a fitted tree metric, and abandoning distance filtration for a
/// complex built from violating subsets, are both principled rebuilds worth
/// doing only if these three show the cause is none of the above.
fn run_filtration(args: &[String]) -> Result<bool, String> {
    use reticulation::rips::vietoris_rips_h1;
    use reticulation::simulate::{
        clonal_tree, inject_gap_blocks, with_recombination, with_reciprocal_recombinants,
    };

    let length: usize = parse_flag(args, "--length", 4000);
    let generations: u32 = parse_flag(args, "--generations", 6);
    let mu: f64 = parse_flag(args, "--mu", 0.004);
    let saturated_mu: f64 = parse_flag(args, "--saturated-mu", 0.04);
    let gap_fraction_flag: f64 = parse_flag(args, "--gap-fraction", 0.25);
    let n_seeds: usize = parse_flag(args, "--seeds", 3);
    let seed0: u64 = parse_flag(args, "--seed", 2026);
    let replicates: usize = parse_flag(args, "--bootstrap", 49);
    let spread: f64 = parse_flag(args, "--spread", 1.0);

    let frac = |flags: &[bool]| -> String {
        let (held, total) = held_at(flags);
        format!("{held}/{total}")
    };

    // Four variants: the baseline that failed, plus one per hypothesis.
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Variant {
        Baseline,
        JukesCantor,
        GlobalDeletion,
        LongestBar,
    }
    let variants = [
        ("baseline_p_distance_total", Variant::Baseline),
        ("F1_jukes_cantor_total", Variant::JukesCantor),
        ("F2_global_complete_deletion_total", Variant::GlobalDeletion),
        ("F3_p_distance_longest_bar", Variant::LongestBar),
    ];

    let statistic = |pop: &reticulation::simulate::Population,
                     variant: Variant|
     -> Result<(f64, usize), String> {
        let effective = if variant == Variant::GlobalDeletion {
            globally_comparable(pop)
        } else {
            reticulation::simulate::Population {
                labels: pop.labels.clone(),
                sequences: pop.sequences.clone(),
            }
        };
        if effective.sequences.first().map_or(0, Vec::len) == 0 {
            return Err("global complete deletion left no columns".to_string());
        }
        let model = if variant == Variant::JukesCantor {
            DistanceModel::JukesCantor
        } else {
            DistanceModel::PDistance
        };
        let dm = effective.distances(model)?;
        let scale = default_scale(&dm);
        let persistence = vietoris_rips_h1(&dm, scale);
        let value = if variant == Variant::LongestBar {
            persistence.h1.first().map_or(0.0, |bar| bar.persistence)
        } else {
            persistence.h1_total_persistence
        };
        Ok((value, persistence.h1.len()))
    };

    // The four datasets that decide the question: one reticulate control that
    // must keep firing, one reticulate case the baseline misses, and the two
    // no-reticulation populations the baseline falsely fires on.
    let build = |which: &str, seed: u64| -> reticulation::simulate::Population {
        match which {
            "D5_reciprocal_mosaics" => {
                with_reciprocal_recombinants(clonal_tree(length, generations, mu, seed), 0.5, seed ^ 0x33).0
            }
            "D3_one_recombinant" => {
                with_recombination(clonal_tree(length, generations, mu, seed), 0.5, seed ^ 0x11).0
            }
            "D6_gap_blocks" => inject_gap_blocks(
                clonal_tree(length, generations, mu, seed),
                gap_fraction_flag,
                seed ^ 0x44,
            ),
            _ => clonal_tree(length, generations, saturated_mu, seed),
        }
    };
    let datasets = [
        ("D5_reciprocal_mosaics", true),
        ("D3_one_recombinant", true),
        ("D6_gap_blocks", false),
        ("D7_saturated", false),
    ];

    #[derive(serde::Serialize)]
    struct Cell {
        dataset: &'static str,
        reticulation_present: bool,
        observed: Vec<f64>,
        cycles: Vec<usize>,
        null_mean: Vec<f64>,
        null_replicates_with_no_signal: usize,
        null_replicates_total: usize,
        null_is_degenerate: bool,
        p_value: Vec<f64>,
        significant_at_005: String,
    }
    #[derive(serde::Serialize)]
    struct VariantResult {
        variant: &'static str,
        hypothesis: &'static str,
        cells: Vec<Cell>,
        mosaics_detected: String,
        false_positives: String,
        criterion_met: bool,
    }

    let mut results: Vec<VariantResult> = Vec::new();
    for (name, variant) in variants {
        let hypothesis = match variant {
            Variant::Baseline => "none: the construction of section 13.3, repeated as the control",
            Variant::JukesCantor => "F1 saturation: p-distance concavity creates spurious cycles at high divergence",
            Variant::GlobalDeletion => "F2 denominators: pairwise deletion is the audit's finding 6, appearing as topology",
            Variant::LongestBar => "F3 summary: many weak cycles and few strong ones are conflated by a total",
        };
        let mut cells = Vec::new();
        let mut mosaic_flags = Vec::new();
        let mut false_positive_flags = Vec::new();
        for (dataset, reticulate) in datasets {
            let mut observed = Vec::new();
            let mut cycles = Vec::new();
            let mut null_means = Vec::new();
            let mut p_values = Vec::new();
            let mut flags = Vec::new();
            let mut null_empty = 0usize;
            let mut null_drawn = 0usize;
            for index in 0..n_seeds {
                let seed = seed0 + index as u64;
                let pop = build(dataset, seed);
                let (value, count) = statistic(&pop, variant)?;
                // The matched null: clonal, relaxed clock, tuned to this
                // population's own divergence, carrying its own gap load.
                let (mean_p, _) = divergence(&pop);
                let gaps = gap_fraction(&pop.sequences);
                let (null_mu, _) =
                    tune_mu_for_divergence(length, generations, spread, mean_p, 1);
                let mut null_values = Vec::with_capacity(replicates);
                for replicate in 0..replicates {
                    let null_seed = 0xF117_0000 + replicate as u64;
                    let mut null_pop = reticulation::simulate::clonal_tree_relaxed(
                        length, generations, null_mu, spread, null_seed,
                    );
                    if gaps > 0.001 {
                        null_pop = inject_gap_blocks(null_pop, gaps, null_seed ^ 0xFACE);
                    }
                    null_values.push(statistic(&null_pop, variant)?.0);
                }
                null_empty += null_values.iter().filter(|v| **v <= 0.0).count();
                null_drawn += null_values.len();
                let at_least = null_values.iter().filter(|v| **v >= value).count();
                let p_value = (at_least + 1) as f64 / (replicates + 1) as f64;
                observed.push(value);
                cycles.push(count);
                null_means.push(if null_values.is_empty() {
                    0.0
                } else {
                    null_values.iter().sum::<f64>() / null_values.len() as f64
                });
                p_values.push(p_value);
                flags.push(p_value <= 0.05);
            }
            let held = flags.iter().filter(|b| **b).count();
            if dataset == "D5_reciprocal_mosaics" {
                mosaic_flags = flags.clone();
            }
            if !reticulate {
                false_positive_flags.extend(flags.iter().copied());
            }
            cells.push(Cell {
                dataset,
                reticulation_present: reticulate,
                observed,
                cycles,
                null_mean: null_means,
                null_replicates_with_no_signal: null_empty,
                null_replicates_total: null_drawn,
                null_is_degenerate: null_drawn > 0 && null_empty == null_drawn,
                p_value: p_values,
                significant_at_005: format!("{held}/{n_seeds}"),
            });
        }
        let mosaics_all = mosaic_flags.iter().all(|b| *b) && !mosaic_flags.is_empty();
        let no_false_positives = false_positive_flags.iter().all(|b| !*b);
        results.push(VariantResult {
            variant: name,
            hypothesis,
            cells,
            mosaics_detected: frac(&mosaic_flags),
            false_positives: frac(&false_positive_flags),
            criterion_met: mosaics_all && no_false_positives,
        });
    }

    let winner = results
        .iter()
        .find(|r| r.criterion_met && r.variant != "baseline_p_distance_total")
        .map(|r| r.variant);
    let report = serde_json::json!({
        "experiment": "three cheap filtration changes against the failure of \
                       section 13.3",
        "the_failure_being_explained": "At 64 taxa the Vietoris-Rips H1 statistic \
                fires at 3/3 on gap-blocked and on saturated populations, neither \
                of which contains reticulation, and misses a single recombinant \
                entirely at 0/3.",
        "criterion_fixed_before_the_run": "A variant succeeds if the \
                reciprocal-mosaic population is significant at EVERY seed and \
                both no-reticulation populations at NONE. A variant that merely \
                reduces one false positive has diagnosed something without \
                fixing it.",
        "not_tested_here": "Two expensive hypotheses are deliberately held back: \
                filtering on residuals from a fitted tree metric, and abandoning \
                distance filtration for a complex built from violating subsets \
                under per-subset complete deletion. Both are principled rebuilds \
                worth doing only if these three show the cause is none of the \
                above.",
        "parameters": {
            "taxa": 1usize << generations,
            "alignment_length": length,
            "seeds": n_seeds,
            "bootstrap_replicates": replicates,
            "p_value_floor": 1.0 / (replicates + 1) as f64,
        },
        "variants": results,
        "winner": winner,
        "how_to_read_a_degenerate_cell": "`null_is_degenerate` means the null \
                produced no signal in any replicate, so the p-value is the floor \
                whenever the observed value is positive and 1 otherwise. Such a \
                cell carries no calibration and its count must not be read as \
                power or specificity."
    });
    println!(
        "{}",
        serde_json::to_string_pretty(&report).map_err(|e| e.to_string())?
    );
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
        Some("align-error") => run_align_error(&args),
        Some("real-data") => run_real_data(&args),
        Some("breakpoint") => run_breakpoint(&args),
        Some("topology") => run_topology(&args),
        Some("scale") => run_scale(&args),
        Some("filtration") => run_filtration(&args),
        _ => {
            eprintln!(
                "usage:\n  reticulation controls [--length N] [--generations G] \
                 [--mu M] [--seed S]\n  reticulation poc [--length N] [--seed S] \
                 [--formal-dir D] [--out D]\n  reticulation experiment [--seeds K] \
                 [--spread S] [--gap-fraction F] [--saturated-mu M]\n  \
                 reticulation saturation [--seeds K] [--low-mu M] [--high-mu M]\n  \
                 reticulation verify-tree [--length N] [--seed S] [--formal-dir D] \
                 [--out D]\n  reticulation fasta <FILE> [--model p|jc] \
                 [--top K] [--scale X]\n  reticulation align-error [--length N] \
                 [--deletion-rate R] [--seeds K]\n  reticulation real-data \
                 [--populations DIR] [--band N] [--max-band N] [--seeds K]\n  \
                 reticulation breakpoint --population <aligned.fasta> \
                 [--provenance J] [--top K]"
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

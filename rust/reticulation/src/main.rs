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

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.first().map(String::as_str) {
        Some("controls") => run_controls(&args),
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

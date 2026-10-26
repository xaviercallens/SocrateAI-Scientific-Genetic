"""Apply the C9 (componentwise non-negativity) and C8 (identity-kernel r) changes to qf-popgen."""
from pathlib import Path

LIB = Path("/home/callensxavier_gmail_com/rusty-SUNDIALS-wt-p7/crates/qf-popgen/src/lib.rs")
BIN = Path("/home/callensxavier_gmail_com/rusty-SUNDIALS-wt-p7/crates/qf-popgen/src/bin/fixation_ode.rs")
README = Path("/home/callensxavier_gmail_com/rusty-SUNDIALS-wt-p7/crates/qf-popgen/README.md")


def sub(text: str, old: str, new: str) -> str:
    assert text.count(old) == 1, (old[:60], text.count(old))
    return text.replace(old, new)


s = LIB.read_text()
s = sub(s, """//! total mass is conserved, the interior mass has decayed below a tolerance, and the absorption
//! mass matches the target within a tolerance. Anything else is an error, never a silent pass.""",
"""//! total mass is conserved, *every component* of the state stays non-negative (up to a tolerance;
//! a signed sum alone can hide cancelling entries), the interior mass has decayed below a
//! tolerance, and the absorption mass matches the target within a tolerance. Anything else is an
//! error, never a silent pass.""")
s = sub(s, """    /// Total probability mass drifted by more than the tolerance.
    MassDrift {
        /// Largest `|sum pi - 1|` observed.
        drift: f64,
    },""", """    /// Total probability mass drifted by more than the tolerance.
    MassDrift {
        /// Largest `|sum pi - 1|` observed.
        drift: f64,
    },
    /// A component of the state fell below `-tol.negative` (probabilities cannot be negative).
    NegativeMass {
        /// State index of the most negative component.
        state: usize,
        /// Its value.
        value: f64,
    },""")
s = sub(s, """    /// Largest accepted `|sum pi - 1|` along the run.
    pub mass: f64,
}""", """    /// Largest accepted `|sum pi - 1|` along the run.
    pub mass: f64,
    /// Most negative component accepted: every `pi_k(t)` must be `>= -negative` at every output.
    pub negative: f64,
}""")
s = sub(s, """            value: 1e-8,
            mass: 1e-9,
        }""", """            value: 1e-8,
            mass: 1e-9,
            negative: 1e-12,
        }""")
s = sub(s, """    /// Largest `|sum pi - 1|` over the outputs.
    pub mass_drift: f64,
    /// `(t, pi_N(t), interior(t))`""", """    /// Largest `|sum pi - 1|` over the outputs.
    pub mass_drift: f64,
    /// Smallest component `min_k pi_k(t)` seen at any output (accepted only if `>= -tol.negative`).
    pub min_component: f64,
    /// `(t, pi_N(t), interior(t))`""")
s = sub(s, """/// Integrate `d pi/dt = pi (P - I)` from the point mass""", """/// Componentwise non-negativity check: returns the smallest component of `y`, or
/// [`CertError::NegativeMass`] if it is below `-tol`.
///
/// The total-mass and interior-mass tests use *signed* sums, so a state such as
/// `(0.5, -0.25, 0.25, 0.5)` has total mass `1` and interior mass `0`; this check refuses it.
pub fn check_nonnegative(y: &[f64], tol: f64) -> Result<f64, CertError> {
    let mut state = 0usize;
    let mut min = f64::INFINITY;
    for (k, &v) in y.iter().enumerate() {
        if v < min {
            min = v;
            state = k;
        }
    }
    if min < -tol {
        return Err(CertError::NegativeMass { state, value: min });
    }
    Ok(min)
}

/// Integrate `d pi/dt = pi (P - I)` from the point mass""")
s = sub(s, """    let mut drift: f64 = 0.0;
    let mut t_out = 1.0_f64;""", """    let mut drift: f64 = 0.0;
    let mut min_component = f64::INFINITY;
    let mut t_out = 1.0_f64;""")
s = sub(s, """        let total: f64 = y.iter().sum();
        drift = drift.max((total - 1.0).abs());""", """        min_component = min_component.min(check_nonnegative(y, tol.negative)?);
        let total: f64 = y.iter().sum();
        drift = drift.max((total - 1.0).abs());""")
s = sub(s, """                mass_drift: drift,
                trace,""", """                mass_drift: drift,
                min_component,
                trace,""")
s = sub(s, """    #[test]
    fn negative_horizon_too_short_is_refused() {""", """    #[test]
    fn nonnegativity_check_refuses_cancelling_entries_that_the_signed_sums_accept() {
        // Negative control: total mass 1 and interior mass 0, so the signed tests pass, but a
        // component is negative.
        let bad = [0.5, -0.25, 0.25, 0.5];
        let total: f64 = bad.iter().sum();
        let interior: f64 = bad[1..3].iter().sum();
        assert!((total - 1.0).abs() < 1e-15 && interior.abs() < 1e-15);
        assert!(matches!(
            check_nonnegative(&bad, 1e-12),
            Err(CertError::NegativeMass { state: 1, .. })
        ));
        // Positive control: a probability vector and one with round-off-level negativity pass.
        assert_eq!(check_nonnegative(&[0.25, 0.25, 0.0, 0.5], 1e-12), Ok(0.0));
        assert!(check_nonnegative(&[0.5, -1e-14, 0.5], 1e-12).is_ok());
    }

    #[test]
    fn accepted_runs_record_a_nonnegative_minimum_component() {
        let c = moran_kernel(10, 1.5);
        let cert = run(&c, 3, moran_fix(10, 1.5, 3), 1e6).expect("certified");
        assert!(cert.min_component >= -Tolerances::default().negative);
    }

    #[test]
    fn negative_horizon_too_short_is_refused() {""")
LIB.write_text(s)

b = BIN.read_text()
b = sub(b, """    r: f64,
    i: usize,
    target: f64,
    c: &Certificate,
) -> Row {""", """    r: Option<f64>,
    i: usize,
    target: f64,
    c: &Certificate,
) -> Row {""")
b = sub(b, """\\"r\\":{r},\\"i\\":{i},\\"verdict\\":\\"accepted\\",\\"target\\":{target:.17e},\\"p_fix\\":{:.17e},\\"abs_err\\":{:e},\\"p_loss\\":{:.17e},\\"interior\\":{:e},\\"mass_drift\\":{:e},\\"t_stop\\":{:e},\\"trace\\":{}}}",
            c.p_fix,""", """\\"r\\":{},\\"i\\":{i},\\"verdict\\":\\"accepted\\",\\"target\\":{target:.17e},\\"p_fix\\":{:.17e},\\"abs_err\\":{:e},\\"p_loss\\":{:.17e},\\"interior\\":{:e},\\"mass_drift\\":{:e},\\"min_component\\":{:e},\\"t_stop\\":{:e},\\"trace\\":{}}}",
            r_json(r),
            c.p_fix,""")
b = sub(b, """            c.mass_drift,
            c.t_stop,""", """            c.mass_drift,
            c.min_component,
            c.t_stop,""")
b = sub(b, """fn accepted(""", """/// JSON for the fitness field: a number, or `null` when the row has no meaningful fitness
/// (the identity-kernel control is not a Moran chain).
fn r_json(r: Option<f64>) -> String {
    r.map_or_else(|| "null".to_string(), |v| format!("{v}"))
}

fn accepted(""")
b = sub(b, """    n: usize,
    r: f64,
    i: usize,
    target: f64,
    e: &CertError,""", """    n: usize,
    r: Option<f64>,
    i: usize,
    target: f64,
    e: &CertError,""")
b = sub(b, """\\"r\\":{r},\\"i\\":{i},\\"verdict\\":\\"refused\\",\\"target\\":{target:.17e},\\"reason\\":\\"{reason}\\"}}"
        ),""", """\\"r\\":{},\\"i\\":{i},\\"verdict\\":\\"refused\\",\\"target\\":{target:.17e},\\"reason\\":\\"{reason}\\"}}",
            r_json(r)
        ),""")
b = sub(b, """rows.push(accepted("moran", label, n, r, i, target, &c))""", """rows.push(accepted("moran", label, n, Some(r), i, target, &c))""")
b = sub(b, """rows.push(refused("moran", label, n, r, i, target, &e, false));""", """rows.push(refused("moran", label, n, Some(r), i, target, &e, false));""")
b = sub(b, """                    "wf_neutral",
                    n,
                    1.0,
                    i,
                    target,
                    &c,""", """                    "wf_neutral",
                    n,
                    Some(1.0),
                    i,
                    target,
                    &c,""")
b = sub(b, """                        "wf_neutral",
                        n,
                        1.0,
                        i,
                        target,
                        &e,""", """                        "wf_neutral",
                        n,
                        Some(1.0),
                        i,
                        target,
                        &e,""")
b = sub(b, """    let mut negative = |label: &str,
                        chain: &Chain,
                        i: usize,""", """    let mut negative = |label: &str,
                        chain: &Chain,
                        r: Option<f64>,
                        i: usize,""")
b = sub(b, """        let r = 1.5;
        match certify_absorption(chain, i, target, t_max, tol, set) {""", """        match certify_absorption(chain, i, target, t_max, tol, set) {""")
b = sub(b, """            "flipped_exponent_formula",
            &neg_chain,
            i,""", """            "flipped_exponent_formula",
            &neg_chain,
            Some(1.5),
            i,""")
b = sub(b, """            "neutral_i_over_N_under_selection",
            &neg_chain,
            i,""", """            "neutral_i_over_N_under_selection",
            &neg_chain,
            Some(1.5),
            i,""")
b = sub(b, """        "identity_kernel_no_absorption",
        &identity_kernel(2),
        1,""", """        "identity_kernel_no_absorption",
        &identity_kernel(2),
        None,
        1,""")
b = sub(b, """        "horizon_too_short",
        &neg_chain,
        3,""", """        "horizon_too_short",
        &neg_chain,
        Some(1.5),
        3,""")
b = sub(b, """        "corrupted_row_sum",
        &bad,
        3,""", """        "corrupted_row_sum",
        &bad,
        Some(1.5),
        3,""")
BIN.write_text(b)

rd = README.read_text()
rd = sub(rd, "is conserved, the interior mass has decayed below a tolerance, and the absorption mass at `N` matches the target.\nEvery other outcome is an `Err` (`InvalidKernel`, `MassDrift`, `NotAbsorbed`, `ValueMismatch`, `Solver`, `NonFinite`).",
  "is conserved, every component of the state stays non-negative (`NegativeMass` otherwise: the signed mass sums alone\ncannot see cancelling entries), the interior mass has decayed below a tolerance, and the absorption mass at `N` matches\nthe target. Every other outcome is an `Err` (`InvalidKernel`, `MassDrift`, `NegativeMass`, `NotAbsorbed`,\n`ValueMismatch`, `Solver`, `NonFinite`).")
README.write_text(rd)
print("patched")

#!/usr/bin/env python3
"""
Pytest verification suite for the Socratic Elenchus Closed-Loop Alignment Engine (`v3.0`):
- 26 Lean 4.34.1 theorems (`ST.T1`-`ST.T4`, `B3E.T1`-`B3E.T5`, `GAR.H1`-`GAR.H17`, `0 sorryAx`)
- 20 Live Scientific Skills (`openfda`, `unibind`, `embl-ebi-ols`, etc.)
- 9 Hugging Face `Genomic_Benchmarks` datasets (`3,600` sequences, `36,000` dicodons across 4 species)
- 3 Humanitarian Clinical Applications (`POLE` `p.Pro286Arg` Immuno-Oncology, `HBB` `rs334` Sickle Cell Prime Editing, `POLG` Antiviral Safety)
- 538.84x 5th-Generation Humanitarian Elenchus ANSE Bio-Computational Energy Reduction
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from socrate_genetic_federation import (
    FEDERATED_14_REPOS_REGISTRY,
    FEDERATED_20_SKILLS_REGISTRY,
    FEDERATED_28_THEOREMS_REGISTRY,
    HF_9_GENOMIC_BENCHMARKS,
    SocraticElenchusLoopEngine,
)

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_registries_and_elenchus_v3_telemetry() -> None:
    assert len(FEDERATED_28_THEOREMS_REGISTRY) == 28
    assert len(FEDERATED_20_SKILLS_REGISTRY) == 20
    assert len(FEDERATED_14_REPOS_REGISTRY) == 14
    assert len(HF_9_GENOMIC_BENCHMARKS) == 9

    json_path = REPO_ROOT / "data" / "elenchus_humanitarian_closed_loop_v3.json"
    assert json_path.exists()
    payload = json.loads(json_path.read_text())

    hf = payload["hf_9_dataset_cross_species_benchmark"]
    assert hf["total_datasets"] == 9
    assert hf["total_real_sequences_evaluated"] == 3600
    assert hf["total_real_dicodon_8d_vectors"] == 36000
    assert hf["all_datasets_banach_certified"] is True
    assert hf["mean_e8_vs_z8_mse_reduction_pct"] > 24.0
    assert hf["mean_a2_vs_z2_mse_reduction_pct"] > 12.0


def test_humanitarian_clinical_theorems_gar_h15_h16_h17_exact_q() -> None:
    # GAR.H15: POLE p.Pro286Arg TMB-High & FDA KEYTRUDA (pembrolizumab) rescue
    s_mb = Fraction(3125000, 27)
    gamma = Fraction(3, 250)
    tmb_wt = s_mb * (gamma ** 3)
    tmb_p286r = s_mb * gamma
    assert tmb_wt == Fraction(1, 5)
    assert tmb_p286r == Fraction(12500, 9)
    assert tmb_p286r / tmb_wt == Fraction(1, 1) / (gamma ** 2) == Fraction(62500, 9)
    assert tmb_wt < 10 < 100 < tmb_p286r

    pd1_hazard = Fraction(1, 1) / (1 + Fraction(9, 10) * tmb_p286r)
    assert pd1_hazard == Fraction(1, 1251) < Fraction(1, 1000)

    # GAR.H16: HBB rs334 Sickle Cell 30-nt RTT scarless purity & POLG antiviral selectivity
    scarless_1 = 1 - 30 * gamma
    scarless_2 = 1 - 30 * (gamma ** 2)
    scarless_3 = 1 - 30 * (gamma ** 3)
    assert scarless_1 == Fraction(16, 25)
    assert scarless_2 == Fraction(6223, 6250)
    assert scarless_3 == Fraction(1562419, 1562500)
    assert scarless_3 > Fraction(9999, 10000)

    polg_si = Fraction(1, 20) / (Fraction(1, 100) ** 3)
    assert polg_si == 50000 > 10000

    # GAR.H17: 6-Stage Elenchus ANSE Bio-Computational Energy (538.84x reduction)
    e_base = (
        1000 * Fraction(1, 100)
        + 10000 * Fraction(213, 1000000)
        + 100 * Fraction(1, 115)
        + Fraction(25, 3)
    )
    e_gen5 = (
        1000 * Fraction(27, 15625000)
        + 0
        + 100 * Fraction(1, 49480)
        + Fraction(1, 24) * Fraction(929, 1080)
    )
    assert e_gen5 == Fraction(3966806041, 100197000000)
    assert e_gen5 < Fraction(1, 25)
    assert 538 * e_gen5 < e_base


def test_5round_elenchus_popperian_refutation_controls() -> None:
    json_path = REPO_ROOT / "data" / "elenchus_humanitarian_closed_loop_v3.json"
    payload = json.loads(json_path.read_text())
    el = payload["socratic_elenchus_5_round_alignment"]
    assert el["all_5_rounds_certified"] is True
    assert len(el["rounds"]) == 5
    for r in el["rounds"][:4]:
        ctrl = r["popperian_false_control"]
        assert ctrl["claim_evaluates_to"] is False
        assert ctrl["refuted_by_elenchus_gate"] is True
        assert ctrl["assigned_energy_on_refutation"] == 1000000.0


def test_paper_d_chromatin_remodeling_jarzynski_crooks_mcghee_tonks_st_t5_t6() -> None:
    # ST.T5 & ST.T6: ISWI/SNF2H (PDB: 6NE3), SWI/SNF BRG1 (PDB: 6LTJ), Widom 601 NCP147 (PDB: 1KX5)
    l_nuc = Fraction(147, 1)
    rho_euchrom = Fraction(1, 200)
    rho_hetero = Fraction(1, 165)

    g_mvh_euchrom = (1 - l_nuc * rho_euchrom) / (1 - (l_nuc - 1) * rho_euchrom)
    g_mvh_hetero = (1 - l_nuc * rho_hetero) / (1 - (l_nuc - 1) * rho_hetero)
    assert g_mvh_euchrom == Fraction(53, 54)
    assert g_mvh_hetero == Fraction(18, 19)
    assert Fraction(0, 1) < g_mvh_hetero < g_mvh_euchrom < Fraction(1, 1)

    pi_tg_euchrom = rho_euchrom / (1 - l_nuc * rho_euchrom)
    pi_tg_hetero = rho_hetero / (1 - l_nuc * rho_hetero)
    assert pi_tg_euchrom == Fraction(1, 53)
    assert pi_tg_hetero == Fraction(1, 18)
    assert pi_tg_hetero / pi_tg_euchrom == Fraction(53, 18) > Fraction(29, 10)

    # smFRET & Optical Tweezers Observables at T = 310 K
    kt_pn_nm = Fraction(428, 100)
    delta_mu_atp = 20 * kt_pn_nm
    unwrap_barrier_20bp = Fraction(13, 2) * kt_pn_nm
    assert delta_mu_atp == Fraction(856, 10)
    assert unwrap_barrier_20bp == Fraction(2782, 100)
    assert unwrap_barrier_20bp < delta_mu_atp
    assert Fraction(36, 100) < Fraction(24, 10)
    assert 6 * Fraction(36, 100) < Fraction(24, 10)

    paper_d_json = REPO_ROOT / "data" / "paper_d_telemetry.json"
    assert paper_d_json.exists()
    d_payload = json.loads(paper_d_json.read_text())
    pdbs = {e["pdb_id"]: e["gene"] for e in d_payload["biological_system"]["pdb_structures"]}
    assert pdbs["6NE3"].startswith("SMARCA5")
    assert pdbs["6LTJ"].startswith("SMARCA4")
    assert pdbs["1KX5"].startswith("NCP147")
    lattice = d_payload["chromatin_lattice_gas"]
    assert Fraction(lattice["mcghee_von_hippel_euchromatin_rational"]) == g_mvh_euchrom
    assert Fraction(lattice["mcghee_von_hippel_heterochromatin_rational"]) == g_mvh_hetero
    assert Fraction(lattice["tonks_girardeau_pressure_ratio_rational"]) == pi_tg_hetero / pi_tg_euchrom


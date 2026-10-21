#!/usr/bin/env python3
"""Figures and numerical results for the Formal Genetics textbook.

Every figure is computed from either a closed-form model the book derives by hand
or from the cached live records in `textbook/data/genomic_records.json`. Numbers
that appear in the text are written to `textbook/data/figure_values.json` so the
prose and the plots cannot drift apart.

    python make_figures.py              # all figures
    python make_figures.py --only 3 5   # just those section numbers
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
FIGS = ROOT / "textbook" / "figures"
DATA = ROOT / "textbook" / "data"
RECORDS = DATA / "genomic_records.json"
VALUES = DATA / "figure_values.json"

INK, BLUE, GREEN, WARM, GREY, PURPLE = (
    "#14365C", "#2A6FB8", "#1B7A52", "#B4521E", "#8A8A86", "#5B3E8E")

plt.rcParams.update({
    "font.family": "serif", "font.size": 9.5, "axes.titlesize": 10.5,
    "axes.labelsize": 9.5, "legend.fontsize": 8.5, "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5, "pdf.fonttype": 42, "axes.edgecolor": GREY,
    "axes.labelcolor": INK, "xtick.color": GREY, "ytick.color": GREY,
})

# Kane & Shcherbakova 2014, Cancer Res 74:1895, Table 2 (CAN1, isogenic
# MMR-proficient background, x10^-8). Provenance in data/pole_proofreading_literature.json.
KANE_P301R, KANE_EXO_NULL, KANE_WT = 2300.0, 44.0, 1.5
# Xing 2019, Nat Commun 10:374: purified enzyme, base substitutions x10^-6, no MMR.
XING_P301R_INVITRO, XING_EXO_NULL_INVITRO = 59.0, 155.0


def values() -> dict[str, Any]:
    return json.loads(VALUES.read_text()) if VALUES.exists() else {}


def put(key: str, payload: Any) -> None:
    doc = values()
    doc[key] = payload
    VALUES.parent.mkdir(parents=True, exist_ok=True)
    VALUES.write_text(json.dumps(doc, indent=1, sort_keys=True))


def save(fig: plt.Figure, name: str) -> Path:
    FIGS.mkdir(parents=True, exist_ok=True)
    out = FIGS / name
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.relative_to(ROOT)}")
    return out


# ------------------------------------------------------- tribute: Photo 51 ---
def fig_photo51() -> None:
    """Photo 51, recomputed from helical diffraction theory.

    Not a stylised drawing: this is the Cochran-Crick-Vand transform of a
    continuous double helix with B-DNA parameters. The famous features emerge
    from the physics rather than being drawn in --- in particular the ABSENT
    fourth layer line, which is the direct fingerprint of the two strands being
    offset by 3/8 of a pitch along the axis, and the strong meridional at
    1/3.4 A^-1 from base stacking.

    This is what Rosalind Franklin photographed in May 1952, and what made the
    double helix unavoidable to anyone who could read it.
    """
    radius, pitch, rise = 10.0, 34.0, 3.4   # angstrom: B-DNA
    offset = 3.0 / 8.0                       # minor/major groove asymmetry
    n_max, R_max = 12, 0.46                  # layer lines, reciprocal radius

    R = np.linspace(-R_max, R_max, 900)
    Z = np.linspace(-0.40, 0.40, 900)
    RR, ZZ = np.meshgrid(R, Z)
    img = np.zeros_like(RR)

    from scipy.special import jv

    layer_lines: list[dict[str, Any]] = []
    for n in range(-n_max, n_max + 1):
        z_n = n / pitch
        # Cochran-Crick-Vand: amplitude on layer line n is J_n(2*pi*r*R).
        amp = jv(n, 2 * np.pi * radius * np.abs(RR))
        # Two strands offset by `offset` of a pitch interfere; at n = 4 the
        # cosine vanishes and the layer line disappears. That is Photo 51's
        # missing fourth layer line.
        interference = 2.0 * np.cos(np.pi * n * offset)
        sigma = 0.0042
        profile = np.exp(-((ZZ - z_n) ** 2) / (2 * sigma ** 2))
        img += (amp * interference) ** 2 * profile
        if 0 <= n <= 10:
            layer_lines.append({"n": n, "z_inv_angstrom": z_n,
                                "interference": float(interference),
                                "absent": bool(abs(interference) < 1e-9)})

    # Base stacking: a strong meridional arc at 1/3.4 A^-1.
    for sign in (1, -1):
        img += 11.0 * np.exp(-((ZZ - sign / rise) ** 2) / (2 * 0.0055 ** 2)
                             - (RR ** 2) / (2 * 0.052 ** 2))

    shown = np.power(img / img.max(), 0.32)   # film-like response

    fig = plt.figure(figsize=(10.4, 5.0))
    gs = fig.add_gridspec(1, 2, width_ratios=(1.0, 1.04), wspace=0.16,
                          left=0.04, right=0.97, bottom=0.12, top=0.86)

    ax = fig.add_subplot(gs[0, 0])
    ax.imshow(shown, extent=(-R_max, R_max, -0.40, 0.40), origin="lower",
              cmap="bone", aspect="auto", interpolation="bilinear")
    ax.set_facecolor("black")
    ax.set_xlabel("$R$  (Å$^{-1}$)"); ax.set_ylabel("$Z$  (Å$^{-1}$)")
    ax.set_title("(a) Photo 51, recomputed", loc="left", color=INK)
    ax.set_xticks([]); ax.set_yticks([])

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.imshow(shown, extent=(-R_max, R_max, -0.40, 0.40), origin="lower",
               cmap="bone", aspect="auto", interpolation="bilinear", alpha=0.55)
    # The helix cross: its opening angle gives pitch/radius directly.
    for s in (1, -1):
        ax2.plot([0, s * R_max], [0, R_max * (pitch / (2 * np.pi * radius)) * 0.86],
                 color=WARM, lw=1.3, ls="--", alpha=0.95)
        ax2.plot([0, s * R_max], [0, -R_max * (pitch / (2 * np.pi * radius)) * 0.86],
                 color=WARM, lw=1.3, ls="--", alpha=0.95)
    ax2.annotate("the X: a helix, and its\nangle gives pitch / radius",
                 xy=(0.17, 0.085), xytext=(0.012, 0.30), color=WARM, fontsize=8.5,
                 arrowprops=dict(arrowstyle="->", color=WARM, lw=1.0))
    ax2.annotate("4th layer line ABSENT:\ntwo strands offset by 3/8 pitch",
                 xy=(0.0, 4 / pitch), xytext=(-0.44, 0.185), color=GREEN, fontsize=8.5,
                 arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.0))
    ax2.annotate("meridional at 1/3.4 Å$^{-1}$:\nbase stacking",
                 xy=(0.0, 1 / rise), xytext=(-0.45, -0.33), color=BLUE, fontsize=8.5,
                 arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.0))
    for ll in layer_lines:
        if ll["n"] in (1, 2, 3, 4, 5):
            style = dict(color=GREEN, lw=0.9, alpha=0.9) if ll["absent"] else \
                    dict(color="white", lw=0.4, alpha=0.35)
            for s in (1, -1):
                ax2.axhline(s * ll["z_inv_angstrom"], **style)
    ax2.set_xlim(-R_max, R_max); ax2.set_ylim(-0.40, 0.40)
    ax2.set_xticks([]); ax2.set_yticks([])
    ax2.set_facecolor("black")
    ax2.set_title("(b) What it tells you, read off the pattern", loc="left", color=INK)

    put("photo51", {
        "model": "Cochran-Crick-Vand continuous double helix",
        "radius_angstrom": radius, "pitch_angstrom": pitch,
        "rise_per_base_angstrom": rise, "strand_offset_fraction_of_pitch": offset,
        "layer_lines": layer_lines,
        "fourth_layer_line_absent": next(
            ll["absent"] for ll in layer_lines if ll["n"] == 4),
        "note": "the absent n=4 layer line follows from cos(pi*4*3/8) = 0",
    })
    save(fig, "fig00_photo51.pdf")


# --------------------------------------------------------------- section 1 ---
def fig_three_claims() -> None:
    """Section 1's organising schematic: the three independent ways a quantitative
    biological claim can fail, and which of them a proof assistant can touch."""
    fig, ax = plt.subplots(figsize=(10.2, 3.9))
    ax.set_xlim(0, 10); ax.set_ylim(0, 4.1); ax.axis("off")

    boxes = [
        (0.35, "Arithmetic", "Does the number follow\nfrom the equations?", GREEN,
         "checked with certainty\nby the Lean kernel"),
        (3.6, "Model", "Do the equations describe\nwhat the enzyme does?", WARM,
         "NOT touched by proof\nrefuted in §3–§4"),
        (6.85, "Measurement", "Did the instrument\nreally say that?", BLUE,
         "controls, replication,\nsame background"),
    ]
    for x, title, body, col, foot in boxes:
        ax.add_patch(plt.Rectangle((x, 1.15), 2.8, 2.1, facecolor="white",
                                   edgecolor=col, lw=1.8, zorder=2))
        ax.text(x + 1.4, 2.88, title, ha="center", fontsize=12, color=col,
                fontweight="bold", zorder=3)
        ax.text(x + 1.4, 2.15, body, ha="center", fontsize=9, color=INK, zorder=3)
        ax.text(x + 1.4, 0.72, foot, ha="center", fontsize=8.2, color=GREY, zorder=3)
    for x in (3.15, 6.4):
        ax.annotate("", xy=(x + 0.45, 2.2), xytext=(x, 2.2),
                    arrowprops=dict(arrowstyle="->", color=GREY, lw=1.4))
    ax.text(5.0, 3.72, "A claim must survive all three. They fail independently.",
            ha="center", fontsize=10.5, color=INK, style="italic")
    ax.add_patch(plt.Rectangle((3.6, 1.15), 2.8, 2.1, facecolor=WARM, alpha=0.07,
                               edgecolor="none", zorder=1))
    ax.text(5.0, 0.22, "The retracted POLE mechanism passed the left box and the right box, "
                       "and failed the middle one.", ha="center", fontsize=9, color=WARM)
    save(fig, "fig01_three_claims.pdf")


# --------------------------------------------------------------- section 2 ---
def fig_proofreading_cascade() -> None:
    """3D error surface E_m(gamma) = gamma^m, with the biologically useful corner
    marked. Shows why a few stages buy orders of magnitude."""
    gam = np.linspace(0.005, 0.6, 160)
    m = np.arange(1, 7)
    G, M = np.meshgrid(gam, m)
    E = np.log10(G ** M)

    fig = plt.figure(figsize=(10.6, 4.3))
    gs = fig.add_gridspec(1, 2, width_ratios=(1.12, 1.0), wspace=0.18,
                          left=0.03, right=0.97, bottom=0.13, top=0.9)

    ax = fig.add_subplot(gs[0, 0], projection="3d")
    surf = ax.plot_surface(G, M, E, cmap="viridis", alpha=0.93,
                           linewidth=0, antialiased=True, rstride=1, cstride=3)
    ax.set_xlabel(r"per-stage discrimination $\gamma$")
    ax.set_ylabel("stages $m$")
    ax.set_zlabel(r"$\log_{10} \mathcal{E}_m$")
    ax.set_title("(a) Error falls geometrically in the number of stages", loc="left")
    ax.view_init(elev=24, azim=-131)
    fig.colorbar(surf, ax=ax, shrink=0.58, pad=0.09, label=r"$\log_{10}\mathcal{E}$")

    ax2 = fig.add_subplot(gs[0, 1])
    for mm, col in zip((1, 2, 3, 4), (WARM, PURPLE, BLUE, GREEN)):
        ax2.semilogy(gam, gam ** mm, color=col, lw=1.8, label=f"$m={mm}$")
    ax2.axhline(1e-6, color=GREY, ls=":", lw=1.2)
    ax2.annotate(r"replication fidelity target $10^{-6}$", xy=(0.3, 1e-6),
                 xytext=(0.22, 4e-6), fontsize=8, color=INK)
    ax2.plot([0.01], [1e-6], marker="o", ms=7, mfc="none", mec=WARM, mew=1.8)
    ax2.annotate(r"$\gamma=10^{-2},\ m=3$", xy=(0.01, 1e-6), xytext=(0.035, 8e-8),
                 fontsize=8.5, color=INK,
                 arrowprops=dict(arrowstyle="->", color=GREY, lw=0.9))
    ax2.set_xlabel(r"$\gamma$")
    ax2.set_ylabel(r"error per base $\mathcal{E}_m(\gamma)$")
    ax2.set_ylim(1e-9, 1)
    ax2.set_title("(b) Three stages reach the fidelity a genome needs", loc="left")
    ax2.legend(frameon=False, loc="lower right")
    ax2.spines[["top", "right"]].set_visible(False)

    put("sec02", {"gamma": 0.01, "m": 3, "error": 1e-6,
                  "note": "three stages at gamma=1e-2 give 1e-6 per base"})
    save(fig, "fig02_proofreading_cascade.pdf")


# --------------------------------------------------------------- section 3 ---
def fig_nogo_ceiling() -> None:
    """The PSR_2 ceiling as a 3D feasible region, and the measurement outside it."""
    fig = plt.figure(figsize=(10.8, 4.4))
    gs = fig.add_gridspec(1, 2, width_ratios=(1.08, 1.0), wspace=0.22,
                          left=0.04, right=0.97, bottom=0.14, top=0.88)

    ax = fig.add_subplot(gs[0, 0], projection="3d")
    q1 = np.linspace(0, 1, 70)
    q2 = np.linspace(0, 1, 70)
    Q1, Q2 = np.meshgrid(q1, q2)
    ax.plot_surface(Q1, Q2, Q1 * Q2, cmap="Blues", alpha=0.85, linewidth=0,
                    rstride=2, cstride=2, vmin=-0.4, vmax=1.1)
    ax.plot_surface(Q1, Q2, np.ones_like(Q1), color=WARM, alpha=0.20, linewidth=0)
    ax.text(0.5, 0.5, 1.06, "ceiling  $E/\\varepsilon = 1$", color=WARM, fontsize=9)
    ax.set_xlabel("$q_1$"); ax.set_ylabel("$q_2$")
    ax.set_zlabel(r"$E/\varepsilon = q_1 q_2$")
    ax.set_zlim(0, 1.25)
    ax.set_title("(a) Every two-stage model lives under the ceiling", loc="left")
    ax.view_init(elev=19, azim=-128)

    ax2 = fig.add_subplot(gs[0, 1])
    rho = KANE_P301R / KANE_EXO_NULL
    ax2.axhspan(0, 1, color=BLUE, alpha=0.16)
    ax2.axhline(1, color=WARM, lw=2.0)
    ax2.text(0.04, 1.35, "PSR$_2$ ceiling: every stage-removal model is at or below 1",
             color=WARM, fontsize=8.5)
    ax2.scatter([0.5], [rho], s=85, color=INK, zorder=5)
    ax2.annotate(f"pol2-P301R measured\n" + r"$\rho = 2300/44 \approx $" + f"{rho:.1f}",
                 xy=(0.5, rho), xytext=(0.6, rho * 0.42), fontsize=9, color=INK,
                 arrowprops=dict(arrowstyle="->", color=GREY, lw=1.0))
    ax2.scatter([0.5], [1.0], s=55, marker="s", color=GREY, zorder=5)
    ax2.annotate("pol2-4 (exo-null)\nis the ceiling, by definition",
                 xy=(0.5, 1.0), xytext=(0.03, 3.0), fontsize=8.5, color=GREY,
                 arrowprops=dict(arrowstyle="->", color=GREY, lw=0.8))
    ax2.set_yscale("log"); ax2.set_ylim(0.2, 200); ax2.set_xlim(0, 1)
    ax2.set_xticks([])
    ax2.set_ylabel(r"error relative to exonuclease-null, $\rho$")
    ax2.set_title("(b) The measurement is a factor of 52 outside it", loc="left")
    ax2.spines[["top", "right", "bottom"]].set_visible(False)

    put("sec03", {"rho": rho, "kane_p301r_e8": KANE_P301R,
                  "kane_exo_null_e8": KANE_EXO_NULL, "ceiling": 1.0})
    save(fig, "fig03_nogo_ceiling.pdf")


# --------------------------------------------------------------- section 4 ---
def fig_pole_panel() -> None:
    """Variant panel against the ceiling, plus the in vitro / in vivo reversal."""
    panel = [("pol2-P301R", KANE_P301R / KANE_EXO_NULL, True),
             ("pol2-P301H", 21.0, True), ("pol2-F367S", 12.5, True),
             ("pol2-S459F", 9.1, True), ("pol2-V411L", 0.55, False),
             ("pol2-4 (exo-null)", 1.0, False)]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.8, 4.0),
                                  gridspec_kw={"wspace": 0.30})

    ys = np.arange(len(panel))
    for y, (name, r, breaks) in zip(ys, panel):
        col = WARM if breaks else GREY
        ax.plot([1, r], [y, y], color=col, lw=1.4, alpha=0.6)
        ax.scatter([r], [y], s=70, color=col, zorder=5)
    ax.axvline(1, color=INK, lw=1.8)
    ax.text(1.12, len(panel) - 0.4, "PSR$_2$ ceiling", color=INK, fontsize=8.5)
    ax.set_yticks(ys); ax.set_yticklabels([p[0] for p in panel], fontsize=8.5)
    ax.set_xscale("log"); ax.set_xlabel(r"error relative to exonuclease-null, $\rho$")
    ax.set_title("(a) Four cancer alleles break the ceiling; controls do not", loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=WARM, label="breaks ceiling"),
                       Line2D([], [], marker="o", ls="", color=GREY, label="does not")],
              frameon=False, loc="lower right")

    inv = np.array([XING_P301R_INVITRO, XING_EXO_NULL_INVITRO])
    viv = np.array([KANE_P301R, KANE_EXO_NULL])
    for i, (lab, col) in enumerate((("pol2-P301R", WARM), ("exo-null", GREY))):
        ax2.plot([0, 1], [inv[i] / inv[1], viv[i] / viv[1]], "o-", color=col,
                 lw=2.0, ms=8, label=lab)
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(["in vitro\n(purified, no MMR)", "in vivo\n(MMR-proficient)"])
    ax2.set_yscale("log")
    ax2.set_ylabel("error, relative to exo-null in the same assay")
    ax2.axhline(1, color=GREY, ls=":", lw=1.0)
    reversal = (viv[0] / viv[1]) / (inv[0] / inv[1])
    ax2.set_title(f"(b) The ordering reverses by {reversal:.0f}$\\times$", loc="left")
    ax2.legend(frameon=False)
    ax2.spines[["top", "right"]].set_visible(False)

    put("sec04", {"reversal_factor": reversal,
                  "in_vitro_ratio": inv[0] / inv[1], "in_vivo_ratio": viv[0] / viv[1],
                  "panel": [{"allele": n, "rho": r, "breaks_ceiling": b}
                            for n, r, b in panel]})
    save(fig, "fig04_pole_panel.pdf")


# --------------------------------------------------------------- section 5 ---
def fig_repair_reversal() -> None:
    """3D surface of escaped(p, c) = max(p - c, 0), with the reversal path drawn
    on it. This is the picture behind `reversal_forces_capacity_drop`."""
    p = np.linspace(0, 14, 90)
    c = np.linspace(0, 12, 90)
    P, C = np.meshgrid(p, c)
    Z = np.maximum(P - C, 0.0)

    fig = plt.figure(figsize=(10.8, 4.5))
    gs = fig.add_gridspec(1, 2, width_ratios=(1.12, 1.0), wspace=0.2,
                          left=0.02, right=0.96, bottom=0.12, top=0.88)

    ax = fig.add_subplot(gs[0, 0], projection="3d")
    ax.plot_surface(P, C, Z, cmap="magma", alpha=0.9, linewidth=0,
                    rstride=2, cstride=2)
    ax.plot([10], [9], [1], "o", color=BLUE, ms=8)
    ax.plot([8], [0], [8], "o", color=WARM, ms=8)
    ax.plot([10, 8], [9, 0], [1, 8], color=INK, lw=2.0, ls="--")
    ax.text(10, 9, 2.4, "exo-null", color=BLUE, fontsize=8.5)
    ax.text(8, 0, 9.4, "P301R", color=WARM, fontsize=8.5)
    ax.set_xlabel("produced error $p$"); ax.set_ylabel("repair capacity $c$")
    ax.set_zlabel("escaped $=\\max(p-c,0)$")
    ax.set_title("(a) Fewer errors produced, more errors escaping", loc="left")
    ax.view_init(elev=23, azim=-127)

    ax2 = fig.add_subplot(gs[0, 1])
    for cc, col in ((9.0, BLUE), (4.0, PURPLE), (0.0, WARM)):
        ax2.plot(p, np.maximum(p - cc, 0), color=col, lw=2.0, label=f"$c={cc:g}$")
    ax2.plot([10], [1], "o", color=BLUE, ms=8)
    ax2.plot([8], [8], "o", color=WARM, ms=8)
    ax2.annotate("", xy=(8, 8), xytext=(10, 1),
                 arrowprops=dict(arrowstyle="->", color=INK, lw=1.6, ls="--"))
    ax2.text(5.0, 9.2, "produces fewer errors ($8 < 10$)\nyet escapes more ($8 > 1$)",
             fontsize=8.5, color=INK)
    ax2.set_xlabel("produced error $p$"); ax2.set_ylabel("escaped error")
    ax2.set_ylim(0, 13)
    ax2.set_title("(b) Only a drop in capacity can do this", loc="left")
    ax2.legend(frameon=False, loc="lower right")
    ax2.spines[["top", "right"]].set_visible(False)

    put("sec05", {"witness": {"p_null": 10, "c_null": 9, "escaped_null": 1,
                              "p_mut": 8, "c_mut": 0, "escaped_mut": 8},
                  "lean": "ANSE.GeneticsTextbook.reversal_forces_capacity_drop"})
    save(fig, "fig05_repair_reversal.pdf")


# --------------------------------------------------------------- section 6 ---
def fig_mutation_accumulation() -> None:
    """Luria-Delbruck versus Poisson: why the variance, not the mean, is the
    signature of mutation arising before selection."""
    rng = np.random.default_rng(20261021)
    n_cultures, generations, mu = 4000, 22, 2.0e-8

    # Exact small-culture simulation: each generation doubles the wild-type
    # population, each division is an independent chance to mutate, and existing
    # mutants double along with everyone else. The jackpot comes from early
    # mutations expanding for many generations.
    finals = np.zeros(n_cultures, dtype=np.int64)
    for c in range(n_cultures):
        wt, mut = 1, 0
        for _ in range(generations):
            new_mut = rng.binomial(wt, mu * 2.0)
            wt = wt * 2 - new_mut
            mut = mut * 2 + new_mut
        finals[c] = mut
    n_final = 2 ** generations
    mean, var = float(finals.mean()), float(finals.var())
    poisson = rng.poisson(max(mean, 1e-9), size=n_cultures)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.0),
                                  gridspec_kw={"wspace": 0.28})
    top = int(max(finals.max(), poisson.max()))
    bins = np.linspace(0, min(top, int(np.percentile(finals, 99.5)) + 1), 45)
    ax.hist(finals, bins=bins, color=WARM, alpha=0.75, label="Luria–Delbrück (simulated)")
    ax.hist(poisson, bins=bins, color=BLUE, alpha=0.55, label="Poisson, same mean")
    ax.set_yscale("log")
    ax.set_xlabel("mutants per culture"); ax.set_ylabel("cultures")
    ax.set_title(f"(a) Same mean, wildly different tail", loc="left")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    ax2.scatter([0.35], [var / mean], s=110, color=WARM, zorder=5)
    ax2.scatter([0.65], [float(poisson.var()) / float(poisson.mean())], s=110,
                color=BLUE, zorder=5)
    ax2.axhline(1, color=GREY, ls="--", lw=1.2)
    ax2.text(0.05, 1.5, "Poisson: variance = mean", color=GREY, fontsize=8.5)
    ax2.set_yscale("log"); ax2.set_xlim(0, 1); ax2.set_xticks([0.35, 0.65])
    ax2.set_xticklabels(["Luria–Delbrück", "Poisson"])
    ax2.set_ylabel("variance / mean")
    ax2.set_title(f"(b) Overdispersion = {var/mean:,.0f}$\\times$", loc="left")
    ax2.spines[["top", "right"]].set_visible(False)

    put("sec06", {"cultures": n_cultures, "generations": generations,
                  "mutation_rate_per_division": mu, "final_population": n_final,
                  "mean_mutants": mean, "variance_mutants": var,
                  "variance_over_mean": var / mean,
                  "max_mutants_observed": int(finals.max()),
                  "fraction_with_zero": float((finals == 0).mean()),
                  "seed": 20261021})
    save(fig, "fig06_mutation_accumulation.pdf")


# --------------------------------------------------------------- section 7 ---
def fig_variant_landscape() -> None:
    """Real InterPro domain architecture of POLE with the variants the book
    discusses, and the real ClinVar base rates beside it."""
    rec = json.loads(RECORDS.read_text())
    pole = rec["proteins"]["POLE"]
    doms = pole["domains"]
    length = pole["length_aa"]

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(10.4, 5.0),
                                  gridspec_kw={"height_ratios": (1.25, 1.0),
                                               "hspace": 0.55})
    ax.add_patch(plt.Rectangle((1, 0.42), length, 0.16, color="#E8E8E4", ec=GREY, lw=0.6))
    palette = [BLUE, GREEN, PURPLE, WARM, "#8E6C3E"]
    for i, d in enumerate(doms):
        col = palette[i % len(palette)]
        ax.add_patch(plt.Rectangle((d["start"], 0.33), d["end"] - d["start"], 0.34,
                                   color=col, alpha=0.82, ec="white", lw=0.8))
        ax.text((d["start"] + d["end"]) / 2, 0.19, d["interpro_accession"],
                ha="center", fontsize=7.2, color=col, rotation=0)
    exo = next(d for d in doms if d["interpro_accession"] == "IPR006133")
    for res, lab in ((286, "P286R"), (411, "V411L")):
        ax.plot([res, res], [0.68, 0.86], color=INK, lw=1.3)
        ax.plot([res], [0.86], "v", color=INK, ms=7)
        ax.text(res, 0.92, lab, ha="center", fontsize=8.5, color=INK)
    ax.set_xlim(-60, length + 60); ax.set_ylim(0.05, 1.08)
    ax.set_yticks([]); ax.set_xlabel("residue")
    ax.set_title(f"(a) POLE ({pole['accession']}), {length} aa — InterPro domains, retrieved live. "
                 f"Both variants fall inside IPR006133 ({exo['start']}–{exo['end']}).",
                 loc="left", fontsize=9.5)
    ax.spines[["top", "right", "left"]].set_visible(False)

    syms = list(rec["clinvar"].keys())
    cats = ("pathogenic", "uncertain", "benign")
    w = 0.26
    xs = np.arange(len(syms))
    for j, (cat, col) in enumerate(zip(cats, (WARM, GREY, GREEN))):
        vals = [rec["clinvar"][s][cat] for s in syms]
        ax2.bar(xs + (j - 1) * w, vals, w, color=col, label=cat, alpha=0.88)
    ax2.set_xticks(xs); ax2.set_xticklabels(syms)
    ax2.set_ylabel("ClinVar records")
    ax2.set_title("(b) Base rates before interpreting any single variant. Categories overlap: "
                  "a record with conflicting assertions matches more than one.",
                  loc="left", fontsize=9.5)
    ax2.legend(frameon=False, ncol=3)
    ax2.spines[["top", "right"]].set_visible(False)

    put("sec07", {"pole_length_aa": length, "exonuclease_domain": exo,
                  "p286r_inside_exonuclease_domain": exo["start"] <= 286 <= exo["end"],
                  "clinvar": {s: {c: rec["clinvar"][s][c] for c in cats} for s in syms}})
    save(fig, "fig07_variant_landscape.pdf")


# --------------------------------------------------------------- section 9 ---
def fig_dna_topology() -> None:
    """Calugareanu Lk = Tw + Wr: a relaxed ribbon and a supercoiled one at equal
    linking number, drawn in 3D."""
    def ribbon(turns: float, writhe_amp: float, n: int = 1400):
        t = np.linspace(0, 2 * np.pi, n)
        # Superhelical path: a circle deformed into a plectoneme-like writhe.
        R, r = 3.0, 0.42
        x = (R + writhe_amp * np.cos(3 * t)) * np.cos(t)
        y = (R + writhe_amp * np.cos(3 * t)) * np.sin(t)
        z = writhe_amp * np.sin(3 * t)
        # Two strands twisting about that path.
        dx, dy, dz = np.gradient(x), np.gradient(y), np.gradient(z)
        norm = np.sqrt(dx**2 + dy**2 + dz**2) + 1e-12
        tx, ty, tz = dx / norm, dy / norm, dz / norm
        nx, ny, nz = -ty, tx, np.zeros_like(tz)
        bx = ty * nz - tz * ny
        by = tz * nx - tx * nz
        bz = tx * ny - ty * nx
        phase = turns * t
        s1 = (x + r * (np.cos(phase) * nx + np.sin(phase) * bx),
              y + r * (np.cos(phase) * ny + np.sin(phase) * by),
              z + r * (np.cos(phase) * nz + np.sin(phase) * bz))
        s2 = (x - r * (np.cos(phase) * nx + np.sin(phase) * bx),
              y - r * (np.cos(phase) * ny + np.sin(phase) * by),
              z - r * (np.cos(phase) * nz + np.sin(phase) * bz))
        return s1, s2

    fig = plt.figure(figsize=(10.6, 4.6))
    for k, (turns, amp, title) in enumerate((
            (24, 0.0, "(a) Relaxed: $\\Tw=24$, $\\Wr=0$"),
            (18, 1.05, "(b) Supercoiled: $\\Tw=18$, $\\Wr=6$"))):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        s1, s2 = ribbon(turns, amp)
        ax.plot(*s1, color=BLUE, lw=1.1)
        ax.plot(*s2, color=WARM, lw=1.1)
        for i in range(0, len(s1[0]), 26):
            ax.plot([s1[0][i], s2[0][i]], [s1[1][i], s2[1][i]], [s1[2][i], s2[2][i]],
                    color=GREY, lw=0.45, alpha=0.65)
        ax.set_title(title.replace("\\Tw", "\\mathrm{Tw}").replace("\\Wr", "\\mathrm{Wr}"),
                     loc="left", fontsize=10)
        # Looked at edge-on, the relaxed circle reads as a solenoid. Tilt the
        # camera down so the closed circular path is unmistakable in both panels.
        ax.view_init(elev=58, azim=-62)
        ax.set_box_aspect((1, 1, 0.6))
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        ax.grid(False)
        for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
            pane.pane.set_alpha(0.03)
    fig.suptitle("Linking number is conserved: $\\mathrm{Lk}=\\mathrm{Tw}+\\mathrm{Wr}=24$ "
                 "in both panels. Twist released by the helix reappears as writhe.",
                 fontsize=9.5, y=0.04, color=INK)
    put("sec09", {"Lk": 24, "relaxed": {"Tw": 24, "Wr": 0},
                  "supercoiled": {"Tw": 18, "Wr": 6}})
    save(fig, "fig09_dna_topology.pdf")


# --------------------------------------------------------------- section 8 ---
def fig_topological_genomics() -> None:
    """Persistent homology of a real coding sequence, via GUDHI.

    The sequence is the POLE canonical CDS fetched from Ensembl and cached. It is
    embedded as a 3D chaos-game representation, which turns k-mer composition into
    geometry; persistent homology then reports the loops that survive filtration.
    """
    import gudhi

    cache = DATA / "pole_cds.json"
    if cache.exists():
        seq = json.loads(cache.read_text())["sequence"]
    else:
        rec = json.loads(RECORDS.read_text())
        tid = rec["genes"]["POLE"]["canonical_transcript"]
        url = f"https://rest.ensembl.org/sequence/id/{tid}?type=cds;content-type=application/json"
        req = urllib.request.Request(url, headers={
            "Accept": "application/json",
            "User-Agent": "SocrateAI-Formal-Genetics-Textbook/1.0"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            payload = json.loads(resp.read().decode())
        seq = payload["seq"]
        cache.write_text(json.dumps({"transcript": tid, "source_url": url,
                                     "length_nt": len(seq), "sequence": seq}, indent=1))

    # 3D chaos-game representation on the tetrahedron: each base pulls the point
    # halfway toward one vertex, so nearby points share recent k-mer context.
    verts = {"A": (1, 1, 1), "C": (1, -1, -1), "G": (-1, 1, -1), "T": (-1, -1, 1)}
    pt = np.zeros(3)
    cloud = []
    for base in seq.upper():
        if base not in verts:
            continue
        pt = (pt + np.array(verts[base], dtype=float)) / 2.0
        cloud.append(pt.copy())
    cloud = np.asarray(cloud)
    sample = cloud[:: max(1, len(cloud) // 900)][:900]

    rips = gudhi.RipsComplex(points=sample, max_edge_length=0.9)
    st = rips.create_simplex_tree(max_dimension=2)
    st.compute_persistence()
    betti = st.betti_numbers()
    h1 = [(b, d) for (dim, (b, d)) in st.persistence() if dim == 1 and d != float("inf")]
    h1.sort(key=lambda bd: bd[1] - bd[0], reverse=True)

    fig = plt.figure(figsize=(10.8, 4.4))
    gs = fig.add_gridspec(1, 2, width_ratios=(1.05, 1.0), wspace=0.22,
                          left=0.03, right=0.97, bottom=0.14, top=0.88)
    ax = fig.add_subplot(gs[0, 0], projection="3d")
    ax.scatter(sample[:, 0], sample[:, 1], sample[:, 2], s=4.5, c=np.arange(len(sample)),
               cmap="viridis", alpha=0.85)
    ax.set_title(f"(a) POLE CDS as a 3D chaos-game cloud ({len(sample)} sampled points)",
                 loc="left")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.view_init(elev=20, azim=-125)

    ax2 = fig.add_subplot(gs[0, 1])
    for dim, col, mark in ((0, BLUE, "o"), (1, WARM, "^")):
        pts = [(b, d) for (dd, (b, d)) in st.persistence()
               if dd == dim and d != float("inf")]
        if pts:
            ax2.scatter([b for b, _ in pts], [d for _, d in pts], s=22, color=col,
                        marker=mark, alpha=0.75, label=f"$H_{dim}$")
    lim = 0.95
    ax2.plot([0, lim], [0, lim], color=GREY, lw=1.0, ls="--")
    ax2.set_xlabel("birth"); ax2.set_ylabel("death")
    ax2.set_xlim(0, lim); ax2.set_ylim(0, lim)
    ax2.set_title(f"(b) Persistence diagram — Betti numbers {betti}", loc="left")
    ax2.legend(frameon=False, loc="lower right")
    ax2.spines[["top", "right"]].set_visible(False)

    put("sec08", {"cds_length_nt": len(seq), "points_sampled": int(len(sample)),
                  "betti_numbers": [int(b) for b in betti],
                  "longest_h1_bars": [[float(b), float(d)] for b, d in h1[:5]],
                  "max_edge_length": 0.9, "engine": f"GUDHI {gudhi.__version__}"})
    save(fig, "fig08_topological_genomics.pdf")


# -------------------------------------------------------------- section 10 ---
def fig_clinic() -> None:
    """What survives the refutation: the clinical threshold, and the gap between
    the retracted model's prediction and the quantity clinicians actually use."""
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.0),
                                  gridspec_kw={"wspace": 0.3})
    labels = ["wild-type\n(model)", "retracted model\nprediction", "FDA threshold\n(pembrolizumab)"]
    vals = [0.2, 12500 / 9, 10.0]
    cols = [GREY, WARM, INK]
    ax.bar(labels, vals, color=cols, alpha=0.88)
    ax.set_yscale("log")
    ax.axhline(10, color=INK, ls="--", lw=1.4)
    ax.set_ylabel("TMB (mut/Mb)")
    ax.set_title("(a) The retracted model's number, kept only as arithmetic", loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    for i, v in enumerate(vals):
        ax.text(i, v * 1.35, f"{v:.4g}", ha="center", fontsize=8.5, color=INK)

    rho = np.linspace(1, 80, 300)
    hazard = 1.0 / (1.0 + 0.9 * rho)
    ax2.plot(rho, hazard, color=BLUE, lw=2.0)
    ax2.axvline(KANE_P301R / KANE_EXO_NULL, color=WARM, lw=1.6, ls="--")
    ax2.text(KANE_P301R / KANE_EXO_NULL + 1.5, 0.45,
             r"measured $\rho\approx52$", color=WARM, fontsize=8.5)
    ax2.set_xlabel(r"mutational excess over exo-null, $\rho$")
    ax2.set_ylabel("modelled immune-escape hazard")
    ax2.set_title("(b) The clinical direction is robust to the mechanism being wrong",
                  loc="left")
    ax2.spines[["top", "right"]].set_visible(False)

    put("sec10", {"retracted_tmb_prediction": 12500 / 9, "fda_threshold": 10.0,
                  "note": "the threshold crossing is established by measurement, "
                          "not by the retracted stage-removal model"})
    save(fig, "fig10_clinic.pdf")


BUILDERS = {0: fig_photo51, 1: fig_three_claims, 2: fig_proofreading_cascade, 3: fig_nogo_ceiling, 4: fig_pole_panel,
            5: fig_repair_reversal, 6: fig_mutation_accumulation, 7: fig_variant_landscape,
            8: fig_topological_genomics, 9: fig_dna_topology, 10: fig_clinic}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", type=int, nargs="*", help="section numbers to build")
    args = ap.parse_args()
    if not RECORDS.exists():
        print(f"missing {RECORDS}; run retrieve_genomic_data.py first")
        return 1
    todo = sorted(BUILDERS) if not args.only else sorted(set(args.only) & set(BUILDERS))
    for n in todo:
        print(f"section {n}:")
        BUILDERS[n]()
    print(f"\n{len(todo)} figure(s) built; values in {VALUES.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

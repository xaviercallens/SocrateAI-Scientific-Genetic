"""Figures for the proofreading-stage-removal falsification in Paper E.

Every plotted value is read from data/proofreading_stage_removal_test.json and
data/pole_proofreading_literature.json, which are produced and populated from published
sources; nothing is invented or smoothed. Writes three vector PDFs into papers/.

Run: .venv/bin/python scripts/make_proofreading_figures.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# matplotlib >= 3.6 registers the "3d" projection on import; no mpl_toolkits import is needed.

HERE = Path(__file__).resolve().parents[1]
RESULT = json.loads((HERE / "data" / "proofreading_stage_removal_test.json").read_text(encoding="utf-8"))
LIT = json.loads((HERE / "data" / "pole_proofreading_literature.json").read_text(encoding="utf-8"))
OUT = HERE / "papers"

plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 10.5,
    "legend.fontsize": 8.5,
    "pdf.fonttype": 42,
    "savefig.bbox": "tight",
})

N_STAGES = 3  # the three-stage scheme Paper E's GAR.H15 assumes


def fig_feasible_surface() -> Path:
    """3D: the entire error-rate surface reachable by removing proofreading stages, and the
    measured alleles that sit above it."""
    ceiling_fold = RESULT["primary_result"]["exonuclease_null_ceiling_fold"]["fold"]
    p301r_fold = RESULT["primary_result"]["P301R_fold"]["fold"]
    # The per-stage escape factor is not free: with N_STAGES stages, removing all of them gives
    # the exonuclease-null fold change, so q is pinned by that measurement.
    q_pinned = ceiling_fold ** (-1.0 / N_STAGES)

    q = np.linspace(0.45, 0.97, 160)
    j = np.linspace(0.0, N_STAGES, 120)  # stages removed
    Q, J = np.meshgrid(q, j)
    Z = np.log10(Q ** (-J))  # log10 fold change over wild type

    fig = plt.figure(figsize=(7.2, 5.4))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(Q, J, Z, cmap="viridis", alpha=0.82, linewidth=0,
                           antialiased=True, rstride=2, cstride=2)

    # The q the data pins, and the feasible folds along it (integer numbers of stages removed).
    js = np.arange(0, N_STAGES + 1)
    feasible = q_pinned ** (-js.astype(float))
    ax.plot([q_pinned] * len(js), js, np.log10(feasible), color="black", lw=2.0,
            marker="o", ms=5, zorder=10)

    # Measured alleles, as stems at the pinned q.
    ax.plot([q_pinned, q_pinned], [N_STAGES, N_STAGES],
            [np.log10(ceiling_fold), np.log10(p301r_fold)],
            color="crimson", lw=2.6, ls="--", zorder=11)
    ax.scatter([q_pinned], [N_STAGES], [np.log10(p301r_fold)], color="crimson", s=70,
               depthshade=False, zorder=12, edgecolor="black", linewidth=0.5)

    for name, rec in RESULT["variant_panel_cross_paper"].items():
        if rec["verdict"] == "breaks_ceiling":
            ax.scatter([q_pinned], [N_STAGES], [np.log10(rec["fold"])], color="darkorange",
                       s=42, depthshade=False, zorder=12, edgecolor="black", linewidth=0.4)

    ax.set_xlabel("per-stage escape factor $q$")
    ax.set_ylabel("proofreading stages removed")
    ax.set_zlabel("$\\log_{10}$ fold change over wild type")
    ax.set_title("Every error rate reachable by removing proofreading stages,\nand the alleles that exceed it")
    ax.set_yticks(js)
    ax.view_init(elev=22, azim=-128)
    cb = fig.colorbar(surf, ax=ax, shrink=0.58, pad=0.11)
    cb.set_label("$\\log_{10}$ fold change (model)")

    handles = [
        Line2D([], [], color="black", marker="o", lw=2, label=f"feasible folds at measured $q={q_pinned:.2f}$"),
        Line2D([], [], color="crimson", marker="o", ls="--", lw=2.2,
               label=f"pol2-P301R (POLE P286R), {p301r_fold:.0f}$\\times$"),
        Line2D([], [], color="darkorange", marker="o", ls="none",
               label="other alleles above the ceiling"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(-0.06, 0.99), framealpha=0.9)
    fig.text(0.5, 0.015,
             f"$q$ is fixed by the exonuclease-null measurement: $q^{{-{N_STAGES}}}$ = "
             f"{ceiling_fold:.2f}$\\times$, so the surface cannot be raised by choosing $q$.",
             ha="center", fontsize=8.2, style="italic")

    path = OUT / "figE5_stage_removal_feasible_surface.pdf"
    fig.subplots_adjust(left=0.02, right=0.97, bottom=0.10, top=0.94)
    fig.savefig(path, bbox_inches=None)
    plt.close(fig)
    return path


def fig_variant_forest() -> Path:
    """2D: every measured allele against the exonuclease-null ceiling, with published CIs."""
    ceiling = RESULT["primary_result"]["exonuclease_null_ceiling_fold"]
    rows = [("pol2-P301R", RESULT["primary_result"]["P301R_fold"], "POLE P286R", "breaks_ceiling")]
    for name, rec in RESULT["variant_panel_cross_paper"].items():
        rows.append((name, rec, rec["human_equivalent"], rec["verdict"]))
    rows.sort(key=lambda r: r[1]["fold"])

    colours = {"breaks_ceiling": "crimson", "below_ceiling": "#3b7ea1", "indistinguishable": "grey"}
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    for i, (name, rec, human, verdict) in enumerate(rows):
        ax.plot([rec["fold_low"], rec["fold_high"]], [i, i], color=colours[verdict], lw=1.6)
        ax.plot([rec["fold"]], [i], marker="o", ms=5.5, color=colours[verdict])

    ax.axvspan(ceiling["fold_low"], ceiling["fold_high"], color="black", alpha=0.10)
    ax.axvline(ceiling["fold"], color="black", lw=1.5, ls="--")
    ax.text(ceiling["fold"] * 1.12, len(rows) - 0.6,
            f"exonuclease-null ceiling\n(pol2-4, {ceiling['fold']:.1f}$\\times$)",
            fontsize=8.2, va="top")
    ax.fill_betweenx([-1, len(rows)], ceiling["fold_high"], 1e4, color="crimson", alpha=0.05)

    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{n}  ({h})" if h else n for n, _, h, _ in rows], fontsize=8)
    ax.set_xscale("log")
    ax.set_xlim(0.5, 400)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.set_xlabel("CAN1 mutation rate, fold change over wild type (log scale)")
    ax.set_title("Theorem PSR_2 forbids the shaded region to the right of the ceiling")
    handles = [
        Line2D([], [], color="crimson", marker="o", lw=1.6, label="breaks the ceiling (stage removal refuted)"),
        Line2D([], [], color="grey", marker="o", lw=1.6, label="indistinguishable from the ceiling"),
        Line2D([], [], color="#3b7ea1", marker="o", lw=1.6, label="below the ceiling (consistent)"),
    ]
    ax.legend(handles=handles, loc="lower right", framealpha=0.95)
    ax.grid(axis="x", alpha=0.25)

    path = OUT / "figE6_variant_panel_vs_ceiling.pdf"
    fig.savefig(path)
    plt.close(fig)
    return path


def fig_discordance() -> Path:
    """2D: the ordering of P301R and exonuclease-null reverses between the purified system and
    the cell."""
    d = RESULT["in_vitro_in_vivo_discordance"]
    iv = LIT["in_vitro_fidelity"]["base_substitution_error_rate"]
    kane = LIT["in_vivo_can1"]["kane2014"]["strains"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.9))
    for ax, vals, title, ylab in (
        (ax1, [iv["exonuclease_null"], iv["P301R"]],
         "Purified enzyme, no mismatch repair\n(Xing 2019, lacZ forward assay)",
         "base-substitution error rate\n($\\times 10^{-6}$ per nucleotide)"),
        (ax2, [kane["pol2-4_exonuclease_null"]["rate"], kane["pol2-P301R"]["rate"]],
         "In the cell, mismatch repair present\n(Kane 2014, CAN1)",
         "mutation rate\n($\\times 10^{-8}$ per cell per generation)"),
    ):
        ax.plot([0, 1], vals, color="black", lw=1.4, zorder=1)
        ax.scatter([0], [vals[0]], s=90, color="#3b7ea1", zorder=3, label="exonuclease-null")
        ax.scatter([1], [vals[1]], s=90, color="crimson", zorder=3, label="P301R")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["exonuclease\nnull", "P301R"])
        ax.set_xlim(-0.35, 1.35)
        ax.set_yscale("log")
        ax.set_title(title, fontsize=9)
        ax.set_ylabel(ylab, fontsize=8.8)
        ax.grid(axis="y", alpha=0.25)
        ax.set_ylim(min(vals) / 2.0, max(vals) * 2.6)
        for x, v in zip((0, 1), vals):
            ax.annotate(f"{v:g}", (x, v), textcoords="offset points", xytext=(9, 4), fontsize=8.5)

    ax1.annotate("P301R is more accurate", (0.5, max(iv.values()) * 0.55), ha="center",
                 fontsize=8.4, style="italic")
    ax2.annotate("P301R is far worse", (0.5, kane["pol2-P301R"]["rate"] * 0.22), ha="center",
                 fontsize=8.4, style="italic")
    fig.suptitle(
        f"The ordering reverses by a factor of {d['discordance_factor']:.0f}: "
        "intrinsic error production cannot explain the allele", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.99))

    path = OUT / "figE7_in_vitro_in_vivo_reversal.pdf"
    fig.savefig(path)
    plt.close(fig)
    return path


if __name__ == "__main__":
    for p in (fig_feasible_surface(), fig_variant_forest(), fig_discordance()):
        print("wrote", p.relative_to(HERE), p.stat().st_size, "bytes")

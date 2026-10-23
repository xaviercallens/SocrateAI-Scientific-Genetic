#!/usr/bin/env python3
"""Fit a surrogate for the four-point margin's null floor, and test it honestly.

Phase 3 of `docs/GENETICA_RETICULATA.md` -- the margin's null distribution -- is
the oldest open item in the programme, and the literature review found no
analytic sampling distribution for it either. Every calibrated answer currently
costs a matched simulation: `real-data` spends a bisection plus five tree
simulations per population, each calibrated statistic spends 99 replicates per
dataset-seed, and two experiments reported counts against nulls that could not
vary because there was no cheap way to know in advance where the null is
degenerate.

This fits the floor as a function of the parameters that generate it, from the
table `reticulation null-sweep` produces.

Three commitments, which are what make a fitted model admissible here at all
(`docs/RETICULATA_NUMERICS_TOPOLOGY_ML.md` section 1):

* **Only the floor is predicted.** The certificate's margin stays an exact
  integer comparison the Lean kernel can check. A wrong surrogate costs power or
  wastes a query; it cannot make a certificate wrong.
* **The model is a scaling law, not a black box.** `log(floor)` is regressed on
  `log` of the length, the divergence and the quadruple count, plus the gap load
  and clock spread. That is interpretable, and it is *falsifiable*: the
  length exponent has an independent prediction. The earlier length-scaling
  diagnostic found the clonal margin shrinking as `L^(-1/2)`, the signature of
  sampling noise, so a fit that recovers an exponent near `-0.5` confirms the law
  rather than merely interpolating the table. An exponent far from it would mean
  the floor is not noise-dominated and the surrogate rests on nothing.
* **Held-out *regions*, not a random split.** The failure mode that matters is
  extrapolation, so the model is tested on the longest alignments and the largest
  taxon count after being trained without them. A random split would hide exactly
  the error a deployment would hit.

Usage:

    python scripts/reticulata_null_surrogate.py \\
      --sweep data/reticulata_phase3/null_sweep.json \\
      --out data/reticulata_phase3/null_surrogate.json
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

# The design matrix's columns, in order. Named so the fitted exponents can be
# read against their predictions rather than inspected as a vector.
FEATURES: Tuple[str, ...] = (
    "intercept",
    "log_length",
    "log_mean_p_distance",
    "log_quadruples",
    "gap_load",
    "clock_spread",
)

# The independent prediction for the length exponent: sampling noise in a
# difference of multinomial counts scales as L^(-1/2).
PREDICTED_LENGTH_EXPONENT: float = -0.5
# How far from it the fit may land before the scaling law is judged unsupported.
LENGTH_EXPONENT_TOLERANCE: float = 0.15


def load_rows(sweep: Path) -> List[Dict[str, float]]:
    """One row per replicate, dropping replicates with no violating quadruple.

    A replicate with no violation has a floor of zero, which has no logarithm and
    is not a floor the surrogate can predict. Those are the degenerate region,
    and they are reported separately rather than imputed.
    """
    payload = json.loads(sweep.read_text())
    rows: List[Dict[str, float]] = []
    for cell in payload["cells"]:
        for index, margin in enumerate(cell["max_margin"]):
            if margin <= 0.0:
                continue
            divergence = cell["mean_p_distance"][index]
            if divergence <= 0.0:
                continue
            rows.append(
                {
                    "taxa": float(cell["taxa"]),
                    "quadruples": float(cell["quadruples"]),
                    "length": float(cell["length"]),
                    "mu": float(cell["mu"]),
                    "gap_load": float(cell["gap_load"]),
                    "clock_spread": float(cell["clock_spread"]),
                    "mean_p_distance": float(divergence),
                    "margin": float(margin),
                }
            )
    return rows


def design(rows: Sequence[Dict[str, float]]) -> Tuple[np.ndarray, np.ndarray]:
    """The design matrix and the response `log(floor)`."""
    matrix = np.column_stack(
        [
            np.ones(len(rows)),
            np.log([r["length"] for r in rows]),
            np.log([r["mean_p_distance"] for r in rows]),
            np.log([r["quadruples"] for r in rows]),
            np.array([r["gap_load"] for r in rows]),
            np.array([r["clock_spread"] for r in rows]),
        ]
    )
    response = np.log([r["margin"] for r in rows])
    return matrix, response


def fit(rows: Sequence[Dict[str, float]]) -> np.ndarray:
    """Ordinary least squares on the log-log design."""
    matrix, response = design(rows)
    coefficients, *_ = np.linalg.lstsq(matrix, response, rcond=None)
    return coefficients


def relative_errors(
    coefficients: np.ndarray, rows: Sequence[Dict[str, float]]
) -> np.ndarray:
    """`|predicted - actual| / actual` per row, in the original units."""
    matrix, response = design(rows)
    predicted = np.exp(matrix @ coefficients)
    actual = np.exp(response)
    return np.abs(predicted - actual) / actual


def summarise(errors: np.ndarray) -> Dict[str, float]:
    """Error distribution, reported rather than reduced to one number."""
    if errors.size == 0:
        return {"rows": 0}
    return {
        "rows": int(errors.size),
        "median": float(np.median(errors)),
        "mean": float(np.mean(errors)),
        "p90": float(np.percentile(errors, 90)),
        "max": float(np.max(errors)),
    }


def held_out_region(
    rows: Sequence[Dict[str, float]], key: str, test_value: float
) -> Tuple[List[Dict[str, float]], List[Dict[str, float]]]:
    """Split by a parameter *region*: train below `test_value`, test at it.

    Not a random split. The error that matters is extrapolation into parameters
    the model has never seen, and a random split would hide it.
    """
    train = [r for r in rows if r[key] < test_value]
    test = [r for r in rows if r[key] == test_value]
    return train, test


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.sweep.read_text())
    rows = load_rows(args.sweep)
    if len(rows) < len(FEATURES) * 10:
        raise RuntimeError(
            f"only {len(rows)} usable replicates for {len(FEATURES)} features; "
            "the sweep is too small to fit anything defensible"
        )

    full = fit(rows)
    in_sample = summarise(relative_errors(full, rows))
    exponents = dict(zip(FEATURES, (float(c) for c in full)))

    # The falsifiable part: does the length exponent recover the noise prediction?
    length_exponent = exponents["log_length"]
    law_supported = (
        abs(length_exponent - PREDICTED_LENGTH_EXPONENT) <= LENGTH_EXPONENT_TOLERANCE
    )

    # Held-out regions, in the two dimensions a deployment would extrapolate in.
    regions: Dict[str, Dict[str, object]] = {}
    for key, label in (("length", "longest alignments"), ("taxa", "largest taxon count")):
        values = sorted({r[key] for r in rows})
        if len(values) < 2:
            continue
        test_value = values[-1]
        train, test = held_out_region(rows, key, test_value)
        if not train or not test:
            continue
        local = fit(train)
        regions[key] = {
            "held_out_value": test_value,
            "description": f"trained without the {label}, tested only on them",
            "train_rows": len(train),
            "test_rows": len(test),
            "in_sample": summarise(relative_errors(local, train)),
            "held_out": summarise(relative_errors(local, test)),
            "length_exponent_without_this_region": float(
                dict(zip(FEATURES, (float(c) for c in local)))["log_length"]
            ),
        }

    degenerate = payload["degenerate_cells"]
    worst_held_out = max(
        (float(r["held_out"]["median"]) for r in regions.values() if r["held_out"]["rows"]),
        default=float("nan"),
    )

    report = {
        "what_this_is": (
            "A surrogate for the four-point margin's null floor, fitted to "
            "`reticulation null-sweep`. Only the floor is predicted; the "
            "certificate's margin stays an exact integer comparison, so a wrong "
            "surrogate costs power or wastes a query and cannot make a "
            "certificate wrong."
        ),
        "model": {
            "form": "log(floor) ~ intercept + b*log(L) + c*log(mean p-distance) "
            "+ d*log(quadruples) + e*gap_load + f*clock_spread",
            "why_a_scaling_law": (
                "Interpretable and falsifiable. The length exponent has an "
                "independent prediction: the earlier length-scaling diagnostic "
                "found the clonal margin shrinking as L^(-1/2), the signature of "
                "sampling noise. A fit recovering an exponent near -0.5 confirms "
                "the law; one far from it would mean the floor is not "
                "noise-dominated and this surrogate rests on nothing."
            ),
            "fitted_exponents": exponents,
            "rows_fitted": len(rows),
        },
        "the_falsifiable_check": {
            "predicted_length_exponent": PREDICTED_LENGTH_EXPONENT,
            "fitted_length_exponent": length_exponent,
            "tolerance": LENGTH_EXPONENT_TOLERANCE,
            "scaling_law_supported": bool(law_supported),
            "reading": (
                "This is the check that distinguishes a principled surrogate from "
                "a curve through a table. It was fixed before the fit was run."
            ),
        },
        "in_sample_relative_error": in_sample,
        "held_out_regions": regions,
        "worst_held_out_median_relative_error": worst_held_out,
        "degenerate_region": {
            "cells": degenerate["count"],
            "why_it_is_reported_not_imputed": (
                "A cell where no replicate produced a violating quadruple has a "
                "floor of zero, which has no logarithm and is not a floor any "
                "model should invent. Knowing which parameters are degenerate "
                "before running an experiment is the most useful thing this table "
                "buys -- it is exactly what was missing when the topology and "
                "breakpoint experiments reported counts against empty nulls."
            ),
            "parameters": degenerate["parameters"],
        },
        "calibrated_range": {
            "why": "A surrogate must refuse to predict outside the parameters it "
            "was fitted on. These are the bounds; anything outside them gets a "
            "simulation, not a prediction.",
            "length": [min(r["length"] for r in rows), max(r["length"] for r in rows)],
            "taxa": [min(r["taxa"] for r in rows), max(r["taxa"] for r in rows)],
            "mean_p_distance": [
                min(r["mean_p_distance"] for r in rows),
                max(r["mean_p_distance"] for r in rows),
            ],
            "gap_load": [min(r["gap_load"] for r in rows), max(r["gap_load"] for r in rows)],
            "clock_spread": [
                min(r["clock_spread"] for r in rows),
                max(r["clock_spread"] for r in rows),
            ],
        },
        "how_to_use_it": (
            "Refuse to predict outside the convex hull of the swept parameters, "
            "and keep the surrogate honest at runtime by simulating a handful of "
            "replicates and checking they fall where it says. A surrogate that "
            "silently extrapolates is the unexamined step this project keeps "
            "finding in its own apparatus."
        ),
    }
    args.out.write_text(json.dumps(report, indent=2) + "\n")

    print(f"rows fitted              : {len(rows)}")
    print(f"length exponent          : {length_exponent:+.3f} "
          f"(predicted {PREDICTED_LENGTH_EXPONENT:+.2f}, "
          f"supported={law_supported})")
    print(f"in-sample median error   : {in_sample['median']:.1%}")
    for key, region in regions.items():
        held = region["held_out"]
        print(
            f"held out {key:>6} = {region['held_out_value']:<7g}: "
            f"median {held['median']:.1%}, p90 {held['p90']:.1%} "
            f"({held['rows']} rows)"
        )
    print(f"degenerate cells         : {degenerate['count']}")
    print(
        "calibrated divergence    : "
        f"{min(r['mean_p_distance'] for r in rows):.5f} .. "
        f"{max(r['mean_p_distance'] for r in rows):.5f}"
    )
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

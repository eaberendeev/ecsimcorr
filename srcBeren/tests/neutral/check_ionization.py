#!/usr/bin/env python3
"""PASS/FAIL for beam ionization of a neutral gas with equal and unequal macro weights.

Expected neutral density: n(t) = n(0) exp(-int nu dt), nu(t) recorded from the beam state (trapezoid rule).
Physical balance: neutrals removed = electrons created = ions created, within sampling noise; with equal
product weights electrons and ions are created in pairs (exact charge).
"""
import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile

CASES = (
    "proton_equal",
    "proton_beam_heavier_x4",
    "proton_beam_lighter_x4",
    "proton_dilute_neutrals",
    "electron_equal",
    "electron_heavy_products",
    "electron_unpaired_products",
)
# Neutral density against n(0) exp(-int nu dt), in standard deviations of the sampling: one reaction removes
# k = w_beam/w_neutral neutral macroparticles on average (integer part plus a Bernoulli remainder).
CURVE_SIGMA = 3.0           # RMS over the samples
FINAL_SIGMA = 5.0           # last sample
BALANCE_SIGMA = 5.0         # products vs removed neutrals
MINIMUM_DECAY = 0.75        # informative: the neutrals must lose at least this fraction


def check_case(directory, name):
    params = json.loads((directory / f"{name}_params.json").read_text())
    with (directory / f"{name}.csv").open(newline="") as stream:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]
    w = params["weights"]
    n0 = rows[0]["neutrals"]
    if not math.isclose(n0, params["neutral_density"], rel_tol=1e-9):
        raise ValueError("initial neutral density disagrees with the configuration")
    integral, expected = 0.0, [n0]
    for a, b in zip(rows, rows[1:]):
        integral += 0.5 * (a["rate"] + b["rate"]) * (b["time"] - a["time"])
        expected.append(n0 * math.exp(-integral))
    ratio = w["beam"] / w["neutrals"]
    whole, remainder = math.floor(ratio), ratio - math.floor(ratio)
    second_moment = (1 - remainder) * whole**2 + remainder * (whole + 1) ** 2   # E[k^2]
    z = []
    for r, e in zip(rows[1:], expected[1:]):
        reactions = (n0 - e) / w["beam"]
        z.append((r["neutrals"] - e) / (w["neutrals"] * math.sqrt(reactions * second_moment)))
    curve = math.sqrt(sum(v * v for v in z) / len(z))
    relative = max(abs(r["neutrals"] / e - 1) for r, e in zip(rows, expected))
    final = rows[-1]
    removed = n0 - final["neutrals"]
    decay = removed / n0
    # One reaction moves w_beam physical particles; each species draws w_beam/w_s macroparticles.
    # Variance of the products-minus-removed difference from the Bernoulli remainders, per reaction.
    reactions = removed / w["beam"]
    def remainder_variance(weight):
        ratio = w["beam"] / weight
        p = ratio - math.floor(ratio)
        return p * (1 - p) * weight**2
    errors = {}
    for product in ("electrons", "ions"):
        sigma = math.sqrt(reactions * (remainder_variance(w[product]) + remainder_variance(w["neutrals"])))
        sigma = max(sigma, 1e-12 * n0)
        errors[product] = (final[product] - removed) / sigma
    paired = w["electrons"] == w["ions"]
    pairs_exact = final["electron_macros"] == final["ion_macros"] and final["electrons"] == final["ions"]
    passed = (
        curve <= CURVE_SIGMA
        and abs(z[-1]) <= FINAL_SIGMA
        and decay >= MINIMUM_DECAY
        and all(abs(z) <= BALANCE_SIGMA for z in errors.values())
        and (pairs_exact or not paired)
    )
    return {
        "case": name,
        "pass": passed,
        "neutral_curve_rms_sigma": curve,
        "neutral_final_sigma": z[-1],
        "neutral_max_relative_error": relative,
        "neutral_decay": decay,
        "electron_balance_sigma": errors["electrons"],
        "ion_balance_sigma": errors["ions"],
        "paired_products": paired,
        "pairs_exact": pairs_exact,
        "weights_over_neutral": {k: v / w["neutrals"] for k, v in w.items()},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=Path.cwd())
    parser.add_argument("--binary", type=Path, help="run the bounded --quick experiment first")
    args = parser.parse_args()
    if args.binary:
        with tempfile.TemporaryDirectory(prefix="ionization-", dir=args.directory) as temporary:
            subprocess.run([str(args.binary.resolve()), "--quick"], cwd=temporary, check=True, timeout=600)
            for path in Path(temporary).iterdir():
                shutil.copyfile(path, args.directory / path.name)
    results = []
    for name in CASES:
        try:
            r = check_case(args.directory, name)
            print(("PASS" if r["pass"] else "FAIL") + f" {name}: neutrals RMS {r['neutral_curve_rms_sigma']:.2f} sigma "
                  f"(limit {CURVE_SIGMA}), final {r['neutral_final_sigma']:+.2f} sigma (limit {FINAL_SIGMA}), max "
                  f"{r['neutral_max_relative_error']:.1%}, decay={r['neutral_decay']:.0%}, electrons {r['electron_balance_sigma']:+.2f} sigma, "
                  f"ions {r['ion_balance_sigma']:+.2f} sigma (limit {BALANCE_SIGMA}), "
                  f"pairs exact={r['pairs_exact']} (required: {r['paired_products']})")
        except (OSError, ValueError, KeyError, ZeroDivisionError) as error:
            r = {"case": name, "pass": False, "error": str(error)}
            print(f"FAIL {name}: {error}")
        results.append(r)
    (args.directory / "ionization_report.json").write_text(json.dumps(results, indent=2) + "\n")
    return 0 if all(r["pass"] for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Compare finite/completed Weibel runs, report growth rates, and draw figures.

The growth window and dispersion relation are inherited from weibel_check.py.
Particle-count-dependent saturation/noise limits of the 500-ppc fixture are
not transferred to reduced-particle experiments.
"""
import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from weibel_check import gamma_theory, GAMMA_FIT_WIN, GAMMA_WB_FIT_WIN, TOL, TOL_ENERGY_CONS


def fit_growth(time, values, window, factor=1.0):
    mask = (time >= window[0]) & (time <= window[1]) & (values > 0)
    if mask.sum() < 5:
        raise ValueError(f"not enough positive samples for growth fit in {window}")
    x, y = time[mask], np.log(values[mask])
    coefficients, covariance = np.polyfit(x, y, 1, cov=True)
    residual = y - np.polyval(coefficients, x)
    variation = np.sum((y - y.mean()) ** 2)
    return {"gamma": float(coefficients[0] / factor),
            "ols_stderr": float(np.sqrt(covariance[0, 0]) / factor),
            "r_squared": float(1 - np.sum(residual ** 2) / variation) if variation else 0,
            "samples": int(mask.sum()), "window": list(window)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdir", type=Path)
    parser.add_argument("--partial", action="store_true", help="analyze only completed runs")
    parser.add_argument("--require-growth", action="store_true", help="require the original linear Weibel growth bands")
    args = parser.parse_args()
    root = args.workdir.resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    theory11 = float(gamma_theory(2 * np.pi * np.sqrt(2) / 9))
    theory_max = float(gamma_theory(1.16))
    arrays, summary, fingerprints, moments = {}, {}, set(), {}
    for name in manifest["variants"]:
        if args.partial and name not in manifest["results"]:
            continue
        result = manifest["results"][name]
        if result["exit_code"] != 0:
            raise ValueError(f"{name} did not finish successfully")
        state = result["initial_state"]
        fingerprints.add((state["fingerprint"], state["particles"]))
        data = np.atleast_1d(np.genfromtxt(root / name / "comparison.csv", delimiter=",", names=True))
        if not all(np.all(np.isfinite(data[column])) for column in data.dtype.names):
            raise ValueError(f"{name}: non-finite observable")
        expected = manifest["steps"] * manifest["dt"]
        if not np.isclose(data["time"][-1], expected, rtol=1e-12) or not np.all(np.diff(data["time"]) > 0):
            raise ValueError(f"{name}: incomplete or nonmonotonic times")
        max_energy_error = float(np.max(np.abs(data["relative_energy_error"])))
        # Same relative threshold as the existing integration energy check.
        if max_energy_error > 1e-3:
            raise ValueError(f"{name}: energy error {max_energy_error} exceeds 1e-3")
        mode_fit = fit_growth(data["time"], data["mode11"], GAMMA_FIT_WIN)
        magnetic_fit = fit_growth(data["time"], data["magnetic"], GAMMA_WB_FIT_WIN, factor=2)
        sample100 = np.argmin(abs(data["time"] - 100))
        energy_table = np.atleast_1d(np.genfromtxt(root / name / "energy.txt", names=True))
        if not all(np.all(np.isfinite(energy_table[column])) for column in energy_table.dtype.names):
            raise ValueError(f"{name}: non-finite energy diagnostic")
        if not np.isclose(energy_table["Time"][-1], expected, rtol=1e-12):
            raise ValueError(f"{name}: incomplete energy diagnostic")
        initial_total = float(energy_table["Electrons"][0] + energy_table["energyFieldE"][0] + energy_table["energyFieldB"][0])
        step_absolute = float(np.max(np.abs(energy_table["energyConserve"])))
        step_residual = step_absolute / initial_total
        if step_absolute > TOL_ENERGY_CONS:
            raise ValueError(f"{name}: absolute step energy residual {step_absolute} exceeds original Weibel tolerance {TOL_ENERGY_CONS}")
        anisotropy = 2 * energy_table["ElectronsZ"] / energy_table["ElectronsXY"]
        if not np.all(np.isfinite(anisotropy)):
            raise ValueError(f"{name}: non-finite temperature anisotropy")
        summary[name] = {
            "gamma11": mode_fit["gamma"], "gamma11_over_theory": mode_fit["gamma"] / theory11,
            "gamma11_fit": mode_fit, "gamma_magnetic_fit": magnetic_fit,
            "anisotropy_initial": float(anisotropy[0]),
            "anisotropy_at_100": float(anisotropy[np.argmin(abs(energy_table["Time"] - 100))]),
            "anisotropy_final": float(anisotropy[-1]),
            "gamma11_in_original_nocorr_band": bool(TOL["nocorr"]["gamma_rel"][0] * theory11 < mode_fit["gamma"] < TOL["nocorr"]["gamma_rel"][1] * theory11),
            "gamma_wb_in_original_nocorr_band": bool(TOL["nocorr"]["gamma_wb"][0] * theory_max < magnetic_fit["gamma"] < TOL["nocorr"]["gamma_wb"][1] * theory_max),
            "max_relative_energy_error": max_energy_error,
            "max_relative_step_residual": step_residual,
            "max_absolute_step_residual": step_absolute,
            "density_rms_at_100": float(data["density_rms"][sample100]),
            "density_rms_initial": float(data["density_rms"][0]),
            "density_sample_time": float(data["time"][sample100]),
            "gauss_rms_final": float(data["gauss_rms"][-1]),
            "gauss_rms_initial": float(data["gauss_rms"][0]),
            "gauss_drift_rms_final": float(data["gauss_drift_rms"][-1]),
            "continuity_rms_final": float(data["continuity_interval_rms"][-1]),
            "magnetic_peak": float(data["magnetic"].max()),
            "magnetic_peak_time": float(data["time"][np.argmax(data["magnetic"])]),
            "magnetic_final": float(data["magnetic"][-1]),
            "electric_early_mean": float(np.mean(energy_table["energyFieldE"][(energy_table["Time"] > 50) & (energy_table["Time"] < 150)])),
            "seconds": result["seconds"],
        }
        arrays[name] = data
        moments[name] = (energy_table["Time"], anisotropy)
    if len(fingerprints) != 1:
        raise ValueError("initial particle fingerprints differ")
    report = {"ppc": manifest["ppc"], "particles": next(iter(fingerprints))[1],
              "final_time": expected, "theory_gamma11": theory11, "theory_gamma_max": theory_max,
              "same_initial_particles": True, "partial": args.partial, "require_growth": args.require_growth, "variants": summary}
    (root / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    columns = ["variant", "gamma11", "gamma11_over_theory", "anisotropy_at_100", "electric_early_mean", "gauss_rms_final", "gauss_drift_rms_final",
               "density_rms_at_100", "max_relative_energy_error", "magnetic_peak", "seconds"]
    with (root / "summary.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for name, values in summary.items():
            writer.writerow({"variant": name, **{key: values[key] for key in columns[1:]}})
    lines = [f"Weibel comparison: {manifest['ppc']} particles/cell, {report['particles']} particles, t={expected:g}.",
             f"Same initial particles: yes. Continuum Gamma11={theory11:.6f}. Fit window: t=72..270.", "",
             "| Variant | Gamma11 | Fit R² | Density RMS at t≈100 | Final Gauss RMS | Max relative energy drift |",
             "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for name, value in summary.items():
        lines.append(f"| {name} | {value['gamma11']:.6f} | {value['gamma11_fit']['r_squared']:.4f} | "
                     f"{100 * value['density_rms_at_100']:.3f}% | {value['gauss_rms_final']:.6g} | "
                     f"{value['max_relative_energy_error']:.3g} |")
    lines += ["", "Density RMS is measured on the 3D physical grid. Gauss includes the neutralizing background.",
              "The E/rho time staggering is unchanged from ECSIM. Lower residuals do not imply exact charge conservation.",
              "OLS fit errors do not replace a seed ensemble or resolution study."]
    if args.partial:
        lines.append("This report contains only the completed subset of runs.")
    (root / "report.md").write_text("\n".join(lines) + "\n")
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    for name, data in arrays.items():
        t = data["time"]
        for axis, column in zip(axes.flat, ("magnetic", "electric", "mode11", "gauss_drift_rms", "density_rms", "relative_energy_error")):
            values = np.abs(data[column]) if column == "relative_energy_error" else data[column]
            axis.semilogy(t, np.where(values > 0, values, np.nan), label=name)
            axis.set_xlabel(r"$t\omega_p$")
            axis.set_title(column.replace("_", " "))
            axis.grid(alpha=0.3)
    axes[0, 2].axvspan(*GAMMA_FIT_WIN, alpha=0.1, color="black")
    first = next(iter(arrays.values()))
    anchor = np.argmin(abs(first["time"] - GAMMA_FIT_WIN[0]))
    reference_time = np.linspace(*GAMMA_FIT_WIN, 50)
    axes[0, 2].semilogy(reference_time,
                        first["mode11"][anchor] * np.exp(theory11 * (reference_time - first["time"][anchor])),
                        "k--", label="continuum growth (arbitrary amplitude)")
    axes[0, 2].legend(fontsize=7)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(f"Weibel, 30x30x5, {manifest['ppc']} particles/cell, dt={manifest['dt']}; identical particles")
    fig.savefig(root / "comparison.png", dpi=160)
    fig.savefig(root / "comparison.pdf")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    names = list(summary)
    gammas = [summary[name]["gamma11"] for name in names]
    errors = [summary[name]["gamma11_fit"]["ols_stderr"] for name in names]
    ax.errorbar(np.arange(len(names)), gammas, yerr=errors, fmt="o", capsize=4)
    ax.axhline(theory11, color="black", linestyle="--", label=f"continuum theory: {theory11:.5f}")
    ax.set_xticks(np.arange(len(names)), names, rotation=20)
    ax.set_ylabel(r"$\Gamma_{11}/\omega_p$")
    ax.set_title("Mode (1,1): log-amplitude fit at t=72..270; bars are OLS fit errors")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.savefig(root / "growth_rates.png", dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for name, (time, anisotropy) in moments.items():
        ax.plot(time, anisotropy, label=name)
    ax.set_xlabel(r"$t\omega_p$")
    ax.set_ylabel(r"$2K_z/(K_x+K_y)$")
    ax.set_title("Anisotropy of global velocity second moments (includes resolved flows)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.savefig(root / "anisotropy.png", dpi=160)
    plt.close(fig)
    print(json.dumps(report, indent=2))
    if args.require_growth:
        failed = [name for name, value in summary.items()
                  if not (value["gamma11_in_original_nocorr_band"] and value["gamma_wb_in_original_nocorr_band"])]
        if failed:
            print("Original linear Weibel growth checks failed: " + ", ".join(failed), file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

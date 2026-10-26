#!/usr/bin/env python3
"""Check a reduced ECCAPFIM Weibel run against the kinetic dispersion relation.

The fit window (40..140 / omega_p) is shorter than in the 500 ppc, five-cell
reference: the larger initial PIC noise reaches saturation earlier here.
Fits, their statistical error and R^2 are all retained in the report.
"""
import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from snapshot_files import instantaneous_snapshots
from weibel_check import gamma_theory

FIT_WINDOW = (40.0, 140.0)
PHASE_WINDOW = (126.0, 270.0)  # Same late-stage window as the reference Weibel checker.
MODES = ((1, 1), (1, 0), (0, 1), (2, 0), (0, 2))


def read_plane(path, nx, ny):
    raw = np.fromfile(path, dtype=np.float32)
    if raw.size < 2 or not np.all(np.isfinite(raw)):
        raise ValueError(f"Invalid or nonfinite field slice: {path}")
    sx, sy = map(int, raw[:2])
    if raw.size != 2 + 3 * sx * sy:
        raise ValueError(f"Truncated field slice: {path}")
    gx, gy = (sx - nx - 1) // 2, (sy - ny - 1) // 2
    if gx < 0 or gy < 0 or sx != nx + 2 * gx + 1 or sy != ny + 2 * gy + 1:
        raise ValueError(f"Inconsistent physical/ghost dimensions in {path}")
    # Beren3D writes both ghost layers and the duplicate periodic end node.
    return raw[2:].reshape(3, sx, sy)[:, gx:gx + nx, gy:gy + ny]


def fit_growth(times, amplitude):
    mask = (times >= FIT_WINDOW[0]) & (times <= FIT_WINDOW[1])
    if mask.sum() < 8 or np.any(amplitude[mask] <= 0):
        raise ValueError("Insufficient positive mode amplitudes in the prescribed fit window")
    x, y = times[mask], np.log(amplitude[mask])
    coeff, cov = np.polyfit(x, y, 1, cov=True)
    residual = y - np.polyval(coeff, x)
    return {
        "gamma": float(coeff[0]), "stderr": float(np.sqrt(cov[0, 0])),
        "r_squared": float(1 - residual.dot(residual) / ((y - y.mean()).dot(y - y.mean()))),
        "first_time": float(x[0]), "last_time": float(x[-1]), "samples": int(mask.sum()),
    }


def analyze(workdir):
    cfg = json.loads((workdir / "system_config.json").read_text())
    pcfg = json.loads((workdir / "particles_config.json").read_text())
    if cfg["Scheme"] != "eccapfim":
        raise ValueError("This checker expects Scheme=eccapfim")
    particles = pcfg["particles"]
    if (len(particles) != 1 or particles[0]["Mass"] != 1 or particles[0]["Charge"] != -1
            or particles[0]["Density"] != 1
            or particles[0]["distribution"][0]["dist_pulse"]["sigma"] != [0.1, 0.1, 1.0]):
        raise ValueError("Dispersion comparison requires the reference 0.1/0.1/1 keV electrons")
    rows = np.atleast_1d(np.genfromtxt(workdir / "energy.txt", names=True))
    if rows.size < 2 or any(not np.all(np.isfinite(rows[n])) for n in rows.dtype.names):
        raise ValueError("Empty or nonfinite energy diagnostics")
    expected_time = cfg["LastTimestep"] * cfg["Dt"]
    total = rows["Electrons"] + rows["energyFieldE"] + rows["energyFieldB"]
    initial_energy = float(total[0])
    if initial_energy <= 0:
        raise ValueError("Nonpositive initial plasma energy")
    outputs = cfg["diagnostics"]["outputs"]
    field_cfg = next(o for o in outputs if o["type"] == "fields_2d")
    directory = workdir / "Fields" / "Diag2D"
    paths = sorted(directory.glob("FieldB_PlaneZ_*_*"))
    prefixes = {p.name.rsplit("_", 1)[0] + "_" for p in paths if "avg" not in p.name}
    if len(prefixes) != 1:
        raise ValueError("Expected exactly one instantaneous Z plane")
    prefix = prefixes.pop()
    snapshots = instantaneous_snapshots(directory, prefix)
    times = np.array(snapshots) * field_cfg["interval"] * cfg["Dt"]
    nx, ny = cfg["NumCellsX"], cfg["NumCellsY"]
    modes = {mode: [] for mode in MODES}
    coefficients = []
    bz_ratios = []
    for snapshot in snapshots:
        plane = read_plane(directory / f"{prefix}{snapshot:04d}", nx, ny)
        ft = np.fft.fft2(plane) / (nx * ny)
        for mode in MODES:
            modes[mode].append(float(np.linalg.norm(ft[:2, mode[0], mode[1]])))
        coefficients.append(ft[:, 1, 1])
        transverse = np.sqrt(np.mean(plane[0] ** 2 + plane[1] ** 2))
        bz_ratios.append(float(np.sqrt(np.mean(plane[2] ** 2)) / transverse) if transverse > 0 else 0)
    modes = {mode: np.array(values) for mode, values in modes.items()}
    fits = {}
    for mode in MODES:
        fit = fit_growth(times, modes[mode])
        k = 2 * np.pi * math.hypot(mode[0] / (nx * cfg["Dx"]), mode[1] / (ny * cfg["Dy"]))
        fit.update(k=k, theory=float(gamma_theory(k)))
        fits[str(mode)] = fit
    fit = fits[str((1, 1))]
    coefficient = np.array(coefficients)
    # Bx and By occupy different Yee faces; remove their spatial phase shifts.
    kx_dx, ky_dy = 2 * np.pi / nx, 2 * np.pi / ny
    bx = coefficient[:, 0] * np.exp(-0.5j * ky_dy)
    by = coefficient[:, 1] * np.exp(-0.5j * kx_dx)
    mask = (times >= FIT_WINDOW[0]) & (times <= FIT_WINDOW[1])
    projected = (-bx + by) / np.sqrt(2)
    early_phase_span = float(np.ptp(np.unwrap(np.angle(projected[mask]))))
    phase_mask = (times >= PHASE_WINDOW[0]) & (times <= PHASE_WINDOW[1])
    if phase_mask.sum() < 5:
        raise ValueError("Insufficient late-stage samples for the aperiodicity check")
    # Initial transverse plasma waves interfere with the unstable mode before
    # it dominates; aperiodicity is measured in the reference late-stage window.
    phase_span = float(np.ptp(np.unwrap(np.angle(projected[phase_mask]))))
    ratio = bx[mask] / by[mask]
    polarization_error = float(np.median(np.abs(ratio + 1)))
    relative_drift = float(np.max(np.abs(total - total[0])) / total[0])
    peak = int(np.argmax(rows["energyFieldB"]))
    anisotropy = 2 * rows["ElectronsZ"] / rows["ElectronsXY"]
    maxwell_residual = float(np.max(rows["nonlinearResidual"]))
    checks = {
        "completed_requested_time": bool(np.isclose(rows["Time"][-1], expected_time)),
        "constant_particle_count": bool(np.all(rows["ElectronsParticles"] == rows["ElectronsParticles"][0])),
        "relative_total_energy_drift_below_1e-8": relative_drift < 1e-8,
        "step_energy_balance_below_1e-9": bool(np.max(rows["energyConserve"]) / initial_energy < 1e-9),
        "nonlinear_residual_below_1e-8": maxwell_residual < 1e-8,
        # 20% includes statistical noise and mesh/time errors for this reduced run.
        # Report both the measured uncertainty and all other modes, without fitting a new window.
        "mode_11_growth_within_20_percent_of_theory": abs(fit["gamma"] / fit["theory"] - 1) < 0.2,
        "mode_11_growth_significant_at_5_sigma": fit["gamma"] > 5 * fit["stderr"],
        "transverse_mode_11_polarization": polarization_error < 0.2,
        "aperiodic_mode_11": phase_span < 0.6,
        "magnetic_energy_above_1_percent_of_initial_kinetic": float(rows["energyFieldB"][peak]) / initial_energy > 0.01,
        "temperature_anisotropy_relaxes": bool(anisotropy[-1] < 0.8 * anisotropy[0]),
        "in_plane_magnetic_field_dominates_at_end": bz_ratios[-1] < 0.3,
    }
    report = {
        "workdir": str(workdir), "grid": [nx, ny, cfg["NumCellsZ"]],
        "particles": int(rows["ElectronsParticles"][0]), "dt": cfg["Dt"],
        "final_time": float(rows["Time"][-1]), "fit_window": FIT_WINDOW,
        "fits": fits, "phase_window": PHASE_WINDOW,
        "phase_samples": int(phase_mask.sum()), "phase_span": phase_span,
        "early_phase_span": early_phase_span, "polarization_error": polarization_error,
        "relative_energy_drift": relative_drift,
        "maximum_step_energy_residual": float(np.max(rows["energyConserve"])),
        "maximum_nonlinear_residual": maxwell_residual,
        "maximum_nonlinear_iterations": int(np.max(rows["nonlinearIterations"])),
        "magnetic_peak": float(rows["energyFieldB"][peak]), "magnetic_peak_time": float(rows["Time"][peak]),
        "anisotropy_initial": float(anisotropy[0]), "anisotropy_final": float(anisotropy[-1]),
        "bz_to_bxy_final": bz_ratios[-1], "checks": checks, "passed": all(checks.values()),
    }
    return report, rows, times, modes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdir", type=Path)
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    try:
        report, rows, times, modes = analyze(args.workdir.resolve())
        directory = args.outdir or args.workdir / "Analysis"
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "weibel_eccapfim_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
        axes[0, 0].semilogy(rows["Time"], rows["energyFieldB"], label="magnetic")
        axes[0, 0].semilogy(rows["Time"], rows["energyFieldE"], label="electric")
        axes[0, 0].set_ylabel("Field energy"); axes[0, 0].legend()
        for mode in MODES:
            axes[0, 1].semilogy(times, modes[mode], label=str(mode))
        fit = report["fits"][str((1, 1))]
        t = np.linspace(*FIT_WINDOW, 60)
        amplitude = np.interp(t[0], times, modes[(1, 1)])
        axes[0, 1].semilogy(t, amplitude * np.exp(fit["theory"] * (t-t[0])), "k--", label="theory (1,1)")
        axes[0, 1].set_ylabel("Magnetic mode amplitude"); axes[0, 1].legend(fontsize=8)
        total = rows["Electrons"] + rows["energyFieldE"] + rows["energyFieldB"]
        axes[1, 0].plot(rows["Time"], (total-total[0])/total[0])
        axes[1, 0].set_ylabel("Relative total energy change")
        axes[1, 1].plot(rows["Time"], 2*rows["ElectronsZ"]/rows["ElectronsXY"])
        axes[1, 1].set_ylabel("T_parallel / T_perp")
        for ax in axes.flat:
            ax.set_xlabel("omega_p t"); ax.grid(alpha=0.3)
        fig.savefig(directory / "weibel_eccapfim.png", dpi=160)
        plt.close(fig)
        print(json.dumps(report, indent=2, allow_nan=False))
        return 0 if report["passed"] else 1
    except (ValueError, OSError, KeyError, StopIteration) as error:
        print(f"Weibel validation failed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

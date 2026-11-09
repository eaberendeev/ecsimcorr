#!/usr/bin/env python3
"""Run a bounded, deterministic 1D3V Weibel benchmark using the production PIC step."""
import argparse
import itertools
import json
import math
import os
from pathlib import Path
import subprocess

import numpy as np
from scipy.optimize import brentq
from scipy.special import erfcx, ndtri


def gamma_theory(k, sigma_cold, anisotropy):
    """Kinetic transverse dispersion, c=omega_pe=1, sigma_cold=sqrt(T_x/m).

    For omega=i*gamma: gamma^2+k^2+1-A+A*sqrt(pi)*eta*erfcx(eta)=0,
    eta=gamma/(sqrt(2)*k*sigma_cold), A=T_y/T_x. Ion response is O(1e-12).
    See Kuldinow & Hara (2025), Eq. (2.2), doi:10.1017/S0022377825000303.
    """
    if anisotropy <= 1 + k * k:
        return 0.0

    def dispersion(gamma):
        eta = gamma / (math.sqrt(2) * k * sigma_cold)
        return gamma**2 + k*k + 1 - anisotropy + anisotropy * math.sqrt(math.pi) * eta * erfcx(eta)

    return brentq(dispersion, 0, 10, xtol=1e-14)


def normal_nodes(count):
    # Equal macro weights: midpoint normal quantiles, with exactly unit variance.
    # Opposite quantiles occur in pairs, so the initial current is zero.
    values = ndtri((np.arange(count) + 0.5) / count)
    return values / np.sqrt(np.mean(values**2))


def prepare(args):
    sigma_cold = args.sigma_hot / math.sqrt(args.anisotropy)
    gamma = gamma_theory(args.k, sigma_cold, args.anisotropy)
    nx = args.cells
    dx = 2 * math.pi / args.k / nx
    steps = round(args.end_time / args.dt)
    if abs(steps * args.dt - args.end_time) > 1e-8:
        raise ValueError("end-time must be a whole number of time steps")
    if args.dt * math.sqrt(dx**-2 + 2) >= 1:
        raise ValueError("requested Dt violates the conservative Maxwell CFL")
    velocities = list(itertools.product(sigma_cold * normal_nodes(args.cold_nodes),
                                       args.sigma_hot * normal_nodes(args.hot_nodes),
                                       sigma_cold * normal_nodes(2)))
    system = {
        "Scheme": "symplectic", "Dt": args.dt, "Dx": dx, "Dy": 1, "Dz": 1,
        "NumCellsX": nx, "NumCellsY": 1, "NumCellsZ": 1,
        "StartFromTime": 0, "StartTimeStep": 0, "LastTimestep": steps,
        "n0": 1e10, "k_particles_reservation": 0, "Collisions": [],
        "Boundary_conditions": [{"periodic": {"face": face}} for face in ("XMIN", "YMIN", "ZMIN")],
        "diagnostics": {"outputs": [{"type": "energy_balance", "interval": max(1, round(1 / args.dt))}]},
    }
    particles = {"particles": [
        {"Name": "Electrons", "Charge": -1, "Mass": 1, "Density": 1,
         "NumPartPerCell": len(velocities), "distribution": []},
        {"Name": "Ions", "Charge": 1, "Mass": 1e12, "Density": 1, "NumPartPerCell": 1, "distribution": []},
    ]}
    config = {
        "system": system, "particles": particles, "velocities": [list(v) for v in velocities],
        "k": args.k, "seed_b": args.seed_b, "seed_gamma": gamma,
        "gamma_theory": gamma, "anisotropy": args.anisotropy,
        "sigma_cold": sigma_cold, "sigma_hot": args.sigma_hot,
        "cold_nodes": args.cold_nodes, "hot_nodes": args.hot_nodes,
        "fit_window": [args.fit_start, args.fit_end], "sample_interval": max(1, round(1 / args.dt)),
    }
    args.directory.mkdir(parents=True, exist_ok=False)
    (args.directory / "weibel.json").write_text(json.dumps(config, indent=2) + "\n")
    return config


def analyze(directory, make_plot=True):
    config = json.loads((directory / "weibel.json").read_text())
    data = np.atleast_1d(np.genfromtxt(directory / "weibel.csv", names=True, delimiter=","))
    if not len(data) or not all(np.isfinite(data[name]).all() for name in data.dtype.names):
        raise ValueError("empty or non-finite Weibel time series")
    t = data["time"]
    amplitude = np.hypot(data["mode_b_real"], data["mode_b_imag"])
    complete = int(data["step"][-1]) == config["system"]["LastTimestep"]
    checks = {
        "completed": complete,
        "energy": bool(np.max(np.abs(data["relative_energy_error"])) < 0.005),
        "gauss": bool(np.max(data["gauss_error"]) < 1e-10),
        "div_b": bool(np.max(data["div_b"]) < 1e-12),
    }
    report = {
        "configuration": {key: config[key] for key in ("k", "anisotropy", "sigma_cold", "sigma_hot", "gamma_theory")},
        "dt": config["system"]["Dt"], "cells": config["system"]["NumCellsX"],
        "particles": len(config["velocities"]) * config["system"]["NumCellsX"],
        "last_time": float(t[-1]), "wall_seconds": float(data["wall_seconds"][-1]),
        "max_relative_energy_error": float(np.max(np.abs(data["relative_energy_error"]))),
        "max_gauss_error": float(np.max(data["gauss_error"])),
        "max_div_b": float(np.max(data["div_b"])),
        "max_amplification": float(amplitude.max() / amplitude[0]),
        "initial_second_moment_ratio": float(data["vy2"][0] / data["vx2"][0]),
        "final_second_moment_ratio": float(data["vy2"][-1] / data["vx2"][-1]),
    }
    if config["gamma_theory"] > 0:
        low, high = config["fit_window"]
        mask = (t >= low) & (t <= high) & (amplitude > 0)
        if mask.sum() < 10:
            raise ValueError("at least ten samples required in the fixed linear fit window")
        x, y = t[mask], np.log(amplitude[mask])
        (growth, intercept), covariance = np.polyfit(x, y, 1, cov=True)
        r2 = 1 - np.sum((y - growth*x - intercept)**2) / np.sum((y - y.mean())**2)
        report.update(gamma_fit=float(growth), gamma_fit_standard_error=float(np.sqrt(covariance[0, 0])),
                      gamma_relative_error=float(abs(growth / config["gamma_theory"] - 1)), fit_r2=float(r2),
                      fit_window=[low, high])
        checks.update(growth=report["gamma_relative_error"] < 0.10, exponential=r2 > 0.98,
                      amplification=report["max_amplification"] > 30)
    else:
        checks["stable_control"] = report["max_amplification"] < 3
    checks = {name: bool(value) for name, value in checks.items()}
    report["checks"] = checks
    report["passed"] = bool(all(checks.values()))
    (directory / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    if make_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
        axes[0, 0].semilogy(t, amplitude, label="PIC |B_k|")
        if config["gamma_theory"] > 0:
            anchor = np.interp(config["fit_window"][0], t, amplitude)
            axes[0, 0].semilogy(t, anchor*np.exp(config["gamma_theory"]*(t-config["fit_window"][0])),
                                "--", label="kinetic theory (anchored)")
            axes[0, 0].axvspan(*config["fit_window"], alpha=0.15)
        axes[0, 0].legend(); axes[0, 0].set_ylabel("Magnetic mode amplitude")
        axes[0, 1].semilogy(t, data["magnetic_energy"], label="magnetic")
        axes[0, 1].semilogy(t, data["electric_energy"], label="electric")
        axes[0, 1].legend(); axes[0, 1].set_ylabel("Field energy")
        axes[1, 0].plot(t, data["relative_energy_error"])
        axes[1, 0].set_ylabel("(H-H0)/H0")
        axes[1, 1].plot(t, data["vy2"] / data["vx2"])
        axes[1, 1].set_ylabel("<v_y²>/<v_x²>")
        for ax in axes.flat:
            ax.set_xlabel("t omega_pe"); ax.grid(alpha=0.25)
        fig.suptitle(f"Symplectic 1D3V Weibel: A={config['anisotropy']:g}, dt={config['system']['Dt']:g}")
        fig.savefig(directory / "weibel.png", dpi=160)
        plt.close(fig)
    return report


def status(directory):
    config = json.loads((directory / "weibel.json").read_text())
    text = (directory / "weibel.csv").read_text()
    # Ignore a final row that is still being written by the simulation.
    lines = text.splitlines()
    if not text.endswith("\n"):
        lines = lines[:-1]
    if len(lines) < 2:
        return {"state": "initializing"}
    row = dict(zip(lines[0].split(","), map(float, lines[-1].split(","))))
    fraction = row["step"] / config["system"]["LastTimestep"]
    return {"state": "completed" if fraction == 1 else "partial", "step": int(row["step"]), "last_step": config["system"]["LastTimestep"],
            "time": row["time"], "fraction_complete": fraction,
            "mode_amplitude": math.hypot(row["mode_b_real"], row["mode_b_imag"]),
            "relative_energy_error": row["relative_energy_error"], "gauss_error": row["gauss_error"],
            "div_b": row["div_b"], "wall_seconds": row["wall_seconds"],
            "estimated_seconds_remaining": row["wall_seconds"]*(1/fraction-1) if fraction else None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    parser.add_argument("--status", action="store_true", help="read progress from a running calculation")
    parser.add_argument("--cells", type=int, default=32)
    parser.add_argument("--cold-nodes", type=int, default=32)
    parser.add_argument("--hot-nodes", type=int, default=8)
    parser.add_argument("--k", type=float, default=1.16)
    parser.add_argument("--anisotropy", type=float, default=10)
    parser.add_argument("--sigma-hot", type=float, default=0.1)
    parser.add_argument("--dt", type=float, default=0.1)
    parser.add_argument("--end-time", type=float, default=160)
    parser.add_argument("--seed-b", type=float, default=1e-5)
    parser.add_argument("--fit-start", type=float, default=60)
    parser.add_argument("--fit-end", type=float, default=110)
    parser.add_argument("--timeout", type=float, default=900)
    args = parser.parse_args()
    args.directory = args.directory.resolve()
    if args.status:
        print(json.dumps(status(args.directory), indent=2))
        return
    for key in ("cells", "cold_nodes", "hot_nodes", "k", "anisotropy", "sigma_hot", "dt", "end_time", "seed_b", "timeout"):
        if not math.isfinite(getattr(args, key)) or getattr(args, key) <= 0:
            parser.error(f"{key} must be finite and positive")
    if args.cold_nodes % 2 or args.hot_nodes % 2:
        parser.error("velocity-node counts must be even for a zero-current quiet start")
    if not args.analyze_only:
        config = prepare(args)
        print(f"Prepared {args.directory}: {config['system']['LastTimestep']} steps, "
              f"{len(config['velocities'])*args.cells} electrons, gamma={config['gamma_theory']:.8g}", flush=True)
        if args.prepare_only:
            return
        if args.binary is None:
            parser.error("--binary is required unless preparing or analyzing only")
        env = dict(os.environ, OMP_NUM_THREADS="2")
        with (args.directory / "run.log").open("w") as log:
            subprocess.run([str(args.binary.resolve())], cwd=args.directory, stdout=log, stderr=subprocess.STDOUT,
                           env=env, timeout=args.timeout, check=True)
    report = analyze(args.directory)
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if report["passed"] else 1)


if __name__ == "__main__":
    main()

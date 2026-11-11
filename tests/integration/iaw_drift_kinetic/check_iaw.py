#!/usr/bin/env python3
"""Checks the ion-acoustic ringdown against kinetic linear theory.

The ion density is sampled by probes at the cell centres of the column (the
grid E_z carries the marker noise of the fast electrons, the ion density much
less). Its mode-1 Fourier coefficient c(t) is fitted, after the first tenth
of the run, to
    c(t) = exp(-gamma t) (P cos(omega t) + Q sin(omega t)),
P, Q complex, and omega is compared with the root of the kinetic dispersion
relation for the species of particles_config.json (iaw_theory.py). Energy
must be balanced to the nonlinear tolerance.

gamma is reported but not checked: with 512 markers per cell the marker noise
of the mode (|c| ~ 0.005-0.01) is comparable with the signal over a damping
time, and at dn/n = 0.03 the resonant electrons are trapped after
1 / (k sqrt(e phi / m_e)) ~ 300 << 1 / gamma. Measured early decay: 0.55 of
the amplitude at t = 800 for dn/n = 0.03, 0.34 for 0.01, theory 0.68. A
damping check needs the xpic statistics (8192 markers per cell, dt = 5).
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parent))
from iaw_theory import ion_acoustic_root  # noqa: E402

OMEGA_TOLERANCE = 0.03
BALANCE_TOLERANCE = 1e-8


def fail(message):
    print(f"iaw_drift_kinetic: {message}", file=sys.stderr)
    return 2


def fit(times, coefficients):
    def basis(omega, gamma):
        decay = np.exp(-gamma * times)
        return np.stack([decay * np.cos(omega * times), decay * np.sin(omega * times)], axis=1)

    def residual(x):
        a = basis(x[0], x[1])
        amplitudes, *_ = np.linalg.lstsq(a, coefficients, rcond=None)
        r = a @ amplitudes - coefficients
        return np.concatenate([r.real, r.imag])

    return residual


def main():
    workdir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    try:
        config = json.loads((workdir / "system_config.json").read_text())
        particles = json.loads((workdir / "particles_config.json").read_text())["particles"]
        lines = (workdir / "Probes" / "probes.txt").read_text().splitlines()
        energy = (workdir / "energy.txt").read_text().splitlines()
    except (OSError, KeyError, ValueError) as error:
        return fail(f"cannot read run data: {error}")

    dt, dz, nz = config["Dt"], config["Dz"], config["NumCellsZ"]
    k = 2 * math.pi / (nz * dz)
    species = [(p["Charge"], p["Mass"], p["Density"], p["distribution"][0]["dist_pulse"]["sigma"][2])
               for p in particles]
    omega_theory, gamma_theory = ion_acoustic_root(k, species, dz)

    header = next(line for line in lines if line.startswith("#")).lstrip("#").split()
    column = {name: index for index, name in enumerate(header)}
    rows = np.array([[float(v) for v in line.split()] for line in lines if line.strip() and not line.startswith("#")])
    if rows.size == 0 or not np.all(np.isfinite(rows)):
        return fail("empty or non-finite probe data")
    steps = np.unique(rows[:, column["Time"]])
    times, coefficients = [], []
    for step in steps:
        sample = rows[rows[:, column["Time"]] == step]
        z, density = sample[:, column["z"]], sample[:, column["dens_Ions"]]
        coefficients.append(2.0 / len(z) * np.sum(density * np.exp(-1j * k * z)))
        times.append(step * dt)
    times, coefficients = np.array(times), np.array(coefficients)
    keep = times >= 0.1 * times[-1]
    solution = least_squares(fit(times[keep], coefficients[keep]), [omega_theory, gamma_theory],
                             x_scale=[omega_theory, gamma_theory])
    omega, gamma = solution.x

    head = energy[0].split()
    table = np.array([[float(v) for v in line.split()] for line in energy[1:] if line.strip()])
    kinetic = sum(table[:, head.index(p["Name"])] for p in particles if p["Name"] in head)
    balance = np.max(np.abs(table[:, head.index("energyConserve")]) / np.maximum(kinetic, 1e-300))

    result = {
        "omega": omega, "omega_theory": omega_theory, "omega_error": abs(omega / omega_theory - 1),
        "gamma": gamma, "gamma_theory": gamma_theory, "gamma_error": abs(gamma / gamma_theory - 1),
        "max_relative_energy_balance": balance, "samples": int(keep.sum()),
    }
    result["gamma_checked"] = False
    result["passed"] = bool(result["omega_error"] < OMEGA_TOLERANCE and balance < BALANCE_TOLERANCE)
    print(json.dumps(result, indent=2))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())

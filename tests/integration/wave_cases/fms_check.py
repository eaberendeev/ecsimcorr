"""Analysis of a fast magnetosonic wave run, common to all schemes.

The self field dB_z is sampled by probes along x; its mode-1 cosine
coefficient y(t) is fitted to A + C cos(omega t + phi) (omega scanned, then
refined), as in xpic, and compared with the kinetic root mapped to the grid
and to the implicit midpoint rule (fms_theory.py; for an explicit scheme with
omega dt << 1 the mapping is the identity). Schemes that conserve energy
exactly must balance it.
"""

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fms_theory import Parameters, omega_discrete, omega_fluid, omega_kinetic  # noqa: E402

MEC2_KEV = 511.0
OMEGA_TOLERANCE = 0.02
BALANCE_TOLERANCE = 1e-8


def fit(t, y, omega0):
    def residual(omega):
        m = np.column_stack([np.ones_like(t), np.cos(omega * t), np.sin(omega * t)])
        coefficients, *_ = np.linalg.lstsq(m, y, rcond=None)
        return np.sum((m @ coefficients - y) ** 2)

    omegas = np.linspace(0.7 * omega0, 1.3 * omega0, 601)
    omega = omegas[np.argmin([residual(w) for w in omegas])]
    lo, hi = omega - 1e-3 * omega0, omega + 1e-3 * omega0
    for _ in range(60):   # golden-section refinement
        a, b = lo + 0.382 * (hi - lo), lo + 0.618 * (hi - lo)
        if residual(a) < residual(b):
            hi = b
        else:
            lo = a
    return 0.5 * (lo + hi)


def analyse(workdir):
    """Frequency fit and energy balance of a run; returns the result dictionary or an error string."""
    try:
        config = json.loads((workdir / "system_config.json").read_text())
        particles = {p["Name"]: p for p in json.loads((workdir / "particles_config.json").read_text())["particles"]}
        lines = (workdir / "Probes" / "probes.txt").read_text().splitlines()
        energy = (workdir / "energy.txt").read_text().splitlines()
    except (OSError, KeyError, ValueError) as error:
        return f"cannot read run data: {error}"

    dt, dx, nx = config["Dt"], config["Dx"], config["NumCellsX"]
    k = 2 * math.pi / (nx * dx)
    # Uniform background field: of the drift-kinetic schemes, the external one of
    # the kinetic schemes, or the initial flux of the symplectic scheme.
    if "DriftKinetic" in config:
        B0 = config["DriftKinetic"]["magnetic"][2]
    elif "ExternalFieldB" in config:
        B0 = next(entry["uniform_field"]["value"][2] for entry in config["ExternalFieldB"]
                  if "uniform_field" in entry)
    else:
        B0 = config["Symplectic"]["initial_magnetic"][2]
    ions, electrons = particles["Ions"], particles["Electrons"]
    p = Parameters(B0, 1.0, ions["Mass"], ions["distribution"][0]["dist_pulse"]["sigma"][0] / MEC2_KEV,
                   electrons["distribution"][0]["dist_pulse"]["sigma"][0] / MEC2_KEV, k, dx, dt)
    omega_theory = omega_discrete(p, omega_kinetic(p, k))

    header = next(line for line in lines if line.startswith("#")).lstrip("#").split()
    column = {name: index for index, name in enumerate(header)}
    rows = np.array([[float(v) for v in line.split()] for line in lines if line.strip() and not line.startswith("#")])
    if rows.size == 0 or not np.all(np.isfinite(rows)):
        return "empty or non-finite probe data"
    times, mode = [], []
    for step in np.unique(rows[:, column["Time"]]):
        sample = rows[rows[:, column["Time"]] == step]
        x, bz = sample[:, column["x"]], sample[:, column["Bnz"]]
        mode.append(2.0 / len(x) * np.sum(bz * np.cos(k * x)) / B0)
        times.append(step * dt)
    times, mode = np.array(times), np.array(mode)
    omega = fit(times, mode, omega_theory)

    head = energy[0].split()
    table = np.array([[float(v) for v in line.split()] for line in energy[1:] if line.strip()])
    kinetic = sum(table[:, head.index(name)] for name in particles if name in head)
    # Energy-conserving schemes report the balance of a step, the others the drift from the start.
    conserving = "energyConserve" in head
    defect = table[:, head.index("energyConserve" if conserving else "energyError")]
    balance = np.max(np.abs(defect) / np.maximum(kinetic, 1e-300))

    result = {"omega": omega, "omega_theory": omega_theory, "omega_fluid": omega_fluid(p, k),
              "omega_error": abs(omega / omega_theory - 1), "initial_mode": mode[0],
              "max_relative_energy_balance": balance, "periods": times[-1] * omega_theory / (2 * math.pi)}
    result["passed"] = bool(result["omega_error"] < OMEGA_TOLERANCE
                            and (balance < BALANCE_TOLERANCE or not conserving))
    return result




def main(name):
    result = analyse(Path(sys.argv[1]) if len(sys.argv) > 1 else Path("."))
    if isinstance(result, str):
        print(f"{name}: {result}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2))
    return 0 if result["passed"] else 1

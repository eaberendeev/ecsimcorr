"""Langmuir wave: Bohm-Gross frequency and Landau damping, for any kinetic scheme.

A periodic column of 64 x 1 x 1 cells, dx = 0.1 c/omega_pe, electrons with
v_te = 0.407 c (84.8 keV) and immobile ions. The wavelength is the box,
k = 2 pi / 6.4, so k lambda_De = 0.4: the kinetic root is omega = 1.285,
gamma = 0.0661 (in omega_pe), a decay of 28 % per period. The large thermal
speed is chosen for the numerics only: it lets one time step dt = 0.04 serve
every scheme, with c dt / dx = 0.4 below the explicit Courant limit and
v dt < dx for the ECSIM cell update up to 6 v_te. The wave is electrostatic
and the particle motion is non-relativistic in the code, so c does not enter.

The electrons get a velocity perturbation u = u0 sin(k x) at uniform density
(quiet start, antithetic pairs), so E = 0 satisfies Gauss's law at t = 0 and a
standing wave E_x ~ cos(k x) sin(omega t) exp(-gamma t) develops. Its
amplitude keeps the bounce frequency sqrt(k E) = 0.05 below gamma: the
damping is linear.

Landau damping is the work of the electrons near v = omega / k = 3.2 v_te,
0.2 % of the distribution: the damping rate is as good as their statistics.
With 8192 markers per cell about a thousand markers are resonant. (At
k lambda_De = 0.3 the resonance is at 3.9 v_te; 2048 markers per cell leave
30 resonant markers there, and the measured rate was 1.7 times the theory
for both ECSIM schemes alike.)

`check` fits the complex mode-1 coefficient of E_x from the probes to
exp(-gamma t) (P cos(omega t) + Q sin(omega t)) and compares omega and gamma
with the root of the kinetic dispersion relation
    k^2 / sinc(k dx / 2)^4 + sum_s (omega_ps^2 / v_ts^2) (1 + zeta Z(zeta)) = 0,
zeta = omega / (sqrt(2) k v_ts), Z the plasma dispersion function; the sinc
factor is that of linear gather and deposit with the finite-difference
divergence (0.2 % here). The time discretization shifts omega by
(omega dt)^2 / 12 = 2e-4, below the tolerance. Schemes that conserve energy
exactly must also balance it; for the others the balance is only reported.
"""

import json
import math
import sys
from pathlib import Path

MEC2_KEV = 511.0
Cells, CellSize, Dt = (64, 1, 1), 0.1, 0.04
ThermalSpeed = 0.4 * Cells[0] * CellSize / (2 * math.pi)   # k lambda_De = 0.4
TemperatureKev = MEC2_KEV * ThermalSpeed**2
IonMass, IonTemperatureKev = 1.0e6, 1.0e-6
FlowAmplitude = 0.006 * ThermalSpeed
ParticlesPerCell = 8192
Periods = 4.0
# The fit over periods 1-4 moves by 0.5 % in omega and 6 % in gamma with the
# window; past four periods the mode falls to the marker noise of its neighbours.
OMEGA_TOLERANCE, GAMMA_TOLERANCE, BALANCE_TOLERANCE = 0.02, 0.20, 1e-6
ENERGY_CONSERVING = ("ecsim", "ecsim_corr", "eccapfim")


def dispersion(omega, k, species, dx):
    """Electrostatic dielectric response times k^2; species = (charge, mass, density, T in keV)."""
    from scipy.special import wofz

    sinc = math.sin(0.5 * k * dx) / (0.5 * k * dx)
    total = k * k / sinc**4
    for charge, mass, density, temperature_kev in species:
        vt = math.sqrt(temperature_kev / MEC2_KEV / mass)
        zeta = omega / (math.sqrt(2.0) * k * vt)
        total += density * charge**2 / mass / vt**2 * (1.0 + zeta * 1j * math.sqrt(math.pi) * wofz(zeta))
    return total


def theory():
    """(omega, gamma) of the Langmuir root, gamma > 0 for damping."""
    from scipy.optimize import fsolve

    k = 2 * math.pi / (Cells[0] * CellSize)
    species = [(-1.0, 1.0, 1.0, TemperatureKev), (1.0, IonMass, 1.0, IonTemperatureKev)]
    guess = math.sqrt(1 + 3 * (k * ThermalSpeed) ** 2)

    def residual(x):
        value = dispersion(x[0] + 1j * x[1], k, species, CellSize)
        return [value.real, value.imag]

    (re, im), _, converged, message = fsolve(residual, [guess, -0.05], xtol=1e-10, full_output=True)
    # The solver stalls at rounding level of the residual, which is of the order k^2 ~ 1.
    if converged != 1 and max(abs(v) for v in residual([re, im])) > 1e-6:
        raise RuntimeError(f"Langmuir dispersion root not found: {message}")
    return re, -im


def configuration(scheme):
    omega, _ = theory()
    box = [Cells[a] * CellSize for a in range(3)]
    system = {
        "Scheme": scheme,
        "Dt": Dt,
        "Dx": CellSize,
        "Dy": CellSize,
        "Dz": CellSize,
        "NumCellsX": Cells[0],
        "NumCellsY": Cells[1],
        "NumCellsZ": Cells[2],
        "StartFromTime": 0,
        "StartTimeStep": 0,
        "LastTimestep": round(Periods * 2 * math.pi / omega / Dt),
        "n0": 1e10,
        "k_particles_reservation": 0,
        "Collisions": [],
        "DampingType": "None",
        "Boundary_conditions": [{"periodic": {"face": face}} for face in ("XMIN", "YMIN", "ZMIN")],
        "diagnostics": {
            "outputs": [
                {"type": "energy_balance", "interval": 1},
                {"type": "console_summary"},
                {"type": "probes", "interval": 1, "species": [],
                 "points": [[(i + 0.5) * CellSize, 0.5 * CellSize, 0.5 * CellSize] for i in range(Cells[0])]},
            ]
        },
    }

    def species(name, charge, mass, temperature, perturbation=None):
        distribution = {
            "type": "initial",
            "sampling": "halton_quiet",
            "stratified_axis": 0,
            "dist_space": {"type": "rectangle", "center": [0.5 * b for b in box],
                           "half_length": [0.5 * b for b in box]},
            "dist_pulse": {"type": "gaussian", "mean": [0, 0, 0], "sigma": [temperature] * 3},
            "density": 1.0,
        }
        if perturbation:
            distribution["perturbation"] = perturbation
        return {"Name": name, "Charge": charge, "Mass": mass, "Density": 1.0,
                "NumPartPerCell": ParticlesPerCell, "distribution": [distribution]}

    particles = {"particles": [
        species("Ions", 1.0, IonMass, IonTemperatureKev),
        species("Electrons", -1.0, 1.0, TemperatureKev, {"amplitude": 0.0, "mode": 1, "flow": FlowAmplitude}),
    ]}
    return system, particles


def write_configuration(scheme):
    system, particles = configuration(scheme)
    for filename, config in (("system_config.json", system), ("particles_config.json", particles)):
        with open(filename, "w", encoding="utf-8") as stream:
            json.dump(config, stream, indent=2)
    wp = math.sqrt(4 * math.pi * system["n0"] * 4.80320427e-10 ** 2 / 9.10938356e-28)
    with open("phys.par", "w", encoding="utf-8") as stream:
        stream.write(f"w_p = {wp}\n1/w_p = {1 / wp}\n")


def check(workdir):
    """Prints the comparison; returns the process exit code."""
    import numpy as np
    from scipy.optimize import least_squares

    workdir = Path(workdir)
    try:
        config = json.loads((workdir / "system_config.json").read_text())
        lines = (workdir / "Probes" / "probes.txt").read_text().splitlines()
        energy = (workdir / "energy.txt").read_text().splitlines()
    except (OSError, ValueError) as error:
        print(f"langmuir: cannot read run data: {error}", file=sys.stderr)
        return 2
    dt, k = config["Dt"], 2 * math.pi / (config["NumCellsX"] * config["Dx"])
    header = next(line for line in lines if line.startswith("#")).lstrip("#").split()
    column = {name: index for index, name in enumerate(header)}
    rows = np.array([[float(v) for v in line.split()] for line in lines if line.strip() and not line.startswith("#")])
    if rows.size == 0 or not np.all(np.isfinite(rows)):
        print("langmuir: empty or non-finite probe data", file=sys.stderr)
        return 2
    times, mode = [], []
    for step in np.unique(rows[:, column["Time"]]):
        sample = rows[rows[:, column["Time"]] == step]
        x, ex = sample[:, column["x"]], sample[:, column["Enx"]]
        mode.append(2.0 / len(x) * np.sum(ex * np.exp(-1j * k * x)))
        times.append(step * dt)
    times, mode = np.array(times), np.array(mode)
    omega_theory, gamma_theory = theory()
    # The first period holds the free-streaming transient of the velocity
    # perturbation, which decays as exp(-(k v_te t)^2 / 2).
    start = np.searchsorted(times, 2 * math.pi / omega_theory)
    t, c = times[start:], mode[start:]

    def residual(x):
        decay = np.exp(-x[1] * t)
        basis = np.stack([decay * np.cos(x[0] * t), decay * np.sin(x[0] * t)], axis=1)
        amplitudes, *_ = np.linalg.lstsq(basis, c, rcond=None)
        r = basis @ amplitudes - c
        return np.concatenate([r.real, r.imag])

    solution = least_squares(residual, [omega_theory, gamma_theory], x_scale=[omega_theory, gamma_theory])
    omega, gamma = float(solution.x[0]), float(solution.x[1])

    head = energy[0].split()
    table = np.array([[float(v) for v in line.split()] for line in energy[1:] if line.strip()])
    kinetic = table[:, head.index("Electrons")]
    # Energy-conserving schemes report the balance of a step, the others the drift from the start.
    defect = "energyConserve" if "energyConserve" in head else "energyError"
    balance = float(np.max(np.abs(table[:, head.index(defect)]) / kinetic))
    result = {
        "scheme": config["Scheme"], "omega": omega, "omega_theory": omega_theory,
        "omega_error": abs(omega / omega_theory - 1), "gamma": gamma, "gamma_theory": gamma_theory,
        "gamma_error": abs(gamma / gamma_theory - 1), "max_mode_amplitude": float(np.max(np.abs(mode))),
        "max_relative_energy_balance": balance,
    }
    result["passed"] = bool(result["omega_error"] < OMEGA_TOLERANCE and result["gamma_error"] < GAMMA_TOLERANCE
                            and (balance < BALANCE_TOLERANCE or config["Scheme"] not in ENERGY_CONSERVING))
    print(json.dumps(result, indent=2))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    print("omega = %.6f, gamma = %.6f, T = %.3f keV" % (*theory(), TemperatureKev))

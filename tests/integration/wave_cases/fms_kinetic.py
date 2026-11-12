"""Fast magnetosonic wave with kinetic electrons and ions, for any kinetic scheme.

The plasma and the wave are those of fms_hybrid: a periodic column of
32 x 1 x 1 cells, dx = 10 c/omega_pe, uniform n = 1, Te = Ti = 1 keV,
m_i / m_e = 100, B = B0 z with B0 = 0.2 and an initial self-field perturbation
dB_z = 0.05 B0 cos(k x), k = 2 pi / Lx; 2.1 periods of the mode. Here both
species are kinetic and the uniform field is external (for the symplectic
scheme it is the initial magnetic flux).

The implicit schemes use Omega_e dt = 0.4: the orbit of the midpoint rule
then deposits 96 % of the electron magnetic moment
(1 / (1 + (Omega_e dt / 2)^2)); with Omega_e dt = 10 the electron diamagnetic
response is lost and the frequency is 4 % low (docs/drift-kinetic.md). The
explicit symplectic scheme needs omega_pe dt < 2 and uses dt = 0.5.

fms_check.py fits dB_z(k, t) and compares the frequency with the hot-ion
kinetic root (fms_theory.py).
"""

import json
import math

B0, Amplitude, IonMass, TemperatureKev = 0.2, 0.05, 100.0, 1.0
Cells, CellSize, ParticlesPerCell = (32, 1, 1), 10.0, 64
Periods = 2.1
TimeStep = {"symplectic": 0.5}   # explicit; the others: Omega_e dt = 0.4
GyrationStep = 0.4


def configuration(scheme):
    k = 2 * math.pi / (Cells[0] * CellSize)
    # Fluid estimate of the frequency, for the length of the run only.
    omega = k * math.sqrt((B0**2 / IonMass + 4 * TemperatureKev / 511.0 / IonMass) / (1 + B0**2 / IonMass))
    dt = TimeStep.get(scheme, GyrationStep / B0)
    steps = round(Periods * 2 * math.pi / omega / dt)
    box = [Cells[a] * CellSize for a in range(3)]
    system = {
        "Scheme": scheme,
        "Dt": dt,
        "Dx": CellSize,
        "Dy": CellSize,
        "Dz": CellSize,
        "NumCellsX": Cells[0],
        "NumCellsY": Cells[1],
        "NumCellsZ": Cells[2],
        "StartFromTime": 0,
        "StartTimeStep": 0,
        "LastTimestep": steps,
        "n0": 1e10,
        "k_particles_reservation": 0,
        "Collisions": [],
        "DampingType": "None",
        "Boundary_conditions": [{"periodic": {"face": face}} for face in ("XMIN", "YMIN", "ZMIN")],
        "InitialFieldB": [{"cosine_field": {"amplitude": [0, 0, Amplitude * B0], "mode": [1, 0, 0]}}],
        "diagnostics": {
            "outputs": [
                {"type": "energy_balance", "interval": max(1, steps // 15000)},
                {"type": "console_summary"},
                # about 3000 samples over the run; a wave period is 14600 time units
                {"type": "probes", "interval": max(1, steps // 3000), "species": [],
                 "points": [[(i + 0.5) * CellSize, 0.5 * CellSize, 0.5 * CellSize] for i in range(Cells[0])]},
            ]
        },
    }
    if scheme == "ecsim_corr":
        # The single-pass scheme is unstable here (Omega_e dt = 0.4): wave_cases/README.md.
        system["CorrectorPasses"] = 1
    if scheme == "symplectic":
        system["Symplectic"] = {"initial_magnetic": [0, 0, B0]}
    else:
        system["ExternalFieldB"] = [{"uniform_field": {"value": [0, 0, B0]}}]

    def species(name, charge, mass):
        return {
            "Name": name,
            "Charge": charge,
            "Mass": mass,
            "Density": 1.0,
            "NumPartPerCell": ParticlesPerCell,
            "distribution": [{
                "type": "initial",
                "sampling": "halton_quiet",
                "stratified_axis": 0,
                "dist_space": {"type": "rectangle", "center": [0.5 * b for b in box],
                               "half_length": [0.5 * b for b in box]},
                "dist_pulse": {"type": "gaussian", "mean": [0, 0, 0], "sigma": [TemperatureKev] * 3},
                "density": 1.0,
            }],
        }

    return system, {"particles": [species("Ions", 1.0, IonMass), species("Electrons", -1.0, 1.0)]}


def write_configuration(scheme):
    system, particles = configuration(scheme)
    for filename, config in (("system_config.json", system), ("particles_config.json", particles)):
        with open(filename, "w", encoding="utf-8") as stream:
            json.dump(config, stream, indent=2)
    wp = math.sqrt(4 * math.pi * system["n0"] * 4.80320427e-10 ** 2 / 9.10938356e-28)
    with open("phys.par", "w", encoding="utf-8") as stream:
        stream.write(f"w_p = {wp}\n1/w_p = {1 / wp}\n")

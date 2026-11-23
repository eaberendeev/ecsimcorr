#!/usr/bin/env python3
"""Fast magnetosonic wave in the hybrid model (kinetic ions, drift-kinetic electrons).

Reduced analogue of xpic hybrid_drift_kinetic_ex2_fms: a periodic column of
32 x 1 x 1 cells, dx = 10, uniform n = 1, Te = Ti = 1 keV, mi / me = 100,
B = B0 z with B0 = 0.2 and an initial self-field perturbation
dB_z = 0.05 B0 cos(k x), k = 2 pi / Lx. dt = 50 gives Omega_i dt = 0.1 and
Omega_e dt = 10; the run covers 2.1 periods of the mode. check_fms.py fits
dB_z(k, t) and compares the frequency with the hot-ion kinetic root
(fms_theory.py), mapped to the grid and to the implicit midpoint rule.
"""

import json
import math

# 4k markers x 615 steps: about 20 s on 8 cores, part of the regular integration set
NumProcs = 8
DirName = "Res_FMS_Hybrid"

B0, Amplitude, IonMass, TemperatureKev = 0.2, 0.05, 100.0, 1.0
Cells, CellSize, ParticlesPerCell = (32, 1, 1), 10.0, 64
Dt = 0.1 / (B0 / IonMass)
_k = 2 * math.pi / (Cells[0] * CellSize)
_omega = _k * math.sqrt((B0**2 / IonMass + 4 * TemperatureKev / 511.0 / IonMass) / (1 + B0**2 / IonMass))
LastTimestep = round(2.1 * 2 * math.pi / _omega / Dt)

system_config = {
    "Scheme": "hybrid_kinetic",
    "Dt": Dt,
    "Dx": CellSize,
    "Dy": CellSize,
    "Dz": CellSize,
    "NumCellsX": Cells[0],
    "NumCellsY": Cells[1],
    "NumCellsZ": Cells[2],
    "StartFromTime": 0,
    "StartTimeStep": 0,
    "LastTimestep": LastTimestep,
    "n0": 1e10,
    "k_particles_reservation": 0,
    "Collisions": [],
    "DampingType": "None",
    "Boundary_conditions": [{"periodic": {"face": face}} for face in ("XMIN", "YMIN", "ZMIN")],
    "DriftKinetic": {"magnetic": [0, 0, B0], "electric": [0, 0, 0]},
    "InitialFieldB": [{"cosine_field": {"amplitude": [0, 0, Amplitude * B0], "mode": [1, 0, 0]}}],
    "DriftKineticPush": {"coord_is_gc": True},
    "Eccapfim": {"coupled_fields": True, "preconditioner": "mass_matrix", "coupled_max_outer": 60},
    "diagnostics": {
        "outputs": [
            {"type": "energy_balance", "interval": 1},
            {"type": "console_summary"},
            {"type": "probes", "interval": 1, "species": [],
             "points": [[(i + 0.5) * CellSize, 0.5 * CellSize, 0.5 * CellSize] for i in range(Cells[0])]},
        ]
    },
}


def species(name, charge, mass, model):
    box = [Cells[a] * CellSize for a in range(3)]
    return {
        "Name": name,
        "Charge": charge,
        "Mass": mass,
        "Density": 1.0,
        "NumPartPerCell": ParticlesPerCell,
        "model": model,
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


particles_config = {"particles": [
    species("Ions", 1.0, IonMass, "kinetic"),
    species("Electrons", -1.0, 1.0, "drift_kinetic"),
]}


def generate_config():
    for filename, config in (("system_config.json", system_config),
                             ("particles_config.json", particles_config)):
        with open(filename, "w", encoding="utf-8") as stream:
            json.dump(config, stream, indent=2)
    wp = math.sqrt(4 * math.pi * system_config["n0"] * 4.80320427e-10 ** 2 / 9.10938356e-28)
    with open("phys.par", "w", encoding="utf-8") as stream:
        stream.write(f"w_p = {wp}\n1/w_p = {1 / wp}\n")
    return DirName


if __name__ == "__main__":
    generate_config()

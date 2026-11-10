#!/usr/bin/env python3
"""Reduced hybrid_ex1 regression: kinetic ions + drift-kinetic electrons.

Periodic 4^3 box, dt=10, Bz=1 (omega_ce dt = omega_pe dt = 10), 10
markers/cell/species.  This is deliberately small: it validates the coupled
E+B solve with NGMRES and the mass-matrix preconditioner, the hybrid current
and the discrete energy theorem without the cost of a production run.
"""

import json
import math

NumProcs = 4
LastTimestep = 5
DirName = "Res_Hybrid_Coupled"

system_config = {
    "Scheme": "hybrid_kinetic",
    "Dt": 10.0,
    "Dx": 1.0,
    "Dy": 1.0,
    "Dz": 1.0,
    "NumCellsX": 4,
    "NumCellsY": 4,
    "NumCellsZ": 4,
    "StartFromTime": 0,
    "StartTimeStep": 0,
    "LastTimestep": LastTimestep,
    "n0": 1e10,
    "k_particles_reservation": 0,
    "Collisions": [],
    "DampingType": "None",
    "Boundary_conditions": [
        {"periodic": {"face": face}} for face in ("XMIN", "YMIN", "ZMIN")
    ],
    "DriftKinetic": {"magnetic": [0, 0, 1], "electric": [0, 0, 0]},
    "Eccapfim": {
        "coupled_fields": True,
        "preconditioner": "mass_matrix",
        "coupled_solver": "ngmres",
        "coupled_max_outer": 30,
    },
    "diagnostics": {
        "outputs": [
            {"type": "energy_balance", "interval": 1},
            {"type": "console_summary"},
        ]
    },
}


def species(name, charge, mass, model):
    return {
        "Name": name,
        "Charge": charge,
        "Mass": mass,
        "Density": 1.0,
        "NumPartPerCell": 10,
        "model": model,
        "distribution": [{
            "type": "initial",
            "dist_space": {
                "type": "rectangle",
                "center": [2, 2, 2],
                "half_length": [1.5, 1.5, 1.5],
            },
            "dist_pulse": {
                "type": "gaussian",
                "mean": [0, 0, 0],
                "sigma": [0.1, 0.1, 0.1],
            },
            "density": 1.0,
        }],
    }


particles_config = {"particles": [
    species("Ions", 1.0, 100.0, "kinetic"),
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

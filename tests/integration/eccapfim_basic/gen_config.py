#!/usr/bin/env python3
"""Small periodic ECCAPFIM plasma; writes configs into the current directory."""

import json
import math

NumProcs = 2
LastTimestep = 8
DirName = "Res_TestInt_Eccapfim"

system_config = {
    "Scheme": "eccapfim",
    "Dt": 0.5,
    "Dx": 0.7,
    "Dy": 1.1,
    "Dz": 0.9,
    "NumCellsX": 4,
    "NumCellsY": 5,
    "NumCellsZ": 3,
    "StartFromTime": 0,
    "StartTimeStep": 0,
    "LastTimestep": LastTimestep,
    "n0": 1e10,
    "k_particles_reservation": 0,
    "Collider": "None",
    "Collisions": [],
    "DampingType": "None",
    "Boundary_conditions": [
        {"periodic": {"face": face}} for face in ("XMIN", "YMIN", "ZMIN")
    ],
    "ExternalFieldB": [{"uniform_field": {"value": [0, 0, 0.2]}}],
    "diagnostics": {
        "outputs": [
            {"type": "energy_balance", "interval": 1},
            {"type": "console_summary"},
            {"type": "recovery", "interval": 8},
        ]
    },
}


def species(name, charge, mass, drift):
    # Both species use the same spatial distribution; injection is disabled.
    return {
        "Name": name,
        "Charge": charge,
        "Mass": mass,
        "Density": 1,
        "NumPartPerCell": 2,
        "distribution": [{
            "type": "initial",
            "dist_space": {
                "type": "rectangle",
                "center": [1.4, 2.75, 1.35],
                "half_length": [1.4, 2.75, 1.35],
            },
            "dist_pulse": {
                "type": "gaussian", "mean": drift, "sigma": [0.01, 0.01, 0.01],
            },
            "density": 1,
        }],
    }


particles_config = {"particles": [
    species("Electrons", -1, 1, [0.04, 0, 0]),
    species("Ions", 1, 1836, [0, 0, 0]),
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

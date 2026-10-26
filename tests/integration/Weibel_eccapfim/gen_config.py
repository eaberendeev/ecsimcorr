#!/usr/bin/env python3
"""Weibel instability, ECCAPFIM, periodic 2D3V plasma.

The x/y domain, mesh and temperatures match Weibel_ecsim. One periodic
cell in z and 16 particles/cell give a reduced physical validation rather
than the 30x30x5, 500 particles/cell paper reproduction.
"""
import json
import math

IntegrationHeavy = True  # Several minutes of nonlinear orbit iterations on four CPU threads.
NumProcs = 4
LastTime = 360.0
Dt = 1.2
LastTimestep = int(round(LastTime / Dt))
DirName = "Res_Weibel_eccapfim_Dx_0.3_np_16_Dt_1.2"

system_config = {
    "Scheme": "eccapfim",
    "Dx": 0.3,
    "Dy": 0.3,
    "Dz": 0.3,
    "Dt": Dt,
    "NumCellsX": 30,
    "NumCellsY": 30,
    "NumCellsZ": 1,
    "StartTimeStep": 0,
    "LastTimestep": LastTimestep,
    "StartFromTime": 0,
    "Collider": "None",
    "Collisions": [],
    "DampingType": "None",
    "n0": 1e13,
    "k_particles_reservation": -1,
    "Boundary_conditions": [
        {"periodic": {"face": face}} for face in ("XMIN", "YMIN", "ZMIN")
    ],
    "diagnostics": {
        "outputs": [
            {"type": "energy_balance", "interval": 1},
            {"type": "console_summary"},
            {"type": "fields_2d", "interval": 5, "fields": ["E", "B"], "planes": {"z": [0.15]}},
            {"type": "density_2d", "interval": 10, "species": "all", "planes": {"z": [0.15]}},
            {"type": "recovery", "interval": 100},
        ],
    },
}

particles_config = {"particles": [{
    "Name": "Electrons",
    "Charge": -1.0,
    "Density": 1.0,
    "Mass": 1.0,
    "NumPartPerCell": 16,
    "distribution": [{
        "type": "initial",
        "dist_space": {
            "type": "rectangle", "center": [4.5, 4.5, 0.15], "half_length": [4.5, 4.5, 0.15],
        },
        "dist_pulse": {"type": "gaussian", "mean": [0, 0, 0], "sigma": [0.1, 0.1, 1.0]},
        "density": 1.0,
    }],
}]}


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

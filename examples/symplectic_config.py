#!/usr/bin/env python3
"""Generate a small periodic electron-positron plasma in a NEW run directory."""
import argparse
import json
from pathlib import Path


def configuration(steps=200):
    system = {
        "Scheme": "symplectic", "Dt": 0.05,
        "NumCellsX": 4, "NumCellsY": 4, "NumCellsZ": 4,
        "Dx": 1.0, "Dy": 1.0, "Dz": 1.0,
        "StartFromTime": 0, "StartTimeStep": 0, "LastTimestep": steps,
        "n0": 1e10, "k_particles_reservation": 0, "Collisions": [],
        "Boundary_conditions": [{"periodic": {"face": face}} for face in ("XMIN", "YMIN", "ZMIN")],
        "co_locate_species": True,
        "Symplectic": {"initial_electric": [0.01, 0, 0], "initial_magnetic": [0, 0, 0.1]},
        "diagnostics": {"outputs": [{"type": "energy_balance", "interval": 1}]},
    }
    distribution = [{
        "type": "initial",
        "dist_space": {"type": "rectangle", "center": [2, 2, 2], "half_length": [2, 2, 2]},
        "dist_pulse": {"type": "gaussian", "mean": [0, 0, 0], "sigma": [0, 0, 0]},
        "density": 1,
    }]
    particles = {"particles": [
        {"Name": name, "Charge": charge, "Mass": 1, "Density": 1, "NumPartPerCell": 8,
         "distribution": distribution}
        for name, charge in (("Electrons", -1), ("Positrons", 1))
    ]}
    return system, particles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--steps", type=int, default=200)
    args = parser.parse_args()
    if args.steps < 0:
        parser.error("--steps must be nonnegative")
    args.directory.mkdir(parents=True, exist_ok=False)
    system, particles = configuration(args.steps)
    for name, config in (("system_config.json", system), ("particles_config.json", particles)):
        (args.directory / name).write_text(json.dumps(config, indent=2) + "\n")


if __name__ == "__main__":
    main()

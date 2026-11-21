#!/usr/bin/env python3
"""Run a small matched outflow problem without rebuilding or deleting results.

Periodic x/y, absorbing z; --reflect replaces z boundaries by the shared
energy-selective electron reflector. Requires NumPy and an existing beren3d.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np


def configurations(scheme, dt, reflecting):
    steps = round(2.4 / dt)
    walls = ({"electron_reflection": [
        {"face": face, "radius": 100, "energy_threshold": 1000}
        for face in ("ZMIN", "ZMAX")
    ]} if reflecting else {"open": [{"face": "ZMIN"}, {"face": "ZMAX"}]})
    config = {
        "Scheme": scheme, "Dt": dt, "Dx": 0.7, "Dy": 1.1, "Dz": 0.9,
        "NumCellsX": 4, "NumCellsY": 4, "NumCellsZ": 8,
        "StartFromTime": 0, "StartTimeStep": 0, "LastTimestep": steps,
        "n0": 1e10, "k_particles_reservation": 0, "Collider": "None",
        "Collisions": [], "DampingType": "None", "co_locate_species": True,
        "Boundary_conditions": [{"periodic": {"face": "XMIN"}},
                                {"periodic": {"face": "YMIN"}}, walls],
        "diagnostics": {"outputs": [
            {"type": "energy_balance", "interval": 1},
            {"type": "boundary_stats", "interval": 1},
            {"type": "recovery", "interval": steps},
        ]},
    }
    particles = []
    for name, charge, mass, drift in (("Electrons", -1, 1, 0.45), ("Ions", 1, 25, 0.1)):
        particles.append({
            "Name": name, "Charge": charge, "Mass": mass, "Density": 1,
            "NumPartPerCell": 8,
            "distribution": [{
                "type": "initial", "density": 1,
                "dist_space": {"type": "rectangle", "center": [1.4, 2.2, 3.6],
                               "half_length": [1.4, 2.2, 3.6]},
                "dist_pulse": {"type": "gaussian", "mean": [0, 0, drift],
                               "sigma": [0.02, 0.02, 0.04]},
            }],
        })
    return config, {"particles": particles}


def read_field(path):
    with path.open("rb") as stream:
        shape = np.fromfile(stream, dtype=np.int32, count=3)
        values = np.fromfile(stream, dtype=np.float64)
    if tuple(shape) != (7, 7, 11) or values.size != 3 * np.prod(shape):
        raise ValueError(f"Unexpected recovery field layout: {path}")
    field = values.reshape(*shape, 3)[1:5, 1:5, :, :]
    if not np.all(np.isfinite(field)):
        raise ValueError(f"Nonfinite field: {path}")
    return field


def run_case(binary, directory, scheme, dt, reflecting):
    config, particles = configurations(scheme, dt, reflecting)
    directory.mkdir()
    for name, value in (("system_config", config), ("particles_config", particles)):
        (directory / f"{name}.json").write_text(json.dumps(value, indent=2) + "\n")
    env = dict(os.environ, OMP_NUM_THREADS="2", VALIDATION_PERIODICITY="1000000")
    start = time.monotonic()
    with (directory / "run.log").open("w") as log:
        subprocess.run([str(binary)], cwd=directory, env=env, stdout=log,
                       stderr=subprocess.STDOUT, timeout=180, check=True)
    rows = np.atleast_1d(np.genfromtxt(directory / "energy.txt", names=True))
    boundary = np.atleast_1d(np.genfromtxt(directory / "boundary.txt", names=True))
    for table in (rows, boundary):
        if len(table) != config["LastTimestep"] + 1 or any(
                not np.all(np.isfinite(table[name])) for name in table.dtype.names):
            raise ValueError(f"Incomplete or nonfinite diagnostics: {directory}")
    if not np.isclose(rows["Time"][-1], 2.4):
        raise ValueError(f"Wrong final time: {directory}")
    initial = rows["Electrons"][0] + rows["Ions"][0]
    total = rows["Electrons"] + rows["Ions"] + rows["energyFieldE"] + rows["energyFieldB"]
    closed = total + np.cumsum(rows["totalLostEnergy"])
    metrics = {"wall_seconds": time.monotonic() - start,
               "initial_kinetic": float(initial),
               "relative_energy_balance": float(np.max(np.abs(closed - closed[0])) / initial),
               "electric_energy": float(rows["energyFieldE"][-1])}
    for name in ("Electrons", "Ions"):
        count = rows[name + "Particles"]
        lost = sum(boundary[name + "_" + face + "LostN"].sum() for face in ("ZMIN", "ZMAX"))
        reflected = sum(boundary[name + "_" + face + "ReflN"].sum() for face in ("ZMIN", "ZMAX"))
        if count[0] != count[-1] + lost or np.any(np.diff(count) > 0):
            raise ValueError(f"Boundary particle accounting failed: {directory}/{name}")
        if lost + reflected == 0:
            raise ValueError(f"No wall interaction: {directory}/{name}")
        metrics[name] = {"initial_count": int(count[0]), "final_count": int(count[-1]),
                         "lost": int(lost), "reflected": int(reflected),
                         "kinetic": float(rows[name][-1])}
    if metrics["relative_energy_balance"] > 1e-7:
        raise ValueError(f"Energy balance failed: {directory}")
    if reflecting and metrics["Electrons"]["lost"] != 0:
        raise ValueError(f"Low energy electrons escaped reflector: {directory}")
    return metrics, read_field(directory / "Recovery/Fields/FieldE.backup")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("output_directory", type=Path, help="New directory; existing results are preserved")
    parser.add_argument("--reflect", action="store_true")
    args = parser.parse_args()
    binary = args.binary.resolve()
    root = args.output_directory.resolve()
    root.mkdir()
    summary = {}
    initial = []
    for dt in (0.2, 0.1, 0.05):
        cases = {}
        for scheme in ("ecsim_corr", "eccapfim"):
            name = f"{scheme}_dt{dt:g}"
            metrics, field = run_case(binary, root / name, scheme, dt, args.reflect)
            initial.append(metrics["initial_kinetic"])
            cases[scheme] = metrics, field
            print(name, json.dumps(metrics), flush=True)
        corr, implicit = cases["ecsim_corr"], cases["eccapfim"]
        comparison = {"electric_difference_energy_over_initial_kinetic":
                      float(0.5 * np.sum((corr[1] - implicit[1]) ** 2) / initial[0])}
        for name in ("Electrons", "Ions"):
            comparison[name + "_count_difference_over_initial"] = abs(
                corr[0][name]["final_count"] - implicit[0][name]["final_count"]) / corr[0][name]["initial_count"]
            comparison[name + "_kinetic_difference_over_initial"] = abs(
                corr[0][name]["kinetic"] - implicit[0][name]["kinetic"]) / initial[0]
        summary[str(dt)] = {"cases": {s: c[0] for s, c in cases.items()}, "comparison": comparison}
        (root / "comparison.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
        print("comparison", dt, json.dumps(comparison), flush=True)
    if max(initial) - min(initial) > 1e-12 * max(initial):
        raise ValueError("Initial particle samples differ")
    # The methods use different response shapes: compare resolved bulk quantities,
    # not equality of individual noisy grid values or identical wall-hit steps.
    fine = summary["0.05"]["comparison"]
    for metric, value in fine.items():
        limit = (0.005 if "count" in metric else
                 1e-4 if metric.startswith("electric") else 0.001)
        if value > limit:
            raise ValueError(f"Fine-step comparison failed: {metric}={value} > {limit}")


if __name__ == "__main__":
    main()

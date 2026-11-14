#!/usr/bin/env python3
"""Benchmark both schemes from identical integer-time particles and fields.

For the recorded comparison STATE_DIRECTORY is eccapfim_dt0.6: its last
recovery is t=120 (interval=200, final step=300, dt=0.6). Each trial performs
30 steps and excludes the first five. Use Release with USE_TIMERS=OFF.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np


def run_trial(binary, directory, cfg, particles):
    directory.mkdir()
    for name, value in (("system_config", cfg), ("particles_config", particles)):
        (directory / f"{name}.json").write_text(json.dumps(value, indent=2))
    env = dict(os.environ, OMP_NUM_THREADS="8", VALIDATION_PERIODICITY="1000000")
    start = time.monotonic()
    with (directory / "run.log").open("w") as log:
        result = subprocess.run([str(binary)], cwd=directory, env=env, stdout=log,
                                stderr=subprocess.STDOUT, timeout=120)
    wall_seconds = time.monotonic() - start
    if result.returncode:
        raise RuntimeError(f"Calculation failed: {directory}")
    rows = np.genfromtxt(directory / "energy.txt", names=True)
    if len(rows) != 31 or any(not np.all(np.isfinite(rows[n])) for n in rows.dtype.names):
        raise ValueError(f"Incomplete or nonfinite diagnostics: {directory}")
    if not np.all(rows["ElectronsParticles"] == rows["ElectronsParticles"][0]):
        raise ValueError(f"Particle count changed: {directory}")
    total = rows["Electrons"] + rows["energyFieldE"] + rows["energyFieldB"]
    drift = float(np.max(np.abs(total - total[0])) / total[0])
    if drift > 1e-7:
        raise ValueError(f"Energy drift exceeds 1e-7: {directory}")
    timing = np.genfromtxt(directory / "Performance/globalFunctions.time", names=True)
    values = timing["Total"][timing["Time"] > 5]
    if len(values) != 25 or not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError(f"Invalid warm-step timings: {directory}")
    return {
        "wall_seconds": wall_seconds,
        "warm_steps": len(values),
        "mean_step_seconds": float(np.mean(values)),
        "median_step_seconds": float(np.median(values)),
        "initial_kinetic": float(rows["Electrons"][0]),
        "relative_energy_drift": drift,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("state_directory", type=Path, help="Common state's configs and Recovery directory")
    parser.add_argument("output_directory", type=Path, help="A new directory; previous results are preserved")
    args = parser.parse_args()
    root = args.output_directory.resolve()
    root.mkdir()
    source, binary = args.state_directory.resolve(), args.binary.resolve()
    shutil.copytree(source / "Recovery", root / "Recovery")
    cfg_base = json.loads((source / "system_config.json").read_text())
    particles = json.loads((source / "particles_config.json").read_text())
    records = []
    for dt in (1.2, 0.6):
        for trial in range(3):
            # Alternate order to reduce bias from cache and background activity.
            order = ("ecsim_corr", "eccapfim") if trial % 2 == 0 else ("eccapfim", "ecsim_corr")
            for scheme in order:
                name = f"{scheme}_dt{dt:g}_trial{trial + 1}"
                cfg = dict(cfg_base, Scheme=scheme, Dt=dt, StartFromTime=120, StartTimeStep=0, LastTimestep=30)
                cfg["diagnostics"] = {"outputs": [{"type": "energy_balance", "interval": 1}]}
                # The shared parent Recovery stays read-only: trials write no recovery output.
                record = run_trial(binary, root / name, cfg, particles)
                record.update(name=name, scheme=scheme, dt=dt, trial=trial + 1, threads=8)
                records.append(record)
                print(json.dumps(record), flush=True)
                (root / "benchmark.json").write_text(json.dumps(records, indent=2))
    initial = [record["initial_kinetic"] for record in records]
    if max(initial) - min(initial) > 1e-11 * max(initial):
        raise ValueError("Benchmark trials did not load the same initial kinetic energy")
    summary = {}
    for dt in (1.2, 0.6):
        row = {}
        for scheme in ("ecsim_corr", "eccapfim"):
            values = [r["mean_step_seconds"] for r in records if r["scheme"] == scheme and r["dt"] == dt]
            row[scheme] = {
                "median_trial_mean": float(np.median(values)),
                "minimum_trial_mean": min(values), "maximum_trial_mean": max(values),
            }
        row["eccapfim_over_corr"] = row["eccapfim"]["median_trial_mean"] / row["ecsim_corr"]["median_trial_mean"]
        summary[str(dt)] = row
    (root / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

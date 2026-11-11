"""Run one desktop workload through the production loop, including collisions."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import time

from gen_config import add_config_arguments, make_config, write_config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_config_arguments(parser)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="new run directory; never overwritten")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    binary = args.binary.resolve()
    if not binary.is_file() or args.threads < 1:
        parser.error("provide an existing binary and positive thread count")
    try:
        system, particles = make_config(args.case, args.scheme, args.cells, args.ppc, args.steps)
    except ValueError as error:
        parser.error(str(error))
    args.output.mkdir(parents=True, exist_ok=False)
    write_config(args.output, system, particles)
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS=str(args.threads), OMP_DYNAMIC="FALSE",
                       OMP_PROC_BIND="close", OMP_PLACES="cores", TIMER_START_STEP="2",
                       TIMER_END_STEP=str(args.steps), TIMERED_SIM_STEPS=str(args.steps),
                       VALIDATION_PERIODICITY=str(args.steps))
    metadata = {
        "case": args.case, "scheme": args.scheme, "cells": args.cells,
        "ppc": system["NumPartPerCell"], "dt": system["Dt"], "steps": args.steps,
        "threads": args.threads, "platform": platform.platform(), "cpu": platform.processor(),
        "binary": str(binary), "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "environment": {key: value for key, value in environment.items()
                        if key.startswith(("OMP_", "TIMER", "USE_MIXED", "COLLISION_", "VALIDATION_"))},
    }
    (args.output / "run.json").write_text(json.dumps(metadata, indent=2) + "\n")
    start = time.perf_counter()
    with (args.output / "stdout.log").open("w") as log:
        subprocess.run([str(binary)], cwd=args.output, env=environment,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    elapsed = time.perf_counter() - start
    profile = args.output / f"profile_it_{args.steps}.json"
    if not profile.is_file():
        raise RuntimeError("No profile: build with USE_TIMERS=ON; see stdout.log")
    events = json.loads(profile.read_text())
    durations = {}
    for event in events:
        if event.get("ph") == "X" and event.get("tid") == 0:
            durations.setdefault(event["name"], []).append(event["dur"] / 1.e6)
    wall_times = [float(value) for value in re.findall(
        r"Step wall time: ([\d.eE+\-]+) s", (args.output / "beren3d.log").read_text())]
    if len(wall_times) != args.steps:
        raise RuntimeError(f"Expected {args.steps} completed steps, got {len(wall_times)}")
    samples = wall_times[1:]
    report = {
        "measured_steps": list(range(2, args.steps + 1)), "step_seconds": samples,
        "mean_step_seconds": statistics.mean(samples), "median_step_seconds": statistics.median(samples),
        "process_seconds_including_initialization": elapsed,
        "master_thread_timers": {name: {"calls": len(values), "total_seconds": sum(values)}
                                 for name, values in durations.items()},
    }
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"{args.case}/{args.scheme}: {statistics.mean(samples):.6f} s/step; {args.output}")


if __name__ == "__main__":
    main()

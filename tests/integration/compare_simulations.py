#!/usr/bin/env python3
"""Compare complete particle/field states from two independently built binaries.

Run this script from the validation checkout (based on master); neither tested
checkout needs to contain these tests. Existing binaries may include uncommitted
changes. No checkout, build, commit or timing comparison is performed.

  python3 tests/integration/compare_simulations.py --baseline /tmp/master/beren3d \
    --candidate /tmp/branch/beren3d --profile smoke --output /tmp/new-results

Use --profile volume50 or dense50 for 50x50x50 cells with 32 or 100 particles/cell.
Large cases are explicit, never part of the default unit suite. Options --cells,
--npc, --steps, --dt and --case select other inputs. --case two_species creates
an electron/ion plasma; --case periodic removes the background magnetic field.
Each run independently checks finite states, conserved periodic particle counts
and relative energy balance <= 1e-8; --continuity also checks ECSIM-Corr charge
continuity (Linf <= 1e-10). Equality with master is a regression check, not a
replacement for these physical invariants. Executables are copied once into the
output directory so a later rebuild cannot change the binary mid-test.
"""
import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import runpy
import shutil
import struct
import subprocess

from run_integration import check_energy_conservation, parse_energy_file

ROOT = Path(__file__).resolve().parents[2]
PROFILES = {
    "smoke": ((8, 6, 4), 16, 3),
    "local": ((16, 12, 8), 64, 4),
    "volume50": ((50, 50, 50), 32, 3),
    "dense50": ((50, 50, 50), 100, 3),
}


def configuration(scheme, profile, case="magnetized"):
    cells, npc, steps = PROFILES[profile]
    template = runpy.run_path(str(ROOT / "tests/integration/Weibel_ecsim/gen_config.py"))
    system = copy.deepcopy(template["system_config"])
    particles = copy.deepcopy(template["particles_config"])
    # Small dt keeps the old and new single-cell trajectory envelopes valid.
    system.update(Scheme=scheme, Dt=0.1, LastTimestep=steps, VerboseStep=False,
                  RecoveryInterval=steps, NumPartPerCell=npc, k_particles_reservation=0)
    system["ExternalFieldB"] = [{"uniform_field": {"value": [0, 0, 0] if case == "periodic" else [0.03, -0.02, 0.1]}}]
    for axis, count in zip("XYZ", cells):
        system["NumCells" + axis] = count
    half = [0.5 * count * system["D" + axis.lower()] for axis, count in zip("XYZ", cells)]
    for species in particles["particles"]:
        species["NumPartPerCell"] = npc
        space = species["distribution"][0]["dist_space"]
        space.update(center=half, half_length=half)
    if case == "two_species":
        ions = copy.deepcopy(particles["particles"][0])
        ions.update(Name="Ions", Charge=1.0, Mass=1836.0)
        ions["distribution"][0]["dist_pulse"]["sigma"] = [0.1, 0.1, 0.1]
        particles["particles"].append(ions)
    outputs = [{"type": "energy_balance", "interval": 1},
               {"type": "recovery", "interval": steps}]
    outputs.append({"type": "density_2d", "interval": steps,
                    "species": "all", "planes": {"z": [half[2]]}})
    system["diagnostics"] = {"outputs": outputs}
    return system, particles


def finite_values(values, label):
    values = list(values)
    if not values or not all(math.isfinite(v) for v in values):
        raise ValueError(f"{label}: empty or non-finite data")
    return values


def read_field(path):
    payload = path.read_bytes()
    if len(payload) < 12:
        raise ValueError(f"{path}: missing field header")
    shape = struct.unpack("=3i", payload[:12])
    count = 3 * math.prod(shape)
    if min(shape) <= 0 or len(payload) != 12 + 8 * count:
        raise ValueError(f"{path}: invalid field dimensions or payload")
    return shape, finite_values(struct.unpack(f"={count}d", payload[12:]), str(path))


def read_plane(path):
    payload = path.read_bytes()
    if len(payload) < 8:
        raise ValueError(f"{path}: missing plane header")
    shape = struct.unpack("=2f", payload[:8])
    if not all(math.isfinite(v) and v > 0 and v.is_integer() for v in shape):
        raise ValueError(f"{path}: invalid plane dimensions")
    count = 3 * math.prod(int(v) for v in shape)
    if len(payload) != 8 + 4 * count:
        raise ValueError(f"{path}: truncated or excess plane payload")
    return shape, finite_values(struct.unpack(f"={count}f", payload[8:]), str(path))


def read_particles(path):
    payload = path.read_bytes()
    if len(payload) < 4:
        raise ValueError(f"{path}: missing particle count")
    count, = struct.unpack("=i", payload[:4])
    # Supported default build: ParticleSimple contains coord and velocity (6 doubles).
    if count <= 0 or len(payload) != 4 + count * 48:
        raise ValueError(f"{path}: invalid count or unsupported Particle layout (expected 6 doubles)")
    records = sorted(struct.iter_unpack("=6d", payload[4:]))
    return count, finite_values((v for row in records for v in row), str(path))


def compare_values(actual, expected, rtol, atol, label):
    if len(actual) != len(expected) or not actual:
        raise ValueError(f"{label}: missing values or different lengths")
    finite_values(actual, label)
    finite_values(expected, label)
    # One scale per physical array/column avoids relative errors at zero.
    scale = max(abs(v) for v in expected)
    difference = max(abs(a - b) for a, b in zip(actual, expected))
    limit = atol + rtol * scale
    if difference > limit:
        raise ValueError(f"{label}: max error {difference:.6g} > {limit:.6g}")
    return {"max_abs_error": difference, "reference_scale": scale, "limit": limit}


def conserved_counts(rows, species):
    counts = {}
    for name in species:
        column = name + "Particles"
        initial = rows[0].get(column, float("nan"))
        if not math.isfinite(initial) or initial <= 0 or int(initial) != initial:
            raise ValueError(f"{name}: missing or invalid initial particle count")
        if any(row.get(column) != initial for row in rows):
            raise ValueError(f"{name}: particle number is not conserved in the periodic box")
        counts[name] = int(initial)
    return counts


def check_continuity(directory, steps, dt, limit=1e-10):
    lines = (directory / "continuity.txt").read_text().splitlines()
    if not lines or lines[0].split() != ["Time", "L2", "Linf"] or len(lines) != steps + 1:
        raise ValueError("incomplete charge continuity history")
    values = []
    for step, line in enumerate(lines[1:], 1):
        row = finite_values(map(float, line.split()), "continuity")
        if len(row) != 3 or not math.isclose(row[0], step * dt, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError("wrong charge continuity timestep")
        if row[1] < 0 or row[2] < 0 or row[2] > limit:
            raise ValueError(f"charge continuity residual {row[2]} exceeds {limit}")
        values.append(row)
    return {"max_L2": max(row[1] for row in values), "max_Linf": max(row[2] for row in values), "limit": limit}


def compare_outputs(candidate, baseline, species, steps, dt, rtol, atol, require_density=False, energy_tol=1e-8):
    results = {}
    tables = []
    for directory in (candidate, baseline):
        rows, error = parse_energy_file(str(directory))
        if error or not rows or len(rows) != steps + 1:
            raise ValueError(f"{directory}: {error or 'incomplete energy history'}")
        for index, row in enumerate(rows):
            if not math.isclose(row["Time"], index * dt, rel_tol=1e-10, abs_tol=1e-12):
                raise ValueError(f"{directory}: wrong diagnostic timestep {index}")
        conserved_counts(rows, species)
        ok, message = check_energy_conservation(str(directory), tol=energy_tol)
        if not ok:
            raise ValueError(f"{directory}: {message}")
        tables.append(rows)
    if tables[0][0].keys() != tables[1][0].keys():
        raise ValueError("energy columns differ")
    for column in tables[0][0]:
        results["energy/" + column] = compare_values(
            [row[column] for row in tables[0]], [row[column] for row in tables[1]],
            rtol, atol, "energy/" + column)
    for field in ("FieldE", "FieldB"):
        relative = Path("Recovery/Fields") / (field + ".backup")
        a_shape, actual = read_field(candidate / relative)
        b_shape, expected = read_field(baseline / relative)
        if a_shape != b_shape:
            raise ValueError(f"{field}: shapes differ")
        results[field] = compare_values(actual, expected, rtol, atol, field)
    for name in species:
        relative = Path("Recovery/Particles") / name / (name + ".backup")
        a_count, actual = read_particles(candidate / relative)
        b_count, expected = read_particles(baseline / relative)
        if a_count != b_count or a_count != tables[0][0][name + "Particles"] or b_count != tables[1][0][name + "Particles"]:
            raise ValueError(f"{name}: particle counts differ from each other or the conserved initial count")
        # Positions and velocities have different units and need separate scales.
        for axis in range(6):
            label = f"particles/{name}/{axis}"
            results[label] = compare_values(actual[axis::6], expected[axis::6], rtol, atol, label)
        results[f"particles/{name}/count"] = a_count
    # Density diagnostics are native float32 plane snapshots (header + data).
    density_files = sorted(path for path in baseline.rglob("Density*") if path.is_file())
    candidate_density_files = sorted(path.relative_to(candidate) for path in candidate.rglob("Density*") if path.is_file())
    if candidate_density_files != [path.relative_to(baseline) for path in density_files]:
        raise ValueError("density snapshot file sets differ")
    if require_density and not density_files:
        raise ValueError("missing density snapshots")
    for path in density_files:
        relative = path.relative_to(baseline)
        a_shape, actual = read_plane(candidate / relative)
        b_shape, expected = read_plane(path)
        if a_shape != b_shape:
            raise ValueError(f"{relative}: density snapshot dimensions differ")
        results[str(relative)] = compare_values(actual, expected, max(rtol, 2e-6), max(atol, 1e-7), str(relative))

    return results


def execute(binary, directory, settings, particles, threads, timeout, mixed_precision=False):
    directory.mkdir(parents=True, exist_ok=False)
    for name, config in (("system", settings), ("particles", particles)):
        (directory / (name + "_config.json")).write_text(json.dumps(config, indent=2))
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE", OMP_PROC_BIND="close", OMP_PLACES="cores")
    # Keep both processes on the requested solver algorithm.
    env["USE_MIXED_PRECISION"] = "1" if mixed_precision else "0"
    env["OMP_THREAD_LIMIT"] = str(threads)
    env["VALIDATION_PERIODICITY"] = "10"
    with (directory / "stdout.log").open("w") as output:
        subprocess.run([str(binary)], cwd=directory, env=env, stdout=output,
                       stderr=subprocess.STDOUT, check=True, timeout=timeout)


def binary_info(path):
    info = {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    for directory in path.parents:
        cache = directory / "CMakeCache.txt"
        if not cache.is_file():
            continue
        settings = {}
        for line in cache.read_text().splitlines():
            if ":" in line and "=" in line and not line.startswith(("#", "//")):
                name, value = line.split("=", 1)
                settings[name.split(":", 1)[0]] = value
        info["build"] = {key: settings[key] for key in
                         ("CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER", "CMAKE_CXX_FLAGS", "CMAKE_CXX_FLAGS_RELEASE",
                          "USE_TIMERS", "MEMORY_TIMERS_LEVEL", "CMAKE_HOME_DIRECTORY") if key in settings}
        source = settings.get("CMAKE_HOME_DIRECTORY")
        if source:
            revision = subprocess.run(["git", "-C", source, "rev-parse", "HEAD"], capture_output=True, text=True)
            diff = subprocess.run(["git", "-C", source, "diff", "HEAD"], capture_output=True)
            if revision.returncode == 0 and diff.returncode == 0:
                info["source_head"] = revision.stdout.strip()
                info["source_diff_sha256"] = hashlib.sha256(diff.stdout).hexdigest()
        break
    return info


def snapshot_binary(binary, destination):
    # Pin the executable for the whole run even if the developer rebuilds later.
    shutil.copy2(binary, destination)
    info = binary_info(binary)
    info["sha256"] = hashlib.sha256(destination.read_bytes()).hexdigest()
    info["snapshot"] = str(destination)
    return info


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profile", choices=PROFILES, default="smoke")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--case", choices=("periodic", "magnetized", "two_species"), default="magnetized")
    parser.add_argument("--schemes", nargs="+", choices=("ecsim", "ecsim_corr"), default=["ecsim", "ecsim_corr"])
    parser.add_argument("--cells", nargs=3, type=int, help="override Nx Ny Nz")
    parser.add_argument("--npc", type=int, help="override particles per cell")
    parser.add_argument("--steps", type=int, help="override steps")
    parser.add_argument("--dt", type=float, default=0.1)
    parser.add_argument("--rtol", type=float, default=1e-8)
    parser.add_argument("--atol", type=float, default=1e-10)
    parser.add_argument("--mixed-precision", action="store_true", help="use the production mixed precision solver instead of double only")
    parser.add_argument("--energy-tol", type=float, default=1e-8)
    parser.add_argument("--continuity", action="store_true", help="also check ECSIM-Corr charge continuity (Linf <= 1e-10)")
    parser.add_argument("--timeout", type=float, default=3600)
    args = parser.parse_args()
    if args.threads <= 0 or args.threads > len(os.sched_getaffinity(0)):
        parser.error("threads must fit the available CPU affinity")
    if any(v is not None and v <= 0 for v in (args.npc, args.steps)) or (args.cells and min(args.cells) <= 0):
        parser.error("cell counts, npc and steps must be positive")
    for value in (args.dt, args.rtol, args.atol, args.energy_tol, args.timeout):
        if not math.isfinite(value) or value <= 0:
            parser.error("dt, tolerances and timeout must be finite and positive")
    args.baseline, args.candidate, args.output = args.baseline.resolve(), args.candidate.resolve(), args.output.resolve()
    for binary in (args.baseline, args.candidate):
        if not binary.is_file() or not os.access(binary, os.X_OK):
            parser.error(f"not an executable file: {binary}")
    cells, npc, steps = PROFILES[args.profile]
    cells, npc, steps = tuple(args.cells or cells), args.npc or npc, args.steps or steps
    PROFILES["requested"] = (cells, npc, steps)
    args.output.mkdir(parents=True, exist_ok=False)
    snapshots = args.output / "binaries"
    snapshots.mkdir()
    binaries = {label: snapshots / label for label in ("baseline", "candidate")}
    information = {label: snapshot_binary(getattr(args, label), binaries[label]) for label in binaries}
    report = {"profile": args.profile, "case": args.case, "cells": cells, "npc": npc, "steps": steps, "dt": args.dt,
              "threads": args.threads, "rtol": args.rtol, "atol": args.atol, "energy_tol": args.energy_tol, "mixed_precision": args.mixed_precision,
              "baseline": information["baseline"], "candidate": information["candidate"], "results": []}
    try:
        for scheme in args.schemes:
            system, particles = configuration(scheme, "requested", args.case)
            system["Dt"] = args.dt
            if args.continuity and scheme == "ecsim_corr":
                system["diagnostics"]["outputs"].append({"type": "charge_continuity", "interval": 1})
            directories = {}
            for label in ("baseline", "candidate"):
                directory = args.output / (scheme + "-" + label)
                execute(binaries[label], directory, system, particles, args.threads, args.timeout, args.mixed_precision)
                directories[label] = directory
            names = [entry["Name"] for entry in particles["particles"]]
            comparison = compare_outputs(directories["candidate"], directories["baseline"], names,
                                         steps, args.dt, args.rtol, args.atol, require_density=True, energy_tol=args.energy_tol)
            if args.continuity and scheme == "ecsim_corr":
                comparison["charge_continuity"] = {label: check_continuity(directory, steps, args.dt)
                                                   for label, directory in directories.items()}
            report["results"].append({"scheme": scheme, "comparison": comparison})
            print(f"{scheme}: numerical comparison PASS", flush=True)
        report["status"] = "passed"
    except Exception as error:
        report.update(status="failed", error=str(error))
        raise
    finally:
        (args.output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

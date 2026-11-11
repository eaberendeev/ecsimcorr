#!/usr/bin/env python3
"""Run physical scenarios against the current build; reference states are optional.

  python3 tests/integration/validate_physics.py --binary build/beren3d
  python3 tests/integration/validate_physics.py --binary master-build/beren3d \
    --save-reference /tmp/master-reference --threads 4
  python3 tests/integration/validate_physics.py --binary build/beren3d \
    --reference /tmp/master-reference --threads 4

The primary checks are conserved particle counts, finite diagnostics, relative
energy balance <= 1e-8 and ECSIM-Corr charge continuity Linf <= 1e-10. The
symplectic scheme keeps energy only bounded, so it checks relative energy error
<= 1e-3, Gauss law drift <= 1e-10 and div B <= 1e-12 instead. Three
scenarios cover unmagnetized, magnetized and electron/ion plasmas. No second
executable or clean working tree is required. CMake target run_physics_validation
builds and checks the current program, including uncommitted source changes.

Reference creation runs the SAME build repeatedly (default 3), measures natural
state differences and adds a configurable 2..10 safety factor (default 5).
A one-ULP noise floor covers a bitwise-repeatable baseline before applying that factor. Reference comparison checks
every E/B value and every particle coordinate/velocity, not executable bytes or
aggregate checksums. References use a JSON manifest and explicit little-endian
float64 arrays, are created separately and must not be committed for large cases.
The optional export adapter currently supports the default six-double Particle.

Profiles volume50/dense50 use 50^3 cells and 32/100 particles per cell. They are
explicit large validations; the default smoke profile uses 3,072 particles/species
and three steps. --cells, --npc, --steps, --dt and --mixed-precision override inputs.
"""
import argparse
import copy
from array import array
import hashlib
import json
import math
import os
from pathlib import Path
import runpy
import shutil
import struct
import subprocess
import sys
import tempfile

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
    field = [0, 0, 0] if case == "periodic" else [0.03, -0.02, 0.1]
    if scheme == "symplectic":
        # The symplectic B is evolved state, not an external field.
        system.pop("ExternalFieldB", None)
        system["Symplectic"] = {"initial_magnetic": field}
    else:
        system["ExternalFieldB"] = [{"uniform_field": {"value": field}}]
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
    outputs = [{"type": "energy_balance", "interval": 1}]
    if scheme == "ecsim_corr":
        outputs.append({"type": "charge_continuity", "interval": 1})
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


def check_symplectic_invariants(rows, energy_tol, gauss_tol, div_tol=1e-12):
    # The symplectic step keeps energy only bounded; Gauss law and div B hold to round-off.
    total = rows[0]["energyTotal"]
    energy = max(abs(row["energyError"]) for row in rows) / total
    gauss = max(row["gaussDrift"] for row in rows)
    divergence = max(row["divB"] for row in rows)
    if not energy <= energy_tol:
        raise ValueError(f"symplectic relative energy error {energy} exceeds {energy_tol}")
    if not gauss <= gauss_tol:
        raise ValueError(f"symplectic Gauss law drift {gauss} exceeds {gauss_tol}")
    if not divergence <= div_tol:
        raise ValueError(f"symplectic div B {divergence} exceeds {div_tol}")
    return {"energy_relative": energy, "energy_limit": energy_tol, "gauss_drift": gauss, "div_b": divergence}


def validate_run(directory, system, particles, energy_tol, continuity_tol=1e-10, symplectic_energy_tol=1e-3):
    steps, dt = system["LastTimestep"], system["Dt"]
    names = [entry["Name"] for entry in particles["particles"]]
    rows, error = parse_energy_file(str(directory))
    if error or not rows or len(rows) != steps + 1:
        raise ValueError(f"{directory}: {error or 'incomplete energy history'}")
    for index, row in enumerate(rows):
        if not math.isclose(row.get("Time", -1), index * dt, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError("wrong diagnostic timestep")
    counts = conserved_counts(rows, names)
    cell_count = math.prod(system["NumCells" + axis] for axis in "XYZ")
    for entry in particles["particles"]:
        expected = cell_count * entry["NumPartPerCell"]
        if counts[entry["Name"]] != expected:
            raise ValueError(f"{entry['Name']}: initial count does not match the prescribed scenario")
    if system["Scheme"] == "symplectic":
        return {"counts": counts, **check_symplectic_invariants(rows, symplectic_energy_tol, continuity_tol)}
    ok, message = check_energy_conservation(str(directory), tol=energy_tol)
    if not ok:
        raise ValueError(message)
    result = {"counts": counts, "energy": message, "energy_limit": energy_tol}
    if system["Scheme"] == "ecsim_corr":
        result["continuity"] = check_continuity(directory, steps, dt, continuity_tol)
    return result


def state_arrays(directory, particles):
    """Optional bridge from application recovery files to explicit state arrays."""
    arrays = {}
    for name in ("FieldE", "FieldB"):
        shape, values = read_field(directory / "Recovery/Fields" / (name + ".backup"))
        arrays[name] = {"shape": [*shape, 3], "values": values, "components": 1}
    for entry in particles["particles"]:
        name = entry["Name"]
        count, values = read_particles(directory / "Recovery/Particles" / name / (name + ".backup"))
        arrays[name] = {"shape": [count, 6], "values": values, "components": 6,
                        "columns": ["x", "y", "z", "vx", "vy", "vz"]}
    return arrays


def component_errors(actual, expected, components):
    if (not actual or len(actual) != len(expected) or components <= 0 or len(actual) % components != 0):
        raise ValueError("state sizes differ or component layout is invalid")
    if not all(math.isfinite(v) for v in actual) or not all(math.isfinite(v) for v in expected):
        raise ValueError("non-finite state")
    return [max(abs(a-b) for a, b in zip(actual[c::components], expected[c::components]))
            for c in range(components)]


def write_payload(path, values, shape):
    if len(values) != math.prod(shape):
        raise ValueError("state shape and payload differ")
    finite_values(values, str(path))
    data = array("d", values)
    if sys.byteorder != "little":
        data.byteswap()
    payload = data.tobytes()
    path.write_bytes(payload)
    return {"file": path.name, "dtype": "f64le", "order": "C", "shape": shape,
            "sha256": hashlib.sha256(payload).hexdigest()}


def read_payload(directory, descriptor):
    name, shape = descriptor["file"], descriptor["shape"]
    if Path(name).name != name or (descriptor["dtype"] != "f64le" or descriptor.get("order") != "C") or not shape or min(shape) <= 0:
        raise ValueError("invalid reference descriptor")
    payload = (directory / name).read_bytes()
    if len(payload) != 8 * math.prod(shape) or hashlib.sha256(payload).hexdigest() != descriptor["sha256"]:
        raise ValueError("reference payload is truncated or changed")
    data = array("d")
    data.frombytes(payload)
    if sys.byteorder != "little":
        data.byteswap()
    return finite_values(data, name)


def reference_input(system, particles, mixed_precision, threads):
    # Outputs do not change the physical initial-value problem.
    physical = {k: v for k, v in system.items() if k not in ("diagnostics", "VerboseStep", "RecoveryInterval")}
    encoded = json.dumps({"system": physical, "particles": particles, "mixed_precision": mixed_precision,
                          "threads": threads}, sort_keys=True).encode()
    return hashlib.sha256(encoded).hexdigest()


def save_reference(directory, arrays, noise, input_id, producer, factor, repeats, scenario=None):
    directory.mkdir(parents=True, exist_ok=False)
    manifest = {"format": "ecsim-state-v1", "input_sha256": input_id, "producer": producer,
                "noise_factor": factor, "calibration_runs": repeats, "scenario": scenario, "arrays": {}}
    for name, item in arrays.items():
        components, values = item["components"], item["values"]
        scales = [max(abs(v) for v in values[c::components]) for c in range(components)]
        limits = [factor * max(n, math.ulp(scale)) for n, scale in zip(noise[name], scales)]
        descriptor = write_payload(directory / (name + ".f64"), values, item["shape"])
        descriptor.update(components=components, columns=item.get("columns"), noise_max_abs=noise[name],
                          scale=scales, abs_limits=limits)
        manifest["arrays"][name] = descriptor
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False))


def compare_reference(directory, arrays, input_id):
    manifest = json.loads((directory / "manifest.json").read_text())
    if manifest["format"] != "ecsim-state-v1" or manifest["input_sha256"] != input_id:
        raise ValueError("reference belongs to another scenario, solver mode or thread count")
    if arrays.keys() != manifest["arrays"].keys():
        raise ValueError("different fields or species in reference")
    result = {}
    for name, item in arrays.items():
        descriptor = manifest["arrays"][name]
        if item["shape"] != descriptor["shape"]:
            raise ValueError(f"{name}: reference shape differs")
        limits = descriptor["abs_limits"]
        if (descriptor["components"] != item["components"] or len(limits) != item["components"] or
                descriptor.get("columns") != item.get("columns") or
                not all(math.isfinite(limit) and limit >= 0 for limit in limits)):
            raise ValueError("invalid reference tolerances or component layout")
        expected = read_payload(directory, descriptor)
        errors = component_errors(item["values"], expected, item["components"])
        if any(error > limit for error, limit in zip(errors, descriptor["abs_limits"])):
            raise ValueError(f"{name}: state errors {errors} exceed calibrated limits {descriptor['abs_limits']}")
        result[name] = {"max_abs_errors": errors, "limits": descriptor["abs_limits"]}
    return result


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
                info["source_head_at_capture"] = revision.stdout.strip()
                info["source_diff_sha256_at_capture"] = hashlib.sha256(diff.stdout).hexdigest()
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
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, help="new directory; default creates a temporary result directory")
    parser.add_argument("--profile", choices=PROFILES, default="smoke")
    parser.add_argument("--case", choices=("all", "periodic", "magnetized", "two_species"), default="all")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--schemes", nargs="+", choices=("ecsim", "ecsim_corr", "symplectic"),
                        default=["ecsim", "ecsim_corr", "symplectic"])
    parser.add_argument("--cells", nargs=3, type=int)
    parser.add_argument("--npc", type=int)
    parser.add_argument("--steps", type=int)
    parser.add_argument("--dt", type=float, default=0.1)
    parser.add_argument("--mixed-precision", action="store_true")
    parser.add_argument("--energy-tol", type=float, default=1e-8)
    parser.add_argument("--continuity-tol", type=float, default=1e-10)
    parser.add_argument("--symplectic-energy-tol", type=float, default=1e-3)
    parser.add_argument("--timeout", type=float, default=3600)
    parser.add_argument("--save-reference", type=Path, help="create a separate state artifact; never overwrites")
    parser.add_argument("--reference", type=Path, help="optionally compare complete states with a saved reference")
    parser.add_argument("--calibration-runs", type=int, default=3)
    parser.add_argument("--noise-factor", type=float, default=5.0)
    args = parser.parse_args()
    if args.threads <= 0 or args.threads > len(os.sched_getaffinity(0)):
        parser.error("threads must fit the available CPU affinity")
    if any(v is not None and v <= 0 for v in (args.npc, args.steps)) or (args.cells and min(args.cells) <= 0):
        parser.error("cell counts, npc and steps must be positive")
    if not all(math.isfinite(v) and v > 0 for v in (args.dt, args.energy_tol, args.continuity_tol,
                                                   args.symplectic_energy_tol, args.timeout)):
        parser.error("dt, energy tolerance and timeout must be finite and positive")
    if args.calibration_runs < 2 or not 2 <= args.noise_factor <= 10:
        parser.error("reference calibration needs >=2 runs and a safety factor in [2,10]")
    args.binary = args.binary.resolve()
    if not args.binary.is_file() or not os.access(args.binary, os.X_OK):
        parser.error("binary must be an executable file")
    cells, npc, steps = PROFILES[args.profile]
    PROFILES["requested"] = (tuple(args.cells or cells), args.npc or npc, args.steps or steps)
    cases = ["periodic", "magnetized", "two_species"] if args.case == "all" else [args.case]
    if args.output:
        args.output = args.output.resolve()
        args.output.mkdir(parents=True, exist_ok=False)
    else:
        args.output = Path(tempfile.mkdtemp(prefix="ecsim-physics-"))
    print("Results:", args.output, flush=True)
    if args.save_reference:
        args.save_reference = args.save_reference.resolve()
        args.save_reference.mkdir(parents=True, exist_ok=False)
    if args.reference:
        args.reference = args.reference.resolve()
    binary = args.output / "program"
    producer = snapshot_binary(args.binary, binary)
    report = {"binary": producer, "profile": args.profile, "threads": args.threads,
              "mixed_precision": args.mixed_precision, "results": []}
    try:
        for case in cases:
            for scheme in args.schemes:
                system, particles = configuration(scheme, "requested", case)
                system["Dt"] = args.dt
                input_id = reference_input(system, particles, args.mixed_precision, args.threads)
                needs_state = bool(args.save_reference or args.reference)
                if needs_state:
                    system["diagnostics"]["outputs"].append({"type": "recovery", "interval": system["LastTimestep"]})
                first = None
                noise = {}
                repeats = args.calibration_runs if args.save_reference else 1
                result = {"case": case, "scheme": scheme, "input_sha256": input_id, "runs": []}
                for repeat in range(repeats):
                    directory = args.output / case / scheme / str(repeat)
                    execute(binary, directory, system, particles, args.threads, args.timeout, args.mixed_precision)
                    result["runs"].append(validate_run(directory, system, particles, args.energy_tol, args.continuity_tol,
                                                      args.symplectic_energy_tol))
                    if needs_state:
                        current = state_arrays(directory, particles)
                        if args.reference:
                            result["reference"] = compare_reference(args.reference / case / scheme, current, input_id)
                        if first is None:
                            first = current
                            noise = {name: [0.0] * item["components"] for name, item in current.items()}
                        else:
                            for name, item in current.items():
                                errors = component_errors(item["values"], first[name]["values"], item["components"])
                                noise[name] = [max(old, new) for old, new in zip(noise[name], errors)]
                if args.save_reference:
                    save_reference(args.save_reference / case / scheme, first, noise, input_id, producer,
                                   args.noise_factor, repeats, {"system": system, "particles": particles,
                                                               "threads": args.threads, "mixed_precision": args.mixed_precision})
                    result["calibration_noise"] = noise
                report["results"].append(result)
                print(f"{case}/{scheme}: physics PASS" + (", reference PASS" if args.reference else ""), flush=True)
        report["status"] = "passed"
    except Exception as error:
        report.update(status="failed", error=str(error))
        raise
    finally:
        (args.output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

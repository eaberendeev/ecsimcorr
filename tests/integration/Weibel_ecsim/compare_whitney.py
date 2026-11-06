#!/usr/bin/env python3
"""Explicit, isolated Weibel comparison. Never removes an existing run directory.

Uses the existing Weibel configuration, changing only particle count, scheme,
predictor, and output frequency. Run the CMake target weibel_whitney first.
"""
import argparse
import copy
from concurrent.futures import ThreadPoolExecutor, as_completed
import importlib.util
import json
import os
from pathlib import Path
from queue import SimpleQueue
import subprocess
import time


VARIANTS = {
    "ecsim": ("ecsim", None),
    "ecsim_corr": ("ecsim_corr", None),
    "point": ("ecsim_whitney", None),
    "ballistic": ("ecsim_whitney_frozen", "ballistic"),
    "old_fields": ("ecsim_whitney_frozen", "old_fields"),
    "half_drift": ("ecsim_whitney_frozen", "half_drift"),
    "centered_fields": ("ecsim_whitney_frozen", "centered_fields"),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--ppc", type=int, default=16)
    parser.add_argument("--steps", type=int)
    parser.add_argument("--profile-step", type=int, default=0, help="write timer profiles for one step (0 disables)")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--jobs", type=int, default=1, help="independent runs in parallel")
    parser.add_argument("--assembly", choices=("reference", "cell"), default="reference")
    parser.add_argument("--variants", nargs="+", choices=VARIANTS, default=list(VARIANTS))
    args = parser.parse_args()
    if min(args.ppc, args.threads, args.jobs) < 1 or (args.steps is not None and args.steps < 1):
        parser.error("counts must be positive")
    binary = args.binary.resolve(strict=True)
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=False)
    config_path = Path(__file__).with_name("gen_config.py")
    spec = importlib.util.spec_from_file_location("weibel_config", config_path)
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    steps = args.steps or source.system_config["LastTimestep"]
    manifest = {
        "source_config": str(config_path), "binary": str(binary),
        "ppc": args.ppc, "original_ppc": source.NumPartPerCell,
        "steps": steps, "dt": source.Dt, "threads": args.threads, "jobs": args.jobs,
        "variants": args.variants, "assembly": args.assembly, "results": {},
    }
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), OMP_PROC_BIND="close",
               VALIDATION_PERIODICITY="1000000", TIMER_START_STEP="2147483647")
    fingerprints = set()
    # Independent OpenMP processes otherwise all bind to the first places.
    # Prefer separate physical cores before using SMT siblings.
    allowed = sorted(os.sched_getaffinity(0))
    physical, seen = [], set()
    for cpu in allowed:
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        identity = ((topology / "physical_package_id").read_text().strip(),
                    (topology / "core_id").read_text().strip())
        if identity not in seen:
            physical.append(cpu)
            seen.add(identity)
    cpus = physical + [cpu for cpu in allowed if cpu not in physical]
    if len(cpus) < args.jobs * args.threads:
        parser.error("not enough available CPUs for jobs*threads")
    places = SimpleQueue()
    for slot in range(args.jobs):
        places.put(cpus[slot * args.threads:(slot + 1) * args.threads])
    def run_one(name):
        assigned = places.get()
        scheme, predictor = VARIANTS[name]
        run = outdir / name
        run.mkdir()
        system = copy.deepcopy(source.system_config)
        particles = copy.deepcopy(source.particles_config)
        system.update(Scheme=scheme, NumPartPerCell=args.ppc, LastTimestep=steps,
                      ComparisonInterval=6)
        if args.profile_step:
            system["ComparisonProfileStep"] = args.profile_step
        system["diagnostics"]["outputs"] = [
            {"type": "energy_balance", "interval": 1},
            {"type": "fields_2d", "interval": 12, "fields": ["E", "B"],
             "planes": {"z": [source.bbox_centerZ]}},
            {"type": "density_2d", "interval": 12, "species": "all",
             "planes": {"z": [source.bbox_centerZ]}},
        ]
        if predictor:
            system["WhitneyOrbit"] = {"predictor": predictor, "assembly": args.assembly}
        elif name in ("point", "corr"):
            system["WhitneyAssembly"] = args.assembly
        for species in particles["particles"]:
            species["NumPartPerCell"] = args.ppc
        for filename, data in (("system_config.json", system), ("particles_config.json", particles)):
            (run / filename).write_text(json.dumps(data, indent=2) + "\n")
        (run / "phys.par").write_text(f"w_p = {source.w_p}\n1/w_p = {1 / source.w_p}\n")
        print(f"START {name}: {steps} steps, ppc={args.ppc}, threads={args.threads}", flush=True)
        start = time.monotonic()
        child_env = dict(env, OMP_PLACES=",".join("{" + str(cpu) + "}" for cpu in assigned))
        with (run / "run.log").open("w") as log:
            result = subprocess.run([str(binary)], cwd=run, env=child_env, stdout=log, stderr=subprocess.STDOUT)
        places.put(assigned)
        elapsed = time.monotonic() - start
        state_path = run / "initial_state.json"
        state = json.loads(state_path.read_text()) if state_path.exists() else None
        print(f"END {name}: exit={result.returncode}, seconds={elapsed:.2f}", flush=True)
        return name, {"exit_code": result.returncode, "seconds": elapsed, "initial_state": state, "cpus": assigned}
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(run_one, name) for name in args.variants]
        for future in as_completed(futures):
            name, result = future.result()
            manifest["results"][name] = result
            (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
            state = result["initial_state"]
            if state:
                fingerprints.add((state["fingerprint"], state["particles"]))
            if result["exit_code"]:
                raise RuntimeError(f"{name} failed; inspect {outdir / name / 'run.log'}")
            if len(fingerprints) != 1:
                raise RuntimeError("initial particle fingerprints differ between schemes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

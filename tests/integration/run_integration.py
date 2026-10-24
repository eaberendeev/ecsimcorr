#!/usr/bin/env python3
"""Integration test runner for beren3d.

Usage: python3 run_integration.py [--test TEST_NAME] [--keep] [--timesteps N]

Each test lives in tests/integration/<name>/gen_config.py. The runner:
1. Runs build.py --config <test>/gen_config.py (builds + generates config + creates workdir)
2. Runs beren3d inside the workdir
3. Collects diagnostics and checks results
4. Reports pass/fail
"""

import os
import shutil
import subprocess
import sys
import re
import glob
import argparse
import json
import math


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
WORKDIR_FILE = os.path.join(PROJECT_ROOT, "workdir.tmp")


def run(cmd, cwd=None, capture=False):
    """Run a shell command, optionally capture stdout/stderr."""
    print(f"  + {cmd}")
    if capture:
        result = subprocess.run(cmd, shell=True, cwd=cwd,
                                capture_output=True, text=True)
    else:
        result = subprocess.run(cmd, shell=True, cwd=cwd)
    if result.returncode != 0:
        print(f"  ERROR: command failed with code {result.returncode}")
        if capture:
            print(f"  stdout:\n{result.stdout}")
            print(f"  stderr:\n{result.stderr}")
        sys.exit(1)
    return result


def parse_energy_file(workdir):
    """Read a complete numeric table; malformed or non-finite data is a failure."""
    fpath = os.path.join(workdir, "energy.txt")
    try:
        with open(fpath) as stream:
            header = stream.readline().split()
            if not header or len(header) != len(set(header)):
                return None, "Empty or duplicate energy columns"
            records = []
            for lineno, line in enumerate(stream, 2):
                if not line.strip():
                    continue
                parts = line.split()
                if len(parts) != len(header):
                    return None, f"Malformed energy row {lineno}: wrong column count"
                values = [float(value) for value in parts]
                if not all(math.isfinite(value) for value in values):
                    return None, f"NaN/Inf in energy row {lineno}"
                records.append(dict(zip(header, values)))
    except (OSError, ValueError) as error:
        return None, f"Cannot read energy data: {error}"
    return records, None


def check_energy_conservation(workdir, tol=1e-3):
    records, error = parse_energy_file(workdir)
    if error or not records:
        return False, error or "Empty energy data"
    if not math.isfinite(tol) or tol <= 0:
        return False, "Energy tolerance must be finite and positive"
    if "energyConserve" not in records[0]:
        return False, "No 'energyConserve' column in energy.txt"
    try:
        with open(os.path.join(workdir, "particles_config.json")) as stream:
            names = [entry["Name"] for entry in json.load(stream)["particles"]]
    except (OSError, KeyError, TypeError, ValueError) as error:
        return False, f"Cannot resolve species energies: {error}"
    if any(name not in records[0] for name in names):
        return False, "Missing species kinetic energy column"

    config_path = os.path.join(workdir, "system_config.json")
    if os.path.isfile(config_path):
        try:
            with open(config_path) as stream:
                config = json.load(stream)
            if "LastTimestep" in config and "Dt" in config:
                interval = next((output.get("interval", 1)
                                 for output in config.get("diagnostics", {}).get("outputs", [])
                                 if output.get("type") == "energy_balance"), 1)
                expected_time = (config["LastTimestep"] // interval) * interval * config["Dt"]
                if not math.isclose(records[-1].get("Time", float("nan")), expected_time,
                                    rel_tol=1e-10, abs_tol=1e-12):
                    return False, "Energy table did not reach the last requested diagnostic timestep"
        except (OSError, ValueError, TypeError, ZeroDivisionError) as error:
            return False, f"Cannot validate diagnostic completion: {error}"

    def total_energy(record):
        return (sum(record[name] for name in names) + record.get("energyFieldE", 0.0)
                + record.get("energyFieldB", 0.0))

    initial_energy = total_energy(records[0])
    max_error = max(abs(record["energyConserve"]) /
                    max(abs(total_energy(record)), abs(initial_energy), 1e-30)
                    for record in records)
    if max_error > tol:
        return False, f"Energy error {max_error:.2e} exceeds tolerance {tol:.2e}"
    return True, f"Energy conserved, max relative error = {max_error:.2e}"


def check_nan_inf(workdir):
    records, error = parse_energy_file(workdir)
    return (False, error or "Empty energy data") if error or not records else (True, "No NaN/Inf")


def run_test_case(test_name, clean_workdir=False, extra_timesteps=None):
    """Run a single integration test case."""
    test_dir = os.path.join(PROJECT_ROOT, "tests", "integration", test_name)
    test_config = os.path.join(test_dir, "gen_config.py")

    if not os.path.isfile(test_config):
        print(f"ERROR: test config not found: {test_config}")
        return False

    print(f"\n{'='*60}")
    print(f"Integration test: {test_name}")
    print(f"{'='*60}")

    # Copy test config to a temp file; apply --timesteps override if needed
    config_to_use = os.path.join(PROJECT_ROOT, f"__test_{test_name}.py")
    shutil.copy(test_config, config_to_use)
    if extra_timesteps:
        with open(config_to_use) as f:
            content = f.read()
        content = re.sub(r'LastTimestep\s*=.*',
                         f'LastTimestep = {extra_timesteps}', content)
        with open(config_to_use, 'w') as f:
            f.write(content)

    build_args = f"--rerun --type Release --timers 0 --config '{config_to_use}'"
    if not os.path.isdir(os.path.join(PROJECT_ROOT, "_build_release")):
        build_args = "--rebuild " + build_args
    # Eigen and AMGCL are vendored in srcBeren/3rd_party, so build.py takes no
    # external library paths.
    run(f"python3 build.py {build_args}", cwd=PROJECT_ROOT)

    if not os.path.isfile(WORKDIR_FILE):
        print("ERROR: workdir.tmp not found after build")
        return False
    with open(WORKDIR_FILE) as f:
        workdir = f.readline().strip()
        numprocs = f.readline().strip() or "1"
    os.remove(WORKDIR_FILE)
    print(f"  Workdir: {workdir}")
    print(f"  Threads (NumProcs from config): {numprocs}")

    binary = os.path.join(workdir, "beren3d")
    if not os.path.isfile(binary):
        print(f"ERROR: binary not found: {binary}")
        return False

    run(f"OMP_NUM_THREADS={numprocs} numactl --interleave=all ./beren3d",
        cwd=workdir, capture=True)

    ok = False
    ok_nan, msg_nan = check_nan_inf(workdir)
    if not ok_nan:
        print(f"  FAIL: {msg_nan}")
    else:
        ok, msg = check_energy_conservation(workdir)
        if ok:
            print(f"  PASS: {msg}")
        else:
            print(f"  FAIL: {msg}")

    check_script = os.path.join(test_dir, "check_weibel.py")
    if ok and os.path.isfile(check_script):
        print(f"  Weibel check (draws figures, checks growth rate...):")
        res = subprocess.run(f"python3 {check_script} {workdir}",
                             shell=True, cwd=PROJECT_ROOT,
                             capture_output=True, text=True)
        print(res.stdout, end="")
        if res.stderr:
            print(res.stderr, end="")
        if res.returncode == 0:
            print("  PASS: Weibel checks all passed")
        else:
            print(f"  FAIL: Weibel check exited with code {res.returncode}")
            ok = False

    if os.path.isfile(config_to_use):
        os.remove(config_to_use)

    if clean_workdir and ok:
        full_path = os.path.join(PROJECT_ROOT, workdir)
        if os.path.isdir(full_path):
            print(f"  Cleaning workdir: {full_path}")
            shutil.rmtree(full_path, ignore_errors=True)

    return ok


def main():
    parser = argparse.ArgumentParser(description="beren3d integration test runner")
    parser.add_argument("--test", default="all",
                        help="Test name (subdirectory of tests/integration/) or 'all'")
    parser.add_argument("--clean", action="store_true",
                        help="Remove workdir after a successful test (default: keep)")
    parser.add_argument("--timesteps", type=int, default=None,
                        help="Override LastTimestep (for debugging)")
    args = parser.parse_args()

    integration_dir = os.path.join(PROJECT_ROOT, "tests", "integration")
    available = sorted(d for d in os.listdir(integration_dir)
                       if os.path.isdir(os.path.join(integration_dir, d))
                       and d != '__pycache__'
                       and os.path.isfile(os.path.join(integration_dir, d, 'gen_config.py')))

    if not available:
        print("No integration tests found in", integration_dir)
        sys.exit(1)

    if args.test == "all":
        # A test config may declare `IntegrationHeavy = True` for production-scale
        # runs that cannot finish in CI; they stay selectable via --test NAME.
        def is_heavy(name):
            cfg_path = os.path.join(integration_dir, name, "gen_config.py")
            with open(cfg_path) as fh:
                return any(line.split("#")[0].strip() == "IntegrationHeavy = True"
                           for line in fh)

        heavy = [t for t in available if is_heavy(t)]
        tests = [t for t in available if t not in heavy]
        if heavy:
            print(f"Skipping heavy integration tests (run explicitly): {', '.join(heavy)}")
    else:
        if args.test not in available:
            print(f"Unknown test '{args.test}'. Available: {', '.join(available)}")
            sys.exit(1)
        tests = [args.test]

    passed = 0
    failed = 0
    for test in tests:
        if run_test_case(test, clean_workdir=args.clean, extra_timesteps=args.timesteps):
            passed += 1
        else:
            failed += 1

    print(f"\n{'='*60}")
    print(f"Results: {passed} passed, {failed} failed")
    print(f"{'='*60}")

    if failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()

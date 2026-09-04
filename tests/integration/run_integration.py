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
    """Parse energy.txt from workdir. Returns list of {key: value} dicts."""
    fpath = os.path.join(workdir, "energy.txt")
    if not os.path.isfile(fpath):
        return None, f"No energy.txt in {workdir}"

    records = []
    with open(fpath) as fh:
        header = fh.readline().strip().split()
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= len(header):
                rec = {}
                for i, h in enumerate(header):
                    try:
                        rec[h] = float(parts[i])
                    except (ValueError, IndexError):
                        rec[h] = parts[i]
                records.append(rec)
    return records, None


def check_energy_conservation(workdir, tol=1e-3):
    """Check energy conservation: |energyConserve| / total_kinetic_energy < tol.
    Reports max fractional energy error over all timesteps."""
    records, err = parse_energy_file(workdir)
    if err:
        return False, err

    if not records:
        return False, "Empty energy data"

    required = ['energyConserve']
    for key in required:
        if key not in records[0]:
            return False, f"No '{key}' column in energy.txt"

    species_keys = [k for k in records[0] if k in ('Electrons', 'Ions', 'Protons', 'Neutrals')]

    max_frac_err = 0.0
    for r in records:
        e_cons = r['energyConserve']
        e_kin = sum(r.get(k, 0.0) for k in species_keys)
        if abs(e_kin) > 1e-30:
            frac = abs(e_cons) / abs(e_kin)
            if frac > max_frac_err:
                max_frac_err = frac

    if max_frac_err > tol:
        return False, f"Energy error {max_frac_err:.2e} exceeds tolerance {tol:.2e}"

    return True, f"Energy conserved, max relative error = {max_frac_err:.2e}"


def check_nan_inf(workdir):
    """Check for NaN or Inf in energy.txt."""
    fpath = os.path.join(workdir, "energy.txt")
    if not os.path.isfile(fpath):
        return False, "No energy.txt"
    with open(fpath) as fh:
        for line in fh:
            if 'nan' in line.lower() or 'inf' in line.lower():
                return False, f"NaN/Inf found: {line.strip()}"
    return True, "No NaN/Inf"


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
    if not os.path.isdir(os.path.join(PROJECT_ROOT, "_build")):
        build_args = "--rebuild " + build_args
    run(f"python3 build.py {build_args} "
        f"--eigen '{eigen_path}' --amgcl '{amgcl_path}'",
        cwd=PROJECT_ROOT)

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
        tests = available
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

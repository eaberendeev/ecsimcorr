#!/usr/bin/env python3
"""Integration test runner for beren3d.

Usage: python3 run_integration.py [--test TEST_NAME] [--keep] [--timesteps N]

Each test lives in tests/integration/<name>/gen_config.py. The runner:
1. Copies the test's gen_config.py to the project root
2. Runs build.py --rebuild --rerun  (builds + generates config + creates workdir)
3. Runs beren3d inside the workdir
4. Collects diagnostics and checks results
5. Reports pass/fail
"""

import json
import os
import shutil
import subprocess
import sys
import re
import glob
import argparse
from pathlib import Path


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

    required = ['energyConserve', 'Ions', 'Electrons']
    for key in required:
        if key not in records[0]:
            return False, f"No '{key}' column in energy.txt"

    max_frac_err = 0.0
    for r in records:
        e_cons = r['energyConserve']
        e_kin = r['Ions'] + r['Electrons']
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


def run_test_case(test_name, keep_workdir=False, extra_timesteps=None):
    """Run a single integration test case."""
    test_dir = os.path.join(PROJECT_ROOT, "tests", "integration", test_name)
    test_config = os.path.join(test_dir, "gen_config.py")

    if not os.path.isfile(test_config):
        print(f"ERROR: test config not found: {test_config}")
        return False

    print(f"\n{'='*60}")
    print(f"Integration test: {test_name}")
    print(f"{'='*60}")

    orig_config = os.path.join(PROJECT_ROOT, "gen_config.py")
    backup_config = os.path.join(PROJECT_ROOT, "gen_config.py.bak")

    if os.path.isfile(orig_config):
        shutil.copy(orig_config, backup_config)

    ok = False
    try:
        shutil.copy(test_config, orig_config)

        if extra_timesteps:
            with open(orig_config) as f:
                content = f.read()
            content = re.sub(r'LastTimestep\s*=\s*\d+',
                             f'LastTimestep = {extra_timesteps}', content)
            with open(orig_config, 'w') as f:
                f.write(content)

        eigen_path = os.environ.get("EIGEN_PATH", os.path.expanduser("~/soft/eigen-3.4.0/"))
        amgcl_path = os.environ.get("AMGCL_PATH", os.path.expanduser("~/soft/amgcl/"))
        build_args = "--rerun --type Release --timers 0"
        if not os.path.isdir(os.path.join(PROJECT_ROOT, "_build")):
            build_args = "--rebuild " + build_args
        run(f"python3 build.py {build_args} "
            f"--eigen '{eigen_path}' --amgcl '{amgcl_path}'",
            cwd=PROJECT_ROOT)

        if not os.path.isfile(WORKDIR_FILE):
            print("ERROR: workdir.tmp not found after build")
            return False
        with open(WORKDIR_FILE) as f:
            workdir = f.read().strip()
        print(f"  Workdir: {workdir}")

        binary = os.path.join(workdir, "beren3d")
        if not os.path.isfile(binary):
            print(f"ERROR: binary not found: {binary}")
            return False

        run(f"OMP_NUM_THREADS=1 numactl --interleave=all ./beren3d",
            cwd=workdir, capture=True)

        ok_nan, msg_nan = check_nan_inf(workdir)
        if not ok_nan:
            print(f"  FAIL: {msg_nan}")
            ok = False
        else:
            ok, msg = check_energy_conservation(workdir)
            if ok:
                print(f"  PASS: {msg}")
            else:
                print(f"  FAIL: {msg}")

        return ok

    finally:
        if os.path.isfile(backup_config):
            shutil.copy(backup_config, orig_config)
            os.remove(backup_config)

        if os.path.isfile(WORKDIR_FILE):
            os.remove(WORKDIR_FILE)

        if not keep_workdir and ok:
            for d in glob.glob(os.path.join(PROJECT_ROOT, "Res_*")):
                shutil.rmtree(d, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(description="beren3d integration test runner")
    parser.add_argument("--test", default="all",
                        help="Test name (subdirectory of tests/integration/) or 'all'")
    parser.add_argument("--keep", action="store_true",
                        help="Keep workdir after test")
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
        if run_test_case(test, keep_workdir=args.keep, extra_timesteps=args.timesteps):
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

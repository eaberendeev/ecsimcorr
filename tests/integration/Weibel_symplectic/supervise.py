#!/usr/bin/env python3
"""Prepare and supervise the full ECSIM-sized Weibel test without deleting runs."""
import argparse
import datetime
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time


def save_status(directory, status):
    status["updated_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    temporary = directory / "status.json.tmp"
    temporary.write_text(json.dumps(status, indent=2) + "\n")
    temporary.replace(directory / "status.json")


def energy_rows(directory):
    path = directory / "energy.txt"
    if not path.exists():
        return []
    text = path.read_text()
    lines = text.splitlines()
    if not text.endswith("\n"):
        lines = lines[:-1]
    if len(lines) < 2:
        return []
    header = lines[0].split()
    rows = []
    for line in lines[1:]:
        values = line.split()
        if len(values) != len(header):
            raise ValueError("malformed energy row")
        row = dict(zip(header, map(float, values)))
        if not all(math.isfinite(v) for v in row.values()):
            raise ValueError("non-finite diagnostic")
        rows.append(row)
    return rows


def stop(process):
    os.killpg(process.pid, signal.SIGTERM)
    try:
        process.wait(timeout=15)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=43200, help="wall-time limit in seconds")
    parser.add_argument("--poll-seconds", type=float, default=30)
    args = parser.parse_args()
    args.binary = args.binary.resolve()
    directory = args.directory.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    config_path = Path(__file__).with_name("gen_config.py")
    spec = importlib.util.spec_from_file_location("full_weibel_config", config_path)
    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)
    os.chdir(directory)
    config.generate_config()
    shutil.copy2(args.binary.resolve(), directory / "beren3d")
    shutil.copy2(config_path, directory / "source_config.py")
    shutil.copy2(config.REFERENCE_CONFIG, directory / "reference_ecsim_config.py")
    binary_hash = hashlib.sha256((directory / "beren3d").read_bytes()).hexdigest()
    system = json.loads((directory / "system_config.json").read_text())
    dt, end = system["Dt"], system["LastTimestep"] * system["Dt"]
    status = {"state": "initializing", "supervisor_pid": os.getpid(), "threads": config.NumProcs,
              "dt": dt, "end_time": end, "last_step": system["LastTimestep"],
              "timeout_seconds": args.timeout, "binary_sha256": binary_hash,
              "reference_config": str(config.REFERENCE_CONFIG)}
    save_status(directory, status)
    started = time.monotonic()
    env = dict(os.environ, OMP_NUM_THREADS=str(config.NumProcs), OMP_DYNAMIC="FALSE")
    with (directory / "run.log").open("w") as log:
        process = subprocess.Popen([str(directory / "beren3d")], stdout=log, stderr=subprocess.STDOUT,
                                   env=env, start_new_session=True)
        status.update(state="running", simulation_pid=process.pid)
        save_status(directory, status)
        reason = None
        while process.poll() is None:
            elapsed = time.monotonic() - started
            try:
                rows = energy_rows(directory)
                if rows:
                    initial = rows[0]["energyTotal"]
                    if initial <= 0:
                        raise ValueError("nonpositive initial total energy")
                    error = max(abs(r["energyError"]) / initial for r in rows)
                    gauss = max(abs(r["gaussDrift"]) for r in rows)
                    div_b = max(abs(r["divB"]) for r in rows)
                    if error > 0.05 or gauss > 1e-8 or div_b > 1e-8:
                        raise ValueError("monitor energy/constraint guard exceeded")
                    last_time = rows[-1]["Time"]
                    last_step = round(last_time / dt)
                    status.update(time=last_time, step=last_step, fraction_complete=last_time/end,
                                  max_relative_energy_error=error, max_gauss_drift=gauss, max_div_b=div_b,
                                  magnetic_energy=rows[-1]["energyFieldB"], electric_energy=rows[-1]["energyFieldE"],
                                  initial_gauss_error=rows[0]["gaussError"],
                                  checkpoint_step=max(0, (last_step-1)//500*500),
                                  estimated_seconds_remaining=elapsed*(end/last_time-1) if last_time else None)
                status["wall_seconds"] = elapsed
                save_status(directory, status)
                print(json.dumps(status), flush=True)
            except (ValueError, OSError) as error:
                reason = str(error)
            if elapsed > args.timeout:
                reason = "wall-time limit reached"
            if reason:
                stop(process)
                break
            time.sleep(args.poll_seconds)
        status["simulation_exit_code"] = process.wait()
    try:
        final_rows = energy_rows(directory)
    except (ValueError, OSError) as error:
        reason = str(error)
        final_rows = []
    if final_rows:
        final = final_rows[-1]
        status.update(time=final["Time"], step=round(final["Time"]/dt), fraction_complete=final["Time"]/end,
                      magnetic_energy=final["energyFieldB"], electric_energy=final["energyFieldE"])
    if reason or status["simulation_exit_code"] != 0:
        status.update(state="failed", reason=reason or "simulation exited unsuccessfully")
        save_status(directory, status)
        return 1
    status.update(state="analyzing", wall_seconds=time.monotonic()-started)
    save_status(directory, status)
    checker = Path(__file__).resolve().with_name("check_weibel.py")
    try:
        with (directory / "analysis.log").open("w") as log:
            result = subprocess.run([sys.executable, str(checker), str(directory), "--quiet"],
                                    stdout=log, stderr=subprocess.STDOUT, timeout=300)
    except (OSError, subprocess.TimeoutExpired) as error:
        status.update(state="failed_analysis", reason=str(error))
        save_status(directory, status)
        return 1
    status.update(state="passed" if result.returncode == 0 else "failed_validation",
                  analysis_exit_code=result.returncode, wall_seconds=time.monotonic()-started)
    save_status(directory, status)
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Checks the reduced hybrid coupled-field regression.

The runner calls this file for every test that has it.  For this case the
physical requirement is a finite, particle-conserving, energy-balanced coupled
step; growth-rate fitting is intentionally left to the full Weibel tests.
"""

import json
import math
import sys
from pathlib import Path


def fail(message):
    print(f"hybrid_coupled: {message}", file=sys.stderr)
    return 2


def main():
    workdir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    try:
        config = json.loads((workdir / "system_config.json").read_text())
        particles = json.loads((workdir / "particles_config.json").read_text())["particles"]
        names = [entry["Name"] for entry in particles]
        header = (workdir / "energy.txt").read_text().splitlines()[0].split()
        rows = []
        for line in (workdir / "energy.txt").read_text().splitlines()[1:]:
            if line.strip():
                values = [float(value) for value in line.split()]
                if len(values) != len(header):
                    return fail("malformed energy row")
                rows.append(dict(zip(header, values)))
    except (OSError, KeyError, TypeError, ValueError, IndexError) as error:
        return fail(f"cannot read run data: {error}")

    if config.get("Scheme") != "hybrid_kinetic":
        return fail("expected Scheme=hybrid_kinetic")
    if not rows or len(rows) < 2:
        return fail("empty energy table")
    required = ["Time", "energyConserve", "nonlinearResidual", "nonlinearIterations",
                "energyFieldE", "energyFieldB", "guidingCenterEnergyResidual"]
    missing = [name for name in required if name not in rows[0]]
    if missing:
        return fail(f"missing columns: {missing}")
    counted_species = [name for name in names if f"{name}Particles" in rows[0]]
    if not counted_species:
        return fail("no species particle-count column in energy table")
    required += [f"{name}Particles" for name in counted_species]
    if any(not all(math.isfinite(row[name]) for name in required) for row in rows):
        return fail("NaN/Inf in energy table")

    expected_time = config["LastTimestep"] * config["Dt"]
    if not math.isclose(rows[-1]["Time"], expected_time, rel_tol=1e-12, abs_tol=1e-12):
        return fail("energy table did not reach the requested time")
    for name in counted_species:
        counts = [row[f"{name}Particles"] for row in rows]
        if any(count != counts[0] for count in counts):
            return fail(f"particle count changed for {name}")

    initial_energy = sum(rows[0][name] for name in counted_species) + rows[0]["energyFieldE"] + rows[0]["energyFieldB"]
    if not (initial_energy > 0 and math.isfinite(initial_energy)):
        return fail("invalid initial energy")
    max_balance = max(abs(row["energyConserve"]) for row in rows) / initial_energy
    max_residual = max(abs(row["nonlinearResidual"]) for row in rows)
    max_theorem = max(abs(row["guidingCenterEnergyResidual"]) for row in rows) / initial_energy
    if max_theorem > 1e-10:
        return fail(f"guiding-centre energy theorem violated: {max_theorem:.3e}")
    max_iterations = max(row["nonlinearIterations"] for row in rows)
    # The scheme conserves energy to the nonlinear tolerance (1e-10).
    if max_balance > 1e-8:
        return fail(f"energy balance too large: {max_balance:.3e}")
    if max_residual > 1e-8:
        return fail(f"nonlinear residual too large: {max_residual:.3e}")
    if max_iterations < 1:
        return fail("coupled solver did not iterate")

    report = {
        "scheme": config["Scheme"],
        "steps": config["LastTimestep"],
        "final_time": rows[-1]["Time"],
        "particles": {name: int(rows[0][f"{name}Particles"]) for name in counted_species},
        "max_relative_energy_balance": max_balance,
        "max_nonlinear_residual": max_residual,
        "max_relative_guiding_center_theorem": max_theorem,
        "max_nonlinear_iterations": int(max_iterations),
        "passed": True,
    }
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

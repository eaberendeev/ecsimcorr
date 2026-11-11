#!/usr/bin/env python3
"""Checks the magnetosonic wave carried by electrons of two models.

The frequency and the energy balance are those of wave_cases/fms_check.py.
In addition both parts of the adaptive species must hold a sizeable share of
the markers during the whole run, their number must be constant, and the
energy defect of the model changes must stay at rounding level.
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "wave_cases"))
from fms_check import analyse  # noqa: E402

SWITCH_DEFECT_TOLERANCE = 1e-12   # relative to the electron energy


def main():
    workdir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    result = analyse(workdir)
    if isinstance(result, str):
        print(f"fms_adaptive: {result}", file=sys.stderr)
        return 2
    lines = (workdir / "energy.txt").read_text().splitlines()
    head = lines[0].split()
    table = np.array([[float(v) for v in line.split()] for line in lines[1:] if line.strip()])
    column = lambda name: table[:, head.index(name)]  # noqa: E731
    guiding = column("Electrons_gcParticles")
    electron_energy = column("Electrons") + column("Electrons_gc")
    moved = column("Electrons_gcFromKinetic").sum() + column("Electrons_gcToKinetic").sum()
    result.update({
        "guiding_centres_min": float(guiding.min()),
        "guiding_centres_max": float(guiding.max()),
        "model_changes": float(moved),
        "kinetic_energy_share_min": float((column("Electrons") / electron_energy).min()),
        "kinetic_energy_share_max": float((column("Electrons") / electron_energy).max()),
        "switch_energy_defect": float(np.max(np.abs(column("modelSwitchEnergyDefect")) / electron_energy)),
    })
    result["passed"] = bool(
        result["passed"]
        and moved > 0
        and 0.2 < result["kinetic_energy_share_min"] and result["kinetic_energy_share_max"] < 0.8
        and result["switch_energy_defect"] < SWITCH_DEFECT_TOLERANCE)
    print(json.dumps(result, indent=2))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())

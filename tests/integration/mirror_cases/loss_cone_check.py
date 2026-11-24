"""Analysis of a loss-cone run (loss_cone.py), common to all schemes.

An isotropic electron at a point where the field is B0 escapes the mirror if
its pitch angle lies in the loss cone, sin^2 < B0 / Bmax: a fraction
1 - sqrt(1 - B0 / Bmax). The cloud is short, so the fraction is averaged over
its length with the on-axis coil field. The run must reproduce it (escaped
markers counted from the electron count of energy.txt) and keep the trapped
electrons: a scheme that loses the magnetic moment scatters them into the
cone, and they keep leaking after the first transit.
"""

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import loss_cone  # noqa: E402

FRACTION_TOLERANCE = 0.1    # relative, against the theory
# Of all electrons, over the second half of the run: slow electrons of the cone
# still leave then (0.35 % at any step, README.md).
LEAK_TOLERANCE = 0.01


def theory():
    centre = 0.5 * loss_cone.Cells[2] * loss_cone.CellSize
    peak = max(loss_cone.coil_field(z0) for z0 in loss_cone.CoilZ)
    starts = centre + np.linspace(-loss_cone.CloudHalfLength, loss_cone.CloudHalfLength, 201)
    return float(np.mean([1 - math.sqrt(1 - loss_cone.coil_field(z) / peak) for z in starts]))


def analyse(workdir):
    try:
        lines = (workdir / "energy.txt").read_text().splitlines()
        config = json.loads((workdir / "system_config.json").read_text())
    except (OSError, ValueError) as error:
        return f"cannot read run data: {error}"
    head = lines[0].split()
    table = np.array([[float(v) for v in line.split()] for line in lines[1:] if line.strip()])
    if table.size == 0 or not np.all(np.isfinite(table)):
        return "empty or non-finite energy data"
    time, count = table[:, head.index("Time")], table[:, head.index("ElectronsParticles")]
    half = np.argmin(np.abs(time - 0.5 * time[-1]))
    escaped = 1 - count[-1] / count[0]
    leak = (count[half] - count[-1]) / count[0]
    expected = theory()
    result = {"scheme": config["Scheme"], "dt": config["Dt"], "markers": int(count[0]),
              "escaped": escaped, "escaped_theory": expected, "escaped_error": abs(escaped / expected - 1),
              "second_half_leak": leak}
    result["passed"] = bool(result["escaped_error"] < FRACTION_TOLERANCE and leak < LEAK_TOLERANCE)
    return result


def main(name):
    result = analyse(Path(sys.argv[1]) if len(sys.argv) > 1 else Path("."))
    if isinstance(result, str):
        print(f"{name}: {result}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2))
    return 0 if result["passed"] else 1

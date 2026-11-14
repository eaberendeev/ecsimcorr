#!/usr/bin/env python3
"""Fast magnetosonic wave, fully kinetic, with the symplectic scheme.

The case and its parameters are in tests/integration/wave_cases/fms_kinetic.py.
"""

import sys
from pathlib import Path

# The integration runner executes a copy of this file from the project root.
_root = next(path for path in Path(__file__).resolve().parents if (path / "srcBeren").is_dir())
sys.path.insert(0, str(_root / "tests" / "integration" / "wave_cases"))
import fms_kinetic  # noqa: E402

# 4k markers x 62k explicit steps: about 3 minutes. Energy is not conserved by
# the explicit scheme: on this grid (dx = 226 Debye lengths) the electrons gain
# 4.5 % over the run; the check reports the drift and tests the frequency.
NumProcs = 2
DirName = "Res_FMS_symplectic"

system_config, particles_config = fms_kinetic.configuration("symplectic")


def generate_config():
    fms_kinetic.write_configuration("symplectic")
    return DirName


if __name__ == "__main__":
    generate_config()

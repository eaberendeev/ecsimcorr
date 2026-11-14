#!/usr/bin/env python3
"""Fast magnetosonic wave, fully kinetic, with the ecsim scheme.

The case and its parameters are in tests/integration/wave_cases/fms_kinetic.py.
"""

import sys
from pathlib import Path

# The integration runner executes a copy of this file from the project root.
_root = next(path for path in Path(__file__).resolve().parents if (path / "srcBeren").is_dir())
sys.path.insert(0, str(_root / "tests" / "integration" / "wave_cases"))
import fms_kinetic  # noqa: E402

# 4k markers x 15k steps: about 4 minutes; more threads do not help a 32-cell column
NumProcs = 2
DirName = "Res_FMS_ecsim"

system_config, particles_config = fms_kinetic.configuration("ecsim")


def generate_config():
    fms_kinetic.write_configuration("ecsim")
    return DirName


if __name__ == "__main__":
    generate_config()

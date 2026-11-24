#!/usr/bin/env python3
"""Loss cone of a magnetic mirror with the ecsim scheme.

The case and its parameters are in tests/integration/mirror_cases/loss_cone.py.
"""

import sys
from pathlib import Path

# The integration runner executes a copy of this file from the project root.
_root = next(path for path in Path(__file__).resolve().parents if (path / "srcBeren").is_dir())
sys.path.insert(0, str(_root / "tests" / "integration" / "mirror_cases"))
import loss_cone  # noqa: E402

# 16k markers x 750 steps
NumProcs = 4
DirName = "Res_LossCone_ecsim"

system_config, particles_config = loss_cone.configuration("ecsim")


def generate_config():
    loss_cone.write_configuration("ecsim")
    return DirName


if __name__ == "__main__":
    generate_config()

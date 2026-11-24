#!/usr/bin/env python3
"""Loss cone of a magnetic mirror with the ecsim_corr scheme.

The case and its parameters are in tests/integration/mirror_cases/loss_cone.py.
One corrector pass, as needed at Omega_e dt > 0.3 (wave_cases/README.md).
"""

import sys
from pathlib import Path

# The integration runner executes a copy of this file from the project root.
_root = next(path for path in Path(__file__).resolve().parents if (path / "srcBeren").is_dir())
sys.path.insert(0, str(_root / "tests" / "integration" / "mirror_cases"))
import loss_cone  # noqa: E402

# 16k markers x 750 steps
NumProcs = 4
DirName = "Res_LossCone_ecsim_corr"

system_config, particles_config = loss_cone.configuration("ecsim_corr", CorrectorPasses=1)


def generate_config():
    loss_cone.write_configuration("ecsim_corr", CorrectorPasses=1)
    return DirName


if __name__ == "__main__":
    generate_config()

#!/usr/bin/env python3
"""Fast magnetosonic wave, fully kinetic, with the ecsim_corr scheme.

The case and its parameters are in tests/integration/wave_cases/fms_kinetic.py.

The scheme runs with one corrector pass ("CorrectorPasses": 1). The single-pass
scheme fails this case: an electron-cyclotron mode at k dx ~ 1.5-2 grows, the
ions lose 18 % of their energy to the electrons and the frequency is 6.7 % off
(tests/integration/wave_cases/README.md).
"""

import sys
from pathlib import Path

# The integration runner executes a copy of this file from the project root.
_root = next(path for path in Path(__file__).resolve().parents if (path / "srcBeren").is_dir())
sys.path.insert(0, str(_root / "tests" / "integration" / "wave_cases"))
import fms_kinetic  # noqa: E402

# 4k markers x 15k steps: about 8 minutes; more threads do not help a 32-cell column
NumProcs = 2
DirName = "Res_FMS_ecsim_corr"

system_config, particles_config = fms_kinetic.configuration("ecsim_corr")


def generate_config():
    fms_kinetic.write_configuration("ecsim_corr")
    return DirName


if __name__ == "__main__":
    generate_config()

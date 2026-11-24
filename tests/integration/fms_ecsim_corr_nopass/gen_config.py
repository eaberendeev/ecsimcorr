#!/usr/bin/env python3
"""Fast magnetosonic wave with the single-pass ecsim_corr ("CorrectorPasses": 0).

The case and its parameters are in tests/integration/wave_cases/fms_kinetic.py.
The single-pass scheme needs Omega_e dt <= 0.35 here (a numerical
electron-cyclotron mode grows at 0.4, wave_cases/README.md); this case runs it
at Omega_e dt = 0.2, twice the steps of fms_ecsim_corr.
"""

import sys
from pathlib import Path

# The integration runner executes a copy of this file from the project root.
_root = next(path for path in Path(__file__).resolve().parents if (path / "srcBeren").is_dir())
sys.path.insert(0, str(_root / "tests" / "integration" / "wave_cases"))
import fms_kinetic  # noqa: E402

# 4k markers x 30k steps: about 16 minutes; more threads do not help a 32-cell column
NumProcs = 2
DirName = "Res_FMS_ecsim_corr_nopass"

system_config, particles_config = fms_kinetic.configuration("ecsim_corr", corrector_passes=0, gyration_step=0.2)


def generate_config():
    fms_kinetic.write_configuration("ecsim_corr", corrector_passes=0, gyration_step=0.2)
    return DirName


if __name__ == "__main__":
    generate_config()

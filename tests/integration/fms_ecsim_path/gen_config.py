#!/usr/bin/env python3
"""Fast magnetosonic wave, fully kinetic, with the ecsim_path scheme.

The case and its parameters are in tests/integration/wave_cases/fms_kinetic.py.
ecsim_path always makes one corrector pass. Its energy balance is not exact:
here (cell = 226 Debye lengths) Omega_e dt = 0.4, omega_pe dt = 2 heats the
electrons through a slow electrostatic mode (the field energy grows 100 times
over the run); this case runs at Omega_e dt = 0.2, omega_pe dt = 1, where the
energy drifts by 1e-3 over the run, twice the steps of fms_ecsim_corr.
"""

import sys
from pathlib import Path

# The integration runner executes a copy of this file from the project root.
_root = next(path for path in Path(__file__).resolve().parents if (path / "srcBeren").is_dir())
sys.path.insert(0, str(_root / "tests" / "integration" / "wave_cases"))
import fms_kinetic  # noqa: E402

# 4k markers x 30k steps: about 16 minutes; more threads do not help a 32-cell column
NumProcs = 2
DirName = "Res_FMS_ecsim_path"

system_config, particles_config = fms_kinetic.configuration("ecsim_path", gyration_step=0.2)


def generate_config():
    fms_kinetic.write_configuration("ecsim_path", gyration_step=0.2)
    return DirName


if __name__ == "__main__":
    generate_config()

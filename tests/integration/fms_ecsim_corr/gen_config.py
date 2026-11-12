#!/usr/bin/env python3
"""Fast magnetosonic wave, fully kinetic, with the ecsim_corr scheme.

The case and its parameters are in tests/integration/wave_cases/fms_kinetic.py.

Known deviation, 2026-11-12: this check FAILS. ECSIM-Corr gives the frequency
6.7 % off the kinetic root (ECSIM on the same input: 0.4 %). Over the run the
ions lose 21 % of their perpendicular energy to the electrons (+18 %) and to
magnetic noise (field energy 8e-4 -> 3.8e-3, mostly B_y at k dx ~ 1), at
exact total energy balance (3e-13). The transfer is there without the wave
too: in the uniform plasma B noise reaches 1.1e-3 in 4000 steps. Here
dx = 10 c/omega_pe = 226 Debye lengths and omega_pe dt = 2. The case is kept
out of `--test all` until the cause is understood; the tolerance is not
loosened.
"""

import sys
from pathlib import Path

# The integration runner executes a copy of this file from the project root.
_root = next(path for path in Path(__file__).resolve().parents if (path / "srcBeren").is_dir())
sys.path.insert(0, str(_root / "tests" / "integration" / "wave_cases"))
import fms_kinetic  # noqa: E402

# 4k markers x 15k steps: about 5 minutes; more threads do not help a 32-cell column
IntegrationHeavy = True
NumProcs = 2
DirName = "Res_FMS_ecsim_corr"

system_config, particles_config = fms_kinetic.configuration("ecsim_corr")


def generate_config():
    fms_kinetic.write_configuration("ecsim_corr")
    return DirName


if __name__ == "__main__":
    generate_config()

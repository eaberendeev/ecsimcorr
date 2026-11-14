#!/usr/bin/env python3
"""Fast magnetosonic wave, fully kinetic, with the ecsim_corr scheme.

The case and its parameters are in tests/integration/wave_cases/fms_kinetic.py.

Known deviation: the single-pass ECSIM-Corr FAILS this check. The frequency is
6.7 % off the kinetic root (ECSIM on the same input: 0.4 %), the ions lose
18 % of their energy and the electrons gain 15 % at an exact total balance.
The cause is an electron-cyclotron instability of the predictor-corrector pair
at Omega_e dt = 0.4; it is removed by the corrector passes of PR #47
("CorrectorPasses": 1, branch fix/ecsim-corr-cyclotron-stability). Until that
change is merged here the case stays out of `--test all`; the tolerance is not
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

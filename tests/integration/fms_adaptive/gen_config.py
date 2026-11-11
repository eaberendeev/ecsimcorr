#!/usr/bin/env python3
"""Fast magnetosonic wave carried by electrons of two models at once.

The case of fms_hybrid with `model: "adaptive"` electrons and the time step
reduced to omega_ce dt = 0.4, where the kinetic orbit of the midpoint rule
still deposits 96 % of the magnetic moment (1 / (1 + (omega_ce dt / 2)^2)) and
both models describe the same plasma. The markers are sorted once, at the
start: omega_ce dt > 0.404 is a guiding centre, so the electrons of the
compressed half of the wave (dB_z > 0.01 B0) are guiding centres and those of
the rarefied half are kinetic. The wave then runs through this plasma for two
periods. check_fms_adaptive.py requires the frequency of the kinetic root and
exact energy balance.

Sorting on every step is not used here on purpose. In this case the wave
would flip whole regions between the models twice per period (40 thousand
changes); every change of a region replaces the response of its electrons to
the fields at the cyclotron frequency, the field noise grows 20 times and
the frequency is off by 3.7 % (docs/drift-kinetic.md). The model change is
meant for markers that cross a border between a magnetized and an
unmagnetized region one by one.
"""

import importlib.util
import json
import math
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "fms_hybrid_config", Path(__file__).resolve().parent.parent / "fms_hybrid" / "gen_config.py")
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)

# 15k steps of 4k markers: about 10 minutes on 8 cores
IntegrationHeavy = True
NumProcs = 8
DirName = "Res_FMS_Adaptive"

GyrationStep = 0.4                       # omega_ce dt in the unperturbed field
ToGuidingCenter, ToKinetic = 1.01, 0.99  # thresholds of omega_ce dt, in units of GyrationStep

system_config = json.loads(json.dumps(_base.system_config))
system_config["Dt"] = GyrationStep / _base.B0
system_config["LastTimestep"] = round(_base.LastTimestep * _base.Dt / system_config["Dt"])
system_config["ModelSwitch"] = {
    "omega_dt_to_guiding_center": ToGuidingCenter * GyrationStep,
    "omega_dt_to_kinetic": ToKinetic * GyrationStep,
    "interval": system_config["LastTimestep"] + 1,   # the initial sorting only
}
# every 5th step is enough for the fit of a period of 7300 steps
for output in system_config["diagnostics"]["outputs"]:
    if output["type"] == "probes":
        output["interval"] = 5

particles_config = json.loads(json.dumps(_base.particles_config))
for entry in particles_config["particles"]:
    if entry["Name"] == "Electrons":
        entry["model"] = "adaptive"


def generate_config():
    for filename, config in (("system_config.json", system_config),
                             ("particles_config.json", particles_config)):
        with open(filename, "w", encoding="utf-8") as stream:
            json.dump(config, stream, indent=2)
    wp = math.sqrt(4 * math.pi * system_config["n0"] * 4.80320427e-10 ** 2 / 9.10938356e-28)
    with open("phys.par", "w", encoding="utf-8") as stream:
        stream.write(f"w_p = {wp}\n1/w_p = {1 / wp}\n")
    return DirName


if __name__ == "__main__":
    generate_config()

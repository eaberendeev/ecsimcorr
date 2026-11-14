"""The existing ECSIM Weibel experiment, with the explicit symplectic time step.

Geometry, species, temperatures, particle count and random initialization are
imported from the ECSIM fixture. The integration is intentionally heavy.
"""
import copy
import importlib.util
import json
from pathlib import Path

PROJECT_ROOT = next(path for path in Path(__file__).resolve().parents if (path / "srcBeren").is_dir())
REFERENCE_CONFIG = PROJECT_ROOT / "tests/integration/Weibel_ecsim/gen_config.py"
spec = importlib.util.spec_from_file_location("weibel_ecsim_reference", REFERENCE_CONFIG)
reference = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reference)

IntegrationHeavy = True
DirName = "Res_Weibel_symplectic"
NumProcs = 4
Dx = Dy = Dz = reference.Dx
NumPartPerCell = reference.NumPartPerCell
LastTime = reference.LastTime
Dt = 0.15
LastTimestep = round(LastTime / Dt)

system_config = copy.deepcopy(reference.system_config)
particles_config = copy.deepcopy(reference.particles_config)
system_config.update(Scheme="symplectic", Dt=Dt, LastTimestep=LastTimestep)
# The reference has only zero external B and an empty coil list. A uniform
# magnetic flux, when present, belongs to the evolved symplectic B state.
system_config.pop("ExternalFieldB", None)
for output in system_config["diagnostics"]["outputs"]:
    if output["type"] in ("fields_2d", "density_2d"):
        physical_interval = output["interval"] * reference.Dt
        output["interval"] = round(physical_interval / Dt)
# Overwritten checkpoints bound disk usage while retaining a restart roughly
# every 75 inverse plasma frequencies. Output interval 0 is never used.
system_config["diagnostics"]["outputs"].append({"type": "recovery", "interval": round(75 / Dt)})
system_config["TimeStepDelayDiag2D"] = round(reference.TimeStepDelayDiag2D * reference.Dt / Dt)


def generate_config():
    system_config["LastTimestep"] = LastTimestep
    for name, config in (("system_config.json", system_config), ("particles_config.json", particles_config)):
        Path(name).write_text(json.dumps(config, indent=2) + "\n")
    Path("phys.par").write_text(f"w_p = {reference.w_p}\n1/w_p = {1/reference.w_p}\n")
    return f"{DirName}_Dx_{Dx}_np_{NumPartPerCell}_Dt_{Dt}"

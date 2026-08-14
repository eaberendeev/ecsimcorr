#!/usr/bin/env python3
# Integration test: Weibel instability — ecsim_corr scheme.
# Anisotropic electron temperature (T⊥=100eV, T∥=1keV), periodic BCs.
# Compares field maps and growth rate with catcode paper Sec 3.4.

import json
import math
import os

system_config = {}

Scheme_name = "ecsim_corr"

BoundaryConditions = []
BoundaryConditions.append({"periodic": {"face": "XMIN"}})
BoundaryConditions.append({"periodic": {"face": "YMIN"}})
BoundaryConditions.append({"periodic": {"face": "ZMIN"}})

Collider = "None"
Collisions = []

StartFromTime = 0
StartTimeStep = 0
LastTime = 750  # end time in 1/w_p (~500 steps at Dt=1.5)

NumProcs = 4
NumAreas = 1

DirName = "Res_Weibel_ecsim_corr"
DEBUG = False

Dx = 0.3
Dy = Dx
Dz = Dx
Dt = 1.0
LastTimestep = int(round(LastTime / Dt + 1))

Tau = 0
NumCellsX = 30
NumCellsY = 30
NumCellsZ = 5



Coils = []

DampingType = "None"
DampCellsX_glob = [0, 0]
DampCellsY_glob = [0, 0]
DampCellsZ_glob = [0, 0]

n0 = 1.e13
k_particles_reservation = -1
NumPartPerCell = 100
BUniform = [0., 0., 0.]

BExternal = []
BExternal.append({"uniform_field": {"value": BUniform}})
BExternal.append({"coils": Coils})
system_config["ExternalFieldB"] = BExternal

w_p = math.sqrt(4.0 * math.pi * n0 * (4.80320427e-10) ** 2 / 9.10938356e-28)

bbox_centerX = 0.5 * Dx * NumCellsX
bbox_centerY = 0.5 * Dy * NumCellsY
bbox_centerZ = 0.5 * Dz * NumCellsZ

electrons = {}
electrons["Name"] = "Electrons"
electrons["Charge"] = -1.0
electrons["Density"] = 1.0
electrons["Mass"] = 1.0
electrons["NumPartPerCell"] = NumPartPerCell

electrons["distribution"] = [
    {
        "type": "injection",
        "dist_space": {
            "type": "rectangle",
            "center": [bbox_centerX, 
                       bbox_centerY, 
                       bbox_centerZ],
            "half_length": [bbox_centerX, 
                            bbox_centerY, 
                            bbox_centerZ],
        },
        "dist_pulse": {
            "type": "gaussian",
            "mean": [0, 0, 0],
            "sigma": [0.1, 0.1, 1.]
        },
        "density": 1.0,
    }
]

particles_config = {"particles": []}
particles_config["particles"].append(electrons)

RecoveryInterval = 20000
TimeStepDelayDiag1D = 1
TimeStepDelayDiag2D = 12

DiagDict = {
    "TimeStepDelayDiag2D": 2,
    "outputs": [
        {"type": "energy_balance", "interval": 1},
        {"type": "boundary_stats", "interval": 1},
        {"type": "console_summary"},
        {"type": "density_2d", "interval": TimeStepDelayDiag2D,
         "species": "all", "planes": {"z": [bbox_centerZ]}},
        {"type": "fields_2d", "interval": TimeStepDelayDiag2D,
         "fields": ["E", "B"], "planes": {"z": [bbox_centerZ]}},
    ],
}

system_config["Scheme"] = Scheme_name
system_config["Dx"] = Dx
system_config["Dy"] = Dy
system_config["Dz"] = Dz
system_config["Dt"] = Dt
system_config["NumCellsX"] = NumCellsX
system_config["NumCellsY"] = NumCellsY
system_config["NumCellsZ"] = NumCellsZ
system_config["DampCellsX_glob"] = DampCellsX_glob
system_config["DampCellsY_glob"] = DampCellsY_glob
system_config["DampCellsZ_glob"] = DampCellsZ_glob
system_config["Coils"] = Coils
system_config["diagnostics"] = DiagDict
system_config["StartTimeStep"] = StartTimeStep
system_config["LastTimestep"] = LastTimestep
system_config["RecoveryInterval"] = RecoveryInterval
system_config["TimeStepDelayDiag1D"] = TimeStepDelayDiag1D
system_config["TimeStepDelayDiag2D"] = TimeStepDelayDiag2D
system_config["Collider"] = Collider
system_config["Collisions"] = Collisions
system_config["DampingType"] = DampingType
system_config["n0"] = n0
system_config["k_particles_reservation"] = k_particles_reservation
system_config["NumPartPerCell"] = NumPartPerCell
system_config["StartFromTime"] = StartFromTime
system_config["Tau"] = Tau
system_config["Boundary_conditions"] = BoundaryConditions


def generate_config():
    with open("system_config.json", "w", encoding="utf-8") as f:
        json.dump(system_config, f, indent=2, ensure_ascii=False)
    with open("particles_config.json", "w", encoding="utf-8") as f:
        json.dump(particles_config, f, indent=2, ensure_ascii=False)
    f = open("phys.par", "w")
    f.write("w_p = " + str(w_p) + "\n")
    f.write("1/w_p = " + str(1.0 / w_p))
    f.close()
    return DirName + "_Dx_" + str(Dx) + "_np_" + str(NumPartPerCell) + "_Dt_" + str(Dt)

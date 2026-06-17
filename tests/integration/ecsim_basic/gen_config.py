#!/usr/bin/env python3
# Integration test: basic ecsim — small grid, few timesteps, no collisions.
# Verifies the simulation runs end-to-end and energy is conserved.

import json
import math

system_config = {}

Scheme_name = "ecsim"

BoundaryConditions = []
BoundaryConditions.append({"open": {"face": "CYLINDER"}})
BoundaryConditions.append({"periodic": {"face": "ZMIN"}})

Collider = "None"
Collisions = []

StartFromTime = 0
StartTimeStep = 0
LastTimestep = 20

NumProcs = 1
NumAreas = 1

DirName = "Res_TestInt_Basic"
DEBUG = False

Dx = 2
Dy = Dx
Dz = Dx
Dt = 2

Tau = 0
NumCellsX = 10
NumCellsY = 10
NumCellsZ = 2

Cyl = {
    "CylinderDomain": {
        "radius": 20.0,
        "center": [10.0, 10.0],
    },
}
system_config.update(Cyl)

Coils = []

DampingType = "None"
DampCellsX_glob = [0, 0]
DampCellsY_glob = [0, 0]
DampCellsZ_glob = [0, 0]

n0 = 1.e13
k_particles_reservation = -1
NumPartPerCell = 4
BUniform = [0., 0., 0.2]

BExternal = []
BExternal.append({"uniform_field": {"value": BUniform}})
BExternal.append({"coils": Coils})
system_config["ExternalFieldB"] = BExternal

ElectronsInjectDensity = 0.5
Tx = 0.01
Ty = 0.01
Tz = 0.01
ratio_x = 1.
ratio_y = 1.

w_p = math.sqrt(4.0 * math.pi * n0 * (4.80320427e-10) ** 2 / 9.10938356e-28)

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
            "type": "cylinder_ring_z",
            "center": [10.0, 10.0, 1.0],
            "r1": 1,
            "r2": 2,
            "half_length": 1.0,
        },
        "dist_pulse": {
            "type": "tangential",
            "rotation_center": [10.0, 10.0, 1.0],
            "mean_speed": 0.0,
            "sigma_speed": 0.0,
            "thermal_sigma": [Tx, Ty, Tz],
        },
        "density": ElectronsInjectDensity,
    }
]

IonsInjectDensity = 0.5

ions = {}
ions["Name"] = "Ions"
ions["Charge"] = 1.0
ions["Density"] = 1.0
ions["Mass"] = 1836.0
ions["NumPartPerCell"] = NumPartPerCell

ions["distribution"] = [
    {
        "type": "injection",
        "dist_space": {
            "type": "cylinder_ring_z",
            "center": [10.0, 10.0, 1.0],
            "r1": 1,
            "r2": 2,
            "half_length": 1.0,
        },
        "dist_pulse": {
            "type": "tangential",
            "rotation_center": [10.0, 10.0, 1.0],
            "mean_speed": 0.0,
            "sigma_speed": 0.0,
            "thermal_sigma": [0.001, 0.001, 0.001],
        },
        "density": IonsInjectDensity,
    }
]

particles_config = {"particles": []}
particles_config["particles"].append(electrons)
particles_config["particles"].append(ions)

RecoveryInterval = 200
TimeStepDelayDiag1D = 1
TimeStepDelayDiag2D = 2

DiagDict = {"outTime3D": [5, 150, 200],
            "zondCoordsLineX": [(0.0, 10.0, 2.0), (0.0, 10.0, -30), (0.0, 20.0, 2.0)],
            "zondCoordsLineY": [(10.0, 0, 2.0), (5.0, 0, 2.0), (20.0, 0, 2.0)],
            "zondCoordsLineZ": [(10.0, 10.0, 0.0), (5.0, 10.0, 0.0), (20.0, 10.0, 0.0)],
            "sliceFieldsPlaneX": [10.0],
            "sliceFieldsPlaneY": [10.0],
            "sliceFieldsPlaneZ": [2.0, 0, 2],
            "sliceRadiationPlaneY": [40, 760],
            "sliceRadiationPlaneZ": [40, -30],
            "radiationDiagRadiuses": [-35.0, -75.0, -15.0],
            "TimeStepDelayDiag2D": 2}

DiagDict["TimeStepDelayDiag2D"] = 2

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

#!/usr/bin/env python3
# Integration test: Weibel instability — ecsim_corr scheme.
# Anisotropic electron temperature (T⊥=100eV, T∥=1keV), periodic BCs.
# Compares field maps and growth rate with catcode paper Sec 3.4.

import json
import math

system_config = {}

Scheme_name = "ecsim_corr"

BoundaryConditions = []
BoundaryConditions.append({"open": {"face": "CYLINDER"}})
BoundaryConditions.append({"open": {"face": "ZMIN"}})
BoundaryConditions.append({"open": {"face": "ZMAX"}})

BoundaryConditions.append({
    "electron_reflection": {
        "face": "ZMIN",
        "radius": 5.0,
        "energy_threshold": 0.1 / 512.,
    }
})

Collider = "None"
Collisions = [
    # {"type": "coulomb", "species1": "Electrons", "species2": "Electrons", "coulomb_log": 15},
    # {"type": "coulomb", "species1": "Ions", "species2": "Ions", "coulomb_log": 15},
    {"type": "coulomb", "species1": "Electrons", "species2": "Ions", "coulomb_log": 15*500},
    # {"type": "neutral_ionization", "charged_species": "Electrons", "neutral_species": "Neutrals",
    #  "electron_product": "Electrons", "ion_product": "Ions", "scheme": "physical_only",
    #  "electron_ionization": True, "proton_ionization": True, "proton_charge_exchange": True},
]

StartFromTime = 0
StartTimeStep = 0
LastTimestep = 30000

NumProcs = 1
NumAreas = 1

DirName = "Res_Weibel_ecsim_corr"
DEBUG = False

Dx = 0.5
Dy = Dx
Dz = Dx
Dt = 1.5

Tau = 4998
NumCellsX = 140
NumCellsY = 140
NumCellsZ = 30



Coils = []

DampingType = "None"
DampCellsX_glob = [0, 0]
DampCellsY_glob = [0, 0]
DampCellsZ_glob = [0, 0]

n0 = 1.e14
k_particles_reservation = -1
NumPartPerCell = 1000
BUniform = [0., 0., 0.2]

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

electron_dist_space = {
    "type": "cylinder_z",
    "center": [
        0.5 * NumCellsX * Dx,
        0.5 * NumCellsY * Dy,
        0.5 * NumCellsZ * Dz,
    ],
    "radius": 10,
    "half_length": 7.5,
}
Te = 0.1 # kev
electron_dist_momentum = {"type": "gaussian", "mean": [0, 0, 0], "sigma": [Te, Te, Te]}

electrons["distribution"] = [
    {
        "type": "injection",
        "dist_space": electron_dist_space,
        "dist_pulse": electron_dist_momentum,
        "density": Dt / Tau,
    },
]

Ti = 0.1 # kev
ion_dist_space = electron_dist_space
ion_dist_momentum = {"type": "gaussian", "mean": [0, 0, 0], "sigma": [Ti, Ti, Ti]}

ions = {
    "Name": "Ions",
    "Charge": 1.0,
    "Density": 1.0,
    "Mass": 100.0,
    "NumPartPerCell": NumPartPerCell,
    "distribution": [
        {
            "type": "injection",
            "dist_space": ion_dist_space,
            "dist_pulse": ion_dist_momentum,
            "density": Dt / Tau,
        }
    ],
}

particles_config = {"particles": []}
particles_config["particles"].append(electrons)
particles_config["particles"].append(ions)

RecoveryInterval = 20000
TimeStepDelayDiag1D = 1
TimeStepDelayDiag2D = 48

DiagDict = {"sliceFieldsPlaneX": [bbox_centerX],
            "sliceFieldsPlaneY": [bbox_centerY],
            "sliceFieldsPlaneZ": [bbox_centerZ],
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

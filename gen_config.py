# Python script for setting parameters
import pprint
import sys
import json

import math

system_config = {}

Scheme_name = "ecsim_corr"

# face: XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX, CYLINDER
BoundaryConditions = []

BoundaryConditions.append({"open": {"face": "CYLINDER"}})
BoundaryConditions.append({"periodic": {"face": "ZMIN"}})
# BoundaryConditions.append({"open": {"face": "ZMAX"}})

# Tx = Ty = Tz = 0.005 # Kev
# BoundaryConditions.append(
#     {
#         "second_emisson": {
#             "face": "ZMIN",
#             "mean": [0, 0, 0],
#             "sigma": [Tx, Ty, Tz],
#         }
#     }
# )
# BoundaryConditions.append(
#     {
#         "second_emisson": {
#             "face": "ZMAX",
#             "mean": [0, 0, 0],
#             "sigma": [Tx, Ty, Tz],
#         }
#     }
# )

Collider = "None" #BinaryCollider"  # "BinaryCollider" # None
#####
StartFromTime = 0

NumProcs = 1  # number of processors
NumAreas = 1  # Number of decomposition region

DirName = "Res_Test"
DEBUG = False

Dx = 2  # step on X
Dy = Dx  # step on Y
Dz = Dx  # step on Z
Dt = 3  # 4*min(Dx,Dy)  # time step

Tau = 4998
NumCellsX = 400  # Number of all cells in computation domain on X
NumCellsY = 400  # NumbeY of all cells in computation domain on Y
# for home usage set 100
NumCellsZ = 5  # NumbeY of all cells in computation domain on Z

Cyl = {
    "radius": 0.5 * Dx * NumCellsX,
    "center": [0.5*NumCellsX * Dx, 0.5*NumCellsY * Dy],
}
system_config["CylinderDomain"] = Cyl

damp = 0
DampingType = "None"  # "CircleXY" #"CircleXY" # CircleXY Rectangle
DampCellsX_glob = [damp, damp]  # Number of Damping layer cells on Z
DampCellsY_glob = [damp, damp]  # Number of Damping layer cells on Y
DampCellsZ_glob = [0, 0]  # Number of Damping layer cells on Y


NumPartPerLine = 1  # Number of particles per line segment cell
NumPartPerCell = 100  # NumPartPerLine**3 # Number of particles per cell
k_particles_reservation = -1.0

MaxTime = 600000  # in 1/w_p
RecTime = 600  #

DiagDelay2D = 6  # in 1 / w_p
DiagDelay1D = 1  # in 1 / w_p
outTime3D = [5, 150, 200]
# DiagDelay3D = 10*DiagDelay2D # in 1 / w_p

DiagDelayEnergy1D = 1

###### PHYSIC CONST
PI = 3.141592653589793
me = 9.10938356e-28  # electron mass
ee = 4.80320427e-10  # electron charge
n0 = 1.0e13  # particles / cm^3
w_p = (4 * PI * n0 * ee * ee / me) ** 0.5
cc = 2.99792458e10  # speed on light cm/sec
MC2 = 511.0

# External radial electric field, 1/r, in w_c / w_p
# system_config["ExternalFieldE"] = [
#     { "uniform_field": { "value": [0.0, 0.0, 0.0] } },
#     { "uniformly_charged_cylinder": { "radius": 25.0, "value": -0.025 } }
#   ]

########
# BUniform = [0, 0, 0.05]  # in w_c / w_p

## KASP
R_coil = 70
I_coil = 3
ncolis = 25
listR = list(R_coil for i in range(ncolis))
listI = list(I_coil * ((-1) ** i) for i in range(-12, 13))
listZ = list(Dz * NumCellsZ * (i + 1) / 2 for i in range(-12, 13))

R_coil = 32
I_coil = 2
ncolis = 0
listR = list(R_coil for i in range(ncolis))
listI = list(I_coil for i in range(ncolis))
listZ = [60.0, 140]

Coils = []

for z, r, i in zip(listZ, listR, listI):
    Coils.append({"z": z, "R": r, "I": i})

system_config["ExternalFieldB"] = [
    {"uniform_field": {"value": [0.0, 0.0, 0.2]}},
    {"coils": Coils},
]

#######################################


LastTimestep = int(round(MaxTime / Dt + 1))
RecoveryInterval = int(round(RecTime / Dt))
StartTimeStep = int(round(StartFromTime / Dt))
TimeStepDelayDiag2D = int(round(max(DiagDelay2D, Dt) / Dt))
TimeStepDelayDiag1D = int(round(max(DiagDelay1D, Dt) / Dt))


#### bounding box without absorbing layer coords

bbox_centerY = 0.5 * Dy * NumCellsY
bbox_centerZ = 0.5 * Dz * NumCellsZ
bbox_centerX = 0.5 * Dx * NumCellsX

bbox_minX = Dx * DampCellsX_glob[0]
bbox_maxX = Dx * (NumCellsX + DampCellsX_glob[0])
bbox_minY = Dy * DampCellsY_glob[0]
bbox_maxY = Dy * (NumCellsY - DampCellsY_glob[1])
bbox_minZ = Dz * DampCellsZ_glob[0]
bbox_maxZ = Dz * (NumCellsZ - DampCellsZ_glob[1])
bbox_lenX = bbox_maxX - bbox_minX
bbox_lenY = bbox_maxY - bbox_minY
bbox_lenZ = bbox_maxZ - bbox_minZ

zondCoordsLineX = [
    (0.0, bbox_centerY, bbox_centerZ),
    (0.0, bbox_centerY, bbox_maxZ - 20 * Dz),
    (0.0, bbox_maxY - 20 * Dy, bbox_centerZ),
]

zondCoordsLineY = [
    (bbox_centerX, 0, bbox_centerZ),
    (bbox_minX + 20 * Dx, 0, bbox_centerZ),
    (bbox_maxX - 20 * Dx, 0, bbox_centerZ),
]

zondCoordsLineZ = [
    (bbox_centerX, bbox_centerY, 0.0),
    (bbox_minX + 20 * Dx, bbox_centerY, 0.0),
    (bbox_maxX - 20 * Dx, bbox_centerY, 0.0),
]

sliceFieldsPlaneX = [bbox_centerX]  # , bbox_minX + 20*Dx, bbox_maxX - 20*Dx]
sliceFieldsPlaneY = [bbox_centerY]  # , bbox_minY + 20*Dy, bbox_maxY - 20*Dy]
sliceFieldsPlaneZ = [bbox_centerZ, 0, Dz]  # , bbox_minZ + 20*Dz, bbox_maxZ - 20*Dz]


########################################
#### Coords for radiation diagnostic
sliceRadiationPlaneY = [bbox_minY + 20 * Dy, bbox_maxY - 20 * Dy]
sliceRadiationPlaneZ = [bbox_minZ + 20 * Dz, bbox_maxZ - 20 * Dz]


radius1 = min(bbox_maxY - 20 * Dy - bbox_centerY, bbox_maxZ - 20 * Dz - bbox_centerZ)
radius2 = min(bbox_maxY - 40 * Dy - bbox_centerY, bbox_maxZ - 40 * Dz - bbox_centerZ)
radius3 = min(bbox_maxY - 10 * Dy - bbox_centerY, bbox_maxZ - 10 * Dz - bbox_centerZ)
radiationDiagRadiuses = [radius1, radius2, radius3]
for radius in radiationDiagRadiuses:
    if radius > 0.5 * bbox_lenY:
        print("Wrong radius for radiation diagnostics!\n")
        exit(-1)
###########################################

# dist_space =  {"type": "rectangle", "center" : [], "half_length" : []}
# dist_space =  {"type": "cylinder_z", "center" : [], "radius" : , "half_length" : }
# dist_space =  {"type": "cylinder_x", "center" : [], "radius" : ,"half_length" : }
# dist_pulse = {"type": "gaussian", "mean" : [], "sigma" : [] }

################
# ELECTRONS
################

electron_dist_space = {
    "type": "cylinder_z",
    "center": [
        0.5 * NumCellsX * Dx,
        0.5 * NumCellsY * Dy,
        0.5 * NumCellsZ * Dz,
    ],
    "radius": 5,
    "half_length": 7.5,
}
electron_dist_space = {
    "type": "cylinder_ring_z",
    "center": [
        0.5 * NumCellsX * Dx,
        0.5 * NumCellsY * Dy,
        0.5 * NumCellsZ * Dz,
    ],
    "r1": 40,
    "r2": 60,
    "half_length": 0.5*Dz * NumCellsZ,
}

center =  [
        0.5 * NumCellsX * Dx,
        0.5 * NumCellsY * Dy,
        0.5 * NumCellsZ * Dz]

Te = 1.0 # kev
electron_dist_momentum = {"type": "gaussian", "mean": [0, 0, 0], "sigma": [Te, Te, Te]}
Te = 0.01 # kev
v_tan = -5.6e-3
electron_dist_momentum = {"type": "tangential", "rotation_center": center, "mean_speed": v_tan, "sigma_speed": 0, "thermal_sigma": [Te, Te, Te]}

electrons = {
    "Name": "Electrons",
    "Charge": -1.0,
    "Density": 1.0,
    "Mass": 1.0,
    "NumPartPerCell": NumPartPerCell,
    "distribution": [
        {
            "type": "injection",
            "dist_space": electron_dist_space,
            "dist_pulse": electron_dist_momentum,
            "density": Dt / Tau,
        },
    ],
}


################
# IONS
################

Ti = 1.0 # kev
ion_dist_space = electron_dist_space
ion_dist_momentum = {"type": "gaussian", "mean": [0, 0, 0], "sigma": [Ti, Ti, Ti]}
ion_dist_momentum = {"type": "tangential", "rotation_center": center, "mean_speed": v_tan, "sigma_speed": 0, "thermal_sigma": [0, 0, 0]}

ions = {
    "Name": "Ions",
    "Charge": 1.0,
    "Density": 1.0,
    "Mass": 1836.0,
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

################
# NEUTRALS
################

neutrals = {
    "Name": "Neutrals",
    "Charge": 0.0,
    "Density": 1.0,
    "Mass": 1837.0,
    "NumPartPerCell": NumPartPerCell,
}

neutrals_dist_space1 = {
    "type": "rectangle",
    "center": [
        -Dx,
        0.5 * NumCellsY * Dy + 25,
        0.5 * NumCellsZ * Dz,
    ],
    "half_length": [
        Dx,
        25,
        0.5 * NumCellsZ * Dz,
    ],
}

neutrals_dist_space2 = {
    "type": "rectangle",
    "center": [
        0.5 * NumCellsX * Dx + 25,
        NumCellsY * Dy + Dy,
        0.5 * NumCellsZ * Dz,
    ],
    "half_length": [
        25,
        Dy,
        0.5 * NumCellsZ * Dz,
    ],
}

Tx = 15  # Kev
vx = (2 * Tx / neutrals["Mass"] / 511.0) ** 0.5
neutrals_dist_momentum1 = {"type": "gaussian", "mean": [vx, 0, 0], "sigma": [0, 0, 0]}
neutrals_dist_momentum2 = {"type": "gaussian", "mean": [0, -vx, 0], "sigma": [0, 0, 0]}

neutrals["distribution"] = [
    {
        "type": "injection_bound",
        "dist_space": neutrals_dist_space1,
        "dist_pulse": neutrals_dist_momentum1,
        "density": 0.1,
    },
    {
        "type": "injection_bound",
        "dist_space": neutrals_dist_space2,
        "dist_pulse": neutrals_dist_momentum2,
        "density": 0.1,
    },
]

##############3

particles_config = {"particles": []}
particles_config["particles"].append(electrons)
# particles_config["particles"].append(ions)
# particles_config["particles"].append(neutrals)

NumOfPartSpecies = len(particles_config["particles"])


#####//////////////////////////////

WorkDir = DirName + "_Dx_" + str(Dx) + "_np_" + str(NumPartPerCell) + "_Dt_" + str(Dt)

if NumCellsX % NumAreas != 0:
    print("***********************************************")
    print("WARNING!!! Domain decomposition is not correct!!")
    print("***********************************************")

###////////////////////////////////
DiagDict = {}
DiagDict["outTime3D"] = outTime3D
DiagDict["zondCoordsLineX"] = zondCoordsLineX
DiagDict["zondCoordsLineY"] = zondCoordsLineY
DiagDict["zondCoordsLineZ"] = zondCoordsLineZ
DiagDict["sliceFieldsPlaneX"] = sliceFieldsPlaneX
DiagDict["sliceFieldsPlaneY"] = sliceFieldsPlaneY
DiagDict["sliceFieldsPlaneZ"] = sliceFieldsPlaneZ
DiagDict["sliceRadiationPlaneY"] = sliceRadiationPlaneY
DiagDict["sliceRadiationPlaneZ"] = sliceRadiationPlaneZ
DiagDict["radiationDiagRadiuses"] = radiationDiagRadiuses
DiagDict["TimeStepDelayDiag2D"] = TimeStepDelayDiag2D

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
# system_config["BUniform"] = BUniform
system_config["Coils"] = Coils
system_config["diagnostics"] = DiagDict
system_config["StartTimeStep"] = StartTimeStep
system_config["LastTimestep"] = LastTimestep
system_config["RecoveryInterval"] = RecoveryInterval
system_config["TimeStepDelayDiag1D"] = TimeStepDelayDiag1D
system_config["TimeStepDelayDiag2D"] = TimeStepDelayDiag2D
system_config["Collider"] = Collider
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
    print(system_config)
    with open("particles_config.json", "w", encoding="utf-8") as f:
        json.dump(particles_config, f, indent=2, ensure_ascii=False)
    print(particles_config)

    f = open("phys.par", "w")
    f.write("w_p = " + str(w_p) + "\n")
    f.write("1/w_p = " + str(1.0 / w_p))
    f.close()
    return WorkDir

if __name__ == "__main__":
    generate_config()

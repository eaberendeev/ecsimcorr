# Python script for setting parameters
import pprint
import sys
import json

import math

sys.path.insert(0, "./utils")
from berenUtils import *

system_config = {}

Scheme_name = "ecsim_corr"


# DampType = ("NONE","DAMP","PML")
# BoundType = ("NONE","PERIODIC","OPEN","OPEN_RADIUS","NEIGHBOUR","MIRROR")

BoundTypeX = ["OPEN_RADIUS", "OPEN_RADIUS"]
BoundTypeY = ["OPEN_RADIUS", "OPEN_RADIUS"]
# BoundTypeX = ["OPEN", "OPEN"]
# BoundTypeY = ["OPEN", "OPEN"]
BoundTypeZ = ["OPEN", "OPEN"]

Collider = "None" #BinaryCollider"  # "BinaryCollider" # None
#####
StartFromTime = 0

NumProcs = 8  # number of processors
NumAreas = 1  # Number of decomposition region

DirName = "Res_Jz_m0.01"
DEBUG = False

Dx = 0.5  # step on X
Dy = Dx  # step on Y
Dz = Dx  # step on Z
Dt = 1.5  # 4*min(Dx,Dy)  # time step

Tau = 4998
NumCellsX_glob = 120  # Number of all cells in computation domain on X
NumCellsY_glob = 120  # NumbeY of all cells in computation domain on Y
# for home usage set 100
NumCellsZ_glob = 100  # NumbeY of all cells in computation domain on Z

damp = 0
DampingType = "None"  # "CircleXY" #"CircleXY" # CircleXY Rectangle
DampCellsX_glob = [damp, damp]  # Number of Damping layer cells on Z
DampCellsY_glob = [damp, damp]  # Number of Damping layer cells on Y
DampCellsZ_glob = [0, 0]  # Number of Damping layer cells on Y


NumPartPerLine = 1  # Number of particles per line segment cell
NumPartPerCell = 1000  # NumPartPerLine**3 # Number of particles per cell
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
# system_config["ExternalFieldE"] = {"uniformly_charged_cylinder" : {"radius": 25, "value" :-0.001 * 25}}

########
BUniform = [0, 0, 0.0]  # in w_c / w_p

## KASP
R_coil = 70
I_coil = 3
ncolis = 25
listR = list(R_coil for i in range(ncolis))
listI = list(I_coil * ((-1) ** i) for i in range(-12, 13))
listZ = list(Dz * NumCellsZ_glob * (i + 1) / 2 for i in range(-12, 13))

R_coil = 32
I_coil = 2
ncolis = 2
listR = list(R_coil for i in range(ncolis))
listI = list(I_coil for i in range(ncolis))
listZ = [60.0, 140]
## Coil parameters.
## number of coils, x-coord coil_1, radius coil_1 (c/w_p), current in coil_1 (e*c/r_e), x-coord colil_2, ...
coils = []
for z, r, i in zip(listZ, listR, listI):
    coils.append(z)
    coils.append(r)
    coils.append(i)

BCoil = [ncolis] + coils
print(BCoil)

#######################################


LastTimestep = int(round(MaxTime / Dt + 1))
RecoveryInterval = int(round(RecTime / Dt))
StartTimeStep = int(round(StartFromTime / Dt))
TimeStepDelayDiag2D = int(round(max(DiagDelay2D, Dt) / Dt))
TimeStepDelayDiag1D = int(round(max(DiagDelay1D, Dt) / Dt))


#### bounding box without absorbing layer coords

bbox_centerY = 0.5 * Dy * NumCellsY_glob
bbox_centerZ = 0.5 * Dz * NumCellsZ_glob
bbox_centerX = 0.5 * Dx * NumCellsX_glob

bbox_minX = Dx * DampCellsX_glob[0]
bbox_maxX = Dx * (NumCellsX_glob + DampCellsX_glob[0])
bbox_minY = Dy * DampCellsY_glob[0]
bbox_maxY = Dy * (NumCellsY_glob - DampCellsY_glob[1])
bbox_minZ = Dz * DampCellsZ_glob[0]
bbox_maxZ = Dz * (NumCellsZ_glob - DampCellsZ_glob[1])
bbox_lenX = bbox_maxX - bbox_minX
bbox_lenY = bbox_maxY - bbox_minY
bbox_lenZ = bbox_maxZ - bbox_minZ

#### ZONDS ##############
numZonds = 15
### First number - number of zonds
zondX = list(bbox_lenX / (numZonds + 1) * (x + 1) + bbox_minX for x in range(numZonds))
zondY = [bbox_maxY - 20 * Dy] * numZonds
zondZ = [bbox_maxZ - 20 * Dz] * numZonds

zondCoords = list(zip(zondX, zondY, zondZ))
# zondCoords = [item for sublist in zondCoords for item in sublist]


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
sliceFieldsPlaneZ = [bbox_centerZ]  # , bbox_minZ + 20*Dz, bbox_maxZ - 20*Dz]


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
        0.5 * NumCellsX_glob * Dx,
        0.5 * NumCellsY_glob * Dy,
        0.5 * NumCellsZ_glob * Dz,
    ],
    "radius": 10,
    "half_length": 15,
}
Te = 1.0 # kev
electron_dist_momentum = {"type": "gaussian", "mean": [0, 0, 0], "sigma": [Te, Te, Te]}

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
            "density": Dt/Tau,
        }
    ],
}


################
# IONS
################

Ti = 1.0 # kev
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
            "density": Dt/Tau,
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
        0.5 * NumCellsY_glob * Dy + 25,
        0.5 * NumCellsZ_glob * Dz,
    ],
    "half_length": [
        Dx,
        25,
        0.5 * NumCellsZ_glob * Dz,
    ],
}

neutrals_dist_space2 = {
    "type": "rectangle",
    "center": [
        0.5 * NumCellsX_glob * Dx + 25,
        NumCellsY_glob * Dy + Dy,
        0.5 * NumCellsZ_glob * Dz,
    ],
    "half_length": [
        25,
        Dy,
        0.5 * NumCellsZ_glob * Dz,
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
particles_config["particles"].append(ions)
# particles_config["particles"].append(neutrals)
with open("particles_config.json", "w", encoding="utf-8") as f:
    json.dump(particles_config, f, indent=2, ensure_ascii=False)
print(particles_config)
NumOfPartSpecies = len(particles_config["particles"])


#####//////////////////////////////

WorkDir = DirName + "_Dx_" + str(Dx) + "_np_" + str(NumPartPerCell) + "_Dt_" + str(Dt)

if DEBUG:
    WorkDir = "Res_Debug3"

if NumCellsX_glob % NumAreas != 0:
    print("***********************************************")
    print("WARNING!!! Domain decomposition is not correct!!")
    print("***********************************************")

###////////////////////////////////
DiagParams = {}
DiagDict = {}
DiagDict["outTime3D"] = outTime3D
DiagDict["zondCoords"] = zondCoords
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
system_config["NumCellsX"] = NumCellsX_glob
system_config["NumCellsY"] = NumCellsY_glob
system_config["NumCellsZ"] = NumCellsZ_glob
system_config["DampingCellsX"] = DampCellsX_glob
system_config["DampingCellsY"] = DampCellsY_glob
system_config["DampingCellsZ"] = DampCellsZ_glob
system_config["BUniform"] = BUniform
system_config["BCoil"] = BCoil
system_config["BoundTypeX"] = BoundTypeX
system_config["BoundTypeY"] = BoundTypeY
system_config["BoundTypeZ"] = BoundTypeZ
system_config["diagnostics"] = DiagDict
with open("system_config.json", "w", encoding="utf-8") as f:
    json.dump(system_config, f, indent=2, ensure_ascii=False)

setParams(DiagParams, "Diagnostics", DiagDict)

SysParams = {}
SysDict = {}

SysDict["Scheme"] = Scheme_name
SysDict["NumProcs"] = NumProcs
SysDict["NumAreas"] = NumAreas

SysDict["Dt"] = Dt
SysDict["Tau"] = Tau

SysDict["DampCellsX_glob"] = DampCellsX_glob
SysDict["DampCellsY_glob"] = DampCellsY_glob
SysDict["DampCellsZ_glob"] = DampCellsZ_glob

SysDict["NumOfPartSpecies"] = NumOfPartSpecies
SysDict["NumPartPerLine "] = NumPartPerLine
SysDict["NumPartPerCell"] = NumPartPerCell

SysDict["LastTimestep"] = LastTimestep
SysDict["RecoveryInterval"] = RecoveryInterval
SysDict["StartTimeStep"] = StartTimeStep
SysDict["TimeStepDelayDiag1D"] = TimeStepDelayDiag1D
SysDict["TimeStepDelayDiag2D"] = TimeStepDelayDiag2D

SysDict["BUniform"] = BUniform
SysDict["BCoil"] = BCoil
SysDict["BoundTypeX"] = BoundTypeX
SysDict["BoundTypeY"] = BoundTypeY
SysDict["BoundTypeZ"] = BoundTypeZ

SysDict["MC2"] = MC2
SysDict["n0"] = n0
SysDict["StartFromTime"] = StartFromTime
SysDict["Collider"] = Collider
SysDict["DampingType"] = DampingType

SysDict["PI"] = PI

SysDict["k_particles_reservation"] = k_particles_reservation

# print(f"Parameters saved to {Particles}")
writeParams("Diagnostics", "Diagnostics.cfg", DiagParams)

setParams(SysParams, "SystemParams", SysDict)
writeParams("SystemParams", "SysParams.cfg", SysParams)


f = open("phys.par", "w")
f.write("w_p = " + str(w_p) + "\n")
f.write("1/w_p = " + str(1.0 / w_p))
f.close()

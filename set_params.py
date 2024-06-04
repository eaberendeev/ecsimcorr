# Python script for setting parameters
import pprint
import sys
import math
sys.path.insert(0, "./Scripts")
from setInitParams import *

DirName = "Res_CircleInjNew2"


DampType = enum("NONE","DAMP","PML")
BoundType = enum("NONE","PERIODIC","OPEN","NEIGHBOUR","MIRROR")

BoundTypeX = [BoundType.PERIODIC, BoundType.PERIODIC]
BoundTypeY = [BoundType.PERIODIC, BoundType.PERIODIC]
BoundTypeZ = [BoundType.PERIODIC, BoundType.PERIODIC]


Queue = "home" # type of queue
Cluster = "icc" # cluster name (spb, nsu, sscc)


#####
StartFromTime = 0

NumProcs = 1 # number of processors
NumAreas = 1 # Number of decomposition region


Dx = 0.5 # step on X
Dy = Dx # step on Y
Dz = Dx # step on Z
Dt = 1.5 #4*min(Dx,Dy)  # time step

PlasmaCellsX_glob = 300 # Number of cells for Plasma on Z

PlasmaCellsY_glob = 1 # Number of cells for Plasma on R 
PlasmaCellsZ_glob = 1 # Number of cells for Plasma on R 

NumCellsY_glob = PlasmaCellsX_glob # NumbeY of all cells in computation domain on R
NumCellsZ_glob = 5 # NumbeY of all cells in computation domain on R

damp = 20
DampCellsX_glob = [damp,damp] # Number of Damping layer cells on Z
DampCellsY_glob = [damp,damp] # Number of Damping layer cells on Y
DampCellsZ_glob = [0,0] # Number of Damping layer cells on Y

NumCellsX_glob = PlasmaCellsX_glob #+ DampCellsX_glob[0]+DampCellsX_glob[1] # Number of all cells in computation domain on Z


NumPartPerLine = 1 # Number of particles per line segment cell 
NumPartPerCell = 2000 #NumPartPerLine**3 # Number of particles per cell
k_particles_reservation = -1.

MaxTime = 130000 # in 1/w_p
RecTime = 600 #

Tau = 6000


DiagDelay2D = 6 # in 1 / w_p
DiagDelay1D = 1 # in 1 / w_p
outTime3D = [5,150,200]
#DiagDelay3D = 10*DiagDelay2D # in 1 / w_p

DiagDelayEnergy1D = 1 

###### PHYSIC CONST
PI = 3.141592653589793
me = 9.10938356e-28 # electron mass
ee = 4.80320427e-10 # electron charge
n0 = 1.e13 # particles / cm^3
w_p = (4*PI*n0*ee*ee/me)**0.5
cc = 2.99792458e10 # speed on light cm/sec 
MC2 = 512.
########
BUniform = [0, 0, 0.2] # in w_c / w_p
R_coil = 70
I_coil = 3
ncolis = 25
listR = list(R_coil for i in range(ncolis))
listI = list(I_coil*((-1)**i) for i in range(-12,13))
listZ = list(Dz*NumCellsZ_glob*(i+1)/2 for i in range(-12,13))
## Coil parameters. 
## number of coils, x-coord coil_1, radius coil_1 (c/w_p), current in coil_1 (e*c/r_e), x-coord colil_2, ...
coils = []
for z,r,i in zip(listZ, listR, listI):
    coils.append(z)
    coils.append(r)
    coils.append(i)

BCoil = [ncolis] + coils
print(BCoil)

#######################################



LastTimestep = int(round(MaxTime/Dt+1))
RecoveryInterval = int(round(RecTime/Dt))
StartTimeStep = int(round(StartFromTime/Dt))
TimeStepDelayDiag2D = int(round(max(DiagDelay2D,Dt)/Dt))
TimeStepDelayDiag1D = int(round(max(DiagDelay1D,Dt)/Dt))


#### bounding box without absorbing layer coords

bbox_centerY = 0.5*Dy*NumCellsY_glob
bbox_centerZ = 0.5*Dz*NumCellsZ_glob
bbox_centerX = 0.5*Dx*NumCellsX_glob

bbox_minX = Dx*DampCellsX_glob[0]
bbox_maxX = Dx*( PlasmaCellsX_glob + DampCellsX_glob[0] )
bbox_minY = Dy*DampCellsY_glob[0]
bbox_maxY = Dy*( NumCellsY_glob - DampCellsY_glob[1] )
bbox_minZ = Dz*DampCellsZ_glob[0]
bbox_maxZ = Dz*( NumCellsZ_glob - DampCellsZ_glob[1] )
bbox_lenX = bbox_maxX - bbox_minX
bbox_lenY = bbox_maxY - bbox_minY
bbox_lenZ = bbox_maxZ - bbox_minZ

#### ZONDS ##############
numZonds = 15
### First number - number of zonds
zondX = list( bbox_lenX / (numZonds + 1) * (x+1) + bbox_minX for x in range(numZonds)  )
zondY = [ bbox_maxY - 20*Dy ] * numZonds
zondZ = [ bbox_maxZ - 20*Dz ] * numZonds

zondCoords = list(zip(zondX,zondY,zondZ) )
#zondCoords = [item for sublist in zondCoords for item in sublist]



zondCoordsLineX = [(0., bbox_centerY, bbox_centerZ),
                   (0., bbox_centerY, bbox_maxZ - 20*Dz),
                   (0., bbox_maxY - 20*Dy, bbox_centerZ) ]

zondCoordsLineY = [(bbox_centerX, 0, bbox_centerZ),
                   (bbox_minX + 20*Dx, 0, bbox_centerZ),
                   (bbox_maxX - 20*Dx, 0, bbox_centerZ) ]

zondCoordsLineZ = [(bbox_centerX, bbox_centerY, 0.),
                   (bbox_minX + 20*Dx, bbox_centerY, 0.),
                   (bbox_maxX - 20*Dx, bbox_centerY, 0.) ]

sliceFieldsPlaneX = [bbox_centerX] #, bbox_minX + 20*Dx, bbox_maxX - 20*Dx]
sliceFieldsPlaneY = [bbox_centerY]#, bbox_minY + 20*Dy, bbox_maxY - 20*Dy]
sliceFieldsPlaneZ = [bbox_centerZ]#, bbox_minZ + 20*Dz, bbox_maxZ - 20*Dz]


########################################
#### Coords for radiation diagnostic
sliceRadiationPlaneY = [bbox_minY + 20*Dy, bbox_maxY - 20*Dy]
sliceRadiationPlaneZ = [bbox_minZ + 20*Dz, bbox_maxZ - 20*Dz]


radius1 = min(bbox_maxY - 20*Dy - bbox_centerY, bbox_maxZ - 20*Dz -bbox_centerZ)
radius2 = min(bbox_maxY - 40*Dy - bbox_centerY, bbox_maxZ - 40*Dz -bbox_centerZ)
radius3 = min(bbox_maxY - 10*Dy - bbox_centerY, bbox_maxZ - 10*Dz -bbox_centerZ)
radiationDiagRadiuses = [ radius1, radius2, radius3]
for radius in radiationDiagRadiuses:
    if radius > 0.5*bbox_lenY:
        print("Wrong radius for radiation diagnostics!\n")
        exit(-1)
###########################################


MaxSizeOfParts=2*NumPartPerCell*PlasmaCellsX_glob*PlasmaCellsY_glob*PlasmaCellsZ_glob/NumProcs+1


NumOfPartSpecies = 0

PartParams = {} # 
 # 
isVelDiag = 0


PName="Electrons"

Exist = True
PartDict = {}
PartDict["Charge"] = -1.0
PartDict["Density"] = 1.
PartDict["Velocity"] = 0.0
PartDict["Mass"] = 1.0
PartDict["Temperature"] = (1./512.)**0.5 #(0.05/512.)**0.5 
PartDict["Px_max"] = 1.e-1 # 
PartDict["Px_min"] = -1.e-1 #
PartDict["WidthY"] = NumCellsY_glob*Dy - 90*Dy
PartDict["WidthZ"] = NumCellsZ_glob*Dz
PartDict["Shift"] = 0.0
PartDict["SmoothMass"] = 0.0
PartDict["BoundResumption"] = 1
InitDist = "StrictUniformCircle"
InitDist = "None"
#InitDist = "UniformCosX_dn_k"

PartDict["DistParams"] = [str(InitDist)]


if Exist:
    NumOfPartSpecies+=1
    setParams(PartParams, PName, PartDict)

PName="Ions"

Exist = True
PartDict = {}
PartDict["Charge"] = 1.0
PartDict["Density"] = 1.
PartDict["Velocity"] = 0.0
PartDict["Mass"] = 100.0
PartDict["Temperature"] = (10./512.)**0.5  
PartDict["Px_max"] = 1.0 # 
PartDict["Px_min"] = -1.0 #
PartDict["WidthY"] = NumCellsY_glob*Dy - 90*Dy
PartDict["WidthZ"] = NumCellsZ_glob*Dz
PartDict["Shift"] = 0.0
PartDict["SmoothMass"] = 0
PartDict["SmoothMassMax"] = 30
PartDict["SmoothMassSize"] = 15
PartDict["BoundResumption"] = 1

InitDist = "StrictUniformCircle"
InitDist = "None"
#InitDist = "UniformCircle"
#InitDist = "Uniform"


PartDict["DistParams"] = [str(InitDist)]


if Exist :
    NumOfPartSpecies+=1
    setParams(PartParams, PName, PartDict)


#####//////////////////////////////

WorkDir = DirName+"_Dx_"+str(Dx)+"_np_"+str(NumPartPerCell )+"_Dt_"+str(Dt)

if PlasmaCellsX_glob % NumAreas != 0:
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


setParams(DiagParams, "Diagnostics", DiagDict)

DefineParams = []
SysParams = []

setConst(SysParams,'NumProcs',[NumProcs])
setConst(SysParams,'NumAreas',[NumAreas])

setConst(SysParams,'Dx',[Dx])
setConst(SysParams,'Dy',[Dy])
setConst(SysParams,'Dz',[Dz])
setConst(SysParams,'Dt',[Dt])
setConst(SysParams,'Tau',[Tau])

setConst(SysParams,'NumCellsX_glob',[NumCellsX_glob])
setConst(SysParams,'NumCellsY_glob',[NumCellsY_glob])
setConst(SysParams,'NumCellsZ_glob',[NumCellsZ_glob])
setConst(SysParams,'DampCellsX_glob',DampCellsX_glob)
setConst(SysParams,'DampCellsY_glob',DampCellsY_glob)
setConst(SysParams,'DampCellsZ_glob',DampCellsZ_glob)
setConst(SysParams,'PlasmaCellsY_glob',[PlasmaCellsY_glob])
setConst(SysParams,'PlasmaCellsX_glob',[PlasmaCellsX_glob])
setConst(SysParams,'PlasmaCellsZ_glob',[PlasmaCellsZ_glob])

setConst(SysParams,'NumOfPartSpecies',[NumOfPartSpecies])
setConst(SysParams,'NumPartPerLine ',[NumPartPerLine ])
setConst(SysParams,'NumPartPerCell',[NumPartPerCell])

setConst(SysParams,'LastTimestep',[LastTimestep])
setConst(SysParams,'RecoveryInterval',[RecoveryInterval])
setConst(SysParams,'StartTimeStep',[StartTimeStep])
setConst(SysParams,'TimeStepDelayDiag1D',[TimeStepDelayDiag1D])
setConst(SysParams,'TimeStepDelayDiag2D',[TimeStepDelayDiag2D])

setConst(SysParams,'BUniform',BUniform)
setConst(SysParams,'BCoil',BCoil)
setConst(SysParams,'BoundTypeX',BoundTypeX)
setConst(SysParams,'BoundTypeY',BoundTypeY)
setConst(SysParams,'BoundTypeZ',BoundTypeZ)

setConst(SysParams,'MC2',[MC2])
setConst(SysParams,'n0',[n0])
setConst(SysParams,'StartFromTime',[StartFromTime])

setConst(SysParams,'MaxSizeOfParts',[MaxSizeOfParts])

setConst(SysParams,'PI',[PI])

setConst(SysParams,'k_particles_reservation',[k_particles_reservation])


writeParams("Particles","PartParams.cfg",PartParams)
writeParams("Diagnostics","Diagnostics.cfg",DiagParams)

writeConst('SysParams.cfg', SysParams)


f = open('phys.par', 'w')
f.write("w_p = " + str(w_p))
f.write("1/w_p = " + str(1./w_p))
f.close()

f = open('workdir.tmp', 'w')
f.write(WorkDir)
f.close()
f = open('queue.tmp', 'w')
f.write(Queue)
f.close()
f = open('cluster.tmp', 'w')
f.write(Cluster)
f.close()
f = open('proc.tmp', 'w')
f.write(str(NumProcs))
f.close()

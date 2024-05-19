# Python script for setting parameters
import pprint
import sys
import math
sys.path.insert(0, "./Scripts")
from setInitParams import *

DirName = "Res_InjectCircle2"


DampType = enum("NONE","DAMP","PML")

DAMP_FIELDS = DampType.DAMP

BoundType = enum("NONE","PERIODIC","OPEN","NEIGHBOUR")

BoundTypeX_glob = [BoundType.OPEN, BoundType.OPEN]
BoundTypeY_glob = [BoundType.OPEN, BoundType.OPEN]
BoundTypeZ_glob = [BoundType.OPEN, BoundType.OPEN]


Queue = "home" # type of queue
Cluster = "icc" # cluster name (spb, nsu, sscc)

####
RECOVERY = 0

DEBUG = 1

SHAPE = 1  # 1 - PIC, 2 - Parabolic e.t.c

IONIZATION = 0 # Fields ionization

PARTICLE_MASS = False #True ### Particle has individual mass
PARTICLE_MPW = False ### Particle has individual macro particle weight

UPD_PARTICLES = 1

UPD_FIELDS = 1 # Update Fields (for Debug)


#####

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

MaxTime = 131 # in 1/w_p
RecTime = 600 #


DiagDelay2D = 60 # in 1 / w_p
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
MC2 = 511.
########
BUniform = [0, 0, 0.2] # in w_c / w_p
B_norm = me*cc*w_p/ee
BCoil_max = 0.4 * 10000 / B_norm # Gauss / B_norm
R_coil = 344*Dx / 2.39
I_coil = BCoil_max*R_coil/(8*PI)*((1+0.25*(344*Dx/R_coil)**2)**1.5)
## Coil parameters. 
## number of coils, x-coord coil_1, radius coil_1 (c/w_p), current in coil_1 (e*c/r_e), x-coord colil_2, ...
BCoil = [0, 50-172 , R_coil, I_coil ,50+172, R_coil, I_coil] 


#######################################



MaxTimeStep = int(round(MaxTime/Dt+1))
RecTimeStep = int(round(RecTime/Dt))
StartTimeStep = int(round(RECOVERY/Dt))
#TimeStepDelayDiag3D = int(round(DiagDelay3D/Dt))
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





PxMax = 800 // NumAreas
PpMax = 800

MaxSizeOfParts=2*NumPartPerCell*PlasmaCellsX_glob*PlasmaCellsY_glob*PlasmaCellsZ_glob/NumProcs+1


NumOfPartSpecies = 0

PartParams = {} # 
 # 
isVelDiag = 0

Vb = 0.3#79 #0.995

PName="Neutrals"

Exist = False
PartDict = {}
PartDict["Charge"] = 0.0
PartDict["Density"] = 0.5
PartDict["Velocity"] = 0.0
PartDict["Mass"] = 4.*1836.0
Pot_I = 0.02459 # Kev
PartDict["Pot_I"] = Pot_I
PartDict["Pot_k"] = 1.
PartDict["Temperature"] = 0;# (0.014/512.)**0.5
PartDict["Px_max"] = 1.0 # 
PartDict["Px_min"] = -1.0 #
PartDict["WidthY"] = PlasmaCellsY_glob*Dy
PartDict["WidthZ"] = PlasmaCellsZ_glob*Dz
PartDict["Shift"] = 0.0
PartDict["SmoothMass"] = 1
InitDist = "StrictUniformCircle"
dn0 = 0.0
k0=3.55
PartDict["DistParams"] = [str(InitDist), dn0, k0]



if Exist :
    NumOfPartSpecies+=1
    setParams(PartParams, PName, PartDict)


PName="Electrons"

Exist = True
PartDict = {}
PartDict["Charge"] = -1.0
PartDict["Density"] = 1.
PartDict["Velocity"] = 0.0
PartDict["Mass"] = 1.0
PartDict["Temperature"] = (1.0/512.)**0.5 
PartDict["Px_max"] = 1.e-1 # 
PartDict["Px_min"] = -1.e-1 #
PartDict["WidthY"] = NumCellsY_glob*Dy - 90*Dy
PartDict["WidthZ"] = NumCellsZ_glob*Dz
PartDict["Shift"] = 0.0
PartDict["SmoothMass"] = 0.0
PartDict["BoundResumption"] = 1
InitDist = "StrictUniformCircle"
InitDist = "Uniform"
#InitDist = "UniformCosX_dn_k"
dn0 = 0.1
k0 = 1./Vb
InitDist = "None"

PartDict["DistParams"] = [str(InitDist), dn0, k0]


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
PartDict["Temperature"] = (10.0/512.)**0.5  
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

dn0 = 0.0
k0=3.55
PartDict["DistParams"] = [str(InitDist), dn0, k0]
Pot_I = 0.05442 # Kev
PartDict["Pot_I"] = Pot_I
PartDict["Pot_k"] = 2.

#SmoothMassMax = 800.
#SmoothMassSize = 0.15*PlasmaCellsZ_glob * Dz

if Exist :
    NumOfPartSpecies+=1
    setParams(PartParams, PName, PartDict)


PName="Ions2"

Exist = False
PartDict = {}
PartDict["Charge"] = 2.0
PartDict["Density"] = 0.5
PartDict["Velocity"] = 0.0
PartDict["Mass"] = 4.*1836.0
PartDict["Temperature"] = 0.0
PartDict["SmoothMass"] = 0.0
PartDict["Px_max"] = 1.0 # 
PartDict["Px_min"] = -1.0 #
PartDict["WidthY"] = PlasmaCellsY_glob*Dy
PartDict["WidthZ"] = PlasmaCellsZ_glob*Dz
PartDict["Shift"] = 0.0
PartDict["BoundResumption"] = 1

InitDist = "None" #"StrictUniformCircle"
dn0 = 0.0
k0=3.55
PartDict["DistParams"] = [str(InitDist), dn0, k0]

SmoothMassMax = 40.
SmoothMassSize = 10. #0.15*PlasmaCellsZ_glob*Dz

if Exist :
    NumOfPartSpecies+=1
    setParams(PartParams, PName, PartDict)



###### GLOBAL BEAMS PARAMETERS #####################
InjType = 0 # 0 - random inject, 1 - periodic strit inject
InjectSmoothTime = 10.001


SourceType = enum("NONE","FOCUSED_GAUSS","UNIFORM")
#### FOCUSED BEAMS 

Lmax = 10. / (cc / w_p) # cm / c/w_p
Larea = 0.5*PlasmaCellsX_glob*Dx
Rmax = 0.1* 2.5 / (cc / w_p)  # cm / c/w_p
Rfocus = 0.05  / (cc / w_p) # cm / c/w_p
Xfocus = 0.5 * Dx * PlasmaCellsX_glob



PName="BeamLeft"

Exist = False
PartDict = {}
PartDict["Charge"] = -1.0
PartDict["Density"] = 0.01
PartDict["Velocity"] = 0.3#Vb
PartDict["Mass"] = 1.0
PartDict["Temperature"] = 0#(5.014/512.)**0.5
PartDict["Px_max"] = 1.0
PartDict["Px_min"] = 0.0
PartDict["WidthY"] = PlasmaCellsY_glob*Dy#Rfocus * 0.5
PartDict["WidthZ"] = PlasmaCellsZ_glob*Dz#Rfocus * 0.5
PartDict["Focus"] = Dx * PlasmaCellsX_glob
PartDict["SourceType"] = SourceType.UNIFORM #FOCUSED_GAUSS
PartDict["SourceAngle"] = 0#PI/10;
PartDict["Shift"] = 0.0
InitDist = "None"
PartDict["DistParams"] = [str(InitDist)]


if Exist :
    NumOfPartSpecies+=1
    setParams(PartParams, PName, PartDict)

PName="BeamRight"

Exist = False
PartDict = {}
PartDict["Charge"] = -1.0
PartDict["Density"] = 0.01
PartDict["Velocity"] = -Vb
PartDict["Mass"] = 1.0
PartDict["Temperature"] = 0.0
PartDict["Px_max"] = 0.0
PartDict["Px_min"] = -10.0
PartDict["WidthY"] = Rfocus
PartDict["WidthZ"] = Rfocus
PartDict["Focus"] = 0.5 * Dx * PlasmaCellsX_glob
PartDict["SourceType"] = SourceType.FOCUSED_GAUSS
PartDict["Shift"] = 0.0
InitDist = "None"
PartDict["DistParams"] = [str(InitDist)]

if Exist :
    NumOfPartSpecies+=1
    setParams(PartParams, PName, PartDict)
    

### GLOBAL LASERS PARAMETERS ##############################
Las_tau = 3.48
Las_w0 = 25.44
Las_vg = (1. - 1. / (Las_w0**2) )**0.5

LasParams = {} # 

DelayLeft =0.
DelayRight = 0.
LasName="LaserLeft"
Exist = False
LasDict = {} #
LasDict["delay"] = 0#DelayLeft*Las_tau
LasDict["tau"] = Las_tau
LasDict["type"] = "Virtual"
LasDict["vg"] = Las_vg
LasDict["a0"] = 0.67
LasDict["sigma0"] = 1.86
LasDict["w0"] = Las_w0
LasDict["focus"] =  [0.5 * Dx * NumCellsX_glob,0.5 * Dy * NumCellsY_glob,0.5 * Dz * NumCellsZ_glob]
LasDict["y0"] =  0.5 * Dy * NumCellsY_glob
LasDict["start"] = [Dx*DampCellsX_glob[0] ,0.5 * Dy * NumCellsY_glob,0.5 * Dz * NumCellsZ_glob]

if Exist:
    setParams(LasParams, LasName, LasDict)

LaserAngle = PI/36.
LasName="LaserRight"
Exist = False
LasDict = {} #
LasDict["delay"] = 0#DelayRight*Las_tau
LasDict["tau"] = Las_tau
LasDict["vg"] = -Las_vg
LasDict["a0"] = 0.8
LasDict["sigma0"] =  5.28
LasDict["type"] = "Virtual"
LasDict["w0"] = Las_w0
LasDict["y0"] =  0.5 * Dy * NumCellsY_glob
LasDict["focus"] =  [0.5 * Dx * NumCellsX_glob,0.5 * Dy * NumCellsY_glob,0.5 * Dz * NumCellsZ_glob]
#LasDict["start_x"] = Dx*PlasmaCellsX_glob + Dx*DampCellsX_glob[0]-Dx
z_start = 0.5 * Dz * NumCellsZ_glob + math.tan(LaserAngle)*(Dx*(PlasmaCellsX_glob + DampCellsX_glob[0]) - 0.5 * Dx * NumCellsX_glob)
LasDict["start"] = [ Dx*(PlasmaCellsX_glob + DampCellsX_glob[0]),0.5 * Dy * NumCellsY_glob,z_start ]

if Exist:
    setParams(LasParams, LasName, LasDict)



#####//////////////////////////////

WorkDir = DirName+"_nx_"+str(PlasmaCellsX_glob)+"_np_"+str(NumPartPerCell )+"_Dt_"+str(Dt)


if SHAPE < 3:
    SHAPE_SIZE = 2
    CELLS_SHIFT = 1
else:
    SHAPE_SIZE = 3
    CELLS_SHIFT = 2    

ADD_NODES = 2 * CELLS_SHIFT + 1


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

setConst(SysParams,'const int','NumProcs',[NumProcs],None)
setConst(SysParams,'const int','NumAreas',[NumAreas],None)

setConst(SysParams,'const double','Dx',[Dx],None)
setConst(SysParams,'const double','Dy',[Dy],None)
setConst(SysParams,'const double','Dz',[Dz],None)
setConst(SysParams,'const double','Dt',[Dt],None)

setConst(SysParams,'const int','NumCellsX_glob',[NumCellsX_glob],None)
setConst(SysParams,'const int','NumCellsY_glob',[NumCellsY_glob],None)
setConst(SysParams,'const int','NumCellsZ_glob',[NumCellsZ_glob],None)
setConst(SysParams,'const int','NumCellsXmax_glob',[NumCellsX_glob + 2 * SHAPE_SIZE - 1],None)
setConst(SysParams,'const int','NumCellsYmax_glob',[NumCellsY_glob + 2 * SHAPE_SIZE - 1],None)
setConst(SysParams,'const int','DampCellsX_glob',DampCellsX_glob,None)
setConst(SysParams,'const int','DampCellsY_glob',DampCellsY_glob,None)
setConst(SysParams,'const int','DampCellsZ_glob',DampCellsZ_glob,None)
setConst(SysParams,'const int','PlasmaCellsY_glob',[PlasmaCellsY_glob],None)
setConst(SysParams,'const int','PlasmaCellsX_glob',[PlasmaCellsX_glob],None)
setConst(SysParams,'const int','PlasmaCellsZ_glob',[PlasmaCellsZ_glob],None)

setConst(SysParams,'const int','NumOfPartSpecies',[NumOfPartSpecies],None)
setConst(SysParams,'const int','NumPartPerLine ',[NumPartPerLine ],None)
setConst(SysParams,'const int','NumPartPerCell',[NumPartPerCell],None)

setConst(SysParams,'const int','MaxTimeStep',[MaxTimeStep],None)
setConst(SysParams,'const int','RecTimeStep',[RecTimeStep],None)
setConst(SysParams,'const int','StartTimeStep',[StartTimeStep],None)
setConst(SysParams,'const int','TimeStepDelayDiag1D',[TimeStepDelayDiag1D],None)
setConst(SysParams,'const int','TimeStepDelayDiag2D',[TimeStepDelayDiag2D],None)

setConst(SysParams,'const double','InjectSmoothTime',[InjectSmoothTime],None)

setConst(SysParams,'const double','BUniform',BUniform,None)
setConst(SysParams,'const double','BCoil',BCoil,None)
setConst(SysParams,'const int','BoundTypeY_glob',BoundTypeY_glob,None)
setConst(SysParams,'const int','BoundTypeZ_glob',BoundTypeZ_glob,None)
setConst(SysParams,'const int','BoundTypeX_glob',BoundTypeX_glob,None)

setConst(SysParams,'const double','MC2',[MC2],None)
setConst(SysParams,'const double','n0',[n0],None)


setConst(SysParams,'const int','isVelDiag',[isVelDiag],None)
setConst(SysParams,'const int','PxMax',[PxMax],None)
setConst(SysParams,'const int','PpMax',[PpMax],None)

setConst(SysParams,'const int','MaxSizeOfParts',[MaxSizeOfParts],None)

setConst(SysParams,'const double','PI',[PI],None)

setConst(SysParams,'const double','Rmax',[Rmax],None)
setConst(SysParams,'const double','Rfocus',[Rfocus],None)
setConst(SysParams,'const double','Lmax',[Lmax],None)
setConst(SysParams,'const double','Larea',[Larea],None)
setConst(SysParams,'const double','SmoothMassMax',[SmoothMassMax],None)
setConst(SysParams,'const double','SmoothMassSize',[SmoothMassSize],None)
setConst(SysParams,'const double','k_particles_reservation',[k_particles_reservation],None)

if DEBUG == 0:
    setConst(DefineParams,'#define','NDEBUG',[' '],None)
setConst(DefineParams,'#define','DEBUG',[DEBUG],None)
if PARTICLE_MASS:
    setConst(DefineParams,'#define','PARTICLE_MASS',[1],None)
if PARTICLE_MPW:
    setConst(DefineParams,'#define','PARTICLE_MPW',[1],None)
setConst(DefineParams,'#define','RECOVERY',[RECOVERY],None)
setConst(DefineParams,'#define','SHAPE',[SHAPE],None)
setConst(DefineParams,'#define','SHAPE_SIZE',[SHAPE_SIZE],None)
setConst(DefineParams,'#define','CELLS_SHIFT',[CELLS_SHIFT],None)
setConst(DefineParams,'#define','ADD_NODES',[ADD_NODES],None)
setConst(DefineParams,'#define','IONIZATION',[IONIZATION],None)

setConst(DefineParams,'#define','DEBUG',[DEBUG],None)
setConst(DefineParams,'#define','UPD_FIELDS',[UPD_FIELDS],None)
setConst(DefineParams,'#define','UPD_PARTICLES',[UPD_PARTICLES],None)
setConst(DefineParams,'#define','DAMP_FIELDS',[DAMP_FIELDS],None)

setConst(DefineParams,'#define','PML',[DampType.PML],None)
setConst(DefineParams,'#define','DAMP',[DampType.DAMP],None)
setConst(DefineParams,'#define','PERIODIC',[BoundType.PERIODIC],None)
setConst(DefineParams,'#define','OPEN',[BoundType.OPEN],None)
setConst(DefineParams,'#define','NEIGHBOUR',[BoundType.NEIGHBOUR],None)
setConst(DefineParams,'#define','SOURCE_FOCUSED_GAUSS',[SourceType.FOCUSED_GAUSS],None)
setConst(DefineParams,'#define','SOURCE_UNIFORM',[SourceType.UNIFORM],None)
setConst(DefineParams,'#define','SOURCE_NONE',[SourceType.NONE],None)



writeParams("Particles","PartParams.cfg",PartParams)
writeParams("Lasers","LasParams.cfg",LasParams)
writeParams("Diagnostics","Diagnostics.cfg",DiagParams)

writeConst('const.h', 'SysParams.cfg', SysParams)
writeDefine('defines.h',DefineParams)



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

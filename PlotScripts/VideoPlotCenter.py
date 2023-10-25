#!/home/eberendeev/anaconda3/bin/python3
#!/usr/bin/python3
#!/home/fano.icmmg/berendeev_e_a/anaconda3/bin/python3
#!/mnt/storage/home/eaberendeev/anaconda3/bin/python3
#!/mnt/storage/home/aaefimova/anaconda3/bin/python3

#from mpi4py import MPI

import time
import multiprocessing
import numpy as np
import os
import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as col
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
from matplotlib import rc
import pprint
import collections
import re


#чтение параметров расчёта
from ReadDataLib import *
from LibPlot import *

def vx_vy_to_vr_va(vx, vy, COS, SIN):
    if (vx.shape != vy.shape) : raise RuntimeError("vx, vy must have the same shape!")

    vr = + COS * vx + SIN * vy
    va = - SIN * vx + COS * vy

    return vr, va


def init_COS_SIN(bx,ex,by,ey):
    x0 = (ex-bx)/2 - 0.5

    X = np.zeros((ey-by, ex-bx), dtype=np.float32) +\
        np.fromiter((x-x0 for x in range(bx, ex, 1)), dtype=float)

    Y = X.transpose()  # no upside down turn here
    R = np.sqrt(X * X + Y * Y)

    return X/R, Y/R

def init_rmap(bx,ex,by,ey):
    r0 = int((ex+bx)/2)
    x0 = (ex-bx)/2
    y0 = (ey-by)/2
    rmap = [[] for r in range(r0)]

    # TODO: traverse only 1/4 of cells
    for x in range(bx, ex, 1):
        for y in range(by, ey, 1):
            r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0)
            if r2 >= r0 * r0:
                continue

            rmap[int(np.sqrt(r2))].append((y, x))

    return [np.asarray(r).transpose() for r in rmap]


def phi_averaged(data, rmap):
    result = np.zeros(len(rmap))

    for i, r in enumerate(rmap):
        result[i] = data[r[0], r[1]].mean()

    return result

#корневой каталог для расчёта (где файл параметров)
WorkDir='../'

Bz0 = 0.2
#словарь с системными параметрами
SystemParameters=ReadParameters()

SystemParameters['2D_Diag_Fields']=['1','1','1','1','1','1'] #ReadParamFile(WorkDir)

FileNameSignNum=3 #len(str(int(int(SystemParameters['Max_Time'][0])/int(SystemParameters['Diagn_Time'][0]))))
print(FileNameSignNum)

def SubPlot(x,y,fig):
	return fig.add_subplot(gs[y,x])


tdelay = int(SystemParameters['TimeDelay2D']) * float( SystemParameters['Dt'])


print(SystemParameters)

stepZ = "0051"
i = 0
imax = 270
iStep = 1
tC = []
Bz = []
PrrE = []
PppE = []
densE = []
PrrI = []
PppI = []
densI = []

while i <= imax:
	TimeStep=str(i).zfill(4)
	#имя файла результирующего
	#проверка на наличие уже построенного графика
	if(True):
		print('TimeStep=',TimeStep)
		MaxTime=int(TimeStep)*tdelay
		MaxTimeDt=round(MaxTime,1)
		#считываем плотности
#		DensDataXY= collections.OrderedDict()
#		for sort in SystemParameters['Particles'].keys():
#			DensDataXY.update(ReadDensFile2D(WorkDir,"Dens","PlaneZ",SystemParameters,sort,TimeStep))

		#PlotIntegMPI(WorkDir)
		FieldsTitles = ["Ex","Ey","Ez","Bx","By","Bz"]
#		FieldDataXY=ReadFieldsFile2D(WorkDir,"PlaneZ_"+stepZ+"_","xy",SystemParameters,TimeStep)
		Name = WorkDir + "Fields/Diag2D/FieldPlaneZ_"+stepZ+"_"
		FieldDataXY = ReadFieldsFile2DNew(Name,FieldsTitles,TimeStep)

		#Name = WorkDir + "Fields/Diag2D/FieldAvgPlaneZ_"
		#FieldDataAvgXY = ReadFieldsFile2DNew(Name,FieldsTitles,TimeStep)

		Jxy = collections.OrderedDict()
		DensXy = collections.OrderedDict()
		Ppp = collections.OrderedDict()
		Prr = collections.OrderedDict()
		for sort in SystemParameters['Particles'].keys():
			Name = WorkDir + "Particles/" + sort + "/Diag2D/CurrentPlaneZ_"+stepZ+"_"
			Jxy[sort] = ReadFieldsFile2DNew(Name,["Jx","Jy","Jz"],TimeStep)
			Name = WorkDir + "Particles/" + sort + "/Diag2D/DensPlaneZ_"+stepZ+"_"
			DensXy[sort] = ReadFieldsFile2DNew(Name,["Dens"],TimeStep)
			Name = WorkDir + "Particles/" + sort + "/Diag2D/PyyPlaneZ_"+stepZ+"_"
			Ppp[sort] = ReadFieldsFile2DNew(Name,["Ppp"],TimeStep)
			Name = WorkDir + "Particles/" + sort + "/Diag2D/PxxPlaneZ_"+stepZ+"_"
			Prr[sort] = ReadFieldsFile2DNew(Name,["Prr"],TimeStep)
		
		ex, ey = Prr["Electrons"].shape[0], Prr["Electrons"].shape[1]
		PrrE.append(Prr["Electrons"][ex//2][ey//2])
		PrrI.append(Prr["Ions"][ex//2][ey//2])
		ex, ey = Ppp["Electrons"].shape[0], Ppp["Electrons"].shape[1]
		PppE.append(Ppp["Electrons"][ex//2][ey//2])
		PppI.append(Ppp["Ions"][ex//2][ey//2])
		ex, ey = DensXy["Electrons"].shape[0], DensXy["Electrons"].shape[1]
		densE.append(DensXy["Electrons"][ex//2][ey//2])
		densI.append(DensXy["Ions"][ex//2][ey//2])

		ex, ey = FieldDataXY["Bz"].shape[0], FieldDataXY["Bz"].shape[1]
		tC.append(MaxTimeDt)
		Bz.append(FieldDataXY["Bz"][ex//2][ey//2])

	i=i+iStep

f = open("Center.txt", 'w')
f.write("time Bz densE densI PppE PppI PrrE PrrI")
f.write("\n")
for t in range(len(tC)):
    f.write(str(tC[t]) + " " + str(Bz[t]) + " " + str(densE[t]) + " " + str(densI[t]) + " " + str(PppE[t]) + " " + str(PppI[t]) + " " + str(PrrE[t]) + " " + str(PrrI[t]) ) 
    f.write("\n")
f.close()

#!/home/eberendeev/anaconda3/bin/python3
#!/usr/bin/python3
#!/mnt/storage/home/eaberendeev/anaconda3/bin/python3

import time
import multiprocessing
import numpy as np
import os
import matplotlib
# matplotlib.use('Qt4Agg')
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
# чтение параметров расчёта
from ReadDataLib import *
from LibPlot import *

FieldsAmp = 0.03 #амплитуда палитры для 3D полей
densAmp = 1.
Bz0 = 0.0

axis = "xy"
sliceNum = "200"
t0 = 851
tMax = 9999 
iStep = 50

# корневой каталог для расчёта (где файл параметров)
WorkDir = '../'
FileNameSignNum = 3 

# словарь с системными параметрами
SystemParameters = ReadParameters()

# SystemParameters['2D_Diag_Fields']=['1','1','1','1','1','1'] #ReadParamFile(WorkDir)
NumOfPartSpecies = int(SystemParameters['NumOfPartSpecies'])

def SubPlot(x,y,fig):
    return fig.add_subplot(gs[y,x])


def vx_vy_to_vr_va(vx, vy, COS, SIN):
    if vx.shape != vy.shape:
        raise RuntimeError("vx, vy must have the same shape!")

    vr = +COS * vx + SIN * vy
    va = -SIN * vx + COS * vy

    return vr, va


def init_COS_SIN(bx, ex, by, ey):
    x0 = (ex - bx) / 2
    COSP = np.zeros((ex - bx, ey - by))
    SINP = np.zeros((ex - bx, ey - by))
    for x in range(0, (ex - bx)):
        for y in range(0, (ey - by)):
            r = ((x - x0) * (x - x0) + (y - x0) * (y - x0)) ** 0.5
            if r != 0:
                cosp = (x - x0) / r
                sinp = (y - x0) / r
            else:
                cosp = 0
                sinp = 0
            COSP[x, y] = cosp
            SINP[x, y] = sinp

    return COSP, SINP


def init_rmap(bx, ex, by, ey):
    r0 = int((ex + bx) / 2)
    x0 = (ex - bx) / 2
    y0 = (ey - by) / 2
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

tdelay = int(SystemParameters['TimeDelay2D']) * float( SystemParameters['Dt'])

try:
      os.mkdir('../Anime/2D')
except OSError:
      print("подкаталог для 2D диагностики существуют")

print(SystemParameters)

if axis == "xz":
	plane = "Y"
elif axis == "yz":
	plane = "X"
elif axis == "xy":
	plane = "Z"
else:
	print("Invalid plane\n")
	exit(0)

grid_x = 0
if "Neutrals" in SystemParameters['Particles'].keys():
	grids_x = -1


for i in range(t0, tMax, iStep):
	TimeStep=str(i).zfill(4)
	#имя файла результирующего
	GraphName='Res'+axis.upper()+"_"+TimeStep+'.png'
	print('TimeStep=',TimeStep)

	fig = plt.figure(figsize=(NumOfPartSpecies*4+34, 42))

	#на основе того, сколько было сортов частиц определить сетку для графиков
 	#первый аргумент - вертикальное количество, 
 	#второй - горизонтальное (по горизонтале число сортов частиц по 2 столбца + 1 столбца для полей
	gs =gridspec.GridSpec(5,1+2*NumOfPartSpecies + grids_x,height_ratios=[0.1,1,1,1,1])

	FieldsTitles = ["Ex","Ey","Ez"]
	Name = WorkDir + "Fields/Diag2D/FieldE_Plane"+plane+"_"+sliceNum+"_"
	FieldE = ReadFieldsFile2DNew(Name,FieldsTitles,TimeStep)
	FieldsTitles = ["Bx","By","Bz"]
	Name = WorkDir + "Fields/Diag2D/FieldB_Plane"+plane+"_"+sliceNum+"_"
	FieldB = ReadFieldsFile2DNew(Name,FieldsTitles,TimeStep)

	J = collections.OrderedDict()
	Dens = collections.OrderedDict()
	Pxx = collections.OrderedDict()
	Pyy = collections.OrderedDict()
	Pzz = collections.OrderedDict()
	for sort in SystemParameters['Particles'].keys():
		Name = WorkDir + "Particles/" + sort + "/Diag2D/Density_Plane"+plane+"_"+sliceNum+"_"
		Dens[sort] = ReadFieldsFile2DNew(Name,["Dens"],TimeStep)
		if sort != "Neutrals":
			Name = WorkDir + "Particles/" + sort + "/Diag2D/Prr_Plane"+plane+"_"+sliceNum+"_"
			Pxx[sort] = ReadFieldsFile2DNew(Name,["Pxx"],TimeStep)
			Name = WorkDir + "Particles/" + sort + "/Diag2D/Ppp_Plane"+plane+"_"+sliceNum+"_"
			Pyy[sort] = ReadFieldsFile2DNew(Name,["Pyy"],TimeStep)
			Name = WorkDir + "Particles/" + sort + "/Diag2D/Pzz_Plane"+plane+"_"+sliceNum+"_"
			Pzz[sort] = ReadFieldsFile2DNew(Name,["Pzz"],TimeStep)			
			Name = WorkDir + "Particles/" + sort + "/Diag2D/Current_Plane"+plane+"_"+sliceNum+"_"
			J[sort] = ReadFieldsFile2DNew(Name,["Jx","Jy","Jz"],TimeStep)

	Plot2Ddens(FieldE["Ex"],axis,-FieldsAmp,FieldsAmp,"Ex", "field", 1, SystemParameters,SubPlot(0,1,fig),fig)
	Plot2Ddens(FieldE["Ey"],axis,-FieldsAmp,FieldsAmp,"Ey", "field", 1, SystemParameters,SubPlot(0,2,fig),fig)
	Plot2Ddens(FieldE["Ez"],axis,-FieldsAmp,FieldsAmp,"Ez", "field", 1, SystemParameters,SubPlot(0,3,fig),fig)
	#Plot2Ddens(FieldB["Bx"],axis,-FieldsAmp+Bz0,FieldsAmp+Bz0,"By", "field", 1, SystemParameters,SubPlot(0,3,fig),fig)
	Plot2Ddens(FieldB["Bz"],axis,-FieldsAmp+Bz0,FieldsAmp+Bz0,"Bz", "field", 1, SystemParameters,SubPlot(0,4,fig),fig)

	gs_x=1
	for sort in SystemParameters['Particles'].keys():
		Plot2Ddens(np.abs(Dens[sort]),axis,0,densAmp*float(SystemParameters['Particles'][sort]['Dens']),
		sort + " density", "dens", 1, SystemParameters,SubPlot(gs_x,1,fig),fig)
		if sort != "Neutrals":
			Plot2Ddens(Pxx[sort],axis,0,0,sort + " Ppp", "dens", 0, SystemParameters,SubPlot(gs_x,2,fig),fig)
			Plot2Ddens(Pyy[sort],axis,0,0,sort + " Prr", "dens", 0, SystemParameters,SubPlot(gs_x,3,fig),fig)
			Plot2Ddens(Pzz[sort],axis,0,0,sort + " Pzz", "dens", 0, SystemParameters,SubPlot(gs_x,4,fig),fig)
			Plot2Ddens(J[sort]["Jx"],axis,0,0,sort + " Jx", "field", 0, SystemParameters,SubPlot(gs_x+1,2,fig),fig)
			Plot2Ddens(J[sort]["Jy"],axis,0,0,sort + " Jy", "field", 0, SystemParameters,SubPlot(gs_x+1,3,fig),fig)
			Plot2Ddens(J[sort]["Jz"],axis,0,0,sort + " Jz", "field", 0, SystemParameters,SubPlot(gs_x+1,4,fig),fig)
		gs_x+=2

	MaxTime=int(TimeStep)*tdelay
	MaxTimeDt=round(MaxTime,1)

	plt.figtext(0.5, 0.96,  r'$t\cdot\omega_p = '+str(MaxTimeDt)+'$',bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
						color='black',fontsize=18,ha='center')

	plt.tight_layout()
	plt.savefig('../Anime/2D/'+GraphName, format='png', dpi=150)
	plt.close(fig)

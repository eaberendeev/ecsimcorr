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

#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#ProcNum=comm.Get_size()
#print('ProcNum=',ProcNum)
##gs = gridspec.GridSpec(3, 4,height_ratios=[1,0.3,1],width_ratios=[1, 0.05,0.5,0.5])

#корневой каталог для расчёта (где файл параметров)
WorkDir='../'

FieldsAmp=0.03#амплитуда палитры для 3D полей
FieldsAmp2D=0.1#амплитуда палитры для 2D полей
densAmp = 1.5
#WindowSize=1000#окно для усреднения эффективности излучения
Bz0 = 0.2
#словарь с системными параметрами
SystemParameters=ReadParameters()

SystemParameters['2D_Diag_Fields']=['1','1','1','1','1','1'] #ReadParamFile(WorkDir)

FileNameSignNum=3 #len(str(int(int(SystemParameters['Max_Time'][0])/int(SystemParameters['Diagn_Time'][0]))))
print(FileNameSignNum)

def SubPlot(x,y,fig):
	return fig.add_subplot(gs[y,x])


tdelay = int(SystemParameters['TimeDelay2D']) * float( SystemParameters['Dt'])

try:
      os.mkdir('../Anime/2D')
except OSError:
      print("подкаталог для 2D диагностики существуют")


print(SystemParameters)

step = "0101"
stepX = "0101"
stepZ = "0003"
i = 1600
imax = 8000 
iStep = 100

while i <= imax:
	TimeStep=str(i).zfill(4)
	#имя файла результирующего
	GraphName='ResR2_'+TimeStep+'.png'
	#проверка на наличие уже построенного графика
	if(True):
		print('TimeStep=',TimeStep)

		fig = plt.figure(figsize=(int(SystemParameters['NumOfPartSpecies'])*4+25, 32))

		#на основе того, сколько было сортов частиц определить сетку для графиков
		#gs =gridspec.GridSpec(4,int(SystemParameters['NumOfPartSpecies'])+3,height_ratios=[0.1,1,1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей
		gs =gridspec.GridSpec(6,6,height_ratios=[0.1,1,1,1,1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей

		#считываем плотности
		DensDataXY= collections.OrderedDict()
		for sort in SystemParameters['Particles'].keys():
			DensDataXY.update(ReadDensFile2D(WorkDir,"Dens","PlaneZ",SystemParameters,sort,TimeStep))

		#PlotIntegMPI(WorkDir)
		FieldsTitles = ["Ex","Ey","Ez","Bx","By","Bz"]
		FieldDataXY=ReadFieldsFile2D(WorkDir,"PlaneZ_"+stepZ+"_","xy",SystemParameters,TimeStep)
		Name = WorkDir + "Fields/Diag2D/FieldPlaneZ_"+stepZ+"_"
		FieldDataXY = ReadFieldsFile2DNew(Name,FieldsTitles,TimeStep)

		Name = WorkDir + "Fields/Diag2D/FieldAvgPlaneZ_"
		FieldDataAvgXY = ReadFieldsFile2DNew(Name,FieldsTitles,TimeStep)

		JxyAvg = collections.OrderedDict()
		DensXyAvg = collections.OrderedDict()
		Ppp = collections.OrderedDict()
		Prr = collections.OrderedDict()
		for sort in SystemParameters['Particles'].keys():
			Name = WorkDir + "Particles/" + sort + "/Diag2D/CurrentPlaneAvgZ"
			JxyAvg[sort] = ReadFieldsFile2DNew(Name,["Jx","Jy","Jz"],TimeStep)
			Name = WorkDir + "Particles/" + sort + "/Diag2D/DensPlaneAvgZ"
			DensXyAvg[sort] = ReadFieldsFile2DNew(Name,["Dens"],TimeStep)
			Name = WorkDir + "Particles/" + sort + "/Diag2D/PyyPlaneAvgZ"
			Ppp[sort] = ReadFieldsFile2DNew(Name,["Ppp"],TimeStep)
			Name = WorkDir + "Particles/" + sort + "/Diag2D/PxxPlaneAvgZ"
			Prr[sort] = ReadFieldsFile2DNew(Name,["Prr"],TimeStep)
			
		bx, by = 0, 0
		ex, ey = FieldDataXY["Ex"].shape[0], FieldDataXY["Ex"].shape[1]
		COS, SIN  = init_COS_SIN(bx,ex,by,ey)
		ER, EA = vx_vy_to_vr_va(FieldDataXY["Ex"], FieldDataXY["Ey"], COS, SIN)
		rmap = init_rmap(bx,ex,by,ey)
		ER1Dr = phi_averaged(ER, rmap)
		EA1Dr = phi_averaged(EA, rmap)
		Bz1Dr = phi_averaged(FieldDataAvgXY["Bz"], rmap)

		Plot2Ddens(ER,"xy",-0.001,0.001,"Er", "field", 0, SystemParameters,SubPlot(-2,1,fig),fig)
		Plot2Ddens(EA,"xy",-0.001,0.001,"Ephi", "field", 0, SystemParameters,SubPlot(-2,2,fig),fig)
		Plot2Dfields(FieldDataAvgXY,"xy","Bz",-FieldsAmp+Bz0,FieldsAmp+Bz0,"Bz avg",SystemParameters,SubPlot(-2,3,fig),fig)

		Plot1D(ER1Dr,0,0,"Er",0,SystemParameters,SubPlot(-1,1,fig))
		Plot1D(EA1Dr,0,0,"Ephi",0,SystemParameters,SubPlot(-1,2,fig))
		Plot1D(Bz1Dr,0,0,"Bz",0,SystemParameters,SubPlot(-1,3,fig))

		gs_x=0
		for sort in SystemParameters['Particles'].keys():
			JR, JA = vx_vy_to_vr_va(JxyAvg[sort]["Jx"], JxyAvg[sort]["Jy"], COS, SIN)
			Plot2Ddens(np.abs(DensXyAvg[sort]["Dens"]),"xy",0,densAmp*float(SystemParameters['Particles'][sort]['Dens']),
			sort + " density", "dens", 1, SystemParameters,SubPlot(gs_x,1,fig),fig)
			Plot2Ddens(JR,"xy",-0.001,0.001,sort + " Jr", "field", 0, SystemParameters,SubPlot(gs_x,2,fig),fig)
			Plot2Ddens(JA,"xy",-0.001,0.001,sort + " Jphi", "field", 0, SystemParameters,SubPlot(gs_x,3,fig),fig)
			Plot2Ddens(Prr[sort]["Prr"],"xy",0,0,sort + " Prr", "dens", 0, SystemParameters,SubPlot(gs_x,4,fig),fig)
			Plot2Ddens(Ppp[sort]["Ppp"],"xy",0,0,sort + " Ppp", "dens", 0, SystemParameters,SubPlot(gs_x,5,fig),fig)

			data1D = phi_averaged(np.abs(DensXyAvg[sort]["Dens"]),rmap)
			Plot1D(data1D,0,densAmp*float(SystemParameters['Particles'][sort]['Dens']),
			sort + " density",0,SystemParameters,SubPlot(gs_x+1,1,fig))
			data1D = phi_averaged(JR,rmap)
			Plot1D(data1D,0,0,sort + " Jr",0,SystemParameters,SubPlot(gs_x+1,2,fig))
			data1D = phi_averaged(JA,rmap)
			Plot1D(data1D,0,0,sort + " Jphi",0,SystemParameters,SubPlot(gs_x+1,3,fig))
			data1D = phi_averaged(Prr[sort]["Prr"],rmap)
			Plot1D(data1D,0,0,sort + " Prr",0,SystemParameters,SubPlot(gs_x+1,4,fig))
			data1D = phi_averaged(Ppp[sort]["Ppp"],rmap)
			Plot1D(data1D,0,0,sort + " Ppp",0,SystemParameters,SubPlot(gs_x+1,5,fig))

			gs_x+=2


		MaxTime=int(TimeStep)*tdelay
		MaxTimeDt=round(MaxTime,1)


		plt.figtext(0.5, 0.96,  r'$t\cdot\omega_p = '+str(MaxTimeDt)+'$',bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
						color='black',fontsize=18,ha='center')

		plt.tight_layout()
		plt.savefig('../Anime/2D/'+GraphName, format='png', dpi=150)
		plt.close(fig)

	i=i+iStep


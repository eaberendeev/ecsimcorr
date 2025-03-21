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
        x0 = (ex-bx)/2
        COSP = np.zeros((ex-bx,ey-by))
        SINP = np.zeros((ex-bx,ey-by))
        for x in range(0,(ex-bx)):
                for y in range(0,(ey-by)):
                        r = ((x-x0)*(x-x0)+(y-x0)*(y-x0))**0.5
                        if r != 0:
                                cosp = (x-x0)/r
                                sinp = (y-x0)/r
                        else:
                                cosp = 0
                                sinp = 0
                        COSP[x,y] = cosp
                        SINP[x,y] = sinp

        return COSP, SINP

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

FieldsAmp=0.1#амплитуда палитры для 3D полей
FieldsAmp2D=0.1#амплитуда палитры для 2D полей
densAmp = 3.
#WindowSize=1000#окно для усреднения эффективности излучения
Bz0 = 0.
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

stepY = "200"
stepX = "070"
stepZ = "200"
i = 0
imax = 2999 
iStep = 50

while i <= imax:
	TimeStep=str(i).zfill(4)
	#имя файла результирующего
	GraphName='ResXY'+TimeStep+'.png'
	#проверка на наличие уже построенного графика
	if(True):
		print('TimeStep=',TimeStep)

		fig = plt.figure(figsize=(int(SystemParameters['NumOfPartSpecies'])*4+34, 42))

		#на основе того, сколько было сортов частиц определить сетку для графиков
		#gs =gridspec.GridSpec(4,int(SystemParameters['NumOfPartSpecies'])+3,height_ratios=[0.1,1,1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей
		gs =gridspec.GridSpec(7,6,height_ratios=[0.1,1,1,1,1,1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей

		#PlotIntegMPI(WorkDir)
		FieldsTitles = ["Ex","Ey","Ez"]
		#FieldDataXY=ReadFieldsFile2D(WorkDir,"PlaneZ_"+stepY+"_","xy",SystemParameters,TimeStep)
		Name = WorkDir + "Fields/Diag2D/FieldE_PlaneZ_"+stepY+"_"
		FieldE = ReadFieldsFile2DNew(Name,FieldsTitles,TimeStep)
		FieldsTitles = ["Bx","By","Bz"]
		#FieldDataXY=ReadFieldsFile2D(WorkDir,"PlaneZ_"+stepY+"_","xy",SystemParameters,TimeStep)
		Name = WorkDir + "Fields/Diag2D/FieldB_PlaneZ_"+stepY+"_"
		FieldB = ReadFieldsFile2DNew(Name,FieldsTitles,TimeStep)
		#Name = WorkDir + "Fields/Diag2D/FieldAvgPlaneZ_"
		#FieldDataAvgXY = ReadFieldsFile2DNew(Name,FieldsTitles,TimeStep)

		Jxz = collections.OrderedDict()
		DensXY = collections.OrderedDict()
		Pxx = collections.OrderedDict()
		Pyy = collections.OrderedDict()
		Pzz = collections.OrderedDict()
		for sort in SystemParameters['Particles'].keys():
			Name = WorkDir + "Particles/" + sort + "/Diag2D/Current_PlaneZ_"+stepY+"_"
			Jxz[sort] = ReadFieldsFile2DNew(Name,["Jx","Jy","Jz"],TimeStep)
			Name = WorkDir + "Particles/" + sort + "/Diag2D/Density_PlaneZ_"+stepY+"_"
			DensXY[sort] = ReadFieldsFile2DNew(Name,["Dens"],TimeStep)
			Name = WorkDir + "Particles/" + sort + "/Diag2D/Prr_PlaneZ_"+stepY+"_"
			Pxx[sort] = ReadFieldsFile2DNew(Name,["Pxx"],TimeStep)
			Name = WorkDir + "Particles/" + sort + "/Diag2D/Ppp_PlaneZ_"+stepY+"_"
			Pyy[sort] = ReadFieldsFile2DNew(Name,["Pyy"],TimeStep)
			Name = WorkDir + "Particles/" + sort + "/Diag2D/Pzz_PlaneZ_"+stepY+"_"
			Pzz[sort] = ReadFieldsFile2DNew(Name,["Pzz"],TimeStep)			

		bx, by = 0, 0
		ex, ey = FieldE["Ex"].shape[0], FieldE["Ex"].shape[1]
		COS, SIN  = init_COS_SIN(bx,ex,by,ey)
		ER, EA = vx_vy_to_vr_va(FieldE["Ex"], FieldE["Ey"], COS, SIN)
		rmap = init_rmap(bx,ex,by,ey)
		ER1Dr = phi_averaged(ER, rmap)
		EA1Dr = phi_averaged(EA, rmap)
		Ez1Dr = phi_averaged(FieldE["Ez"], rmap)
		#Br1Dr = phi_averaged(FieldB["Br"], rmap)
		Bz1Dr = phi_averaged(FieldB["Bz"], rmap)
			
		Plot2Ddens(ER,"xy",-FieldsAmp,FieldsAmp," Er", "field", 1, SystemParameters,SubPlot(0,1,fig),fig)
		Plot2Ddens(EA,"xy",-FieldsAmp,FieldsAmp," Ephi", "field", 1, SystemParameters,SubPlot(1,1,fig),fig)
		Plot2Ddens(FieldE["Ez"],"xy",-FieldsAmp,FieldsAmp," Ez", "field", 1, SystemParameters,SubPlot(2,1,fig),fig)
		#Plot2Ddens(FieldB["Bx"],"xy",-FieldsAmp,FieldsAmp," Br", "field", 1, SystemParameters,SubPlot(3,1,fig),fig)
		Plot2Ddens(FieldB["Bz"],"xy",-FieldsAmp,FieldsAmp," Bz", "field", 1, SystemParameters,SubPlot(3,1,fig),fig)

		Plot1D(ER1Dr,0,0,"Er",0,SystemParameters,SubPlot(0,2,fig))
		Plot1D(EA1Dr,0,0,"Ephi",0,SystemParameters,SubPlot(1,2,fig))
		Plot1D(Ez1Dr,0,0,"Br",0,SystemParameters,SubPlot(2,2,fig))
		Plot1D(Bz1Dr,0,0,"Bz",0,SystemParameters,SubPlot(3,2,fig))

		gs_y=3
		for sort in SystemParameters['Particles'].keys():
			Plot2Ddens(np.abs(DensXY[sort]),"xy",0,densAmp*float(SystemParameters['Particles'][sort]['Dens']),
			sort + " density", "dens", 1, SystemParameters,SubPlot(0,gs_y,fig),fig)
			data1D = phi_averaged(np.abs(DensXY[sort]),rmap)
			Plot1D(data1D,0,0,"Dens",0,SystemParameters,SubPlot(0,gs_y+1,fig))
			JR, JA = vx_vy_to_vr_va(Jxz[sort]["Jx"], Jxz[sort]["Jy"], COS, SIN)
			Plot2Ddens(JR,"xy",0,0,sort + " Jx", "field", 0, SystemParameters,SubPlot(1,gs_y,fig),fig)
			Plot2Ddens(JA,"xy",0,0,sort + " Jy", "field", 0, SystemParameters,SubPlot(2,gs_y,fig),fig)
			Plot2Ddens(Jxz[sort]["Jz"],"xy",0,0,sort + " Jz", "field", 0, SystemParameters,SubPlot(3,gs_y,fig),fig)
			data1D = phi_averaged(JR,rmap)

			Plot1D(data1D,0,0,"Jr",0,SystemParameters,SubPlot(1,gs_y+1,fig))
			data1D = phi_averaged(JA,rmap)

			Plot1D(data1D,0,0,"Jphi",0,SystemParameters,SubPlot(2,gs_y+1,fig))
			data1D = phi_averaged(Jxz[sort]["Jz"],rmap)
			Plot1D(data1D,0,0,"Jz",0,SystemParameters,SubPlot(3,gs_y+1,fig))
			Plot2Ddens(Pzz[sort],"xy",0,0,sort + " Pzz", "dens", 0, SystemParameters,SubPlot(4,gs_y,fig),fig)
			gs_y+=2


		MaxTime=int(TimeStep)*tdelay
		MaxTimeDt=round(MaxTime,1)


		plt.figtext(0.5, 0.96,  r'$t\cdot\omega_p = '+str(MaxTimeDt)+'$',bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
						color='black',fontsize=18,ha='center')

		plt.tight_layout()
		plt.savefig('../Anime/2D/'+GraphName, format='png', dpi=150)
		plt.close(fig)

	i=i+iStep


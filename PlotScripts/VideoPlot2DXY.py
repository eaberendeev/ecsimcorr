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


#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#ProcNum=comm.Get_size()
#print('ProcNum=',ProcNum)
##gs = gridspec.GridSpec(3, 4,height_ratios=[1,0.3,1],width_ratios=[1, 0.05,0.5,0.5])

#корневой каталог для расчёта (где файл параметров)
WorkDir='../'

FieldsAmp=0.1#амплитуда палитры для 3D полей
FieldsAmp2D=0.1#амплитуда палитры для 2D полей
densAmp = 3.
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

stepY = "005"
stepX = "0176"
stepZ = "0035"
i = 0
imax = 2999 
iStep = 1

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
		gs =gridspec.GridSpec(6,6,height_ratios=[0.1,1,1,1,1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей

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
			Name = WorkDir + "Particles/" + sort + "/Diag2D/Pxx_PlaneZ_"+stepY+"_"
			Pxx[sort] = ReadFieldsFile2DNew(Name,["Pxx"],TimeStep)
			Name = WorkDir + "Particles/" + sort + "/Diag2D/Pyy_PlaneZ_"+stepY+"_"
			Pyy[sort] = ReadFieldsFile2DNew(Name,["Pyy"],TimeStep)
			Name = WorkDir + "Particles/" + sort + "/Diag2D/Pzz_PlaneZ_"+stepY+"_"
			Pzz[sort] = ReadFieldsFile2DNew(Name,["Pzz"],TimeStep)			

			
		Plot2Ddens(FieldE["Ex"],"xy",-FieldsAmp,FieldsAmp," Ex", "field", 1, SystemParameters,SubPlot(0,1,fig),fig)
		Plot2Ddens(FieldE["Ey"],"xy",-FieldsAmp,FieldsAmp," Ey", "field", 1, SystemParameters,SubPlot(1,1,fig),fig)
		Plot2Ddens(FieldB["Bz"],"xy",-FieldsAmp,FieldsAmp," Bz", "field", 1, SystemParameters,SubPlot(2,1,fig),fig)

		FieldE1D = np.mean(FieldE["Ex"], axis=1)
		FieldE1D = np.diag(FieldE["Ex"])
		Plot1D(FieldE1D,0,0,"Ex",0,SystemParameters,SubPlot(0,2,fig))
		FieldE1D = np.mean(FieldE["Ey"], axis=1)
		FieldE1D = np.diag(FieldE["Ey"])
		Plot1D(FieldE1D,0,0,"Ey",0,SystemParameters,SubPlot(1,2,fig))
		FieldE1D = np.mean(FieldB["Bz"], axis=1)
		FieldE1D = np.diag(FieldB["Bz"])
		Plot1D(FieldE1D,0,0,"Bz",0,SystemParameters,SubPlot(2,2,fig))
		#Plot2Dfields(FieldE,"xz","Ex",-FieldsAmp,FieldsAmp,"Ex",SystemParameters,SubPlot(0,1,fig),fig)
		#Plot2Dfields(FieldE,"xz","Ey",-FieldsAmp,FieldsAmp,"Ey",SystemParameters,SubPlot(1,1,fig),fig)
		#Plot2Dfields(FieldDataXY,"xy","Ez",-FieldsAmp,FieldsAmp,"Ez",SystemParameters,SubPlot(-1,1,fig),fig)
		#Plot2Dfields(FieldDataXY,"xy","Bx",-FieldsAmp,FieldsAmp,"Bx",SystemParameters,SubPlot(-3,2,fig),fig)
		#Plot2Dfields(FieldDataXY,"xy","By",-FieldsAmp,FieldsAmp,"By",SystemParameters,SubPlot(-2,2,fig),fig)
		#Plot2Dfields(FieldDataXY,"xy","Bz",-FieldsAmp+Bz0,FieldsAmp+Bz0,"Bz",SystemParameters,SubPlot(-1,2,fig),fig)
		#Plot2Dfields(FieldDataAvgXY,"xy","Bx",-FieldsAmp,FieldsAmp,"Bx avg",SystemParameters,SubPlot(-3,3,fig),fig)
		#Plot2Dfields(FieldDataAvgXY,"xy","By",-FieldsAmp,FieldsAmp,"By avg",SystemParameters,SubPlot(-2,3,fig),fig)
		#Plot2Dfields(FieldDataXZ,"xz","Bz",-FieldsAmp+Bz0,FieldsAmp+Bz0,"Bz",SystemParameters,SubPlot(2,1,fig),fig)


		gs_y=2
		for sort in SystemParameters['Particles'].keys():
			Plot2Ddens(np.abs(DensXY[sort]),"xy",0,densAmp*float(SystemParameters['Particles'][sort]['Dens']),
			sort + " density", "dens", 1, SystemParameters,SubPlot(3,gs_y,fig),fig)
			FieldE1D = np.mean(np.abs(DensXY[sort]), axis=1)
			FieldE1D = np.diag(np.abs(DensXY[sort]))
			Plot1D(FieldE1D,0,0,"Dens",0,SystemParameters,SubPlot(3,gs_y+1,fig))
			#Plot2Ddens(Pxx[sort],"xz",0,0,sort + " Ppp", "dens", 0, SystemParameters,SubPlot(0,gs_y,fig),fig)
			#Plot2Ddens(Pyy[sort],"xz",0,0,sort + " Prr", "dens", 0, SystemParameters,SubPlot(1,gs_y,fig),fig)
			#Plot2Ddens(Pzz[sort],"xz",0,0,sort + " Pzz", "dens", 0, SystemParameters,SubPlot(2,gs_y,fig),fig)
			gs_y+=1
		for sort in SystemParameters['Particles'].keys():
			Plot2Ddens(Jxz[sort]["Jx"],"xz",0,0,sort + " Jx", "field", 0, SystemParameters,SubPlot(0,gs_y,fig),fig)
			Plot2Ddens(Jxz[sort]["Jy"],"xz",0,0,sort + " Jy", "field", 0, SystemParameters,SubPlot(1,gs_y,fig),fig)
			Plot2Ddens(Jxz[sort]["Jz"],"xz",0,0,sort + " Jz", "field", 0, SystemParameters,SubPlot(2,gs_y,fig),fig)
			FieldE1D = np.mean(Jxz[sort]["Jx"], axis=1)
			FieldE1D = np.diag(Jxz[sort]["Jx"])
			Plot1D(FieldE1D,0,0,"Jx",0,SystemParameters,SubPlot(0,gs_y+1,fig))
			FieldE1D = np.mean(Jxz[sort]["Jy"], axis=1)
			FieldE1D = np.diag(Jxz[sort]["Jy"])
			Plot1D(FieldE1D,0,0,"Jy",0,SystemParameters,SubPlot(1,gs_y+1,fig))
			FieldE1D = np.mean(Jxz[sort]["Jz"], axis=1)
			FieldE1D = np.diag(Jxz[sort]["Jz"])
			Plot1D(FieldE1D,0,0,"Jz",0,SystemParameters,SubPlot(2,gs_y+1,fig))
			gs_y+=1


		MaxTime=int(TimeStep)*tdelay
		MaxTimeDt=round(MaxTime,1)


		plt.figtext(0.5, 0.96,  r'$t\cdot\omega_p = '+str(MaxTimeDt)+'$',bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
						color='black',fontsize=18,ha='center')

		plt.tight_layout()
		plt.savefig('../Anime/2D/'+GraphName, format='png', dpi=150)
		plt.close(fig)

	i=i+iStep


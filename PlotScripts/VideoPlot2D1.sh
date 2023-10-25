#!/usr/bin/python3
#!/home/eberendeev/anaconda3/bin/python3
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

FieldsAmp=0.025#амплитуда палитры для 3D полей
FieldsAmp2D=0.1#амплитуда палитры для 2D полей
#WindowSize=1000#окно для усреднения эффективности излучения

#словарь с системными параметрами
SystemParameters=ReadParameters()

SystemParameters['2D_Diag_Fields']=['1','1','1','1','1','1'] #ReadParamFile(WorkDir)

FileNameSignNum=3 #len(str(int(int(SystemParameters['Max_Time'][0])/int(SystemParameters['Diagn_Time'][0]))))
print(FileNameSignNum)

def SubPlot(x,y,fig):
	return fig.add_subplot(gs[y,x])


tdelay = int(SystemParameters['TimeDelay2D']) * float( SystemParameters['Dt'])
#TimeStep='Last'#если это выбрано, то будет построен последний насчитанный файл

#if TimeStep=='Last':
	#files=os.listdir(WorkDir+'Fields/DiagXY/')#получили список файлов в каталоге с полями (они есть всегда)
	##print(files)
	#TimeStep=files[-1][-FileNameSignNum:]
try:
      os.mkdir('../Anime/2D')
except OSError:
      print("подкаталог для 2D диагностики существуют")

Fieldfiles=os.listdir(WorkDir+'Fields/Diag2D/')#получили список файлов в каталоге с полями (они есть всегда)
Graphfiles=os.listdir(WorkDir+'Anime/2D/')#получили список файлов в каталоге с графиками
pattern = "Field2Dxy"
Files = []
for f in Fieldfiles:
	if re.match(pattern,f):
		Files.append(f)


i =1 #len(Files)-1 #-rank #StartNum+  rank
print(i,Files)

print(SystemParameters)

      
while i >=0:
	TimeStep='100' #Files[i][-FileNameSignNum:]
	#имя файла результирующего
	GraphName='1D'+TimeStep+'.png'

	print('TimeStep=',TimeStep)

	fig = plt.figure(figsize=(16, 8))
		#fig1 = plt.figure(figsize=16, 10))

		#на основе того, сколько было сортов частиц определить сетку для графиков
	gs =gridspec.GridSpec(2,2,height_ratios=[0.1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей

		#PlotIntegMPI(WorkDir)
	FieldDataXY=ReadFieldsFile2D(WorkDir,"xy_201_","xy",SystemParameters,TimeStep)
	FieldDataXZ=ReadFieldsFile2D(WorkDir,"xz_201_","xz",SystemParameters,TimeStep)


	FieldDataC=ReadFieldsFile2D(WorkDir,"circle_","xy",SystemParameters,TimeStep)
	ax = SubPlot(0,1,fig)
	#ylist = np.swapaxes(FieldDataXY["Ex"],0,1)
	ylist = FieldDataXY["Ex"]
	ylist1 = ylist[:][550]
	#ax.set_ylim(-0.006,0.006)
	ax.plot(ylist1,linestyle=':',lw=2.4, zorder=4)
	ylist1 = ylist[:][750]
	#ax.set_ylim(-0.006,0.006)
	ax.plot(ylist1,linestyle=':',lw=2.4, zorder=4)
	ylist1 = ylist[:][650]
	#ax.set_ylim(-0.006,0.006)
	ax.plot(ylist1,linestyle=':',lw=2.4, zorder=4)		
	for k,item in enumerate(ylist1):
		print(k,item)

#	ylist1 = ylist[:][201-131]
#	ax.plot(ylist1,linestyle=':',lw=2.4, zorder=4)
#	ylist1 = ylist[:][201-132]
#	ax.plot(ylist1,linestyle=':',lw=2.4, zorder=4)
#	ylist2 = ylist[:][201+130]
#	ax.plot(ylist2,linestyle=':',lw=2.4, zorder=4)
                #axEff.legend()
	#Plot2Dfields(FieldDataC,"xz","Ex",-FieldsAmp,FieldsAmp,"Ex",SystemParameters,SubPlot(0,1,fig),fig)
	#Plot2Dfields(FieldDataC,"xz","By",-FieldsAmp,FieldsAmp,"By",SystemParameters,SubPlot(1,1,fig),fig)
	#Plot2Dfields(FieldDataC,"xz","Bz",-FieldsAmp,FieldsAmp,"Bz",SystemParameters,SubPlot(2,1,fig),fig)
	#Plot2Dfields2(FieldDataC,"xz","Bz",-FieldsAmp,FieldsAmp,"Pointing",SystemParameters,SubPlot(3,1,fig),fig)

	plt.tight_layout()
	plt.savefig('../Anime/2D/'+GraphName, format='png', dpi=150)
	plt.close(fig)



	i=i-1 #ProcNum
#plt.show()

#EOF

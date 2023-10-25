#!/home/fano.icmmg/berendeev_e_a/anaconda3/bin/python3
#!/usr/bin/python3
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

SystemParameters['3D_Diag_Fields']=['1','1','1','1','1','1'] #ReadParamFile(WorkDir)

FileNameSignNum=3 #len(str(int(int(SystemParameters['Max_Time'][0])/int(SystemParameters['Diagn_Time'][0]))))
print(FileNameSignNum)

def SubPlot(x,y,fig):
	return fig.add_subplot(gs[y,x])

delay = int(SystemParameters['TimeDelay3D'])//int(SystemParameters['TimeDelay2D'])
tdelay = int(SystemParameters['TimeDelay2D']) * float( SystemParameters['Dt'])
#TimeStep='Last'#если это выбрано, то будет построен последний насчитанный файл

#if TimeStep=='Last':
	#files=os.listdir(WorkDir+'Fields/DiagXY/')#получили список файлов в каталоге с полями (они есть всегда)
	##print(files)
	#TimeStep=files[-1][-FileNameSignNum:]

Fieldfiles=os.listdir(WorkDir+'Fields/Diag3D/')#получили список файлов в каталоге с полями (они есть всегда)
Graphfiles=os.listdir(WorkDir+'Anime/')#получили список файлов в каталоге с графиками


i =delay*len(Fieldfiles)-1 #-rank #StartNum+  rank
print(SystemParameters)

try:
	for sort in SystemParameters['Particles'].keys():
		os.mkdir('../Anime/'+sort)
except OSError:
      print("Подкаталоги для частиц существуют")


try:
      os.mkdir('../Anime/Fields')
except OSError:
      print("Подкаталог для полей существуют")

      
while i >=0:
	TimeStep=Fieldfiles[i//delay][-FileNameSignNum:]
	#имя файла результирующего
	GraphName='Res'+TimeStep+'.png'
	#проверка на наличие уже построенного графика
	if(GraphName not in Graphfiles):
		print('TimeStep=',TimeStep)

		fig = plt.figure(figsize=(int(SystemParameters['NumOfPartSpecies'])*4+8, 12))
		#fig1 = plt.figure(figsize=16, 10))

		#на основе того, сколько было сортов частиц определить сетку для графиков
		gs =gridspec.GridSpec(4,int(SystemParameters['NumOfPartSpecies'])+3,height_ratios=[0.1,1,1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей
		#gs1 =gridspec.GridSpec(2,1,height_ratios=[0.1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей

			
		#считываем плотности
		DensData= collections.OrderedDict()
		for sort in SystemParameters['Particles'].keys():
			DensData.update(ReadDensFile(WorkDir,SystemParameters,sort,TimeStep))
#
		PhaseData= collections.OrderedDict()
		for sort in SystemParameters['Particles'].keys():
			PhaseData.update(ReadPhaseFile(WorkDir,SystemParameters,sort,TimeStep))

		#PlotIntegMPI(WorkDir)
		FieldData=ReadFieldsFile(WorkDir,SystemParameters,TimeStep)
		#pprint.pprint(FieldData)
		#pprint.pprint(FieldData.keys())

		#3D плотности

		
		gs_x=0
		for sort in SystemParameters['Particles'].keys():
			Plot3Ddens(DensData[sort],"xz",sort,0,2*float(SystemParameters['Particles'][sort]['Dens']),sort,SystemParameters,SubPlot(gs_x,1,fig),fig)
			Plot3Ddens(DensData[sort],"yz",sort,0,2*float(SystemParameters['Particles'][sort]['Dens']),sort,SystemParameters,SubPlot(gs_x,2,fig),fig)
			Plot2DPhase(PhaseData,sort,0,2*float(SystemParameters['Particles'][sort]['Dens']),sort,SystemParameters,SubPlot(gs_x,3,fig),fig)
			#Plot1Ddens(DensData,sort,'long',0,4*float(SystemParameters['Particles'][sort]['Dens']),'',SystemParameters,fig.add_subplot(gs[gs_y,gs_x], label="1"),fig)
			
			gs_x+=1
		#gs_x=0
		#for sort in SystemParameters['Particles'].keys():
		#	gs_x+=1
		#gs_x=0
		#gs_y=2
		#for sort in SystemParameters['Particles'].keys():
		#	gs_x+=1
		#gs_x=0
		#gs_y=3
#		for sort in SystemParameters['Particles'].keys():
#			Plot1Ddens(DensData,sort,'long',320,0,4*float(SystemParameters['Particles'][sort]['Dens']),'',SystemParameters,fig.add_subplot(gs[gs_y,gs_x], label="1"),fig)
#			gs_x+=1

		#если есть правый пучок, то фазовый портрет его vx совмещаем с портретом от левого пучка
		#if 'RIGHTBEAM' in SystemParameters['Particles'].keys():
	#		PhaseData['LEFTBEAM']['vx']=PhaseData['LEFTBEAM']['vx']+PhaseData['RIGHTBEAM']['vx']
		#gs_x=0#поля начинаются со второго справа столбика
		#gs_y=1#первая строка
	#	for sort in SystemParameters['Particles'].keys():
	#		if sort!='RIGHTBEAM':#
#		Plot2Dphase(PhaseData,'BeamLeft','vx',0,0,'BeamLeft',SystemParameters,SubPlot(gs_x,2,fig),fig)
	#		Plot2Dphase(PhaseData,sort,'vy',0,0,sort,SystemParameters,SubPlot(gs_x,3,fig),fig)
	#		gs_x+=1


	#	gs_x=0
	#	for sort in SystemParameters['Particles'].keys():
	#		Plot1DField(FieldData,'Ex',int(Ny*0.5),-FieldsAmp2D,FieldsAmp2D,'',SystemParameters,fig.add_subplot(gs[-1,gs_x], label="2", frame_on=False),fig)
	#		gs_x+=1
	
	#	gs_x=-3 #поля начинаются со второго справа столбика
	#	gs_y=1 #первая строка
	#	for sort in FieldData.keys():
	#		if sort!= 'Bx':
     #                           Plot2Dfields(FieldData,sort,-FieldsAmp,FieldsAmp,sort,SystemParameters,SubPlot(gs_x,gs_y,fig),fig)
     #                           gs_y+=1
	#		if gs_y==4:
	#			gs_x+=1
	#			gs_y=1
	#		if sort== 'Bx':
	#			Plot2Dfields(FieldData,sort,float(SystemParameters['UniformB'][0])-FieldsAmp,float(SystemParameters['UniformB'][0])+FieldsAmp,sort,SystemParameters,SubPlot(gs_x,gs_y,fig),fig)
	#			gs_y+=1

		Plot3Dfields(FieldData,"xz","Ex",-FieldsAmp,FieldsAmp,"Ex",SystemParameters,SubPlot(-3,1,fig),fig)
		Plot3Dfields(FieldData,"xz","Ey",-FieldsAmp,FieldsAmp,"Ey",SystemParameters,SubPlot(-2,1,fig),fig)
		Plot3Dfields(FieldData,"xz","Ez",-FieldsAmp,FieldsAmp,"Ez",SystemParameters,SubPlot(-1,1,fig),fig)
		Plot3Dfields(FieldData,"xy","Ex",-FieldsAmp,FieldsAmp,"Ex",SystemParameters,SubPlot(-3,2,fig),fig)
		Plot3Dfields(FieldData,"xy","Ey",-FieldsAmp,FieldsAmp,"Ey",SystemParameters,SubPlot(-2,2,fig),fig)
		Plot3Dfields(FieldData,"xy","Ez",-FieldsAmp,FieldsAmp,"Ez",SystemParameters,SubPlot(-1,2,fig),fig)

		Plot3Dfields(FieldData,"xz","Bx",float(SystemParameters['UniformB'][0])-FieldsAmp,float(SystemParameters['UniformB'][0])+FieldsAmp,"Bx",SystemParameters,SubPlot(-3,3,fig),fig)
		Plot3Dfields(FieldData,"xz","By",-FieldsAmp,FieldsAmp,"By",SystemParameters,SubPlot(-2,3,fig),fig)
		Plot3Dfields(FieldData,"xz","Bz",-FieldsAmp,FieldsAmp,"Bz",SystemParameters,SubPlot(-1,3,fig),fig)

		MaxTime=int(TimeStep)*tdelay
		MaxTimeDt=round(MaxTime,1)



		plt.figtext(0.5, 0.96,  r'$t\cdot\omega_p = '+str(MaxTimeDt)+'$',bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
						color='black',fontsize=18,ha='center')

		#axEff.legend()

		plt.tight_layout()
		plt.savefig('../Anime/'+GraphName, format='png', dpi=150)
		plt.close(fig)



		#for sort in SystemParameters['Particles'].keys():
		#	fig = plt.figure(figsize=(12, 8))
		#	gs =gridspec.GridSpec(2,1,height_ratios=[0.1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей
		#	Plot2Ddens(DensData,sort,0,2*float(SystemParameters['Particles'][sort]['Dens']),sort,SystemParameters,SubPlot(0,1,fig),fig)
		#	plt.figtext(0.5, 0.96,  r'$t\cdot\omega_p = '+str(MaxTimeDt)+'$',bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
		#				color='black',fontsize=18,ha='center')
		#	plt.tight_layout()
		#	plt.savefig('../Anime/'+sort+'/Density'+TimeStep+'.png', format='png', dpi=150)
		#	plt.close(fig)
		
		#for sort in FieldData.keys():
		#	fig = plt.figure(figsize=(12, 8))
		#	gs =gridspec.GridSpec(2,1,height_ratios=[0.1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей
		#	if sort!= 'Bx':
        #                        Plot2Dfields(FieldData,sort,-FieldsAmp,FieldsAmp,sort,SystemParameters,SubPlot(0,1,fig),fig)
		#	if sort== 'Bx':
		#		Plot2Dfields(FieldData,sort,float(SystemParameters['UniformB'][0])-FieldsAmp,float(SystemParameters['UniformB'][0])+FieldsAmp,sort,SystemParameters,SubPlot(0,1,fig),fig)
		#	
		#	plt.figtext(0.5, 0.96,  r'$t\cdot\omega_p = '+str(MaxTimeDt)+'$',bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
		#				color='black',fontsize=18,ha='center')
		#	plt.tight_layout()
		#	plt.savefig('../Anime/Fields/'+sort+TimeStep+'.png', format='png', dpi=150)
		#	plt.close(fig)
	print(i,int(SystemParameters['TimeDelay3D'])//int(SystemParameters['TimeDelay2D']))
	i=i-delay #ProcNum
#plt.show()

#EOF

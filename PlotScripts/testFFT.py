# gratings.py
import numpy as np
import matplotlib.pyplot as plt

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

#корневой каталог для расчёта (где файл параметров)
WorkDir='../'

#словарь с системными параметрами
SystemParameters=ReadParameters()

SystemParameters['2D_Diag_Fields']=['1','1','1','1','1','1'] #ReadParamFile(WorkDir)

def SubPlot(x,y,fig):
    return fig.add_subplot(gs[y,x])


tdelay = int(SystemParameters['TimeDelay2D']) * float( SystemParameters['Dt'])

try:
      os.mkdir('../Anime/2D')
except OSError:
      print("подкаталог для 2D диагностики существуют")

Fieldfiles=os.listdir(WorkDir+'Fields/Diag2D/')#получили список файлов в каталоге с полями (они есть всегда)
Graphfiles=os.listdir(WorkDir+'Anime/2D/')#получили список файлов в каталоге с графиками


print(SystemParameters)

step = "016"
stepX = "016"
stepZ = "003"
t = 0 
tmax = 310
def write_fft_time(data):
    pass

fftTime = {}
harm_keys = []
for i in range(4):
    for j in range(4):
        harm_keys.append((i,j))
for key in harm_keys:
    fftTime[key] = []

print(fftTime)
while t <= tmax:
    TimeStep=str(t).zfill(3) #Files[i][-FileNameSignNum:]
    #имя файла результирующего
    GraphName='Res'+TimeStep+'.png'
    #проверка на наличие уже построенного графика
    if(True):
        print('TimeStep=',TimeStep)

        fig = plt.figure(figsize=(28, 32))

        #на основе того, сколько было сортов частиц определить сетку для графиков
        gs =gridspec.GridSpec(4,int(SystemParameters['NumOfPartSpecies'])+3,height_ratios=[0.1,1,1,1]) #первый аргумент - вертикальное количество, второй - горизонтальное (по горизонтале число сортов частиц + 2 столбца для полей

        FieldDataXY=ReadFieldsFile2D(WorkDir,"AvgPlaneZ_","xy",SystemParameters,TimeStep)
        #FieldDataXZ=ReadFieldsFile2D(WorkDir,"PlaneY_"+step+"_","xz",SystemParameters,TimeStep)
        #FieldDataYZ=ReadFieldsFile2D(WorkDir,"PlaneX_"+stepX+"_","yz",SystemParameters,TimeStep)

        FieldDataXY_BzFFT = np.fft.fft2(FieldDataXY["Bx"])
        for i in range(4):
            for j in range(4):
                fftTime[(i,j)].append(abs(FieldDataXY_BzFFT[i][j]))        

        # Plot2Dfields(FieldDataXZ,"xz","Bx",-FieldsAmp,FieldsAmp,"Ex",SystemParameters,SubPlot(-3,1,fig),fig)
        # Plot2Dfields(FieldDataXZ,"xz","By",-FieldsAmp,FieldsAmp,"Ey",SystemParameters,SubPlot(-2,1,fig),fig)
        # Plot2Dfields(FieldDataXZ,"xz","Bz",-FieldsAmp,FieldsAmp,"Ez",SystemParameters,SubPlot(-1,1,fig),fig)
        # Plot2Dfields(FieldDataXY,"xy","Bx",-FieldsAmp,FieldsAmp,"Ex",SystemParameters,SubPlot(-3,2,fig),fig)
        # Plot2Dfields(FieldDataXY,"xy","By",-FieldsAmp,FieldsAmp,"Ey",SystemParameters,SubPlot(-2,2,fig),fig)
        # Plot2Dfields(FieldDataXY,"xy","Bz",-FieldsAmp,FieldsAmp,"Ez",SystemParameters,SubPlot(-1,2,fig),fig)
        # Plot2Dfields(FieldDataYZ,"yz","Bx",-FieldsAmp,FieldsAmp,"Ex",SystemParameters,SubPlot(-3,3,fig),fig)
        # Plot2Dfields(FieldDataYZ,"yz","By",-FieldsAmp,FieldsAmp,"Ey",SystemParameters,SubPlot(-2,3,fig),fig)
        # Plot2Dfields(FieldDataYZ,"yz","Bz",-FieldsAmp,FieldsAmp,"Ez",SystemParameters,SubPlot(-1,3,fig),fig)

        # MaxTime=int(TimeStep)*tdelay
        # MaxTimeDt=round(MaxTime,1)

        # plt.figtext(0.5, 0.96,  r'$t\cdot\omega_p = '+str(MaxTimeDt)+'$',bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
        #                color='black',fontsize=18,ha='center')

        # plt.tight_layout()
        # plt.savefig('../Anime/2D/'+GraphName, format='png', dpi=150)
        # plt.close(fig)

    t=t+1

f = open("harmsCorrBx.txt", 'w')
for i in range(4):
    for j in range(4):
        f.write("("+str(i)+","+str(j) + ")  ")
f.write("\n")
for t in range(len(fftTime[(0,0)])):
    for i in range(4):
        for j in range(4):
            f.write(str(fftTime[(i,j)][t]) + " ") 
    f.write("\n")

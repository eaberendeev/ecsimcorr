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

WorkDir='../'

Bz0 = 0.2
#словарь с системными параметрами
SystemParameters=ReadParameters()

def SubPlot(x,y,fig):
	return fig.add_subplot(gs[y,x])

tdelay = int(SystemParameters['TimeDelay2D']) * float( SystemParameters['Dt'])

try:
      os.mkdir('../Anime')
except OSError:
      print("подкаталог Bz диагностики существуют")


print(SystemParameters)

i = 0
imax = 4300
iStep = 20
tBz = []
Bz = []
while i <= imax:
	TimeStep=str(i).zfill(4)
	print(TimeStep)
	Name = WorkDir + "Fields/Diag2D/FieldAvgPlaneZ_"
	FieldsTitles = ["Ex","Ey","Ez","Bx","By","Bz"]
	FieldDataAvgXY = ReadFieldsFile2DNew(Name,FieldsTitles,TimeStep)
	ex, ey = FieldDataAvgXY["Bz"].shape[0], FieldDataAvgXY["Bz"].shape[1]
	MaxTime=int(TimeStep)*tdelay
	MaxTimeDt=round(MaxTime,1)
	tBz.append(MaxTimeDt)
	Bz.append(FieldDataAvgXY["Bz"][ex//2][ey//2])
	i=i+iStep

f = open("Bz.txt", 'w')
f.write("time Bz")
f.write("\n")
for t in range(len(tBz)):
    f.write(str(tBz[t]) + " " + str(Bz[t])) 
    f.write("\n")
f.close()

fig = plt.figure(figsize=(12, 8))

gs =gridspec.GridSpec(1,1)

Plot1Dxy(tBz,Bz,0,0,"Bz",0,SubPlot(0,0,fig))
plt.tight_layout()
plt.savefig('../Anime/Bz.png', format='png', dpi=150)


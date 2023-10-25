
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as col
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
from matplotlib import rc
import struct
#f = open('rotJz', 'rb')

from ReadDataLib import *

#опеределенеи цветовых гамм
#для полей
cdict = {'red':	    ((0.0, 0.0, 0),
					 (0.35, 0.5859375, 0.5859375),
					 (0.4, 0, 0),#blue
					 (0.45, 0, 0),#0096ff
					 (0.495, 1, 1),#white
					 (0.505, 1, 1),#white
					 (0.55, 1, 1),#ff9600
					 (0.65, 1, 1),
                     (1.0, 0.5859375, 0.5859375)),
         'green':   ((0.0, 0.0, 0.0),
         		     (0.35, 0, 0),
         		     (0.4, 0, 0),#blue
         		     (0.45, 0.5859375, 0.5859375),#0096ff
         		     (0.495, 1, 1),#white
					 (0.505, 1, 1),#white
         		     (0.55, 0.5859375, 0.5859375),#ff9600
         		     (0.65, 0, 0),
                     (1.0, 0, 0)),
         'blue':    ((0.0, 0.5859375, 0.5859375),
                     (0.35, 1, 1),
                     (0.4, 1, 1),#blue
                     (0.45, 1, 1),#0096ff
                     (0.495, 1, 1),#white
					 (0.505, 1, 1),#white
                     (0.55, 0, 0),#ff9600
                     (0.65, 0.5859375, 0.5859375),
                     (1.0, 0, 0))}
#фазовая
phasecdict = {'red':(
(0.0, 1, 1),#white
(0.15, 0, 0),#0096ff
(0.35, 0, 0),#blue
(0.55, 1, 1),#ff9600
(0.75, 1, 1),
(1.0, 0.5859375, 0.5859375)),

'green':   (
(0.0, 1, 1),#white
(0.15, 0.5859375, 0.5859375),#0096ff
(0.35, 0, 0),#blue
(0.55, 0.5859375, 0.5859375),#ff9600
(0.75, 0, 0),
(1.0, 0, 0)),

'blue':    (
(0.0, 1, 1),#white
(0.15, 1, 1),#0096ff
(0.35, 1, 1),#blue
(0.55, 0, 0),#ff9600
(0.75, 0.5859375, 0.5859375),
(1.0, 0, 0))}

def Plot1Dxy(x,data,Vmin,Vmax,title,scale,ax):
	if(scale == 1):
		ax.set_ylim(Vmin,Vmax)
	ax.set_title(title)

	ax.plot(x,data, color="g",linestyle=':',lw=2.4, zorder=4)

def Plot1D(data,Vmin,Vmax,title,scale,SystemParameters,ax):
	Dx = float(SystemParameters['Dx'])
	xlist = np.arange(0,data.size)
	if(scale == 1):
		ax.set_ylim(Vmin,Vmax)
	ax.set_title(title)

	ax.plot(Dx*xlist,data, color="g",linestyle=':',lw=2.4, zorder=4)


def Plot1Ddens(DensData,Sort,Direction,Vmin,Vmax,title,SystemParameters,ax,fig):
	#Direction - 'long' или 'trans'
	#Coord - координата на которой делать срез Игрик для long и Икс для trans
	#if(Direction=='long'):
	#	Data1D=DensData[Sort][:][Coord]
	#	ax.set_xlim(int(SystemParameters['X_DampCellsLeft']), int(SystemParameters['Nx'])-int(SystemParameters['X_DampCellsRight']))
	#if(Direction=='trans'):#похоже что не работает
	#	print('!!',Direction)
	Coord=int(SystemParameters['NumCellsX'])
	#Data1D=DensData[Sort].np.swapaxes(Dens[PartSort],0,1)
	Data2D=np.swapaxes(DensData[Sort],0,1)
#	Data2D=DensData[Sort]
	Data1D=Data2D[:][Coord//2+2]


#        Dens[PartSort]=data.reshape(( Nr, Nz))

	#	ax.set_xlim(0, int(SystemParameters['Ny']))
	xlist = np.arange(0,int(SystemParameters['NumCellsY']))
	ylist = [Data1D[x] for x in xlist]
#       print(ylist)
	ax.plot(xlist,ylist, color="g",linestyle=':',lw=2.4, zorder=4)
	Data1D=Data2D[:][Coord//2-2]


#        Dens[PartSort]=data.reshape(( Nr, Nz))

	#	ax.set_xlim(0, int(SystemParameters['Ny']))
	xlist = np.arange(0,int(SystemParameters['NumCellsY']))
	ylist = [Data1D[x] for x in xlist]
#       print(ylist)
	ax.plot(xlist,ylist, color="r",linestyle=':',lw=2.4, zorder=4)
	Data1D=Data2D[:][2*Coord//3]


#        Dens[PartSort]=data.reshape(( Nr, Nz))

	#	ax.set_xlim(0, int(SystemParameters['Ny']))
	xlist = np.arange(0,int(SystemParameters['NumCellsY']))
	ylist = [Data1D[x] for x in xlist]
#       print(ylist)
	ax.plot(xlist,ylist, color="b",linestyle=':',lw=2.4, zorder=4)

	#ax.plot(Data1D, color="g")
	#ax.set_xlim(0,50)
	ax.set_ylim(0,Vmax)
	ax.tick_params(axis='y', colors="g")
	ax.set_xlabel('$x$')
	#ax.set_ylabel('dens, r = '+str(Coord*0.025), color="g")
	ax.yaxis.grid() # horizontal lines

	return 0

def Plot1DField(FieldData,Type,Coord,Vmin,Vmax,title,SystemParameters,ax,fig):
	#Direction - 'long' или 'trans'
	#Coord - координата на которой делать срез Игрик для long и Икс для trans
	Data1D=FieldData[Type][:][Coord]
	ax.set_xlim(int(SystemParameters['X_DampCellsLeft']), int(SystemParameters['Nx'])-int(SystemParameters['X_DampCellsRight']))
	ax.plot(Data1D, color="C0",lw=0.4,label=Type)
	ax.set_ylim(Vmin,Vmax)
	ax.yaxis.grid() # horizontal lines

	#ax.set_ylabel(Type,color="C0")
	ax.yaxis.tick_right()
	ax.tick_params(
			axis='x',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom='off',      # ticks along the bottom edge are off
			top='off',         # ticks along the top edge are off
			labelbottom='off') # labels along the bottom edge are off
	ax.yaxis.set_label_position('right')
	ax.tick_params(axis='y', colors="C0")
	ax.legend()
	return 0


def Plot2DRad(DensData,plane,title,SystemParameters,ax,fig):
	cm = col.LinearSegmentedColormap('phase',phasecdict,N=1024,gamma=1)
	sizeX = int(SystemParameters['NumCellsX'])*float(SystemParameters['Dx'])
	sizeY = int(SystemParameters['NumCellsY'])*float(SystemParameters['Dy'])
	sizeZ = int(SystemParameters['NumCellsZ'])*float(SystemParameters['Dz'])
	dampSizeX = int(SystemParameters['DampCellsX'])*float(SystemParameters['Dx'])
	dampSizeY = int(SystemParameters['DampCellsY'])*float(SystemParameters['Dy'])
	dampSizeZ = int(SystemParameters['DampCellsZ'])*float(SystemParameters['Dz'])
	shape = DensData.shape

	if plane == "xy":
		ax.set_xlabel('$x$')
		ax.set_ylabel('$y$')
		ax.set_xlim(0,sizeX-2*dampSizeX)
		ax.set_ylim(0,sizeY-2*dampSizeX)
		ex = [-dampSizeX,sizeX-dampSizeX,-dampSizeY,sizeY - dampSizeY ]	
		data = np.swapaxes(DensData,0,1)
		im=ax.imshow(data, cmap=cm,origin='lower',extent=ex,aspect='auto')
	elif plane == "xz":
		ax.set_xlabel('$x$')
		ax.set_ylabel('$z$')
		ax.set_xlim(0,sizeX-2*dampSizeX)
		ax.set_ylim(0,sizeZ-2*dampSizeZ)
		ex = [-dampSizeX,sizeX-dampSizeX,-dampSizeZ,sizeZ - dampSizeZ ]		
		data = np.swapaxes(DensData,0,1)
		im=ax.imshow(data, cmap=cm,origin='lower',extent=ex,aspect='auto')
	elif plane == "circle":
		ax.set_xlabel('$x$')
		ax.set_ylabel('$angle, graduses$')
		ax.set_xlim(0,sizeX-2*dampSizeX)
		ax.set_ylim(0,360)
		ex = [-dampSizeX,sizeX-dampSizeX,0,360 ]	
		data = np.swapaxes(DensData,0,1)
		im=ax.imshow(data, cmap=cm,origin='lower',extent=ex,aspect='auto')
	else :
		return

#	im=ax.imshow(DensData[Sort], vmin=Vmin, vmax=Vmax, cmap='Greys',origin='lower',aspect='auto')

	ax.set_title(title)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("bottom", size="8%", pad=0.6)#свойства colorbar'а
	cbar =fig.colorbar(im, cax=cax, orientation='horizontal')
	#ax.set_xlim(int(SystemParameters['Z_DampCellsLeft'][0]), int(SystemParameters['Nz'][0])-int(SystemParameters['Z_DampCellsRight'][0]))
#	ax.set_ylim(-1,25)
	#ax.set_xlim(0,sizeZ-2*dampSizeZ)
	#xlist = np.arange(0,200,0.1)
	#ylist = [6.7 for x in xlist]
#	print(ylist)
#	ax.plot(xlist,ylist, color="g",linestyle=':',lw=2.4, zorder=4)
	tick_locator = ticker.MaxNLocator(nbins=5)
	cbar.locator = tick_locator
	cbar.formatter.set_powerlimits((0, 0))
	locator=ticker.MaxNLocator(prune=None, nbins=5)
	cbar.ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_major_locator(locator)
	locator=ticker.MaxNLocator(prune=None, nbins=5)
	ax.xaxis.set_major_locator(locator)
	cbar.update_ticks()

def Plot2Ddens(DensData,plane,Vmin,Vmax,title,cmap,scale,SystemParameters,ax,fig):
	if cmap == 'dens':
		cm = col.LinearSegmentedColormap('phase',phasecdict,N=1024,gamma=1)
	elif cmap == 'field':
		cm = col.LinearSegmentedColormap('fields',cdict,N=1024,gamma=1)
	else:
		print("Invalid colormap name!\n")
		exit(0)
	sizeX = int(SystemParameters['NumCellsX'])*float(SystemParameters['Dx'])
	sizeY = int(SystemParameters['NumCellsY'])*float(SystemParameters['Dy'])
	sizeZ = int(SystemParameters['NumCellsZ'])*float(SystemParameters['Dz'])
	dampSizeX = 0 #int(SystemParameters['DampCellsX'])*float(SystemParameters['Dx'])
	shape = DensData.shape

	if plane == "xy":
		ax.set_xlabel('$x$')
		ax.set_ylabel('$y$')
		#ax.set_xlim(150,sizeX-2*dampSizeX-150)
		ex = [-dampSizeX,sizeX-dampSizeX,0,sizeY]	
	
	if plane == "xz":
		ax.set_xlabel('$x$')
		ax.set_ylabel('$z$')
		#ax.set_xlim(150,sizeX-2*dampSizeX-150)		
		ex = [-dampSizeX,sizeX-dampSizeX,0,sizeZ]

	if plane == "yz":
		ax.set_xlabel('$y$')
		ax.set_ylabel('$z$')
		ax.set_xlim(0,sizeY)		
		ex = [0,sizeY,0,sizeZ]

	data = np.swapaxes(DensData,0,1)

	if scale != 0:
		im=ax.imshow(data, vmin=Vmin, vmax=Vmax,  cmap=cm,origin='lower',extent=ex,aspect='auto')
	else:
		im=ax.imshow(data, cmap=cm,origin='lower',extent=ex,aspect='auto')

	ax.set_title(title)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("bottom", size="8%", pad=0.6)#свойства colorbar'а
	cbar =fig.colorbar(im, cax=cax, orientation='horizontal')

	tick_locator = ticker.MaxNLocator(nbins=5)
	cbar.locator = tick_locator
	cbar.formatter.set_powerlimits((0, 0))
	locator=ticker.MaxNLocator(prune=None, nbins=5)
	cbar.ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_major_locator(locator)
	locator=ticker.MaxNLocator(prune=None, nbins=5)
	ax.xaxis.set_major_locator(locator)
	cbar.update_ticks()

def Plot3Ddens(DensData,plane,Sort,Vmin,Vmax,title,SystemParameters,ax,fig):
	cm = col.LinearSegmentedColormap('phase',phasecdict,N=1024,gamma=1)
	sizeX = int(SystemParameters['NumCellsX'])*float(SystemParameters['Dx'])
	sizeY = int(SystemParameters['NumCellsY'])*float(SystemParameters['Dy'])
	sizeZ = int(SystemParameters['NumCellsZ'])*float(SystemParameters['Dz'])
	dampSizeX = int(SystemParameters['DampCellsX'])*float(SystemParameters['Dx'])
	shape = DensData.shape

	if plane == "xy":
		ax.set_xlabel('$x$')
		ax.set_ylabel('$y$')
		ax.set_xlim(0,sizeX-2*dampSizeX)
		ex = [-dampSizeX,sizeX-dampSizeX,0,sizeY]
		data = np.zeros((shape[0],shape[1]) )
		for i in range(shape[0]):
			for j in range(shape[1]):
				data[i][j] = DensData[i][j][shape[2]//2]	
	
	if plane == "xz":
		ax.set_xlabel('$x$')
		ax.set_ylabel('$z$')
		ex = [-dampSizeX,sizeX-dampSizeX,0,sizeZ]

		data = np.zeros( (shape[0],shape[2]) )

		for i in range(shape[0]):
			for j in range(shape[2]):
				data[i][j] = DensData[i][shape[1]//2][j]
	

	if plane == "yz":
		ex = [0,sizeY,0,sizeZ]
		ax.set_xlabel('$y$')
		ax.set_ylabel('$z$')	
		data = np.zeros( (shape[1],shape[2]) )
		for i in range(shape[1]):
			for j in range(shape[2]):
				data[i][j] = DensData[shape[0]//2][i][j]		

	data = np.swapaxes(data,0,1)
	#if plane == "yz":
	#	data = np.swapaxes(DensData,0,1)
	#	data = np.swapaxes(data,1,2)
	
	#data = data[:][:][shape[2]//2 ]




	im=ax.imshow(data, vmin=Vmin, vmax=Vmax,  cmap=cm,origin='lower',extent=ex,aspect='auto')
#	im=ax.imshow(DensData[Sort], vmin=Vmin, vmax=Vmax, cmap='Greys',origin='lower',aspect='auto')

	ax.set_title(title)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("bottom", size="8%", pad=0.6)#свойства colorbar'а
	cbar =fig.colorbar(im, cax=cax, orientation='horizontal')
	#ax.set_xlim(int(SystemParameters['Z_DampCellsLeft'][0]), int(SystemParameters['Nz'][0])-int(SystemParameters['Z_DampCellsRight'][0]))
#	ax.set_ylim(-1,25)
	#ax.set_xlim(0,sizeZ-2*dampSizeZ)
	xlist = np.arange(0,200,0.1)
	ylist = [6.7 for x in xlist]
#	print(ylist)
#	ax.plot(xlist,ylist, color="g",linestyle=':',lw=2.4, zorder=4)
	tick_locator = ticker.MaxNLocator(nbins=5)
	cbar.locator = tick_locator
	cbar.formatter.set_powerlimits((0, 0))
	locator=ticker.MaxNLocator(prune=None, nbins=5)
	cbar.ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_major_locator(locator)
	locator=ticker.MaxNLocator(prune=None, nbins=5)
	ax.xaxis.set_major_locator(locator)
	cbar.update_ticks()

def Plot2DPhase(DensData,Sort,Vmin,Vmax,title,SystemParameters,ax,fig):
	cm = col.LinearSegmentedColormap('phase',phasecdict,N=1024,gamma=1)
	sizeX = int(SystemParameters['NumCellsX'])*float(SystemParameters['Dx'])
	px_min = float(SystemParameters['Particles'][Sort]['Px_min'])
	px_max = float(SystemParameters['Particles'][Sort]['Px_max'])
	#sizeY = int(SystemParameters['NumCellsY'])*float(SystemParameters['Dy'])
	dampSizeX = 0 # int(SystemParameters['DampCellsX'])*float(SystemParameters['Dx'])
	im=ax.imshow(DensData[Sort], cmap=cm,origin='lower',extent=[0,sizeX-2*dampSizeX,px_min,px_max],aspect='auto')
#	im=ax.imshow(DensData[Sort], vmin=Vmin, vmax=Vmax, cmap='Greys',origin='lower',aspect='auto')
	ax.set_xlabel('$x$')
	ax.set_ylabel('$p$')
	ax.set_title(title)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("bottom", size="8%", pad=0.6)#свойства colorbar'а
	cbar =fig.colorbar(im, cax=cax, orientation='horizontal')
	#ax.set_xlim(int(SystemParameters['Z_DampCellsLeft'][0]), int(SystemParameters['Nz'][0])-int(SystemParameters['Z_DampCellsRight'][0]))
#	ax.set_ylim(-1,25)
	#ax.set_xlim(0,sizeZ-2*dampSizeZ)
	xlist = np.arange(0,200,0.1)
	ylist = [6.7 for x in xlist]
#	print(ylist)
#	ax.plot(xlist,ylist, color="g",linestyle=':',lw=2.4, zorder=4)
	tick_locator = ticker.MaxNLocator(nbins=5)
	cbar.locator = tick_locator
	cbar.formatter.set_powerlimits((0, 0))
	locator=ticker.MaxNLocator(prune=None, nbins=5)
	cbar.ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_major_locator(locator)
	locator=ticker.MaxNLocator(prune=None, nbins=5)
	ax.xaxis.set_major_locator(locator)
	cbar.update_ticks()
	
	
def Plot2Dphase(PhaseData,Sort,MomType,Vmin,Vmax,title,SystemParameters,ax,fig):
	cm = col.LinearSegmentedColormap('phase',phasecdict,N=1024,gamma=1)
#	im=ax.imshow(PhaseData[Sort][MomType], vmin=2*int(SystemParameters['PartNumPerCell'][0]), vmax=Vmax, cmap=cm,origin='lower',aspect='auto')
	im=ax.imshow(PhaseData[Sort][MomType], cmap=cm,origin='lower',aspect='auto')
	ax.set_xlabel('$x$')
	#y=[0,int(SystemParameters['PhaseDiagNum'][0])]
	if(MomType=='vz'):
		ax.set_ylabel('$v_x$')
		labels = [SystemParameters['Particles'][Sort]['PhaseMinX_Y'][0], SystemParameters['Particles'][Sort]['PhaseMaxX_Y'][0]]
		if(float(SystemParameters['Particles'][Sort]['PhaseMinX_Y'][0])<0):
			dv=(float(SystemParameters['Particles'][Sort]['PhaseMaxX_Y'][0])-float(SystemParameters['Particles'][Sort]['PhaseMinX_Y'][0]))/float(SystemParameters['PhaseDiagNum'][0])
			labels.append('0')
			y.append(-float(SystemParameters['Particles'][Sort]['PhaseMinX_Y'][0])/dv)
		if(Sort[-4:]=='BEAM'):
			vb=float(SystemParameters['Particles'][Sort]['Velocity'][0])
			dv=(float(SystemParameters['Particles'][Sort]['PhaseMaxX_Y'][0])-float(SystemParameters['Particles'][Sort]['PhaseMinX_Y'][0]))/float(SystemParameters['PhaseDiagNum'][0])
			vby=vb/dv-float(SystemParameters['Particles'][Sort]['PhaseMinX_Y'][0])/dv
			#print('!!!',vb,dv,vby)
			labels.append(str(round(vb,2)))
			y.append(vby)
			#x1, y1 = [0, int((SystemParameters['PhaseDiagNum'][0]))], [vby, vby]
			#plt.plot(x1, y1, '--')
	if(MomType=='vy'):
		ax.set_ylabel('$v_y$')
		labels = [SystemParameters['Particles'][Sort]['PhaseMinX_Y'][1], SystemParameters['Particles'][Sort]['PhaseMaxX_Y'][1]]

		if(float(SystemParameters['Particles'][Sort]['PhaseMinX_Y'][1])<0):
			dv=(float(SystemParameters['Particles'][Sort]['PhaseMaxX_Y'][1])-float(SystemParameters['Particles'][Sort]['PhaseMinX_Y'][1]))/float(SystemParameters['PhaseDiagNum'][0])
			labels.append('0')
			y.append(-float(SystemParameters['Particles'][Sort]['PhaseMinX_Y'][1])/dv)
	
	#plt.yticks(y, labels)
	ax.yaxis.grid() # horizontal lines
	#coordConvert=float(SystemParameters['PhaseDiagCoordNum'][0])/SystemParameters['Nx']
	#ax.set_xlim(float(SystemParameters['X_DampCellsLeft'][0])*coordConvert, float(SystemParameters['PhaseDiagCoordNum'][0])-float(SystemParameters['X_DampCellsRight'][0])*coordConvert)
	divider = make_axes_locatable(ax)
	ax.set_ylim(0, 4)
	cax = divider.append_axes("bottom", size="2%", pad=0.6)#свойства colorbar'а
	cbar =fig.colorbar(im, cax=cax, orientation='horizontal')
	#ax.set_xlim(int(SystemParameters['X_DampCellsLeft'][0]), int(SystemParameters['Nx'])-int(SystemParameters['X_DampCellsRight'][0]))
	#ax.set_ylim(int(SystemParameters['Y_DampCellsBottom'][0]), int(SystemParameters['Ny'])-int(SystemParameters['Y_DampCellsTop'][0]))
	tick_locator = ticker.MaxNLocator(nbins=5)
	cbar.locator = tick_locator
	cbar.formatter.set_powerlimits((0, 0))
	cbar.ax.xaxis.set_ticks_position('bottom')
	cbar.update_ticks()

def Plot2Dfields(FieldData,plane,Type,Vmin,Vmax,title,SystemParameters,ax,fig):

	sizeX = int(SystemParameters['NumCellsX'])*float(SystemParameters['Dx'])
	sizeY = int(SystemParameters['NumCellsY'])*float(SystemParameters['Dy'])
	sizeZ = int(SystemParameters['NumCellsZ'])*float(SystemParameters['Dz'])
	dampSizeX = 0 #int(SystemParameters['DampCellsX'])*float(SystemParameters['Dx'])
	cm = col.LinearSegmentedColormap('fields',cdict,N=1024,gamma=1)
	
	data = FieldData[Type]

	if plane == "xy":
		ax.set_xlabel('$x$')
		ax.set_ylabel('$y$')
		ex = [-dampSizeX,sizeX-dampSizeX,0,sizeY]
		#ax.set_xlim(150,sizeX-2*dampSizeX-150)
	
	if plane == "xz":
		ax.set_xlabel('$x$')
		ax.set_ylabel('$z$')
		ex = [-dampSizeX,sizeX-dampSizeX,0,sizeZ]
		#ax.set_xlim(150,sizeX-2*dampSizeX-150)
	if plane == "yz":
		ax.set_xlabel('$y$')
		ax.set_ylabel('$z$')
		ex = [0,sizeY,0,sizeZ]
		#ax.set_xlim(0,sizeX-2*dampSizeX)

	data = np.swapaxes(data,0,1)

	im=ax.imshow(data, vmin=Vmin, vmax=Vmax, cmap=cm,origin='lower',extent=ex,aspect='auto')
	#ax.xaxis.set_ticks(x)
	#ax.yaxis.set_ticks(y)
	locator=ticker.MaxNLocator(prune=None, nbins=5)
	ax.yaxis.set_major_locator(locator)
	locator=ticker.MaxNLocator(prune=None, nbins=5)
	ax.xaxis.set_major_locator(locator)
	#ax.set_title(title)
#	ax.set_ylim(-1,25)
	#ax.set_xlim(0,sizeZ-2*dampSizeZ)
	#ax.set_xlabel('$x$')
	#ax.set_ylabel('$y$')
	ax.set_title(title)
	#ax.set_xlim(int(SystemParameters['X_DampCellsLeft'][0]), int(SystemParameters['Nx'][0])-int(SystemParameters['X_DampCellsRight'][0]))
	#ax.set_xlim(int(SystemParameters['Z_DampCellsLeft'][0]), int(SystemParameters['Nz'][0])-int(SystemParameters['Z_DampCellsRight'][0]))
	#ax.set_ylim(-10, int(SystemParameters['Nr'][0])-300)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("bottom", size="8%", pad=0.6)#свойства colorbar'а
	cbar =fig.colorbar(im, cax=cax, orientation='horizontal')
	tick_locator = ticker.MaxNLocator(nbins=3)
	cbar.locator = tick_locator
		#cb = fig.colorbar(im,cax=cax, orientation='horizontal')#,ticks=[-1e-6,-0.5e-6,0,0.5e-6,1e-6]
	#cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
	#cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])  # horizontal colorbar
	#cb = ax.colorbar(im,cax=cax, orientation='horizontal')#,ticks=[-1e-6,-0.5e-6,0,0.5e-6,1e-6]
	#cb.outline.set_linewidth(0.2)#ширина border'а у colorbar'а
	cbar.formatter.set_powerlimits((0, 0))
	#cb.format=ticker.FuncFormatter(fmt)	
	cbar.ax.xaxis.set_ticks_position('bottom')
	cbar.update_ticks()


def Plot2Dfields2(FieldData,plane,Type,Vmin,Vmax,title,SystemParameters,ax,fig):

	sizeX = int(SystemParameters['NumCellsX'])*float(SystemParameters['Dx'])
	sizeY = int(SystemParameters['NumCellsY'])*float(SystemParameters['Dy'])
	sizeZ = int(SystemParameters['NumCellsZ'])*float(SystemParameters['Dz'])
	dampSizeX = 0 # int(SystemParameters['DampCellsX'])*float(SystemParameters['Dx'])
	#cm = col.LinearSegmentedColormap('phase',cdict,N=1024,gamma=1)
	cm = col.LinearSegmentedColormap('phase',phasecdict,N=1024,gamma=1)

	
	data = FieldData#[Type]
	data1 = FieldData[Type]


	ex = [-dampSizeX,sizeX-dampSizeX,0,360]



	for i in range(data["Ex"].shape[0]):
		for j in range(data["Ex"].shape[1]):
			ry = 13*np.cos(2*3.1415926*j/data["Ex"].shape[1])
			rz = 13*np.sin(2*3.1415926*j/data["Ex"].shape[1])
			Bp = -rz/13*data["By"][i][j] + ry/13*data["Bz"][i][j]
			Ep = -rz/13*data["Ey"][i][j] + ry/13*data["Ez"][i][j]
			data1[i][j] = abs(data["Ex"][i][j]*Bp - data["Bx"][i][j]*Ep)

	data1 = np.swapaxes(data1,0,1)
	ax.set_xlim(0,sizeX-2*dampSizeX)
	im=ax.imshow(data1,  cmap=cm,origin='lower',extent=ex,aspect='auto')
	#ax.xaxis.set_ticks(x)
	#ax.yaxis.set_ticks(y)
	locator=ticker.MaxNLocator(prune=None, nbins=5)
	ax.yaxis.set_major_locator(locator)
	locator=ticker.MaxNLocator(prune=None, nbins=5)
	ax.xaxis.set_major_locator(locator)
	#ax.set_title(title)
#	ax.set_ylim(-1,25)
	#ax.set_xlim(0,sizeZ-2*dampSizeZ)
	#ax.set_xlabel('$x$')
	#ax.set_ylabel('$y$')
	ax.set_title(title)
	#ax.set_xlim(int(SystemParameters['X_DampCellsLeft'][0]), int(SystemParameters['Nx'][0])-int(SystemParameters['X_DampCellsRight'][0]))
	#ax.set_xlim(int(SystemParameters['Z_DampCellsLeft'][0]), int(SystemParameters['Nz'][0])-int(SystemParameters['Z_DampCellsRight'][0]))
	#ax.set_ylim(-10, int(SystemParameters['Nr'][0])-300)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("bottom", size="8%", pad=0.6)#свойства colorbar'а
	cbar =fig.colorbar(im, cax=cax, orientation='horizontal')
	tick_locator = ticker.MaxNLocator(nbins=3)
	cbar.locator = tick_locator
		#cb = fig.colorbar(im,cax=cax, orientation='horizontal')#,ticks=[-1e-6,-0.5e-6,0,0.5e-6,1e-6]
	#cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
	#cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])  # horizontal colorbar
	#cb = ax.colorbar(im,cax=cax, orientation='horizontal')#,ticks=[-1e-6,-0.5e-6,0,0.5e-6,1e-6]
	#cb.outline.set_linewidth(0.2)#ширина border'а у colorbar'а
	cbar.formatter.set_powerlimits((0, 0))
	#cb.format=ticker.FuncFormatter(fmt)	
	cbar.ax.xaxis.set_ticks_position('bottom')
	cbar.update_ticks()


	
def Plot3Dfields(FieldData,plane,Type,Vmin,Vmax,title,SystemParameters,ax,fig):
	#x = np.arange(SystemParameters['Nx'])*float(SystemParameters['Cell_Dx'][0])
	#y = np.arange(SystemParameters['Ny'])*float(SystemParameters['Cell_Dx'][0])
	#X, Y = np.meshgrid(x,y)
	#print(X,Y)
	sizeX = int(SystemParameters['NumCellsX'])*float(SystemParameters['Dx'])
	sizeY = int(SystemParameters['NumCellsY'])*float(SystemParameters['Dy'])
	sizeZ = int(SystemParameters['NumCellsZ'])*float(SystemParameters['Dz'])
	dampSizeX = int(SystemParameters['DampCellsX'])*float(SystemParameters['Dx'])
	cm = col.LinearSegmentedColormap('fields',cdict,N=1024,gamma=1)
#	im=ax.imshow(FieldData[Type], vmin=Vmin, vmax=Vmax, cmap=cm,origin='lower',aspect='auto',extent=[0,72,0,64])
	
	shape = FieldData[Type].shape

	if plane == "xy":
		ax.set_xlabel('$x$')
		ax.set_ylabel('$y$')
		#ax.set_xlim(0,sizeX-2*dampSizeX)
		ex = [-dampSizeX,sizeX-dampSizeX,0,sizeY]
		data = np.zeros((shape[0],shape[1]) )
		for i in range(shape[0]):
			for j in range(shape[1]):
				data[i][j] = FieldData[Type][i][j][shape[2]//2]	
	
	if plane == "xz":
		ax.set_xlabel('$x$')
		ax.set_ylabel('$z$')
		ex = [-dampSizeX,sizeX-dampSizeX,0,sizeZ]

		data = np.zeros( (shape[0],shape[2]) )

		for i in range(shape[0]):
			for j in range(shape[2]):
				data[i][j] = FieldData[Type][i][shape[1]//2][j]
	

	if plane == "yz":
		ex = [0,sizeY,0,sizeZ]
		ax.set_xlabel('$y$')
		ax.set_ylabel('$z$')	
		data = np.zeros( (shape[1],shape[2]) )
		for i in range(shape[1]):
			for j in range(shape[2]):
				data[i][j] = FieldData[Type][shape[0]//2][i][j]		
	#print(data)

	data = np.swapaxes(data,0,1)


	im=ax.imshow(data, vmin=Vmin, vmax=Vmax, cmap=cm,origin='lower',extent=ex,aspect='auto')
	#ax.xaxis.set_ticks(x)
	#ax.yaxis.set_ticks(y)
	locator=ticker.MaxNLocator(prune=None, nbins=5)
	ax.yaxis.set_major_locator(locator)
	locator=ticker.MaxNLocator(prune=None, nbins=5)
	ax.xaxis.set_major_locator(locator)
	#ax.set_title(title)
#	ax.set_ylim(-1,25)
	#ax.set_xlim(0,sizeZ-2*dampSizeZ)
	#ax.set_xlabel('$x$')
	#ax.set_ylabel('$y$')
	ax.set_title(title)
	#ax.set_xlim(int(SystemParameters['X_DampCellsLeft'][0]), int(SystemParameters['Nx'][0])-int(SystemParameters['X_DampCellsRight'][0]))
	#ax.set_xlim(int(SystemParameters['Z_DampCellsLeft'][0]), int(SystemParameters['Nz'][0])-int(SystemParameters['Z_DampCellsRight'][0]))
	#ax.set_ylim(-10, int(SystemParameters['Nr'][0])-300)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("bottom", size="8%", pad=0.6)#свойства colorbar'а
	cbar =fig.colorbar(im, cax=cax, orientation='horizontal')
	tick_locator = ticker.MaxNLocator(nbins=3)
	cbar.locator = tick_locator
		#cb = fig.colorbar(im,cax=cax, orientation='horizontal')#,ticks=[-1e-6,-0.5e-6,0,0.5e-6,1e-6]
	#cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
	#cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])  # horizontal colorbar
	#cb = ax.colorbar(im,cax=cax, orientation='horizontal')#,ticks=[-1e-6,-0.5e-6,0,0.5e-6,1e-6]
	#cb.outline.set_linewidth(0.2)#ширина border'а у colorbar'а
	cbar.formatter.set_powerlimits((0, 0))
	#cb.format=ticker.FuncFormatter(fmt)	
	cbar.ax.xaxis.set_ticks_position('bottom')
	cbar.update_ticks()


def PlotIntegMPI(IntegEnergy,SystemParameters,ax):
	time=np.arange(0,len(IntegEnergy['AllFields']))*float(SystemParameters['dt'][0])
	ax.plot(time,IntegEnergy['AllFields']+IntegEnergy['0']['IONS_Chg']+IntegEnergy['1']['ELECTRONS_Chg']+IntegEnergy['2']['LEFTBEAM_Chg'],label=SystemParameters['FolderName'][0])
	ax.set_xlabel('$t/\omega_p$')


import os
import numpy as np
import struct
import collections

def ReadParameters():
  SystemParameters = {}
  ParticlesParameters={}
  
  with open('../SysParams.cfg', 'r') as f1:
    for line in f1:
      t = line.split()
      if t[0]=="NumOfPartSpecies": SystemParameters['NumOfPartSpecies'] = t[1]
      if t[0]=="NumCellsX_glob": SystemParameters['NumCellsX'] = t[1]
      if t[0]=="NumCellsY_glob": SystemParameters['NumCellsY'] = t[1]
      if t[0]=="NumCellsZ_glob": SystemParameters['NumCellsZ'] = t[1]
      if t[0]=="DampCellsX_glob": SystemParameters['DampCellsX'] = t[1].replace(',', '')
      if t[0]=="DampCellsY_glob": SystemParameters['DampCellsY'] = t[1].replace(',', '')
      if t[0]=="DampCellsZ_glob": SystemParameters['DampCellsZ'] = t[1].replace(',', '')
      if t[0]=="Dx": SystemParameters['Dx'] = t[1]
      if t[0]=="Dy": SystemParameters['Dy'] = t[1]
      if t[0]=="Dz": SystemParameters['Dz'] = t[1]
      if t[0]=="Dt": SystemParameters['Dt'] = t[1]
      if t[0]=="TimeStepDelayDiag3D": SystemParameters['TimeDelay3D'] = t[1]
      if t[0]=="TimeStepDelayDiag2D": SystemParameters['TimeDelay2D'] = t[1]
      if t[0]=="TimeStepDelayDiag1D": SystemParameters['TimeDelay1D'] = t[1]
      if t[0]=="BUniform": SystemParameters['UniformB'] = [t[1].replace(',', '')]
  
  with open('../PartParams.cfg', 'r') as f2:
    for line in f2:
      t = line.split()
      if t[0]=="Particles" : 
        ParticlesParameters[t[1]] = {}
        sort = t[1]
      if t[0]=="Density": ParticlesParameters[sort]['Dens'] = t[1]      
      if t[0]=="Px_max": ParticlesParameters[sort]['Px_max'] = t[1]      
      if t[0]=="Px_min": ParticlesParameters[sort]['Px_min'] = t[1]  

  SystemParameters['Particles'] = ParticlesParameters
  return SystemParameters

def ReadFieldsFile(WorkDir,SystemParameters,TimeStep):
	FieldsTitles=['Ex','Ey','Ez','Bx','By','Bz']
	#открыли файл
	FieldsName=WorkDir+'Fields/Diag3D/'+'Field3D'+TimeStep
	file = open(FieldsName, 'rb')
	#прочли параметры файла
	buf = file.read(3*4)
	Prop=struct.unpack("fff", buf[:3*4])
	Nx=int(Prop[0])
	Ny=int(Prop[1])
	Nz=int(Prop[2])
	size=str(int(Nx)*int(Ny)*int(Nz))
	print(Nx,Ny,Nz,size,SystemParameters['3D_Diag_Fields'])
	#SystemParameters['Nr'][0]=Nr
	#SystemParameters['Nz'][0]=Nz
	#составить формат файла
	ftype="("+str(3)+")f4"

	for f in range(0,6):
		if(SystemParameters['3D_Diag_Fields'][f]=='1'):
			ftype+=",("+size+")f4"
	print(ftype)
	RawData = np.fromfile(FieldsName, dtype=ftype)
	FieldData=[]
	for j in range(0, 6):
		if(SystemParameters['3D_Diag_Fields'][j]=='1'):
			data=RawData[0][j+1]
			FieldData.append(data.reshape(( Nx,Ny,Nz)))
#	print(FieldData[2][10][1000])
	result=collections.OrderedDict()
	for f in range(0,6):
		if(SystemParameters['3D_Diag_Fields'][f]=='1'):
			result[FieldsTitles[f]]=FieldData[f]
			#result[FieldsTitles[f]]=np.swapaxes(result[FieldsTitles[f]],0,1)
#pprint.pprint(result)
	return result

def ReadDensFile(WorkDir,SystemParameters,PartSort,TimeStep):
	#открыли файл
	DensName=WorkDir+'Particles/'+PartSort+'/Diag3D/'+'Dens3D'+TimeStep
	#в первых двух флоатах записаны размеры сетки (на самом деле уже не актуально)
	file = open(DensName, 'rb')
	buf = file.read(3*4)
	Prop=struct.unpack("fff", buf[:3*4])
	Nx=int(Prop[0])
	Ny=int(Prop[1])
	Nz=int(Prop[2])
	size=str(int(Nx)*int(Ny)*int(Nz))
	#print('dens',Nx,Ny,size)
	ftypeDens="("+str(3)+")f4,("+size+")f4"
	RawData = np.fromfile(DensName, dtype=ftypeDens)
	data=RawData[0][1]
	Dens={}
	Dens[PartSort]=data.reshape(( Nx, Ny, Nz))
	#Dens[PartSort]=np.swapaxes(Dens[PartSort],0,1)

	return Dens

def ReadRadFile2D(WorkDir,name,SystemParameters,TimeStep):
	#открыли файл
	DensName=WorkDir+'Fields/Diag2D/'+'Rad'+name+TimeStep
	#в первых двух флоатах записаны размеры сетки (на самом деле уже не актуально)
	file = open(DensName, 'rb')
	buf = file.read(2*4)
	Prop=struct.unpack("ff", buf[:2*4])
	Nx=int(Prop[0])
	Ny=int(Prop[1])
	#Nz=int(Prop[2])
	size=str(int(Nx)*int(Ny))
	#print('dens',Nx,Ny,size)
	ftypeDens="("+str(2)+")f4,("+size+")f4"
	RawData = np.fromfile(DensName, dtype=ftypeDens)
	data=RawData[0][1]
	#Dens={}
	Dens=data.reshape(( Nx, Ny))
	#Dens[PartSort]=np.swapaxes(Dens[PartSort],0,1)

	return Dens


def ReadFieldsFile2DNew(Name,FieldsTitles,TimeStep):
	N = len(FieldsTitles)
	FieldsName=Name+TimeStep
	file = open(FieldsName, 'rb')
	buf = file.read(2*4)
	Prop=struct.unpack("ff", buf[:2*4])
	Nx=int(Prop[0])
	Ny=int(Prop[1])
	size=str(int(Nx)*int(Ny))
	print(Nx,Ny,size,FieldsName)

	ftype="("+str(2)+")f4"

	for f in range(0,N):
		ftype+=",("+size+")f4"
	print(ftype)
	RawData = np.fromfile(FieldsName, dtype=ftype)
	FieldData=[]
	for j in range(0, N):
		data=RawData[0][j+1]
		FieldData.append(data.reshape(( Nx,Ny)))
	result=collections.OrderedDict()
	for f in range(0,N):
		result[FieldsTitles[f]]=FieldData[f]
	if N > 1 :
		return result
	else:
		return FieldData[0]

def ReadFieldsFile2D(WorkDir,name,axes,SystemParameters,TimeStep):
	FieldsTitles=['Ex','Ey','Ez','Bx','By','Bz']
	#открыли файл
	FieldsName=WorkDir+'Fields/Diag2D/'+'Field'+name+TimeStep
	file = open(FieldsName, 'rb')
	#прочли параметры файла
	buf = file.read(2*4)
	Prop=struct.unpack("ff", buf[:2*4])
	Nx=int(Prop[0])
	Ny=int(Prop[1])
	size=str(int(Nx)*int(Ny))
	print(Nx,Ny,size,SystemParameters['2D_Diag_Fields'])
	#SystemParameters['Nr'][0]=Nr
	#SystemParameters['Nz'][0]=Nz
	#составить формат файла
	ftype="("+str(2)+")f4"

	for f in range(0,6):
		if(SystemParameters['2D_Diag_Fields'][f]=='1'):
			ftype+=",("+size+")f4"
	print(ftype)
	RawData = np.fromfile(FieldsName, dtype=ftype)
	FieldData=[]
	for j in range(0, 6):
		if(SystemParameters['2D_Diag_Fields'][j]=='1'):
			data=RawData[0][j+1]
			FieldData.append(data.reshape(( Nx,Ny)))
#	print(FieldData[2][10][1000])
	result=collections.OrderedDict()
	for f in range(0,6):
		if(SystemParameters['2D_Diag_Fields'][f]=='1'):
			result[FieldsTitles[f]]=FieldData[f]
			#result[FieldsTitles[f]]=np.swapaxes(result[FieldsTitles[f]],0,1)
#pprint.pprint(result)
	return result

def ReadDensFile2D(WorkDir,dataType,axes,SystemParameters,PartSort,TimeStep):
	#открыли файл
	DensName=WorkDir+'Particles/'+PartSort+'/Diag2D/'+dataType+axes+TimeStep
	#в первых двух флоатах записаны размеры сетки (на самом деле уже не актуально)
	file = open(DensName, 'rb')
	buf = file.read(2*4)
	Prop=struct.unpack("ff", buf[:2*4])
	Nx=int(Prop[0])
	Ny=int(Prop[1])
	size=str(int(Nx)*int(Ny) )
	#print('dens',Nx,Ny,size)
	ftypeDens="("+str(2)+")f4,("+size+")f4"
	RawData = np.fromfile(DensName, dtype=ftypeDens)
	data=RawData[0][1]
	Dens={}
	Dens[PartSort]=data.reshape(( Nx, Ny))
	#Dens[PartSort]=np.swapaxes(Dens[PartSort],0,1)
	if dataType == "Lambda":
		print(Dens[PartSort][10,10])
	return Dens

def ReadPhaseFile(WorkDir,SystemParameters,PartSort,TimeStep):
	#открыли файл
	DensName=WorkDir+'Particles/'+PartSort+'/Diag2D/'+'Phase2D'+TimeStep
	#в первых двух флоатах записаны размеры сетки (на самом деле уже не актуально)
	file = open(DensName, 'rb')
	buf = file.read(2*4)
	Prop=struct.unpack("ff", buf[:2*4])
	Nr=int(Prop[0])
	Nz=int(Prop[1])
	size=str(int(Nr)*int(Nz))
	#print('dens',Nx,Ny,size)
	ftypeDens="("+str(2)+")f4,("+size+")f4"
	RawData = np.fromfile(DensName, dtype=ftypeDens)
	data=RawData[0][1]
	Dens={}
	Dens[PartSort]=data.reshape(( Nr, Nz))
	Dens[PartSort]=np.swapaxes(Dens[PartSort],0,1)

	return Dens

def ReadPhaseFileOld(WorkDir,SystemParameters,PartSort,TimeStep):
	#открыли файл
	PhaseName=WorkDir+'Particles/'+PartSort+'/Diag2D/'+'Phase2D'+TimeStep
	file = open(PhaseName, 'rb')
	#прочли параметры файла
	buf = file.read(2*4)
	Prop=struct.unpack("ff", buf[:2*4])
	#Npart=int(Prop[0])#общее число частиц, если вдруг вздумаем как-то нормировать
	Npx=int(Prop[0])#число отрезков по скорости
	Ncx=int(Prop[1])#число отрезков по координате (никак не связано с расчётной сеткой)
	size=str(int(Npx)*int(Ncx))
	ftypePhase="("+str(2)+")f4,("+size+")f4"
	##("+size+")f4"

	RawData = np.fromfile(PhaseName, dtype=ftypePhase)
	data=RawData[0][1]
	Phase_x=data.reshape(( Ncx, Npx))
	#data=RawData[0][2]
	#Phase_y=data.reshape(( Npx, Ncx))
	#Phase_x=np.swapaxes(Phase_x,0,1)

	Phase={}
	Phase[PartSort]={}
	Phase[PartSort]['vx']=Phase_x
	#Phase[PartSort]['vy']=Phase_y
	return Phase

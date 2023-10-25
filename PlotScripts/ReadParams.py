#!/usr/bin/python3

ParticlesParameters={}
with open('../SysParams.cfg', 'r') as f1:
    for line in f1:
      t=line.split()
      print(t)
      if t[0]=="NumOfPartSpecies": SystemParameters['NumOfPartSpecies'] = t[1]
      if t[0]=="NumCellsZ_glob": print("afsfagvae") #SystemParameters['NumCellsZ'] = t[1]
      if t[0]=="NumCellsZ_glob": SystemParameters['NumCellsZ'] = t[1]
      if t[0]=="NumCellsR_glob": SystemParameters['NumCellsR'] = t[1]
      if t[0]=="DampCellsZ_glob": SystemParameters['DampCellsZ'] = t[1]
      if t[0]=="DampCellsR_glob": SystemParameters['DampCellsR'] = t[1]
      if t[0]=="Dz": SystemParameters['Dz'] = t[1]
      if t[0]=="Dr": SystemParameters['Dr'] = t[1]
      if t[0]=="Dt": SystemParameters['Dt'] = t[1]
      if t[0]=="TimeDelay": SystemParameters['TimeDelay'] = t[1]
      if t[0]=="TimeDelay1D": SystemParameters['TimeDelay1D'] = t[1]
      if t[0]=="BUniform": 
          SystemParameters['BUniform'] = t[1]
      
with open('../PartParams.cfg', 'r') as f2:
    for line in f2:
      t = line.split()
      if t[0]=="Particles": 
          ParticlesParameters[t[1]]={}
          sort = t[1]
      if t[0]=="Density": ParticlesParameters[sort]['Density'] = t[1]      

SystemParameters['Particles']=ParticlesParameters


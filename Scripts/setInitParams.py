
def setConst(Param,Type,Name,Value,Dim):
	Param.append({})
	Param[len(Param)-1]['Name']=Name
	Param[len(Param)-1]['Type']=Type
	Param[len(Param)-1]['Value']=Value

def setParams(PartParam, PName, PartDict):
       PartParam[PName] = {}
       for p in PartDict.keys():
               PartParam[PName][p] = PartDict[p]

    

def writeParams(ptype, filename,Param):
	f = open(filename, 'w')
	for p in Param.keys():
		tempParam = Param[p]
		f.write('####################################\n')
		f.write( ptype + ' ' + p + '\n')
		for k in tempParam.keys():
				val=str(tempParam[k])
				val=val.replace('[', '')
				val=val.replace("'", '')
				val=val.replace('"', '')
				val=val.replace(',', '')
				val=val.replace(']', '')
				f.write( k +' ' +  str(val)+'\n')

	f.close()

def writeConst(cppFile, cfgFile, Param):
	f1 = open(cppFile, 'w')
	f2 = open(cfgFile, 'w')

	for p in range(0,len(Param)):
		tempParam=Param[p]
		if(len(tempParam['Value'])!=1):
			val=str((tempParam['Value']))
			val=val.replace('[', '{')
			val=val.replace(']', '}')
			f2.write(tempParam['Name']+' '+ str((val.replace('{','')).replace('}',''))+"\n")
			tempParam['Name']=tempParam['Name']+"["+str(len(tempParam['Value']))+"]"
			f1.write(tempParam['Type']+' '+tempParam['Name'] + " = "+ str(val)+";\n")
		else:
			f1.write(tempParam['Type']+' '+tempParam['Name'] + " = "+ str(tempParam['Value'][0])+";\n")
			f2.write(tempParam['Name'] + " "+ str(tempParam['Value'][0])+"\n")

	f1.close()
	f2.close()
		
def writeDefine(filename,Param):
	f = open(filename, 'w')
	for p in range(0,len(Param)):
		f.write(Param[p]['Type']+' '+Param[p]['Name'] + " "+ str(Param[p]['Value'][0])+"\n")
	f.close()

def enum(*sequential, **named):
	enums = dict(zip(sequential, range(len(sequential))), **named)
	return type('Enum', (), enums)
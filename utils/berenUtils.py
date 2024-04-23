import sys

def setParams(PartParam, PName, PartDict):
    PartParam[PName] = {}
    for p in PartDict.keys():
        PartParam[PName][p] = PartDict[p]


def writeParams(ptype, filename, Param):
    f = open(filename, 'w')
    for p in Param.keys():
        tempParam = Param[p]
        f.write('####################################\n')
        f.write(ptype + ' ' + p + '\n')
        for k in tempParam.keys():
            val = str(tempParam[k])
            val = val.replace('[', '')
            val = val.replace("'", '')
            val = val.replace('"', '')
            val = val.replace(',', '')
            val = val.replace(']', '')
            f.write(k + ' ' + str(val)+'\n')
    f.close()


def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

def option():
    if len(sys.argv) > 1:
        return sys.argv[1]
    else:
        return ""
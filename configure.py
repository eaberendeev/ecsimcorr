from set_params import *
from utils.berenUtils import *
import os
import sys
import shutil

### Name folder which will be used for simulation
execName = "beren3d"

### Path to Eigen library
EigenPath = "~/soft/eigen-3.4.0/" #
AmgclPath = "~/soft/amgcl/" #
#EigenPath = "/home/berendeev/bpi/Progs/eigen-3.4.0/"
#AmgclPath = "/home/berendeev/bpi/Progs/amgcl/"
BuildType = "Release"
platform = "nix"

CurrentDir = os.getcwd()

print("SIMULATION WILL BE EXECUTE IN " + WorkDir)

SourceDir = CurrentDir + "/srcBeren"
PlottingDir = CurrentDir + "/PlotScripts"

#WorkDir = CurrentDir + "/" + WorkDirName  # WORK DIRECTORY

BuildDir = CurrentDir + "/_build"

if option(1) == "rebuild":
    shutil.rmtree(BuildDir, ignore_errors=True)

os.makedirs(BuildDir, exist_ok=True)

os.chdir(BuildDir)
os.system("cmake -DPATH_TO_EIGEN=" + EigenPath +
        " -DPATH_TO_AMGCL=" + AmgclPath +
          " -DCMAKE_BUILD_TYPE=" + BuildType +
          " " + SourceDir)
os.system("cmake --build . -j8")
os.system("cmake --install .")
os.chdir(CurrentDir)

if option(1) == "rerun":
    shutil.rmtree(WorkDir, ignore_errors=True)
    os.system("rm -rf " + WorkDir)

try:
    os.makedirs(WorkDir)
except:
    print("********** RUN FAILED! **********\n")
    print("Work directory exist! Please change WorkDir parameters or delete work directory before run\n")
    os.remove(CurrentDir + "/workdir.tmp")
    sys.exit()

try:
    print(BuildDir + "/bin/"+execName)
    shutil.copy(BuildDir + "/bin/"+execName, WorkDir)
except:
    print("********** COMPILATION FAILED! **********\n")
    sys.exit()

os.system("rm " + BuildDir + "/bin/"+execName)

try:
    shutil.move(CurrentDir + "/SysParams.cfg", WorkDir)
    shutil.move(CurrentDir + "/PartParams.cfg", WorkDir)
    shutil.move(CurrentDir + "/Diagnostics.cfg", WorkDir)
    shutil.move(CurrentDir + "/phys.par", WorkDir)
except:
    print("********** configs is not used! **********\n")

shutil.copytree(SourceDir, WorkDir+"/srcBeren")
shutil.copytree(PlottingDir, WorkDir+"/PlotScripts")

shutil.copy(CurrentDir + "/start.sh", WorkDir)
shutil.copy(CurrentDir + "/configure.py", WorkDir)
shutil.copy(CurrentDir + "/set_params.py", WorkDir)


f = open('workdir.tmp', 'w')
f.write(WorkDir)
f.close()
f = open('proc.tmp', 'w')
f.write(str(NumProcs))
f.close()
f = open('name.tmp', 'w')
f.write(execName)
f.close()

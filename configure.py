from set_params import *
from utils.berenUtils import *
import os
import sys
import shutil

### Name folder which will be used for simulation
CurrentSimulation = "testCircle2"
print("USING " + CurrentSimulation)

#WorkDirName = "Res_tst"
print("SIMULATION WILL BE EXECUTE IN " + WorkDir)

### Path to Eigen library
EigenPath = "~/soft/eigen-3.4.0/" #
BuildType = "Release"
platform = "nix"

CurrentDir = os.getcwd()

SourceDir = CurrentDir + "/srcBeren"
PlottingDir = CurrentDir + "/PlotScripts"

#WorkDir = CurrentDir + "/" + WorkDirName  # WORK DIRECTORY

BuildDir = CurrentDir + "/_build"

if option() == "rebuild":
    shutil.rmtree(BuildDir, ignore_errors=True)

os.makedirs(BuildDir, exist_ok=True)

os.remove(CurrentDir + "/cluster.tmp")
os.remove(CurrentDir + "/queue.tmp")
os.remove(CurrentDir + "/proc.tmp")

os.chdir(BuildDir)
os.system("cmake -DPATH_TO_EIGEN=" + EigenPath +
          " -DCMAKE_BUILD_TYPE=" + BuildType +
          " -DCURRENT_SIMULATION=" + CurrentSimulation +
          " " + SourceDir)
os.system("cmake --build . -j8")
os.system("cmake --install .")
os.chdir(CurrentDir)

if option() == "rerun":
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
    shutil.copy(BuildDir + "/bin/beren3d.exe", WorkDir)
except:
    print("********** COMPILATION FAILED! **********\n")
    sys.exit()

os.system("rm " + BuildDir + "/bin/beren3d.exe")
shutil.copytree(SourceDir, WorkDir+"/srcBeren")
shutil.copytree(PlottingDir, WorkDir+"/PlotScripts")
shutil.move(CurrentDir + "/SysParams.cfg", WorkDir)
shutil.move(CurrentDir + "/PartParams.cfg", WorkDir)
shutil.move(CurrentDir + "/Diagnostics.cfg", WorkDir)
shutil.copy(CurrentDir + "/start.sh", WorkDir)

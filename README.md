###
Plasma simulation program.
We use particle in cell model with ECSim correction method

## Dependency
Eugen library - (headers)
https://eigen.tuxfamily.org/

## ./configure.py 
script for compile and create work directory. Call set_params.py
In this file you set variant of simulation and path to Eigen library
## ./start.sh - for run configure.py and execute program
For start work set paths in configure.py and parameters in set_params.py and run ./start.sh

## ./set_params.py
Set main simulation parameters and write it to *.cfg files

## ./Scripts and ./utils 
## TO DO: merge this files
util scripts for write parameters to *.cfg files 

## ./PlotScripts
Visualization python files

## ./srcBeren
Source code
./srcBeren/constants - parameters for simulation bound. ## TO DO: move this to parameters
./srcBeren/diagnostics - sources for output diagnostic (fields, densities at all)
./srcBeren/fields - sources for fields storage and manipulation. ## TO DO: move storage to simulation class
./srcBeren/fields - sources for fields storage and manipulation. ## TO DO: move storage to simulation class
./srcBeren/particles - sources for particles storage and manipulation. ## TO DO: move storage to simulation class

17.05.2024 comment - this folder use for perdormance testing!
./srcBeren/testCircle - simulation circle injection particles to pseudo 2D cilinder - 

17.05.2024 comment - this folder use for future implement MPI
./srcBeren/testCircle2 - simulation circle injection particles to pseudo 2D cilinder

./srcBeren/testCollision - simulation test collision problem

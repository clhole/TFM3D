# Simulation and Evaluation of 3D Traction Force Microscopy:

This Software can be used to simulate and evaluate 3D traction force microscopy (TFM) data. 
The integrated approach simulates data in ANSYS and then evaluates the generated data with a 
range of different algorithms in MATLAB. Nearly every parameter can be changed at will. This will 
help to find guidelines for optimal data collection and evaluation. The software is written in a 
modular way, thus it is also possible to analyze experimental data, or to combine a variety of algorithms 
to calculate the displacement and traction. 

A deeper description of software and background will hopefully soon be published (and linked here).

# Software Requirements:
You need MATLAB as well as ANSYS (if you want to simulate data, and also one of the 
traction reconstruction algorithms needs ANSYS). The code was written on Windows, and 
not tested on Apple or Linux machines.

# File Organization:
The main file is TFM3D.m, which also contains documentation about how to use it. 
The subfolders then contain a variety of functions used for the different modules. 
There are also some external functions from the MATLAB file exchange, as well as 
from the Franck Lab: https://github.com/FranckLab/FIDVC. All credit for the FIDVC 
algorithm goes to them. We copied the files to this GitHub page to make sure that our code 
implementation is compatible also if they choose to update their code.

# Additional Files needed:
Some of the files were too large for GitHub, 
you'll find them here: https://doi.org/10.3929/ethz-b-000289734. (at some point, in the mean time here:https://www.research-collection.ethz.ch/handle/20.500.11850/289734)

# How to start:
If you run the code for the first time, you will need to 
generate some data first with ANSYS, after that you can also only use some submodules.

# Questions:
If you have any questions or recommendations, 
please don't hesitate to contact us: http://www.biomech.ethz.ch/research/jess-snedeker.html

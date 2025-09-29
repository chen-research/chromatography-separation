This repository stores the codes for the paper of Optimal Parameter Design of High Resolution for Preparative Chromatography by You Sun, Chen Xu, and Ye Zhang.

Python Version: 3.8.8
TensorFlow Version: 2.8.0
SciPy Version: 1.6.2

Train_FNN_GsPowerOpt.py 
is the code for performing any one of Experiment 1-8 (transport-dispersive model with Langmuir isotherm) in the paper. The simulated data, which is generated in another code file and should be imported to Train_FNN_GsPowerOpt.py, determines the index of experiment to be performed.

predefined_functions.py 
defines the classes and functions to be called by Train_FNN_GsPowerOpt.py.


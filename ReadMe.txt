This repository stores the codes for the paper of Optimal Parameter Design of High Resolution for Preparative Chromatography by You Sun, Chen Xu, and Ye Zhang.

Python Version: 3.8.8
TensorFlow Version: 2.8.0
SciPy Version: 1.6.2

Train_FNN_GsPowerOpt.py 
is the code for performing any one of Experiments 1-8 (transport-dispersive model with Langmuir isotherm) in the paper. The simulated data, which is generated in another code file and should be imported to Train_FNN_GsPowerOpt.py, determines the index of experiment to be performed.

predefined_functions.py 
defines the classes and functions to be called by Train_FNN_GsPowerOpt.py.


MATLAB: R2023a (9.14)
- OS: Windows 11 Pro (x64)
- CPU: 24 physical cores / 32 logical cores  
- Toolboxes: Parallel Computing Toolbox

Files in dataGenerate
are the codes for generating the dataset (x, r), used to train the Feedforward Neural Network (FNN) in Experiments 1–8. The code simulates outlet concentrations of a column by solving a partial differential equation (PDE) using the finite volume approximation method, based on a large set of randomly generated x-samples. These concentrations are then used to construct the response curve r.

Files in generate_response
are the code for performing Experiments 1-8 in the paper. The code is used for generating and plotting chromatographic outlet response with known parameters.


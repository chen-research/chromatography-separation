############# - Load packages
from sklearn.model_selection import train_test_split
from predefined_functions import *
import random
import tensorflow as tf
import pandas as pd
import numpy as np
from sklearn.metrics import r2_score
import os
#from sklearn.neural_network import MLPRegressor
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import shift
print(tf.__version__)
#import datetime
import gc
from sklearn.metrics import mean_squared_error as mse
from scipy.signal import find_peaks, peak_widths
#To limit TensorFlow to CPU
cpus = tf.config.experimental.list_physical_devices('CPU') 
#tf.config.experimental.set_visible_devices(gpus[0], 'GPU')
tf.config.experimental.set_visible_devices(cpus[0], 'CPU')
data_path = "K:/data/" #Replace with path to your simulated data of (x,r)


############### - Load the simulated data and preprocess for training the FNN
density = np.array(pd.read_csv(data_path+"Out6-time.csv")) #r
param_inj = np.array(pd.read_csv(data_path+"params6-time.csv"))  #[y, h]

#scale the density by 1,000,000
density = density*1000000

#Set the sample size
sample_size = density.shape[0]

#Check for duplicates
print('param shape before dropping duplicates',param_inj.shape, 
      'after duplicates are dropped', pd.DataFrame(param_inj).drop_duplicates().shape)
print('density shape before dropping duplicates',density.shape, 
      'after duplicates are dropped', pd.DataFrame(density).drop_duplicates().shape)

#Set the smaple size for training, validation, and testing
tr_size = int(sample_size*0.8)
val_size = int(sample_size*0.1)
te_size = int(sample_size*0.1)
train_val_size = val_size+tr_size

#Get the train, validation, and test data
train_x, test_x, train_y, test_y = train_test_split(param_inj, 
                                                    density, 
                                                    test_size=te_size/sample_size 
                                                   )

train_x, val_x, train_y, val_y = train_test_split(train_x, 
                                                   train_y, 
                                                   test_size=val_size/train_val_size 
                                                   )

train_m = np.mean(train_x,axis=0)
train_sd = np.std(train_x,axis=0)


############### - Plot a simulated response r - for Figure 3/6.
i = 210
plt.figure()
plt.plot(train_y[i])
plt.show()
print(train_x[i])

#Build fw-FNN (with the optimal hyper-parameters) for training bw-FNN
#Only the training dataset is used for training this fw-FNN


###################--Build and train FNN F
fw_FNN = build_fwmodel(loss_func = 'L2', #Loss function to be used, "L1L2":error is L1 norm, bias and weights are L2 norms
                       activ='tanh', #The activation function
                       hid_layers= (20,100,200), #(50,100,200,400,800), #(1st hidden layer nodes, 2nd hidden layer nodes)
                       bias_regu_cosnt= 0.01, #The regularization coeff. for bias terms
                       w_regu_const= 0.01, #The regularization coeff. for weights
                       output_size = 800
         )

early_stop = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=20) 
history = fw_FNN.fit((train_x-train_m)/train_sd, 
                     train_y, 
                     epochs=2000, 
                     batch_size=1000, 
                     shuffle=True, 
                     verbose=2, 
                     callbacks=[early_stop],
                     validation_data=((val_x-train_m)/train_sd, val_y)
                    )

#fw-FNN Train and Test Performance - Test_R2 is for Table 1
train_pred = fw_FNN.predict((train_x-train_m)/train_sd)
train_r2 = compute_r2(train_y,train_pred,low_percentile_to_remove=0.5)     
test_pred = fw_FNN.predict((test_x-train_m)/train_sd)
test_r2 = compute_r2(test_y,test_pred,low_percentile_to_remove=0.5)                                     
print('fw-FNN:', 
      "Train_R2:",round(train_r2,3),
      "Test_R2:",round(test_r2,3),
      "Train_MSE:", round(history.history['MSE'][-1],3),
     ) 


#################### GS-PowerOpt
def objective(inputs):
    """
    Inputs.
    --------
    inputs: np.2darray, each row of params is a vector of:  
                        the adsorption parameters and 2 projection profile parameters
                        
    Outputs.
    --------
    objective_to_maximize, np.1darray, objective_to_maximize[n] is the fitness of params[n]
    """
    
    #rescale the inputs to the desired range
    params = inputs.copy()
    params[:,0:-2] = 10 + 40*(1+np.sin(inputs[:,0:-2])) #adsorption params
    params[:,-2:] = 5 + 10*(1+np.sin(inputs[:,-2:])) #injection profile
    
    density = fw_FNN.predict((params-train_m)/train_sd,verbose=0)
    batch_size = density.shape[0]
    threshold = 0.1
    objective_to_maximize = np.ones(batch_size, dtype=float)*(-10)

    for row in range(batch_size):
        response = density[row]
        peaks, _ = find_peaks(response, height = 10)  #indices of all peaks
                                                      #If the height of a peak is less than 10, it will be ignored.
        if len(peaks)<2:
            continue
        
        # sort peaks to get the two highest ones
        peak_heights = response[peaks].copy()  #peak heights
        sorted_peak_indices = np.argsort(peak_heights)
        two_tops = peaks[sorted_peak_indices[-2:]]
        t_R1 = np.min(two_tops)
        t_R2 = np.max(two_tops)
        
        # get the two widths
        widths, _, _, _ = peak_widths(response, two_tops, rel_height=0.99)
        w1, w2 = widths
        
        # compute Rs
        Rs = 2*(t_R2-t_R1)/(w1+w2)
        objective_to_maximize[row] = Rs
    return objective_to_maximize


# Hyper-parameters for GS-PowerOpt
dim_num = param_inj.shape[1]
learning_rates = [0.01]
smooth_param_candidates = [0.1]
generation_num = 2000
pop_size = 100
final_fit = 0.0
final_params = {}
initial_guess = np.random.uniform(low=-0.1, high=0.1, size=dim_num)

# Perform the experiment under each value of the hyper-parameters, and select the optimal one
for N in [1,5,8]:
    for lr in learning_rates:
        for sp in smooth_param_candidates:
            print("N=",N, "lr=",lr, "sp=",sp)
            mu_times = [] #number of iterations taken to achieve the best mu
            mu_fit = []
            mu_list = []
            
            mu_log = power_gs(
                            init_mu=initial_guess, 
                            dim=dim_num, 
                            fitness_fcn=objective,
                            init_lr=lr, 
                            sigma=sp,
                            power=N,
                            sga_sample_size=pop_size,
                            total_step_limit=generation_num
            )
            #print("N:",N, "init_lr:", lr, "smoothing_param:",sp)
            best_index, best_sol, best_fit = evaluate_gs(np.array(mu_log), objective, verbose=False)
            mu_times.append(best_index) #number of iterations taken to achieve the best mu
            mu_fit.append(best_fit)
            mu_list.append(best_sol)
                
            #print(" ")
            if best_fit>final_fit:
                final_solution = best_sol.copy()
                final_index = best_index
                final_fit = best_fit
                final_params["N"] = N
                final_params["lr"] = lr
                final_params["smooth_param"] = sp
                print("best fit so far:", round(final_fit,3))
                print(final_params)
                
                Final_param = best_sol.copy()
                Final_param[0:-2] = 10 + 40*(1+np.sin(best_sol[0:-2])) #adsorption params
                Final_param[-2:] = 5 + 10*(1+np.sin(best_sol[-2:])) #injection profile
                print("best sol", Final_param)
                
                plt.figure()
                pred_den = fw_FNN.predict( ((Final_param-train_m)/train_sd).reshape([1,dim_num]) )
                plt.plot(pred_den[0])
                plt.show()

# Print mu* (Final Solution) - For Table 2/3
print(" ")
print("Final Solution:", Final_param.round(3))
print("Fitness:", round(final_fit,3))
print("Time to best:", final_index)
print(final_params)


#Plot the final solution
plt.figure(figsize=(6,5))
pred_den = fw_FNN.predict( ((Final_param-train_m)/train_sd).reshape([1,dim_num]) )/1000000
plt.xlabel("time",size=15)
plt.tick_params(labelsize=15)
plt.plot(pred_den[0])
plt.savefig("./Graphs/FF-output-pgs-time-dataGenerate4.eps",bbox_inches="tight")
plt.show()
Final_param
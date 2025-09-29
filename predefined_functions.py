import pandas as pd
import numpy as np
from tensorflow import keras
import matplotlib.pyplot as plt
import tensorflow as tf
from sklearn.metrics import r2_score

#Define the function for training a neural net model
def build_fwmodel(loss_func = 'L1', #Loss function to be used, "L1L2":error is L1 norm, bias and weights are L2 norms
                activ = 'tanh', #The activation function
                hid_layers = (25,36), #(1st hidden layer nodes, 2nd hidden layer nodes)
                bias_regu_cosnt = 0.01, #The regularization coeff. for bias terms
                w_regu_const = 0.01, #The regularization coeff. for weights
                show_process = 0, #If equals 1, show the training process, if equals 0, do not show. 
                output_size = 800
                 ):
    """
    This function returns an un-trained fw-FNN.
    """   
    #Build the model structure
    model = keras.Sequential()
    if loss_func == 'L1': 
        error_loss = 'MAE'
        for node_num in hid_layers:
            model.add(keras.layers.Dense(node_num, activation=activ, 
                                     kernel_initializer='glorot_uniform', 
                                     bias_initializer='glorot_uniform',
                                     kernel_regularizer=tf.keras.regularizers.l1(w_regu_const),
                                     bias_regularizer=tf.keras.regularizers.l1(bias_regu_cosnt))) 
    elif loss_func == 'L2':
        error_loss = 'MSE'
        for node_num in hid_layers:
            model.add(keras.layers.Dense(node_num, activation=activ, 
                                     kernel_initializer='glorot_uniform', 
                                     bias_initializer='glorot_uniform',
                                     kernel_regularizer=tf.keras.regularizers.l2(w_regu_const),
                                     bias_regularizer=tf.keras.regularizers.l2(bias_regu_cosnt)))
             
    model.add(keras.layers.Dense(output_size,activation='relu'))   
    model.compile(optimizer='adam', loss=error_loss, metrics=['MSE'])
    return model


def compute_r2(y_true,y_pred,low_percentile_to_remove=1):
    """
    Compute the average of r2 of each paired rows in (y_true,y_pred)
    """
    r2_scores = np.array([r2_score(y_true[i], y_pred[i]) for i in range(len(y_true))])
    threshold = np.percentile(r2_scores, low_percentile_to_remove)  # 1st percentile
    r2_scores_filtered = r2_scores[r2_scores > threshold]
    average_r2 = r2_scores_filtered.mean()
    return average_r2



def power_gs(init_mu, dim, fitness_fcn, 
             init_lr=0.001, sigma=1.0, power=2,
             sga_sample_size=1000,
             total_step_limit=10000, 
             verbose=False):
    """
    Powered Gaussian Smooth with a Baseline Algorithm, for optimization.
    
    Inputs.
    ---------
    init_generation:np.2darray, each row is a candidate solution.
    fitness_fcn, the fitness function.
    elite_ratio:float, the portion of best candidates in each generation used to produce a Gaussian distribution.
    update_weight:float, mu = (1-update_weight)*mu + update_weight*updated_mu
    generation_num:number of generations to be generated.
    verbose:binary, whether to show the training process.
    
    Outputs.
    ---------
    best_solution:np.1darray, the solution with the largest fitness value among all the generations.
    best_fitness:float, fitness of the best solution.
    """
    mu = init_mu
    #logs
    mu_log = []
    for k in range(total_step_limit):
        generation = mu+np.random.normal(loc=0, scale=sigma, size=(sga_sample_size,dim))
        sol = np.append(generation, mu.reshape((1,-1)), axis=0)
        stacked_fvalues = fitness_fcn(sol)
        fitness = stacked_fvalues[0:-1]
        mu_fit = stacked_fvalues[-1]
        ############ - Update mu
        alpha = init_lr #learning rate
        #v = (fitness/mu_fit)**power
        v = (fitness)**power
        gradient = np.mean( (generation-mu)*v.reshape((-1,1)), axis=0 ) #E[(X-mu)*f(X)]
        gradient = gradient/(np.sqrt(np.sum(gradient**2))) #normalize the gradient
        mu = mu + alpha*gradient
        #Print to Screen to check training statistics
        if k%100==0:
            print("k:",k)
            print("mu_fit:",mu_fit)
            print("mu_norm:", np.sqrt(np.sum(mu**2)))
        mu_log.append(mu)   
        
    return mu_log

def evaluate_gs(solution_log, f, verbose=False):
    """
    Output the best solution in solution_log.
    
    Inputs.
    --------
    solution_log:list of solutions, each row in solution_log is a solution.
    f: the fitness function.

    Outputs.
    --------
    best_sol:1darray, the best solution in solution_log.
    best_fit:float, the fitness value of the best solution.
    
    """
    fitness = f(solution_log)
    best_fit_index = np.argmax(fitness)
    best_sol = solution_log[best_fit_index]
    best_fit = f(best_sol.reshape([1,-1]))[0]
    if verbose:
        print("best fitness:",round(best_fit,3))
        print("Steps taken to reach the best fit:",best_fit_index)
        print("best solution:",best_sol.round(3))
    
    return [best_fit_index, best_sol, best_fit]
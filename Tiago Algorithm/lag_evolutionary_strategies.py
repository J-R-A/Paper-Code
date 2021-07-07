#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 03:14:56 2021

@author: tjpc
"""
import numpy as np

from evo_strategies import evolutionary_strategies
from utils_lagged_regression import get_data_matrix

##################
# TODO: pass params as dicts and pass all params
##################
def cost_function_lag(y,eval_params,fixed_params):
    normalize_lag = 100
    eval_params = np.hstack([np.zeros(fixed_params['padding']), eval_params, np.zeros(fixed_params['padding'])])
    eval_params = eval_params*normalize_lag
   
    _, X, time = get_data_matrix(fixed_params['spiketimes'],fixed_params['neurons'],
                                     fixed_params['trials'],fixed_params['tmin'],
                                     fixed_params['tmax'], eval_params, 
                                     fixed_params['sigma'], fixed_params['binsize'],  
                                     type_kernel = 'real')
    
    
    X_1 = X[ fixed_params['padding'] : - fixed_params['padding'] ,:,:]
    X_reg = X_1.reshape(len(X_1[:,0,0])*len(X_1[0,:,0]),len(X_1[0,0,:]),order='F')
    
    pred_y = np.matmul(X_reg,fixed_params['regressors']) #change this line to include 
    error_est = ((y-pred_y)**2).mean()
    error_reg_smooth = fixed_params['reg_smoothess']*np.mean(np.diff(np.diff(eval_params))**2)
    
    # print(np.abs(eval_params).mean())
    # print(error_est)  
    # print(error_reg_smooth)
    # print('')
    
    return error_est + error_reg_smooth


def compute_optimal_lag_ES(y, lag_vector, regressors, spiketimes, neurons, trials,
                           tmin, tmax, kern_off_bins, sigma, binsize,lrate =0.5, std_perturbations= 0.5,
                           n_perturbations = 20, max_iterations= 50, anneling_lrate = 1, 
                           anneling_std = 1, reg_smoothess = 5, minimize = True):
    
    dict_fparams = {'spiketimes'      : spiketimes,
                   'neurons'         : neurons,
                   'trials'          : trials,
                   'tmin'            : tmin,
                   'tmax'            : tmax,
                   'padding'         : kern_off_bins,
                   'regressors'      : regressors,
                   'sigma'           : sigma,
                   'binsize'         : binsize,
                   'reg_smoothess'   : reg_smoothess
                   }
    
    
    # X_reg, X, time = get_data_matrix(dict_fparams['spiketimes'],dict_fparams['neurons'],
    #                               dict_fparams['trials'],dict_fparams['tmin'],
    #                               dict_fparams['tmax'], lag_vector, 
    #                               dict_fparams['sigma'], dict_fparams['binsize'],  
    #                               type_kernel = 'real')
    

    
    results = evolutionary_strategies(y, cost_function_lag, lag_vector, lrate, std_perturbations,
                                      n_perturbations, max_iterations, anneling_lrate= anneling_lrate, 
                                      anneling_std = anneling_std, minimize = minimize, fixed_params = dict_fparams)
    
    return results[2]




def cost_function_regressor(y,eval_params,fixed_params):
    eval_params = eval_params
    
    pred_y = np.matmul(fixed_params['X_reg'],eval_params) #change this line to include 
    error_est = ((y-pred_y)**2).mean()
    # error_reg_smooth = fixed_params['lambda_r']*np.sum(eval_params)**2)
    
    # print('error inside cost function',error_est)
    # print('error', error_est)  
    # print(error_reg_smooth)
    # print('')
    # print('params',eval_params[:3])
    return error_est #+ error_reg_smooth



def compute_optimal_regressors_ES(y, lag_vector, regressors, spiketimes, neurons, trials,
                                  tmin, tmax, sigma, binsize, lrate =0.5, std_perturbations= 0.5,
                                  n_perturbations = 20, max_iterations= 50, anneling_lrate = 1, 
                                  anneling_std = 1,lambda_r = 0, minimize = True):
    
    X_reg, X, time = get_data_matrix(spiketimes,neurons,
                                  trials,tmin,
                                  tmax, lag_vector, 
                                  sigma, binsize,  
                                  type_kernel = 'real')
    
    dict_fparams = {'X_reg'      : X_reg,
                    'lambda_r'        : lambda_r,
                    }

    
    
    pred_y = np.matmul(X_reg,regressors)
    print('Initial error:', ((y-pred_y)**2).mean())    

    results = evolutionary_strategies(y, cost_function_regressor, regressors, lrate, std_perturbations,
                                      n_perturbations, max_iterations, anneling_lrate= anneling_lrate, 
                                      anneling_std = anneling_std, minimize = minimize, fixed_params = dict_fparams,
                                      mirror_perturbations = True)
    
    return results[2],results[0]



def cost_function_lag_regressors(y,eval_params,fixed_params):
    normalize_lag = 100

    lag_vector_initial = eval_params[0][:fixed_params['size_lag']]
    regressors_initial = eval_params[0][fixed_params['size_lag']:] #code this properly
    perturbation_lag = eval_params[1][:fixed_params['size_lag']]
    perturbation_regressors = eval_params[1][fixed_params['size_lag']:]
    
    lag_vector = lag_vector_initial + perturbation_lag
    regressors = regressors_initial + perturbation_regressors
    
    lag_vector = lag_vector*normalize_lag
    regressors = regressors
    
    X_reg, X, time = get_data_matrix(fixed_params['spiketimes'],fixed_params['neurons'],
                                     fixed_params['trials'],fixed_params['tmin'],
                                     fixed_params['tmax'], lag_vector, 
                                     fixed_params['sigma'], fixed_params['binsize'],  
                                     type_kernel = 'real')
    
    
    pred_y = np.matmul(X_reg,regressors) #change this line to include 
    
    error_est = ((y-pred_y)**2).mean()
    error_reg_regressors = fixed_params['lambda_r']*np.mean(regressors**2)
    error_reg_smooth = fixed_params['reg_smoothess']*np.mean(np.diff(np.diff(lag_vector))**2)
    # print(np.abs(regressors_initial).mean())

    # print(np.abs(regressors).mean())
    # print('est error: ',error_est)    
    # print('reg regressors error: ',error_reg_regressors)
    # print('reg smooth error: ',error_reg_smooth)
    # print('')
    # print(fixed_params['size_lag'])
    return error_est + error_reg_regressors + error_reg_smooth




def compute_optimal_lag_regressors_ES(y, lag_vector, regressors, spiketimes, neurons, 
                                      trials, tmin, tmax, sigma, binsize, lrate =0.5, 
                                      std_perturbations= 0.5, n_perturbations = 20, 
                                      max_iterations= 50, anneling_lrate = 1, 
                                      anneling_std = 1, lambda_r = 0, reg_smoothess = 0, 
                                      minimize = True):
    
    dict_fparams = {'spiketimes'     : spiketimes,
                    'neurons'         : neurons,
                    'trials'          : trials,
                    'tmin'            : tmin,
                    'tmax'            : tmax,
                    'lambda_r'        : lambda_r,
                    'sigma'           : sigma,
                    'binsize'         : binsize,
                    'reg_smoothess'   : reg_smoothess,
                    'size_lag'        : len(lag_vector),
                    'size_regressors' : len(regressors)
                    }
    
    
    # X_reg, X, time = get_data_matrix(dict_fparams['spiketimes'],dict_fparams['neurons'],
    #                                   dict_fparams['trials'],dict_fparams['tmin'],
    #                                   dict_fparams['tmax'], lag_vector, 
    #                                   dict_fparams['sigma'], dict_fparams['binsize'],  
    #                                   type_kernel = 'real')
  # print('X mean',X.mean())
    # print('vel mean', y.mean())
    # pred_y = np.matmul(X_reg,regressors)
    # error_est = ((y-pred_y)**2).sum()
    # error_reg_lag = dict_fparams['lambda_r']*np.sum(regressors**2)
    # error_reg_smooth = dict_fparams['reg_smoothess']*np.sum(lag_vector**2)
    # print(np.abs(regressors).mean())
    # # print(error_est, ((y-pred_y)**2).std())    
    # # print(error_reg_lag)
    # # print(error_reg_smooth, fixed_params['reg_smoothess']*np.std(eval_params**2))
    # # print('')
    # error_est + error_reg_lag + error_reg_smooth
    # print(error_est)
    # print(anneling_lrate)
    # print(anneling_std)
    # print(minimize)
    
    optim_params = np.concatenate((lag_vector,regressors))
    results = evolutionary_strategies(y, cost_function_lag_regressors, optim_params, lrate, std_perturbations,
                                      n_perturbations, max_iterations,  anneling_lrate= anneling_lrate, 
                                      anneling_std = anneling_std, minimize = minimize , fixed_params = dict_fparams, pass_perturbations = True)
    

    
    return results[2]
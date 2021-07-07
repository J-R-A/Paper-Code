#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 02:10:47 2021

@author: tjpc
"""
import numpy as np
 

from glmnet import ElasticNet

# import scipy 
# MAD = scipy.stats.median_absolute_deviation 
#custom imports
from utils_lagged_regression import get_data_matrix
from utils_lagged_regression import compute_optimal_lag_approx, apply_approx_latency
from lag_evolutionary_strategies import compute_optimal_lag_ES




def approximated_lagged_regression(speed_window, spiketimes, neurons, trials,binsize, tmin, 
                                   tmax,  niter = 5, sigma = 150, 
                                   lambda_b = 0.0, lambda_l = 0.5, alpha = 0.5, betas = []):
    
    # TODO: Change nomenculature betas -> regressors
    # TODO: implement stopping condition based on relative error and abs error
    # TODO: Change latency to lag (no idea why i choose latency in the first place...)


    ####################################################################
    # Algorithm parameters
    # sigma -> size gaussian kernel
    # lambda_b -> regularization on regressors
    # lambda_l -> regularization on lag function
    # niter -> niters algo
    # alpha -> reg parameter measures passage from rigdge to lasso reg
    #          0 -> ridge, 1 -> lasso
    #################################################################
    nbins = int((tmax - tmin) / binsize)
    
    betas_vectors = np.zeros((niter,len(np.unique(neurons))))
    latency_vectors =  np.zeros((niter,nbins+1))
    latency_vector = np.zeros(nbins+1)

    regression_model = ElasticNet(alpha = alpha, lambda_path=np.array([0]) ,fit_intercept = False, verbose = False) # check this part
    
    errs = np.zeros(2*niter)    
    X_reg, X, dX, time = get_data_matrix(spiketimes,neurons,
                                        trials,tmin,tmax,latency_vector,sigma, binsize)
    print('x_reg', X_reg.shape)
    print('x', X.shape)

    y = speed_window.reshape((nbins+1)*len(np.unique(trials)), order = 'F')
    use_reg = False
    factor_normalize_error = len(y)
    if len(betas)== 0:
        use_reg = True
    for i in range(niter):
        print('iteration num:',str(i+1))
        if use_reg:
        # Computing regressors
            fit = regression_model.fit(X_reg.copy(), y.copy())
            betas = fit.coef_
            pred_y = regression_model.predict(X_reg)
        else:
            print(X_reg.shape)
            pred_y = np.matmul(X_reg, betas)

        error1 = ((y-pred_y)**2).sum()+lambda_b*np.sum(betas**2)+lambda_l*np.sum(latency_vector**2)
        
        # Computing lag function
        latency_vector = compute_optimal_lag_approx(X, dX, speed_window, betas,0, lambda_l)
        X_reg= apply_approx_latency(X, dX,latency_vector)
        # pred_y = regression_model.predict(X_reg)    
        pred_y = np.matmul(X_reg, betas)

        error2 = ((y-pred_y)**2).sum()+lambda_b*np.sum(betas**2)+lambda_l*np.sum(latency_vector**2)
        
        # saving results
        errs[2*i]= error1/factor_normalize_error
        errs[2*i+1]= error2/factor_normalize_error
        betas_vectors[i,:] = betas
        latency_vectors[i,:] = latency_vector
    return betas_vectors, latency_vectors, errs



def compute_regressors(speed_window, lag_vector, spiketimes, neurons, trials,binsize, tmin, 
                     tmax, sigma = 150):
    ####################################################################
    # Algorithm parameters
    # sigma -> size gaussian kernel
    # lambda_b -> regularization on regressors
    # lambda_l -> regularization on lag function
    # niter -> niters algo
    # alpha -> reg parameter measures passage from rigdge to lasso reg
    #          0 -> ridge, 1 -> lasso
    #################################################################
    nbins = int((tmax - tmin) / binsize)

        
    regression_model = ElasticNet(alpha = 0, lambda_path=np.array([0]) ,fit_intercept = False, verbose = False) # check this part
    
    X_reg, X, time = get_data_matrix(spiketimes,neurons,
                                     trials,tmin,tmax,lag_vector,sigma, 
                                     binsize, type_kernel = 'real')

    
    y = speed_window.reshape((nbins+1)*len(np.unique(trials)), order = 'F')
    # Computing regressors
    fit = regression_model.fit(X_reg.copy(), y.copy())
    betas = fit.coef_
    pred_y = regression_model.predict(X_reg)

    return betas

def generate_dataset(speed_window, spiketimes, neurons, trials,binsize, tmin, 
                     tmax, sigma = 150, lag_magnitute = 150, nolag = False):
    

    ####################################################################
    # Algorithm parameters
    # sigma -> size gaussian kernel
    # lambda_b -> regularization on regressors
    # lambda_l -> regularization on lag function
    # niter -> niters algo
    # alpha -> reg parameter measures passage from rigdge to lasso reg
    #          0 -> ridge, 1 -> lasso
    ###################################################################
    nbins = int((tmax - tmin) / binsize)
    latency_vector  = generate_latency_vector(nbins, lag_magnitute)

    if nolag:
        latency_vector  = latency_vector*0
        
    regression_model = ElasticNet(alpha = 0, lambda_path=np.array([0]) ,fit_intercept = False, verbose = False) # check this part
    
    X_reg, X, time = get_data_matrix(spiketimes,neurons,
                                     trials,tmin,tmax,latency_vector,sigma, 
                                     binsize, type_kernel = 'real')

    
    y = speed_window.reshape((nbins+1)*len(np.unique(trials)), order = 'F')
    # Computing regressors
    # print(X_reg.shape,y.shape)
    fit = regression_model.fit(X_reg.copy(), y.copy())
    betas = fit.coef_
    pred_y = regression_model.predict(X_reg)

    return betas, latency_vector, pred_y, pred_y.reshape((nbins+1),len(np.unique(trials)), order = 'F'), X_reg, X, time

def generate_latency_vector(nbins, lag_magnitute):
    
    t = np.linspace(0,1,nbins+1)
    
    # y = np.exp(-10*(0.5-t)**2)
    # return (np.sin(19*t)+np.sin(11*t)+ np.sin(2*np.pi*t)-0.5*t)*y*5

    y = np.exp(-10*(0.4-t)**2)
    return (np.sin(19*t)+np.sin(11*t)+ np.sin(2*np.pi*t)-t+0.5)*y*lag_magnitute



if __name__ == '__main__':
    import matplotlib.pyplot as plt
    n = 96
    plt.plot(generate_latency_vector(n))

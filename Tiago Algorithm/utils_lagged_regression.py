#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 21:19:57 2020

@author: tjpc
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy as spy
import numpy.random as npr


def gaussian_filter_shifted(x,t_center,t_shift, sigma):
    x = x-t_center
    n_const = 1/(sigma*(2*np.pi)**0.5)
    a_exp =-1/2*((x-t_shift)/sigma)**2
    return n_const*np.exp(a_exp)

    
def gaussian_derivative_filter_shifted(x,t_center,t_shift, sigma):
    x = x-t_center
    n_const = -(x-t_shift)/(sigma**3*(2*np.pi)**0.5)
    a_exp =-1/2*((x-t_shift)/sigma)**2
    return n_const*np.exp(a_exp)


def gaussian_filter_shifted_approx(x,t_center,t_shift, sigma):
    gaussian_k = gaussian_filter_shifted(x,t_center,0, sigma)
    dgaussian_k = gaussian_derivative_filter_shifted(x,t_center,0, sigma)
    
    return gaussian_k + t_shift*dgaussian_k


def convolution_spikes_variable_kernel(x, spikes, latency_vector, 
                                       sigma, type_kernel = 'approximate'):
    y = np.zeros_like(x)
    if len(spikes)!= 0:
        latency_spikes = np.interp(spikes,x,latency_vector)
        
        for spike,shift in zip(spikes, latency_spikes):
            if type_kernel == 'approximate':
                y += gaussian_filter_shifted_approx(x,spike,shift, sigma)
            else:
                y += gaussian_filter_shifted(x,spike,shift, sigma)
    return y

def convolution_spikes_gaussian(x, spikes, latency_vector, 
                                       sigma, type_kernel = 'approximate'):
    y = np.zeros_like(x)
    if len(spikes)!= 0:
        for spike in spikes:
                y += gaussian_filter_shifted(x,spike,0, sigma)
    return y

def convolution_spikes_gaussian_derivative(x, spikes, latency_vector, 
                                       sigma, type_kernel = 'approximate'):
    y = np.zeros_like(x)
    if len(spikes)!= 0:
        for spike in spikes:
                y += gaussian_derivative_filter_shifted(x,spike,0, sigma)
    return y


def compute_optimal_lag_approx(X, dX, y, weights,a0, lambdau):
    
    #This is incorrect i need to do this computation for every time step
    
    r = np.einsum('kji,i->kj', X, weights)+a0
    dr =  np.einsum('kji,i->kj', dX, weights)
    em = (y-r)
    return -np.sum(em*dr,1)/(np.sum(dr**2,1)+lambdau)

    
def get_data_matrix(sorted_spiketimes_windowed,sorted_neurons_id_windowed,
                    sorted_trial_id_windowed, TMIN, TMAX,latency_vector,
                    sigma = 100, BINSIZE = 10, type_kernel = 'approximate'):
    
    # sigma = sigma #/ BINSIZE
    
    
    n_neurons = len(np.unique(sorted_neurons_id_windowed))
    n_trials  = len(np.unique(sorted_trial_id_windowed))
    time = np.linspace(TMIN,TMAX,int(TMAX/BINSIZE)+1)
    print(time.shape)
    # print(n_neurons)
    if type_kernel != 'approximate':
        X = np.zeros((len(time),n_trials,n_neurons))
    
        for trial in np.unique(sorted_trial_id_windowed):
            k = 0
            for neuron in np.unique(sorted_neurons_id_windowed):
                cond = (sorted_neurons_id_windowed==neuron) & (sorted_trial_id_windowed == trial)
                # print(neuron,trial, cond.sum() )
                # if cond.sum() > 50:
                spiketimes = sorted_spiketimes_windowed[cond]
                C = convolution_spikes_variable_kernel(time, spiketimes, latency_vector, 
                                                        sigma, 
                                                        type_kernel = 'real')
                
                X[:, trial, k] = C
                k+=1
        X = X*1000 # convert from spikes per ms to spikes per s
        data = X.reshape(len(time)*n_trials,n_neurons,order='F')
        return data, X, time
    
    else:
        X = np.zeros((len(time),n_trials,n_neurons))
        dX = np.zeros((len(time),n_trials,n_neurons))
        # print(X.shape)
        for trial in np.unique(sorted_trial_id_windowed):
            k = 0
            for neuron in np.unique(sorted_neurons_id_windowed):
                cond = (sorted_neurons_id_windowed==neuron) & (sorted_trial_id_windowed == trial)
                # print(neuron,trial, cond.sum() )
                # if cond.sum() > 50:
                spiketimes = sorted_spiketimes_windowed[cond].copy()
                
                C = convolution_spikes_gaussian(time, spiketimes, latency_vector, 
                                                        sigma)
                
                X[:, trial, k] = C
                
                dC = convolution_spikes_gaussian_derivative(time, spiketimes, latency_vector, 
                                                         sigma)
                
                dX[:, trial, k] = dC
                k+=1
                
        X = X*1000 # convert from spikes per ms to spikes per s
        dX = dX*1000 # convert from spikes per ms to spikes per s

        data = X - dX*latency_vector[:,None,None]
        data = data.reshape(len(time)*n_trials,n_neurons,order='F')
        return data, X, dX, time

def apply_approx_latency(X, dX,latency_vector):
    dims = X.shape
    data = X - dX*latency_vector[:,None,None]
    return data.reshape(dims[0]*dims[1],dims[2],order='F')

def get_windowed_data(tmin,tmax, BINSIZE, sorted_trials, sorted_spiketimes, 
                      sorted_neurons_id,sorted_trial_id, ind_sound_onset, data):
    
    time_window = np.linspace(tmin,tmax,int((tmax-tmin)/BINSIZE)+1)/1000
    
    
    speed_window = np.zeros((len(time_window),len(sorted_trials)))
    
    sorted_spiketimes_windowed = sorted_spiketimes.copy()
    cond = np.zeros(len(sorted_spiketimes)).astype(bool)
    for i in range(len(sorted_trials)):
        trial = sorted_trials[i]
        # print(ind_sound_onset[trial]-all_data['sound_onset'][trial]*1000)
        cond_trial =  (sorted_trial_id == i)
        cond_temp = (((sorted_spiketimes-ind_sound_onset[trial])>=tmin) & ((sorted_spiketimes-ind_sound_onset[trial])<=tmax) & (sorted_trial_id == i))
        cond = cond | cond_temp
        
        sorted_spiketimes_windowed[cond_trial] = sorted_spiketimes_windowed[cond_trial] - ind_sound_onset[trial] - tmin
        
    
        temp_speed  = data['speed'][trial]
        temp_time   = data['time'][trial] - ind_sound_onset[trial]/1000
        
        speed_window[:,i] = np.interp(time_window,temp_time,temp_speed)
        
        
        
    sorted_spiketimes_windowed = sorted_spiketimes_windowed[cond].copy()
    sorted_neurons_id_windowed = sorted_neurons_id[cond].copy()
    sorted_trial_id_windowed  = sorted_trial_id[cond].copy()
    
    
    y = speed_window.reshape(len(time_window)*len(sorted_trials), order = 'F')
    return y, speed_window, sorted_spiketimes_windowed, sorted_neurons_id_windowed, sorted_trial_id_windowed



def resample_speed(speeds, times, n_trials, t_min, t_max, bin_size):
    
    speeds_tr = speeds.reshape(int(len(speeds)/n_trials), n_trials, order = 'F')
    times_tr = times.reshape(int(len(speeds)/n_trials), n_trials, order = 'F')



    time_window = np.linspace(t_min,t_max,int((t_max-t_min)/bin_size)+1)
    speed_window = np.zeros((len(time_window),n_trials))

    for i in range(n_trials):
    
    
    
        temp_speed  = speeds_tr[:,i]
        temp_time   = times_tr[:,i]
        
        speed_window[:,i] = np.interp(time_window,temp_time,temp_speed)
    


    speeds_rs = speed_window.reshape(len(time_window)*n_trials, order = 'F')
    
    return speeds_rs
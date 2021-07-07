#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 18:24:48 2021

@author: JRA
"""
 # Import Stuff

import pickle
import numpy as np
import matplotlib.pyplot as plt


from glmnet import ElasticNet
from assistance_functions_neuron_warping import get_inds_alignment_points, get_sorted_data_for_alignment, plot_raw_rasters_with_alignment
from assistance_functions_neuron_warping import plot_error_alignment, get_speed_time_good_trials
from utils_lagged_regression import convolution_spikes_variable_kernel, get_data_matrix, get_windowed_data, resample_speed
from utils_lagged_regression import compute_optimal_lag_approx, apply_approx_latency, gaussian_filter_shifted
from lagged_regression_functions import generate_dataset
from lag_evolutionary_strategies import compute_optimal_lag_ES

import warnings
warnings.filterwarnings("ignore")
import timeit
import sys

sys.argv[1:]



# load data
Name =  'bM76_DRS_Day1_2017-02-23_19-34-57'

#machine_path = '/nfs/tank/renart/users/joao.afonso/'
machine_path = '/Volumes/'

folder_path = machine_path + 'Cerebro/Recording Data/' + Name + '/Predict from Neurons/'
del_group = 1
int_del = 4
sigma_conv = 140 #ms


load_name_1 = 'LagsAlgoTrials Delay_'+str(del_group)+'_'+str(int_del)+' Stdv_'+str(sigma_conv/1000)+'.npz'

load_name_2 = 'LagsAlgoResults Delay_'+str(del_group)+'_'+str(int_del)+' Stdv_'+str(sigma_conv/1000)+'.npz'

#loadName = 'TimeWarp Data Session_13 IntCond_1 Events_1  4.npz'
load_path_1 = folder_path + load_name_1
load_path_2 = folder_path + load_name_2

data_1 = np.load(load_path_1)
data_2 = np.load(load_path_2)

sorted_spikes = np.array(data_1['spikes']) * 1000 # in ms
sorted_neurons = np.array(data_1['n_ids'])
sorted_trials = np.array(data_1['t_ids'])


# get number of neurons and trials
n_trials = np.max(sorted_trials)
n_neurons = np.max(sorted_neurons)


# python indexes

sorted_neurons = sorted_neurons - 1
sorted_trials = sorted_trials - 1 





times = np.array(data_1['times']) * 1000 # in ms
speeds = np.array(data_1['speeds'])




sigma_conv = 140 #ms
t_max_speed = np.round(np.max(times))
t_min_speed = np.min(times)
kern_off = sigma_conv * 4
t_min_neurs = t_min_speed
t_max_neurs = t_max_speed + kern_off * 2 


bin_size = 20  # ms
n_bins_neurs = int((t_max_neurs - t_min_neurs) / bin_size)
n_bins_speed = int((t_max_speed - t_min_speed) / bin_size)
kern_off_bins = int(kern_off / bin_size)



n_neurons = len(np.unique(sorted_neurons))
n_trials  = len(np.unique(sorted_trials))
time = np.linspace(t_min_neurs,t_max_neurs,int(t_max_neurs/bin_size)+1)

X = np.zeros((len(time),n_trials,n_neurons))

trial = 0
neuron = 0

cond = (sorted_neurons==neuron) & (sorted_trials == trial)
spiketimes = sorted_spikes[cond]


lags_0 = np.zeros(n_bins_neurs+1)

all_lags = data_2['lags'] 
lags_t = all_lags[5,:]
#lags_1 =  np.hstack([np.zeros(kern_off_bins), lags_t, np.zeros(kern_off_bins )])

lags_1 =  np.hstack([lags_0[:134] + 100, lags_0[134:] - 100])


y_0 = np.zeros_like(time)
latency_spikes_0 = np.interp(spiketimes,time,lags_0)
for spike,shift in zip(spiketimes, latency_spikes_0):
    y_0 += gaussian_filter_shifted(time,spike,shift, sigma_conv)
    
    
    

y_1 = np.zeros_like(time)
latency_spikes_1 = np.interp(spiketimes,time,lags_1)
for spike,shift in zip(spiketimes, latency_spikes_1):
    y_1 += gaussian_filter_shifted(time,spike,shift, sigma_conv)    
    

plt.figure()
plt.subplot(2,1,1)
plt.plot(y_0[28:-28]*1000)
plt.plot(y_1[28:-28]*1000)
plt.xlim([0,208])
plt.xlabel('Samples (0.02 s)')
plt.ylabel('FR spikes/s')

plt.subplot(2,1,2)
plt.plot(lags_1[28:-28]/1000)
plt.xlim([0,208])
plt.ylim([-0.12,0.12])
plt.ylabel('lag s)')
plt.xlabel('Samples (0.02 s)')
plt.show()


def gaussian_filter_shifted(x,t_center,t_shift, sigma):
    x = x-t_center
    n_const = 1/(sigma*(2*np.pi)**0.5)
    a_exp =-1/2*((x-t_shift)/sigma)**2
    return n_const*np.exp(a_exp)










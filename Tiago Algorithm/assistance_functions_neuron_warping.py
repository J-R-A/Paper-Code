#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 22:35:08 2020

@author: tjpc
"""
import numpy as np
import matplotlib.pyplot as plt


def get_inds_alignment_points(data,t_end_trial = 5):
    ind_end_trial    = np.zeros(len(data['speed'])).astype(int)
    ind_sound_onset  = np.zeros(len(data['speed'])).astype(int)
    ind_sound_offset = np.zeros(len(data['speed'])).astype(int)
    ind_begin_accel  = np.zeros(len(data['speed'])).astype(int)
    ind_max_speed    = np.zeros(len(data['speed'])).astype(int)
    
    for i in range(len(data['speed'])):
        ind_end_trial[i]    = np.argmin(np.abs(data['time'][i]-t_end_trial)).astype(int)
        ind_sound_onset[i]  = np.argmin(np.abs(data['time'][i]-data['sound_onset'][i])).astype(int)
        ind_sound_offset[i] = np.argmin(np.abs(data['time'][i]-data['sound_offset'][i])).astype(int)
        ind_begin_accel[i]  = np.argmin(data['speed'][i][ind_sound_onset[i]:ind_sound_offset[i]])+ ind_sound_onset[i]
        ind_max_speed[i]    = np.argmax(data['speed'][i][ind_sound_offset[i]:ind_end_trial[i]])+ind_sound_offset[i]
    return ind_end_trial, ind_sound_onset, ind_sound_offset, ind_begin_accel, ind_max_speed



def get_sorted_data_for_alignment(all_data, TMAX, TMIN, NMAX_NEURONS, KEEP_TYPE_TRIALS = np.arange(11)):

    cond_sound   = all_data['trial_type'] == 1
    cond_correct = all_data['trial_outcome'] == 1
    cond_sound_start = all_data['i_sound_start'].copy()
    
    conds = cond_sound & cond_correct 
    REMOVE_TYPE_TRIALS = np.arange(11).astype(int)
    
    for KEEP_TYPE_TRIAL in KEEP_TYPE_TRIALS:
        REMOVE_TYPE_TRIALS = REMOVE_TYPE_TRIALS[REMOVE_TYPE_TRIALS!=KEEP_TYPE_TRIAL]
    
    if REMOVE_TYPE_TRIALS != []:
        for remove_trial in REMOVE_TYPE_TRIALS:
            conds_st = cond_sound_start!= remove_trial
            conds = conds & conds_st
    for i in range(len(conds)):
        if conds[i]:
            if (all_data['time'][i][-1]<TMAX/1000) | (all_data['speed'][i][min(TMAX,len(all_data['speed'][i])-1)]>5):
                conds[i] = False
    
    good_trials = np.arange(len(all_data['speed']))[conds]
    # =============================================================================
    # fig = plt.figure(figsize = (10,10))
    # axs = [fig.add_subplot(3,4,i+1) for i in range(11)]
    # for i in good_trials:
    #     axs[cond_sound_start[i]].plot(all_data['time'][i]*1000,all_data['speed'][i], c = 'b')
    # 
    # [axs[i].set_xlim([0,TMAX]) for i in range(11)]
    # =============================================================================

    
    # sort trials by sound start ...
    sorted_neurons_id   = []
    sorted_trials       = []
    sorted_spiketimes   = []
    sorted_trial_id     = []
    sorted_trial_cond   = []

    i_sound_start = cond_sound_start[good_trials]
    for trial_sound_start in np.unique(i_sound_start):
        for i in range(len(good_trials)):
            if i_sound_start[i] == trial_sound_start:
                sorted_neurons_id.extend(list(all_data['neurons_id'][all_data['neurons_trial_id'] == good_trials[i]]))
                sorted_trials.extend([good_trials[i]])
                sorted_spiketimes.extend(list(all_data['neurons_timestamps'][all_data['neurons_trial_id'] == good_trials[i]]))
                sorted_trial_id.extend(list(i*np.ones((all_data['neurons_trial_id'] == good_trials[i]).sum())))
                sorted_trial_cond.extend(list(all_data['neurons_i_sound_start'][all_data['neurons_trial_id'] == good_trials[i]]))
                
    sorted_neurons_id = np.array(sorted_neurons_id)  
    sorted_trials = np.array(sorted_trials).astype(int)  
    sorted_spiketimes = np.array(sorted_spiketimes)*1000    #convert to ms
    sorted_trial_id    = np.array(sorted_trial_id).astype(int)
    sorted_trial_cond = np.array(sorted_trial_cond).astype(int)
    # crop the spike times for a given interval 
    c_time_interval = (sorted_spiketimes > TMIN) & (sorted_spiketimes < TMAX)
    
    sorted_neurons_id = sorted_neurons_id[c_time_interval]
    sorted_trial_id   = sorted_trial_id[c_time_interval]
    sorted_spiketimes = sorted_spiketimes[c_time_interval]
    sorted_trial_cond = sorted_trial_cond[c_time_interval]
    # Remove certain neurons
    if NMAX_NEURONS!= None:
        c_neurons = sorted_neurons_id < NMAX_NEURONS
        
        sorted_neurons_id = sorted_neurons_id[c_neurons]
        sorted_trial_id   = sorted_trial_id[c_neurons]
        sorted_spiketimes = sorted_spiketimes[c_neurons]
        sorted_trial_cond = sorted_trial_cond[c_neurons]

    return sorted_neurons_id, sorted_trials, sorted_trial_id, sorted_spiketimes, sorted_trial_cond

def get_sorted_data_around_event(sorted_neurons_id, sorted_trials, sorted_trial_id, sorted_spiketimes, sorted_trial_cond, align_inds):
    pass

    # return sorted_neurons_id, sorted_trials, sorted_trial_id, sorted_spiketimes, sorted_trial_cond



def plot_raw_rasters_with_alignment(data, sorted_trials,  ind_sound_onset, ind_sound_offset, ind_begin_accel, ind_max_speed):
    raster_kws = dict(s=4, c='k', lw=0, alpha = 0.4)
    align_typess = [0,1,2,3]
    colors = ['r','b','g','b']
    colors = ['C0','C2','C1','C3']
    for align_types in align_typess:
        align_types = [align_types]
        fig = plt.figure(figsize=(35, 20))
        axes = []
        
        for i in range(20):
            axes.append(fig.add_subplot(4,5,i+1))
            # raw data
            int_Data_raw = data.select_neurons(i)
            spks = int_Data_raw.spiketimes
            trials = int_Data_raw.trials
            axes[i].scatter(spks,trials,**raster_kws ) 
        
            for align_type in align_types:
                if align_type == 0:
                    align_sound_onset  = ind_sound_onset[sorted_trials]
                    axes[i].scatter(align_sound_onset, range(align_sound_onset.size), c = colors[align_type], alpha=0.5, lw=0)
        
                if align_type == 1:
                    align_sound_offset = ind_sound_offset[sorted_trials]
                    axes[i].scatter(align_sound_offset, range(align_sound_offset.size), c = colors[align_type], alpha=0.5, lw=0)
        
                if align_type == 2:
                    align_begin_accel  = ind_begin_accel[sorted_trials]
                    axes[i].scatter(align_begin_accel, range(align_begin_accel.size), c = colors[align_type], alpha=0.5, lw=0)
        
                if align_type == 3:
                    align_peak_speed   = ind_max_speed[sorted_trials]
                    axes[i].scatter(align_peak_speed, range(align_peak_speed.size), c = colors[align_type], alpha=0.5, lw=0)
        
            if i >= 15:
                axes[i].set_xlabel('time (ms)', fontsize = 20)
            if (i % 5)== 0:
                axes[i].set_ylabel('trials', fontsize = 20)
            # axes[i].set_title('Neuron ind:' + str(i))
            

def get_speed_time_good_trials(data, t_end_trial, good_trials, ind_end_trial):
    speeds = np.zeros((len(good_trials),t_end_trial))
    t_types = np.zeros(len(good_trials))

    count = 0
    for i in range(len(data['speed'])):
        if np.isin(i, good_trials):
            # ind_end_trial    = np.argmin(np.abs(data['time'][i]-t_end_trial)).astype(int)
            # print(len(data['speed'][i][:ind_end_trial[i]]))
            speeds[count,:] = data['speed'][i][:t_end_trial]
            t_types[count] = data['i_sound_start'][i]

            count += 1
    return speeds, t_types


def plot_error_alignment(binned,time = None):
    figs = plt.figure(figsize= (14,10))
    
    axs = [figs.add_subplot(7,7,i+1) for i in range(49)]
    
    for j in range(49):
        # for i in range(binned_cond.shape[0]):
        #     axs[j].plot(binned_cond[i,:,j], alpha = 0.05, c ='b')
            
        mean_j = binned[:,:,j].mean(0)
        std_j = binned[:,:,j].std(0)
        m_std = binned[:,:,j].mean(0) - std_j
        M_std = binned[:,:,j].mean(0) + std_j
        
        mean_j = np.quantile(binned[:,:,j], 0.5, axis=0)
        m_std  = np.quantile(binned[:,:,j], 0.1, axis=0)
        M_std  = np.quantile(binned[:,:,j], 0.9, axis=0)
        axs[j].fill_between(np.arange(len(mean_j))/4.1*10,m_std, M_std, alpha = 0.4,  color='C1')
        axs[j].plot(np.arange(len(mean_j))/4.1*10,binned[:,:,j].mean(0), c = 'C0')
        if j<= 6:
            axs[j].set_ylim([0,8])
        else:
            axs[j].set_ylim([0,5])
    plt.tight_layout()
        

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 01:42:36 2020

@author: tjpc
"""


from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import json
import pickle

path ='/home/tjpc/Desktop/work_projects/Joao_project/joao_data/'

raw_data_path = 'Behaviour_and_Neural_Data/'


sessions = ['bM67_DRS_Day4_2017-02-19_18-59-10.mat',
            'bM69_DRS_Day1_2017-02-14_13-21-36.mat',
            'bM76_DRS_Day1_2017-02-23_19-34-57.mat',
            'bM76_DRS_Day2_2017-02-24_18-49-30.mat']


# sessions = ['bM69_DRS_Day1_2017-02-14_13-21-36.mat',
#             'bM76_DRS_Day1_2017-02-23_19-34-57.mat',
#             'bM76_DRS_Day2_2017-02-24_18-49-30.mat']

conds_by_type = [1,5,10,15,20,25,30,35,40,45,50]
# fig = plt.figure(figsize = (12,14))
# axs = [fig.add_subplot(4,3,i+1) for i in range(len(conds_by_type))]


session_counter = 0
for session in sessions:
    
    raw_data = loadmat(path + raw_data_path + session)
    print('Start extracting data session ' + str(session_counter))
    print('')
    print('')
    neurons = raw_data['neurons']
    behaviour = raw_data['organizedTrials']
    data = {}
    data['speed'] = []
    data['time'] = []
    data['time_interval'] = []
    data['length'] = []
    data['sound_onset'] = []
    data['sound_offset'] = []
    data['area_start'] = []
    data['area_end'] = []
    data['trial_type'] = []
    data['trial_outcome'] = []
    data['sound_start'] = []
    data['i_sound_start'] = []
    
    data['neurons_timestamps'] = []
    data['neurons_lengthstamps'] = []
    data['neurons_trial_type'] = []
    data['neurons_trial_outcome'] = []
    data['neurons_i_sound_start'] = []
    data['neurons_id'] = []
    data['neurons_trial_id'] = []

    for i in range(behaviour[0].shape[0]):
        print('start extracting trial ' + str(i))
        data['speed'].append(behaviour[0][i][0])
        data['time'].append(behaviour[1][i][0][:-1]-behaviour[1][i][0][:-1][0])
        data['time_interval'].append([behaviour[1][i][0].min(), behaviour[1][i][0].max()])
        data['length'].append(behaviour[2][i][0][:-1]-behaviour[2][i][0][:-1][0])
        data['sound_onset'].append(behaviour[4][i][0][0]-behaviour[1][i][0][:-1][0])
        data['sound_offset'].append(behaviour[4][i][0][1]-behaviour[1][i][0][:-1][0])
        data['area_start'].append(behaviour[5][i][0][0]-behaviour[1][i][0][:-1][0])
        data['area_end'].append(behaviour[5][i][0][1]-behaviour[1][i][0][:-1][0])
        data['trial_type'].append(behaviour[6][i][0][0])
        data['trial_outcome'].append(behaviour[7][i][0][0])
        data['sound_start'].append(behaviour[-1][i][0,0])
        data['i_sound_start'].append(np.where(conds_by_type == behaviour[-1][i][0,0])[0][0])
        
        t = data['time_interval'][i]

        for j in range(neurons.shape[0]):
            inds = ((neurons[j,0] >= t[0]) & (neurons[j,0] < t[1])).reshape(-1)
            placeholder = np.ones(sum(inds)).astype(int)
            data['neurons_timestamps'].extend(neurons[j,0].reshape(-1)[inds] - t[0])
            data['neurons_i_sound_start'].extend(data['i_sound_start'][i]*placeholder)
            data['neurons_id'].extend(j*placeholder)
            data['neurons_trial_type'].extend(data['trial_type'][i]*placeholder)
            data['neurons_trial_outcome'].extend(data['trial_outcome'][i]*placeholder)
            data['neurons_trial_id'].extend(i*placeholder)
    data['neurons_timestamps']    = np.array(data['neurons_timestamps'])
    data['neurons_trial_type']    = np.array(data['neurons_trial_type'])
    data['neurons_trial_outcome'] = np.array(data['neurons_trial_outcome'])
    data['neurons_i_sound_start'] = np.array(data['neurons_i_sound_start'])
    data['neurons_id']            = np.array(data['neurons_id'])
    data['neurons_trial_id']      = np.array(data['neurons_trial_id'])
    
    data['sound_onset']   = np.array(data['sound_onset'])
    data['sound_offset']  = np.array(data['sound_offset'])
    data['area_start']    = np.array(data['area_start'])
    data['area_end']      = np.array(data['area_end'])
    data['trial_type']    = np.array(data['trial_type'])
    data['trial_outcome'] = np.array(data['trial_outcome'])
    data['sound_start']   = np.array(data['sound_start'])
    data['i_sound_start'] = np.array(data['i_sound_start'])
    
    spath = path  + raw_data_path + 'session_'+str(session_counter)
    pickle.dump(data, open( spath+'.p', "wb" ) )

    # with open(spath, 'w') as outfile:
    #     json.dump(data, outfile)
    session_counter += 1


# import pickle
# import numpy as np
# import matplotlib.pyplot as plt


# from glmnet import ElasticNet

# # import scipy
# # MAD = scipy.stats.median_absolute_deviation
# #custom imports
# from assistance_functions_neuron_warping import get_inds_alignment_points, get_sorted_data_for_alignment, plot_raw_rasters_with_alignment
# from assistance_functions_neuron_warping import plot_error_alignment, get_speed_time_good_trials
# from utils_lagged_regression import convolution_spikes_variable_kernel, get_data_matrix, get_windowed_data
# from utils_lagged_regression import compute_optimal_lag_approx, apply_approx_latency
# from lagged_regression_functions import generate_dataset
# from lag_evolutionary_strategies import compute_optimal_lag_ES

# import warnings
# warnings.filterwarnings("ignore")
# import timeit

# #PATH DATA
# path ='/home/tjpc/Desktop/work_projects/Joao_project/joao_data/'

# raw_data_path = 'Behaviour_and_Neural_Data/'
# sessions = ["session_0.p","session_1.p","session_2.p","session_3.p"]
# all_data = pickle.load( open(path + raw_data_path + sessions[2], "rb" ))

# ##############################################################################################
# # organizing data
# # you don't need to change this
# NMAX_NEURONS = None
# REMOVE_TYPE_TRIALS = [5]
# TMIN = 0
# TMAX = 4100 #ms
# sorted_neurons_id, sorted_trials, sorted_trial_id, sorted_spiketimes, sorted_trial_cond = get_sorted_data_for_alignment(all_data, TMAX, TMIN, NMAX_NEURONS, REMOVE_TYPE_TRIALS)
# ind_end_trial, ind_sound_onset, ind_sound_offset, ind_begin_accel, ind_max_speed = get_inds_alignment_points(all_data,TMAX/1000)

# conds = np.unique(sorted_trial_cond)


# #########################################################################
# #WINDOW TRIALS AROUND SOUND ON SET
# # Other alignment points can be given by chaging ind_sound_onset
# # in the get_windowed_data function

# sigma_conv = 160 #ms
# tmin = -500
# # tmax    =  3300
# tmax    =  1500
# TMIN = 0
# TMAX = tmax-tmin#ms

# BINSIZE = 20  # ms
# NBINS = int((TMAX - TMIN) / BINSIZE)


# y,speed_window,sorted_spiketimes_windowed,sorted_neurons_id_windowed,sorted_trial_id_windowed = get_windowed_data(tmin,tmax, BINSIZE, sorted_trials, sorted_spiketimes,
#                       sorted_neurons_id,sorted_trial_id, ind_sound_onset, all_data)

# # Everything from this line upwards can be replaced by your functions in matlab.
# # In order to do it you need:
# # y -> the velocity vector of the animal with all trials concatenated.

# # speed_window -> a velocity matrix without the trials concatenated. You don't need this one, but it's nice for plotting.

# # sorted_spiketimes_windowed -> A vector with all the spike times (in ms) in a certain window of time. Note that the spiketimes have to be relative to the start of the window
# # for example if the window of time is -500 ms before the sound onset to 1500 ms after for a given trial a spike happens 499 ms before the sound onset the value
# # of the timestamp as to be 1 ms.

# # sorted_trial_id_windowed -> A index vector with the length = n_of_spikes, with the corresponding trial index (without missing values)

# # sorted_neurons_id_windowed -> A index vector with the length = n_of_spikes, with the corresponding neuron index (without missing values)
# #####################################################################################
# #  Generate dataset
# # you can delete this part just replace pred_y with the y (the real speed) from the previous section, in the compute_optimal_lag_ES function.
# betas, latency_vector, pred_y, vel, _,_,_= generate_dataset(speed_window, sorted_spiketimes_windowed, sorted_neurons_id_windowed,
#                                                         sorted_trial_id_windowed,BINSIZE, TMIN,
#                                                         TMAX, sigma = sigma_conv, nolag = False)



###################################################   MY VERSION   ###############################################################

# Import Stuff

import pickle
import numpy as np
import matplotlib.pyplot as plt


from glmnet import ElasticNet
from assistance_functions_neuron_warping import get_inds_alignment_points, get_sorted_data_for_alignment, plot_raw_rasters_with_alignment
from assistance_functions_neuron_warping import plot_error_alignment, get_speed_time_good_trials
from utils_lagged_regression import convolution_spikes_variable_kernel, get_data_matrix, get_windowed_data, resample_speed
from utils_lagged_regression import compute_optimal_lag_approx, apply_approx_latency
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
load_name = 'LagsAlgoTrials Delay_'+str(del_group)+'_'+str(int_del)+'.npz'
save_name = 'LagsAlgoResults Delay_'+str(del_group)+'_'+str(int_del)+'.npz'

#loadName = 'TimeWarp Data Session_13 IntCond_1 Events_1  4.npz'
load_path = folder_path + load_name
save_path = folder_path + save_name

data = np.load(load_path)
sorted_spikes = np.array(data['spikes']) * 1000 # in ms
sorted_neurons = np.array(data['n_ids'])
sorted_trials = np.array(data['t_ids'])


# get number of neurons and trials
n_trials = np.max(sorted_trials)
n_neurons = np.max(sorted_neurons)


# python indexes

sorted_neurons = sorted_neurons - 1
sorted_trials = sorted_trials - 1 





times = np.array(data['times']) * 1000 # in ms
speeds = np.array(data['speeds'])




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


# transform speeds

speed_rs = resample_speed(speeds, times, n_trials, t_min_speed, t_max_speed, bin_size)
    

#####################################################################################

# These hyper parameters are just suggestions that seem to work. I would start playing around with
# n_iter and max_iterations tell you how good is your estimation and increase the computational time
# The lrate and std_perturbations seem to work around these values but feel free to explore other regions of the parameter space.
# I wouldn't mess with the annelings in the begginings
# reg_smoothess seems small but it actually as a good impact due to the scalling of the lag function.

n_iter = 6                                 # number of cycles computing the betas and the lag function


lrate = 0.025                            # Learning rate for the ES algorithm
std_perturbations = 0.15                    # Standard deviation of the perturbations right now the rewards are
                                            # not being standardized so this parameter doesn't change the actual learning rate
n_perturbations = 15                        # N of random perturbations in the ES
max_iterations = 100                        # Num of iterations in the ES
anneling_lrate = 1                          # Anneling rate of the learning rate
anneling_std = 1                            # Anneling rate of the Standard deviation of the perturbations
reg_smoothess = 0.009                       # regularization parameter of the lag function5
                                            # This controls how smooth the lag function is. It has also a minor effect on it's magnitute


# This is the standart glmnet model in python. You have to change a line of code in the package (i can help you)
# you can installing it using: pip install glmnet you need some requirements and a python version of 3.8xx
# You can change the regularization parameters of the glm net.
# Alpha is how much lasso vs ridge you want
# lambda_path is how big is your regularization.
regression_model = ElasticNet(alpha = 0, lambda_path=np.array([0]) ,fit_intercept = False, verbose = False) # check this part


# lag_vectors -> saves the lag vectors of all the iterations in a matrix iteration x time
# regressors_vectors -> saves the regressors of all the iterations in a matrix iteration x beta index
# neural_signals -> saves the convolved neural signals in a tensor $n_iter x time x trials x neurons
lag_vectors = np.zeros((n_iter, n_bins_speed + 1))
regressors_vectors = np.zeros((n_iter,n_neurons))
neural_signals = np.zeros((n_iter,n_bins_speed + 1,n_trials,n_neurons))



dummy_lag_vector = np.zeros(n_bins_neurs+1)



start_0 = timeit.default_timer()

for i in range(n_iter):
    start_1 = timeit.default_timer()
    print('starting iteration num', str(i))
    # print(lag_vector)
    _, X, time = get_data_matrix(sorted_spikes, sorted_neurons,sorted_trials, t_min_neurs, t_max_neurs, dummy_lag_vector, sigma_conv,bin_size, type_kernel = 'real')

    
    X_1 = X[ kern_off_bins : - kern_off_bins ,:,:]
    X_reg = X_1.reshape(len(X_1[:,0,0])*len(X_1[0,:,0]),len(X_1[0,0,:]),order='F')

    fit = regression_model.fit(X_reg.copy(), speed_rs.copy())
    regressors = fit.coef_

    initial_lag_vector = np.zeros(n_bins_speed + 1)
    # replace pred_y with the y (the real speed) from the previous section, in the compute_optimal_lag_ES function.


    lag_vector = compute_optimal_lag_ES(speed_rs, initial_lag_vector, regressors, sorted_spikes,
                                        sorted_neurons,sorted_trials,
                                        t_min_neurs, t_max_neurs, kern_off_bins, sigma_conv, bin_size, lrate = lrate,
                                        std_perturbations= std_perturbations, n_perturbations = n_perturbations,
                                        max_iterations= max_iterations, anneling_lrate = anneling_lrate, anneling_std = anneling_std,
                                        reg_smoothess = reg_smoothess, minimize = True)

    lag_vector = lag_vector*100 # converts the lag back to ms the algorithm uses a scalling factor of 100 in the lag.
    dummy_lag_vector =  np.hstack([np.zeros(kern_off_bins), lag_vector, np.zeros(kern_off_bins )])
    
    
    lag_vectors[i,:] = lag_vector.copy()
    regressors_vectors[i,:] = regressors.copy()
    neural_signals[i,:,:,:] = X_1.copy()
    
    print('time it took inner iter: ', str(timeit.default_timer() - start_1))


print('time it took outer iter: ', str(timeit.default_timer() - start_0))



#np.savez(save_path, coefs = regressors_vectors, neurons = neural_signals, lags = lag_vectors, speeds = speed_rs)



    






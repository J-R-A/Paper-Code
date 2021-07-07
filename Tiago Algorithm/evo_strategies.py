#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 17:27:23 2020

@author: tjpc
"""

import numpy as np
import numpy.random as npr

#####################################################
# TODO: Implement parallel properly...
# TODO: Code pass_perturbations as requiring a list of divisions
# TODO: Implement stopping conditionsfitness[1,i]
# TODO: Implement linear anneling
# TODO: Implement adaptative leaning rate
import matplotlib.pyplot as plt
#####################################################

def evolutionary_strategies(data, cost_function, optim_params, lrate, std_perturbations,
                            n_perturbations, max_iterations, reltol = 1E-8 ,
                            abstol = 1E-8, anneling_lrate = 1, anneling_std = 1,
                            paralel = False, minimize = True,
                            min_lrate_anneling = 1E-6,
                            min_std_anneling = 1E-6,
                            fixed_params = None, mirror_perturbations = True,
                            pass_perturbations = False):
    # print('anneling_lrate',anneling_lrate)
    # print('anneling_std',anneling_std)
    # print('minimize',minimize)
    # print(n_perturbations)
    # print('std_perturbations', std_perturbations)
    # print(mirror_perturbations)
    pararm_dims = optim_params.shape
    fitness_mean_hist = np.zeros(max_iterations)
    fitness_std_hist = np.zeros(max_iterations)




    for i in range(max_iterations):

        lrate = _apply_anneling(lrate, anneling_lrate, min_lrate_anneling) # replicate this
        std_perturbations = _apply_anneling(std_perturbations, anneling_std, min_std_anneling) # replicate this
        # print('std_perturbations aneel', std_perturbations)
        # print(lrate)
        perturbation_vector = _create_perturbation(std_perturbations, pararm_dims, n_perturbations)

        # print('std_perturbations', perturbation_vector.std())

        if paralel:
            fitness = _run_parallel(data, cost_function, optim_params, fixed_params,
                                    perturbation_vector, n_perturbations,
                                    mirror_perturbations, pass_perturbations)
        else:
           fitness = _run_serial(data, cost_function, optim_params, fixed_params,
                                 perturbation_vector, n_perturbations,
                                 mirror_perturbations, pass_perturbations)

# =============================================================================
#         To debugg and later remove
#        if (i % max(max_iterations//20,1)) == 0:
#            if pass_perturbations:
#                fig = plt.figure()
#                ax1 = fig.add_subplot(1,2,1)
#                ax2 = fig.add_subplot(1,2,2)

#                ax1.plot(optim_params[:fixed_params['size_lag']])
#                ax2.plot(optim_params[fixed_params['size_lag']:])

#            else:
                #pass
#                plt.plot(optim_params.copy())
#                plt.ylim([-2.4,2.4])

# =============================================================================

        if pass_perturbations:
            optim_params = _update_optim_params_different_lr(optim_params, lrate,
                                                             fitness, perturbation_vector,
                                                             n_perturbations,
                                                             std_perturbations,
                                                             minimize, mirror_perturbations,
                                                             fixed_params['size_lag'])
        else:
            optim_params = _update_optim_params(optim_params, lrate, fitness, perturbation_vector,
                                            n_perturbations, std_perturbations, minimize, mirror_perturbations)
        fitness_mean_hist[i]=fitness.mean()
        fitness_std_hist[i]=fitness.std()


        print('current generation', str(i+1),' of ',str(max_iterations),
              ' with fitness: ', str(round(fitness_mean_hist[i],2)),'(',str(round(np.log10(fitness_mean_hist[i]),1)),')',
              ' and  std: ', str(round(fitness_std_hist[i],2)),'(',str(round(np.log10(fitness_std_hist[i]),1)),')')
        print('Learning rate:',  str(round(lrate,6)), ', std perturbations: ', str(round(std_perturbations,6)))
        print('')


# =============================================================================
#         To debugg and later remove
        if (i % max(max_iterations//20,1)) == 0:
            if pass_perturbations:
                ax1.plot(optim_params[:fixed_params['size_lag']])
                ax2.plot(optim_params[fixed_params['size_lag']:])
                ax1.set_ylim([-2.4,2.4])
                ax2.set_ylim([-15,5])
                ax1.set_title('lag')
                ax2.set_title('regressors')

                plt.show()
                plt.close()
            else:

                plt.plot(optim_params.copy())
                plt.ylim([-2.4,2.4])


                plt.show()
                plt.close()
# # =============================================================================
        # print('optimize params', optim_params)


        # print('params: ', i, optim_params)
        if _stopping_conditions(fitness_mean_hist[:i+1],reltol, abstol):        # implement stopping conditions
            break
    #implement stoppage conditions mensage
    stoppage_conds = []
    return fitness_mean_hist, fitness_std_hist, optim_params, stoppage_conds

def _update_optim_params(optim_params, lrate, fitness, perturbation_vector,
                         n_perturbations, std_perturbations, minimize,
                         mirror_perturbations):

    gradient = _compute_gradient(fitness,perturbation_vector, n_perturbations,
                                 std_perturbations, mirror_perturbations)
    # print('test',np.abs(optim_params[101:]).mean())
    # print('test grad',np.abs(lrate * gradient[101:]).mean())
    if minimize:
        # print(lrate)
        # print(optim_params[:3])
        # print((gradient)[:3].mean())
        gamma = 1#np.abs(gradient).max()/10
        print('max_update: ', lrate/gamma*gradient.max())
        return optim_params - lrate * gradient/gamma
    else:
        return optim_params + lrate * gradient

def _update_optim_params_different_lr(optim_params, lrate, fitness, perturbation_vector,
                                      n_perturbations, std_perturbations, minimize,
                                      mirror_perturbations, split_params):

    gradient = _compute_gradient(fitness,perturbation_vector, n_perturbations,
                                 std_perturbations, mirror_perturbations)
    # print('test',np.abs(optim_params[101:]).mean())
    # print('test grad',np.abs(lrate * gradient[101:]).mean())
    if minimize:
        # print(lrate)
        # print(optim_params[:3])
        # print((gradient)[:3].mean())


        gamma1 = 1#np.abs(gradient[:split_params]).max()
        gamma2 = 1#np.abs(gradient[split_params:]).max()

        print('max_update lag: ', lrate*np.abs(gradient[:split_params]).max())
        print('max_update reg: ', lrate*np.abs(gradient[:split_params]).max())

        # print('max_update betas: ', lrate/gamma1*gradient.max())

        optim_params[:split_params] = optim_params[:split_params]  - lrate * gradient[:split_params]
        optim_params[split_params:] = optim_params[split_params:]  - lrate * gradient[split_params:]

        return optim_params
    else:
        optim_params[:split_params] = optim_params[:split_params]  + lrate * gradient[:split_params]
        optim_params[split_params:] = optim_params[split_params:]  + lrate * gradient[split_params:]
        return optim_params


def _compute_gradient(fitness,perturbation_vector, n_perturbations,
                      std_perturbations, mirror_perturbations):
    # check order os multiplication
    # print(fitness[0,:4])
    # print(fitness[1,:4])
    print('mean fitness MIRROR:', np.abs(fitness[0,:]-fitness[1,:]).mean(), (np.abs(fitness[0,:]).std()+np.abs(fitness[1,:]).std())/2)

    # fitness = (fitness-fitness.mean())/fitness.std()
    print('mean fitness:', np.abs(fitness).mean())

    if mirror_perturbations:
        fitness = fitness[0,:] - fitness[1,:]
        # fitness = (fitness-fitness.mean())/fitness.std()
        # # print('print fittness',fitness.mean(), fitness.std())
        # # print(std_perturbations)
        # # print(np.matmul(perturbation_vector.T,fitness))
        return 1/(2*n_perturbations*std_perturbations**2) * np.matmul(perturbation_vector.T,fitness)
    else:
        return 1/(2*n_perturbations*std_perturbations**2) * np.matmul(perturbation_vector.T,fitness)

# _apply_anneling(lrate, anneling_lrate, min_lrate_anneling, i)
def _apply_anneling(rate, anneling, min_value):
    if anneling == 'linear':
        return anneling
    else:
        return max([anneling*rate,min_value])


def _create_perturbation(sigma, param_dims, n_perturbations):
    return sigma*npr.randn(n_perturbations, *param_dims)

def _stopping_conditions(fitness,reltol, abstol):
    # implement stopping conditionsfitness[1,i]
    return False



def _run_parallel(data, cost_function, optim_params, fixed_params, perturbation_vector, n_perturbations, pass_perturbations):
    return np.zeros((n_perturbations))


def _run_serial(data, cost_function, optim_params, fixed_params, perturbation_vector,
                n_perturbations, mirror_perturbations, pass_perturbations):
    if mirror_perturbations:
        fitness = np.zeros((2,n_perturbations))
        for i in range(n_perturbations):
            if pass_perturbations:
                eval_params = [optim_params.copy(), perturbation_vector[i,:]]
                fitness[0,i] = cost_function(data, eval_params, fixed_params)
                eval_params = [optim_params.copy(), - perturbation_vector[i,:]]
                fitness[1,i] = cost_function(data, eval_params, fixed_params)

            else:
                # print(perturbation_vector.max(), perturbation_vector.min())
                eval_params = optim_params.copy() + perturbation_vector[i,:]
                fitness[0,i] = cost_function(data, eval_params, fixed_params)
                # print('check + :',fitness[0,i])
                eval_params = optim_params.copy() - perturbation_vector[i,:]
                fitness[1,i] = cost_function(data, eval_params, fixed_params)
                # print('check - :',fitness[1,i])
                # print('')
        return fitness
    else:
        fitness = np.zeros((n_perturbations))
        for i in range(n_perturbations):
            if pass_perturbations:
                eval_params = [optim_params.copy(), perturbation_vector[i,:]]
                fitness[i] = cost_function(data, eval_params, fixed_params)
            else:
                eval_params = optim_params.copy() + perturbation_vector[i,:]
                fitness[i] = cost_function(data, eval_params, fixed_params)
        return fitness

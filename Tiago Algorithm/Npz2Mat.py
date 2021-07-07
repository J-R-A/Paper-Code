import numpy as np
import scipy.io
import re

Name =  'bM76_DRS_Day1_2017-02-23_19-34-57'
folderPath = '/Volumes/Cerebro/Recording Data/'+Name+'/Predict from Neurons/'
del_group = 1
int_del = 4
stdv = 0.14


file_path = '/Volumes/Cerebro/Recording Data/'+Name+'/Predict from Neurons/'
load_name = 'LagsAlgoResults Delay_'+str(del_group)+'_'+str(int_del)+' Stdv_'+str(stdv)+'.npz'
save_name = 'LagsAlgoResults Delay_'+str(del_group)+'_'+str(int_del)+' Stdv_'+str(stdv)+'.mat'
save_name_sup = 'LagsAlgoResults_Sup Delay_'+str(del_group)+'_'+str(int_del)+' Stdv_'+str(stdv)+'.mat'

load_path = file_path + load_name
save_path = file_path + save_name
save_path_sup = file_path + save_name_sup




with open('/Users/JRA/Downloads/output_158-0.txt') as f:
    lines = f.readlines()

fitting_data = {'mean_fitness_mirror':[],'mean_fitness':[],'max_update':[],'current_iter':[],'learning_rate':[]}

for l in lines:

    if l.startswith('mean fitness MIRROR:'):
        fitting_data['mean_fitness_mirror'].append([float(i) for i in re.findall("\d+\.\d+", l)])


    if l.startswith('mean fitness:'):
        fitting_data['mean_fitness'].append(float(re.findall("\d+\.\d+", l)[0]))


    if l.startswith('max_update:'):
        fitting_data['max_update'].append(float(re.findall("\d+\.\d+", l)[0]))

    if l.startswith('current generation'):
        fitting_data['current_iter'].append(int(re.findall(r'\d+', l)[0]))



    if l.startswith('Learning rate:'):
        fitting_data['learning_rate'] = [float(i) for i in re.findall("\d+\.\d+", l)]





file = np.load(load_path)
scipy.io.savemat(save_path, file)
scipy.io.savemat(save_path_sup, fitting_data)

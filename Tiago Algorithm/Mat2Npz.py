import numpy as np
import scipy.io


Name =  'bM76_DRS_Day1_2017-02-23_19-34-57'
folderPath = '/Volumes/Cerebro/Recording Data/'+Name+'/Predict from Neurons/'
del_group = 1
stdv = 0.14

for int_del in range(1,6):
    print(int_del)
    loadName = 'LagsAlgoTrials Delay_'+str(del_group)+'_'+str(int_del)+' Stdv_'+str(stdv)+'.mat'
    saveName = 'LagsAlgoTrials Delay_'+str(del_group)+'_'+str(int_del)+' Stdv_'+str(stdv)
    loadPath = folderPath + loadName
    savePath = folderPath + saveName

    # load matfile
    matFile = scipy.io.loadmat(loadPath)

    # saved arrays
    spikes = matFile['mergedSpikes'][:,0]
    n_ids = matFile['mergedNeuronsId'][:,0]
    t_ids = matFile['mergedTrialsId'][:,0]
    times = matFile['mergedTimes'][:,0]
    speeds = matFile['mergedSpeeds'][:,0]

    # save npz file

    np.savez(savePath, spikes = spikes, n_ids = n_ids, t_ids = t_ids, times = times, speeds = speeds)

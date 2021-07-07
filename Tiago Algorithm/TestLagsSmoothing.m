

machine = 1;
session = 13;
delGroup = 1;
intDel = 4;
stdv = 0.14;

if machine == 1
    machinePath = '/Volumes/';
elseif machine == 2
    machinePath = '/nfs/tank/renart/users/joao.afonso/';
end

load([machinePath,'Cerebro/Recording Data/SessionsMap.mat'])

sess = 1;
for m = 1:length(sessionsMap(:,1))
    for s = 1:length(sessionsMap{m,2}(1,:))
        Names{sess,1} = sessionsMap{m,2}{1,s};
        sess = sess+1;
    end
end


Name =  Names{session,1};



load([machinePath,'Cerebro/Recording Data/',Name,'/Predict from Neurons/LagsAlgoTrials Delay_',num2str(delGroup),'_',num2str(intDel),' Stdv_',num2str(stdv),'.mat'])


 
sorted_spikes = mergedSpikes;
sorted_neurons = mergedNeuronsId;
sorted_trials = mergedTrialsId;
times = mergedTimes;
speeds = mergedSpeeds;

n_trials = max(sorted_trials);
n_neurons = max(sorted_neurons);


t_max_speed = round(max(times),2);
t_min_speed = min(times);
kern_off = stdv * 4;
t_min_neurs = t_min_speed;
t_max_neurs = t_max_speed + kern_off * 2; 


bin_size = 0.02;
n_bins_neurs = (t_max_neurs - t_min_neurs) / bin_size;
n_bins_speed = (t_max_speed - t_min_speed) / bin_size;
kern_off_bins = round(kern_off / bin_size);

time = linspace(t_min_neurs,t_max_neurs,n_bins_neurs+1);
X = zeros(length(time),n_trials,n_neurons);

lags = ones(size(time)) * 0.1;

for t = 1:n_trials
    k = 1;
    for n = 1:n_neurons
        cond_map = (sorted_neurons == n) & (sorted_trials == t);
        spike_times = sorted_spikes(cond_map);
        
        conv = convolution_spikes_variable_kernel(time, spike_times, lags, stdv);
        
        X(:, t, k) = conv;
        k = k+1;
        
    end
end


x_final = X(kern_off_bins+1:end-kern_off_bins,:,:);









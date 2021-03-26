function PredictSpeedClusterPaper(session,event,intCond,stdv,elnetVal,spdSmt,machine)
clearvars -except session event intCond stdv elnetVal spdSmt machine 

procTime = tic;


if machine == 1
    machinePath = '/Volumes/';
elseif machine == 2
    machinePath = '/nfs/tank/renart/users/joao.afonso/';
end

load([machinePath,'Cerebro/Recording Data/SessionsMap.mat'])



% Extract all session Names from sessionsMap cell

sess = 1;
for m = 1:length(sessionsMap(:,1))
    for s = 1:length(sessionsMap{m,2}(1,:))
        Names{sess,1} = sessionsMap{m,2}{1,s};
        sess = sess+1;
    end
end


% Select name of the session to use

Name =  Names{session,1};
disp(['Processing session: ',  Name])


% Create features matrix 
 
[Speeds,pMats,trialInfo,trialEvents,Distances,paramInfo] = MakeFeaturesMatrixSpdPredictPaper(Name,event,intCond,stdv,spdSmt,machinePath);

[Laps, dummyLaps, allPredMat, allTrials, trialEvents, allDists, coefs, fitInfo] = PredictSpeedTimePaper(pMats,trialInfo,trialEvents,Speeds,Distances,elnetVal);


paramInfo.elnet = elnetVal;
spdSmt = paramInfo.spdSmt;

folderPath=[machinePath,'Cerebro/Recording Data/',Name,'/Predict from Neurons/'];
filename= [folderPath,'SpdPredResultsTimePaper Cond_',num2str(intCond),' TSeg_',num2str(event),' spdSmt_', num2str(spdSmt),' Sigma_',num2str(stdv),' Alpha_',num2str(elnetVal),'.mat'];
save(filename,'Laps', 'dummyLaps', 'allPredMat', 'allTrials', 'trialEvents', 'allDists', 'coefs', 'fitInfo','paramInfo','-v7.3')

%clearvars -except s sessionsMap  machinePath intCond stdv kWdthS event shuffleTrials elnet Names
%end

disp(['The entire process took: ',num2str(toc(procTime)/60),' minutes'])

end
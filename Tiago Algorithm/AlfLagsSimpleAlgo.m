clearvars -except session tSelectMethod intCond stdv deformation spdSmt p1 p2 machine  delGroup intDel


procTime = tic;


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
% 
% Names = Names(intSess,1);


%for s = 1:length(Names(:,1))


Name =  Names{session,1};
disp(['Processing session: ',  Name])




[Speeds,pMats,Distances,Times,trialInfo,trialEvents,paramInfo] = PredSpdSelectDataLagsTrials(Name, intCond, delGroup,tSelectMethod, spdSmt, stdv, p1, p2, machinePath);



clearvars -except pMats trialInfo trialEvents Speeds Distances Times stdv deformation  n_resamples intDel

%%% Predicts Speed from neural activity


%%% Set random number seed (for repeatability)

nSeed=1; %
rng(nSeed,'twister');



% Algorithm parameters

%n_resamples = 100; % nr of folds in which to perform the cross validation process to determine the regularization parameters
nInFolds = 3;


%nParams = [100 1000];
%lParams = [5000 10000];

nParams = [1 10 100 1000 10000 100000 ];
lParams = [5000:1000:10000 15000:5000:100000 200000:100000:1000000];



n_iter = 30; % nr of iterations to determine best lags and coefficients
d = intDel;

rsp = 1;

%%%%%%%%%%%%%%%% PREDICTION %%%%%%%%%%%%%%%%




if stdv >= 0.1;
    sampling = 0.1 / 5;
else
    sampling = stdv / 5;
end

tlim = round([Times{1,d}{3,1}(1) Times{1,d}{3,1}(end)] - Times{1,d}{3,1}(1),4);
smplSize = length(tlim(1):sampling:tlim(2));




[kernel,diffKernel] = GenerateKernelNew(sampling,stdv);

% Create merged data matrices


[tSpikes,tSpds,tTrialInfo,tEventTrials,tDsts,tTms,tTmsSpk] = CreateMergedMatrices(pMats{1,d},trialInfo{1,d},trialEvents{1,d},Speeds{1,d},Distances{1,d},Times{1,d},stdv,sampling);

% outfoldtraining data
trPMats = tSpikes;
trSpeeds = tSpds;


% further partition training set in n folds for cross validation of the ridge parameters

ICV = cvpartition(length(trSpeeds(:,1))','KFold',nInFolds);

for icv = 1:ICV.NumTestSets% for each validation set
    
    % get map of training and test trials
    icvIdx = ICV.training(icv);
    valIdx = ICV.test(icv);
    
    
    % infold training data
    icvPMats = trPMats(icvIdx,:);
    icvSpeeds = cell2mat(trSpeeds(icvIdx,1));
    
    % validation data
    valPMats = trPMats(valIdx,:);
    valSpeeds = cell2mat(trSpeeds(valIdx,1));
    
    
    [icvSmtSpikes, ~] = SmoothSpikeVecsNew(icvPMats,stdv,kernel,diffKernel,tTms,tTmsSpk,sampling);
    [valSmtSpikes, ~] = SmoothSpikeVecsNew(valPMats,stdv,kernel,diffKernel,tTms,tTmsSpk,sampling);
    
    
    % Run a grid search for the ridge parameteres in each fold
    for lN = 1:length(nParams)
        lags = zeros(smplSize,1);
        icvLagSpikes = icvSmtSpikes;
        
        % Fit ridge regression with specified lambda value
        itCoefs{rsp,icv}{1,lN} = ridge(icvSpeeds,icvLagSpikes,nParams(lN),0);
        
        itCoefs{rsp,icv}{2,lN} = ridge(icvSpeeds,icvLagSpikes,nParams(lN),1);
        
        valLagSpikes = valSmtSpikes;
        
        [RMSEVal{rsp}(lN,icv),~,~,~,~,~] = TotalErrorCalculatorNew(valSpeeds,valLagSpikes,itCoefs{rsp,icv}{1,lN},itCoefs{rsp,icv}{2,lN},nParams(lN),lags,1);
        
    end
    
    disp(['Finished inner fold ',num2str(icv),' in ',num2str(nInFolds),' !!!'])
    
    allIcvIdx{rsp,icv} = icvIdx;
    allValIdx{rsp,icv} = valIdx;
    
end


disp('Starting algorithm in entire training set with validated reg params')
% Select best params combination across validation sets

meanE = mean(RMSEVal{rsp},2);
[~,lNIdx] = min(meanE);

trSpeeds = cell2mat(trSpeeds);

[trSmtSpikes, trSmtSpikesDiff] = SmoothSpikeVecsNew(trPMats,stdv,kernel,diffKernel,tTms,tTmsSpk,sampling);

for lL = 1:length(lParams)
    lags = zeros(smplSize,1);
    trLagSpikes = trSmtSpikes;
    
    for iter = 1:n_iter
        
        allLagsTr{rsp}{lL}{iter} = lags;
        
        % Fit ridge regression with specified lambda value
        itCoefsTr{rsp}{lL}{1,iter} = ridge(trSpeeds,trLagSpikes,nParams(lNIdx),0);
        itCoefsTr{rsp}{lL}{2,iter} = ridge(trSpeeds,trLagSpikes,nParams(lNIdx),1);
        
        if iter > 1
            betasDiffTr{rsp}{lL}(1,iter-1) = sqrt(sum((itCoefsTr{rsp}{lL}{1,iter}(2:end) - itCoefsTr{rsp}{lL}{1,iter-1}(2:end)).^2) / length(itCoefsTr{rsp}{lL}{1,iter-1}(2:end)));
            betasDiffTr{rsp}{lL}(2,iter-1) = sqrt(sum((itCoefsTr{rsp}{lL}{2,iter} - itCoefsTr{rsp}{lL}{2,iter-1}).^2) / length(itCoefsTr{rsp}{lL}{2,iter-1}));
        else
            betasDiffTr{rsp}{lL}(1,iter) = 0;
            betasDiffTr{rsp}{lL}(2,iter) = 0;
        end
        
        
        %%% Calculate Error %%%
        
        [RMSETr{rsp}{lL}(1,iter),totErrorTr{rsp}{lL}(1,iter),predSqrdErrorTr{rsp}{lL}(1,iter),lagsCompTr{rsp}{lL}(1,iter),betasCompTr{rsp}{lL}(1,iter),PredsTr{rsp}{lL}{1,iter}] = TotalErrorCalculatorNew(trSpeeds,trLagSpikes,itCoefsTr{rsp}{lL}{1,iter},itCoefsTr{rsp}{lL}{2,iter},nParams(lNIdx),lags,lParams(lL));
        
        
        [lags] = ComputeLagsNew(trSpeeds,trSmtSpikes,trSmtSpikesDiff,itCoefsTr{rsp}{lL}{1,iter},lParams(lL),smplSize);
        
        
        
        [trLagSpikes] = LagsTransformation(trSmtSpikes,trSmtSpikesDiff,lags);
        
        [RMSETr{rsp}{lL}(2,iter),totErrorTr{rsp}{lL}(2,iter),predSqrdErrorTr{rsp}{lL}(2,iter),lagsCompTr{rsp}{lL}(2,iter),betasCompTr{rsp}{lL}(2,iter),PredsTr{rsp}{lL}{2,iter}] = TotalErrorCalculatorNew(trSpeeds,trLagSpikes,itCoefsTr{rsp}{lL}{1,iter},itCoefsTr{rsp}{lL}{2,iter},nParams(lNIdx),lags,lParams(lL));
        
        if sum(abs(lags) > stdv*deformation) > 0
            break
        end
        
    end
end


for reg = 1:length(RMSETr{rsp}(1,:))
    if length(RMSETr{rsp}{reg}(1,:)) == 1
        regPerfs(1,reg) = RMSETr{rsp}{reg}(1,1);
    else
        regPerfs(1,reg) = RMSETr{rsp}{reg}(2,end);
    end
end

[~,trLLIdx] = min(regPerfs);

% store selected lags

finalLags{rsp} =  allLagsTr{rsp}{trLLIdx}{end};

finalPredsTr{rsp,1} = PredsTr{rsp}{trLLIdx}{1,end};
finalPredsTr{rsp,2} = PredsTr{rsp}{trLLIdx}{2,end};

allRegParams{rsp} = [nParams(lNIdx) lParams(trLLIdx) ; lNIdx  trLLIdx];
    
   
% Get data from Alf Algo
predsNoLagOrig = PredsTr{1,1}{1,1}{1,1};
predsLagOrig = finalPredsTr{rsp,2};


% Get data from Tiago Algo
data = load([machinePath,'Cerebro/Recording Data/',Name,'/Predict from Neurons/LagsAlgoResults Delay_',num2str(delGroup),'_',num2str(intDel),' Stdv_',num2str(stdv),'.mat']);
sup_data = load([machinePath,'Cerebro/Recording Data/',Name,'/Predict from Neurons/LagsAlgoResults_Sup Delay_',num2str(delGroup),'_',num2str(intDel),' Stdv_',num2str(stdv),'.mat']);

coefs = data.coefs;
lags = data.lags;
neurons = data.neurons;
speeds = data.speeds;
model_perf = sup_data.mean_fitness;



for t = 1:35

allTrialsNoLag{t,1} =  reshape(neurons(1,:,t,:),211,96);
allTrialsLag{t,1} =  reshape(neurons(6,:,t,:),211,96);

end

NoLagN = cell2mat(allTrialsNoLag);
LagN = cell2mat(allTrialsLag);

betasNoLag = coefs(1,:)';
predsNoLag = NoLagN * betasNoLag;

betasLag = coefs(6,:)';
predsLag = LagN * betasLag;




% Mean lagged an non lagged predictions fot both algorithms

mSpd = mean(reshape(speeds,211,35)');
mPredO = mean(reshape(predsNoLagOrig,211,35)');
mPredOLag = mean(reshape(predsLagOrig,211,35)');
mPred = mean(reshape(predsNoLag,211,35)');
mPredLag = mean(reshape(predsLag,211,35)');

xLinesA = [15.5, 36
           36, 43.3;
           43.3, 75; 
           75, 122; 
           122 176];

clrAMap = ['r','b','r','b','r'];       

xLinesT = [13, 54.8;
           54.8, 103;
           154 190.5];
       
clrTMap = ['r','b','r'];              

subplot(2,2,1)
plot(mSpd);
hold on
plot(mPredO);
hold on
plot(mPredOLag);
title('Alf Lag')

for l = 1:length(xLinesA(:,1))
    patch([xLinesA(l,1) xLinesA(l,2) xLinesA(l,2) xLinesA(l,1)],[-10 -10 60 60],clrAMap(l),'EdgeColor','none','FaceAlpha',.05)
    hold on
end
title('old Algorithm')
ylabel('Speed cm/s')
xlabel('Samples (0.02 s)')


subplot(2,1,1)
plot(mSpd);
hold on
plot(mPred);
hold on
plot(mPredLag);


for l = 1:length(xLinesT(:,1))
    patch([xLinesT(l,1) xLinesT(l,2) xLinesT(l,2) xLinesT(l,1)],[-10 -10 60 60],clrTMap(l),'EdgeColor','none','FaceAlpha',.05)
    hold on
end
ylabel('Speed cm/s')
xlabel('Samples (0.02 s)')
title('New Algorithm')


subplot(2,2,3)
plot(finalLags{1,1}*-1)
line([0 198],[0 0],'Color','r')

for l = 1:length(xLinesA(:,1))
    patch([xLinesA(l,1) xLinesA(l,2) xLinesA(l,2) xLinesA(l,1)],[-0.15 -0.15 0.15 0.15],clrAMap(l),'EdgeColor','none','FaceAlpha',.05)
    hold on
end

ylabel('lag s')
xlabel('Samples (0.02 s)')

subplot(2,1,2)
plot((lags(6,:)*-1)/1000)
line([0 211],[0 0],'Color','r')

for l = 1:length(xLinesT(:,1))
    patch([xLinesT(l,1) xLinesT(l,2) xLinesT(l,2) xLinesT(l,1)],[-0.15 -0.15 0.15 0.15],clrTMap(l),'EdgeColor','none','FaceAlpha',.05)
    hold on
end

ylabel('lag s')
xlabel('Samples (0.02 s)')





% Plot evolution of predictions and lags Tiago algo



subplot(2,1,1)

for it = 1:3


for t = 1:35
allTrialsLag{t,1} =  reshape(neurons(it,:,t,:),211,96);

end

LagN = cell2mat(allTrialsLag);

betasLag = coefs(it,:)';
predsLag = LagN * betasLag;





plot(mean(reshape(predsLag,198,35)'));
hold on

end

% Mean lagged an non lagged predictions fot both algorithms











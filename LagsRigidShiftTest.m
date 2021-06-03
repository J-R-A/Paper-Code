function LagsRigidShiftTest(session, intCond, stdv, deformation, intDel, machine)

clearvars -except session intCond stdv deformation  machine intDel


spdSmt = 500;
p1 = 0.25;
p2 = 0.5;
tSelectMethod = 1;
delGroup = 1;

if stdv >= 0.1;
    sampling = 0.1 / 5;
else
    sampling = stdv / 5;
end


if machine == 1
    machinePath = '/Volumes/';
elseif machine == 2
    machinePath = '/nfs/tank/renart/users/joao.afonso/';
end

load([machinePath,'Cerebro/Recording Data/SessionsMap.mat'])


sess = 1;
for m = 1:length(sessionsMap(:,1))
    for s1 = 1:length(sessionsMap{m,2}(1,:))
        Names{sess,1} = sessionsMap{m,2}{1,s1};
        sess = sess+1;
    end
end


Name =  Names{session,1};

AlfonsoRenart1990

%%% Get data from each session selected trials


[Speeds,pMats,Distances,Times,trialInfo,trialEvents,~] = PredSpdSelectDataLagsTrials(Name, intCond, delGroup,tSelectMethod, spdSmt, stdv, p1, p2, machinePath);



disp(['Processing session: ',  Name, ' stdv ',num2str(stdv),' Deformation ',num2str(deformation),' group ', num2str(intDel)])


% Create merged data matrices

[tSpikes,tSpds,~,~,~,tTms,tTmsSpk] = CreateMergedMatrices(pMats{1,intDel},trialInfo{1,intDel},trialEvents{1,intDel},Speeds{1,intDel},Distances{1,intDel},Times{1,intDel},stdv,sampling);


load([machinePath,'Cerebro/Recording Data/',Name,'/Predict from Neurons/SpdPredLagsBootCheck Trials tsMethod_',num2str(tSelectMethod),' Cond_',num2str(intCond),' Delay_',num2str(delGroup),'_',num2str(intDel),' P1_',num2str(p1),' P2_',num2str(p2),' Sigma_',num2str(stdv),' Deform_',num2str(deformation),'.mat'],'testData');

trTrials = testData.testTrials{1,1};
teTrials = testData.testTrials{1,2};
trPreds = testData.PredsTr;
tePreds = testData.Preds;
lags = testData.finalLags;
betas = testData.betas;
regParams = testData.regParams;



% length(intLags(1,:))
for rs = 1:length(lags(1,:))
    
    
    if mod(rs,10) == 0
        disp(['Processing sample ',num2str(rs),' of ',num2str(length(lags(1,:)))])
    end
    
    rsTrialsTr = trTrials{rs};
    rsTrialsTe = teTrials{rs};
    Lags(:,rs) = lags{rs};
    
    [matSpikesLagTr] = SmoothRigid(tSpikes(rsTrialsTr,:),stdv,tTms(rsTrialsTr,1),tTmsSpk(rsTrialsTr,:),sampling,lags{rs});
    [matSpikesNoLagTr] = SmoothRigid(tSpikes(rsTrialsTr,:),stdv,tTms(rsTrialsTr,1),tTmsSpk(rsTrialsTr,:),sampling,zeros(length(lags{rs}),1));
    
    [matSpikesLagTe] = SmoothRigid(tSpikes(rsTrialsTe,:),stdv,tTms(rsTrialsTe,1),tTmsSpk(rsTrialsTe,:),sampling,lags{rs});
    [matSpikesNoLagTe] = SmoothRigid(tSpikes(rsTrialsTe,:),stdv,tTms(rsTrialsTe,1),tTmsSpk(rsTrialsTe,:),sampling,zeros(length(lags{rs}),1));
    
    
    
    lIdx = regParams{rs}(2,2);
    fBetasLag = betas{rs}{lIdx}{1,end};
    fBetasNoLag = betas{rs}{lIdx}{1,1};
    
    PredsLagTr = matSpikesLagTr * fBetasLag(2:end) + fBetasLag(1);
    PredsNoLagTr = matSpikesNoLagTr * fBetasNoLag(2:end) + fBetasNoLag(1);
    
    PredsLagTe = matSpikesLagTe * fBetasLag(2:end) + fBetasLag(1);
    PredsNoLagTe = matSpikesNoLagTe * fBetasNoLag(2:end) + fBetasNoLag(1);
    
    origPredsLagTr = trPreds{rs,2};
    origPredsLagTe = tePreds{2}{rs};
    origPredsNoLagTe = tePreds{1}{rs};
    
    trSpds = cell2mat(tSpds(rsTrialsTr));
    teSpds = cell2mat(tSpds(rsTrialsTe));
    
    
    
    if rs <= 10
        
        exBetasLag{rs} = fBetasLag;
        exBetasNoLag{rs} = fBetasNoLag;
        
        exPredsLagTr(:,rs) = PredsLagTr;
        exPredsNoLagTr(:,rs) = PredsNoLagTr;
        exPredsLagTe(:,rs) = PredsLagTe;
        exPredsNoLagTe(:,rs) = PredsNoLagTe;
        
        exOrigPredsLagTr(:,rs) = origPredsLagTr;
        exOrigPredsLagTe(:,rs) = origPredsLagTe;
        exOrigPredsNoLagTe(:,rs) = origPredsNoLagTe;
        
        exSpdsTr(:,rs) = trSpds;
        exSpdsTe(:,rs) = teSpds;
        
    end
    
    RMSETr(rs,1) = sqrt(sum((trSpds - PredsNoLagTr).^2) / length(trSpds(:,1)));
    RMSETr(rs,2) = sqrt(sum((trSpds - PredsLagTr).^2) / length(trSpds(:,1)));
    
    RMSETr(rs,3) = sqrt(sum((trSpds - origPredsLagTr).^2) / length(trSpds(:,1)));
   
    RMSETe(rs,1) = sqrt(sum((teSpds - PredsNoLagTe).^2) / length(teSpds(:,1)));
    RMSETe(rs,2) = sqrt(sum((teSpds - PredsLagTe).^2) / length(teSpds(:,1)));
    
    RMSETe(rs,3) = sqrt(sum((teSpds - origPredsNoLagTe).^2) / length(teSpds(:,1)));
    RMSETe(rs,4) = sqrt(sum((teSpds - origPredsLagTe).^2) / length(teSpds(:,1)));
    
end
    

RMSE.Tr = RMSETr;
RMSE.Te = RMSETe;

Betas.Lag = exBetasLag;
Betas.NoLag = exBetasNoLag;

Preds.LagTr = exPredsLagTr;
Preds.NoLagTr = exPredsNoLagTr;
Preds.LagTe = exPredsLagTe;
Preds.NoLagTe = exPredsNoLagTe;

Preds.Orig.LagTr =  exOrigPredsLagTr;
Preds.Orig.LagTe = exOrigPredsLagTe;
Preds.Orig.NoLagTe = exOrigPredsNoLagTe;

Spds.Tr = exSpdsTr;
Spds.Te = exSpdsTe;

           
disp('Saving Results...')

folderPath=[machinePath,'Cerebro/Recording Data/',Name,'/Predict from Neurons/'];

filename_1 = [folderPath,'RigidShiftResults Trials tsMethod_',num2str(tSelectMethod),' Cond_',num2str(intCond),' Delay_',num2str(delGroup),'_',num2str(intDel),' P1_',num2str(p1),' P2_',num2str(p2),' Sigma_',num2str(stdv),' Deform_',num2str(deformation),'.mat'];

save(filename_1, 'RMSE','Betas', 'Preds','Spds','Lags','-v7.3')








clearvars -except Name intCond  spdSmt  p1 p2 machinePath stdv sampling tSelectMethod delGroup

load([machinePath,'Cerebro/Recording Data/',Name,'/OrganizedSessionInfo.mat'],'validSession','organizedTrials','timestampsPerTrial')
load([machinePath,'Cerebro/Recording Data/',Name,'/Neurons.mat']);           % load spike times
if tSelectMethod == 0
    ref = 1;
    load([machinePath,'Cerebro/Recording Data/',Name,'/Predict from Neurons/SimilarDelayTrials ref_',num2str(ref),'.mat'],'delayTrials')
elseif tSelectMethod == 1
    load([machinePath,'Cerebro/Recording Data/',Name,'/Predict from Neurons/SimilarDelayTrialsAuto Group_',num2str(delGroup), ' Method_',num2str(tSelectMethod),'.mat'],'delayTrials')
elseif tSelectMethod == 2
    load([machinePath,'Cerebro/Recording Data/',Name,'/Predict from Neurons/SimilarDelayTrialsAuto Group_',num2str(delGroup), ' Method_',num2str(tSelectMethod),'.mat'],'delayTrials')
end


mouseNumber=str2num(Name(3:4));  % define rewarded and and non rewarded sound based on animal number
rwddSound=mod(mouseNumber,2)+1;
nRwddSound=2-mod(mouseNumber,2);

conditions=[rwddSound 1; rwddSound 0 ;nRwddSound 1;nRwddSound 2]; % Behavioral conditions to analyse

delays=unique(cell2mat(validSession(9,:)));% Sound Starts
delays=delays(delays>0);

trialConditions=cell2mat(validSession(7:9,:)); % Type of trial and outcome for all session trials
unknownDelays = find(trialConditions(3,:) == 0);
timestampsPerTrial(unknownDelays,:) = [];
trialConditions(:,unknownDelays)=[];
organizedTrials(:,unknownDelays)=[];
neurons = Neurons(:,2); % spiketimes


% Categorize trials in 1: hit ; 2: miss 3: correct rejection; 4: false alarms
for t = 1:length(trialConditions(1,:))
    if trialConditions(1,t) == conditions(1,1) && trialConditions(2,t) == conditions(1,2)
        trialTypes(1,t) = 1;
    elseif trialConditions(1,t) == conditions(2,1) && trialConditions(2,t) == conditions(2,2)
        trialTypes(1,t) = 2;
    elseif trialConditions(1,t) == conditions(3,1) && trialConditions(2,t) == conditions(3,2)
        trialTypes(1,t) = 3;
    elseif trialConditions(1,t) == conditions(4,1) && trialConditions(2,t) == conditions(4,2)
        trialTypes(1,t) = 4;
    end
    delayPos(1,t) = find(delays == trialConditions(3,t));
    trialsInfo(t,:) = [trialConditions(1,t) t trialConditions(3,t)];
end


organizedTrials = organizedTrials(:,trialTypes == intCond);
delayPos = delayPos(trialTypes == intCond);
trialsInfo = trialsInfo(trialTypes == intCond,:);
timestampsPerTrial = timestampsPerTrial(trialTypes == intCond,:);
trialTypes = trialTypes(trialTypes == intCond);



% Smooth speeds and dists
for t = 1:length(organizedTrials(1,:))
    spn = spdSmt/length(organizedTrials{1,t});
    organizedTrials{1,t} = smooth(organizedTrials{1,t},spn,'lowess');
    organizedTrials{1,t}(organizedTrials{1,t} < 0) = 0;
    organizedTrials{3,t} = smooth(organizedTrials{3,t},spn,'lowess');
end



% intDels = 11;
% dropTrials{1} = [8,11,12];


for d = 1:length(delayTrials(:,1))
   
    intTrials = organizedTrials(:,ismember(delayPos, delayTrials{d,1}));
    intTInfo{1,d} = trialsInfo(ismember(delayPos, delayTrials{d,1}),:);
    intEventsT = timestampsPerTrial(ismember(delayPos, delayTrials{d,1}),:);
    intTrials(:,delayTrials{d,2}) = [];
    intTInfo{1,d}(delayTrials{d,2},:) = [];
    intEventsT(delayTrials{d,2},:) = [];
    
    [spikesFr{1,d}, speeds{1,d}, dists{1,d}, trialTs{1,d},intEvents{1,d}] = TrialSelectorLags(intTrials,intEventsT,neurons,stdv,p1,p2);
end


for d = 1:length(intTInfo(1,:))
    for t = 1:length(intTInfo{1,d}(:,1))
        cellInfo{1,d}{t,1} = intTInfo{1,d}(t,:); 
        cellEvents{1,d}{t,1} = intEvents{1,d}(t,:); 
        
    end
end

pMats = spikesFr;
Speeds = speeds;
Distances = dists;
Times = trialTs;
trialInfo = cellInfo;
trialEvents = cellEvents;


for d = 1:length(pMats(1,:))
    a = 1;
    
    for t = 1:length(pMats{1,d}(:,1))
         times = Times{1,d}{t,1};
         allTimes{1,d}{t,1} = times - times(1,1);
        for n = 1:length(pMats{1,d}(1,:))
            
            spikes = pMats{1,d}{t,n};
%             spikes_d = zeros(length(pMats{1,d}{t,n}),1);
            %spikesFilter = spikes>= times(1,1) & spikes<=times(end,1);
%             spikes_d(spikesFilter) = 1;
            
            %allSpikes{1,d}{a,1} = (spikes(spikesFilter) - times(1,1))';
            allSpikes{1,d}{a,1} = spikes' - (times(1,1)-3.5*stdv);
            allNeuronsId{1,d}{a,1} = ones(length(allSpikes{1,d}{a}),1) * n ;
            allTrialsId{1,d}{a,1} = ones(length(allSpikes{1,d}{a}),1) * t ;
            a = a+1;
        end
    end
    
    mergedSpikes = cell2mat(allSpikes{1,d});
    mergedNeuronsId = cell2mat(allNeuronsId{1,d});
    mergedTrialsId = cell2mat(allTrialsId{1,d});
    mergedTimes = cell2mat(allTimes{1,d});
    mergedSpeeds = cell2mat(Speeds{1,d});
    
     disp(['Saving data for session ',Name,' trials group ',num2str(d)])
     save([machinePath,'Cerebro/Recording Data/',Name,'/Predict from Neurons/LagsAlgoTrials Delay_',num2str(delGroup),'_',num2str(d),' Stdv_',num2str(stdv),'.mat'],'mergedSpikes','mergedNeuronsId','mergedTrialsId','mergedTimes','mergedSpeeds')
    
end




% t = 1;
% d = 4;
% 
% spd = speeds{1,d}{t,1};
% tms = Times{1,d}{t,1} - Times{1,d}{t,1}(1,1);
% %neur = 
% 
% ev = trialEvents{1,d}{t,1};
% 
% [~,sSidx] = min(abs(tms - ev(2)));
% [~,aSidx] = min(abs(tms - ev(4)));
% mSidx = find(spd <= 5,1,'last');
% 
% subplot(1,2,1)
% plot(tms,spd)
% hold on
% plot(tms([sSidx,aSidx,mSidx]),spd([sSidx,aSidx,mSidx]),'r*')
% 
% subplot(1,2,2)
% plot(intTrials{2,t}(1:end-1)-intTrials{2,t}(1,1), intTrials{1,t})



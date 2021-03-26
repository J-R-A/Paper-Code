function [SpeedMat,FeatMat,tInfoMat,tEvMat,distsMat,paramInfo] = MakeFeaturesMatrixSpdPredictPaper(Name,event,intCond,stdv,spdSmt,machinePath)

clearvars -except Name event intCond stdv  machinePath  spdSmt gauss


% Load beahvioral session info
load([machinePath,'Cerebro/Recording Data/',Name,'/OrganizedSessionInfo.mat'],'validSession','organizedTrials','timestampsPerTrial')
% Load Neurons info
load([machinePath,'Cerebro/Recording Data/',Name,'/Neurons.mat']);           % load spike times


 % Define rewarded and and non rewarded sound based on animal number
mouseNumber=str2num(Name(3:4)); 
rwddSound=mod(mouseNumber,2)+1;
nRwddSound=2-mod(mouseNumber,2);


 % Map of behavioral conditions to analyse
conditions=[rwddSound 1; rwddSound 0 ;nRwddSound 1;nRwddSound 2];

% Retrieve valid sound start locations
delays=unique(cell2mat(validSession(9,:)));
delays=delays(delays>0);

% Trial type, trial outcome and Sound start location of all session trials
trialConditions=cell2mat(validSession(7:9,:)); 

% Get rid of trials with invalid sound start location 0
unknownDelays = find(trialConditions(3,:) == 0);
timestampsPerTrial(unknownDelays,:) = [];
trialConditions(:,unknownDelays)=[];
organizedTrials(:,unknownDelays)=[];


% Spike times of all Neurons
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
    
    % Sound start position indixes of every trial
    sSIdxs(1,t) = find(delays == trialConditions(3,t));
    
    % Trial info matrix
    trialsInfo(t,:) = [trialConditions(1,t) t trialConditions(3,t)];
end


% Smooth speeds and distances
for t = 1:length(organizedTrials(1,:))
    spn = spdSmt/length(organizedTrials{1,t});
    organizedTrials{1,t} = smooth(organizedTrials{1,t},spn,'lowess');
    organizedTrials{3,t} = smooth(organizedTrials{3,t},spn,'lowess');
end


% Select time window, within each trial, from which we want to use data 

timeWindow = [timestampsPerTrial(:,floor(event(1))) timestampsPerTrial(:,floor(event(2)))];

%{
if mod(event(1),1) ~= 0
    
    if event(1) == 3.5
        for t = 1:length(organizedTrials(1,:))
            
            asP = find((organizedTrials{3,t} - organizedTrials{3,t}(1,1)) >= 95,1,'first');
            beforeArea(t,1) = organizedTrials{2,t}(asP);
            
        end
    end
    timeWindow(:,1) = beforeArea;
end


if mod(event(2),1) ~= 0
    if event(2) == 4.5
        for t = 1:length(organizedTrials(1,:))
            
            if length(intCond) == 1
                
                asP = find((organizedTrials{3,t} - organizedTrials{3,t}(1,1)) >= 104,1,'first');
                aeP = find((organizedTrials{3,t} - organizedTrials{3,t}(1,1)) >= 119,1,'first');
                
                stopTime = find(organizedTrials{1,t}(asP : aeP) < 5,1,'first');
                
                
                if isempty(stopTime) && trialTypes(1,t) ~= 1
                    stopTimes(t,1) = timeWindow(t,2);
                else
                    stopTimes(t,1) = organizedTrials{2,t}(asP + stopTime) + 0.2;
                end
                
            elseif length(intCond) ==2
                
                asP = find((organizedTrials{3,t} - organizedTrials{3,t}(1,1)) >= 104,1,'first');
                
                stopTimes(t,1) = organizedTrials{2,t}(asP) + 0.2;
                
            end
            
        end
        timeWindow(:,2) = stopTimes;
    end
    
end
%}



% Calculate sampling of smoothed fr based on the kernel width 
if stdv >= 0.1;
    dt = 0.1 / 5;
else
    dt = stdv / 5;
end





% Prealocate variables
speeds = cell(length(timeWindow(:,1)),1);
dists = cell(length(timeWindow(:,1)),1);
spikesFr = cell(length(neurons(:,1)),length(timeWindow(:,1)));
smtNeurs = cell(length(neurons(:,1)),1);


% Generate gaussian kernel
[kernel,~] = GenerateKernelNew(dt,stdv);


% Smooth neurons by convolving with gaussian kernel 
for n = 1:length(neurons(:,1))
    
    if mod(n,10) == 0
        disp(['Neurons smoothing for session ', Name,' ',num2str(round((n/length(neurons(:,1))*100))),'% completed'])
    end
   
    [smtNeurs{n,1},sessionTS] = GaussianSmoothingPaper(neurons{n,1},stdv,dt,[timeWindow(1,1) timeWindow(end,2)],kernel);
end



% Organize and transform data acording to fr sampling  
for t = 1:length(timeWindow(:,1))
   for n = 1:length(smtNeurs(:,1))
        
        %[smtSpikes,timeStamps] = GaussianSmoothing(neurons{n},stdv,dt,timeWindow(t,:));
        %spikesFr{n,t} = smtSpikes';
        
        spikesFr{n,t} = smtNeurs{n,1}(sessionTS >= timeWindow(t,1) & sessionTS <= timeWindow(t,2))';
        trialTS =  sessionTS(sessionTS >= timeWindow(t,1) & sessionTS <= timeWindow(t,2));
        
        
        if n == 1;
            tempSpds = organizedTrials{1,t};
            tempSpds(end+1) = tempSpds(end);
            tempDst  = organizedTrials{3,t};
            tempTs = organizedTrials{2,t};
            speeds{t,1} = zeros(length(trialTS(1,:)),1);
            dists{t,1} = zeros(length(trialTS(1,:)),1);
            for i = 1:length(trialTS(1,:))
                [~,idx] = min(abs(tempTs-trialTS(i)));
                speeds{t,1}(i,1)= tempSpds(idx);
                dists{t,1}(i,1)= tempDst(idx);
            end
            dists{t,1} = dists{t,1} - tempDst(1);
        end
    end
end


% prealocate cells for later data storage
orgData = cell(4,11);
orgSpeeds = cell(4,11);
orgInfo = cell(4,11);
orgEvents = cell(4,11);
orgDists = cell(4,11);

% organize trials by condition and delay
for t = 1:length(spikesFr(1,:))
    [r,~]=size(orgData{trialTypes(1,t),sSIdxs(1,t)});
    orgData{trialTypes(1,t),sSIdxs(1,t)}(r+1,:) = spikesFr(:,t);
    orgSpeeds{trialTypes(1,t),sSIdxs(1,t)}(r+1,1) = speeds(t,1);
    orgInfo{trialTypes(1,t),sSIdxs(1,t)}{r+1,:} = trialsInfo(t,:);
    orgEvents{trialTypes(1,t),sSIdxs(1,t)}{r+1,:} = timestampsPerTrial(t,:);
    orgDists{trialTypes(1,t),sSIdxs(1,t)}{r+1,:} = dists{t,1};
end



% concatenate trials from the conditions we are interested in: intCond
for d = 1:length(orgData(1,:))
    Features{1,d} = vertcat(orgData{intCond,d});
    Speeds{1,d} = vertcat(orgSpeeds{intCond,d});
    trialInfo{1,d} = vertcat(orgInfo{intCond,d});
    trialEvents{1,d} = vertcat(orgEvents{intCond,d});
    Distances{1,d} = vertcat(orgDists{intCond,d});
end


% Concatenate all delays in a last feature and shuffle to avoid delay speciric effects


trialsNum = sum(cell2mat(cellfun(@(x) length(x(:,1)),Features,'UniformOutput',0)));
trialsOrder = randperm(trialsNum);

Features{1,d+1} = vertcat(Features{:});
Features{1,end} = Features{1,end}(trialsOrder,:); 

Speeds{1,d+1} = vertcat(Speeds{:});
Speeds{1,end} = Speeds{1,end}(trialsOrder,:); 

trialInfo{1,d+1} = vertcat(trialInfo{:});
trialInfo{1,end} = trialInfo{1,end}(trialsOrder,:);

trialEvents{1,d+1} = vertcat(trialEvents{:});
trialEvents{1,end} = trialEvents{1,end}(trialsOrder,:);

Distances{1,d+1} = vertcat(Distances{:});
Distances{1,end} = Distances{1,end}(trialsOrder,:);




for d = 1:length(Features(1,:))
    
    FeatMat{1,d} = cell2mat(Features{1,d});
    SpeedMat{1,d} = cell2mat(Speeds{1,d});
    distsMat{1,d} = cell2mat(Distances{1,d});
    
end





% Create map between each datapoint and that trial Info


for d=1:length(trialInfo(1,:))
        for t = 1:length(trialInfo{1,d}(:,1))
            
            tempTType{1,d}{t,1} = repmat(trialInfo{1,d}{t,1}(1,1),length(Features{1,d}{t,1}(:,1)),1);
            tempTNumber{1,d}{t,1} = repmat(trialInfo{1,d}{t,1}(1,2),length(Features{1,d}{t,1}(:,1)),1);
            tempTDelay{1,d}{t,1} = repmat(trialInfo{1,d}{t,1}(1,3),length(Features{1,d}{t,1}(:,1)),1);
            tempTEvent{1,d}{t,1} = repmat(trialEvents{1,d}{t,1}(1,:),length(Features{1,d}{t,1}(:,1)),1);
            
        end
        
        tInfoMat{1,d}(:,1) = cell2mat(tempTType{1,d});
        tInfoMat{1,d}(:,2) = cell2mat(tempTNumber{1,d});
        tInfoMat{1,d}(:,3) = cell2mat(tempTDelay{1,d});
        tEvMat{1,d} = cell2mat(tempTEvent{1,d});
end







paramInfo.Name = Name;
paramInfo.trialPart = event;
paramInfo.rwddSound = rwddSound;
paramInfo.intCond = intCond;
paramInfo.spdSmt = spdSmt;
paramInfo.GaussWidth = stdv;    



% [~,ia,~]=intersect(origDelays, delays);  % If not using trials for all delays this will place the data
% in the correct columns


end







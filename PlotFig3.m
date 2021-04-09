clearvars -except toPlot
close all


%%% Script to plot the different panels of the paper second figure


%% Load and process needed Data

load('/Volumes/Cerebro/Recording Data/SessionsMap.mat')


if ~exist('toPlot','var')
    
    %%% Load panel A data
    
    
    m = 5;
    s = 1;
    
    
    Name = sessionsMap{m,2}{1,s} % Name of the session
    mouseNumber = str2num(Name(3:4)); % number of the mouse being loaded
    rwddSound=mod(mouseNumber,2)+1; % define rewarded and and non rewarded sound based on animal number
    nRwddSound=2-mod(mouseNumber,2);
    
    
    
    % load raw neural data
    disp('Loading Raw Neurons')
    load(['/Volumes/Cerebro/Recording Data/',Name,'/Neurons.mat'])
    
    
    % load session behavioral data
    disp('Loading Behavior data')
    load(['/Volumes/Cerebro/Recording Data/',Name,'/OrganizedSessionInfo.mat'],'validSession','organizedTrials')
    
    
    
    
    cond1_filter = cell2mat(validSession(7,:)) == rwddSound & cell2mat(validSession(8,:)) == 1;
    cond3_filter = cell2mat(validSession(7,:)) == nRwddSound & cell2mat(validSession(8,:)) == 1;
    
    cond1_trials = organizedTrials(:,cond1_filter);
    cond3_trials = organizedTrials(:,cond3_filter);
    
    for n = 1:length(Neurons(:,1))
        
        for t = 1:length(cond1_trials(1,:))
            
            t_filter = Neurons{n,2} >= cond1_trials{4,t}(1,1) &  Neurons{n,2} <= cond1_trials{4,t}(1,2);
            cond1_spikes{t,n} = Neurons{n,2}(t_filter);
            cond1_nSpikes(t,n) = sum(t_filter);
            cond1_tTimes(t,n) = cond1_trials{4,t}(1,2) - cond1_trials{4,t}(1,1);
        end
        
        for t = 1:length(cond3_trials(1,:))
            
            t_filter = Neurons{n,2} >= cond3_trials{4,t}(1,1) &  Neurons{n,2} <= cond3_trials{4,t}(1,2);
            cond3_spikes{t,n} = Neurons{n,2}(t_filter);
            cond3_nSpikes(t,n) = sum(t_filter);
            cond3_tTimes(t,n) = cond3_trials{4,t}(1,2) - cond3_trials{4,t}(1,1);
        end
    end
    
    toPlot.A{1,1}{1,1} = cond1_spikes;
    toPlot.A{1,1}{1,2} = cond1_nSpikes;
    toPlot.A{1,1}{1,3} = cond1_tTimes;
    toPlot.A{1,1}{1,4} = cond1_trials;
    
    toPlot.A{1,2}{1,1} = cond3_spikes;
    toPlot.A{1,2}{1,2} = cond3_nSpikes;
    toPlot.A{1,2}{1,3} = cond3_tTimes;
    toPlot.A{1,2}{1,4} = cond3_trials;
    
    clearvars -except toPlot sessionsMap
    
    
    
    
    
    
    %%% Load panel B data
    
    intCond = [1 3];
    intCondO = [1 3];
    event=[1 6];
    eventO = [2 4];
    stdv = 0.14;
    alpha = 0.5;
    spdSmt = 500;
    
    
    m = 5;
    % :length(intSess{m,2}(1,:))
    ss = 1;
    
    
    
    Name = sessionsMap{m,2}{1,ss};
    mouseNumber = str2num(Name(3:4)); % number of the mouse being loaded
    rwddSound=mod(mouseNumber,2)+1; % define rewarded and and non rewarded sound based on animal number
    nRwddSound=2-mod(mouseNumber,2);
    sID = [rwddSound nRwddSound];
    sID(sID==nRwddSound) = 0;
    sID(sID == rwddSound) = 1;
    
    disp(['Processing session: ',num2str(ss),' of mice ',num2str(Name(1:4))])
    
    disp('Loading data')
    
    
    for c = 1:length(intCond(1,:))
        load(['/Volumes/Cerebro/Recording Data/',Name,'/Predict from Neurons/SpdPredResultsTimePaper Cond_',num2str(intCond(1,c)),' TSeg_',num2str(event),' spdSmt_',num2str(spdSmt),' Sigma_', num2str(stdv),' Alpha_',num2str(alpha),'.mat'],...
            'allPredMat','Laps','dummyLaps','allTrials','coefs','allDists');
        
        for d = 1:length(allPredMat(1,:))
            
            LapsS{c,d} = Laps{1,d};
            coefsS{c,d} = coefs{1,d};
            %dummyLapsS{c,1} = dummyLaps{1,12};
            dummyLapsS{c,d} = (allPredMat{1,d} * coefsS{c,d}(2:end)) + coefsS{c,d}(1);
            outTrialsS{c,d} = allTrials{1,d};
            outDistsS{c,d} = allDists{1,d};
        end
        
        clear Laps dummyLaps allTrials coefs allPredMat
        
    end
    
    
    
    load(['/Volumes/Cerebro/Recording Data/',Name,'/Predict from Neurons/OutFeaturesTimePaper Cond_',num2str(intCondO),' TSeg_',num2str(event),' spdSmt_',num2str(spdSmt),' Sigma_', num2str(stdv),' Alpha_',num2str(alpha),'.mat'],...
        'allPredMat','allTrials','allDists');
    
    load(['/Volumes/Cerebro/Recording Data/',Name,'/Predict from Neurons/OutPredResultsTimePaper Cond_',num2str(intCondO),' TSeg_',num2str(eventO),' spdSmt_',num2str(spdSmt),' Sigma_', num2str(stdv),' Alpha_',num2str(alpha),'.mat'],...
        'coefs','Laps','dummyLaps');
    
    
    for d = 1:length(allPredMat(1,:))
        
        pMatsO = allPredMat{1,d};
        coefsO{1,d} = coefs{1,d};
        
        pp = exp(-(coefsO{1,d}(1) + pMatsO * coefsO{1,d}(2:end)));
        dummyLapsOT = 1 ./ (1 + pp);
          
        for c = 1:length(intCondO(1,:))
            tType = sID(c);
            
            %         if sID(c) == 2
            %             tType = 0;
            %         end
            
            
            dummyLapsO{c,d} = dummyLapsOT(allTrials{1,d}(:,1) == tType,:);
            outTrialsO{c,d} = allTrials{1,d}(allTrials{1,d}(:,1) == tType,:);
            LapsO{c,d} = allTrials{1,d}(allTrials{1,d}(:,1) == tType,1);
            outDistsO{c,d} = allDists{1,d}(allTrials{1,d}(:,1) == tType,1);
        end
        
    end
    
    
    clear Laps dummyLaps outTrials coefs allPredMat
    
    
    load(['/Volumes/Cerebro/Recording Data/',Name,'/OrganizedSessionInfo.mat'],'validSession','organizedTrials','timestampsPerTrial')
    
    %load(['/Volumes/Cerebro/Recording Data/',Name,'/Predict from Neurons/Neur-OutPredTime Cond_',num2str(intCondO),' TSeg_',num2str(eventO),' KWdth_',num2str(kWdthS),' BinSz_',num2str(binSize),' SB_',num2str(SB),'.mat'],'Distances')
    
    % Separate speed values by trial and delay type
    
    
    
    for d = 1:length(LapsS(1,:))
        for c = 1:length(LapsS(:,1));
            
            
            %outTrialsS{c,d} = cell2mat(outTrialsS{c,d});
            
            trialDelaysS = unique(outTrialsS{c,d}(:,3));
            trialTypesS = unique(outTrialsS{c,d}(:,1));
            trialNumsS = unique(outTrialsS{c,d}(:,2));
            
            trialDelaysO = unique(outTrialsO{c,d}(:,3));
            trialTypesO = unique(outTrialsO{c,d}(:,1));
            trialNumsO = unique(outTrialsO{c,d}(:,2));
            
            orgSpeeds = cell(length(trialTypesS(:,1)),length(trialDelaysS(:,1)));
            orgDummySpeeds = cell(length(trialTypesS(:,1)),length(trialDelaysS(:,1)));
            orgTrialsSpeeds = cell(length(trialTypesS(:,1)),length(trialDelaysS(:,1)));
            orgDistsSpds = cell(length(trialTypesO(:,1)),length(trialDelaysO(:,1)));
            
            orgOuts = cell(length(trialTypesO(:,1)),length(trialDelaysO(:,1)));
            orgDummyOuts = cell(length(trialTypesO(:,1)),length(trialDelaysO(:,1)));
            orgTrialsOuts = cell(length(trialTypesO(:,1)),length(trialDelaysO(:,1)));
            orgDistsOuts = cell(length(trialTypesO(:,1)),length(trialDelaysO(:,1)));
            
            
            for tn = 1:length(trialNumsS(:,1))
                
                speeds = { LapsS{c,d}(outTrialsS{c,d}(:,2) == trialNumsS(tn,1),1)};
                dSpeeds ={ dummyLapsS{c,d}(outTrialsS{c,d}(:,2) == trialNumsS(tn,1),1)};
                distSpds = { outDistsS{c,d}(outTrialsS{c,d}(:,2) == trialNumsS(tn,1),1)};
                tDelayS = unique(outTrialsS{c,d}(outTrialsS{c,d}(:,2) == trialNumsS(tn,1),3));
                tTypeS = unique(outTrialsS{c,d}(outTrialsS{c,d}(:,2) == trialNumsS(tn,1),1));
                
                Outs = { LapsO{c,d}(outTrialsO{c,d}(:,2) == trialNumsO(tn,1),1)};
                dOuts = { dummyLapsO{c,d}(outTrialsO{c,d}(:,2) == trialNumsO(tn,1),1)};
                distOuts = { outDistsO{c,d}(outTrialsO{c,d}(:,2) == trialNumsO(tn,1),1)};
                tDelayO = unique(outTrialsO{c,d}(outTrialsO{c,d}(:,2) == trialNumsO(tn,1),3));
                tTypeO = unique(outTrialsO{c,d}(outTrialsO{c,d}(:,2) == trialNumsO(tn,1),1));
                
                orgSpeeds{tTypeS == trialTypesS,tDelayS == trialDelaysS} = [orgSpeeds{tTypeS == trialTypesS,tDelayS == trialDelaysS}; speeds];
                orgDummySpeeds{tTypeS == trialTypesS,tDelayS == trialDelaysS} = [orgDummySpeeds{tTypeS == trialTypesS,tDelayS == trialDelaysS}; dSpeeds];
                orgDistsSpds{tTypeS == trialTypesS,tDelayS == trialDelaysS} = [orgDistsSpds{tTypeS == trialTypesS,tDelayS == trialDelaysS}; distSpds];
                orgTrialsSpeeds{tTypeS == trialTypesS,tDelayS == trialDelaysS} = [orgTrialsSpeeds{tTypeS == trialTypesS,tDelayS == trialDelaysS}; trialNumsS(tn,1)];
                
                orgOuts{tTypeO == trialTypesO,tDelayO == trialDelaysO} = [orgOuts{tTypeO == trialTypesO,tDelayO == trialDelaysO}; Outs];
                orgDummyOuts{tTypeO == trialTypesO,tDelayO == trialDelaysO} = [orgDummyOuts{tTypeO == trialTypesO,tDelayO == trialDelaysO}; dOuts];
                orgDistsOuts{tTypeO == trialTypesO,tDelayO == trialDelaysO} = [orgDistsOuts{tTypeO == trialTypesO,tDelayO == trialDelaysO}; distOuts];
                orgTrialsOuts{tTypeO == trialTypesO,tDelayO == trialDelaysO} = [orgTrialsOuts{tTypeO == trialTypesO,tDelayO == trialDelaysO}; trialNumsO(tn,1)];
            end
            
            %         for tt = 1:length(trialTypesS(:,1))
            %             orgSpeeds(tt,length(trialDelaysS)+1) = { LapsS{c,d}(outTrialsS{c,d}(:,1) == trialTypesS(tt,1),1)};
            %             orgDummySpeeds(tt,length(trialDelaysS)+1) = { dummyLapsS{c,d}(outTrialsS{c,d}(:,1) == trialTypesS(tt,1),1)};
            %             orgDistsSpds(tt,length(trialDelaysS)+1) = { outDistsS{c,d}(outTrialsS{c,d}(:,1) == trialTypesS(tt,1),1)};
            %
            %             orgOuts(tt,length(trialDelaysO)+1) = { LapsO{c,d}(outTrialsO{c,d}(:,1) == trialTypesO(tt,1),1)};
            %             orgDummyOuts(tt,length(trialDelaysO)+1) = { dummyLapsO{c,d}(outTrialsO{c,d}(:,1) == trialTypesO(tt,1),1)};
            %             orgDistsOuts(tt,length(trialDelaysO)+1) = { outDistsO{c,d}(outTrialsO{c,d}(:,1) == trialTypesO(tt,1),1)};
            %         end
            
            
            if d < length(LapsS(1,:))
                finalSpeeds(c,d) = orgSpeeds;
                finalDummySpeeds(c,d) = orgDummySpeeds;
                finalDistsSpeeds(c,d) = orgDistsSpds;
                
                finalOuts(c,d) = orgOuts;
                finalDummyOuts(c,d) = orgDummyOuts;
                finalDistsOuts(c,d) = orgDistsOuts;
                
            else
                finalSpeeds{c,d}(1,:) = orgSpeeds(1:end);
                finalDummySpeeds{c,d}(1,:) = orgDummySpeeds(1:end);
                finalDistsSpeeds{c,d}(1,:) = orgDistsSpds(1:end);
                
                finalOuts{c,d}(1,:) = orgOuts(1:end);
                finalDummyOuts{c,d}(1,:) = orgDummyOuts(1:end);
                finalDistsOuts{c,d}(1,:) = orgDistsOuts(1:end);
            end
            
            %         finalDistsSpeeds{c,d} = orgDistsSpds{1:end-1};
            %
            %         finalOuts{c,d} = orgOuts{1:end-1};
            %         finalDummyOuts{c,d} = orgDummyOuts{1:end-1};
            %         finalDistsOuts{c,d} = orgDistsOuts{1:end-1};
            
            
        end
    end

toPlot.B{1,1}{1,1} = finalSpeeds;
toPlot.B{1,1}{1,2} = finalDummySpeeds;
toPlot.B{1,1}{1,3} = finalDistsSpeeds;

toPlot.B{1,2}{1,1} = finalOuts;
toPlot.B{1,2}{1,2} = finalDummyOuts;
toPlot.B{1,2}{1,3} = finalDistsOuts;

clearvars -except toPlot sessionsMap


%%% Load panel C and D data

s  = 1;

for m = 1:5
    % :length(intSess{m,2}(1,:))
    for ss = 1:length(sessionsMap{m,2}(1,:))
        
        intCond = [1 3];
        intCondO = [1 3];
        event=[1 6];
        eventO = [2 4];
        stdv = 0.14 ;
        alpha = 0.5;
        spdSmt = 500;
        
        
        Name = sessionsMap{m,2}{1,ss};
        mouseNumber = str2num(Name(3:4)); % number of the mouse being loaded
        rwddSound=mod(mouseNumber,2)+1; % define rewarded and and non rewarded sound based on animal number
        nRwddSound=2-mod(mouseNumber,2);
        sID = [rwddSound nRwddSound];
        sID(sID == nRwddSound) = 0;
        sID(sID == rwddSound) = 1;
        
        
        
        disp(['Processing session: ',num2str(ss),' of mice ',num2str(Name(1:4))])
        
        disp('Loading data')
        
        
        for c = 1:length(intCond(1,:))
            
            load(['/Volumes/Cerebro/Recording Data/',Name,'/Predict from Neurons/SpdPredResultsTimePaper Cond_',num2str(intCond(1,c)),' TSeg_',num2str(event),' spdSmt_',num2str(spdSmt),' Sigma_', num2str(stdv),' Alpha_',num2str(alpha),'.mat'],...
                'allPredMat','Laps','dummyLaps','allTrials','coefs');
            
            LapsS{c,1} = Laps{1,12};
            coefsS{s,c} = coefs{1,12};
            dummyLapsS{c,1} = dummyLaps{1,12};
            %dummyLapsS{c,1} = (allPredMat{1,12} * coefsS{c,1}(2:end)) + coefsS{c,1}(1);
            outTrialsS{c,1} = allTrials{1,12};
            
            
            clear Laps dummyLaps allTrials coefs allPredMat
            
        end
        
        
         
        
        load(['/Volumes/Cerebro/Recording Data/',Name,'/Predict from Neurons/OutPredResultsTimePaper Cond_',num2str(intCondO),' TSeg_',num2str(eventO),' spdSmt_',num2str(spdSmt),' Sigma_', num2str(stdv),' Alpha_',num2str(alpha),'.mat'],...
            'Laps','dummyLaps','allTrials','coefs','allDists');
        
        for c = 1:length(intCond(1,:))
            
            tType = sID(c);
            
%             if sID(c) == 2
%                 tType = 0;
%             end
        
            
            
            LapsO{c,1} = Laps{1,12}(allTrials{1,12}(:,1) == tType,1);
            
            %     if c == 2
            %         LapsO{c,1} = LapsO{c,1} - 2;
            %     end
            
            %dummyLaps{1,12} = ones(length(dummyLaps{1,12}(:,1)),1) - dummyLaps{1,12}(:,1);
            dummyLapsO{c,1} = dummyLaps{1,12}(allTrials{1,12}(:,1) == tType,1);
            outTrialsO{c,1} = allTrials{1,12}(allTrials{1,12}(:,1) == tType,:);
            outDistsO{c,1} =  allDists{1,12}(allTrials{1,12}(:,1) == tType,1);
        end
        
        coefsO{s,1} = coefs{1,12};
        
        
        clear Laps dummyLaps allTrials coefs allPredMat
        
        
        load(['/Volumes/Cerebro/Recording Data/',Name,'/OrganizedSessionInfo.mat'],'validSession','organizedTrials','timestampsPerTrial')
        
        %load(['/Volumes/Cerebro/Recording Data/',Name,'/Predict from Neurons/Neur-OutPredTime Cond_',num2str(intCondO),' TSeg_',num2str(eventO),' KWdth_',num2str(kWdthS),' BinSz_',num2str(binSize),' SB_',num2str(SB),'.mat'],'Distances')
        
        % Separate speed values by trial and delay type
        
        
        for d = 1
            for c = 1:length(LapsS(:,1));
                
                %outTrialsS{c,d} = cell2mat(outTrialsS{c,d});
              
                trialDelaysS = unique(outTrialsS{c,d}(:,3));
                trialTypesS = unique(outTrialsS{c,d}(:,1));
                trialNumsS = unique(outTrialsS{c,d}(:,2));
                
                trialDelaysO = unique(outTrialsO{c,d}(:,3));
                trialTypesO = unique(outTrialsO{c,d}(:,1));
                trialNumsO = unique(outTrialsO{c,d}(:,2));
                
                orgSpeeds = cell(length(trialTypesS(:,1)),length(trialDelaysS(:,1)));
                orgDummySpeeds = cell(length(trialTypesS(:,1)),length(trialDelaysS(:,1)));
                orgTrialsSpeeds = cell(length(trialTypesS(:,1)),length(trialDelaysS(:,1)));
                 
                orgOuts = cell(length(trialTypesO(:,1)),length(trialDelaysO(:,1)));
                orgDummyOuts = cell(length(trialTypesO(:,1)),length(trialDelaysO(:,1)));
                orgTrialsOuts = cell(length(trialTypesO(:,1)),length(trialDelaysO(:,1)));
                orgDistsOuts = cell(length(trialTypesO(:,1)),length(trialDelaysO(:,1)));
                
                
                for tn = 1:length(trialNumsS(:,1))
                    
                    speeds = { LapsS{c,d}(outTrialsS{c,d}(:,2) == trialNumsS(tn,1),1)};
                    dSpeeds ={ dummyLapsS{c,d}(outTrialsS{c,d}(:,2) == trialNumsS(tn,1),1)};
                    tDelayS = unique(outTrialsS{c,d}(outTrialsS{c,d}(:,2) == trialNumsS(tn,1),3));
                    tTypeS = unique(outTrialsS{c,d}(outTrialsS{c,d}(:,2) == trialNumsS(tn,1),1));
                    
                    Outs = { LapsO{c,d}(outTrialsO{c,d}(:,2) == trialNumsO(tn,1),1)};
                    dOuts = { dummyLapsO{c,d}(outTrialsO{c,d}(:,2) == trialNumsO(tn,1),1)};
                    distOuts = { outDistsO{c,d}(outTrialsO{c,d}(:,2) == trialNumsO(tn,1),1)};
                    tDelayO = unique(outTrialsO{c,d}(outTrialsO{c,d}(:,2) == trialNumsO(tn,1),3));
                    tTypeO = unique(outTrialsO{c,d}(outTrialsO{c,d}(:,2) == trialNumsO(tn,1),1));
                    
                    orgSpeeds{tTypeS == trialTypesS,tDelayS == trialDelaysS} = [orgSpeeds{tTypeS == trialTypesS,tDelayS == trialDelaysS}; speeds];
                    orgDummySpeeds{tTypeS == trialTypesS,tDelayS == trialDelaysS} = [orgDummySpeeds{tTypeS == trialTypesS,tDelayS == trialDelaysS}; dSpeeds];
                    orgTrialsSpeeds{tTypeS == trialTypesS,tDelayS == trialDelaysS} = [orgTrialsSpeeds{tTypeS == trialTypesS,tDelayS == trialDelaysS}; trialNumsS(tn,1)];
                    
                    
                    orgOuts{tTypeO == trialTypesO,tDelayO == trialDelaysO} = [orgOuts{tTypeO == trialTypesO,tDelayO == trialDelaysO}; Outs];
                    orgDummyOuts{tTypeO == trialTypesO,tDelayO == trialDelaysO} = [orgDummyOuts{tTypeO == trialTypesO,tDelayO == trialDelaysO}; dOuts];
                    orgDistsOuts{tTypeO == trialTypesO,tDelayO == trialDelaysO} = [orgDistsOuts{tTypeO == trialTypesO,tDelayO == trialDelaysO}; distOuts];
                    orgTrialsOuts{tTypeO == trialTypesO,tDelayO == trialDelaysO} = [orgTrialsOuts{tTypeO == trialTypesO,tDelayO == trialDelaysO}; trialNumsO(tn,1)];
                end
                
                for tt = 1:length(trialTypesS(:,1))
                    orgSpeeds(tt,length(trialDelaysS)+1) = { LapsS{c,d}(outTrialsS{c,d}(:,1) == trialTypesS(tt,1),1)};
                    orgDummySpeeds(tt,length(trialDelaysS)+1) = { dummyLapsS{c,d}(outTrialsS{c,d}(:,1) == trialTypesS(tt,1),1)};
                    
                    orgOuts(tt,length(trialDelaysO)+1) = { LapsO{c,d}(outTrialsO{c,d}(:,1) == trialTypesO(tt,1),1)};
                    orgDummyOuts(tt,length(trialDelaysO)+1) = { dummyLapsO{c,d}(outTrialsO{c,d}(:,1) == trialTypesO(tt,1),1)};
                    orgDistsOuts(tt,length(trialDelaysO)+1) = { outDistsO{c,d}(outTrialsO{c,d}(:,1) == trialTypesO(tt,1),1)};
                end
                
                finalSpeeds(c,:) = orgSpeeds(1:end);
                finalDummySpeeds(c,:) = orgDummySpeeds(1:end);
                 
                finalOuts(c,:) = orgOuts(1:end);
                finalDummyOuts(c,:) = orgDummyOuts(1:end);
                finalDistsOuts(c,:) = orgDistsOuts(1:end);
                
                
                for dd = 1:length(finalSpeeds(1,:))
                    
                    if dd < length(finalSpeeds(1,:))
                        
                        spd = cell2mat(finalSpeeds{c,dd});
                        dSpd = cell2mat(finalDummySpeeds{c,dd});
                        
                        sRSS{c,dd} = sum((spd -  dSpd).^2);
                        sTSS{c,dd} = sum((spd - mean(spd)).^2);
                        sR2{s}(c,dd) = 1- (sRSS{c,dd} / sTSS{c,dd});
                        
                    else
                       
                        spd = finalSpeeds{c,dd};
                        dSpd = finalDummySpeeds{c,dd};
                        
                        sRSS{c,dd} = sum((spd -  dSpd).^2);
                        sTSS{c,dd} = sum((spd - mean(spd)).^2);
                        sR2{s}(c,dd) = 1- (sRSS{c,dd} / sTSS{c,dd});
                        
                    end
                    
                end
                
            end
            
            for dd = 1:length(finalOuts(1,:))
                thresholds=0:0.001:1;
                
                if dd < length(finalOuts(1,:))
                    
                    outT1=cell2mat(finalOuts{1,dd});
                    outT2=cell2mat(finalOuts{2,dd});
                    out = [outT1;outT2];
                    
                    dOutT1 = cell2mat(finalDummyOuts{1,dd});
                    dOutT2 = cell2mat(finalDummyOuts{2,dd});
                    dOut = [dOutT1;dOutT2];
                    
                    
                    
                    for j=1:length(thresholds(1,:))
                        
                        predStopsALL(j,:)=(sign(dOut-thresholds(j))+1)/2;
                        AllPerfAll{2,dd}(j)=length(find(out'.*predStopsALL(j,:)==1))/length(find(out==1));
                        AllPerfAll{3,dd}(j)=length(find((1-out)'.*predStopsALL(j,:)==1))/length(find(out==0));
                        AllPerfAll{1,dd}= abs(trapz([AllPerfAll{3,dd}(:)],[AllPerfAll{2,dd}(:)]));
                        
                    end
                    
                    perfAll(s,dd)=AllPerfAll{1,dd};
                    
                else
                    
                    outT1=finalOuts{1,dd};
                    outT2=finalOuts{2,dd};
                    out = [outT1;outT2];
                    
                    dOutT1 = finalDummyOuts{1,dd};
                    dOutT2 = finalDummyOuts{2,dd};
                    dOut = [dOutT1;dOutT2];
                    
                    
                    
                    for j=1:length(thresholds(1,:))
                        
                        predStopsALL(j,:)=(sign(dOut-thresholds(j))+1)/2;
                        AllPerfAll{2,dd}(j)=length(find(out'.*predStopsALL(j,:)==1))/length(find(out==1));
                        AllPerfAll{3,dd}(j)=length(find((1-out)'.*predStopsALL(j,:)==1))/length(find(out==0));
                        AllPerfAll{1,dd}= abs(trapz([AllPerfAll{3,dd}(:)],[AllPerfAll{2,dd}(:)]));
                        
                    end
                    
                    perfAll(s,dd)=AllPerfAll{1,dd};
                end
                
                clear predStopsALL AllPerfAll
            end
        end
        
        
        load(['/Volumes/Cerebro/Recording Data/',Name,'/Neurons.mat'])
        nShank{s} = cell2mat(Neurons(:,3));
        
        
        s = s+1;
        clearvars -except m ss perfAll sR2 dR2   toPlot sessionsMap s coefsS coefsO nShank
    end
end


toPlot.D{1,1} = sR2;
toPlot.D{1,2} = perfAll;

clearvars -except toPlot sessionsMap


%%% Load panel E data

s  = 1;

for m = 1:5
    % :length(intSess{m,2}(1,:))
    for ss = 1:length(sessionsMap{m,2}(1,:))
        
        intCond = [1 3];
        intCondO = [1 3];
        event=[1 6];
        eventO = [2 4];
        stdv = 0.14 ;
        alpha = 0.5;
        spdSmt = 500;
        
        
        Name = sessionsMap{m,2}{1,ss};
        mouseNumber = str2num(Name(3:4)); % number of the mouse being loaded
        rwddSound=mod(mouseNumber,2)+1; % define rewarded and and non rewarded sound based on animal number
        nRwddSound=2-mod(mouseNumber,2);
        sID = [rwddSound nRwddSound];
        sID(sID == nRwddSound) = 0;
        sID(sID == rwddSound) = 1;
        
        
        
        disp(['Processing session: ',num2str(ss),' of mice ',num2str(Name(1:4))])
        
        disp('Loading data')
        
        
        for c = 1:length(intCond(1,:))
            
            load(['/Volumes/Cerebro/Recording Data/',Name,'/Predict from Neurons/SpdPredResultsTimePaper Cond_',num2str(intCond(1,c)),' TSeg_',num2str(event),' spdSmt_',num2str(spdSmt),' Sigma_', num2str(stdv),' Alpha_',num2str(alpha),'.mat'],...
                'coefs','allPredMat');
            
            
            for d = 1:length(coefs(1,:))
                coefsS{s,d}{c} = coefs{1,d};
                featsS{s,d}{c} = allPredMat{1,d};
            end
            
            clear  coefs allPredMat
            
        end
        
        
        
        
        load(['/Volumes/Cerebro/Recording Data/',Name,'/Predict from Neurons/OutPredResultsTimePaper Cond_',num2str(intCondO),' TSeg_',num2str(eventO),' spdSmt_',num2str(spdSmt),' Sigma_', num2str(stdv),' Alpha_',num2str(alpha),'.mat'],...
            'coefs','allPredMat');
        
        
        for d = 1:length(coefs(1,:))
            coefsO{s,d}{1} = coefs{1,d};
            featsO{s,d}{1} = allPredMat{1,d};
        end
        
        
        clear  coefs  allPredMat
        
        
        
        load(['/Volumes/Cerebro/Recording Data/',Name,'/Neurons.mat'])
        nShank{s} = cell2mat(Neurons(:,3));
        
        
        s = s+1;
        clearvars -except m ss  sessionsMap s coefsS coefsO featsO featsS nShank toPlot
        
    end
end

toPlot.E{1,1} = coefsS; 
toPlot.E{1,2} = coefsO;
toPlot.E{1,3} = featsS;
toPlot.E{1,4} = featsO;
toPlot.E{1,5} = nShank;

clearvars -except toPlot sessionsMap


end






%%% Define size and position of figure and panels

figSize = [18.3 14.8];


panelA1Size = [3.5 2.5];
panelBSize = [2 2];
panelCSize = [3.5 2.5];
panelDSize = [2 2];
panelESize = [2 2];
panelGSize = [1.8 1.8];
panelFSize = [2.8 1.2];


panelA1Pos = [1 figSize(2)-panelA1Size(2)-0.5];
panelB1Pos = [1+panelA1Size(1)+1.5 figSize(2)-panelBSize(2)-0.7];
panelB2Pos = [1+panelA1Size(1)+panelBSize(1)+2.5 figSize(2)-panelBSize(2)-0.7];
panelCPos = [1+panelA1Size(1)+panelBSize(1)*2+4 figSize(2)-panelCSize(2)-0.6];
panelD1Pos = [1 figSize(2)-panelA1Size(2)-panelDSize(2)-2.2];
panelEPos = [1 figSize(2)-panelA1Size(2)-panelDSize(2)-panelESize(2)-3];
panelGPos = [1+panelDSize(1)+panelFSize(1)+3.6 figSize(2)-panelA1Size(2)-panelGSize(2)-2.2];
panelF1Pos = [1+panelDSize(1)+2 figSize(2)-panelA1Size(2)-panelFSize(2)-2.2];
panelF2Pos = [1+panelDSize(1)+2 figSize(2)-panelA1Size(2)-panelFSize(2)*2-2.25];
panelF3Pos = [1+panelDSize(1)+2 figSize(2)-panelA1Size(2)-panelFSize(2)*3-2.3];
panelF4Pos = [1+panelDSize(1)+2 figSize(2)-panelA1Size(2)-panelFSize(2)*4-2.35];


lAPos = [-0.8 panelA1Size(2)-0.25];
lBPos = [-1 panelBSize(2)-0.15];
lCPos = [-0.8 panelCSize(2)-0.25];
lDPos = [-0.8 panelDSize(2)-0.25];
lEPos = [-0.8 panelESize(2)-0.25];
lFPos = [-0.8 panelFSize(2)-0.25];
lGPos = [-1.2 panelGSize(2)-0.25];


% Define colors to use

spdColor = [0, 119, 187] / 255;
hitsColor = [0, 153, 136] / 255;
faColor = [204, 51, 17] / 255;
greyColor = [187, 187, 187] / 255;
lightGreyColor = [220, 220, 220] / 255;
soundColor = [238,51,119] / 255;




f3 = figure;
f3.Units = 'centimeters';
f3.Position = [20, 15, figSize(1), figSize(2)];




%% Panel A, Raw Data

cond1_spikes = toPlot.A{1,1}{1,1};
cond1_nSpikes = toPlot.A{1,1}{1,2};
cond1_tTimes = toPlot.A{1,1}{1,3};
cond1_trials = toPlot.A{1,1}{1,4};

cond3_spikes = toPlot.A{1,2}{1,1};
cond3_nSpikes = toPlot.A{1,2}{1,2};
cond3_tTimes =  toPlot.A{1,2}{1,3};
cond3_trials = toPlot.A{1,2}{1,4};


cond1_mFR = cond1_nSpikes ./ cond1_tTimes;
cond1_mFRSum = sum(cond1_mFR,2); 
c1_Info = [cond1_mFRSum cond1_tTimes(:,1)];



[~,sorted_c1Spks] = sort(cond1_mFRSum,'descend');


cond3_mFR = cond3_nSpikes ./ cond3_tTimes;
cond3_mFRSum = sum(cond3_mFR,2); 
[~,sorted_c3Spks] = sort(cond3_mFRSum,'descend');


t1_select = 4;
t3_select = 2;
spdSmt = 500;


axes
box on

spd = cond1_trials{1,sorted_c1Spks(t1_select)};
spn = spdSmt/length(spd);
spd = smooth(spd,spn,'lowess');
dist = cond1_trials{3,sorted_c1Spks(t1_select)};
dist = smooth(dist,spn,'lowess');
maxSpd = max(spd);
minSpd = min(spd);
spdOf = maxSpd - 20;

sTimes = cond1_trials{5,sorted_c1Spks(t1_select)};
aTimes = cond1_trials{6,sorted_c1Spks(t1_select)}; 
tms = cond1_trials{2,sorted_c1Spks(t1_select)};



for i = 1:2
    [~,idxA(i)] = min(abs((tms - aTimes(i))));
    [~,idxS(i)] = min(abs((tms - sTimes(i))));
end

sDists = round(dist(idxS) - dist(1));
aDists = round(dist(idxA) - dist(1));
sTimes = round(sTimes - tms(1),1);
aTimes = round(aTimes - tms(1),1);




nNeurs = length(cond1_spikes(1,:));

for s=1:length(cond1_spikes(sorted_c1Spks(t1_select),:))
    spikes = cond1_spikes{sorted_c1Spks(t1_select),s};
    for j=1:length(spikes(:,1))
        line([spikes(j,1) spikes(j,1)] - tms(1),[nNeurs-(s-1)-0.15 nNeurs-s+0.15],'Color',hitsColor)
    end
end

hold on


plot(tms(2:end)-tms(1),spd - spdOf,'Color',spdColor)
line([sTimes(1) sTimes(1)],[min(spd)-spdOf-5 nNeurs+2],'Color',[soundColor 0.3])
line([sTimes(2) sTimes(2)],[min(spd)-spdOf-5 nNeurs+2],'Color',[soundColor 0.3])
line([aTimes(1) aTimes(1)],[min(spd)-spdOf-5 nNeurs+2],'Color',greyColor)
line([aTimes(2) aTimes(2)],[min(spd)-spdOf-5 nNeurs+2],'Color',greyColor)

text(sTimes(1)-1.5,min(spd)-spdOf-30,'Time','FontSize',6,'Color',greyColor - 0.15)

text(sTimes(1)-0.2,min(spd)-spdOf-30,[num2str(sTimes(1)),' s'],'FontSize',6,'Color',greyColor - 0.15)
text(aTimes(1)-0.2,min(spd)-spdOf-30,[num2str(aTimes(1)),' s'],'FontSize',6,'Color',greyColor - 0.15)
text(aTimes(2)-0.2,min(spd)-spdOf-30,[num2str(aTimes(2)),' s'],'FontSize',6,'Color',greyColor - 0.15)


text(sTimes(1)-1.5,min(spd)-spdOf-40,'Posit','FontSize',6,'Color',greyColor - 0.15)

text(sTimes(1)-0.2,min(spd)-spdOf-40,[num2str(sDists(1)),' cm'],'FontSize',6,'Color',greyColor - 0.15)
text(aTimes(1)-0.2,min(spd)-spdOf-40,[num2str(aDists(1)),' cm'],'FontSize',6,'Color',greyColor - 0.15)
text(aTimes(2)-0.2,min(spd)-spdOf-40,[num2str(aDists(2)),' cm'],'FontSize',6,'Color',greyColor - 0.15)


%xlabel('Time(s) / Pos(cm)')

xlim([0 tms(end) - tms(1)])
ylim([min(spd)-spdOf-5 nNeurs+2])

axA(1) = gca;
axA(1).YTick = [(minSpd-spdOf) + ((maxSpd-spdOf) - (minSpd-spdOf))/2,round((nNeurs/2)+10)];
axA(1).YTickLabel = {'Speed','Neurons'};
axA(1).YTickLabelRotation = 90;
axA(1).XTick = [sTimes(1)+(sTimes(2)-sTimes(1))/2 aTimes(1)+(aTimes(2)-aTimes(1))/2];
axA(1).XTickLabel = {'Sound','Area'};
axA(1).FontSize = 6;
axA(1).Units = 'centimeters';
axA(1).Position = [panelA1Pos(1,1), panelA1Pos(1,2), panelA1Size(1), panelA1Size(2)];

lA = text(lAPos(1),lAPos(2),'a','Units','centimeters');
lA.FontSize = 9;
lA.FontWeight = 'bold';




%{
axes
box on

spd = cond3_trials{1,sorted_c3Spks(t3_select)};
spn = spdSmt/length(spd);
spd = smooth(spd,spn,'lowess');
dist = cond3_trials{3,sorted_c3Spks(t3_select)};
dist = smooth(dist,spn,'lowess');
maxSpd = max(spd);
minSpd = min(spd);
spdOf = maxSpd - 20;


tms = cond3_trials{2,sorted_c3Spks(t3_select)};
xTms = linspace(0,round(tms(end)-tms(1),2),4);

for i = 1:length(xTms)
    [~,idx(i)] = min(abs((tms-tms(1) - xTms(i))));
end

xDists = round(dist(idx) - dist(1));

for  i = 1:length(xTms)
   xLbls{i} = [num2str(xTms(i)),'/',num2str(xDists(i))];  
    
end



for s=1:length(cond3_spikes(sorted_c3Spks(t3_select),:))
    spikes = cond3_spikes{sorted_c3Spks(t3_select),s};
    for j=1:length(spikes(:,1))
        line([spikes(j,1) spikes(j,1)] - tms(1),[nNeurs-(s-1)-0.15 nNeurs-s+0.15],'Color',faColor)
    end
end

hold on

plot(tms(2:end)-tms(1),spd - spdOf,'Color',spdColor)
text(0.1,min(spd)-spdOf,num2str(round(minSpd)),'FontSize',6,'Color',spdColor)
text(0.1,max(spd)-spdOf,num2str(round(maxSpd)),'FontSize',6,'Color',spdColor)
xlabel('Time(s) / Pos(cm)')

xlim([0 tms(end) - tms(1)])
ylim([min(spd)-spdOf-5 nNeurs])

axA(2) = gca;
axA(2).YTick = [];
axA(2).XTick = xTms;
axA(2).XTickLabel = xLbls;
axA(2).Units = 'centimeters';
axA(2).FontSize = 6;
axA(2).Position = [panelA2Pos(1,1), panelA2Pos(1,2), panelA2Size(1), panelA2Size(2)];
%}

%% Panel B, Example Trials


finalSpeeds = toPlot.B{1,1}{1,1};
finalDummySpeeds = toPlot.B{1,1}{1,2};
finalDistsSpeeds = toPlot.B{1,1}{1,3};

finalOuts = toPlot.B{1,2}{1,1};
finalDummyOuts = toPlot.B{1,2}{1,2};
finalDistsOuts = toPlot.B{1,2}{1,3};

for c = 1:length(finalSpeeds(:,1))
    for d = 1:length(finalSpeeds(1,:))-1
        for t = 1:length(finalSpeeds{c,d}(:,1))
            spd = finalSpeeds{c,d}{t,1};
            pSpd = finalDummySpeeds{c,d}{t,1};
            
            out = finalOuts{c,d}{t,1};
            pOut = finalDummyOuts{c,d}{t,1};
            
            allErr{c,d}(t,1) = sqrt(sum((spd - pSpd).^2)/length(spd));
            allErr{c,d}(t,3) = sqrt(sum((out - pOut).^2)/length(spd));
        end
        meanErr{c}(d,:) = min(allErr{1,d});
    end
end

sS = 0:5:50;
sS(1) = 1;


% Example Spd

% Stop Trials


d = 7;
edges = 0:3:130;

delErr1 = allErr{1,d}(:,1);
[~,sTrials1] = sort(delErr1);
t1 = sTrials1(1);
spd1 = finalSpeeds{1,d}{t1,1};
pSpd1 = finalDummySpeeds{1,d}{t1,1};
dst1 = finalDistsSpeeds{1,d}{t1,1};

clear mSpd1 mPSpd1

for e = 2:length(edges)
    mSpd1(e-1) = mean(spd1(dst1 >= edges(e-1) & dst1 < edges(e)));
    mPSpd1(e-1) = mean(pSpd1(dst1 >= edges(e-1) & dst1 < edges(e)));
end

mSpd1(mSpd1<0) = 0;
mSpd1 = smooth(mSpd1,0.07,'lowess');

mPSpd1(mPSpd1<0) = 0;
mPSpd1 = smooth(mPSpd1,0.07,'lowess');

% No Stop Trials

d = 7;
edges = 0:3:130;

delErr2 = allErr{2,d}(:,1);
[~,sTrials2] = sort(delErr2);
t2 = sTrials2(1);
spd2 = finalSpeeds{2,d}{t2,1};
pSpd2 = finalDummySpeeds{2,d}{t2,1};
dst2 = finalDistsSpeeds{2,d}{t2,1};

clear mSpd mPSpd

for e = 2:length(edges)
    mSpd2(e-1) = mean(spd2(dst2 >= edges(e-1) & dst2 < edges(e)));
    mPSpd2(e-1) = mean(pSpd2(dst2 >= edges(e-1) & dst2 < edges(e)));
end

mSpd2(mSpd2<0) = 0;
mSpd2 = smooth(mSpd2,0.07,'lowess');

mPSpd2(mPSpd2<0) = 0;
mPSpd2 = smooth(mPSpd2,0.07,'lowess');





axes
box on



plot(edges(2:end),mSpd1,'color',hitsColor,'LineWidth',1)
hold on
plot(edges(2:end),mPSpd1,'color',[hitsColor 0.4],'LineWidth',1)
hold on
plot(edges(2:end),mSpd2,'color',faColor,'LineWidth',1)
hold on
plot(edges(2:end),mPSpd2,'color',[faColor 0.4],'LineWidth',1)

maxSpd = max([mSpd1;mPSpd1;mSpd2;mPSpd2]);

line([sS(d) sS(d)],[0 maxSpd+5],'Color',[soundColor 0.3])
line([sS(d) sS(d)]+9,[0 maxSpd+5],'Color',[soundColor 0.3])
line([103 103],[0 maxSpd+5],'Color',greyColor)
line([117 117],[0 maxSpd+5],'Color',greyColor)

ylim([0 maxSpd+5])
xlim([0 edges(end)])

ylabel('Speed')
xlabel('Position (cm)')
title('Speed')

axB(1) = gca;
axB(1).YTick = [0 round(maxSpd/2) round(maxSpd)];
axB(1).XTick = [sS(d)+4.5 110];
axB(1).XTickLabel = {'Sound','Area'};
axB(1).Units = 'centimeters';
axB(1).FontSize = 6;
axB(1).Position = [panelB1Pos(1,1), panelB1Pos(1,2), panelBSize(1), panelBSize(2)];

lB = text(lBPos(1),lBPos(2),'b','Units','centimeters');
lB.FontSize = 9;
lB.FontWeight = 'bold';



% Example Out

d = 10;
edges = 0:3:130;

delErr3 = allErr{2,d}(:,3);
[~,sTrials3] = sort(delErr3);
t3 = sTrials3(1);

delErr1 = allErr{1,d}(:,3);
[~,sTrials1] = sort(delErr1);
t1 = sTrials1(1);



pOut1 = finalDummyOuts{1,12}{1,d}{t1,1};
pOut3 = finalDummyOuts{2,12}{1,d}{t3,1};

dst1 = finalDistsSpeeds{1,d}{t1,1};
dst3 = finalDistsSpeeds{2,d}{t3,1};

clear  mPout1 mPout3

for e = 2:length(edges)
    mPOut1(e-1) = mean(pOut1(dst1 >= edges(e-1) & dst1 < edges(e)));
    mPOut3(e-1) = mean(pOut3(dst3 >= edges(e-1) & dst3 < edges(e)));
end


mPOut1 = smooth(mPOut1,0.07,'lowess');
mPOut3 = smooth(mPOut3,0.07,'lowess');



axes
box on



plot(edges(2:end),mPOut1,'color',[hitsColor 0.4],'LineWidth',1)
hold on
plot(edges(2:end),mPOut3,'color',[faColor 0.4],'LineWidth',1)

line([sS(d) sS(d)],[0 1],'Color',[soundColor 0.3])
line([sS(d) sS(d)]+9,[0 1],'Color',[soundColor 0.3])
line([103 103],[0 1],'Color',greyColor)
line([117 117],[0 1],'Color',greyColor)

line([sS(d) edges(end)],[1 1],'color',hitsColor)
line([sS(d) edges(end)],[0 0],'color',faColor)

ylabel('ST Prob')
xlabel('Position (cm)')
title('Trial ID')

xlim([0 edges(end)])
ylim([-0.1 1.1])


axB(2) = gca;
axB(2).YTick = [0 0.5 1];
axB(2).XTick = [sS(d)+4.5 110];
axB(2).XTickLabel = {'Sound','Area'};
axB(2).Units = 'centimeters';
axB(2).FontSize = 6;
axB(2).Position = [panelB2Pos(1,1), panelB2Pos(1,2), panelBSize(1), panelBSize(2)];

%% Panel C, mean Preds
axes
box on

edges = 0:3:130;

% Speeds Stop Trials

dS = 7;

for t = 1:length(finalSpeeds{1,dS}(:,1))
      spd1 = finalSpeeds{1,dS}{t,1};
      pSpd1 = finalDummySpeeds{1,dS}{t,1};
      dst1 = finalDistsSpeeds{1,dS}{t,1};
      
      for e = 2:length(edges)
          mSpd1(e-1) = mean(spd1(dst1 >= edges(e-1) & dst1 < edges(e)));
          mPSpd1(e-1) = mean(pSpd1(dst1 >= edges(e-1) & dst1 < edges(e)));
      end 

      
    mSpd1(~isnan(mSpd1)) =  smooth(mSpd1(~isnan(mSpd1)),0.08,'lowess');  
    mPSpd1(~isnan(mPSpd1)) =  smooth(mPSpd1(~isnan(mPSpd1)),0.08,'lowess');  
    
    meanTrialSpds1(t,:) = mSpd1;
    meanTrialPSpds1(t,:) = mPSpd1;  
      
   
%     meanSpds(t,:) = smooth(mSpd,0.07,'lowess');
%     meanPSpds(t,:) = smooth(mPSpd,0.07,'lowess'); 
    clear mSpd1 mPSpd1
end


mPSpd1 = mean(meanTrialPSpds1);
mSpd1 = mean(meanTrialSpds1);
xPatchS1 = edges(~isnan(mPSpd1));
mPSpd1 = mPSpd1(~isnan(mPSpd1));
mSpd1 = mSpd1(~isnan(mSpd1));
errPSpd1 = std(meanTrialPSpds1) / sqrt(length(meanTrialPSpds1(:,1)));
errPSpd1 = errPSpd1(~isnan(errPSpd1));
errSpd1 = std(meanTrialSpds1) / sqrt(length(meanTrialSpds1(:,1)));
errSpd1 = errSpd1(~isnan(errSpd1));

errPSupS1 = mPSpd1 + errPSpd1;
errPInfS1 = mPSpd1 - errPSpd1;

errSupS1 = mSpd1 + errSpd1;
errInfS1 = mSpd1 - errSpd1;

% Speeds No Stop Trials

dS = 7;

for t = 1:length(finalSpeeds{2,dS}(:,1))
      spd2 = finalSpeeds{2,dS}{t,1};
      pSpd2 = finalDummySpeeds{2,dS}{t,1};
      dst2 = finalDistsSpeeds{2,dS}{t,1};
      
      for e = 2:length(edges)
          mSpd2(e-1) = mean(spd2(dst2 >= edges(e-1) & dst2 < edges(e)));
          mPSpd2(e-1) = mean(pSpd2(dst2 >= edges(e-1) & dst2 < edges(e)));
      end 

      
    mSpd2(~isnan(mSpd2)) =  smooth(mSpd2(~isnan(mSpd2)),0.08,'lowess');  
    mPSpd2(~isnan(mPSpd2)) =  smooth(mPSpd2(~isnan(mPSpd2)),0.08,'lowess');  
    
    meanTrialSpds2(t,:) = mSpd2;
    meanTrialPSpds2(t,:) = mPSpd2;  
      
   
%     meanSpds(t,:) = smooth(mSpd,0.07,'lowess');
%     meanPSpds(t,:) = smooth(mPSpd,0.07,'lowess'); 
    clear mSpd2 mPSpd2
end


mPSpd2 = mean(meanTrialPSpds2);
mSpd2 = mean(meanTrialSpds2);
xPatchS2 = edges(~isnan(mPSpd2));
mPSpd2 = mPSpd2(~isnan(mPSpd2));
mSpd2 = mSpd2(~isnan(mSpd2));
errPSpd2 = std(meanTrialPSpds2) / sqrt(length(meanTrialPSpds2(:,1)));
errPSpd2 = errPSpd2(~isnan(errPSpd2));
errSpd2 = std(meanTrialSpds2) / sqrt(length(meanTrialSpds2(:,1)));
errSpd2 = errSpd2(~isnan(errSpd2));

errPSupS2 = mPSpd2 + errPSpd2;
errPInfS2 = mPSpd2 - errPSpd2;

errSupS2 = mSpd2 + errSpd2;
errInfS2 = mSpd2 - errSpd2;


% Out

dO = 7;

for t = 1:length(finalOuts{1,dO}(:,1))
      out1 = finalOuts{1,12}{1,dO}{t,1};
      pOut1 = finalDummyOuts{1,12}{1,dO}{t,1};
      dst = finalDistsOuts{1,dO}{t,1};
      
      for e = 2:length(edges)
          mOut1(e-1) = mean(out1(dst >= edges(e-1) & dst < edges(e)));
          mPOut1(e-1) = mean(pOut1(dst >= edges(e-1) & dst < edges(e)));
      end 

      
    mOut1(~isnan(mOut1)) =  smooth(mOut1(~isnan(mOut1)),0.08,'lowess');  
    mPOut1(~isnan(mPOut1)) =  smooth(mPOut1(~isnan(mPOut1)),0.08,'lowess');  
    
    meanTrialOut1(t,:) = mOut1;
    meanTrialPOut1(t,:) = mPOut1;  
      
  
    clear mOut1 mPOut1
end

for t = 1:length(finalOuts{2,dO}(:,1))
      out2 = finalOuts{2,12}{1,dO}{t,1};
      pOut2 = finalDummyOuts{2,12}{1,dO}{t,1};
      dst = finalDistsOuts{2,dO}{t,1};
      
      for e = 2:length(edges)
          mOut2(e-1) = mean(out2(dst >= edges(e-1) & dst < edges(e)));
          mPOut2(e-1) = mean(pOut2(dst >= edges(e-1) & dst < edges(e)));
      end 

      
    mOut2(~isnan(mOut2)) =  smooth(mOut2(~isnan(mOut2)),0.08,'lowess');  
    mPOut2(~isnan(mPOut2)) =  smooth(mPOut2(~isnan(mPOut2)),0.08,'lowess');  
    
    meanTrialOut2(t,:) = mOut2;
    meanTrialPOut2(t,:) = mPOut2;  
      
  
    clear mOut2 mPOut2
end

mPOut1 = mean(meanTrialPOut1);
xPatchO1 = edges(~isnan(mPOut1));
mPOut1 = mPOut1(~isnan(mPOut1));
errPOut1 = std(meanTrialPOut1) / sqrt(length(meanTrialPOut1(:,1)));
errPOut1 = errPOut1(~isnan(errPOut1));
errSupO1 = mPOut1 + errPOut1;
errInfO1 = mPOut1 - errPOut1;

mPOut2 = mean(meanTrialPOut2);
xPatchO2 = edges(~isnan(mPOut2));
mPOut2 = mPOut2(~isnan(mPOut2));
errPOut2 = std(meanTrialPOut2) / sqrt(length(meanTrialPOut2(:,1)));
errPOut2 = errPOut2(~isnan(errPOut2));
errSupO2 = mPOut2 + errPOut2;
errInfO2 = mPOut2 - errPOut2;


hold on
[axC,hLine11,hLine12] = plotyy(xPatchS1,mSpd1,sS(dO):130,ones(1,length(sS(dO):130)));
hLine11.Color = hitsColor;
hLine12.Color = hitsColor;
line([sS(dO) 130],[0 0],'Color',faColor,'parent',axC(2))
%line([0 130],[0 130],'Color',greyColor - 0.4,'parent',axC(1))

patch([xPatchS1 fliplr(xPatchS1)],[errPInfS1 fliplr(errPSupS1)],hitsColor,'EdgeColor','none','FaceAlpha',0.4,'parent',axC(1))
patch([xPatchS1 fliplr(xPatchS1)],[errInfS1 fliplr(errSupS1)],hitsColor,'EdgeColor','none','parent',axC(1))
patch([xPatchO1 fliplr(xPatchO1)],[errInfO1 fliplr(errSupO1)],[1 1 1],'EdgeColor',hitsColor,'EdgeAlpha',0.4,'FaceAlpha',0.4,'parent',axC(2))

text(sS(dO),100,'Sound','Color',soundColor,'FontSize',6)
text(103,100,'Area','Color',greyColor,'FontSize',6)


axC(1).XLim = [0 130];
axC(1).YLim = [-11 105];
axC(1).XTick = [0 65 130];
axC(1).YTick = [0 round(max(errPSupS2)/2) round(max(errPSupS2))];
axC(1).Units = 'centimeters';
axC(1).FontSize = 6;
axC(1).Position = [panelCPos(1,1), panelCPos(1,2), panelCSize(1), panelCSize(2)];


axC(2).YLim = [-0.1 1.1];
axC(2).XLim = [0 130];
axC(2).XTick = [];
axC(2).YTick = [0 0.5 1];
axC(2).Units = 'centimeters';
axC(2).FontSize = 6;
axC(2).Position = [panelCPos(1,1), panelCPos(1,2), panelCSize(1), panelCSize(2)];
hold on

patch([xPatchS2 fliplr(xPatchS2)],[errPInfS2 fliplr(errPSupS2)],faColor,'EdgeColor','none','FaceAlpha',0.4,'parent',axC(1))
patch([xPatchS2 fliplr(xPatchS2)],[errInfS2 fliplr(errSupS2)],faColor,'EdgeColor','none','parent',axC(1))
patch([xPatchO2 fliplr(xPatchO2)],[errInfO2 fliplr(errSupO2)],[1 1 1],'EdgeColor',faColor,'EdgeAlpha',0.4,'FaceAlpha',0.4,'parent',axC(2))


line([sS(dO) sS(dO)+9],[-5 -5],'Color',soundColor,'parent',axC(1))
line([103 117],[-5 -5],'Color',greyColor,'parent',axC(1))

xlabel(axC(1),'Position cm');
ylabel(axC(1),'Speed');
ylabel(axC(2),'ST Prob');

title('All Predictions')

lC = text(lCPos(1),lCPos(2),'c','Units','centimeters');
lC.FontSize = 9;
lC.FontWeight = 'bold';

%% Panel D, Preds performance

sR2 = toPlot.D{1,1};
oAUC = toPlot.D{1,2};

for s = 1:length(sR2(1,:))
    s1R2(s,:) = sR2{1,s}(1,:);
    s1R2(s,s1R2(s,:) < 0 ) = 0;
    s2R2(s,:) = sR2{1,s}(2,:);
    s2R2(s,s2R2(s,:) < 0 ) = 0;
end

intDels = [3 6 11 12];

mS1 = mean(s1R2(:,intDels));
mS2 = mean(s2R2(:,intDels));
eS1 = std(s1R2(:,intDels)) / sqrt(length(s1R2(:,1)));
eS2 = std(s2R2(:,intDels)) / sqrt(length(s2R2(:,1)));

mO = mean(oAUC(:,intDels));
eO = std(oAUC(:,intDels)) / sqrt(length(oAUC(:,1)));

errsS = [eS1' eS2'];
mnsS = [mS1' mS2'];

errsO = eO;
mnsO = mO;


% All Preds

axes 
box on

[axD, ~, ~]=plotyy([0,1],[0,1],[0,1],[0,1]); 
cla(axD(1)); 
cla(axD(2));


xS = [1 5 9 13];
xNS = [2 6 10 14];
xO = [3 7 11 15];

hold(axD(1),'on')
bar(axD(1),xS,mnsS(:,1),0.2,'FaceColor',hitsColor);
errorbar(axD(1),xS,mnsS(:,1),errsS(:,1),'.','Color','k')
bar(axD(1),xNS,mnsS(:,2),0.2,'FaceColor',faColor);
errorbar(axD(1),xNS,mnsS(:,2),errsS(:,2),'.','Color','k')
hold(axD(1),'off')
ylim(axD(1),[0 1.05])
xlim(axD(1),[-1 17])
ylabel(axD(1),'R2')
title(axD(1),'Model Perfs')

axD(1).YTick = [0 0.5 1];
axD(1).XTick = [2 6 10 14];
axD(1).XTickLabel = {'10','25','50','All'};
axD(1).YColor = 'k';
axD(1).Units = 'centimeters';
axD(1).FontSize = 6;
axD(1).Position = [panelD1Pos(1,1), panelD1Pos(1,2), panelDSize(1), panelDSize(2)];

hold(axD(2),'on')
bar(axD(2),xO,mnsO,0.2,'FaceColor',spdColor,'EdgeColor',greyColor)
errorbar(axD(2),xO,mnsO,errsO,'.','Color','k')
hold(axD(2),'off')
ylim(axD(2),[0 1.05])
xlim(axD(2),[-1 17])
ylabel(axD(2),'AUC')

axD(2).YTick = [0 0.5 1];
axD(2).YColor = greyColor;
axD(2).Units = 'centimeters';
axD(2).FontSize = 6;


lD = text(lDPos(1),lDPos(2),'d','Units','centimeters');
lD.FontSize = 9;
lD.FontWeight = 'bold';





%% Panel E, Plot Distribution of encodings

coefsS = toPlot.E{1,1};
coefsO = toPlot.E{1,2};
featsS = toPlot.E{1,3};
featsO = toPlot.E{1,4};
nShank = toPlot.E{1,5};


% Calculate correlations between predictions and var explained by each neuron in each prediction
d=12;

for s = 1:length(coefsS(:,1))
    
    % get coefs
    cS1{s} = coefsS{s,d}{1,1}(2:end)';
    cS2{s} = coefsS{s,d}{1,2}(2:end)';
    cO{s} = coefsO{s,d}{1,1}(2:end)';
    
    % get frs
    frS1{s} = featsS{s,d}{1,1};
    frS2{s} = featsS{s,d}{1,2}; 
    frO{s} = featsO{s,d}{1,1};
    
    % calculate correlations between coefs vectors
    c1 = corrcoef(cS1{s},cS2{s});
    c2 = corrcoef(cS1{s},cO{s});
    c3 = corrcoef(cS2{s},cO{s});
    
    corrs(s,:) = [c1(1,2), c2(1,2), c3(1,2)];
   
    % calculate feature importance
    
    
    % Stop trial speeds
    allFrS1 = frS1{s};
    allCoefS1 = cS1{s};
    allImpS1 = var(allFrS1 .* repmat(allCoefS1,length(allFrS1(:,1)),1));
    
    featImp{s}(:,1) =  allImpS1 ./ repmat(sum(allImpS1),1,length(allImpS1));
    
    % No Stop trial speeds
    allFrS2 = frS2{s};
    allCoefS2 = cS2{s};
    allImpS2 = var(allFrS2 .* repmat(allCoefS2,length(allFrS2(:,1)),1));
    
    featImp{s}(:,2) =  allImpS2 ./ repmat(sum(allImpS2),1,length(allImpS2));
    
    % Trial outcomes
    allFrO = frO{s};
    allCoefO = cO{s};
    allImpO = var(allFrO .* repmat(allCoefO,length(allFrO(:,1)),1));
    
    featImp{s}(:,3) =  allImpO ./ repmat(sum(allImpO),1,length(allImpO));
   
end


% calculate cumsums of feat imp and # of neurs untill treshold
tresh = 0.8;
for s = 1:16
    
    featImpS1{s} = cumsum(sort(featImp{s}(:,1),'descend'));
    featImpS2{s} = cumsum(sort(featImp{s}(:,2),'descend'));
    featImpO{s} = cumsum(sort(featImp{s}(:,3),'descend'));
    
    [~,treshNum(s,1)] = min(abs(featImpS1{s} - tresh));
    [~,treshNum(s,2)] = min(abs(featImpS2{s} - tresh));
    [~,treshNum(s,3)] = min(abs(featImpO{s} - tresh));
    
end


% Plot neurons variance explained in all predictions 

axes
box on

[axE h1 h2]=plotyy([0,1],[0,1],[0,1],[0,1]); 
cla(axE(1)); 
cla(axE(2));


hold(axE(1),'on')
plot(axE(1),1:length(featImpS1{s}),featImpS1{s},'Color',hitsColor)
plot(axE(1),1:length(featImpS1{s}),featImpS2{s},'Color',faColor)
plot(axE(1),1:length(featImpS1{s}),featImpO{s},'Color',spdColor)
line([0 length(featImpS1{s})],[0.8 0.8],'Color','r','LineStyle','--','parent',axE(1))
hold(axE(1),'off')
ylim(axE(1),[0 1.05])
xlim(axE(1),[0 100])
ylabel(axE(1),'Expl Var')
axE(1).XTick = [];
axE(1).YColor = 'k';
axE(1).Units = 'centimeters';
axE(1).FontSize = 6;
%axE(1).Position = [panelE2Pos(1,1), panelE2Pos(1,2), panelESize(1), panelESize(2)];


hold(axE(2),'on')
bar(axE(2),1,mean(treshNum(:,1)),'FaceColor',hitsColor,'EdgeColor',greyColor);
bar(axE(2),2.5,mean(treshNum(:,2)),'FaceColor',faColor,'EdgeColor',greyColor);
bar(axE(2),4,mean(treshNum(:,3)),'FaceColor',spdColor,'EdgeColor',greyColor);
errorbar([1 2.5 4],mean(treshNum),std(treshNum)/sqrt(length(treshNum(:,1))),'.','Color','k','parent',axE(2))
text(0.1,-4.8,'Neurons  /','FontSize',6,'parent',axE(2))
text(3,-4.8,'Models','FontSize',6,'Color',greyColor,'parent',axE(2))
hold(axE(2),'off')
ylim(axE(2),[0 40])
xlim(axE(2),[-0.5 5])
ylabel(axE(2),'# Neurons to 0.8')
axE(2).YTick = [0 12 24];
axE(2).YColor = greyColor;
axE(2).XColor = 'k';
axE(2).XTick = [];
axE(2).Units = 'centimeters';
axE(2).FontSize = 6;
axE(2).Position = [panelEPos(1,1), panelEPos(1,2), panelESize(1), panelESize(2)];

lE = text(lEPos(1),lEPos(2),'e','Units','centimeters');
lE.FontSize = 9;
lE.FontWeight = 'bold';


%% Panel G, Correlation between encoders

%Plot correlation between encoders

mSizeScat = 10;
mSizeErr = 2;


axes
box on

scatter(ones(1,length(corrs)),corrs(:,1),mSizeScat,hitsColor,'filled','MarkerFaceAlpha',0.2)
hold on
errorbar(1,mean(corrs(:,1)),std(corrs(:,1))/sqrt(length(corrs(:,1))),'o','Color','k','MarkerFaceColor',hitsColor,'MarkerEdgeColor',hitsColor,'MarkerSize',mSizeErr)
hold on
scatter(ones(1,length(corrs))*2,corrs(:,2),mSizeScat,faColor,'filled','MarkerFaceAlpha',0.2)
hold on
errorbar(2,mean(corrs(:,2)),std(corrs(:,2))/sqrt(length(corrs(:,2))),'o','Color','k','MarkerFaceColor',faColor,'MarkerEdgeColor',faColor,'MarkerSize',mSizeErr)
hold on
scatter(ones(1,length(corrs))*3,corrs(:,3),mSizeScat,spdColor,'filled','MarkerFaceAlpha',0.2)
hold on
errorbar(3,mean(corrs(:,3)),std(corrs(:,3))/sqrt(length(corrs(:,3))),'o','Color','k','MarkerFaceColor',spdColor,'MarkerEdgeColor',spdColor,'MarkerSize',mSizeErr)
xlim([0.5 3.5])
ylim([-1 1])
ylabel('Corr Coef')

line([0.5 3.5],[0 0],'Color',greyColor,'LineStyle','--')

axF = gca;
axF.Box = 'on';
axF.XTick = 2;
axF.XTickLabel = 'Models';
axF.Units = 'centimeters';
axF.FontSize = 6;
axF.Position = [panelGPos(1,1), panelGPos(1,2), panelGSize(1), panelGSize(2)];

lG = text(lGPos(1),lGPos(2),'g','Units','centimeters');
lG.FontSize = 9;
lG.FontWeight = 'bold';


%% Panel F, Fraction of Variance per shanks

%Group var exp per shank

shanks = 1:6;

for s = 1:16
    for sk = shanks
        
         [s1Fi, s1FiIdx] = sort(featImp{s}(:,1),'descend');
         s1Shanks = nShank{s}(s1FiIdx);
         s1Fi = s1Fi(1:treshNum(s,1));
         s1Shanks = s1Shanks(1:treshNum(s,1));
         
         [s2Fi, s2FiIdx] = sort(featImp{s}(:,2),'descend');
         s2Shanks = nShank{s}(s2FiIdx);
         s2Fi = s2Fi(1:treshNum(s,2));
         s2Shanks = s2Shanks(1:treshNum(s,2));
         
         [oFi, oFiIdx] = sort(featImp{s}(:,3),'descend');
         oShanks = nShank{s}(oFiIdx);
         oFi = oFi(1:treshNum(s,3));
         oShanks = oShanks(1:treshNum(s,3));
        
        
        nSkImps{1}{s,sk} = s1Fi(s1Shanks == sk,1);
        nSkImps{2}{s,sk} = s2Fi(s2Shanks == sk,1);
        nSkImps{3}{s,sk} = oFi(oShanks == sk,1);
        
        skImp{1}(s,sk) =  mean(nSkImps{1}{s,sk});
        skImp{2}(s,sk) =  mean(nSkImps{2}{s,sk});
        skImp{3}(s,sk) =  mean(nSkImps{3}{s,sk});
        
    end
    
end


% Example session

s = 13;
filt = nShank{s} ~= 7;
mSizeScat = 5;
mSizeErr = 1.5;


axes
box on

scatter(nShank{s}(filt)-0.25,featImp{s}(filt,1),mSizeScat,hitsColor,'filled','MarkerFaceAlpha',0.2) 
hold on
scatter(nShank{s}(filt),featImp{s}(filt,2),mSizeScat,faColor,'filled','MarkerFaceAlpha',0.2)
hold on
scatter(nShank{s}(filt)+0.25,featImp{s}(filt,3),mSizeScat,spdColor,'filled','MarkerFaceAlpha',0.2)

for i = 1:6
errorbar(i-0.25,mean(featImp{s}(nShank{s}==i,1)),std(featImp{s}(nShank{s}==i,1))/sqrt(sum(nShank{s}==i)),'o','MarkerFaceColor',hitsColor,'MarkerEdgeColor',hitsColor,'Color','k','MarkerSize',mSizeErr)
hold on
errorbar(i,mean(featImp{s}(nShank{s}==i,2)),std(featImp{s}(nShank{s}==i,2))/sqrt(sum(nShank{s}==i)),'o','MarkerFaceColor',faColor,'MarkerEdgeColor',faColor,'Color','k','MarkerSize',mSizeErr)
hold on
errorbar(i+0.25,mean(featImp{s}(nShank{s}==i,3)),std(featImp{s}(nShank{s}==i,3))/sqrt(sum(nShank{s}==i)),'o','MarkerFaceColor',spdColor,'MarkerEdgeColor',spdColor,'Color','k','MarkerSize',mSizeErr)
end
xlim([0 7])
ylim([-0.02 0.15])



axG = gca;
axG(1).Box = 'on';
axG(1).XTick = [];
axG(1).Units = 'centimeters';
axG(1).FontSize = 6;
axG(1).Position = [panelF1Pos(1,1), panelF1Pos(1,2), panelFSize(1), panelFSize(2)];

lF = text(lFPos(1),lFPos(2),'f','Units','centimeters');
lF.FontSize = 9;
lF.FontWeight = 'bold';

% All Sessions

allSnksS.S1 = cell2mat(nSkImps{1,1}(:,1));
allSnksS.S2 = cell2mat(nSkImps{1,1}(:,2));
allSnksS.S3 = cell2mat(nSkImps{1,1}(:,3));
allSnksS.S4 = cell2mat(nSkImps{1,1}(:,4));
allSnksS.S5 = cell2mat(nSkImps{1,1}(:,5));
allSnksS.S6 = cell2mat(nSkImps{1,1}(:,6));

offSet = 0;
offSet2 = 0;

allSnksNS.S1 = cell2mat(nSkImps{1,2}(:,1))+offSet;
allSnksNS.S2 = cell2mat(nSkImps{1,2}(:,2))+offSet;
allSnksNS.S3 = cell2mat(nSkImps{1,2}(:,3))+offSet;
allSnksNS.S4 = cell2mat(nSkImps{1,2}(:,4))+offSet;
allSnksNS.S5 = cell2mat(nSkImps{1,2}(:,5))+offSet;
allSnksNS.S6 = cell2mat(nSkImps{1,2}(:,6))+offSet;

allSnksO.S1 = cell2mat(nSkImps{1,3}(:,1))+ offSet2;
allSnksO.S2 = cell2mat(nSkImps{1,3}(:,2))+ offSet2;
allSnksO.S3 = cell2mat(nSkImps{1,3}(:,3))+ offSet2;
allSnksO.S4 = cell2mat(nSkImps{1,3}(:,4))+ offSet2;
allSnksO.S5 = cell2mat(nSkImps{1,3}(:,5))+ offSet2;
allSnksO.S6 = cell2mat(nSkImps{1,3}(:,6))+ offSet2;


% hold on
% violinplot(allSnksNS,[],'ViolinColor',faColor);
% hold on
% violinplot(allSnksO,[],'ViolinColor',spdColor);




axes
box on

violinplot(allSnksS,[],'ViolinColor',hitsColor,'ShowMean',true,'width',0.3,'ShowData',false);
xlim([0 7])
ylim([-0.05 0.5])
yl = ylabel('Expl Var');
yl.Position = [-1.3 -0.1];

axG(2) = gca;
axG(2).Box = 'on';
axG(2).XTick = 1:6;
axG(2).Units = 'centimeters';
axG(2).FontSize = 6;
axG(2).Position = [panelF2Pos(1,1), panelF2Pos(1,2), panelFSize(1), panelFSize(2)];



axes
box on

violinplot(allSnksNS,[],'ViolinColor',faColor,'ShowMean',true,'width',0.3,'ShowData',false);
xlim([0 7])
ylim([-0.05 0.3])
axG(3) = gca;
axG(3).Box = 'on';
axG(3).XTick = 1:6;
axG(3).Units = 'centimeters';
axG(3).FontSize = 6;
axG(3).Position = [panelF3Pos(1,1), panelF3Pos(1,2), panelFSize(1), panelFSize(2)];


axes
box on

violinplot(allSnksO,[],'ViolinColor',spdColor,'ShowMean',true,'width',0.3,'ShowData',false);
xlim([0 7])
ylim([-0.05 0.3])
xlabel('Shanks')

axG(4) = gca;
axG(4).Box = 'on';
axG(4).XTick = 1:6;
axG(4).XTickLabel = {'1','2','3','4','5','6'};
axG(4).Units = 'centimeters';
axG(4).FontSize = 6;
axG(4).Position = [panelF4Pos(1,1), panelF4Pos(1,2), panelFSize(1), panelFSize(2)];
 

%{



s = 13;

axes
box on
plot(featImpS1{s},'Color',hitsColor)
hold on
plot(featImpS2{s},'Color',faColor)
hold on
plot(featImpO{s},'Color',spdColor)
hold on
text(48,0.7,'0.8 var expl','FontSize',5)





line([0 length(featImpS1{s})],[0.8 0.8],'Color','r','LineStyle','--')
ylim([0 1.05])
xlim([0 length(featImpS1{s})])
xlabel('Sorted Neurons')
ylabel('Fract Var Explained')
title('Example Session')

axE(1) = gca;
axE(1).XTick = [0 ceil(length(featImpS1{s})/2) length(featImpS1{s})];
axE(1).Units = 'centimeters';
axE(1).FontSize = 6;
axE(1).Position = [panelE1Pos(1,1), panelE1Pos(1,2), panelESize(1), panelESize(2)];

lE = text(lEPos(1),lEPos(2),'e','Units','centimeters');
lE.FontSize = 9;
lE.FontWeight = 'bold';

% All sessions 

axes
box on

bar(1,mean(treshNum(:,1)),'FaceColor',hitsColor);
hold on
bar(2.5,mean(treshNum(:,2)),'FaceColor',faColor);
hold on
bar(4,mean(treshNum(:,3)),'FaceColor',spdColor);

hold on 
errorbar([1 2.5 4],mean(treshNum),std(treshNum)/sqrt(length(treshNum(:,1))),'.','Color','k')

xlim([0 5])
ylabel('# neurons to treshold')
title('All Sessions')

axE(2) = gca;
axE(2).XTick = [1 2.5 4];
axE(2).XTickLabel = {'Sp.S','Sp.NS','T.Id'};
axE(2).Units = 'centimeters';
axE(2).FontSize = 6;
axE(2).Position = [panelE2Pos(1,1), panelE2Pos(1,2), panelESize(1), panelESize(2)];

%}

disp('Saving fig')
tic
filename = '/Volumes/GoogleDrive/O meu disco/paper/RawFig';
print(f3,filename,'-depsc','-r300')
disp('Finished saving fig')
toc



clear all
close all

%% Load and process needed Data

load('/Volumes/Cerebro/Recording Data/SessionsMap.mat')

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
            
            %             LapsS{c,1} = Laps{1,12};
            
            for d = 1:length(Laps(1,:))
                coefsS{s,d}{c} = coefs{1,d};
            end
            %             dummyLapsS{c,1} = dummyLaps{1,12};
            %             %dummyLapsS{c,1} = (allPredMat{1,12} * coefsS{c,1}(2:end)) + coefsS{c,1}(1);
            %             outTrialsS{c,1} = allTrials{1,12};
            %
            
            clear Laps dummyLaps allTrials coefs allPredMat
            
        end
        
        
        
        
        load(['/Volumes/Cerebro/Recording Data/',Name,'/Predict from Neurons/OutPredResultsTimePaper Cond_',num2str(intCondO),' TSeg_',num2str(eventO),' spdSmt_',num2str(spdSmt),' Sigma_', num2str(stdv),' Alpha_',num2str(alpha),'.mat'],...
            'Laps','dummyLaps','allTrials','coefs','allDists');
        
        
        
        %             if sID(c) == 2
        %                 tType = 0;
        %             end
        
        
        
        %             LapsO{c,1} = Laps{1,12}(allTrials{1,12}(:,1) == tType,1);
        
        %     if c == 2
        %         LapsO{c,1} = LapsO{c,1} - 2;
        %     end
        
        %dummyLaps{1,12} = ones(length(dummyLaps{1,12}(:,1)),1) - dummyLaps{1,12}(:,1);
        %             dummyLapsO{c,1} = dummyLaps{1,12}(allTrials{1,12}(:,1) == tType,1);
        %             outTrialsO{c,1} = allTrials{1,12}(allTrials{1,12}(:,1) == tType,:);
        %             outDistsO{c,1} =  allDists{1,12}(allTrials{1,12}(:,1) == tType,1);
        
        for d = 1:length(Laps(1,:))
            coefsO{s,d}{c} = coefs{1,d};
        end
        
        
        clear Laps dummyLaps allTrials coefs allPredMat
        
        
        %load(['/Volumes/Cerebro/Recording Data/',Name,'/OrganizedSessionInfo.mat'],'validSession','organizedTrials','timestampsPerTrial')
        
        %load(['/Volumes/Cerebro/Recording Data/',Name,'/Predict from Neurons/Neur-OutPredTime Cond_',num2str(intCondO),' TSeg_',num2str(eventO),' KWdth_',num2str(kWdthS),' BinSz_',num2str(binSize),' SB_',num2str(SB),'.mat'],'Distances')
        
        % Separate speed values by trial and delay type
        
        s = s+1;
        clearvars -except m ss  sessionsMap s coefsS coefsO
    end
end
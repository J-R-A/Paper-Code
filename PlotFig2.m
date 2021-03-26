clearvars -except toPlot
close all


%%% Script to plot the different panels of the paper second figure


%% Load and process needed Data

load('/Volumes/Cerebro/Recording Data/SessionsMap.mat')

%%% Panel B

if ~exist('toPlot','var')
    
    disp('Loading Panel B Data')
    
    
    
    conditions=[1 3]; % Trial conditions to plot
    
    mice = {2; 5}; % Mice from which to pick neurons
    days = {1 ; 2}; % Days of the recording sessions, from the picked mice, from which to pick neurons
    units = {{12} ; {14}}; % Number of the units to plot
    soundStarts = {{[3 6 11]} ; {[3 6 11]}}; % Sound Starts of the trials to plot
    intSess = [mice days units soundStarts]; % Structure with all the info
    
    
    a=1;
    for m = 1:length(intSess(:,1)) % for all the picked mice
        for s = 1:length(intSess{m,2}(1,:)) % for all the selected sessions
            
            Name = sessionsMap{intSess{m,1},2}{1,intSess{m,2}(1,s)} % Name of the session
            mouseNumber = Name(3:4); % number of the mouse being loaded
            sessDate = sessionsMap{intSess{m,1},3}{1,intSess{m,2}(1,s)}; % Date of the session being loaded
            
            % load session neural data
            if sessionsMap{intSess{m,1},4}{1,intSess{m,2}(1,s)} == 1
                load(['/Volumes/Cerebro/Recording Data/',Name,'/Engaged/OrganizedNeuralData.mat'])
            elseif sessionsMap{intSess{m,1},4}{1,intSess{m,2}(1,s)} == 0
                load(['/Volumes/Cerebro/Recording Data/',Name,'/OrganizedNeuralData.mat'])
            end
            
            load(['/Volumes/Cerebro/Recording Data/',Name,'/Neurons.mat'])
            
            % load session behavioral data
            load(['/Volumes/Cerebro/Recording Data/',Name,'/OrganizedSessionInfo.mat'])
            load(['/Volumes/GoogleDrive/O meu disco/Data/DelayedResponseSimple/bMouse',num2str(mouseNumber),'/SessionsInfo/',sessDate,'.mat'])
            
            
            psthBinSize = 4; % Psth Bin size in 0.5 cm units
            numBins = length(meanFRCondition{1,1}{1,1}{1,1}(1,:));
            psthBins = round(linspace(1,numBins,numBins/psthBinSize));
            
            
            for n = 1 : length(intSess{m,3}{s,1}(1,:))
                neur = intSess{m,3}{s,1}(1,n);
                for d = 1:length(intSess{m,4}{s,1}(1,:))
                    del = intSess{m,4}{s,1}(1,d);
                    for c = 1 : length(conditions(1,:))
                        if ~isempty(meanFRCondition{1,del}{1,conditions(c)})
                            for b = 1:length(psthBins(1,:))-1
                                %toPlot{a,1}{1,d}{1,c}(1,b) = mean(meanFRCondition{1,del}{1,conditions(c)}{neur,1}(psthBins(b):psthBins(b+1)));
                                toPlot.B{a,1}{1,d}{1,c}(:,b) = mean(histogramMatricesBinned{1,del}{1,conditions(c)}{neur,1}(:,psthBins(b):psthBins(b+1)),2);
                            end
                        end
                    end
                end
                toPlot.B{a,2} = intSess{m,4}{s,1};
                toPlot.B{a,3} = round(SessionInfo.SessionParameters.AreasMatrix(1,:));
                toPlot.B{a,4} = psthBins;
                a=a+1;
            end
            clearvars -except m s sessionsMap intSess  conditions toPlot a
        end
    end
    
    
    
    
    %%% Panel C
    
    clearvars -except toPlot sessionsMap
    disp('Loading Panel C Data')
    
    
    %%% Get pValues data
    
    for m = 1:length(sessionsMap(:,1))
        for s = 1:length(sessionsMap{m,2}(1,:))
            Name = sessionsMap{m,2}{1,s};
            
            disp(['Processing session: ',num2str(s),' of mice ',num2str(Name(1:4))])
            
            if sessionsMap{m,4}{1,s} == 0
                load(['/Volumes/Cerebro/Recording Data/',Name,'/pValuesSelectivityC1  3_3 Segments.mat'])
            elseif sessionsMap{m,4}{1,s} == 1
                load(['/Volumes/Cerebro/Recording Data/',Name,'/Engaged/pValuesSelectivityC1  3_3 Segments.mat'])
            end
            
            toPlot.C{1,1}{s,m} = [pValuesFRConDiff{1,1}(:,2) pValuesFRConDiff{1,2}(:,2) pValuesFRConDiff{1,3}(:,2)];
            
            clearvars -except m s sig sessionsMap toPlot
            
        end
    end
    
    
    
    
    %%% Get FR data fro example psth
    
    pM = 1; % Mouse Number
    pS = 1; % Session Number
    pN = 9; % Neurons Number
    pD = 6; % Delay Number
    
    
    Name = sessionsMap{pM,2}{1,pS};
    mouseNumber = Name(3:4);
    sessDate = sessionsMap{pM,3}{1,pS};
    if sessionsMap{pM,4}{1,pS} == 1
        load(['/Volumes/Cerebro/Recording Data/',Name,'/Engaged/OrganizedNeuralData.mat'])
    elseif sessionsMap{pM,4}{1,pS} == 0
        load(['/Volumes/Cerebro/Recording Data/',Name,'/OrganizedNeuralData.mat'])
    end
    load(['/Volumes/Cerebro/Recording Data/',Name,'/OrganizedSessionInfo.mat'],'validSession')
    load(['/Volumes/GoogleDrive/O meu disco/Data/DelayedResponseSimple/bMouse',num2str(mouseNumber),'/SessionsInfo/',sessDate,'.mat'])
    
    toPlot.C{1,3} = meanFRCondition;
    
    area = round(SessionInfo.SessionParameters.AreasMatrix(1,:)); % area start location in cms
    area = area(1)*2; % area start location in .5cm bins
    conditions = [1 3]; % conditions
    soundStarts = unique(cell2mat(validSession(9,:))); % unique sound start locations
    soundStarts = soundStarts(soundStarts > 0); % valid unique soiund start locations
    soundEnds = (soundStarts+8)*2; % sound end location location in .5cm bins
    boundaries = 10;% Distance in 0.5cm bins from first and lasted divisions edges to sound end and area start
    nrSegments = 3; % Number of delay divisions to test for difference between conditions
    
    
    
    
    % Mean FRS of selected neuron, per condition and delay, in each of the 3 delay segments
    
    for s = 1:nrSegments
        for n = 1:length(pN)
            for c = 1:length(conditions(1,:))
                for d = 1 : length(soundStarts(1,:))
                    segments = round(linspace(soundEnds(d)+boundaries,area-boundaries,nrSegments+1));
                    meanFRS{n,1}{1,conditions(c)}{d,s} = sum(histogramMatrices{1,d}{1,conditions(c)}{pN(n),1}(:,11+segments(1,s):11+segments(1,s+1)),2)./(allNeuronsPerCondition{1,d}{3,conditions(c)}(:,11+segments(1,s+1))-allNeuronsPerCondition{1,d}{3,conditions(c)}(:,11+segments(1,s)));
                end
            end
        end
    end
    
    
    %  Mean FRS of selected neuron, per condition, in each of the 3 delay segments, with all delays concatenated
    
    for n = 1:length(meanFRS(:,1))
        for c = 1:length(conditions(1,:))
            for s = 1:nrSegments
                toPlot.C{1,2}{n,1}{1,conditions(c)}{1,s} = vertcat(meanFRS{n,1}{1,conditions(c)}{:,s});
            end
        end
    end
    
    
    
    
    
    
    
    %%% Panel D
    clearvars -except toPlot sessionsMap
    
    disp('Loading Panel D Data')
    
    
    conditions=[1 3]; % Trial conditions to plot
    
    mice = {1; 2}; % Mice from which to pick neurons
    days = {1;1}; % Days of the recording sessions, from the picked mice, from which to pick neurons
    units = {25;18}; % Number of the units to plot
    taskEv = {'SON';'AON'}; % Task Events for which to plot rasters and histograms
    intSess = [mice days units taskEv]; % Structure with all the info
    
    
    
    
    
    a=1;
    for m = 1:length(intSess(:,1)) % for all the picked mice
        
        
        Name = sessionsMap{intSess{m,1},2}{1,intSess{m,2}} % Name of the session
        
        % load session neural data
        if sessionsMap{intSess{m,1},4}{1,intSess{m,2}} == 1
            load(['/Volumes/Cerebro/Recording Data/',Name,'/Engaged/OrganizedNeuralDataTime',intSess{m,4},'.mat'])
        elseif sessionsMap{intSess{m,1},4}{1,intSess{m,2}} == 0
            load(['/Volumes/Cerebro/Recording Data/',Name,'/OrganizedNeuralDataTime',intSess{m,4},'.mat'])
        end
        
        neur = intSess{m,3};
        
        if strcmp(intSess{m,4},'RON')
            toPlot.D{m,1} = cell2mat({cat(1, allNeuronsPerCondition{12}{2,conditions(1)}{neur,:})});
            toPlot.D{m,4}{1,1} = allNeuronsPerCondition{12}{2,conditions(1)}(neur,:);
            toPlot.D{m,5} = length(allNeuronsPerCondition{12}{2,conditions(1)}(1,:));
        else
            for c = 1 : length(conditions(1,:))
                toPlot.D{m,c} = histogramMatricesBinned{1,12}{1,conditions(c)}{neur,1};
                toPlot.D{m,4}{c,1} = allNeuronsPerCondition{12}{2,conditions(c)}(neur,:);
                toPlot.D{m,5}(1,c) = length(toPlot.D{m,c}(1,:));
            end
        end
        toPlot.D{m,3} = intSess{m,4};
        toPlot.D{m,6} = allNeuronsPerCondition{12}{3,conditions(1)}(1,1:end-1);
        
        
        clearvars -except m  sessionsMap intSess  conditions toPlot pValues
    end
    
    
    
    
    %%% Panel E
    clearvars -except toPlot sessionsMap
    disp('Loading Panel E Data')
    
    
    
    % Get data about neurons' response
    
    a=1;
    for m = 1:length(sessionsMap(:,1))
        for s = 1:length(sessionsMap{m,2}(1,:))
            Name = sessionsMap{m,2}{1,s};
            
            if sessionsMap{m,4}{1,s} == 0
                load(['/Volumes/Cerebro/Recording Data/',Name,'/AllNeurPVals.mat'])
            elseif sessionsMap{m,4}{1,s} == 1
                load(['/Volumes/Cerebro/Recording Data/',Name,'/Engaged/AllNeurPVals.mat'])
            end
            
            allUnitsTemp{a,1} = neurPatterns{1,1};
            allUnitsTemp{a,2} = neurPatterns{1,3};
            a=a+1;
            clearvars -except m s sessionsMap a allUnitsTemp toPlot
            
        end
    end
    
    toPlot.E{1,1}{1,1} = vertcat(allUnitsTemp{:,1});
    toPlot.E{1,1}{1,2} = vertcat(allUnitsTemp{:,2});
    
    
    
    % Get data for example neuron
    
    mE = 2; % Mouse to get example neuron
    sE = 1; % Session to get example neuron
    nE = 10; % Number of example neuron
    
    
    Name =  sessionsMap{mE,2}{1,sE}; % Name of session
    
    if  sessionsMap{mE,4}{1,sE} == 1
        load(['/Volumes/Cerebro/Recording Data/',Name,'/Engaged/OrganizedNeuralDataTimeRON.mat'])
    elseif  sessionsMap{mE,4}{1,sE} == 0
        load(['/Volumes/Cerebro/Recording Data/',Name,'/OrganizedNeuralDataTimeRON.mat'])
    end
    
    
    cond = 1;
    
    toPlot.E{1,2}{1,cond} = histogramMatricesBinned{1,12}{1,cond}{nE,1};
    toPlot.E{1,2}{1,3} = allNeuronsPerCondition{12}{2,cond}(nE,:);
    toPlot.E{1,2}{1,4} = length(allNeuronsPerCondition{12}{2,cond}(1,:));
    
    
    toPlot.E{1,2}{1,2} = 'RON';
    toPlot.E{1,2}{1,5} = allNeuronsPerCondition{12}{3,cond}(1,1:end-1);
    
    clearvars -except toPlot sessionsMap 
end




%%% Plot Pannels





%%% Define size and position of panels

% Figure Size
figSize = [18.3 14.8]; % [width lenght]


% Panels Sizes

panelASize = [6.6 4];
panelBSize = [2 1.3];
panelC1Size = [2 1.3];
panelC2Size = [3 2.5];
panelDSize = [2 1.3];
panelE1Size = [1.8 1.1];
panelE2Size = [2.2 1.3];
panelE3Size = [2 1.3];


% panelD1Size = [2.5 2.5];
% panelD2Size = [2.5 2.5];
% panelD3Size = [2.5 2.5];
% panelD4Size = [2.5 2.5];
% panelE1Size = [2.5 2.5];
% panelE2Size = [2.5 2.5];
% panelE3Size = [2.5 2.5];

% Panels Positions

panelAPos = [1 figSize(2)-panelASize(2)-0.1];

panelB1Pos = [1.5+panelASize(1)   figSize(2)-panelBSize(2)-0.1];
panelB2Pos = [1.5+panelASize(1)   figSize(2)-panelBSize(2)*2-0.25];
panelB3Pos = [1.5+panelASize(1)   figSize(2)-panelBSize(2)*3-0.4];
panelB4Pos = [(1.5 + 0.5)+panelASize(1)+panelBSize(1)   figSize(2)-panelBSize(2)-0.1];
panelB5Pos = [(1.5 + 0.5)+panelASize(1)+panelBSize(1)   figSize(2)-panelBSize(2)*2-0.25];
panelB6Pos = [(1.5 + 0.5)+panelASize(1)+panelBSize(1)   figSize(2)-panelBSize(2)*3-0.4];

% panelB3Pos = [1*(2.5 + 1)+panelASize(1)+panelBSize(1)*2   figSize(2)-panelBSize(2)-0.1];
% panelB4Pos = [1*2.5+panelASize(1)   figSize(2)-panelBSize(2)*2-0.25];
% panelB5Pos = [1*(2.5 + 0.5)+panelASize(1)+panelBSize(1)   figSize(2)-panelBSize(2)*2-0.25];
% panelB6Pos = [1*(2.5 + 1)+panelASize(1)+panelBSize(1)*2   figSize(2)-panelBSize(2)*2-0.25];
% panelB7Pos = [1*2.5+panelASize(1)   figSize(2)-panelBSize(2)*3-0.4];
% panelB8Pos = [1*(2.5 + 0.5)+panelASize(1)+panelBSize(1)   figSize(2)-panelBSize(2)*3-0.4];
% panelB9Pos = [1*(2.5 + 1)+panelASize(1)+panelBSize(1)*2   figSize(2)-panelBSize(2)*3-0.4];

panelBPos = [panelB1Pos;...
    panelB2Pos;...
    panelB3Pos;...
    panelB4Pos;...
    panelB5Pos;...
    panelB6Pos];
panelBPos(:,2) = panelBPos(:,2) - 0.3;

panelC1Pos = [(1.5 + 1.8)+panelASize(1)+panelBSize(1)+panelBSize(1)   figSize(2)-panelC1Size(2)-0.1];
panelC2Pos = [(1.5 + 1.8)+panelASize(1)+panelBSize(1)+panelBSize(1)-0.5   figSize(2)-panelC1Size(2)-panelC2Size(2)-0.5];
panelCPos = [panelC1Pos;...
    panelC2Pos];
panelCPos(:,2) = panelCPos(:,2) - 0.3;





panelD1Pos = [1 figSize(2)-panelASize(2)-panelDSize(2)-0.1-1.5];
panelD2Pos = [1 figSize(2)-panelASize(2)-panelDSize(2)*2-0.1-1.5*1.1];
panelD3Pos = [1 figSize(2)-panelASize(2)-panelDSize(2)*3-0.1-1.5*1.18];
panelD4Pos = [1*1.5+panelDSize(1) figSize(2)-panelASize(2)-panelDSize(2)-0.1-1.5];
panelD5Pos = [1*1.5+panelDSize(1) figSize(2)-panelASize(2)-panelDSize(2)*2-0.1-1.5*1.1];
panelD6Pos = [1*1.5+panelDSize(1) figSize(2)-panelASize(2)-panelDSize(2)*3-0.1-1.5*1.18];
panelDPos = [panelD1Pos;...
    panelD2Pos;...
    panelD3Pos;...
    panelD4Pos;...
    panelD5Pos;...
    panelD6Pos];

panelDPos(:,2) = panelDPos(:,2) - 0.7;

% panelD1Pos = [1 figSize(2)-panelASize(2)-panelD1Size(2)-0.1-0.8];
% panelD2Pos = [1*1.5+panelD1Size(1) figSize(2)-panelASize(2)-panelD1Size(2)-0.1-0.8];
% panelD3Pos = [1 figSize(2)-panelASize(2)-panelD1Size(2)-panelD3Size(2)-0.1-0.8*1.5];
% panelD4Pos = [1*1.5+panelD3Size(1) figSize(2)-panelASize(2)-panelD1Size(2)-panelD3Size(2)-0.1-0.8*1.5];


panelE1Pos = [1*2.7+panelDSize(1)+panelE1Size(1) figSize(2)-panelASize(2)-panelE1Size(2)-0.1-1.2];
panelE2Pos = [1*2.1+panelDSize(1)+panelE2Size(1) figSize(2)-panelASize(2)-panelE1Size(2)-panelE2Size(2)-0.1-1.2*1.6];
panelE3Pos = [1*2.2+panelDSize(1)+panelE2Size(1) figSize(2)-panelASize(2)-panelE1Size(2)-panelE2Size(2)-panelE3Size(2)-0.1-1.2*1.9];
%panelE2Pos = [1*4+panelD3Size(1)+panelE1Size(1)+panelE2Size(1) figSize(2)-panelASize(2)-panelD1Size(2)-panelD3Size(2)-0.1-0.9];
% panelE3Pos = [1*5+panelD3Size(1)+panelE1Size(1)+panelE2Size(1)+panelE3Size(1) figSize(2)-panelASize(2)-panelD1Size(2)-panelD3Size(2)-0.1-0.9];
panelEPos = [panelE1Pos;...
             panelE2Pos;...
             panelE3Pos]; 

panelEPos(:,2) = panelEPos(:,2) - 0.7;         

% Panel letters positions

lAPos = [-0.8 panelASize(2)-0.25];
lBPos = [-0.8 panelBSize(1,2)-0.25];
lCPos = [-1 panelC1Size(1,2)-0.25];
lDPos = [-0.8 panelDSize(1,2)-0.25];
lEPos = [-0.8 panelE1Size(1,2)-0.25];


% Define colors to use

hitsColor = [0, 153, 136] / 255;
faColor = [204, 51, 17] / 255;
greyColor = [187, 187, 187] / 255;
lightGreyColor = [220, 220, 220] / 255;
soundColor = [238,51,119] / 255;


f2 = figure;
f2.Units = 'centimeters';
f2.Position = [20, 15, figSize(1), figSize(2)];

%% Place holder axes for Panel A, Recordings scheme and histology

hold on
plot(1:10,1:10,'visible','off')
axA = gca;
axA.Units = 'centimeters';
axA.XTick = [];
axA.YTick = [];
axA.Visible = 'off';
axA.Position = [panelAPos(1), panelAPos(2), panelASize(1), panelASize(2)];
hold on
lA = text(lAPos(1),lAPos(2),'a','Units','centimeters');
lA.FontSize = 9;
lA.FontWeight = 'bold';



%% Panel B, Space PSTHS

% Plot

cColors = [hitsColor;faColor];
conditions=[1 3];
soundStarts = [1 5 10 15 20 25 30 35 40 45 50] ;
offSetCounts = 11;
smoothing=0.15;

% :length(toPlot(:,1))

p = 1;
for n = 1:2
    for d = 1 : length(toPlot.B{n,2}(1,:))
        
        ssIdx = toPlot.B{n,2}(1,d);
        area = toPlot.B{n,3};
        sStartPos = soundStarts(ssIdx)*2+offSetCounts;
        sEndPos = soundStarts(ssIdx)*2+8*2+offSetCounts;
        aStartPos = area(1)*2+offSetCounts;
        aEndPos = area(2)*2+offSetCounts;
        
        axes
        for c = 1:length(conditions)
            
            FRs = toPlot.B{n,1}{1,d}{1,c};
            
            smtFrs = zeros(size(FRs));
            for t = 1:length(FRs(:,1))
                smtFRs(t,:) = smooth(FRs(t,:),smoothing,'lowess');
            end
            
            mFR = mean(smtFRs);
            eFR = std(smtFRs) / sqrt(length(smtFRs(:,1)));
            supErr = mFR + eFR;
            maxFrs(c)=max(supErr);
            infErr = mFR - eFR;
            xBins = toPlot.B{n,4}(1:end-1);
            
            xPatch = [xBins fliplr(xBins)];
            yPatch = [supErr fliplr(infErr)];
            
            frPatchB(p,c) = patch(xPatch,yPatch,cColors(c,:),'LineStyle','none');
            alpha(frPatchB(p,c),0.4)
            hold on
        end
        
        maxFR = round(max(maxFrs));
        yMax= round(maxFR + maxFR*0.1);
        
        ylim([0 yMax])
        xlim([1 xBins(end)])
        
        
        line([sStartPos sStartPos],[0 yMax],'Color',[soundColor 0.3]);
        line([sEndPos sEndPos],[0 yMax],'Color',[soundColor 0.3]);
        
        line([aStartPos aStartPos],[0 yMax],'Color',greyColor);
        line([aEndPos aEndPos],[0 yMax],'Color',greyColor);
        
        text(sStartPos+2,yMax*0.1,['cm ' num2str(soundStarts(toPlot.B{n,2}(d)))],'FontSize',6)
        
        
        if p == 2
            ylabel('Firing Rate (Hz)','Fontsize',6)
        end
        
        box on
        axB(p) = gca;
        
        if mod(p,3) ~= 0
            axB(p).XTick = [];
        else
            axB(p).XTick = [sStartPos+(sEndPos-sStartPos)/2 aStartPos+(aEndPos-aStartPos)/2];
            axB(p).XTickLabel = {'Sound','Area'};
        end
        
        axB(p).TickLength = [0 0];
        axB(p).YTick = [0 round(maxFR/2)  maxFR];
        axB(p).FontSize = 6;
        axB(p).Units = 'centimeters';
        axB(p).Position = [panelBPos(p,1), panelBPos(p,2), panelBSize(1), panelBSize(2)];
        
        if p == 1
            lB = text(lBPos(1),lBPos(2),'b','Units','centimeters');
            lB.FontSize = 9;
            lB.FontWeight = 'bold';
        end
        
        p = p+1;
    end
end

lgB = legend([frPatchB(1,1) frPatchB(1,2)],{'Stop Trials','No Stop Trials'},'Location','southoutside');
lgB.Units = 'centimeters';
lgB.Position = [panelBPos(3,1)+1.6, panelBPos(3,2) - 0.8, 0, 0.5];
lgB.Box = 'off';
lgB.Orientation = 'horizontal';
lgB.ItemTokenSize = [10,10];

%% Panel C, Conditions difference during delay

%%%PanelC1

cColors = [hitsColor;faColor];
conditions=[1 3];
smoothing=0.1;
nrSegments = 3; % Number of delay divisions to test for difference between conditions

pM = 1; % Mouse Number
pS = 1; % Session Number
pN = 9; % Neurons Number
pD = 6; % Delay Number

Sound = [25 33]*2+11;
Area = [100 115]*2+11;
bWidth = (((Area(1)-10) - (Sound(2)+10))/3) ;
b1 = Sound(2)+10;
b2 = b1 + bWidth;
b3 = b2 + bWidth;
b4 = b3 + bWidth;
bs = [b1 b2 b3 b4];

axes
for c = 1:length(conditions(1,:))
    hold on
    pl(c) = plot(smooth(toPlot.C{1,3}{1,pD}{1,conditions(c)}{pN,1},smoothing,'lowess'),'Color',[cColors(c,:) 0.6],'LineWidth',1.5);
    
    for b = 1:nrSegments
        
        line([bs(b) bs(b+1)],[mean(toPlot.C{1,2}{1,1}{1,conditions(c)}{1,b}) mean(toPlot.C{1,2}{1,1}{1,conditions(c)}{1,b})],'Color',cColors(c,:),'LineWidth',0.5)
        
        mFR = mean(toPlot.C{1,2}{1,1}{1,conditions(c)}{1,b});
        stdFR = std(toPlot.C{1,2}{1,1}{1,conditions(c)}{1,b})/sqrt(numel(toPlot.C{1,2}{1,1}{1,conditions(c)}{1,b}));
        errBCX = bs(b)+(bs(b+1) - bs(b))/2;
        errBC(c,b) = errorbar(errBCX,mFR,stdFR,'Color',cColors(c,:),'LineWidth',1.5);
        errBC(c,b).Marker = '.';
        errBC(c,b).MarkerSize = 5;
        errBC(c,b).MarkerFaceColor = cColors(c,:);
        errBC(c,b).MarkerEdgeColor = cColors(c,:);
        errBC(c,b).Color = 'k';
        
        text(bs(b)+8,35,['S',num2str(b)],'Fontsize',4)
        
    end
    maxsP(1,c) = mFR + stdFR;
end
maxP = round(max(maxsP));

line([Sound(1) Sound(1)],[0 40],'Color',[soundColor 0.3])
line([Sound(2) Sound(2)],[0 40],'Color',[soundColor 0.3])

line([Area(1)+4 Area(1)+4],[0 40],'Color',greyColor)
line([Area(2)+4 Area(2)+4],[0 40],'Color',greyColor)

ylim([0 40])
xlim([11 length(toPlot.C{1,3}{1,pD}{1,conditions(c)}{pN,1})-11])

ylabel('Firing Rate (Hz)','Fontsize',6)

box on
axC(1) = gca;


axC(1).XTick = [Sound(1)+(Sound(2)-Sound(1))/2 Area(1)+(Area(2)-Area(1))/2];
axC(1).XTickLabel = {'Sound','Area'};


axC(1).TickLength = [0 0];
axC(1).YTick = [0 round(maxP/2)  maxP];
axC(1).FontSize = 6;
axC(1).Units = 'centimeters';
axC(1).Position = [panelCPos(1,1), panelCPos(1,2), panelC1Size(1), panelC1Size(2)];

lC = text(lCPos(1),lCPos(2),'c','Units','centimeters');
lC.FontSize = 9;
lC.FontWeight = 'bold';

lgC1 = legend([pl(1)  pl(2)],{'STrials ','NSTrials'},'Location','eastoutside','Fontsize',5);
lgC1.Units = 'centimeters';
lgC1.Position = [panelCPos(1,1)+panelC1Size(1)+0.6, panelCPos(1,2)+0.4, 0, 0.5];
lgC1.Box = 'off';
lgC1.Orientation = 'vertical';
lgC1.ItemTokenSize = [5,5];




%%%PanelC2

rgbW = [245 245 245]/255;

cColors = [hitsColor;faColor];

Sound = [25 33];
Area = [100 114];
bWidth = ((Area(1)-5) - (Sound(2)+5))/3 ;
b1 = Sound(2)+5+bWidth/2;
b2 = b1 + bWidth;
b3 = b2 + bWidth;
allSig = vertcat(toPlot.C{1,1}{:});



% selectivity categories
patterns = [1 1 1;-1 -1 -1;1 1 0;-1 -1 0;0 1 1;0 -1 -1;1 0 1;-1 0 -1;1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1;0 0 0];

% preallocate for selectivity categories counter
patternMatches = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
a=1;

% group and count neurons in each selectivity category
for p = 1:length(patterns(:,1))
    for n = 1:length(allSig(:,1))
        if isequal(allSig(n,:),patterns(p,:))
            patternMatches(1,p) =  patternMatches(1,p)+1;
            newAllSig(a,:) = allSig(n,:);
            a=a+1;
        end
    end
end


neurPat = patternMatches(1:14);
neurPat = reshape(neurPat,2,length(neurPat)/2)';
neurPat = neurPat / length(allSig(:,1));
sumNeurPat = sum(neurPat,2);
textNeurPats = flipud(sumNeurPat);
propPat = neurPat ./ repmat(sumNeurPat,1,2);
propPat = reshape(propPat',length(propPat)*2,1);

schemeXBase = 0;
schemeYBase = 0.2;
barBase = 0.15 + 0.01;

Xcrds = {[0.1 0.15 0.15 0.1] + schemeXBase;...
    [0.05 0.1 0.1 0.05] + schemeXBase;...
    [0 0.05 0.05 0] + schemeXBase;...
    [0 0.05 0.05 0;0.1 0.15 0.15 0.1] + schemeXBase;...
    [0.05 0.15 0.15 0.05] + schemeXBase;...
    [0 0.1 0.1 0] + schemeXBase;...
    [0 0.15 0.15 0] + schemeXBase};


Ycrds = {[6 6 7 7] + schemeYBase;...
    [5 5 6 6];...
    [4 4 5 5];...
    [3 3 4 4;3 3 4 4] + schemeYBase;...
    [2 2 3 3] + schemeYBase;...
    [1 1 2 2] + schemeYBase;...
    [0 0 1 1] + schemeYBase};


axes

b = barh((0.5:6.5) + schemeYBase,neurPat + barBase,'BaseValue',barBase);
b(1).FaceColor = hitsColor;
b(1).EdgeColor = hitsColor;
b(2).FaceColor = faColor;
b(2).EdgeColor = faColor;

text(0.2,7.8,['Total: ',num2str(round(sum(sum(neurPat)),2))],'Fontsize',5)


hold on



for i = 1:length(Xcrds(:,1))
    for j = 1:length(Xcrds{i,1}(:,1))
        pX = Xcrds{i,1}(j,:);
        pY = Ycrds{i,1}(j,:);
        
        patch(pX,pY,greyColor,'EdgeColor','w','FaceColor',lightGreyColor)
        
        if length(Xcrds{i,1}(:,1)) == 1
            text((pX(1) + (pX(2) - pX(1))/2) - 0.018 ,(pY(2) + (pY(3) - pY(2))/2) - 0.08,num2str(round(textNeurPats(i),2)),'Fontsize',5)
        else
            text(0.075 - 0.018, (pY(2) + (pY(3) - pY(2))/2) - 0.08 ,num2str(round(sumNeurPat(i),2)),'Fontsize',5)
        end
        
        
    end
end

segDist = (Xcrds{end,1}(1,2) - Xcrds{end,1}(1,1)) / 6;
seg1 = Xcrds{end,1}(1,1) + segDist + 0.01;
seg2 = Xcrds{end,1}(1,1) + segDist*3 + 0.01;
seg3 = Xcrds{end,1}(1,1) + segDist*5 + 0.01;

text(seg1 - 0.02,7.8,'S1','Fontsize',5,'FontWeight','bold')
text(seg2 - 0.02,7.8,'S2','Fontsize',5,'FontWeight','bold')
text(seg3 - 0.02,7.8,'S3','Fontsize',5,'FontWeight','bold')


maxX = max(max(neurPat + barBase));
xlim([-0.01 maxX+0.01])
ylim([0 8.5])

xlC2 = xlabel('Prop Neurons','Fontsize',6);
xlC2.Units = 'centimeters';
xlC2.Position = [2.4,-0.4304];

box on
axC(2) = gca;
axC(2).FontSize = 6;
halfDist = (Xcrds{end,1}(1,2) - Xcrds{end,1}(1,1)) / 2;
axC(2).XTick = [Xcrds{end,1}(1,1) + halfDist [0 0.065 0.13]+barBase];
axC(2).XTickLabel ={'Delay Period', '0', '0.07','0.14'};
axC(2).YTick = [];
axC(2).TickLength = [0 0];
axC(2).Units = 'centimeters';
axC(2).Position = [panelCPos(2,1), panelCPos(2,2), panelC2Size(1), panelC2Size(2)];


lgC2 = legend([b(1)  b(2)],{'ST Preference','NST Preference'},'Location','southoutside','Fontsize',5);
lgC2.Units = 'centimeters';
lgC2.Position = [panelCPos(2,1)+0.9, panelCPos(2,2) - 1.2, 0, 0.5];
lgC2.Box = 'off';
lgC2.Orientation = 'vertical';
lgC2.ItemTokenSize = [5,5];





%% Panel D, Time Rasters and PSTHS

cColors = [hitsColor;faColor];
conditions=[1 3];
smoothing=0.1;
Events={'Sound Start','Area Start'};
neur = 1:length(toPlot.D(:,1));


p = 1;

for n = 1:2
    
    
    % Plot Rasters
    
    for s=1:length(toPlot.D{n,5}(1,:))
        spikes{1,conditions(s)}=toPlot.D{n,4}{s,1};
        ticks{1,conditions(s)}=0:length(spikes{conditions(s)});
        rasterMaxs(1,s)=toPlot.D{n,5}(1,s);
    end
    
    
    for c = 1:length(conditions)
        axes
        for s=2:length(spikes{conditions(c)})
            for j=1:length(spikes{conditions(c)}{s})
                line([spikes{conditions(c)}{s}(j) spikes{conditions(c)}{s}(j)],[ticks{conditions(c)}(s-1) ticks{conditions(c)}(s)],'Color',cColors(c,:))
            end
        end
        
        rasterMax=rasterMaxs(1,1);
        line([0 0],[0  rasterMax],'Color',greyColor,'LineStyle',':','LineWidth',1);
        ylim([0 rasterMax])
        xlim([-1 1])
        
        ylabel('Trials','Fontsize',6)
        
        if c == 1
            title(Events{1,n},'Fontsize',6)
        end
        
        
        
        box on
        axD(p) = gca;
        
        axD(p).XTick = [];
        axD(p).YTick = [];
        
        %         axC(p).XTick = [sStartPos+(sEndPos-sStartPos)/2 aStartPos+(aEndPos-aStartPos)/2];
        %         axC(p).XTickLabel = {'Sound','Area'};
        %     end
        
        axD(p).TickLength = [0 0];
        %axC(p).YTick = [0 round(maxFR/2)  maxFR];
        axD(p).FontSize = 6;
        axD(p).Units = 'centimeters';
        axD(p).Position = [panelDPos(p,1), panelDPos(p,2), panelDSize(1), panelDSize(2)];
        
        
        if p == 1
            lD = text(lDPos(1),lDPos(2),'d','Units','centimeters');
            lD.FontSize = 9;
            lD.FontWeight = 'bold';
        end
        
        
        p = p+1;
    end
    
    
    % Plot PSTHS
    
    axes
    for c = 1:length(conditions)
        FRs = toPlot.D{n,c};
        
        smtFRs = zeros(size(FRs));
        for t = 1:length(FRs(:,1))
            smtFRs(t,:) = smooth(FRs(t,:),smoothing,'lowess');
        end
        
        mFR = mean(smtFRs);
        eFR = std(smtFRs) / sqrt(length(smtFRs(:,1)));
        supErr = mFR + eFR;
        maxFrs(c)=max(supErr);
        infErr = mFR - eFR;
        xBins = toPlot.D{n,6};
        
        xPatch = [xBins fliplr(xBins)];
        yPatch = [supErr fliplr(infErr)];
        
        frPatchD(p,c) = patch(xPatch,yPatch,cColors(c,:),'LineStyle','none');
        alpha(frPatchD(p,c),0.4)
        hold on
        
        
    end
    
    maxFR = round(max(maxFrs));
    yMax= round(maxFR + maxFR*0.1);
    
    line([0 0],[0  yMax],'Color',greyColor,'LineStyle',':','LineWidth',1);
    
    ylim([0 yMax])
    xlim([-1.1 1.1])
    
    if n == 1
        ylabel('Firing Rate (Hz)','Fontsize',6)
    end
    
    xlabel('Time (s)','Fontsize',6)
    
    
    box on
    axD(p) = gca;
    
    axD(p).XTick = [-1 0 1];
    axD(p).YTick = [ 0 round(maxFR/2) maxFR];
    
    
    %         axC(p).XTick = [sStartPos+(sEndPos-sStartPos)/2 aStartPos+(aEndPos-aStartPos)/2];
    %         axC(p).XTickLabel = {'Sound','Area'};
    %     end
    
    axD(p).TickLength = [0 0];
    %axC(p).YTick = [0 round(maxFR/2)  maxFR];
    axD(p).FontSize = 6;
    axD(p).Units = 'centimeters';
    axD(p).Position = [panelDPos(p,1), panelDPos(p,2), panelDSize(1), panelDSize(2)];
          
    p = p+1;
end

lgD = legend([frPatchD(3,1) frPatchD(3,2)],{'Stop Trials','No Stop Trials'},'Location','southoutside');
lgD.Units = 'centimeters';
lgD.Position = [panelDPos(3,1)+1.6, panelDPos(3,2) - 1.1, 0, 0.5];
lgD.Box = 'off';
lgD.Orientation = 'horizontal';
lgD.ItemTokenSize = [10,10];




%% Panel E, Neurons FR modulation by task events

%%%PanelE1, example PSTH

cColors = [hitsColor;faColor];
conditions=[1 3];
smoothing=0.05;
Events={'Sound Start','Area Start'};

% Calculate FRs for segments being compared

before=[0 0.15]; % Segment of time before the event
after=[0.05 0.2]; % Segment of time after the event
xBins = toPlot.E{1,2}{1,5};
allTrials=toPlot.E{1,2}{1,1}; 

% get indexes of segment edges

[~,minBefore] = min(abs(xBins-(before(1)-before(2))));
[~,maxBefore] = min(abs(xBins-before(1)));
[~,minAfter] = min(abs(xBins-after(1)));
[~,maxAfter] = min(abs(xBins-(after(1)+after(2))));


segBefore=allTrials(:,minBefore:maxBefore-1);
segAfter=allTrials(:,minAfter:maxAfter-1);

meanBefore = mean(mean(segBefore,2));
errBefore = std(mean(segBefore,2)) / sqrt(length(segBefore(:,1)));
meanAfter = mean(mean(segAfter,2));
errAfter = std(mean(segAfter,2)) / sqrt(length(segAfter(:,1)));

mSegs = [meanBefore; meanAfter];
eSegs = [errBefore; errAfter];

segXEdges = [before(2)*-1 0;after(1) after(2)];
segCenters = [segXEdges(1,1) + (segXEdges(1,2) -  segXEdges(1,1))/2 ; segXEdges(2,1) + (segXEdges(2,2) -  segXEdges(2,1))/2];

% Plot PSTHS

FRs = toPlot.E{1,2}{1,1};

smtFRs = zeros(size(FRs));
for t = 1:length(FRs(:,1))
    smtFRs(t,:) = smooth(FRs(t,:),smoothing,'lowess');
end

mFR = mean(smtFRs);
eFR = std(smtFRs) / sqrt(length(smtFRs(:,1)));
supErr = mFR + eFR;
maxFrs=max(supErr);
infErr = mFR - eFR;
minFrs=min(infErr);
xBins = toPlot.E{1,2}{1,5};

axes
hold on

xPatch = [xBins fliplr(xBins)];
yPatch = [supErr fliplr(infErr)];

frPatchD = patch(xPatch,yPatch,cColors(1,:),'LineStyle','none');
alpha(frPatchD,0.4)

for e = 1:length(mSegs(:,1))
    
    line(segXEdges(e,:),[mSegs(e) mSegs(e)],'Color',greyColor,'LineWidth',0.5)
    
    errBE(e) = errorbar(segCenters(e),mSegs(e),eSegs(e),'Color',cColors(1,:),'LineWidth',1.5);
    errBE(e).Marker = '.';
    errBE(e).MarkerSize = 5;
    errBE(e).MarkerFaceColor = cColors(1,:);
    errBE(e).MarkerEdgeColor = cColors(1,:);
    errBE(e).Color = 'k';
    
end


maxFR = round(max(maxFrs));
minFR = round(min(minFrs));

yMax= round(maxFR + maxFR*0.1);
yMin= round(minFR - minFR*0.1);


line([0 0],[yMin  yMax],'Color',greyColor,'LineStyle',':','LineWidth',1);


ylim([yMin+4 yMax])
xlim([-0.32 0.32])


ylabel('FR (Hz)','Fontsize',6)


xlabel('Time (s)','Fontsize',6)

box on
axE(1) = gca;

axE(1).XTick = [-0.3 0 0.3];
axE(1).YTick = [];


%         axC(p).XTick = [sStartPos+(sEndPos-sStartPos)/2 aStartPos+(aEndPos-aStartPos)/2];
%         axC(p).XTickLabel = {'Sound','Area'};
%     end

axE(1).TickLength = [0 0];
%axC(p).YTick = [0 round(maxFR/2)  maxFR];
axE(1).FontSize = 6;
axE(1).Units = 'centimeters';
axE(1).Position = [panelEPos(1,1), panelEPos(1,2), panelE1Size(1), panelE1Size(2)];


lE = text(lEPos(1),lEPos(2),'e','Units','centimeters');
lE.FontSize = 9;
lE.FontWeight = 'bold';



%%%PanelE2, Response per event

allUnits = toPlot.E{1,1};

stopEx = sum(allUnits{1,1} > 0);
stopInh = sum(allUnits{1,1} < 0);
noStopEx = sum(allUnits{1,2} > 0);
noStopEx = [noStopEx(1:4) 0 noStopEx(5)];
noStopInh = sum(allUnits{1,2} < 0);
noStopInh = [noStopInh(1:4) 0 noStopInh(5)];
excited = [stopEx' noStopEx' ]/length(allUnits{1,1}(:,1));
inhibited = [stopInh' noStopInh' ]/length(allUnits{1,1}(:,1));


axes
hold on

bEx = barh(excited);
bEx(1).FaceColor = hitsColor;
bEx(1).EdgeColor = hitsColor;
bEx(2).FaceColor = faColor;
bEx(2).EdgeColor = faColor;

bInh = barh(inhibited*-1);
bInh(1).FaceColor = hitsColor;
bInh(1).EdgeColor = hitsColor;
bInh(2).FaceColor = faColor;
bInh(2).EdgeColor = faColor;

text(-0.16,7.5,'Inhibited','Fontsize',5);
text(0.04,7.5,'Excited','Fontsize',5);

ylim([0 8.5])
xlim([-0.2 0.2])


yTL = {'TON' 'SON' 'SOFF' 'AON' 'RWD' 'AOFF'};
yTLYs = [1 2 3 4 5 6];
yTLX = -0.19;

for i = yTLYs;
    text(yTLX,i,yTL{i},'Fontsize',4,'FontWeight','bold');
end


box on
axE(2) = gca;
axE(2).FontSize = 6;
axE(2).XTick = [-0.14 0 0.14];
axE(2).XTickLabel = {'0.14', '0', '0.14'};
axE(2).YTick = [];
axE(2).TickLength = [0 0];
axE(2).Units = 'centimeters';
axE(2).Position = [panelEPos(2,1), panelEPos(2,2), panelE2Size(1), panelE2Size(2)];

%%%PanelE3, #events a neuron responds to


allResponsesStop = sum(abs(allUnits{1,1}),2);
allResponsesNoStop = sum(abs(allUnits{1,2}),2); 
 
sCounts = histcounts(allResponsesStop(allResponsesStop>0));
nsCounts = histcounts( allResponsesNoStop( allResponsesNoStop>0));

totProps = [sCounts' nsCounts'] / length(allUnits{1,1}(:,1));

axes
hold on

bE3 = barh(totProps);
bE3(1).FaceColor = hitsColor;
bE3(1).EdgeColor = hitsColor;
bE3(2).FaceColor = faColor;
bE3(2).EdgeColor = faColor;

text(0.05,5,['Tot Mod Neur: ', num2str(0.67)],'Fontsize',5);


xlim([0 0.3])
ylim([0 6])

ylabel('#Events','Fontsize',6)
xlabel('Prop Neurons','Fontsize',6)
 
box on
axE(3) = gca;
axE(3).FontSize = 6;
axE(3).XTick = [0 0.15 0.3];
axE(3).TickLength = [0 0];
axE(3).Units = 'centimeters';
axE(3).Position = [panelEPos(3,1), panelEPos(3,2), panelE3Size(1), panelE3Size(2)];


disp('Saving fig')
tic
filename = '/Volumes/GoogleDrive/O meu disco/paper/RawFig';
print(f2,filename,'-depsc','-r300')
disp('Finished saving fig')
toc
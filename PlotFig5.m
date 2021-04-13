clearvars -except toPlot
close all


%%% Script to plot the different panels of the paper fifth figure


%% Load and process needed Data

load('/Volumes/Cerebro/Recording Data/SessionsMap.mat')


if ~exist('toPlot','var')
    
    %%% Load panel A data
    
    sessions = [1 2 13 14]; % Sessions to load data from
    stdv = 0.14; % Stdv of kernel used to smooth data
    def = 0.5; % Parameter of allowed deformation
    tSelectMethod = 1; % Method used to select similar trials: 0 - manual; 1 - based on correlations; 2 - based on vector similarity
    intDels = [3 5 7 9 11]; % Delays on which the analysis was run
    intCond = 1; % Trial's condition: 1-hits; 2-misses; 3-correct rejections; 4-false alarms
    p1 = 0.25; %Cut trials p1 seconds after the higher area time across trials
    p2 = 0.5;  %Cut trials p2 seconds before the lowest moving after stopping time across trials
    
    
    for s = 1:length(sessions)
        for d = 1:length(intDels)
            
            disp(['loading delay ', num2str(intDels(d)),' from session ', num2str(sessions(s))])
            [Name] = SessionNumber2Name(sessions(s));
            load(['/Volumes/Cerebro/Recording Data/',Name,'/Predict from Neurons/SpdPredLagsBootCheck Trials tsMethod_',num2str(tSelectMethod),' Cond_',num2str(intCond),' Delay_',num2str(intDels(d)),' P1_',num2str(p1),' P2_',num2str(p2),' Sigma_',num2str(stdv),' Deform_',num2str(def),'.mat'],'testData','origData');
            
            
            allOrig{s,d} = origData;
            allTest{s,d} = testData;
            
            clear origData testData
        end
    end
    toPlot.A{1,1} = allOrig;
    toPlot.A{1,2} = allTest;
    clear allOrig allTest
    
end

%%% Define size and position of figure and panels

figSize = [18.3 14.8];


panelASize = [3.5 2.5];
% panelBSize = [2 2];
% panelCSize = [3.5 2.5];
% panelDSize = [2 2];
% panelESize = [2 2];
% panelGSize = [1.8 1.8];
% panelFSize = [2.8 1.2];


panelA1Pos = [1 figSize(2)-panelASize(2)-0.5];
% panelB1Pos = [1+panelA1Size(1)+1.5 figSize(2)-panelBSize(2)-0.7];
% panelB2Pos = [1+panelA1Size(1)+panelBSize(1)+2.5 figSize(2)-panelBSize(2)-0.7];
% panelCPos = [1+panelA1Size(1)+panelBSize(1)*2+4 figSize(2)-panelCSize(2)-0.6];
% panelD1Pos = [1 figSize(2)-panelA1Size(2)-panelDSize(2)-2.2];
% panelEPos = [1 figSize(2)-panelA1Size(2)-panelDSize(2)-panelESize(2)-3];
% panelGPos = [1+panelDSize(1)+panelFSize(1)+3.6 figSize(2)-panelA1Size(2)-panelGSize(2)-2.2];
% panelF1Pos = [1+panelDSize(1)+2 figSize(2)-panelA1Size(2)-panelFSize(2)-2.2];
% panelF2Pos = [1+panelDSize(1)+2 figSize(2)-panelA1Size(2)-panelFSize(2)*2-2.25];
% panelF3Pos = [1+panelDSize(1)+2 figSize(2)-panelA1Size(2)-panelFSize(2)*3-2.3];
% panelF4Pos = [1+panelDSize(1)+2 figSize(2)-panelA1Size(2)-panelFSize(2)*4-2.35];


lAPos = [-0.8 panelASize(2)-0.25];
% lBPos = [-1 panelBSize(2)-0.15];
% lCPos = [-0.8 panelCSize(2)-0.25];
% lDPos = [-0.8 panelDSize(2)-0.25];
% lEPos = [-0.8 panelESize(2)-0.25];
% lFPos = [-0.8 panelFSize(2)-0.25];
% lGPos = [-1.2 panelGSize(2)-0.25];


% Define colors to use

spdColor = [0, 119, 187] / 255;
hitsColor = [0, 153, 136] / 255;
faColor = [204, 51, 17] / 255;
greyColor = [187, 187, 187] / 255;
lightGreyColor = [220, 220, 220] / 255;
soundColor = [238,51,119] / 255;




% f3 = figure;
% f3.Units = 'centimeters';
% f3.Position = [20, 15, figSize(1), figSize(2)];


s = 3;
d = 5;

errLag = toPlot.A{1,2}{s,d}.RMSETest{1,1};
errNoLag = toPlot.A{1,2}{s,d}.RMSETest{1,2};
diffErr = errLag - errNoLag;
%[~,it] = max(diffErr);
it = 41;
testTrials = toPlot.A{1,2}{s,d}.testTrials{1,2}{1,it};
spds = [toPlot.A{1,1}{s,d}.Speeds{testTrials(1)} ; toPlot.A{1,1}{s,d}.Speeds{testTrials(2)}];



predsNoLag = toPlot.A{1,2}{s,d}.Preds{1,1}{1,it};
predsLag = toPlot.A{1,2}{s,d}.Preds{1,2}{1,it};

plot(spds,'b')
hold on
plot(predsNoLag,'r')
hold on
plot(predsLag,'g')







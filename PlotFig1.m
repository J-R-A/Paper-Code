clearvars -except allValidLaps lapsToPlot
close all


%%% Script to plot the different panels of the paper first figure

% Load performances and speeds data 

load('/Volumes/GoogleDrive/O meu disco/PhD Thesis/Chapter1/Performances.mat')

if ~exist('allValidLaps','var')
disp('Loading Heavy Stuff')    
load('/Volumes/GoogleDrive/O meu disco/PhD Thesis/Chapter1/allValidLaps.mat');
load('/Volumes/GoogleDrive/O meu disco/PhD Thesis/Chapter1/lapsToPlot.mat')
end


%%% Define size and position of panels

% Figure Size
figSize = [18.3 14.8]; % [width lenght]


% Panels Sizes

panelASize = [7.4 4];
panelBSize = [3 3];
panelCSize = [3 3];
panelD1Size = [2.5 2.5];
panelD2Size = [2.5 2.5];
panelD3Size = [2.5 2.5];
panelD4Size = [2.5 2.5];
panelE1Size = [2.5 2.5];
panelE2Size = [2.5 2.5];
panelE3Size = [2.5 2.5];

% Panels Positions

panelAPos = [1 figSize(2)-panelASize(2)-0.1];
panelBPos = [1*2.5+panelASize(1)   figSize(2)-panelBSize(2)-0.1];
panelCPos = [1*3.5+panelASize(1)+panelBSize(1) figSize(2)-panelCSize(2)-0.1];
panelD1Pos = [1 figSize(2)-panelASize(2)-panelD1Size(2)-0.1-0.8];
panelD2Pos = [1*1.5+panelD1Size(1) figSize(2)-panelASize(2)-panelD1Size(2)-0.1-0.8];
panelD3Pos = [1 figSize(2)-panelASize(2)-panelD1Size(2)-panelD3Size(2)-0.1-0.8*1.5];
panelD4Pos = [1*1.5+panelD3Size(1) figSize(2)-panelASize(2)-panelD1Size(2)-panelD3Size(2)-0.1-0.8*1.5];
panelE1Pos = [1*3+panelD3Size(1)+panelE1Size(1) figSize(2)-panelASize(2)-panelD1Size(2)-panelD3Size(2)-0.1-0.9];
panelE2Pos = [1*4+panelD3Size(1)+panelE1Size(1)+panelE2Size(1) figSize(2)-panelASize(2)-panelD1Size(2)-panelD3Size(2)-0.1-0.9];
panelE3Pos = [1*5+panelD3Size(1)+panelE1Size(1)+panelE2Size(1)+panelE3Size(1) figSize(2)-panelASize(2)-panelD1Size(2)-panelD3Size(2)-0.1-0.9];



% Panel letters positions

lAPos = [-0.8 panelASize(2)-0.25];
lBPos = [-0.8 panelBSize(2)-0.25];
lCPos = [-0.8 panelCSize(2)-0.25];
lDPos = [-0.8 panelD1Size(2)-0.25];
lEPos = [-0.8 panelE1Size(2)*2-0.25];






% Define colors to use

hitsColor = [0, 153, 136] / 255;
faColor = [204, 51, 17] / 255;
greyColor = [187, 187, 187] / 255;
soundColor = [238,51,119] / 255;


p1 = figure;
p1.Units = 'centimeters';
p1.Position = [20, 15, figSize(1), figSize(2)];

%% Place holder axes for Panel A, task scheme

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



%% Panel B, performance in specific delays

intDels = [3 5 7 9 11];
sSLocs = {'10','20','30','40','50'};
jitters = linspace(-0.2,0.2,8);


axes
for d = 1:length(intDels)
    for m = 1:length(delHits)
        
        % Plot Hit Rate
        mouseHitsMean = mean(delHits{1,m}(:,intDels(d)));
        miceDelHits(:,d) = delHits{1,m}(:,intDels(d));
        mouseHitsError = std(delHits{1,m}(:,intDels(d))) / sqrt(length(delHits{1,m}(:,intDels(d))));
        hitsE = errorbar(d + jitters(m),mouseHitsMean,mouseHitsError);
        hitsE.Marker = '.';
        hitsE.MarkerSize = 9;
        hitsE.MarkerFaceColor = hitsColor;
        hitsE.MarkerEdgeColor = hitsColor;
        hitsE.Color = 'k';
        hold on
         
        
        % Plot FA Rate
        mouseFAMean = mean(delFalseAlarms{1,m}(:,intDels(d)));
        miceDelFA(:,d) = delFalseAlarms{1,m}(:,intDels(d));
        mouseFAError = std(delFalseAlarms{1,m}(:,intDels(d))) / sqrt(length(delFalseAlarms{1,m}(:,intDels(d))));
        faE = errorbar(d + jitters(m),mouseFAMean,mouseFAError);
        faE.Marker = '.';
        faE.MarkerSize = 9;
        faE.MarkerFaceColor = faColor;
        faE.MarkerEdgeColor = faColor;
        faE.Color = 'k';
        hold on
    end
%     miceHitsMean = mean(mean(miceDelHits));
%     line([d+jitters(1)-0.05 d+jitters(end)+0.05],[miceHitsMean miceHitsMean],'Color',hitsColor)
%     
%     miceFAMean = mean(mean(miceDelFA));
%     line([d+jitters(1) d+jitters(end)],[miceFAMean miceFAMean],'Color',faColor)
    ylim([-0.1 1.1])
end
xlim([0.5 5.5])
line([0.5 5.5],[0.5 0.5],'Color',greyColor,'LineStyle','--','LineWidth',3)


axB = gca;
axB.XTick = [1 2 3 4 5];
axB.XTickLabel = sSLocs;
axB.TickLength = [0 0];
axB.YTick = [0  0.5 1];
axB.FontSize = 6;
axB.Units = 'centimeters';
axB.Position = [panelBPos(1), panelBPos(2), panelBSize(1), panelBSize(2)]; 
xlabel('Sound Start (cm)','Fontsize',6)
ylabel('Proportion of Trials','Fontsize',6)
lgB = legend({'Hit Rate','FA Rate'},'Location','best');
lgB.Units = 'centimeters';
lgB.Position = [panelBPos(1) + 1.2, panelBPos(2) - 1.1, 0, 0.5];
lgB.Orientation = 'horizontal';
lgB.ItemTokenSize = [10,10];
legend('boxoff')

lB = text(lBPos(1),lBPos(2),'b','Units','centimeters');
lB.FontSize = 9;
lB.FontWeight = 'bold';


%% Panel C, Duration of delay periods
intDels = [3 5 7 9 11];
sSLocs = {'10','20','30','40','50'};
jitters = linspace(-0.3,0.3,8);
yReSet = 4.5;

axes
for d = 1:length(intDels)
    for m = 1:length(rwddDelTimes)
        
        % Plot Hit trials times
        mouseHitsTimes = cell2mat(rwddDelTimes{1,m}(:,intDels(d)));
        mouseHitsTimes50(m) = prctile(cell2mat(rwddDelTimes{1,m}(:,intDels(d))),50);
        mouseHitsTimes25 = prctile(cell2mat(rwddDelTimes{1,m}(:,intDels(d))),25);
        mouseHitsTimes75 = prctile(cell2mat(rwddDelTimes{1,m}(:,intDels(d))),75);
        line([d + jitters(m) d + jitters(m)],[mouseHitsTimes25 mouseHitsTimes75],'color',hitsColor,'LineWidth',2)
        hold on
        
      
       mouseFATimes = cell2mat(rwddDelTimes{1,m}(:,intDels(d)));
       mouseFATimes50(m) = prctile(cell2mat(nRwddDelTimes{1,m}(:,intDels(d))),50);
       mouseFATimes25 = prctile(cell2mat(nRwddDelTimes{1,m}(:,intDels(d))),25);
       mouseFATimes75 = prctile(cell2mat(nRwddDelTimes{1,m}(:,intDels(d))),75);
       line([d + jitters(m) d + jitters(m)],[mouseFATimes25 mouseFATimes75] + yReSet,'color',faColor,'LineWidth',2)
       hold on



    end
    patch([d-0.5, d+0.5, d+0.5, d-0.5],[0, 0, median(mouseHitsTimes50), median(mouseHitsTimes50)],hitsColor,'FaceColor','none','EdgeColor',hitsColor,'EdgeAlpha',.4,'LineWidth',2)
    patch([d-0.5, d+0.5, d+0.5, d-0.5],[0, 0, median(mouseFATimes50), median(mouseFATimes50)] + yReSet,faColor,'FaceColor','none','EdgeColor',faColor,'EdgeAlpha',.4,'LineWidth',2)
end

box on 
ylim([0 9])
line([0 6],[yReSet yReSet],'Color','k')
text(3,8,'No Stop Trials','Fontsize',6,'Color',faColor)
text(3.5,3.5,'Stop Trials','Fontsize',6,'Color',hitsColor)


axC = gca;
axC.XTick = [1 2 3 4 5];
axC.XTickLabel = sSLocs;
axC.TickLength = [0 0];
axC.YTick = [0  2  4  4.5  6.5  8.5];
axC.YTickLabel = {'0','2','4','0','2','4'};
axC.FontSize = 6;
axC.Units = 'centimeters';
axC.Position = [panelCPos(1), panelCPos(2), panelCSize(1), panelCSize(2)];
xlabel('Sound Start (cm)','Fontsize',6)
ylabel('Time (s)','Fontsize',6)


lC = text(lBPos(1),lCPos(2),'c','Units','centimeters');
lC.FontSize = 9;
lC.FontWeight = 'bold';


%% Panel D, Speed Profiles   

%%% mean and single trials in example delay of example session

intDels = [3 5 7 9 11];
sSLocs = {'10','20','30','40','50'};
jitters = linspace(-0.2,0.2,8);


% Numbers of the animals that correpond to each row of allValidLaps, the order is the same 
miceNumbers = [59,61,67,69,71,73,74,76];
mToPlt = 8; % mouse to plot 
ssToPlt = [10 20 30 40 50]; % Sound starts to plot

if miceNumbers(1,mToPlt) == 59 || miceNumbers(1,mToPlt) == 61
    rwddSound = 1;
    nRwddSound = 2;
else
    rwddSound = 1+mod(miceNumbers(1,mToPlt),2);
    nRwddSound = 2-mod(miceNumbers(1,mToPlt),2);
end

sToPlt = 10; % session to plot 
sess1 = allValidLaps{mToPlt,sToPlt}; % data from session to plot
toPlot1 = lapsToPlot{mToPlt,sToPlt}; % altered speed data from session to plot
sLapPos = [1:length(toPlot1(1,:))]*0.5; % Values in cm for x axis  
areaStart =round(sess1{3,1}(1,sess1{2,1}/2048 == sess1{5,1}(1,1))); % Area start position
areaEnd =  round(sess1{3,1}(1,sess1{2,1}/2048 == sess1{5,1}(1,2))); % Area End position
intDel = 3; %

axes
hold on
rwdLaps = toPlot1(cell2mat(sess1(9,:)) == ssToPlt(1,intDel) & cell2mat(sess1(7,:)) == rwddSound & cell2mat(sess1(8,:)) == 1,:); 
nRwdLaps = toPlot1(cell2mat(sess1(9,:)) == ssToPlt(1,intDel) & cell2mat(sess1(7,:)) == nRwddSound & cell2mat(sess1(8,:)) == 1,:);
rwdP = plot(sLapPos,rwdLaps','Color',[hitsColor 0.2],'LineWidth',2);
rwdPM = plot(sLapPos,mean(rwdLaps),'Color',hitsColor,'LineWidth',4);
nRwdP = plot(sLapPos,nRwdLaps','Color',[faColor 0.2],'LineWidth',2);
nRwdPM = plot(sLapPos,mean(nRwdLaps),'Color',faColor,'LineWidth',4);
line([ssToPlt(1,intDel) ssToPlt(1,intDel)],[0 95],'Color',[soundColor 0.3],'LineStyle','--','LineWidth',3)
line([ssToPlt(1,intDel)+8 ssToPlt(1,intDel)+8],[0 95],'Color',[soundColor 0.3],'LineStyle','--','LineWidth',3)
line([areaStart areaStart],[0 95],'Color',greyColor,'LineWidth',3)
line([areaEnd areaEnd],[0 95],'Color',greyColor,'LineWidth',3)
text(65,5,['bM' num2str(76)],'FontSize',6,'FontWeight','bold')
xlim([1 sLapPos(end)])
yMax = round(max(max([rwdLaps' nRwdLaps'])));
ylim([0 yMax])
ylD = ylabel('Speed (cm/s)','Fontsize',6);
ylD.Units = 'centimeters';
ylD.Position = [lDPos(1)+0.25,lDPos(2)-2.5];

box on 
axD1 = gca;
axD1.XTick = [ssToPlt(1,intDel)+4 areaStart+(areaEnd - areaStart)/2];
axD1.XTickLabel = {'Sound','Area'};
axD1.TickLength = [0 0];
axD1.YTick = [0 yMax/2  yMax];
axD1.FontSize = 6;
axD1.Units = 'centimeters';
axD1.Position = [panelD1Pos(1), panelD1Pos(2), panelD1Size(1), panelD1Size(2)];

lD1 = text(lDPos(1),lDPos(2),'d','Units','centimeters');
lD1.FontSize = 9;
lD1.FontWeight = 'bold';

%%% means and SEMs in specific delays od specific sessions



% organize speed info to plot 

intMice = [4,6,8];
ssToPlt = [10 20 30 40 50]; % Sound starts to plot

for m = 1:length(intMice) 
    for s = 1:length(lapsToPlot(1,:))
        for id = 1:length(ssToPlt)
            sess3 = allValidLaps{intMice(m),s};
            toPlot3 = lapsToPlot{intMice(m),s};
            
            if miceNumbers(1,intMice(m)) == 59 || miceNumbers(1,intMice(m)) == 61
                rwddSound = 1;
                nRwddSound = 2;
            else
                rwddSound = 1+mod(miceNumbers(1,m),2);
                nRwddSound = 2-mod(miceNumbers(1,m),2);
            end
            
            
            
            % Correct rewarded trial laps
            allSpdsRwdd2C{1,m}{s,id} = toPlot3(cell2mat(sess3(9,:)) == ssToPlt(1,id) & cell2mat(sess3(7,:)) == rwddSound & cell2mat(sess3(8,:)) == 1,:);
            % Wrong non rewarded trial laps
            allSpdsNRwdd2C{1,m}{s,id} = toPlot3(cell2mat(sess3(9,:)) == ssToPlt(1,id) & cell2mat(sess3(7,:)) == nRwddSound & cell2mat(sess3(8,:)) == 1,:);
            % Correct rewarded trial laps
            allSpdsRwdd2W{1,m}{s,id} = toPlot3(cell2mat(sess3(9,:)) == ssToPlt(1,id) & cell2mat(sess3(7,:)) == rwddSound & cell2mat(sess3(8,:)) == 0,:);
            % Wrong non rewarded trial laps
            allSpdsNRwdd2W{1,m}{s,id} = toPlot3(cell2mat(sess3(9,:)) == ssToPlt(1,id) & cell2mat(sess3(7,:)) == nRwddSound & cell2mat(sess3(8,:)) == 2,:);
            
            
            
        end
        areaStarts(1,s) = round(sess3{3,1}(1,sess3{2,1}/2048 == sess3{5,1}(1,1)));
        areaEnds(1,s) =  round(sess3{3,1}(1,sess3{2,1}/2048 == sess3{5,1}(1,2)));
    end
    lapLens2{1,m} = cellfun('length',allSpdsRwdd2C{1,m});
    areaStart(1,m) = round(mean(areaStarts));
    areaEnd(1,m) = round(mean(areaEnds));
end
    


for m = 1:length(allSpdsRwdd2C(1,:))
    for ss = 1:length(allSpdsRwdd2C{1,m}(1,:))
        for s = 1:length(allSpdsRwdd2C{1,m}(:,ss))
            cutRwddC2{1,m}{s,ss} = allSpdsRwdd2C{1,m}{s,ss}(:,1:min(lapLens2{1,m}(:,ss)));
            meanRwddC2{1,m}{ss}(s,:) = mean(cutRwddC2{1,m}{s,ss});
            cutNRwddC2{1,m}{s,ss} = allSpdsNRwdd2C{1,m}{s,ss}(:,1:min(lapLens2{1,m}(:,ss)));
            meanNRwddC2{1,m}{ss}(s,:) = mean(cutNRwddC2{1,m}{s,ss});
            allXLapsPos2{1,m}{1,ss} = (1:min(lapLens2{1,m}(:,ss)))*0.5;
        end
    end
end



% mice 4

m = 1;
plotMean = 1;
axes
hold on

for l = 1:length(ssToPlt(1,:))
    line([ssToPlt(1,l) ssToPlt(1,l)],[0 95],'Color',[soundColor 0.3],'LineStyle','--','LineWidth',3);
end

line([areaStart(1,m) areaStart(1,m)],[0 95],'Color',greyColor,'LineWidth',3)
line([areaEnd(1,m) areaEnd(1,m)],[0 95],'Color',greyColor,'LineWidth',3)


for ss = 1:length(cutRwddC2{1,m}(1,:))
    
    xLapsPos2 = allXLapsPos2{1,m}{1,ss};
    
    if plotMean == 0
        allSpdsRwddCMat2 = vertcat(cutRwddC2{1,m}{:,ss});
        allSpdsNRwddCMat2 = vertcat(cutNRwddC2{1,m}{:,ss});
    elseif plotMean == 1
        allSpdsRwddCMat2 = meanRwddC2{1,m}{ss};
        allSpdsNRwddCMat2 = meanNRwddC2{1,m}{ss};
    end
    
    mRwddCToPlt2 = mean(allSpdsRwddCMat2);
    mNRwddCToPlt2 = mean(allSpdsNRwddCMat2);
    
    sampleNRwddC2 = length(allSpdsRwddCMat2(:,1));
    sampleNNRwddC2 = length(allSpdsNRwddCMat2(:,1));
    
    seUPRwddCToPlt2 = mRwddCToPlt2 + std(allSpdsRwddCMat2,0,1)/sqrt(sampleNRwddC2);
    seDNRwddCToPlt2 = mRwddCToPlt2 - std(allSpdsRwddCMat2,0,1)/sqrt(sampleNRwddC2);
    
    seUPNRwddCToPlt2 = mNRwddCToPlt2 + std( allSpdsNRwddCMat2,0,1)/sqrt(sampleNNRwddC2);
    seDNNRwddCToPlt2 = mNRwddCToPlt2 - std( allSpdsNRwddCMat2,0,1)/sqrt(sampleNNRwddC2);
    
    yMaxs(ss) = max([seUPRwddCToPlt2 seUPNRwddCToPlt2]);
    
    
    xPatch2 = [xLapsPos2 fliplr( xLapsPos2)];
    yPatchRwddC2 = [seUPRwddCToPlt2 fliplr(seDNRwddCToPlt2)];
    yPatchNRwddC2 = [seUPNRwddCToPlt2 fliplr(seDNNRwddCToPlt2)];
    
    prwdC2 = patch(xPatch2,yPatchRwddC2,hitsColor,'LineStyle','none');
    pnrwdC2 = patch(xPatch2,yPatchNRwddC2,faColor,'LineStyle','none');
    
    alpha(prwdC2,0.4)
    alpha(pnrwdC2,0.4)
end
text(65,5,['bM' num2str(69)],'FontSize',6,'FontWeight','bold')
yMax = round(max(yMaxs))+5;
ylim([0 yMax])
xlim([1 xLapsPos2(end)])

box on 
axD2 = gca;
axD2.XTick = [25+4 areaStart(1,m)+(areaEnd(1,m) - areaStart(1,m))/2];
axD2.XTickLabel = {'Sounds','Area'};
axD2.TickLength = [0 0];
axD2.YTick = [0 yMax/2  yMax];
axD2.FontSize = 6;
axD2.Units = 'centimeters';
axD2.Position = [panelD2Pos(1), panelD2Pos(2), panelD2Size(1), panelD2Size(2)];



% mice 6

m = 2;
plotMean = 1;
axes
hold on

for l = 1:length(ssToPlt(1,:))
    line([ssToPlt(1,l) ssToPlt(1,l)],[0 95],'Color',[soundColor 0.3],'LineStyle','--','LineWidth',3);
end

line([areaStart(1,m) areaStart(1,m)],[0 95],'Color',greyColor,'LineWidth',3)
line([areaEnd(1,m) areaEnd(1,m)],[0 95],'Color',greyColor,'LineWidth',3)


for ss = 1:length(cutRwddC2{1,m}(1,:))
    
    xLapsPos2 = allXLapsPos2{1,m}{1,ss};
    
    if plotMean == 0
        allSpdsRwddCMat2 = vertcat(cutRwddC2{1,m}{:,ss});
        allSpdsNRwddCMat2 = vertcat(cutNRwddC2{1,m}{:,ss});
    elseif plotMean == 1
        allSpdsRwddCMat2 = meanRwddC2{1,m}{ss};
        allSpdsNRwddCMat2 = meanNRwddC2{1,m}{ss};
    end
    
    mRwddCToPlt2 = mean(allSpdsRwddCMat2);
    mNRwddCToPlt2 = mean(allSpdsNRwddCMat2);
    
    sampleNRwddC2 = length(allSpdsRwddCMat2(:,1));
    sampleNNRwddC2 = length(allSpdsNRwddCMat2(:,1));
    
    seUPRwddCToPlt2 = mRwddCToPlt2 + std(allSpdsRwddCMat2,0,1)/sqrt(sampleNRwddC2);
    seDNRwddCToPlt2 = mRwddCToPlt2 - std(allSpdsRwddCMat2,0,1)/sqrt(sampleNRwddC2);
    
    seUPNRwddCToPlt2 = mNRwddCToPlt2 + std( allSpdsNRwddCMat2,0,1)/sqrt(sampleNNRwddC2);
    seDNNRwddCToPlt2 = mNRwddCToPlt2 - std( allSpdsNRwddCMat2,0,1)/sqrt(sampleNNRwddC2);
    
    yMaxs(ss) = max([seUPRwddCToPlt2 seUPNRwddCToPlt2]);
    
    
    xPatch2 = [xLapsPos2 fliplr( xLapsPos2)];
    yPatchRwddC2 = [seUPRwddCToPlt2 fliplr(seDNRwddCToPlt2)];
    yPatchNRwddC2 = [seUPNRwddCToPlt2 fliplr(seDNNRwddCToPlt2)];
    
    prwdC2 = patch(xPatch2,yPatchRwddC2,hitsColor,'LineStyle','none');
    pnrwdC2 = patch(xPatch2,yPatchNRwddC2,faColor,'LineStyle','none');
    
    alpha(prwdC2,0.4)
    alpha(pnrwdC2,0.4)
end
text(65,5,['bM' num2str(73)],'FontSize',6,'FontWeight','bold')
yMax = round(max(yMaxs))+5;
ylim([0 yMax])
xlim([1 xLapsPos2(end)])

box on 
axD3 = gca;
axD3.XTick = [25+4 areaStart(1,m)+(areaEnd(1,m) - areaStart(1,m))/2];
axD3.XTickLabel = {'Sounds','Area'};
axD3.TickLength = [0 0];
axD3.YTick = [0 yMax/2  yMax];
axD3.FontSize = 6;
axD3.Units = 'centimeters';
axD3.Position = [panelD3Pos(1), panelD3Pos(2), panelD3Size(1), panelD3Size(2)];



% mice 8

m = 3;
plotMean = 1;
axes
hold on

for l = 1:length(ssToPlt(1,:))
    line([ssToPlt(1,l) ssToPlt(1,l)],[0 95],'Color',[soundColor 0.3],'LineStyle','--','LineWidth',3);
end

line([areaStart(1,m) areaStart(1,m)],[0 95],'Color',greyColor,'LineWidth',3)
line([areaEnd(1,m) areaEnd(1,m)],[0 95],'Color',greyColor,'LineWidth',3)


for ss = 1:length(cutRwddC2{1,m}(1,:))
    
    xLapsPos2 = allXLapsPos2{1,m}{1,ss};
    
    if plotMean == 0
        allSpdsRwddCMat2 = vertcat(cutRwddC2{1,m}{:,ss});
        allSpdsNRwddCMat2 = vertcat(cutNRwddC2{1,m}{:,ss});
    elseif plotMean == 1
        allSpdsRwddCMat2 = meanRwddC2{1,m}{ss};
        allSpdsNRwddCMat2 = meanNRwddC2{1,m}{ss};
    end
    
    mRwddCToPlt2 = mean(allSpdsRwddCMat2);
    mNRwddCToPlt2 = mean(allSpdsNRwddCMat2);
    
    sampleNRwddC2 = length(allSpdsRwddCMat2(:,1));
    sampleNNRwddC2 = length(allSpdsNRwddCMat2(:,1));
    
    seUPRwddCToPlt2 = mRwddCToPlt2 + std(allSpdsRwddCMat2,0,1)/sqrt(sampleNRwddC2);
    seDNRwddCToPlt2 = mRwddCToPlt2 - std(allSpdsRwddCMat2,0,1)/sqrt(sampleNRwddC2);
    
    seUPNRwddCToPlt2 = mNRwddCToPlt2 + std( allSpdsNRwddCMat2,0,1)/sqrt(sampleNNRwddC2);
    seDNNRwddCToPlt2 = mNRwddCToPlt2 - std( allSpdsNRwddCMat2,0,1)/sqrt(sampleNNRwddC2);
    
    yMaxs(ss) = max([seUPRwddCToPlt2 seUPNRwddCToPlt2]);
    
    
    xPatch2 = [xLapsPos2 fliplr( xLapsPos2)];
    yPatchRwddC2 = [seUPRwddCToPlt2 fliplr(seDNRwddCToPlt2)];
    yPatchNRwddC2 = [seUPNRwddCToPlt2 fliplr(seDNNRwddCToPlt2)];
    
    prwdC2 = patch(xPatch2,yPatchRwddC2,hitsColor,'LineStyle','none');
    pnrwdC2 = patch(xPatch2,yPatchNRwddC2,faColor,'LineStyle','none');
    
    alpha(prwdC2,0.4)
    alpha(pnrwdC2,0.4)
end
text(65,5,['bM' num2str(76)],'FontSize',6,'FontWeight','bold')
yMax = round(max(yMaxs))+5;
ylim([0 yMax])
xlim([1 xLapsPos2(end)])

box on 
axD4 = gca;
axD4.XTick = [25+4 areaStart(1,m)+(areaEnd(1,m) - areaStart(1,m))/2];
axD4.XTickLabel = {'Sounds','Area'};
axD4.TickLength = [0 0];
axD4.YTick = [0 yMax/2  yMax];
axD4.FontSize = 6;
axD4.Units = 'centimeters';
axD4.Position = [panelD4Pos(1), panelD4Pos(2), panelD4Size(1), panelD4Size(2)];

lgD = legend([rwdPM nRwdPM],{'Stop Trials','No Stop Trials'},'Location','southoutside');
lgD.Units = 'centimeters';
lgD.Position = [panelD3Pos(1)+1.6, panelD4Pos(2) - 0.8, 0, 0.5];
lgD.Box = 'off';
lgD.Orientation = 'horizontal';
lgD.ItemTokenSize = [10,10];

%% Panel E, Behavioral Controls 

%%% catch trials


sSLocs = {'70','80','90'};
jitters = linspace(-0.2,0.2,8);


axes
for d = 1:3
    for m = 1:5
        
        % Plot Hit Rate
        mouseHitsMean = nanmean(catchDelHits{1,m}(:,d));
        %miceDelHits(:,d) = catchDelHits{1,m}(:,d);
        mouseHitsError = nanstd(catchDelHits{1,m}(:,d)) / sqrt(length(catchDelHits{1,m}(:,d)));
        hitsE = errorbar(d + jitters(m),mouseHitsMean,mouseHitsError);
        hitsE.Marker = '.';
        hitsE.MarkerSize = 9;
        hitsE.MarkerFaceColor = hitsColor;
        hitsE.MarkerEdgeColor = hitsColor;
        hitsE.Color = 'k';
        hold on
         
        
        % Plot FA Rate
        mouseFAMean = nanmean(catchDelFalseAlarms{1,m}(:,d));
        %miceDelFA(:,d) = catchDelFalseAlarms{1,m}(:,d);
        mouseFAError = nanstd(catchDelFalseAlarms{1,m}(:,d)) / sqrt(length(catchDelFalseAlarms{1,m}(:,d)));
        faE = errorbar(d + jitters(m),mouseFAMean,mouseFAError);
        faE.Marker = '.';
        faE.MarkerSize = 9;
        faE.MarkerFaceColor = faColor;
        faE.MarkerEdgeColor = faColor;
        faE.Color = 'k';
        
        hold on
    end
%     miceHitsMean = mean(mean(miceDelHits));
%     line([d+jitters(1)-0.05 d+jitters(end)+0.05],[miceHitsMean miceHitsMean],'Color',hitsColor)
%     
%     miceFAMean = mean(mean(miceDelFA));
%     line([d+jitters(1) d+jitters(end)],[miceFAMean miceFAMean],'Color',faColor)
    ylim([-0.1 1.1])
end
xlim([0.5 3.5])
line([0.5 3.5],[0.5 0.5],'Color',greyColor,'LineStyle','--','LineWidth',3)


axE1 = gca;
axE1.XTick = [1 2 3];
axE1.XTickLabel = sSLocs;
axE1.TickLength = [0 0];
axE1.YTick = [0  0.5 1];
axE1.FontSize = 6;
axE1.Units = 'centimeters';
axE1.Position = [panelE1Pos(1), panelE1Pos(2), panelE1Size(1), panelE1Size(2)]; 
xlabel('Sound Start (cm)','Fontsize',6)
ylabel('Proportion of Trials','Fontsize',6)
lgE1 = legend({'Hit Rate','FA Rate'},'Location','best');
lgE1.Units = 'centimeters';
lgE1.Position = [panelE1Pos(1) + 1.2, panelE1Pos(2) - 1.1, 0, 0.5];
lgE1.Orientation = 'horizontal';
lgE1.ItemTokenSize = [10,10];
legend('boxoff')

lE1 = text(lEPos(1),lEPos(2),'e','Units','centimeters');
lE1.FontSize = 9;
lE1.FontWeight = 'bold';


%%% Short Trials

jitters = linspace(-0.2,0.2,6);

axes
hold on

% Plot Hit Rate
for s = 1:length(shortHits(1,:))
    for m = 1:length(shortHits(:,1))
        if s == 3
            plot(s + jitters(m),shortHits(m,s),'Marker','v','MarkerSize',3,'MarkerFaceColor',hitsColor,'MarkerEdgeColor','none')
        else
            plot(s + jitters(m),shortHits(m,s),'Marker','.','MarkerSize',9,'Color',hitsColor)
        end
    end
end




% Plot FA Rate
for s = 1:length(shortFalseAlarms(1,:))
    for m = 1:length(shortFalseAlarms(:,1))
        if s == 3
            plot(s + jitters(m),shortFalseAlarms(m,s),'Marker','v','MarkerSize',3,'MarkerFaceColor',faColor,'MarkerEdgeColor','none')
        else
            plot(s + jitters(m),shortFalseAlarms(m,s),'Marker','.','MarkerSize',9,'Color',faColor)
        end
    end
end

pE2 = patch([2.5 3.5 3.5 2.5],[-0.05 -0.05 1.05 1.05],greyColor);
pE2.LineStyle = ':';
pE2.FaceColor = 'none';
pE2.EdgeAlpha = 0.4;
pE2.LineWidth = 2;


ylim([-0.1 1.1])
xlim([0.5 5.5])
line([0.5 5.5],[0.5 0.5],'Color',greyColor,'LineStyle','--','LineWidth',3)

box on
axE2 = gca;
axE2.XTick = [1 2 3 4 5];
axE2.TickLength = [0 0];
axE2.YTick = [0  0.5 1];
axE2.FontSize = 6;
axE2.Units = 'centimeters';
axE2.Position = [panelE2Pos(1), panelE2Pos(2), panelE2Size(1), panelE2Size(2)]; 
xlabel('Sessions','Fontsize',6)
ylabel('Proportion of Trials','Fontsize',6)


%%% No Area Trials

jitters = linspace(-0.2,0.2,6);


axes
hold on

% Plot Hit Rate
for s = 1:length(noAreaHits(1,:))
    for m = 1:length(noAreaHits(:,1))
        if s == 3
            scatter(s + jitters(m),noAreaHits(m,s),10,'v','MarkerFaceColor',hitsColor,'MarkerEdgeColor','none','MarkerFaceAlpha',0.3);
        else
            plot(s + jitters(m),noAreaHits(m,s),'Marker','.','MarkerSize',9,'Color',hitsColor)
        end
    end
end



% Plot FA Rate
for s = 1:length(noAreaFalseAlarms(1,:))
    for m = 1:length(noAreaFalseAlarms(:,1))
        if s == 3
            scatter(s + jitters(m),noAreaFalseAlarms(m,s),10,'v','MarkerFaceColor',faColor,'MarkerEdgeColor','none','MarkerFaceAlpha',0.3);
        else
            plot(s + jitters(m),noAreaFalseAlarms(m,s),'Marker','.','MarkerSize',9,'Color',faColor)
        end
    end
end

pE3 = patch([2.5 3.5 3.5 2.5],[-0.05 -0.05 1.05 1.05],greyColor);
pE3.LineStyle = ':';
pE3.FaceColor = 'none';
pE3.EdgeAlpha = 0.4;
pE3.LineWidth = 2;



ylim([-0.1 1.1])
xlim([0.5 5.5])
line([0.5 5.5],[0.5 0.5],'Color',greyColor,'LineStyle','--','LineWidth',3)

box on
axE3 = gca;
axE3.XTick = [1 2 3 4 5];
axE3.TickLength = [0 0];
axE3.YTick = [0  0.5 1];
axE3.FontSize = 6;
axE3.Units = 'centimeters';
axE3.Position = [panelE3Pos(1), panelE3Pos(2), panelE3Size(1), panelE3Size(2)]; 
xlabel('Sessions','Fontsize',6)
ylabel('Proportion of Trials','Fontsize',6)









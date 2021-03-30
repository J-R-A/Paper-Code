clear all
close all

%% Load and process needed Data

load('/Volumes/Cerebro/Recording Data/SessionsMap.mat')

spdColor = [0, 119, 187] / 255;
hitsColor = [0, 153, 136] / 255;
faColor = [204, 51, 17] / 255;
greyColor = [187, 187, 187] / 255;
lightGreyColor = [220, 220, 220] / 255;
soundColor = [238,51,119] / 255;

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
        clearvars -except m ss  sessionsMap s coefsS coefsO featsO featsS nShank
    end
end


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
    
%     for n = 1:length(frS1{s}(1,:))
%         
%         % Stop trial speeds
%         nFrS1 = frS1{s}(:,n);
%         nCoefS1 = cS1{s}(n);
%         fImpS1 = sum(abs(nFrS1 * nCoefS1));
%         
%         
%         % No Stop trial speeds
%         nFrS2 = frS2{s}(:,n);
%         nCoefS2 = cS2{s}(n);
%         fImpS2 = sum(abs(nFrS2 * nCoefS2));
%         
%         
%         % Trial outcomes
%         nFrO = frO{s}(:,n);
%         nCoefO = cO{s}(n);
%         fImpO = sum(abs(nFrO * nCoefO));
%         
%         
%         featImp{s}(n,1) =  fImpS1 / allImpS1;
%         featImp{s}(n,2) =  fImpS2 / allImpS2;
%         featImp{s}(n,3) =  fImpO / allImpO;
%         
%         
%     end
    
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


% plot cumsum of feature importance of selected session

s = 13;



plot(featImpS1{s},'Color',hitsColor,'LineWidth',2)
hold on
plot(featImpS2{s},'Color',faColor,'LineWidth',2)
hold on
plot(featImpO{s},'Color',spdColor,'LineWidth',2)

line([0 length(featImpS1{s})],[0.8 0.8],'Color','r','LineStyle','--')
ylim([0 1.05])
xlim([0 length(featImpS1{s})])
xlabel('Sorted Neurons')
ylabel('Fraction Variance Explained')


% plot mean +- SEM #neurs untill treshold
figure
bar(1,mean(treshNum(:,1)),'FaceColor',hitsColor);
hold on
bar(2,mean(treshNum(:,2)),'FaceColor',faColor);
hold on
bar(3,mean(treshNum(:,3)),'FaceColor',spdColor);

hold on 
errorbar(mean(treshNum),std(treshNum)/sqrt(length(treshNum(:,1))),'.','Color','k')

ylabel('# neurons until Treshold')
ax = gca;
ax.XTick = [1 2 3];
ax.XTickLabel = {'Speed Stop','Speed No Stop','Trial Id'};

% Scatters and correlations of coefs of different models for example session 

s = 2;

subplot(2,2,1) 
scatter(cS1{s},cS2{s},[],hitsColor,'filled')
xlim([-0.5 0.5])
ylim([-0.5 0.5])
text(0.15,-0.4,['CorrCoef: ',num2str(round(corrs(s,1),2))])
title('Stop Speed vs No Stop Speed')
xlabel('coefs')
ylabel('coefs')
subplot(2,2,2) 
scatter(cS1{s},cO{s},[],faColor,'filled')
xlim([-0.5 0.5])
ylim([-0.5 0.5])
text(0.15,-0.4,['CorrCoef: ',num2str(round(corrs(s,2),2))])
title('Stop Speed vs Trial Id')
xlabel('coefs')
ylabel('coefs')
subplot(2,2,4)
scatter(cS2{s},cO{s},[],spdColor,'filled')
xlim([-0.5 0.5])
ylim([-0.5 0.5])
text(0.15,-0.4,['CorrCoef: ',num2str(round(corrs(s,3),2))])
title('No Stop Speed vs Trial Id')
xlabel('coefs')
ylabel('coefs')


subplot(2,2,3)
scatter(ones(1,length(corrs)),corrs(:,1),[],hitsColor,'filled','MarkerFaceAlpha',0.2)
hold on
plot(1,mean(corrs(:,1)),'o','MarkerFaceColor',hitsColor,'MarkerEdgeColor',hitsColor)
hold on
errorbar(1,mean(corrs(:,1)),std(corrs(:,1))/sqrt(length(corrs(:,1))),'.','Color','k')
hold on
scatter(ones(1,length(corrs))*2,corrs(:,2),[],faColor,'filled','MarkerFaceAlpha',0.2)
hold on
plot(2,mean(corrs(:,2)),'o','MarkerFaceColor',faColor,'MarkerEdgeColor',faColor)
hold on
errorbar(2,mean(corrs(:,2)),std(corrs(:,2))/sqrt(length(corrs(:,2))),'.','Color','k')
hold on
scatter(ones(1,length(corrs))*3,corrs(:,3),[],spdColor,'filled','MarkerFaceAlpha',0.2)
hold on
plot(3,mean(corrs(:,3)),'o','MarkerFaceColor',spdColor,'MarkerEdgeColor',spdColor)
hold on
errorbar(3,mean(corrs(:,3)),std(corrs(:,3))/sqrt(length(corrs(:,3))),'.','Color','k')
hold on
xlim([0.5 3.5])
ylim([-0.6 0.8])
ylabel('Corr Coef')

ax = gca;
ax.XTick = [1 2 3];
ax.XTickLabel = {'S1vsS2', 'S1vsTId','S2vsTId'};

% Exploration of location per shank of neurons coefs 

% Calculate shank importance

shanks = 1:6;

for s = 1:16
    
%     % Stop trial speeds
%     allFrS1 = frS1{s}(:,nShank{s} ~= 7);
%     allCoefS1 = cS1{s}(nShank{s} ~= 7);
%     allImpS1 = sum(allFrS1 * abs(allCoefS1));
%     
%     % No Stop trial speeds
%     allFrS2 = frS2{s}(:,nShank{s} ~= 7);
%     allCoefS2 = cS2{s}(nShank{s} ~= 7);
%     allImpS2 = sum(allFrS2 * abs(allCoefS2));
%     
%     % Trial outcomes
%     allFrO = frO{s}(:,nShank{s} ~= 7);
%     allCoefO = cO{s}(nShank{s} ~= 7);
%     allImpO = sum(allFrO * abs(allCoefO));
    
    for sk = shanks
        
         nSkImps{1}{s,sk} = featImp{s}(nShank{s} == sk,1);
         nSkImps{2}{s,sk} = featImp{s}(nShank{s} == sk,2);
         nSkImps{3}{s,sk} = featImp{s}(nShank{s} == sk,3);
        
         skImp{1}(s,sk) =  mean(nSkImps{1}{s,sk});
         skImp{2}(s,sk) =  mean(nSkImps{2}{s,sk});
         skImp{3}(s,sk) =  mean(nSkImps{3}{s,sk});
        
        
        
%          skFrS1 = frS1{s}(:,nShank{s} == sk);
%          skCoefS1 = cS1{s}(nShank{s} == sk);
%          skImpS1 = sum(skFrS1 * abs(skCoefS1));
%          
%          skFrS2 = frS2{s}(:,nShank{s} == sk);
%          skCoefS2 = cS2{s}(nShank{s} == sk);
%          skImpS2 = sum(skFrS2 * abs(skCoefS2));
%          
%          skFrO = frO{s}(:,nShank{s} == sk);
%          skCoefO = cO{s}(nShank{s} == sk);
%          skImpO = sum(skFrO * abs(skCoefO));
%          
%          
%         skImp{1}(s,sk) =  skImpS1 / allImpS1;
%         skImp{2}(s,sk) =  skImpS2 / allImpS2;
%         skImp{3}(s,sk) =  skImpO / allImpO;
    end
    
end

% Plot example session 

s = 1;

filt = nShank{s} ~= 7;
subplot(2,3,1)
scatter(nShank{s}(filt),featImp{s}(filt,1),[],hitsColor,'filled','MarkerFaceAlpha',0.15) 
hold on
for i = 1:6
errorbar(i,mean(featImp{s}(nShank{s}==i,1)),std(featImp{s}(nShank{s}==i,1))/sqrt(sum(nShank{s}==i)),'o','MarkerFaceColor',hitsColor,'MarkerEdgeColor',hitsColor,'Color','k')
hold on
end
ylim([-0.05 inf])
xlim([0.5 6.5])
xlabel('Shanks')
ylabel('Fraction Var Exp')
title('Stop Speed')


subplot(2,3,2)
scatter(nShank{s}(filt),featImp{s}(filt,2),[],faColor,'filled','MarkerFaceAlpha',0.15) 
hold on
for i = 1:6
errorbar(i,mean(featImp{s}(nShank{s}==i,2)),std(featImp{s}(nShank{s}==i,2))/sqrt(sum(nShank{s}==i)),'o','MarkerFaceColor',faColor,'MarkerEdgeColor',faColor,'Color','k')
hold on
end
ylim([-0.05 inf])
xlim([0.5 6.5])
xlabel('Shanks')
ylabel('Fraction Var Exp')
title('No Stop Speed')

subplot(2,3,3)
scatter(nShank{s}(filt),featImp{s}(filt,3),[],spdColor,'filled','MarkerFaceAlpha',0.15) 
hold on
for i = 1:6
errorbar(i,mean(featImp{s}(nShank{s}==i,3)),std(featImp{s}(nShank{s}==i,3))/sqrt(sum(nShank{s}==i)),'o','MarkerFaceColor',spdColor,'MarkerEdgeColor',spdColor,'Color','k')
hold on
end
ylim([-0.05 inf])
xlim([0.5 6.5])
xlabel('Shanks')
ylabel('Fraction Var Exp')
title('Trial Id')


jitters = linspace(-1.8,1.8,16);
xShanks = linspace(1,30,6);

subplot(2,3,4)
for i = 1:6
    for s = 1:16
        errorbar(xShanks(i)+jitters(s),mean(featImp{s}(nShank{s}==i,1)),std(featImp{s}(nShank{s}==i,1))/sqrt(sum(nShank{s}==i)),'o','MarkerFaceColor',hitsColor,'MarkerEdgeColor',hitsColor,'Color','k')
        hold on
    end 
end
xlim([xShanks(1)-5 xShanks(end)+5])
ylim([-0.005 inf])
xlabel('Shanks')
ylabel('Fraction Var Exp')
ax = gca;
ax.XTick = xShanks;
ax.XTickLabel = {'1','2','3','4','5','6'};

subplot(2,3,5)
for i = 1:6
    for s = 1:16
        errorbar(xShanks(i)+jitters(s),mean(featImp{s}(nShank{s}==i,2)),std(featImp{s}(nShank{s}==i,2))/sqrt(sum(nShank{s}==i)),'o','MarkerFaceColor',faColor,'MarkerEdgeColor',faColor,'Color','k')
        hold on
    end 
end
xlim([xShanks(1)-5 xShanks(end)+5])
ylim([-0.005 inf])
xlabel('Shanks')
ylabel('Fraction Var Exp')
ax = gca;
ax.XTick = xShanks;
ax.XTickLabel = {'1','2','3','4','5','6'};

subplot(2,3,6)
for i = 1:6
    for s = 1:16
        errorbar(xShanks(i)+jitters(s),mean(featImp{s}(nShank{s}==i,3)),std(featImp{s}(nShank{s}==i,3))/sqrt(sum(nShank{s}==i)),'o','MarkerFaceColor',spdColor,'MarkerEdgeColor',spdColor,'Color','k')
        hold on
    end 
end
xlim([xShanks(1)-5 xShanks(end)+5])
ylim([-0.005 inf])
xlabel('Shanks')
ylabel('Fraction Var Exp')

ax = gca;
ax.XTick = xShanks;
ax.XTickLabel = {'1','2','3','4','5','6'};



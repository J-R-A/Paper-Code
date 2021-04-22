clear all



sessions = [1 2 13 14]; % Sessions to load data from
stdv = 0.14; % Stdv of kernel used to smooth data
def = 0.5; % Parameter of allowed deformation
tSelectMethod = 1; % Method used to select similar trials: 0 - manual; 1 - based on correlations; 2 - based on vector similarity
intDels = [1 2 3 4]; % Delays on which the analysis was run
delGroup = [1];
intCond = 1; % Trial's condition: 1-hits; 2-misses; 3-correct rejections; 4-false alarms
p1 = 0.25; %Cut trials p1 seconds after the higher area time across trials
p2 = 0.5;  %Cut trials p2 seconds before the lowest moving after stopping time across trials



% Load data

for s = 1:length(sessions)
    for d = 1:length(intDels)
        
        disp(['loading delay ', num2str(intDels(d)),' from session ', num2str(sessions(s))])
        [Name] = SessionNumber2Name(sessions(s));
        load(['/Volumes/Cerebro/Recording Data/',Name,'/Predict from Neurons/SpdPredLagsBootCheck Trials tsMethod_',num2str(tSelectMethod),' Cond_',num2str(intCond),' Delay_',num2str(delGroup),'_',num2str(intDels(d)),' P1_',num2str(p1),' P2_',num2str(p2),' Sigma_',num2str(stdv),' Deform_',num2str(def),'.mat'],'testData','origData','valData');
        
        allOrig{s,d} = origData;
        allTest{s,d} = testData;
        allVal{s,d} = valData;
      
    end
end


% Calculate R2s


for s = 1:length(sessions)
    for iD = 1:length(intDels)
        allSpds = allOrig{s,iD}.Speeds;
        allPredsLag = allTest{s,iD}.Preds{1,2};
        for rs = 1:length(allPredsLag(1,:))
            teTrials = allTest{s,iD}.testTrials{1,2}{1,rs};
            teSpeeds = cell2mat(allSpds(teTrials));
            tePredsLag = allPredsLag{1,rs};
            
            RSS = sum((teSpeeds - tePredsLag).^2);
            TSS = sum((teSpeeds - mean(teSpeeds)).^2);
            R2{s,iD}(rs) = 1 - (RSS / TSS);
        end
    end
end



%%% Plot no lag performance results for all sessions and delay groups

% Plot R2

s = 3;
iD = 3;



for s = 1:length(sessions)
for g = 1:length(intDels)

end
end





subplot(2,2,1)
plot(R2{s,iD})
hold on
line([0 length(R2{s,iD})],[0 0],'Color','r')
line([0 length(R2{s,iD})],[1 1],'Color','r')
ylim([-inf 1.5])
xlabel('resamples')
ylabel('R2')

subplot(2,2,2)
plot(R2{s,iD},'.')
ylim([0 1])
xlabel('resamples')
ylabel('R2')

subplot(2,2,3)
gProp = sum(R2{s,iD} >= 0.4);
bProp = sum(R2{s,iD} < 0.4);
pie([gProp bProp])
legend('R2>0.4','R2<0.4','Location','southoutside')



% Get speeds and preds from individual trials used in the test groups of all resamples

allSpds = allOrig{s,iD}.Speeds;
allPredsNoLag = allTest{s,iD}.Preds{1,2};
allTeTrials = allTest{s,iD}.testTrials{1,2};
allTP = cell(1,11);
a = 1;
for rs = 1:length(allPredsLag)
    teTrials = allTeTrials{rs};
    for t = 1:2
        spds = allSpds{teTrials(t),1};
        tLen = length(spds);
        
        % separate the 2 trials in the test group
        if t == 1
        preds = allPredsNoLag{rs}(1:tLen);
        else
        preds = allPredsNoLag{rs}(tLen+1:end);    
        end
        
 
        if isempty(allTP{teTrials(t)})
            allTP{teTrials(t)}(1,:) = preds;
        else
            allTP{teTrials(t)} =[allTP{teTrials(t)}; preds'];
        end
       
        
    end
end


% Plot speeds and preds from individual trials used in the test groups of all resamples

figure
for  i = 1:11
subplot(3,4,i)
plot(allTP{i}')
title(['Trial ',num2str(i)])
end

figure
for  i = 1:11
subplot(3,4,i)
plot(allSpds{i}')
title(['Trial ',num2str(i)])
end

% plot proportion of resamples with R2 above and below a given treshold

sb = 1;
perfTresh = 0.4;
for s = 1:length(sessions)
    for iD = 1:length(intDels)
        gProp = sum(R2{s,iD} >= perfTresh);
        bProp = sum(R2{s,iD} < perfTresh);
        subplot(length(sessions),length(intDels),sb)
        pie([gProp bProp])
        
        if s == 1
            title(['Del',num2str(intDels(iD))])
        end
       
        sb = sb+1;
    end
end




%%%%% Several plots not yet commented 


sb = 1;
for s = 1:3
    for iD = 1:5
        gProp = sum(R2{s,iD} >= 0.4);
        bProp = sum(R2{s,iD} < 0.4);
        subplot(3,5,sb)
        pie([gProp bProp])
        sb = sb+1;
    end
end


for s = 1:length(sessions)
    for iD = 1:length(intDels)
        b = 1;
        g = 1;
        m = 1;
        vb = 1;
        for rs = 1:length(allPredsLag(1,:))
            if R2{s,iD}(rs) >= 0.4
            gLags{s,iD}(g,:) = allTest{s,iD}.finalLags{1,rs}*-1;
            g = g+1;
            elseif R2{s,iD}(rs) >= 0 && R2{s,iD}(rs) < 0.4
            mLags{s,iD}(b,:) = allTest{s,iD}.finalLags{1,rs}*-1;
            m = m+1;
            elseif R2{s,iD}(rs) >= -10 && R2{s,iD}(rs) < 0
            bLags{s,iD}(b,:) = allTest{s,iD}.finalLags{1,rs}*-1;
            b = b+1;
            elseif R2{s,iD}(rs) < -10 
            vbLags{s,iD}(b,:) = allTest{s,iD}.finalLags{1,rs}*-1;
            vb = b+1;
            end
        end 
    end
end


s = 3;
figure
sb = 1;
for iD = 1:length(intDels)
    subplot(5,4,sb)
    plot(mean(gLags{s,iD}))
    xl
    
    subplot(5,4,sb+1)
    plot(mean(mLags{s,iD}))
    
    subplot(5,4,sb+2)
    plot(mean(bLags{s,iD}))
    
    subplot(5,4,sb+3)
    plot(mean(vbLags{s,iD}))
    
    
    sb = sb + 4;
end






iD = 5;
s = 3;


lagErr = allTest{s,iD}.RMSETest{1,2};
nLagErr = allTest{s,iD}.RMSETest{1,1};

diffErr = nLagErr - lagErr;

for iD = 1:5
   lagErr = allTest{s,iD}.RMSETest{1,2};
   nLagErr = allTest{s,iD}.RMSETest{1,1};
   subplot(2,3,iD)
  
   plot(nLagErr)
   hold on
   plot(lagErr)
  
end

figure
for iD = 1:5
   lagErr = allTest{s,iD}.RMSETest{1,2};
   nLagErr = allTest{s,iD}.RMSETest{1,1};
   diffErr = nLagErr - lagErr;
   subplot(2,3,iD)
   plot(diffErr)
end



s = 4;
for iD = 1:5
   lagErr = allTest{s,iD}.RMSETest{1,2};
   
   for rs = 1:length(allTest{s,iD}.testTrials{1,1}(1,:))
       trTrials = allTest{s,iD}.testTrials{1,1}{1,rs};
       uTrials(rs) = length(unique(trTrials));
   end
   
   subplot(2,3,iD)
   scatter(uTrials,lagErr)
end

s = 3;

for iD = 1:5
    sumPreds = zeros(11,11);
    trialPreds = zeros(11,11);
    
    
    
    for rs = 1:length(allTest{s,iD}.testTrials{1,1}(1,:))
        teTrials = allTest{s,iD}.testTrials{1,2}{1,rs};
        lagErr = allTest{s,iD}.RMSETest{1,2}(rs);
        trialPreds(teTrials(1),teTrials(2)) = trialPreds(teTrials(1),teTrials(2)) + 1;
        sumPreds(teTrials(1),teTrials(2)) = meanPreds(teTrials(1),teTrials(2)) + lagErr;
    end
    
    sumPreds(trialPreds == 0) = nan;
    trialPreds(trialPreds == 0) = nan;
    meanPreds = sumPreds ./ trialPreds;
    
    subplot(2,3,iD)
    imagesc(meanPreds)
    colorbar;
end



% Inspect individual resamples

s = 3;
iD = 5;
rs = 92;

allSpds = allOrig{s,iD}.Speeds;
allPredsLag = allTest{s,iD}.Preds{1,2};
allPredsNoLag = allTest{s,iD}.Preds{1,1};
trTrials = allTest{s,iD}.testTrials{1,1}{1,rs};
teTrials = allTest{s,iD}.testTrials{1,2}{1,rs};
trSpeeds = cell2mat(allSpds(trTrials));
teSpeeds = cell2mat(allSpds(teTrials));
tePredsLag = allPredsLag{1,rs};
tePredsNoLag = allPredsNoLag{rs};
lagErr = allTest{s,iD}.RMSETest{1,2}(rs);
nLagErr = allTest{s,iD}.RMSETest{1,1}(rs);


figure

plot(teSpeeds(180:end))
hold on
plot(tePredsNoLag(180:end),'k')
hold on
plot(tePredsLag(180:end),'r')
title(['No:',num2str(nLagErr),';yes:',num2str(lagErr)])

figure
plot(allTest{s,iD}.finalLags{1,rs}*-1)


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
        
        s = s+1;
        clearvars -except m ss  sessionsMap s coefsS coefsO featsO featsS
    end
end


d=12;

for s = 1:length(coefsS(:,1))
    
    % get coefs
    cS1{s} = coefsS{s,d}{1,1}(2:end);
    cS2{s} = coefsS{s,d}{1,2}(2:end);
    cO{s} = coefsO{s,d}{1,1}(2:end);
    
    % get frs
    frS1{s} = featsS{s,d}{1,1};
    frS2{s} = featsS{s,d}{1,2}; 
    frO{s} = featsO{s,d}{1,1};
    
    % calculate correlations between coefs vectors
    c1 = corrcoef(cS1{s},cS2{s});
    c2 = corrcoef(cS1{s},cO{s});
    c3 = corrcoef(cS2{s},cO{s});
    
    cS1S2(s) = c1(1,2);
    cS1O(s) = c2(1,2);
    cS2O(s) = c3(1,2);
    
    % calculate feature importance
    
    
    % Stop trial speeds
    allFrS1 = frS1{s};
    allCoefS1 = cS1{s};
    allImpS1 = sum(allFrS1 * abs(allCoefS1));
    
    % No Stop trial speeds
    allFrS2 = frS2{s};
    allCoefS2 = cS2{s};
    allImpS2 = sum(allFrS2 * abs(allCoefS2));
    
    % Trial outcomes
    allFrO = frO{s};
    allCoefO = cO{s};
    allImpO = sum(allFrO * abs(allCoefO));
    
    for n = 1:length(frS1{s}(1,:))
        
        % Stop trial speeds
        nFrS1 = frS1{s}(:,n);
        nCoefS1 = cS1{s}(n);
        fImpS1 = sum(abs(nFrS1 * nCoefS1));
        
        
        % No Stop trial speeds
        nFrS2 = frS2{s}(:,n);
        nCoefS2 = cS2{s}(n);
        fImpS2 = sum(abs(nFrS2 * nCoefS2));
        
        
        % Trial outcomes
        nFrO = frO{s}(:,n);
        nCoefO = cO{s}(n);
        fImpO = sum(abs(nFrO * nCoefO));
        
        
        featImp{s}(n,1) =  fImpS1 / allImpS1;
        featImp{s}(n,2) =  fImpS2 / allImpS2;
        featImp{s}(n,3) =  fImpO / allImpO;
        
        
    end
    
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



plot(featImpS1{s})
hold on
plot(featImpS2{s})
hold on
plot(featImpO{s})

line([0 length(featImpS1{s})],[0.8 0.8],'Color','r','LineStyle','--')


% plot mean +- SEM #neurs untill treshold

bar(mean(treshNum))
hold on 
errorbar(mean(treshNum),std(treshNum)/sqrt(length(treshNum(:,1))))







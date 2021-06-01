clear all

sessions = [1 13]; % Sessions to load data from
stdvs = [0.14 0.2 0.24 0.3]; % Stdv of kernel used to smooth data
deformations = [0.5 0.75 1]; % Parameter of allowed deformation
tSelectMethod = 1; % Method used to select similar trials: 0 - manual; 1 - based on correlations; 2 - based on vector similarity
delGroup = 1;
intDels = [1 2 3 4]; % Delays on which the analysis was run
intCond = 1; % Trial's condition: 1-hits; 2-misses; 3-correct rejections; 4-false alarms
p1 = 0.25; %Cut trials p1 seconds after the higher area time across trials
p2 = 0.5;  %Cut trials p2 seconds before the lowest moving after stopping time across trials

testFields = {'RMSE','totError','totErrorTest','predSqrdError','predSqrdErrorTest','lagsComp',...
                    'betasComp','betasCompTest','betas','betasDiff','regParams'};

%********************  Load and process data ******************

for s = 1:length(sessions)
    for st = 1:length(stdvs)
        for d = 1:length(intDels)
            for def = 1:length(deformations)
                
                disp(['loading deformation ',num2str(deformations(def)),' from delay ', num2str(intDels(d)),' of kernel width ',num2str(stdvs(st)),' session ', num2str(sessions(s))])
                [Name] = SessionNumber2Name(sessions(s));
                load(['/Volumes/Cerebro/Recording Data/',Name,'/Predict from Neurons/SpdPredLagsBootCheck Trials tsMethod_',num2str(tSelectMethod),' Cond_',num2str(intCond),' Delay_',num2str(delGroup),'_',num2str(intDels(1,d)),' P1_',num2str(p1),' P2_',num2str(p2),' Sigma_',num2str(stdvs(st)),' Deform_',num2str(deformations(def)),'.mat'],'testData','origData');
                
                allOrig{s,st}{d,def} = origData;
                
                testData = rmfield(testData,testFields);
                allTest{s,st}{d,def} = testData;
                
            end
        end
    end
end


% Calculate test R2


for s = 1:length(allTest(:,1))
    for st = 1:length(allTest(1,:))
        for d = 1:length(allTest{s,st}(:,1))
            for def = 1:length(allTest{s,st}(1,:))
                
                allSpds = allOrig{s,st}{d,def}.Speeds;
                tePredsLag{s,st}{d,def} = allTest{s,st}{d,def}.Preds{1,2};
                tePredsNoLag{s,st}{d,def} = allTest{s,st}{d,def}.Preds{1,1};
                for rs = 1:length(tePredsLag{s,st}{d,def}(1,:))
                    teTrials = allTest{s,st}{d,def}.testTrials{1,2}{1,rs};
                    teSpeeds{s,st}{d,def}{rs} = cell2mat(allSpds(teTrials));
                    
                    tePNoLag = tePredsNoLag{s,st}{d,def}{1,rs};
                    RSSNoLag = sum((teSpeeds{s,st}{d,def}{rs} - tePNoLag).^2);
                    TSSNoLag = sum((teSpeeds{s,st}{d,def}{rs} - mean(teSpeeds{s,st}{d,def}{rs})).^2);
                    teR2{s,st}{d,def}(rs,1) = 1 - (RSSNoLag / TSSNoLag);
                    
                    tePLag = tePredsLag{s,st}{d,def}{1,rs};
                    RSSLag = sum((teSpeeds{s,st}{d,def}{rs} - tePLag).^2);
                    TSSLag = sum((teSpeeds{s,st}{d,def}{rs} - mean(teSpeeds{s,st}{d,def}{rs})).^2);
                    teR2{s,st}{d,def}(rs,2) = 1 - (RSSLag / TSSLag);
                end
                
            end
        end
    end
end




% Get time of events


for s = 1:length(allOrig(:,1))
    for st = 1:length(allOrig(1,:))
        for d = 1:length(allOrig{s,st}(:,1))
            for def = 1:length(allOrig{s,st}(1,:))
                
                events = allOrig{s,st}{d,def}.eventTrials;
                spds = allOrig{s,st}{d,def}.Speeds;
                tms = allOrig{s,st}{d,def}.Times;
                dsts = allOrig{s,st}{d,def}.Dists;
                for  t = 1:length(events(:,1))
                    tTms =  tms{t,1} -  tms{t,1}(1);
                    tDsts = dsts{t,1} - dsts{t,1}(1);
                    
                    soundTimes{s,st}{d,def}(t) = events{t,1}(1,2);
                    [~,sIdx] = min(abs(tTms - soundTimes{s,st}{d,def}(t)));
                    
                    areaTimes{s,st}{d,def}(t) = events{t,1}(1,4);
                    [~,aIdx] =  min(abs(tTms - areaTimes{s,st}{d,def}(t)));
                    
                    movIdx = find(spds{t,1} <=5,1,'last');
                    movTimes{s,st}{d,def}(t) =tTms(movIdx);
                    
                    [~,maxIdx] = max(spds{t,1}(sIdx+10:aIdx));
                    maxSTimes{s,st}{d,def}(t) =tTms(maxIdx+sIdx+10);
                end
                
            end
        end
    end
end



% Get train predictions and calculate train R2s

for s = 1:length(allOrig(:,1))
    for st = 1:length(allOrig(1,:))
        for d = 1:length(allOrig{s,st}(:,1))
            for def = 1:length(allOrig{s,st}(1,:))
                
                allSpds = allOrig{s,st}{d,def}.Speeds;
                trPredsLag{s,st}{d,def} = allTest{s,st}{d,def}.PredsTr(:,2)';
                trPredsNoLag{s,st}{d,def} = allTest{s,st}{d,def}.PredsTr(:,1)';
                for rs = 1:length(trPredsLag{s,st}{d,def}(1,:))
                    trTrials = allTest{s,st}{d,def}.testTrials{1,1}{1,rs};
                    trSpds{s,st}{d,def}{1,rs} = cell2mat(allSpds(trTrials));
                    
                    trPLag = trPredsLag{s,st}{d,def}{1,rs};
                    RSSLag = sum((trSpds{s,st}{d,def}{1,rs} - trPLag).^2);
                    trRMSE{s,st}{d,def}(rs,2) = sqrt(RSSLag / length(trPLag));
                    TSSLag = sum((trSpds{s,st}{d,def}{1,rs} - mean(trSpds{s,st}{d,def}{1,rs})).^2);
                    trR2{s,st}{d,def}(rs,2) = 1 - (RSSLag / TSSLag);
                    
                    trPNoLag = trPredsNoLag{s,st}{d,def}{1,rs};
                    RSSNoLag = sum((trSpds{s,st}{d,def}{1,rs} - trPNoLag).^2);
                    trRMSE{s,st}{d,def}(rs,1) = sqrt(RSSNoLag / length(trPNoLag));
                    TSSNoLag = sum((trSpds{s,st}{d,def}{1,rs} - mean(trSpds{s,st}{d,def}{1,rs})).^2);
                    trR2{s,st}{d,def}(rs,1) = 1 - (RSSNoLag / TSSNoLag);
                    
                end
            end
        end
    end
end











% Compare out of sample performance for different kernel widths across sessions and deformagtions 


jitters = [-0.2 0 0.2];
sb = 1;
clrs = {'r','b','g'};
figure
for s = 1:length(teR2(:,1))
    for d = 1:length(teR2{s,1}(:,1))
        subplot(length(teR2(:,1)),length(teR2{s,1}(:,1)),sb)
        for st = 1:length(teR2(1,:))
            for def = 1:length(teR2{s,st}(1,:))
                
                mR2 = mean(teR2{s,st}{d,def}(teR2{s,st}{d,def}(:,2)>0,2));
                eR2 = std(teR2{s,st}{d,def}(teR2{s,st}{d,def}(:,2)>0,2)) / sqrt(length(teR2{s,st}{d,def}(teR2{s,st}{d,def}(:,2)>0,2)));
                
                errorbar(st+jitters(def),mR2,eR2,clrs{def})
                hold on
                
            end
        end
        ax = gca;
        ax.XTick = 1:length(teR2(1,:));
        ax.XTickLabel = {'0.14','0,2','0.24','0.3'};
        
        ylim([0.65 0.85])
        xlim([0 4.8])
        
        if d == 1
           ylabel(['Session: ', num2str(sessions(s))]) 
        end
        
        if s == 1
           title(['Delay : ', num2str(d)]) 
        end
        
        sb = sb+1;
    end
end




% Plot lag gains in all sessions and delays across smoothings and deformations


jitters = [-0.2 0 0.2];
sb = 1;
clrs = {'r','b','g'};
figure
for s = 1:length(teR2(:,1))
    for d = 1:length(teR2{s,1}(:,1))
        subplot(length(teR2(:,2)),length(teR2{s,1}(:,1)),sb)
        
        for st = 1:length(teR2(1,:))
            for def = 1:length(teR2{s,st}(1, :))
                diffPerf =allTest{s,st}{d,def}.RMSETest{1,1} - allTest{s,st}{d,def}.RMSETest{1,2};
                diffImprove = diffPerf ./ allTest{s,st}{d,def}.RMSETest{1,1} * 100;
                
                mDI = mean(diffImprove);
                eDI = std(diffImprove) / sqrt(length(diffImprove));
                
                
                errorbar(st+jitters(def),mDI,eDI,clrs{def})
                hold on
                
            end
        end
        %line([1 length(diffImprove)],[0 0],'Color','r')
        sb = sb+1;
        xlim([0 4.8])
    end
end





% Plot mean lag functions, in all sessions and delays, fordifferent deformations


figure
sb = 1;
st = 1;


for st = 1:length(allTest(1,:))
    figure
    sb = 1;
    for s = 1:length(allTest(:,1))
        for d = 1:length(allTest{s,st}(:,1))
            subplot(length(allTest(:,1)),length(allTest{s,st}(:,1)),sb)
            hold on
            
            for def = 1:length(allTest{s,st}(1,:))
                
                sS = soundTimes{s,st}{d,def};
                mSS = mean(sS);
                stdSS = std(sS);
                
                aS = areaTimes{s,st}{d,def};
                mAS = mean(aS);
                stdAS = std(aS);
                
                mT = movTimes{s,st}{d,def};
                mMT = mean(mT);
                stdMT = std(mT);
                
                maxS = maxSTimes{s,st}{d,def};
                mMax = mean(maxS);
                stdMax = std(maxS);
                
                
                
                lagF = cell2mat(allTest{s,st}{d,def}.finalLags)' *-1;
                meanLF = mean(lagF);
                meanLFs{s,st} = meanLF;
                errLF = std(lagF) / sqrt(length(lagF(:,1)));
                tTime = allOrig{s,st}{d,def}.Times{1,1} - allOrig{s,st}{d,def}.Times{1,1}(1,1);
                
                
                errorbar(tTime, meanLF, errLF)
            end
            
            
            line([0 tTime(end)],[0 0],'Color','r')
            stdv = stdvs(st);
            line([mSS mSS],[-max(deformations)*stdv max(deformations)*stdv],'Color','b')
            line([mAS mAS],[-max(deformations)*stdv max(deformations)*stdv],'Color','b')
            line([mMT mMT],[-max(deformations)*stdv max(deformations)*stdv],'Color','b')
            line([mMax mMax],[-max(deformations)*stdv max(deformations)*stdv],'Color','b')
            
            %         patch([mSS-stdSS mSS+stdSS mSS+stdSS mSS-stdSS],[-stdv*def -stdv*def stdv*def stdv*def],'c','EdgeColor','none','FaceAlpha',.3)
            %
            %         patch([mAS-stdAS mAS+stdAS mAS+stdAS mAS-stdAS],[-stdv*def -stdv*def stdv*def stdv*def],'c','EdgeColor','none','FaceAlpha',.3)
            %
            %         patch([mMT-stdMT mMT+stdMT mMT+stdMT mMT-stdMT],[-stdv*def -stdv*def stdv*def stdv*def],'c','EdgeColor','none','FaceAlpha',.3)
            
            if s == 1
                title(['Group: ', num2str(st)])
            end
            
            ylim([-max(deformations)*stdv max(deformations)*stdv])
            xlim([0 tTime(end)])
            ylabel(['Session: ',num2str(sessions(s))])
            
            ax = gca;
            ax.XTick = [mSS mMax mAS mMT];
            ax.XTickLabel = {'1','2','3','4'};
            
            
            sb = sb + 1;
        end
    end
end


%%% Look at test predictions, lag and no lag

% get test predictions for all sessions



for s = 1:length(allOrig(:,1))
    for st = 1:length(allOrig(1,:))
        for d = 1:length(allOrig{s,st}(:,1))
            
            allSpds = allOrig{s,st}{d,1}.Speeds;
            for rs = 1:length(allTest{s,st}{d,def}.testTrials{1,2}(1,:))
                teTrials = allTest{s,st}{d,def}.testTrials{1,2}{1,rs};
                teSpds{s,st}{d}{rs} = cell2mat(allSpds(teTrials));
            end
            
            for def = 1:length(allOrig{s,st}(1,:))
                tePredsLag{s,st}{d,def} = allTest{s,st}{d,def}.Preds{1,2};
                tePredsNoLag{s,st}{d,def} = allTest{s,st}{d,def}.Preds{1,1};
            end
        end
    end
end



s = 2;
st = 1;
d = 3;



smplSize = 1;
sample = randperm(length(tePredsLag{s,st}{d,1}(1,:)),smplSize);
gSpds = teSpds{s,st}{d};
gNoLags = tePredsNoLag{s,st}{d,1};

hold on
plot(gSpds{sample},'b')
plot(gNoLags{sample},'r')

for def = 1:length(tePredsLag{s,st}(1,:))
    gLags = tePredsLag{s,st}{d,def};
    plot(gLags{sample})
    hold on
end
 
xlim([0 length(gSpds{sample})])
ylim([0 max(gSpds{sample})+5])
    
    

%%% Look at train predictions, lag and no lag


% Plot all example resample trials


s = 1;
d = 2;
rsp = randi(100);


% Plot lag funtions and predictions for all smoothings and deformations of example session, group
% and resample
for st = 1:length(trSpds(1,:))
    
    
    spds = trSpds{s,st}{d,1}{1,rsp};
    pNoLag = trPredsNoLag{s,st}{d,1}{1,rsp};
    
    tms = allOrig{s,st}{d,1}.Times{1,1} - allOrig{s,st}{d,1}.Times{1,1}(1,1);
    tNum = length(spds) / length(tms);
    
    totData = 1:length(spds);
    nSamp = 16;
    
    totSamp = linspace(totData(1),totData(end), nSamp*tNum);
    
    
    % trialsSamp = reshape(totSamp,nSamp,tNum)';
    % trialsSamp = trialsSamp - repmat(trialsSamp(:,1),1,length(trialsSamp(1,:)));
    % trialsTime = trialsSamp * 0.02;
    % trialsTime = round(reshape(trialsTime',1,nSamp*tNum),1);
    
    
    figure('Units','centimeters','Position',[10,15,15,7])
    
    subplot(1,2,2)
    hold on
    plot(spds)
    plot(pNoLag)
    
    for def = 1:length(trSpds{s,st}(1,:))
        subplot(1,2,1)
        hold on
        lags = cell2mat(allTest{s,st}{d,def}.finalLags)' *-1;
        plot(lags(rsp,:))
        ax = gca;
        ax.XTick = 5:50:205;
        title(['Kernel Width',num2str(stdvs(st))])
        xlim([0 length(lags)])
        
        
        subplot(1,2,2)
        hold on
        pLag = trPredsLag{s,st}{d,def}{1,rsp};
        plot(pLag)
        xlim([2845 3081])
        ax = gca;
        ax.XTick = 2850:50:3050;
        ax.XTickLabel = {'5','55','105','155','205'};
        %     ax.XTickLabel = trialsTime;
    end
    
    
end




% Plot speed and predictions for all smoothings and deformations of example trial(defined by the xlim)
% and resample


for st = 1:length(trSpds(1,:))
    
    
    subplot(2,2,st)
    spds = trSpds{s,st}{d,1}{1,rsp};
    pNoLag = trPredsNoLag{s,st}{d,1}{1,rsp};
    
    hold on
    plot(spds,'b')
    plot(pNoLag,'r')
    
    for def = 1:3
     pLag = trPredsLag{s,st}{d,def}{1,rsp};
     plot(pLag)
    end
    
    xlim([2845 3081])
    title(['Kernel Width',num2str(stdvs(st))])
    
    ax = gca;
    ax.XTick = 2850:50:3050;
    ax.XTickLabel = {'0.1','1.1','2.1','3.1','4.1'};
    xlabel('Time(s)')
    ylabel('Speed')
    
end



% Plot lag functions, calculated on the sample from where the trial above was chosen, for the different smoothings deformations  

for st = 1:length(trSpds(1,:))
    
    
    subplot(2,2,st)
    
    for def = 1:3
     hold on
        lags = cell2mat(allTest{s,st}{d,def}.finalLags)' *-1;
        plot(lags(rsp,:))
    end
    
    title(['Kernel Width',num2str(stdvs(st))])
    
    ax = gca;
    ax.XTick = 5:50:205;
    ax.XTickLabel = {'0.1','1.1','2.1','3.1','4.1'};
    xlabel('Time(s)')
    ylabel('Lag')
    
end









%%% Plot individual trials of example session, smoothing, deformation ,group and sample 



s = 1;
st = 1;
rsp = randi(100);
d = 3;

tLen = length(allOrig{s,st}{d,1}.Speeds{1,1});
spds = trSpds{s,st}{d,1}{1,rsp};
spds = reshape(spds,tLen,length(spds)/tLen);

pNoLag = trPredsNoLag{s,st}{d,1}{1,rsp};
pNoLag = reshape(pNoLag,tLen,length(pNoLag)/tLen);



% Plot all speeds and non lagged speed predictions 

figure
plot(spds,'b')
hold on
plot(pNoLag,'r')


% Plot mean speeds and non lagged speed predictions 

figure
plot(mean(spds,2),'b')
hold on
plot(mean(pNoLag,2),'r')


% Plot mean accelerations and non lagged accelaration predictions 

figure
plot(diff(mean(spds,2)),'b')
hold on
plot(diff(mean(pNoLag,2)),'r')
line([0 length(pNoLag(:,1))],[0 0],'Color','r')


% Plot lags calculated for the sample used above and bellow

figure
lags = cell2mat(allTest{s,st}{d,1}.finalLags)' *-1;
plot(lags(rsp,:))





% Plot real and predicted accellerations in all trials  

figure
for i=1:length(pNoLag(1,:))
    
    subplot(5,7,i)
    plot(diff(spds(:,i)),'b')
    hold on
    plot(diff(pNoLag(:,i)),'r')
    line([0 length(pNoLag(:,1))],[0 0],'Color','k')
    line([6 6],[-3 3],'Color','k')
    line([27 27],[-3 3],'Color','k')
    axis tight
end



% Plot real and predicted speeds in all trials  
figure
for i=1:length(pNoLag(1,:))
    
    subplot(5,7,i)
    plot(spds(:,i),'b')
    hold on
    plot(pNoLag(:,i),'r')
    line([0 length(pNoLag(:,1))],[0 0],'Color','k')
    line([6 6],[0 70],'Color','k')
    line([27 27],[0 70],'Color','k')
    axis tight
end





%%% Test stuff with kernels



deformations = [0.5 0.75 1];
stdv = 0.14; % Stdv of kernel used to smooth data
sampling = 0.1 / 5;





%%% gaussian Kernel

n_const = 1/(stdv*(2*pi)^0.5);
a_exp =-1/2*((kernRange - 0.03)/stdv).^2;
kernel1 = n_const*exp(a_exp);

%%% gaussian derivative
n_const = -(kernRange)/(stdv^3*(2*pi)^0.5);
diffKernel = n_const.*exp(a_exp);




[kernel,diffKernel] = GenerateKernelNew(sampling,stdv);
tms = [0:length(kernel)-1]*sampling;

for d = 1:length(deformations)
    subplot(2,2,d)
    plot(tms,kernel)
    hold on
    plot(tms,kernel-diffKernel*stdv*deformations(d))
    hold on
    plot(tms+stdv*deformations(d),kernel)
end




for st = 1:4
    for def = 1:3
        lags = cell2mat(allTest{s,st}{d,def}.finalLags)';
        lags1{st,def} = lags(rsp,:);
    end
end

dummy = zeros(1,length(lags1{st,def}));


for i = 1:length(dummy)
    spk = rand(1);
    
    if i<length(dummy)*0.25
        
        if spk >= 0.7
            dummy(i) = 1;
        end
        
    elseif i>=length(dummy)*0.25 && i<=length(dummy)*0.75
        if spk >= 0.4
            dummy(i) = 1;
        end
        
    elseif i>length(dummy)*0.75
        
        if spk >= 0.9
            dummy(i) = 1;
        end
    end
   
end


kernRange = (0-5*stdv-sampling):sampling:(0+5*stdv+sampling);

n_const = 1/(stdv*(2*pi)^0.5);
a_exp =-1/2*((kernRange)/stdv).^2;
kernel = n_const*exp(a_exp);
n_const = -(kernRange)/(stdv^3*(2*pi)^0.5);
diffKernel = n_const.*exp(a_exp);

tms = [0:length(dummy)-1]*sampling;

for i = 1:length(deformations)
    
    n_const = 1/(stdv*(2*pi)^0.5);
    a_exp =-1/2*((kernRange - stdv*deformations(i))/stdv).^2;
    kernel1 = n_const*exp(a_exp);
    
    subplot(2,2,i)
    plot(tms,conv(dummy,kernel,'same'))
    hold on
    plot(tms,conv(dummy,kernel-diffKernel*stdv*deformations(i),'same'))
    hold on
    hold on
    plot(tms,conv(dummy,kernel1,'same'))
    
end



plot(conv(dummy,kernel))
hold on
plot(conv(dummy,kernel-diffKernel*stdv*deformations(1)))
hold on
plot(conv(dummy,kernel-diffKernel*stdv*deformations(2)))
hold on
plot(conv(dummy,kernel-diffKernel*stdv*deformations(3)))
hold on



stdv = 0.3;
sampling = 0.1 / 5;

kernRange = (0-5*stdv-sampling):sampling:(0+5*stdv+sampling);

n_const = 1/(stdv*(2*pi)^0.5);
a_exp =-1/2*((kernRange)/stdv).^2;
kernel = n_const*exp(a_exp);
n_const = -(kernRange)/(stdv^3*(2*pi)^0.5);
diffKernel = n_const.*exp(a_exp);


data = conv(dummy,kernel,'same');
diffData = conv(dummy,diffKernel,'same');


pNoLag = trPredsNoLag{1,st}{2,1}{1,rsp}(2845 : 3081);




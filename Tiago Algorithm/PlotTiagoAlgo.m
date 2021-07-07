clear all

sessions = 13; % Sessions to load data from
stdvs = 0.2; % Stdv of kernel used to smooth data
deformations = [0.5 0.75 1]; % Parameter of allowed deformation
tSelectMethod = 1; % Method used to select similar trials: 0 - manual; 1 - based on correlations; 2 - based on vector similarity
delGroup = 1;
intDels = 4; % Delays on which the analysis was run
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


data = load('/Volumes/GoogleDrive/O meu disco/Paper/MFunc/Paper-Code-1/Tiago Algorithm/test.mat');
lags_1 = data.lags;


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









lagF = cell2mat(allTest{1,1}{1,1}.finalLags)' *-1;
mLagF = mean(lagF);

plot((0:197)*0.02,mLagF)

hold on

plot((0:197)*0.02,(t_lags.lags_1(5,:))/1000)


figure
sb = 1;
st = 1;


for st = 1:length(allTest(1,:))
    figure
    sb = 1;
    for s = 1:length(allTest(:,1))
        for d = 1:length(allTest{s,st}(:,1))
            
            for def = 1%:length(allTest{s,st}(1,:))
                
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
                hold on
            end
            
            hold on
            
            plot(tTime,(lags_1(3,:)*-1)/1000)
            
            
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
                title('Group: 4')
            end
            
            %ylim([-max(deformations)*stdv max(deformations)*stdv])
            ylim([-0.15 0.15])
            xlim([0 tTime(end)])
            ylabel('Lags (s)')
            
            ax = gca;
            ax.XTick = [mSS mMax mAS mMT];
            ax.XTickLabel = {'1','2','3','4'};
            
            
            sb = sb + 1;
        end
    end
end



%  

data = load('/Volumes/GoogleDrive/O meu disco/Paper/MFunc/Paper-Code-1/Tiago Algorithm/test.mat');

coefs = data.coefs;
lags = data.lags;
neurons = data.neurons;
speeds = data.speeds;


for t = 1:35
allTrialsNoLag{t,1} =  reshape(neurons(1,:,t,:),198,96);
allTrialsLag{t,1} =  reshape(neurons(3,:,t,:),198,96);
end

NoLagN = cell2mat(allTrialsNoLag);
LagN = cell2mat(allTrialsLag);

betasNoLag = coefs(1,:)';
predsNoLag = NoLagN * betasNoLag;

betasLag = coefs(3,:)';
predsLag = LagN * betasLag;

plot(speeds)
hold on
plot(predsNoLag)
hold on
plot(predsLag)




for i = 1:3
    
n = reshape(neurons(i,:,1,:),198,96);    
plot(n(:,1))
hold on
end







edit clear all

sessions = [1 13]; % Sessions to load data from
stdvs = [0.14 0.2]; % Stdv of kernel used to smooth data
deformations = [0.5 0.75 1]; % Parameter of allowed deformation
tSelectMethod = 1; % Method used to select similar trials: 0 - manual; 1 - based on correlations; 2 - based on vector similarity
delGroup = 1;
intDels = [1 2 3 4]; % Delays on which the analysis was run
intCond = 1; % Trial's condition: 1-hits; 2-misses; 3-correct rejections; 4-false alarms
p1 = 0.25; %Cut trials p1 seconds after the higher area time across trials
p2 = 0.5;  %Cut trials p2 seconds before the lowest moving after stopping time across trials


for s = 1:length(sessions)
    for st = 1:length(stdvs)
        for d = 1:length(intDels)
            for def = 1:length(deformations)
                
                disp(['loading deformation ',num2str(deformations(def)),' from delay ', num2str(intDels(d)),' of kernel width ',num2str(stdvs(st)),' session ', num2str(sessions(s))])
                [Name] = SessionNumber2Name(sessions(s));
                load(['/Volumes/Cerebro/Recording Data/',Name,'/Predict from Neurons/RigidShiftResults Trials tsMethod_',num2str(tSelectMethod),' Cond_',num2str(intCond),' Delay_',num2str(delGroup),'_',num2str(intDels(1,d)),' P1_',num2str(p1),' P2_',num2str(p2),' Sigma_',num2str(stdvs(st)),' Deform_',num2str(deformations(def)),'.mat']);
                
                allRMSE{s,st}{d,def} = RMSE;
                allBetas{s,st}{d,def} = Betas;
                allPreds{s,st}{d,def} = Preds;
                allSpds{s,st}{d,def} = Spds;
                allLags{s,st}{d,def} = Lags;
                
            end
        end
    end
end





% Plot test performance and lag no lag difference



clrs = {'r','b','c'};
jitters = [-0.25 0 0.25];
sub=1;
for s = 1:length(sessions)
    for st = 1:length(stdvs)
        subplot(2,2,sub)
        
        for d = 1:length(intDels)
            for def = 1:length(deformations)
                
                lagRMSETe = allRMSE{s,st}{d,def}.Te(:,2);
                noLagRMSETe = allRMSE{s,st}{d,def}.Te(:,1);
                
                mLagE = mean(lagRMSETe);
                mNoLagE = mean(noLagRMSETe);
                
                eLagE = std(lagRMSETe)/sqrt(length(lagRMSETe));
                eNoLagE = std(noLagRMSETe)/sqrt(length(noLagRMSETe));
                
                
                if def == 1
                    
                    patch([d+jitters(1) d+jitters(3) d+jitters(3) d+jitters(1)],[mNoLagE-eNoLagE mNoLagE-eNoLagE mNoLagE+eNoLagE mNoLagE+eNoLagE],'y','EdgeColor','none')
                    hold on
                    
                end
                
                errorbar(d+jitters(def),mLagE,eLagE,'o','MarkerFaceColor',clrs{def},'MarkerEdgeColor','none','Color','k')
                hold on
               
                
            end
        end
        ylim([-Inf 11.5])
        xlabel('delay groups')
        ylabel('RMSE')
        title(['kWidth: ',num2str(stdvs(st))])
        ax = gca;
        ax.XTick = [1 2 3 4];
        sub = sub+1;
    end
end



% Plot  lag - nolag difference as percentage of total error




clrs = {'r','b','c'};
jitters = [-0.25 0 0.25];
sub=1;
for s = 1:length(sessions)
    for st = 1:length(stdvs)
        subplot(2,2,sub)
        
        for d = 1:length(intDels)
            for def = 1:length(deformations)
                
                lagRMSETe = allRMSE{s,st}{d,def}.Te(:,2);
                noLagRMSETe = allRMSE{s,st}{d,def}.Te(:,1);
                errDiff = ((noLagRMSETe - lagRMSETe)./noLagRMSETe)*100;
                
                
                mErrDiff = mean(errDiff);
                
                
                eErrDiff = std(errDiff)/sqrt(length(errDiff));
                
                errorbar(d+jitters(def),mErrDiff,eErrDiff,'o','MarkerFaceColor',clrs{def},'MarkerEdgeColor','none','Color','k')
                hold on
               
                
            end
        end
        ylim([0 mErrDiff+eErrDiff+2])
        xlabel('delay groups')
        ylabel('Error Diff %')
        title(['kWidth: ',num2str(stdvs(st))])
        sub = sub+1;
    end
end





% Plot rigid lag lag - approx lag difference 




clrs = {'r','b','c'};
jitters = [-0.25 0 0.25];
sub=1;
for s = 1:length(sessions)
    for st = 1:length(stdvs)
        subplot(2,2,sub)
        
        for d = 1:length(intDels)
            for def = 3
                
                lagRMSETe = allRMSE{s,st}{d,def}.Te(:,2);
                origLagRMSETe = allRMSE{s,st}{d,def}.Te(:,4);
                
                mLagE = mean(lagRMSETe);
                mOrigLagE = mean(origLagRMSETe);
                
                eLagE = std(lagRMSETe)/sqrt(length(lagRMSETe));
                eOrigLagE = std(origLagRMSETe)/sqrt(length(origLagRMSETe));
                
                errorbar([d+jitters(1) d+jitters(3)],[mLagE mOrigLagE],[eLagE eOrigLagE],'o')
                hold on 
               
                
            end
        end
        %ylim([0 mErrDiff+eErrDiff+2])
        xlabel('delay groups')
        ylabel('RMSE')
        title(['kWidth: ',num2str(stdvs(st))])
        ax = gca;
        ax.XTick = [1 2 3 4];
        sub = sub+1;
    end
end


rsp = 4;


figure
plot(allSpds{1,2}{1,3}.Te(:,rsp))
hold on
plot(allPreds{1,2}{1,3}.NoLagTe(:,rsp))
hold on
plot(allPreds{1,2}{1,3}.LagTe(:,rsp))
hold on
plot(allPreds{1,2}{1,3}.Orig.LagTe(:,rsp))



figure

plot(allSpds{1,2}{1,3}.Te(:,rsp) - allPreds{1,2}{1,3}.NoLagTe(:,rsp))
hold on
plot(allSpds{1,2}{1,3}.Te(:,rsp) - allPreds{1,2}{1,3}.LagTe(:,rsp))
hold on
plot(allSpds{1,2}{1,3}.Te(:,rsp) - allPreds{1,2}{1,3}.Orig.LagTe(:,rsp))







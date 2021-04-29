clear all



sessions = [1 13]; % Sessions to load data from
stdv = 0.14; % Stdv of kernel used to smooth data
deformations = [0.5 0.75 1]; % Parameter of allowed deformation
tSelectMethod = 1; % Method used to select similar trials: 0 - manual; 1 - based on correlations; 2 - based on vector similarity
delGroup = [1 3];
intDels = {[1 2 3 4];[1 2 3 4 5]}; % Delays on which the analysis was run
intCond = 1; % Trial's condition: 1-hits; 2-misses; 3-correct rejections; 4-false alarms
p1 = 0.25; %Cut trials p1 seconds after the higher area time across trials
p2 = 0.5;  %Cut trials p2 seconds before the lowest moving after stopping time across trials


%********************  Load and process data ******************

for s = 1:length(sessions)
    for g = 1:length(delGroup)
        for d = 1:length(intDels{g})
            for def = 1:length(deformations)
                
                disp(['loading deformation ',num2str(deformations(def)),' from delay ', num2str(intDels{g}(d)),' of group ',num2str(delGroup(g)),' session ', num2str(sessions(s))])
                [Name] = SessionNumber2Name(sessions(s));
                load(['/Volumes/Cerebro/Recording Data/',Name,'/Predict from Neurons/SpdPredLagsBootCheck Trials tsMethod_',num2str(tSelectMethod),' Cond_',num2str(intCond),' Delay_',num2str(delGroup(g)),'_',num2str(intDels{g}(1,d)),' P1_',num2str(p1),' P2_',num2str(p2),' Sigma_',num2str(stdv),' Deform_',num2str(deformations(def)),'.mat'],'testData','origData','valData');
                
                allOrig{s,g}{d,def} = origData;
                allTest{s,g}{d,def} = testData;
                allVal{s,g}{d,def} = valData;
            end
        end
    end
end
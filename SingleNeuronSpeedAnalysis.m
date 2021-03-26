clearvars -except  data simpleData
close all

processRaw = 1; % if process == 1: process raw data again.  If process == 0 load pre-processed data.

if ~exist('data','var') && ~exist('simpleData','var')
    
    if processRaw == 1
        
        load('/Volumes/Cerebro/Recording Data/SessionsMap.mat')
        
        conditions=[1 3]; % Trial conditions to get
        
        mice = {1; 5}; % Mice from which to pick neurons
        days = {1 ; 1}; % Days of the recording sessions, from the picked mice, from which to pick neurons
        intSess = [mice days]; % Structure with all the info
        
        cutPoint = 102; 
        spdSmt = 500; % Speed smoothing factor
        stdv = 0.14;
        
        sess = 1;
        tT = tic;
        for m = 1:length(sessionsMap(:,1)) % for all the picked mice
            for s = 1:length(sessionsMap{m,2}(1,:)) % for all the selected sessions
                
                Name = sessionsMap{m,2}{1,s} % Name of the session
                mouseNumber = str2num(Name(3:4)); % number of the mouse being loaded
                
                
                
                % load raw neural data
                disp('Loading Raw Neurons')
                load(['/Volumes/Cerebro/Recording Data/',Name,'/Neurons.mat'])
                
                data{sess,1} = Neurons(:,2:3);
                
                
                % load session behavioral data
                disp('Loading Behavior data')
                load(['/Volumes/Cerebro/Recording Data/',Name,'/OrganizedSessionInfo.mat'],'validSession','organizedTrials')
                
                
                rwddSound=mod(mouseNumber,2)+1; % define rewarded and and non rewarded sound based on animal number
                nRwddSound=2-mod(mouseNumber,2);
                
                conditions=[rwddSound 1; rwddSound 0 ;nRwddSound 1;nRwddSound 2]; % Behavioral conditions to analyse
                
                % Add sound start info to organized trials and get rid of souns start 0 trials bug on data
                organizedTrials(9,:) = validSession(9,:);
                unknownDelays = find(cell2mat(organizedTrials(9,:)) == 0);
                organizedTrials(:,unknownDelays)=[];
                
                disp('Smoothing speeds and dists')
                
                % Smooth speeds and dists
                for t = 1:length(organizedTrials(1,:))
                    spn = spdSmt/length(organizedTrials{1,t});
                    
                    organizedTrials{1,t} = smooth(organizedTrials{1,t},spn,'lowess');
                    organizedTrials{1,t}(organizedTrials{1,t} < 0) = 0;
                    organizedTrials{3,t} = smooth(organizedTrials{3,t},spn,'lowess');
                    
                end
                
                disp('Organizing data by delay and condition')
                % Organize trials by condition and delay
                delays=unique(cell2mat(organizedTrials(9,:)));% Sound Starts
                
                for c = 1:length(conditions(:,1))
                    for d = 1:length(delays)
                        filter = cell2mat(organizedTrials(9,:)) == delays(d) & cell2mat(organizedTrials(7,:)) == conditions(c,1) & cell2mat(organizedTrials(8,:)) == conditions(c,2);
                        if sum(filter) == 0
                            data{sess,2}{c,d} = [];
                            data{sess,3}{c,d} = {};
                            data{sess,4}{c,d} = {};
                        else
                            data{sess,2}{c,d} = organizedTrials(:,filter);
                            % Select and organize each neuron spike times according to the trial they belong to
                            for t = 1:length(data{sess,2}{c,d}(1,:))
                                
                                ts = data{sess,2}{c,d}{4,t}(1,1);
                                
                                if c == 1 || c == 4
                                    dist = data{sess,2}{c,d}{3,t} - data{sess,2}{c,d}{3,t}(1,1);
                                    [~,idx] = min(abs(dist - cutPoint));
                                    te = data{sess,2}{c,d}{2,t}(idx);
                                else
                                    te = data{sess,2}{c,d}{4,t}(1,2);
                                end
                                
                                for n = 1:length(data{sess,1}(:,1))
                                    [smtSpikes,smtTS] = GaussianSmoothNeuronsPaper(data{sess,1}{n,1},stdv,[ts te]);
                                    data{sess,3}{c,d}{n,t} = smtSpikes';
                                    if n == 1
                                        spd = data{sess,2}{c,d}{1,t};
                                        tm = data{sess,2}{c,d}{2,t}(2:end);
                                        dowSpds = zeros(1,length(smtTS));
                                        for i = 1:length(smtTS)
                                            [~,idx] = min(abs(tm - smtTS(i)));
                                            dowSpds(i) = spd(idx);
                                        end
                                        data{sess,4}{c,d}{n,t} = dowSpds;
                                    end
                                end
                                
                            end
                        end
                    end
                    data{sess,2}{c,d+1} = horzcat(data{sess,2}{c,:});
                    data{sess,3}{c,d+1} = horzcat(data{sess,3}{c,:});
                    data{sess,4}{c,d+1} = horzcat(data{sess,4}{c,:});
                end
                
                
                clearvars -except m s  intSess sessionsMap data spdSmt stdv tT sess cutPoint
                disp(['Finished processing session ', num2str(sess),' of 16, elapsed time: ',num2str(toc(tT)/60),' minutes.'])
                sess = sess+1;
            end
        end
        
        simpleSpd = vertcat(data{:,4});
        simpleSpd = simpleSpd(:,12);
        
        simpleNeurs = vertcat(data{:,3});
        simpleNeurs = simpleNeurs(:,12);
        
        simpleData = [simpleNeurs simpleSpd];
        
        save('/Volumes/GoogleDrive/O meu disco/paper/SNSpeedsSimple.mat','simpleData')
        save('/Volumes/GoogleDrive/O meu disco/paper/SNSpeedsComplete.mat','data','-v7.3')
        
    elseif processRaw == 0
        
        disp('Loading pre-processed data')
        load('/Volumes/GoogleDrive/O meu disco/paper/SNSpeedsSimple.mat')
    end
end




cond = [1 3];

% trial correlations

for s = 1:16
    for c = 1:length(cond)
        idx = 4*(s-1)+cond(c);
        
        spds = simpleData{idx,2};
        frs = simpleData{idx,1};
        neursPerSess(s) = length(frs(:,1));
        
        for t = 1:length(spds)
            tSpd = spds{1,t};
            tempFr = cellfun(@transpose,frs(:,t),'UniformOutput',0)';
            tFr = cell2mat(tempFr);
            totMatTrial = [tSpd' tFr];
            crcf = corrcoef(totMatTrial);
            rTrial{s,c}(t,:) = crcf(2:end,1);
        end
        
    end
end


% total correlations 

for s = 1:16
    for c = 1:length(cond)
        idx = 4*(s-1)+cond(c);
    
        spd = cell2mat(simpleData{idx,2})';
        fr = cell2mat(simpleData{idx,1})';
        
        totMat = [spd fr];
        crcf = corrcoef(totMat);
        r{s,c}(:,1) = crcf(2:end,1);
        allSpds{s,c} = spd;
        allFrs{s,c} = fr;
        
    end
end





% histogram of total correlations by session

figure
bins = -1:0.05:1;
clrs = {'g','r'};

f = figure;
f.Units = 'centimeters';
f.Position = [10,10,12,12];
for s = 1:length(r(:,1))
    subplot(4,4,s)
    hold on
    for c = 1:length(r(1,:))
        h = histogram(r{s,c},bins);
        h.FaceColor = clrs{c};
        h.EdgeColor = clrs{c};
        xlim([-1 1])
    end
    ax = gca;
    ax.FontSize = 8;
end



% histogram of total across sessions
f = figure;
f.Units = 'centimeters';
f.Position = [10,10,7,7];
for c = 1:length(r(1,:))
    hold on
    allR(:,c) = cell2mat(r(:,c));
    h = histogram(allR(:,c),bins);
    h.FaceColor = clrs{c};
    h.EdgeColor = clrs{c};
    xlim([-1 1])
end
xlabel('Corr Coef')
ylabel('# of neurons')
ax = gca;
ax.FontSize = 8;



% Select example neurons

[extremeCorr,extremeCorrIdx] = max(allR); % max positive corr for both cond
[extremeACorr,extremeACorrIdx] = min(allR); % max negative corr for both cond


extremeCorrs = [extremeCorr' extremeACorr']; 
extremeCorrsIdx = [extremeCorrIdx' extremeACorrIdx'];

cumNPerSess = cumsum(neursPerSess);

f = figure;
f.Units = 'centimeters';
f.Position = [10,10,8,8];
sp = 1;
for  i = 1:length(extremeCorrsIdx(:,1))
    for j = 1:length(extremeCorrsIdx(1,:))
        
        nIdx = extremeCorrsIdx(i,j);
        baseIdx = find(cumNPerSess > nIdx,1,'first');
        
        if baseIdx > 1
            base = cumNPerSess(baseIdx-1);
            neur = nIdx - base;
        else
            neur = nIdx;
        end
        
        if j == 1
            [tCorr,tNum] = max(rTrial{baseIdx,i}(:,neur));
        else
            [tCorr,tNum] = min(rTrial{baseIdx,i}(:,neur));
        end
        
        idx = 4*(baseIdx-1)+cond(i);
         
        subplot(2,2,sp)
        fr = simpleData{idx,1}{neur,tNum};
        spd = simpleData{idx,2}{1,tNum};
        plot(fr,clrs{i})
        hold on
        plot(spd,'b')
        title(num2str(extremeCorrs(i,j)))
        xlim([0 length(spd)])
        
        sp = sp+1;
        
        ax = gca;
        ax.XTick = [];
        ylabel('Speed / FR')

    end
end



% scatter both conditions corr coef
f = figure;
f.Units = 'centimeters';
f.Position = [10,10,6,6];

scatter(allR(:,1),allR(:,2),5,allR(:,1).*allR(:,2),'filled')
colormap jet
line([-1 1],[0 0],'Color','r','LineStyle','--')
line([0 0],[-1 1],'Color','r','LineStyle','--')
xlim([-1 1])
ylim([-1 1])
xlabel('Stop Trials Corr Coef')
ylabel('No Stop Trials Corr Coef')
ax = gca;
ax.FontSize = 6;






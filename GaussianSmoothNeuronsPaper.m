function [smtSpikes,timeStamps] = GaussianSmoothNeuronsPaper(spikeTimes,stdv,tlim)
clearvars -except spikeTimes stdv tlim

% Create vector of spike counts
sampling = stdv/5;
tlimOff = 3.5;
tlim = round(tlim,4);
extTLim = [tlim(1)-tlimOff*stdv tlim(2)+tlimOff*stdv];
trialVec = extTLim(1) :sampling: extTLim(2);

spikes=zeros(1,length(trialVec));
spkt=spikeTimes(spikeTimes>extTLim(1,1) & spikeTimes<=extTLim(1,2));

for iSpk=1:length(spkt) % go over spikes in the spike train of cell i
    tSpk=spkt(iSpk);
    [~,currInd] = min(abs(trialVec-tSpk));
    spikes(currInd) = spikes(currInd) + 1;
end


% Create smoothing kernel

[kernel,~] = GenerateKernelNew(sampling,stdv);
kernel = fliplr(kernel);



% Convolve vector of spike counts with kernel


padding = zeros(1,length(-8*stdv:sampling:-sampling));
spikesVec = [padding spikes padding];
timeVec = [padding trialVec padding];


spkS = length(padding) + 1;
spkE = spkS + length(spikes) - 1;
tStart = tlim(1);
tEnd = tlim(end);
smtSpikes=zeros(1,length(spikesVec)); %initialize output
rangeInd= (length(kernel)-1)/2; %number of dt's that fit into 6 stdv (rounded upwards)

for  i = 1:length(spikesVec)
    if i >= spkS && i <= spkE
        %if timeVec(i) >= tStart && timeVec(i) <= tEnd
            smtSpikes(i) = spikesVec(1,i-rangeInd:i+rangeInd) * kernel';
        %end
    end
end

[~,sIdx] = min(abs(timeVec - tStart));
[~,eIdx] = min(abs(timeVec - tEnd));

% if length(sIdx:eIdx) < length(tS)
%     eIdx = eIdx + 1;
% elseif length(sIdx:eIdx) > length(tS)
%     eIdx = eIdx - 1;
% end

smtSpikes = smtSpikes(sIdx : eIdx)';
timeStamps = timeVec(sIdx : eIdx)';
end








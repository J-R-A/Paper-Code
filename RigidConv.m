function [smtSpikes1,timeStamps] = RigidConv(spikes,tS,tsSpk,stdv,sampling,lags)

clearvars -except spikes tS tsSpk stdv sampling lags

spikesVec = [zeros(1,length(-8*stdv:sampling:-sampling)) spikes' zeros(1,length(sampling:sampling:+8*stdv))];
timeVec = [zeros(1,length(-8*stdv:sampling:-sampling)) tsSpk zeros(1,length(sampling:sampling:+8*stdv))];

spkS = length(-8*stdv:sampling:-sampling) + 1;
spkE = spkS + length(spikes) - 1;
tStart = tS(1);
tEnd = tS(end);
smtSpikes=zeros(1,length(spikesVec));

n_const = 1/(stdv*(2*pi)^0.5);
kernRange = (0-5*stdv-sampling):sampling:(0+5*stdv+sampling);
idx = 1;

for  i = 1:length(spikesVec)
    if i >= spkS && i <= spkE
        
        if timeVec(i) >= tStart && timeVec(i) <= tEnd
            %disp('Using altered lag')
            l = lags(idx);
            idx = idx+1;
        else
            l = 0;
        end
        
        
        a_exp =-1/2*((kernRange - l)/stdv).^2;
        kernel = fliplr(n_const*exp(a_exp));
        rangeInd= (length(kernel)-1)/2;
        
        smtSpikes(i) = spikesVec(1,i-rangeInd:i+rangeInd) * kernel';
        
    end
end

[~,sIdx] = min(abs(timeVec - tStart));
[~,eIdx] = min(abs(timeVec - tEnd));

if length(sIdx:eIdx) < length(tS)
    eIdx = eIdx + 1;
elseif length(sIdx:eIdx) > length(tS)
    eIdx = eIdx - 1;
end

smtSpikes1 = smtSpikes(sIdx : eIdx)';
timeStamps = timeVec(sIdx : eIdx)';


end
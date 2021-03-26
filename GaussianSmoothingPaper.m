function [smtSpikes,timeStamps] = GaussianSmoothingPaper(spikeTimes,sigma,sampling,tlim,kernel)


% Invert kernel for conv
kernel = fliplr(kernel);



% Extend tlims to avoid edge effects
extTlim(1) = tlim(1) - sigma*3.5;
extTlim(2) = tlim(2) + sigma*3.5;

%times where time series will be evaluated
tvec=extTlim(1)-8*sigma:sampling:extTlim(2)+8*sigma; 





csign=zeros(1,length(tvec)); %initialize output
spkt=spikeTimes(spikeTimes>=extTlim(1) & spikeTimes<=extTlim(2));
rangeInd= (length(kernel)-1)/2;

for iSpk=1:length(spkt) % go over spikes in the spike train of cell i
    tSpk=spkt(iSpk);
    
    [~,currInd] = min(abs(tvec-tSpk));
    
    csign(1,currInd-rangeInd:currInd+rangeInd)=csign(1,currInd-rangeInd:currInd+rangeInd)+kernel;
end

smtSpikes = csign(tvec>=tlim(1) & tvec<=tlim(2));
timeStamps=tvec(tvec>=tlim(1) & tvec<=tlim(2));

end








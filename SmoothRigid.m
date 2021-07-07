function [matSpikes] = SmoothRigid(tSpikes,stdv,tTms,tTmsSpk,sampling,lags)

%%% intermediate function to transform predictors matrices before feeding
%%% them to glmnet speed prediction functions

clearvars -except tSpikes stdv tTms tTmsSpk sampling lags


for n = 1:length(tSpikes(1,:))
    %disp(['Processing Neuron: ',num2str(n),' of ',num2str(length(tSpikes(1,:)))])
    for t = 1:length(tSpikes(:,1))
        %disp(['Processing trial ',num2str(t),' of ',num2str(length(tSpikes(:,1))),' from neuron ',num2str(n),' of ',num2str(length(tSpikes(1,:)))])
        
        [smtSpikes{t,n},~] = RigidConv(tSpikes{t,n},tTms{t,1},tTmsSpk{t,n},stdv,sampling,lags);
    end
end

matSpikes = cell2mat(smtSpikes);
end


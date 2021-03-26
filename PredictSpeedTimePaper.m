function [Laps, dummyLaps, allPredMat, allTrials, trialEvents, allDists, coefs, fitInfo] = PredictSpeedTimePaper(pMats,trialInfo,trialEvents,Speeds,Distances,elnetVal)
clearvars -except pMats trialInfo trialEvents Speeds Distances paramInfo shuffleTrials elnetVal

%%% Predicts Speed from neural activity




% Creates the final predictor and target matrices by essencially performing three basic operations:
% 1. intermingles trials of both conditions so that the model always has a balanced training set to work with
% 2. shuffles the order of the predictors relative to the targets if one
% wants to perform one of the controls
%3. Adds a extra predictor column which has the trials from all delays.





%%% Set random number seed (for repeatability)

nSeed=1; %
rng(nSeed,'twister');


%%% Set parameters and upload them into the options object if necessary

nOutFolds = 10; % folds of the outer nested cross-validation loop
nfolds=10;
options = glmnetSet;
options.alpha = elnetVal; % set regularization parameter alpha for glm fits



%%%%%%%%%%%%%%%% PREDICTION %%%%%%%%%%%%%%%%


% Prealocating variables


disp('Starting model fitting')

del=1:length(pMats(1,:));
for d=del
    
    
    predictorsMatrix=pMats{1,d}; % FR predictors
    out=Speeds{1,d}; % Target Speeds
    delTrials=trialInfo{1,d}; % Trial type info
    dsts = Distances{1,d};
    
    % just saving the actual data being used for sanity check...
    allPredMat{1,d}=predictorsMatrix;
    allTrials{1,d}=delTrials;
    allDists{1,d} = dsts;
    
    % Dimensions of pred matrix
    [nTrls,Npred]=size(predictorsMatrix);
    
    
    sz=size(out); % Make sure outcome is column vector
    if sz(2)>sz(1)
        out=out';
    end
    
    %%% Initialize prediction arrays
    predsSpeeds=zeros(nTrls,1);
    
    %%% Start outer Nested Cross Validation Loop to assess prediction
    nFold=floor(nTrls/nOutFolds);
    for icv=1:nOutFolds
        
        %%% Select trial indices for training set (iTrain) and test (iTest) set of this fold
        if icv<nOutFolds
            iTest=((icv-1)*nFold)+1:icv*nFold;
        else
            iTest=((icv-1)*nFold)+1:nTrls;
        end
        iTrain=1:nTrls;
        iTrain(iTest)=[];
        
        %%% Fit model using Cross Validation for current Training Set. This
        %%% means a best value of lambda will be found in each fold. We will
        %%% test the accuracy of the model (lambda,coefs) to predict the
        %%% remaining (test) data in each fold, and average across folds.
        
        
        cvFitSpeeds=cvglmnet(predictorsMatrix(iTrain,:),out(iTrain,:),'gaussian',options,[],nfolds,[],0); % Fits the model parameters and hyperparameters using "nfolds" cross validation
        predsSpeeds(iTest)=cvglmnetPredict(cvFitSpeeds,predictorsMatrix(iTest,:),'lambda_1se','response'); % uses the fitted model to predict the data left out in the nOutFolds cross validation
        
        % Saves the data for each nOutFold prediction
        %                 cvData{1,d}{1,icv}=cvFitSpeeds;
        %                 cvData{1,d}{2,icv}=predictorsMatrix(iTrain,:);
        %                 cvData{1,d}{3,icv}=out(iTrain,:);
        %                 cvData{1,d}{4,icv}=predictorsMatrix(iTest,:);
        %                 cvData{1,d}{5,icv}=out(iTest,:);
        %                 cvData{1,d}{6,icv}=predsSpeeds(iTest);
        
        
        disp(['Outer cross validation fold ' num2str(icv) ' from delay ' num2str(d) ' done!!'])
    end
    
    
    % for each delay saves all the predicted and real speeds
    dummyLaps{1,d}(:,1)=predsSpeeds;
    Laps{1,d}(:,1)=out;
    
    
    
    %%%%%%%%%%%%%%%% COEFFICIENTS %%%%%%%%%%%%%%%%
    %
    % % Now we run a single cross-validation sweep **ON THE WHOLE DATA** to get the coefficients
    % % (but of course, since we're using the whole data, we don't use this for doing prediction, because we would be overfitting).
    %
    %
    
    cvFitBothAll=cvglmnet(predictorsMatrix,out,'gaussian',options); % run cross val to find optimal lambda
    predsTrainSpeeds=cvglmnetPredict(cvFitBothAll,predictorsMatrix,'lambda_1se','response');
    
    
    iLmin=find(cvFitBothAll.glmnet_fit.lambda==cvFitBothAll.lambda_1se); % find out which of the lambdas used corresponds to lambda_min
    coefsBoth=cvFitBothAll.glmnet_fit.beta(:,iLmin); % from the glmnet fit inside the cvglmnet fit extract the coefficients at lambda=lambda_min
    coefsBoth=[cvFitBothAll.glmnet_fit.a0(iLmin) ; coefsBoth]; % include the offset in the predictor array
    
    % For lambda plot
    lamsBothAll=cvFitBothAll.lambda;
    cvMErrBothAll=cvFitBothAll.cvm;
    
    
    % Save the data from fitting the entire data set
    
    coefs{1,d}(:,1)=coefsBoth;
    fitInfo{1,d}=cvFitBothAll;
    
    
end
end
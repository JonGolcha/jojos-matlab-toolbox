function [TE] = proc_TE(DATA, SRATE)
% [TE] = proc_TE(DATA, SRATE)
%
% Computes transfer entropy. 
% (this function is a wrapper on top of (Montalto et al.)'s MUTE toolbox
%
% Inputs:
% - DATA    	--> 2D matrix with the data (P channels by N samples)
%                       OR 3D matrix (P channels, N samples, T trials).
% - SRATE    	--> sampling frequency
%
% Outputs:
% - PHASE    	--> 2D or 3D matrix of instantaneous phase values 
%                   (P channels by N samples)
%
%
% History
% Last version:  10/10/2016
% J.Chatel-Goldman, jonas.chatel.goldman(at)gmail.com

%% assess inputs
if(nargin < 3)  
    error('srate and filter specs must be provided'),
end
P = size(DATA, 1);
N = size(DATA, 2);
if ndims(DATA)==3
    T = size(DATA, 3);
end


%% Defining the experiment parameters
channels             = 1:5;
samplingRate         = 1;
pointsToDiscard      = 4500;
listRealization      = dir([dataDir [dataFileName '*' dataLabel '*' dataExtension]]);
autoPairwiseTarDriv  = [1 1 1 1 1 1 1 1];
handPairwiseTarDriv  = [0 0 0 0 0 0 0 0];

% PAY ATTENTION: if you are able to run the parallel session you can set numProcessors > 1
numProcessors	= 4;


%% Parameters for binnue algorithm
idTargets             = [ones(1,numSeries)*2 ones(1,numSeries)*6 ones(1,numSeries)*17];
idDrivers             = [[1 3:76] [1:63 65:76] [1:73 75 76]];
idOtherLagZero        = [3,5,8];
modelOrder            = 8;
multi_bivAnalysis     = 'multiv';
numQuantLevels        = 6;
entropyFun            = @evaluateNonUniformEntropy;
preProcessingFun      = @quantization;
secondTermCaseVect    = caseVect;%[1 1];
numSurrogates         = 100;
alphaPercentile       = 0.05;
******** Set the following fields together *******
genCondTermFun        = @generateConditionalTerm;%@generateConditionalTerm_selectionVar;%@generateCondTermLagZero
usePresent            = 0;
scalpConduction       = 0;
% **************************************************


%% Computing statistical methods...
[output1,params1]               = parametersAndMethods(listRealization,samplingRate,pointsToDiscard,channels,autoPairwiseTarDriv,...
                                      handPairwiseTarDriv,resultDir,dataDir,copyDir,numProcessors,...
                                      'linue',[],[],[],5,'multiv',5,5,'bayesian',@linearEntropy,[1 0],[1 1],@generateConditionalTerm,0,...
                                      'linnue',[],[],[],5,'multiv',@evaluateLinearNonUniformEntropy,[1 1],100,0.05,@generateConditionalTerm,0,...
                                      'binue',[],[],[],5,'multiv',6,@conditionalEntropy,@quantization,[1 0],[1 1],100,0.05,20,...
                                      @generateConditionalTerm,0,...
                                      'binnue',[],[],[],5,'multiv',6,@evaluateNonUniformEntropy,@quantization,[1 1],100,0.05,@generateConditionalTerm,0,0,...
                                      'neunetue',[],[],[],5,[1 1],'multiv',[],[],{@sigmoid @identity},30,0,4000,2/3,15,...
                                      valThreshold,@resilientBackPropagation,1.1,0.9,1,numHiddenNodes,100,20,0.05,@generateCondTerm,1,...
                                      'neunetnue',[],[],[],[],5,[1 0],[1 1],'multiv',[],[],{@sigmoid @identity},30,0,4000,threshold,2/3,15,...
                                      valThreshold,@resilientBackPropagation,1.1,0.9,1,numHiddenNodes,@generateConditionalTerm,1,...
                                      'nnue',[],[],[],5,'multiv',[1 1],100,'maximum',10,nnMexa64Path,mutePath,0.05,10,@generateConditionalTerm,0,...
                                      'nnnue',[],[],[],5,'multiv',[1 1],100,'maximum',10,@nearNeiConditionalMutualInformation,...
                                      @evalNearNeiTestSurrogates2rand,nnMexa64Path,mutePath,0.05,@generateConditionalTerm,0);
    

return,
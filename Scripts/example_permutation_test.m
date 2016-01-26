% This script shows how to perform a simple paired permutation test on two 
% vectors of observations (simulated). Results are then displayed  in the 
% form of a bar plot of surrogate and real values.
%
% History:
% --------
% First version: 2014-05-15
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com


clear all
close all
addpath(genpath([fileparts(pwd) '\JOJOS_MATLAB_TOOLBOX'])); % add all subfolders below main functions folder

% --- set some parameters
alpha   = 0.01;     % significance level
nPerm   = 10000;    % number of permutations (i.e., size of surrogate data)
nObs    = 100;     	% number of observations
mean_1  = 10;       % mean of first vector of observations
mean_2  = 12.5;       % mean of second vector of observations
var_1   = 5;        % variance of first vector of observations
var_2   = 5;        % variance of second vector of observations


% --- simulate two vectors of paired obervation
vec_OBS_1   = mean_1 + var_1*randn(1,nObs);
vec_OBS_2   = mean_2 + var_2*randn(1,nObs);

% --- perform permutation tests
if alpha < (1/nPerm)
   disp(['Not enough permutations for this significance level (' num2str(alpha) ')']);
   alpha = 1/nPerm;
   disp(['--> Significance level set to minimum possible (' num2str(alpha) ')']);
end
schuffle_index = uniqueShuffle2(nPerm, nObs, 1); % permutation matrix used to shuffle values between observation sets
diff_PERM = zeros(1,nPerm);     
% compute OBSERVATION (reference/real value) 
diff_OBS = mean(vec_OBS_2 - vec_OBS_1);  % here: average paired differences (CAN BE A SIMPLE DIFFERENCE, OR ANY CALCULATION!!!)
% compute PERMUTATION (surrogate values) 
vec_PERM_1 = zeros(1,nObs);
vec_PERM_2 = zeros(1,nObs);
for perm_ix = 1:nPerm
    % draw specific permutation using 'schuffle_index' logical indexes
	vec_PERM_1(schuffle_index(perm_ix,:))   = vec_OBS_1(schuffle_index(perm_ix,:));
    vec_PERM_1(~schuffle_index(perm_ix,:))  = vec_OBS_2(~schuffle_index(perm_ix,:));
    vec_PERM_2(~schuffle_index(perm_ix,:))  = vec_OBS_1(~schuffle_index(perm_ix,:));
    vec_PERM_2(schuffle_index(perm_ix,:))   = vec_OBS_2(schuffle_index(perm_ix,:));
    % compute surrogate surrogate value for this specific permutation
    diff_PERM(perm_ix) =  mean(vec_PERM_2 - vec_PERM_1);
end

% compute statistical significance and plot results
p_val = (1/nPerm)*(1+length(find(diff_PERM > diff_OBS)));
h_val = p_val < alpha;  
dispBootstrap(diff_OBS, diff_PERM)      
        

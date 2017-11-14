function [h_out,p_out] = permSignif(OBS,SURRO,ALPHA)
% function  [h_out,p_out] = permSignif(OBS,SURRO,ALPHA)
%
% ************************************************************************
% Get significance level from observation (real) values and permutation 
% (surrogate) values.
%
% TODO: ADD POSSIBILITY TO CHANGE THE EFFECT DIRECTION (OBS> or OBS< SURRO)
%
% Input:
% ------
% OBS       --> scalar or vector (1,nbTests) of observation values
% SURRO     --> vector (1,nPerms) or matrix (nbPerms,nbTests) of surrogate values  
% ALPHA   	--> (optional) significance level (must be >= 1/nPerms)
%               default: ALPHA=.05
%
% Output:
% -------
% h_out, p_out --> scalars (or vectors) with null hyp. rejection, and associated p-value
%
% Example:
% --------
% SURRO = randn(1000,10);
% OBS = [.25:.25:2.5];
% [h_out,p_out] = permSignif(OBS,SURRO),
%
% see also: dispBootstrap()
% 
% History:
% --------
% *** 2016-05-15
% Created by J.Chatel-Goldman, jonas.chatel.goldman(at)gmail.com
%



% -- assess inputs
if(size(OBS,1) == 1)
    OBS = OBS';
end
if(size(SURRO,1) == 1)
    SURRO = SURRO';
end
nbTests = size(OBS,1);
nbPerms  = size(SURRO,1);
if (nargin<3) || (isempty(ALPHA))
    ALPHA = .05;
end
if ALPHA < 1/nbPerms
    error('not enough permutations for this significance level.');
end

% -- Process significance level from permutation values
p_out 	= zeros(1,nbTests);
SURRO   = sort(SURRO,1,'ascend');
for test_ix = 1:nbTests
    p_out(test_ix) = (1/nbPerms)*length(find(SURRO(:,test_ix) > OBS(test_ix)));
    if(p_out(test_ix)==0)
        p_out(test_ix) = 1/nbPerms;     % case when no permutation gives lower p value than observation: put to minimum sensitivity
    end;
end    
h_out = (p_out <= ALPHA);



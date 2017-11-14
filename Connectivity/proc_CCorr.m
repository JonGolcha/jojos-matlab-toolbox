function [RHO,PVAL] = proc_CCorr(X)
% [RHO PVAL] = proc_CCorr(X)
%
% Wrapper on top of functions included in the Circular Statistics Toolbox 
% for Matlab. Allows to compute CCOR for all possible combinations in 2D
% data. Also takes 3D data with multipe trials in input.
%
% Inputs:
% - X   --> 2D matrix with the data (P channels by N samples)
%                       OR 3D matrix (P channels, N samples, T trials).
%
% Outputs:
% - RHO 	--> correlation coefficients
% - PVAL    --> p-values
%
%
% History
% Last version:  07/10/2016
% Created by J.Chatel-Goldman, jonas.chatel.goldman(at)gmail.com
% WITH CODE FROM Circular Statistics Toolbox for Matlab
% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

nb_chans = size(X,1);
if ismatrix(X)
    RHO = zeros(nb_chans,nb_chans);
    PVAL= zeros(nb_chans,nb_chans);
    for i = 1:nb_chans-1
        for j = i+1:nb_chans
            [RHO(i,j), PVAL(i,j)] = circ_corrcc_internal(X(i,:), X(j,:));
        end
    end
        
elseif ndims(X)==3
    nb_trials = size(X,3);
    RHO = zeros(nb_chans,nb_chans,nb_trials);
    PVAL= zeros(nb_chans,nb_chans,nb_trials);
    for t = 1:nb_trials
        for i = 1:nb_chans-1
            for j = i+1:nb_chans
                [RHO(i,j,t), PVAL(i,j,t)] = circ_corrcc_internal(squeeze(X(i,:,t)), squeeze(X(j,:,t)));
            end
        end
    end
end
end


function [rho,pval] = circ_corrcc_internal(alpha1, alpha2)
    % compute mean directions
    n = length(alpha1);
    alpha1_bar = circ_mean_internal(alpha1);
    alpha2_bar = circ_mean_internal(alpha2);

    % compute correlation coeffcient from p. 176
    num = sum(sin(alpha1 - alpha1_bar) .* sin(alpha2 - alpha2_bar));
    den = sqrt(sum(sin(alpha1 - alpha1_bar).^2) .* sum(sin(alpha2 - alpha2_bar).^2));
    rho = num / den;	

    % compute pvalue
    l20 = mean(sin(alpha1 - alpha1_bar).^2);
    l02 = mean(sin(alpha2 - alpha2_bar).^2);
    l22 = mean((sin(alpha1 - alpha1_bar).^2) .* (sin(alpha2 - alpha2_bar).^2));

    ts = sqrt((n * l20 * l02)/l22) * rho;
    pval = 2 * (1 - normcdf(abs(ts)));
end


function [mu] = circ_mean_internal(alpha, w)
    % check vector size
    if size(alpha,2) > size(alpha,1)
        alpha = alpha';
    end

    if nargin<2
      % if no specific weighting has been specified
      % assume no binning has taken place
        w = ones(size(alpha));
    else
      if size(w,2) > size(w,1)
        w = w';
      end 
    end

    % compute weighted sum of cos and sin of angles
    r = w'*exp(1i*alpha);

    % obtain mean by
    mu = angle(r);
end

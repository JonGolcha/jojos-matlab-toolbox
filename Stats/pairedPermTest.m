function [h_out,p_out] = pairedPermTest(x,y,alpha,nbPerms)
% function  [h_out,p_out] = pairedPermTest(x,y,alpha,nbPerms)
%
% ************************************************************************
% Perform paired-sample statistical testing using randomization procedure.
% Null hypothesis: differences x(i)-y(i) drawn from distribution with mean 
% O and unknown variance.
%
% Input:
% ------
% x,y           --> vectors of paired observations with same lengths (N,1)
%                   can also be matrices of size (N,M), in which case M 
%                   separate tests are performed along each column of x,y 
%                   and a vector of results is returned.
% alpha     	--> (optional) significance level, (must be >= 2^-N)
% nbPerms     	--> (optional) performs either systematic or random data
%                   permutation:
%                   'all': compute all permutations (2^N), or
%                	scalar: compute a random subset only, must be: (1/alpha
%                	< scalar < 2^N)
%
% Output:
% -------
% h_out, p_out --> scalars with null hyp. rejection, and associated p-value
%
% History:
% --------
% *** 2013-06-12
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com
%
% ADDITIONAL EXPLANATION NOTE (from EEGlab wiki!)
% We do not compute the mean condition difference, but its t-value (the 
% mean difference divided by the standard deviation of the difference and 
% multiplied by the square root of the number of observations less one).
% The result is equivalent to using the mean difference. The advantage is
% that when we have more conditions, we can use the comparable ANOVA
% measure. Computing the probability density distribution of the t-test or
% ANOVA is only a "trick" to be able to obtain a difference measure across
% all subjects and conditions. It has nothing to do with relying on a
% parametric t-test or ANOVA model, which assume underlying gaussian value
% distributions. Note that only measures that are well-behaved (e.g., are
% not strongly non-linearly related) should be processed by this kind of
% non-parametric testing.


% -- assess inputs
if(size(x,1) == 1)
    x = x';
    y = y';
end
N       = size(x,1);
nbTests = size(x,2);
if (nargin<3) || (isempty(alpha))
    alpha = .05;
end
if (nargin<4) || (isempty(nbPerms))
    nbPerms = 'all';
end
if length(x(:)) ~= length(y(:))
   error('input data x and y must have equal length.');
end
if alpha < 2^-N 
    error('not enough possible permutations for this significance level.');
end

if ~strcmp(nbPerms,'all') && ((nbPerms<(1/alpha))||(nbPerms>2^N))
    error('invalid requested number of permutations.');
end


% -- create shuffling matrix
if strcmp(nbPerms,'all')
    nbPerms = 2^N;
end
permMat = uniqueShuffle2(nbPerms,N,1);  % permutation matrix used to shuffle observations


% -- Perform paired t-test for OBSERVATION SET
[h_obs,p_obs,ci,stats] = ttest(x,y);
tval_obs = stats.tstat;

% -- Perform paired t-test for PERMUTATION SET
perm_set1   = zeros(N,nbTests);
perm_set2   = zeros(N,nbTests);
tval_perm	= zeros(nbPerms,nbTests);
for perm_ix = 1:nbPerms
    % user feedback
    if ~mod(100*perm_ix/nbPerms,10)
        disp(['paired-sample statistical testing... ' int2str(100*perm_ix/nbPerms) '%']);
    end
    
    % draw specific permutation using 'permMat' logical indexes
    if(nbTests == 1)
        perm_set1(permMat(perm_ix,:))  = x(permMat(perm_ix,:));
        perm_set1(~permMat(perm_ix,:)) = y(~permMat(perm_ix,:));
        perm_set2(~permMat(perm_ix,:)) = x(~permMat(perm_ix,:));
        perm_set2(permMat(perm_ix,:))  = y(permMat(perm_ix,:));
    else
        perm_set1(permMat(perm_ix,:),:)  = x(permMat(perm_ix,:),:);
        perm_set1(~permMat(perm_ix,:),:) = y(~permMat(perm_ix,:),:);
        perm_set2(~permMat(perm_ix,:),:) = x(~permMat(perm_ix,:),:);
        perm_set2(permMat(perm_ix,:),:)  = y(permMat(perm_ix,:),:);
    end
    % compute t-test on this permutation
    [h_temp,p_temp,ci,stats] = ttest(perm_set1,perm_set2);
    tval_perm(perm_ix,:) = stats.tstat;
end

% -- Process significance level from permutation values
p_out       = zeros(1,nbTests);
tval_perm   = sort(tval_perm,1,'ascend');
for test_ix = 1:nbTests
    p_out(test_ix) = (1/nbPerms)*length(find(tval_perm(:,test_ix) < tval_obs(test_ix)));
    if(p_out(test_ix)==0)
        p_out(test_ix) = 1/nbPerms;     % case when no permutation gives lower p value than observation: put to minimum sensitivity
    end;
end    
h_out = (p_out <= alpha);





function SCHUFFLE_INDEX = uniqueShuffle2(N_SHUFFLE, NB_WIN, BINARY_OUTPUT)
    % function  SCHUFFLE_INDEX = uniqueShuffle2(N_SHUFFLE, NB_WIN, BINARY_OUTPUT)
    %
    % ************************************************************************
    % This function return a random permutation matrix of size 
    % (N_SHUFFLE,NB_WIN) with only unique rows filled with random indexes 
    % in range [1->NB_WIN].
    % Main purpose is to compute statistical permutation tests between two 
    % vectors of identical length NB_WIN.
    %
    % Input:
    % ------
    % N_SHUFFLE     --> number of shuffling requested
    % NB_WIN        --> scalar with length of vectors to shuffle 
    % BINARY_OUTPUT --> (optional) if set to 1, outputs a binary vector of shuffling indexes with logical values
    %
    % Output:
    % -------
    % SCHUFFLE_INDEX    --> matrix of shuffling indexes with decimal values (N_SHUFFLE,NB_WIN)
    %
    % History:
    % --------
    % *** 2012-01-14
    % Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com


    if (nargin<3)||isempty(BINARY_OUTPUT)
       BINARY_OUTPUT = 0; 
    end

    if(N_SHUFFLE > (2^NB_WIN))
        error(['Requested number of permutations exceed max possible perm with ' int2str(NB_WIN) ' time windows']);
    % handle the particular case where all possible different permutations must be found
    elseif (N_SHUFFLE == (2^NB_WIN))
        SCHUFFLE_INDEX = dec2bin(0:(2^NB_WIN)-1) == '1'; 
    else 
        SCHUFFLE_INDEX      = zeros(N_SHUFFLE,NB_WIN);
        SCHUFFLE_INDEX(1,:) = randperm(NB_WIN); 
        for p_ix = 2:N_SHUFFLE
        %         disp(['progress: ' int2str(100*p_ix/NB_WIN) '%'])
            allUnique = false;
            while(~allUnique) 
                % create random vector with unique values
                SCHUFFLE_INDEX(p_ix,:) = randperm(NB_WIN); 

                % reiterate until this vector is unique within SCHUFFLE_INDEX
                sameShuff_prev = false(1,NB_WIN);
                sameShuff_prev(1:p_ix-1) = true;
                for win_ix = 1:NB_WIN
                    sameShuff_curr = (SCHUFFLE_INDEX(sameShuff_prev,win_ix) == SCHUFFLE_INDEX(p_ix,win_ix));
                    if ~any(sameShuff_curr)
                        allUnique = true;
                        break;
                    end
                    sameShuff_prev = sameShuff_curr;
                end
            end
        end
    end

    if BINARY_OUTPUT
        SCHUFFLE_INDEX = SCHUFFLE_INDEX > NB_WIN/2; % set half of random shuffling values to zero, and half to one.
    end
end


end



%% TRASH
% trying to consider NaN values
% if find(isnan(x)) || find(isnan(y))
%     try 
%         [h_obs,p_obs,ci,stats] = ttest2(x,y);
%     catch
%         error('data has NaN values, therefore we need the statistics and machine learning toolbox!');
%     end
% else
%     [h_obs,p_obs,ci,stats] = ttest(x,y);
% end
% tval_obs = stats.tstat;
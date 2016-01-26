function [h_sep,p_sep,h_all,p_all] = pairedPermTest2(X,alpha,nbPerms)
% function  [h_sep,p_sep,h_all,p_all] = pairedPermTest2(X,alpha,nbPerms)
%
% ************************************************************************
% Perform paired-sample statistical testing using randomization procedure.
% Considers G groups X{i} of dependent (paired) observations, that can
% be randomized within groups but not across groups. 
% Each group X{i} may correspond e.g., to a subject dyad with measures 
% X{i}(:,1) and X{i}(:,2) obtained e.g., in two conditions to be 
% contrasted, possibly with different number of observations.
%
% We test two Null hypothesis:
% 1) differences within groups X{i}(:,1)-X{i}(:,2) drawn from distribution 
% with mean O and unknown variance.
% 2) global difference X{1:G}(:,1)-X{1:G}(:,2) drawn from distribution with 
%   mean O and unknown variance.
%
% NOTE: This is an upgraded version of pairedPermTest.m:
% - handles non equal number of paired observation within groups.
% - performs global computation on multiple groups of paired observations.
%
% Input:
% ------
% X             --> cells with 1 to G groups of observation
%                   each group is a matrix of size(N,2,nb_tests) 
%                   IMPORTANT NOTE 1: each group may have different numbers 
%                   of observations, in which case smallest vector MUST BE 
%                   FILLED WITH NaNs and "Statistics and Machine Learning 
%                   Toolbox" MUST BE INSTALLED. 
%                   IMPORTANT NOTE 1: multiple separate tests can be 
%                   performed at once (nb_tests >1) and a vector of results is returned. 
%                   within groups, same permutation is performed for all separate tests 
% alpha     	--> (optional) significance level, (must be >= 2^-N)
% nbPerms     	--> (optional) performs either systematic or random X permutation
%                   'all': compute all permutations (2^N), or
%                	scalar: compute a random subset only, must be: (1/alpha < scalar < 2^N) 
%
% Output:
% -------
% h_sep, p_sep --> vectors with null hyp. rejection, and associated p-value (one scalar for each group and separate test)
% h_all, p_all --> scalars with null hyp. rejection, and associated p-value (all groups taken into consideration)
%
% History:
% --------
% *** 2015-11-03
% J.Chatel-Goldman, jonas.chatel.goldman(at)gmail.com


% -- assess inputs
if ~iscell(X)
    error('input must be a cell')
end
nb_groups = size(X,2);
if (nargin<3) || (isempty(alpha))
    alpha = .05;
end
if (nargin<3) || (isempty(nbPerms))
    nbPerms = 'all';
end
N = zeros(2,nb_groups); % sample size for each group
nb_tests = size(X{1},3);
sameSampleSize = true;
for group_ix = 1:nb_groups
    N(1,group_ix) = length(find(~isnan(X{group_ix}(:,1))));
    N(2,group_ix) = length(find(~isnan(X{group_ix}(:,2))));
    if alpha < 2^-(min(N(1,group_ix),N(2,group_ix)))    % not cool, but still feasable
        warning(['min sample size is too small for this significance level (not enough permutation are possible), group #' int2str(group_ix)]);
    end
    if alpha < 2^-(max(N(1,group_ix),N(2,group_ix)))    % not feasable
        error(['max sample size is too small for this significance level (not enough permutation are possible), group #' int2str(group_ix)]);
    end
    if N(1,group_ix) ~= N(2,group_ix)
        sameSampleSize = false;
    end
end
N_max = max(max(N)); % this is the maximum sample size across all groups
N_all = sum(max(N));
if ~strcmp(nbPerms,'all') && ((nbPerms<(1/alpha))||(nbPerms>2^N_max))
    error('invalid requested number of permutations.');
end


% -- Create shuffling matrix for each single group and for all groups together
permMat = cell(1,nb_groups);
for group_ix = 1:nb_groups
    curr_N = max(N(1,group_ix),N(2,group_ix));
    if strcmp(nbPerms,'all')    % if not specified as function input
        nbPerms = 2^curr_N;
    end
    permMat{group_ix} = uniqueShuffle2(nbPerms,curr_N,1);  % permutation matrix used to shuffle observations
end
permMat_all = uniqueShuffle2(nbPerms,N_all,1);  % permutation matrix used to shuffle observations


% -- Perform paired t-test for OBSERVATION SET (each group and global)
% the result h_obs is 1 if the test rejects the null hypothesis at the 5% 
% significance level, and 0 otherwise.
h_obs_sep = zeros(nb_groups,nb_tests);
t_obs_sep = zeros(nb_groups,nb_tests);
x_all = [];
y_all = [];
for group_ix = 1:nb_groups
    x = squeeze(X{group_ix}(:,1,:));
    y = squeeze(X{group_ix}(:,2,:));
    x_all = cat(1,x_all,x);
    y_all = cat(1,y_all,y);
	if sameSampleSize 
        [h_obs_sep(group_ix,:),p,ci,stats] = ttest(x,y); 
     	t_obs_sep(group_ix,:) = stats.tstat ;
    else % different sample sizes
        try
            [h_obs_sep(group_ix,:),p,ci,stats] = ttest2(x,y); 
            t_obs_sep(group_ix,:) = stats.tstat ;
        catch
            error('Statistics and Machine Learning Toolbox must be installed to perform paired test with different sample sizes.')
        end
    end
end
if sameSampleSize
    [h_obs_all,p,ci,stats] = ttest(x_all,y_all); 
    t_obs_all = stats.tstat ;
else
    [h_obs_all,p,ci,stats] = ttest2(x_all,y_all); 
    t_obs_all = stats.tstat ;
end


% -- Perform paired t-test for PERMUTATION SETS (each group and global)
% first each group separate
t_perm_sep = zeros(nbPerms,nb_groups,nb_tests);
for group_ix = 1:nb_groups
    curr_N      = max(N(1,group_ix),N(2,group_ix));
    perm_set1   = zeros(curr_N,nb_tests);
    perm_set2   = zeros(curr_N,nb_tests);
    x = squeeze(X{group_ix}(:,1,:));
    y = squeeze(X{group_ix}(:,2,:));
    for perm_ix = 1:nbPerms
        % draw specific permutation using 'permMat' logical indexes
        perm_set1(permMat{group_ix}(perm_ix,:),:)  = x(permMat{group_ix}(perm_ix,:),:);
        perm_set1(~permMat{group_ix}(perm_ix,:),:) = y(~permMat{group_ix}(perm_ix,:),:);
        perm_set2(~permMat{group_ix}(perm_ix,:),:) = x(~permMat{group_ix}(perm_ix,:),:);
        perm_set2(permMat{group_ix}(perm_ix,:),:)  = y(permMat{group_ix}(perm_ix,:),:);
        % compute t-test on this specific permutation
        if sameSampleSize
            [h_temp,p,ci,stats] = ttest(perm_set1,perm_set2); 
            t_perm_sep(perm_ix,group_ix,:) = stats.tstat ;
        else
            [h_temp,p,ci,stats] = ttest2(perm_set1,perm_set2); 
            t_perm_sep(perm_ix,group_ix,:) = stats.tstat ;
        end 
    end
end

% compute permutation on all groups together
t_perm_all = zeros(nbPerms,nb_tests);
perm_set1   = zeros(N_all,nb_tests);
perm_set2   = zeros(N_all,nb_tests);
for perm_ix = 1:nbPerms
    % draw specific permutation using 'permMat' logical indexes
    perm_set1(permMat_all(perm_ix,:),:)  = x_all(permMat_all(perm_ix,:),:);
    perm_set1(~permMat_all(perm_ix,:),:) = y_all(~permMat_all(perm_ix,:),:);
    perm_set2(~permMat_all(perm_ix,:),:) = x_all(~permMat_all(perm_ix,:),:);
    perm_set2(permMat_all(perm_ix,:),:)  = y_all(permMat_all(perm_ix,:),:);
    % compute t-test on this specific permutation
    if sameSampleSize
            [h_temp,p,ci,stats] = ttest(perm_set1,perm_set2); 
            t_perm_all(perm_ix,:) = stats.tstat ;
    else
        [h_temp,p,ci,stats] = ttest2(perm_set1,perm_set2); 
     	 t_perm_all(perm_ix,:) = stats.tstat ;
    end 
end


% -- Process significance level from permutation values
p_sep   = zeros(nb_groups,nb_tests);
h_sep   = zeros(nb_groups,nb_tests);
t_perm_sep  = sort(t_perm_sep,1,'ascend');
t_perm_all  = sort(t_perm_all,1,'ascend');
for group_ix = 1:nb_groups
    p_sep(group_ix,:) = (1/nbPerms)*sum((squeeze(t_perm_sep(:,group_ix,:)) < repmat(t_obs_sep(group_ix,:),[nbPerms 1])),1);
    if(p_sep(group_ix,:)==0)
        p_sep(group_ix,:) = 1/nbPerms;     % case when no permutation gives lower t value than observation
    end;
    h_sep(group_ix) = (p_sep(group_ix) <= alpha);
end    
p_all = (1/nbPerms)*sum(t_perm_all < repmat(t_obs_all,[nbPerms 1]),1);
if p_all==0
    p_all = 1/nbPerms;     % case when no permutation gives lower t value than observation
end;
h_all = (p_all <= alpha);


end


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
        SCHUFFLE_INDEX = dec2bin(0:(2^NB_WIN)-1) == '1'; % LOGICAL ONLY (DECIMAL CASE FOR ALL PERMUTATIONS IS NOT TREATED YET)
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
        if BINARY_OUTPUT
            SCHUFFLE_INDEX = SCHUFFLE_INDEX > NB_WIN/2; % set half of random shuffling values to zero, and half to one.
        end
    end
end


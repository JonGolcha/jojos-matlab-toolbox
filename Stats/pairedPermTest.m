function [h_out,p_out] = pairedPermTest(x,y,alpha,nbPerms)
% function  [h,p] = pairedPermTest(x,y,alpha,nbPerms)
%
% ************************************************************************
% Perform paired-sample statistical testing using randomization procedure.
% Null hypothesis: differences x(i)-y(i) drawn from distribution with mean 
% O and unknown variance.
%
% Input:
% ------
% x,y           --> vectors of paired observations with same lengths (1,N)
%                   can also be matrices of size (M,N), in which case M 
%                   separate tests are performed along each column of x,y 
%                   and a vector of results is returned.
% alpha     	--> (optional) significance level, (must be >= 2^-N)
% nbPerms     	--> (optional) performs either systematic or random data permutation
%                   'all': compute all permutations (2^N), or
%                	scalar: compute a random subset only, must be: (1/alpha < scalar < 2^N) 
%
% Output:
% -------
% h_out, p_out --> scalars with null hyp. rejection, and associated p-value
%
% History:
% --------
% *** 2013-06-12
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com


% -- assess inputs
if(size(x,2) == 1)
    x = x';
    y = y';
end
N       = size(x,2);
nbTests = size(x,1);
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
[h_obs,p_obs] = ttest(x,y);

% -- Perform paired t-test for PERMUTATION SET
perm_set1   = zeros(N,nbTests);
perm_set2   = zeros(N,nbTests);
p_perm      = zeros(nbPerms,nbTests);
for perm_ix = 1:nbPerms
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
    [h_temp,p_perm(perm_ix,:)] = ttest(perm_set1,perm_set2);
end

% -- Process significance level from permutation values
p_out   = zeros(1,nbTests);
p_perm  = sort(p_perm,1,'ascend');
for test_ix = 1:nbTests
    p_out(test_ix) = (1/nbPerms)*length(find(p_perm(:,test_ix) < p_obs(test_ix)));
    if(p_out(test_ix)==0)
        p_out(test_ix) = 1/nbPerms;     % case when no permutation gives lower p value than observation
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
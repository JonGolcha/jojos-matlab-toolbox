function [S1,S2] = uniqueShuffle2(N_SHUFFLE, NB_WIN, BINARY_OUTPUT)
% function  [S1,S2] = uniqueShuffle2(N_SHUFFLE, NB_WIN, BINARY_OUTPUT)
%
% ************************************************************************
% This function return two random permutation matrices of size 
% (N_SHUFFLE,NB_WIN) with only unique rows filled with random indexes 
% in range [1->NB_WIN] and no identical indexes for each corresponding row.
% e.g., case S1 = [ 1 2 3 ; 1 3 2] and S2 = [ 3 2 1 ; 2 1 3] is forbidden
% due to value 2 in first row.
%
% Main purpose is to create surrogate data from two vectors of identical 
% length NB_WIN. These vectors can be for example trial indexes in some
% experiment.
%
% 
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
% *** 2012-01-14, update on 2016-05-14 to add a second output index
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com



if (nargin<3)||isempty(BINARY_OUTPUT)
   BINARY_OUTPUT = 0; 
end

% create the first shuffle index
if(N_SHUFFLE > (2^NB_WIN))
    error(['Requested number of permutations exceed max possible perm with ' int2str(NB_WIN) ' time windows']);
% handle the particular case where all possible different permutations must be found
elseif (N_SHUFFLE == (2^NB_WIN))
    S1 = dec2bin(0:(2^NB_WIN)-1);   
else 
    S1      = zeros(N_SHUFFLE,NB_WIN);
    S1(1,:) = randperm(NB_WIN); 
    for p_ix = 2:N_SHUFFLE
    %         disp(['progress: ' int2str(100*p_ix/NB_WIN) '%'])
        allUnique = false;
        while(~allUnique) 
            % create random vector with unique values
            S1(p_ix,:) = randperm(NB_WIN); 

            % reiterate until this vector is unique within S1
            sameShuff_prev = false(1,NB_WIN);
            sameShuff_prev(1:p_ix-1) = true;
            for win_ix = 1:NB_WIN
                sameShuff_curr = (S1(sameShuff_prev,win_ix) == S1(p_ix,win_ix));
                if ~any(sameShuff_curr)
                    allUnique = true;
                    break;
                end
                sameShuff_prev = sameShuff_curr;
            end
        end
    end
end

% create the second shuffle index, for which each corresponding row must
% have different corresponding indexes, 
allUnique = false;
while(~allUnique) 
    % create vector of N_SHUFFLE random index offsets
    offset = randi(NB_WIN,[N_SHUFFLE 1]);
    % add it to S1 to create S2
    S2 = mod(S1+repmat(offset,[1 NB_WIN]),NB_WIN)+1;
    % check that all rows are unique (extremly low probability, but still...)
    C = unique(S2,'rows');
    if(size(C,1)==N_SHUFFLE)
        allUnique = true;
    end
end

if BINARY_OUTPUT
    S1 = S1 > NB_WIN/2; % set half of random shuffling values to zero, and half to one.
    S2 = S2 > NB_WIN/2; % set half of random shuffling values to zero, and half to one.
end

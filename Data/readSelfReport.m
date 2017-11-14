function ITEMS = readSelfReport(FILENAME, NB_ITEMS)
% function ITEMS = readSelfReport(FILENAME, NB_ITEMS)
% -------------------------------------------------------------------------
% Analysis of dual-EEG dual-empathy experiment data (Gispa-Lab 2012)
% -------------------------------------------------------------------------
%
% This function reads subject answers to self-report questionnaires. 
%
% Inputs:
% - FILENAME  	--> Name of .txt file with answers, with datapath included.
% - NB_ITEMS    --> (optional) Number of items for each question
%
% Outputs:
% - ITEMS    	--> Matrix with subject answers (N_trials,N_questions)
%
% History:
% --- 2013-05-07
% Created by Jonas Chatel-Goldman @ GIPSA-Lab 


% Check input file .txt extension
[pathstr, name, ext] = fileparts(FILENAME);
if ~strcmpi(ext,'.txt')
    error('[Function ''readSelfReport'']: input file should be a simple ''.txt'' file');
end
if (nargin < 2) || isempty(NB_ITEMS)
    NB_ITEMS = 1;
end

%  Open and read text file
[fid, message] = fopen(FILENAME);
if fid == -1
   error(['[readSelfReport] ' FILENAME ' -> ' message]), 
end
C = textscan(fid, '%d %d','HeaderLines',1);    	% save item index C{1} and subject answer  C{2}
fclose(fid);    % close file

% Put answer to each question type in a specific column of output ITEMS
questions   = C{1,1};
answers     = C{1,2};


if (NB_ITEMS == 1)
    ITEMS = answers;
else
    try
        ITEMS = zeros(length(questions)/NB_ITEMS,NB_ITEMS);
        for item_ix = 1:NB_ITEMS
            ITEMS(:,item_ix)  = answers(questions==item_ix);
        end
    catch
        error('[readSelfReport] Did you specify the right number of items?...'); 
    end
end


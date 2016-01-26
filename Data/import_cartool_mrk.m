function MASK_MRK = import_cartool_mrk(MRK_FILENAME, DATA_LENGTH)

% -------------------------------------------------------------------------
% Analysis of dual-EEG dual-empathy experiment data (Gispa-Lab 2012)
% -------------------------------------------------------------------------
%
% This function imports a Cartool .mrk marker file into a logical mask. 
%
% Inputs:
% - MRK_FILENAME    --> File name of Cartool marker file, with datapath 
%                       included (file extension must be '.mrk')
% - DATA_LENGTH     --> Total number of samples 
%
% Outputs:
% - MASK_MRK        --> Mask: binary line vector of size 1 x DATA_LENGTH
%                       '0' -> noise    ;   '1' -> artifact free
%
% History:
% --- 2012-07-09
% Created by Jonas Chatel-Goldman @ GIPSA-Lab 


%% Import mrk file
% Check input file .mrk extension
[pathstr, name, ext] = fileparts(MRK_FILENAME);
if ~strcmpi(ext,'.mrk')
    error('[Function ''import_cartool_mrk'']: input file should be a cartool ''.mrk'' file');
end

%  Open and read text file
fid = fopen(MRK_FILENAME);
C = textscan(fid, '%*s', 1);      	% read first line (and do nothing)
C = textscan(fid, '%u32 %u32 %s'); 	% save markers begin and end times
% C{1} are mrk begin times ; C{2} are mrk end times ;  C{3} are mrk labels
fclose(fid);

% Create mask from marker times (mask = 0 within marker artifacts periods)
nb_mrk = length(C{1}); % number of markers
MASK_MRK = ones(1,DATA_LENGTH);
for mrk_ix = 1:nb_mrk  
    if C{1}(mrk_ix)==0  C{1}(mrk_ix)=1;     end
    if C{2}(mrk_ix)==0  C{2}(mrk_ix)=1;     end
    
    if(strcmp(cell2mat(C{3}(mrk_ix)),'"Artifact"'))
        MASK_MRK(C{1}(mrk_ix):C{2}(mrk_ix)) = 0; 
    end
end
MASK_MRK = logical(MASK_MRK);



function saveCellToFile(CELL, Filename, Path)
% Simple function to save a cell array of strings in a path with a given filename
% Output format is ascii text file.
% Useful when having different number of characters in each cell part
% Created by Jonas Chatel-Goldman @ GIPSA-Lab 


if(nargin < 3)
    Filepath = Filename;            % If a path is not given, saves on the current directory
else
    Filepath = [Path Filename];
    if ~exist(Path,'dir')
    	mkdir(Path);
    end
end
isUI = false;

total_lines = size(CELL, 1);        % Total of lines in the file
[fileId errorMsg] = fopen(Filepath, 'w+t');                 % Open/create the file to write ('w+') in text mode ('t')
if(fileId == -1)
   disp(errorMsg);
end
    
for line = 1 : total_lines   	% For each lines
    fprintf(fileId, '%s', cell2mat(CELL(line,:)));        % Write the line
    fprintf(fileId, '\n');                      % Line feed    
    
    if isUI && ~mod(line,500) 
        disp(['Writing file ' Filename ': ' int2str(100*line/total_lines) '%']); 
    end;
end          
fclose(fileId);                                 % Close the file
function [X] = vec2vec(VEC1,VEC2, SPACING)
% ************************************************************************
% creates a 2-D matrix X(i,j) with all values from vec1(i) to vec2(i).
% defaults linear spacing can be set by user. 
% /!\ it is assumed here that :
% 1) input vectors have values with equal spacing
% 2) vec2(i) > vec1(i) for all i
%
% INPUTS
% - VEC1       --> vector of size i
% - VEC2       --> vector of size i
% - SPACING    --> (optionnal) spacing for creating vectors (default is 1)
%
% OUTPUTS
% - X       --> 2-D matrix os size (i,j)
%
% HISTORY 
% First version: 18-11-2013
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

 
if nargin<3 || isempty(SPACING)
    SPACING = 1;
end
if length(VEC1) ~= length(VEC2)
    error('[vec2vec] input vectors must have the same length!'); 
end

i = length(VEC1);
j = length(VEC1(1):SPACING:VEC2(1));
X = zeros(i,j);
try
    for ix = 1:i
        X(ix,:) = VEC1(ix):SPACING:VEC2(ix);
    end
catch
   error('[vec2vec] input vectors have values with unequal spacing!'); 
end

end


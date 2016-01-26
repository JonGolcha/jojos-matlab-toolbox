function [Y] = concatDim3(X,DIM)
% ************************************************************************
% Concatenates 3rd dimension on the rows or on the colums
%
% INPUTS
% - X       --> 3-D matrix 
% - DIM     --> scalar with dimension index: {'1'} for rows, '2' for columns
%
% OUTPUTS
% - Y       --> 2-D concatenated matrix 
% 
% EXAMPLE
%  a = [1 2 3 ; 4 5 6]; a(:,:,2) = [7 8 9 ; 10 11 12],
%  b = concatDim3(a,1), c = concatDim3(a,2),
% reverse operation: a2 = reshape(b',[s1 s2 s3]),
%
% HISTORY 
% First version: 18-11-2013
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com

if nargin<2
    DIM = 1;
end

if (DIM == 1)
    Y = reshape(permute(X,[2 1 3]),size(X,2),[])';
elseif (DIM == 2)
    Y = reshape(X,size(X,1),[]);
end

end


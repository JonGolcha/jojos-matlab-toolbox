function [Y]= normalizeMat(X)
% This function outputs Y, a normalized version of X (must be square) along 
% rows and columns, i.e., Y has unitary diagonal.
%
%
% Input     -->     X (N * N)
% Output    -->     Y (N * N)
%
% HISTORY 
% last modified 2014-07-18
% J.Chatel-Goldman @OIST, jonas.chatel.goldman(at)gmail.com


Y = diag(diag(X.^(-1/2))) * X * diag(diag(X.^(-1/2)));


end


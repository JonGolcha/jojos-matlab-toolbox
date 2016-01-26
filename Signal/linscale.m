function [Y] = linscale(X,MIN,MAX)
% ************************************************************************
% map X linearly so that its minimum value is MIN and its maximum value is
% MAX.
%
% INPUTS
% - X       --> vector of size N
% - MIN   	--> scalar
% - MAX     --> scalar
%
% OUTPUTS
% - Y       --> vector of size N
%
% HISTORY 
% First version: 04-09-2014
% Created by J.Chatel-Goldman @OIST, jonas.chatel.goldman(at)gmail.com

 


X_delta = max(X)-min(X);
Y_delta = MAX - MIN;

Y = (X-min(X))*(Y_delta/X_delta)+MIN;

end


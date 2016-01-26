function [] = plotSigDev(X, X_devs, col_X, col_X_dev, doFill)
% function [] = plotSigDev(X, X_devs, col_X, col_X_dev, doFill)
% *************************************************************************
% Display a set of signals (e.g. means) and their deviation (e.,g.,
% percentiles) in the form of filled region or dashed lines.
%
% INPUTS
% - X           --> matrix of main signals to plot (samples, variables).
% - X_devs      --> (optional) matrix with deviation values (samples, variables, 2)
%                   In the last dimension, hold respectively inferior and
%                   superior vectors.
% - col_X       --> (optional) RGB matrix with colors for main signals (variables x 3).
% - col_X_dev	--> (optional) RGB matrix with colors for deviation signals (variables x 3).
% - doFill  	--> (optional) show deviation as filled area '1' or dashed line '0'.  
%
%
% HISTORY 
% First version: Jonas Chatel-Goldman @ OIST, 13/11/2014


% assess inputs
% TODO: correct for the case where input in just one signal 
% TODO: add the dashed line case


[T, N] = size(X);
showDev = true;
if nargin<2 || isempty(X_devs)
    showDev = false;
end
if nargin<3 || isempty(col_X)
    col_X = lines(N);
end
if nargin<4 || isempty(col_X_dev)
    col_X_dev = col_X;
end
if nargin<5 || isempty(doFill)
    doFill = true;
end

% draw main signal
set(gca,'NextPlot','replacechildren'),
set(gca,'ColorOrder',col_X_dev),
plot(X,'LineWidth',2),


% draw deviation
if showDev
    for sig_ix = 1:N
        upper = squeeze(X_devs(:,sig_ix,2))';
        lower = squeeze(X_devs(:,sig_ix,1))';
        color = col_X_dev(sig_ix,:); % color of the filled area
        edge = clamp(col_X_dev(sig_ix,:)-50/255,[0 1]); % color around the edge of the filled area
        transparency = .1;
        filled=[upper,fliplr(lower)];
        xpoints=[(1:T),fliplr(1:T)];
        hold on,
        fillhandle=fill(xpoints',filled',color);%plot the data
        set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency,'EdgeAlpha',transparency);%set edge color
        hold off,
    end
end



function [y] = clamp(x, range)
    y = x;
    y(x<range(1)) = range(1);
    y(x>range(2)) = range(2);
end

end


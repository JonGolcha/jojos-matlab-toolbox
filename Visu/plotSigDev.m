function [h] = plotSigDev(X, Y, Y_devs, col_X, col_X_dev, doFill)
% function [h] = plotSigDev(X, Y, Y_devs, col_X, col_X_dev, doFill)
% *************************************************************************
% Display a set of signals (e.g. means) and their deviation (e.,g.,
% percentiles) in the form of filled region or dashed lines.
%
% INPUTS
% - X           --> (optional) vector to plot Y versus X (samples, 1).
% - Y           --> matrix of main signals to plot (samples, variables).
% - Y_devs      --> (optional) matrix with deviation values (samples, variables, 2)
%                   In the last dimension, hold respectively inferior and
%                   superior vectors.
% - col_X       --> (optional) RGB matrix with colors for main signals (variables x 3).
% - col_X_dev	--> (optional) RGB matrix with colors for deviation signals (variables x 3).
% - doFill  	--> (optional) show deviation as filled area '1' or dashed line '0'.  
%
% EXAMPLES
% % just one signal:
% Y = randn(100,1), Y_devs = [Y-1,Y+1], figure, plotSigDev([], Y, Y_devs,[],[],1),
% % multiple signals:
% X = -50:49; Y = randn(100,3)+repmat([-2 0 2],[100 1]), Y_devs = cat(3,Y-1,Y+1), figure, plotSigDev(X, Y, Y_devs,[],[],1),
% 
% HISTORY 
% First version: Jonas Chatel-Goldman @ OIST, 13/11/2014
% updated 17/06/2016


% assess inputs
if isrow(Y)
    Y=Y';
end
[T, N] = size(Y);
showDev = true;
if isempty(X)
    X = (1:T)';
end
if nargin<3 || isempty(Y_devs)
    showDev = false;
end
if nargin<4 || isempty(col_X)
    col_X = lines(N);
end
if nargin<5 || isempty(col_X_dev)
    col_X_dev = col_X;
end
if nargin<6 || isempty(doFill)
    doFill = true;
end

% draw main signal
set(gca,'NextPlot','replacechildren'),
set(gca,'ColorOrder',col_X_dev),
h = plot(X,Y,'LineWidth',1);


% draw deviation
if showDev
    for sig_ix = 1:N
        if N > 1
         	upper = squeeze(Y_devs(:,sig_ix,2))';
            lower = squeeze(Y_devs(:,sig_ix,1))';  
        else  % case only one signal to plot
            upper = Y_devs(:,2)';
            lower = Y_devs(:,1)';
        end
        color = col_X_dev(sig_ix,:); % color of the filled area
        if doFill && abs(sum(upper-lower))>10^-12 % bug when 'upper' is almost equal to 'lower'
            faceColor = color;
            transparency = .1;
            lineStyle = '-';
        else
            faceColor = 'none';
            transparency = 1;
            lineStyle = '--';
        end
        edge = clamp(col_X_dev(sig_ix,:)-50/255,[0 1]); % color around the edge of the filled area
        
        filled=[upper,fliplr(lower)];
        xpoints=[X',fliplr(X')];
        hold on,
        fillhandle=fill(xpoints',filled',color);%plot the data
        set(fillhandle,'FaceColor',faceColor,'EdgeColor',edge,'FaceAlpha',transparency,'EdgeAlpha',transparency,'LineStyle',lineStyle);%set edge color
        hold off,
    end
end



function [y] = clamp(x, range)
    y = x;
    y(x<range(1)) = range(1);
    y(x>range(2)) = range(2);
end

end


function bar_2sides(Y1, Y2, X, YLABELS, MAX_DISP)    
% function bar_2sides(Y1, Y2, X)      
% 
% This function plots corresponding horizontal bars on two sides of the y 
% axis ('tornado' mode). /!\ data must be positive only, or the separation
% between top and bottom display will not be possible.
%
% *** Inputs *** 
% - Y1  	--> draws one bar for each element in DATA_LEFT 
%                   (can be a 2d array for group bars)
% - Y2      --> draws one bar for each element in DATA_RIGHT 
%                   (can be a 2d array for group bars)
% - X    	--> (optional) values for x axis
% - YLABELS --> (optional) cell array with 2 strings 
% - MAX_DISP--> (optional) maximum display for Y1 and Y2 [max_Y1, max_Y2]
%
%  *** History *** 
% Last version: 09/08/2014
% J.Chatel-Goldman @OIST, jonas.chatel.goldman(at)gmail.com


%% assess inputs + various init
nb_tick_Y1 = 12;
nb_tick_Y2 = 20;

if nargin < 3
    X = 1:size(Y1,1);
end
if (nargin < 5) || (isempty(MAX_DISP))
    LIMS_Y1 = 1.2*[-max(max(Y1)), max(max(Y1))];
    LIMS_Y2 = 1.2*[-max(max(Y2)), max(max(Y2))];
else
  	LIMS_Y1 = [-MAX_DISP(1), MAX_DISP(1)];
    LIMS_Y2 = [-MAX_DISP(2), MAX_DISP(2)];
end

YTick1 = (0:(LIMS_Y1(2)-LIMS_Y1(1))/nb_tick_Y1:LIMS_Y1(2));
YTick2 = (0:(LIMS_Y2(2)-LIMS_Y2(1))/nb_tick_Y2:LIMS_Y2(2));

% This uses the "plotyy(x1,y1,x2,y2,fun1,fun2)" variant of plotyy
[ax,h1,h2] = plotyy(X,Y1,X,Y2, @(X,Y1) bar(X,Y1), @(X,Y2) bar(X,Y2));
set(ax(1),'Ylim',LIMS_Y1),
set(ax(2),'Ylim',LIMS_Y2),
set(ax(1),'Ylim',LIMS_Y1,'YTick',YTick1,'Ygrid','on'),
set(ax(2),'Ylim',LIMS_Y2,'YTick',YTick2,'Ygrid','on','YDir','reverse'),    
line(get(gca,'Xlim'),[0 0],'linewidth',5,'color','k'),

if nargin > 3
    set(get(ax(1),'Ylabel'),'String',['                                                                           ' YLABELS{1} ]),
    set(get(ax(2),'Ylabel'),'String',[ YLABELS{2} '                                                                           ']),
end


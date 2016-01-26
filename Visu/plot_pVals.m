function [h_plot] = plot_pVals(P_VALS,ALPHA,X_LABELS,LEGENDS,MIN_P,LOG_SCALE)
% function [h_plot] = plot_pVals(P_VALS,ALPHA,X_LABELS,LEGENDS,MIN_P,LOG_SCALE)
%
% ************************************************************************
% Display synthesis of multiple p-values, typically obtained from a 
% permutation test.
%
% Input:
% ------
% P_VALS        --> vector with p_values to display. 
%                   can be a matrix, in which case each column is displayed 
%                   with a separate abscissa.
% ALPHA     	--> (optional) significance level. Default: 0.05 
% X_LABELS     	--> (optional) cell with string 
% LEGENDS       --> (optional) cell with string 
% MIN_P         --> (optional) minimum p-value (in general 1/nPerms)
% LOG_SCALE     --> (optional) log-scale y axis (default: true) 
%
% Output:
% -------
% h_plot        --> plot handle
%
% History:
% --------
% *** 2015-11-14
% Created by J.Chatel-Goldman, jonas.chatel.goldman(at)gmail.com


% assess inputs
N = size(P_VALS,1);
if (nargin<2) || (isempty(ALPHA))
    ALPHA = .05;
end
if (nargin<3) || (isempty(X_LABELS))
    X_LABELS = [];
end
if (nargin<4) || (isempty(LEGENDS))
    LEGENDS = {strcat('obs #',int2str([1:N]'))};
end
if (nargin<5) || (isempty(MIN_P))
    MIN_P = 0;
end
if (nargin<6) || (isempty(LOG_SCALE))
    LOG_SCALE = true;
end


strSize     = 8;
markerSize  = 10;
markersVec  = {'o','s','d','p','x','*','+','v','<','>'};
x_val   = 1:size(P_VALS,2);
if LOG_SCALE
    y_val   = -log(P_VALS);
else
    y_val   = P_VALS;
end
colmap  = flipud(summer(6));
colMarker = [1 1 .9];
h_plot = plot(0,0,'linestyle','none');
hasbehavior(h_plot(1),'legend',false); 
set(gca,'ColorOrder',colmap(1:5,:));
hold on,
for i=1:N
    marker_ix = mod(i,10);
    if marker_ix == 0
        marker_ix = 1;
    end
    plot(x_val,y_val(i,:),'linestyle','none','marker',markersVec{marker_ix},...
        'MarkerFaceColor',colMarker,'MarkerSize',markerSize,'linewidth',1),
end
legend(LEGENDS,'FontSize',8,'Location','BestOutside'),
set(gca,'XTick',x_val,'XTickLabel',X_LABELS,'FontSize',strSize),
xlim([min(x_val)-.5 max(x_val)+.5])
% ylim([0 max(max(y_val))+1])    
ylim([0 8]) 
line(get(gca,'Xlim'),[-log(ALPHA) -log(ALPHA)],'color','r','linewidth',1,'linestyle','--'),
text(x_val(end)-.5,-log(ALPHA)+.1, ['p-val:' num2str(ALPHA)],'Color','r','FontSize',strSize);
line(get(gca,'Xlim'),[-log(0.01) -log(0.01)],'color',[255 102 0]./255,'linewidth',1,'linestyle','--'),
text(x_val(end)-.5,-log(.01)+.1, ['p-val: 0.01'],'color',[255 102 0]./255,'FontSize',strSize);
if MIN_P ~= 0
    line(get(gca,'Xlim'),[-log(MIN_P) -log(MIN_P)],'color','k','linewidth',1),
    text(x_val(end)-.5,-log(MIN_P)+.1, ['MIN p-val:' num2str(MIN_P)],'Color','k','FontSize',strSize);
end
if LOG_SCALE
    ylabel('p-value significance  ( -ln(p) )','FontSize',strSize),
else
    ylabel('p-value significance (p)','FontSize',strSize),
end










end


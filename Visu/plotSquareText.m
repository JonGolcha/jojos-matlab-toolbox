function plotSquareText(x,y,str,colEdge,colBgd)    
% function plotSquareText(x,y,str,col)    
% 
% This function given y=f(x) data using text annotation given in str.
% /!\ THIS FUNCTION DOES NOT SUPPORT ZOOMING IN/OUT AFTER!!!
%
% *** Inputs *** 
% - x       --> vector (1,N)
% - y       --> vector (1,N)
% - str   	--> (optional) cell of strings (1,N)
% - colEdge     --> (optional) color edge, [R,G,B] in range 0-255
% - colBgd     --> (optional) color background, [R,G,B] in range 0-255
%
%  *** History *** 
% Last version: 30/04/2015
% Created by J.Chatel-Goldman @GIPSA Lab, jonas.chatel.goldman(at)gmail.com


% default colors
if(nargin < 5)||(isempty(colBgd)) 
    colBgd  = [150 200 255];  
end
if(nargin < 4)||(isempty(colEdge)) 
    colEdge  = [0 0 255];  
end
colEdge = colEdge./255;
colBgd = colBgd./255;
% default cell of strings
if(nargin < 3)||(isempty(str)) 
    str = strsplit(int2str(1:length(x)))';  
    % strcat('mode #',strsplit(int2str(1:length(x)))'); % this is how to cat text if needed
end

plot(x,y,'s','markersize',16,'MarkerEdgeColor',colEdge, 'MarkerFaceColor',colBgd); 
for i = 1:length(x)
    text(x(i)-.1,y(i),str(i),'FontSize',12,'FontWeight','Bold'),  
end




% old way (bad as zooming change coordinates...)
% plot(x,y,' w'); % plot nothing (just to get normalized figure units)
% [xf,yf] = ds2nfu(x,y); % Convert to normalized figure units 
% for i = 1:length(x)
%     annotation('textbox',[xf(i)-.01 yf(i)-.01 .02 .02],'String',str(i),...
%         'HorizontalAlignment','center','VerticalAlignment','middle',...
%         'EdgeColor',colEdge,'BackgroundColor',colBgd,'FaceAlpha',.5,...
%         'FontWeight','Bold');
% end



% for a check of text annotation position:
% h = plot(x,y,'s','markersize',16,'MarkerEdgeColor','b', 'MarkerFaceColor',[150 200 255]./255); 

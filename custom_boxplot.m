function [f,b] = custom_boxplot(data,med_col,box_col,mrk_col,mrk_type,lwd,resc,x_lab,y_lab,fsize,fsize2,varargin)

% function for custom boxplots - especially to modify colors
%
% Inputs:
% data          data to be plotted, with one column per box
% med_col       colors of boxplot median line
% box_col       colors of boxplot patches
% mrk_col       colors of markers placed on median lines
% mrk_type      type of median marker - 'o' (circle),'d' (diamond), ...
%               (see "markertype" options for matlab "scatter" function
% lwd           linewidth
% resc          factor for rescaling median markers
% x_lab         x-axis label
% y_lab         y_axis label
% fsize         size of label fonts
% fsize2        size of axis fonts
%
% Outputs:
% f             figure handle
% b             boxplot handle
%
% Example use:
% data = [(2+randn(50,1)) randn(50,1) 3+randn(50,1) -1*randn(50,1) randn(50,1)];
% [b] = custom_boxplot(data);
%
% Author: 
% Frantisek Vasa (fv247@cam.ac.uk) - February 2016

nd = size(data,2);

% if ~(nd == length(box_col))
%     error('data size does not match');
% end

% default settings for parameters
if nargin < 11; fsize2 = 12; end
if nargin < 10; fsize = 12; end
if nargin < 9; y_lab = ''; end
if nargin < 8; x_lab = repmat({'box'},[nd 1]); end
if nargin < 7; resc = 1; end
if nargin < 6; lwd = 2; end
if nargin < 5; mrk_type = 'o'; end
if nargin < 4; mrk_col = repmat([1 0 0],[nd 1]); end % red
if nargin < 3; box_col = repmat([1 1 1],[nd 1]); end % black
if nargin < 2; med_col = [1 0 0]; end % red

f = figure
b = boxplot(data,'labels',x_lab,'symbol','k.','widths',0.5); % setting labels as blank is the easiest way to "remove" them

% change colors of boxplot patches
box_col = flipud(box_col);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    set(h(j),'Color','w');
    patch(get(h(j),'XData'),get(h(j),'YData'),box_col(j,:),'FaceAlpha',1);
end

% set linewidth on boxplot
set(b,{'linew'},{lwd},'Color',rgb('black'))

% change color of median lines
m = findobj(gca,'Tag','Median');
hold on;
for i = 1:1:nd
    plot(get(m(i),'XData')',get(m(i),'YData')','LineWidth',lwd,'Color',med_col); % median line
end
hold off;

% add circular marker on boxplot (separately from median line as ordering of the objects is nontrivial)
hold on;
for i = 1:1:nd
    s = scatter(i,nanmedian(data(:,i)),750*resc,mrk_type);
    set(s,'LineWidth',lwd,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',mrk_col(i,:));
end
hold off;

% set y-label
ylabel(y_lab,'FontName','Arial','FontSize',fsize);
t = findobj(gca,'Type','text'); set(t,'FontSize',fsize,'VerticalAlignment','top');
set(gca,'FontSize',fsize2);
set(gcf,'Color','w'); 

end



function [fig,ax] = pattern_summary(x_coords, y_coords, theta_sulc)
%PATTERN_SUMMARY  Plot summary of a given set of patterns, superimposed on approximation of sulcus.
%
% Syntax:
%   [fig,ax] = plot.pattern_summary(x_coords, y_coords, theta_sulc);
%
% Inputs:
%   x_coords - X-coordinates of each pattern focus (mm)
%   y_coords - Y-coordinates of each pattern focus (mm)
%
% Output:
%   fig - Figure handle
%   ax  - Axes handle
%
% See also: Content

if nargin < 1
    x_coords = [0,0,2,0,-2];
end

if nargin < 2
    y_coords = [0,2,0,-2,0];
end

if nargin < 3
    theta_sulc = -pi/12;
end

fig = figure(...
    'Name', 'Pattern Summary', ...
    'Color', 'w'); 
ax = axes(fig, 'NextPlot','add','XColor','none','YColor','none'); 
scatter(ax,x_coords,y_coords, ...
    'Color', 'r', ...
    'Marker', 'x', ...
    'LineWidth', 2.5, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'none', ...
    'DisplayName','Pattern Foci', ...
    'SizeData', 64);

A = 1.5*max([x_coords(:); y_coords(:)]);
l = line(ax,[A*cos(theta_sulc),A*cos(theta_sulc+pi)],[A*sin(theta_sulc)-1.5, A*sin(theta_sulc+pi)-1.5],...
    'Color','k','LineStyle',':','LineWidth',2.5, ...
    'DisplayName','CS (approx)');
l.Annotation.LegendInformation.IconDisplayStyle = 'on';

xx = [0.5*A*cos(pi/4),0.5*A*cos(pi/4)+1];
yy = [0.5*A*sin(pi/4),0.5*A*sin(pi/4)+1];
lx = line(ax, [xx(1), xx(2)],[yy(1), yy(1)], 'Color', [0.65 0.65 0.65], 'LineWidth', 2,'DisplayName','1-mm');
lx.Annotation.LegendInformation.IconDisplayStyle = 'on';

line(ax, [xx(1), xx(1)],[yy(1), yy(2)], 'Color', [0.65 0.65 0.65], 'LineWidth', 2,'DisplayName','1-mm');
leg = legend(ax, 'TextColor','Black', 'FontName', 'Tahoma','Location','Northwest');
end
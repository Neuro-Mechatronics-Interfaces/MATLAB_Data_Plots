function [fig, h, u] = raincloud(X, options)
%RAINCLOUD Plots a combination of half-violin, boxplot, and raw
%
% Adapted from original, which is at:
% https://github.com/RainCloudPlots/RainCloudPlots/blob/master/tutorial_matlab/raincloud_plot.m
%
% (Updated to use MATLAB `arguments` block)
% ---------------------------- INPUT ----------------------------
%
% X - vector of data to be plotted, required.
%
% --------------------- OPTIONAL ARGUMENTS ----------------------
%
% color             - color vector for rainclouds (default gray, i.e. = [.5 .5 .5])
% band_width         - band_width of smoothing kernel (default = 1)
% density_type       - choice of density algo ('ks' or 'rath'). Default = 'ks'
% box_on             - logical to turn box plots on/off (default = 0)
% box_dodge          - logical to turn on/off box plot dodging (default = 0)
% box_dodge_amount    - mutiplicative value to increase dodge amount (default = 0)
% alpha             - scalar positive value to increase cloud alpha (defalut = 1)
% dot_dodge_amount    - scalar value to increase dot dodge amounts (defalut =0.6)
% box_col_match       - logical to set it so that boxes match the colour of clouds (default = 0)
% line_width         - scalar value to set global line width (default = 2)
% lwr_bnd        - mutiplicative value to increase spacing at bottom of plot(default = 1)
% bxcl           - color of box outline
% bxfacecl       - color of box face
%
% ---------------------------- OUTPUT ----------------------------
% fig - Figure handle containing graphics children
% h   - Graphics object handles to change more stuff
% u   - parameter from kernel density estimatae
%
% ------------------------ EXAMPLE USAGE -------------------------
%
% fig = plot.raincloud('X', myData, 'box_on', 1, 'color', [0.5 0.5 0.5])
%
% Based on https://micahallen.org/2018/03/15/introducing-raincloud-plots/
% Inspired by https://m.xkcd.com/1967/
% v1 - Written by Tom Marshall. www.tomrmarshall.com
% v2 - Updated inputs to be more flexible - Micah Allen 12/08/2018
%
% Thanks to Jacob Bellmund for some improvements
% v3 - Max Murphy @ NML (2023-07-11)

arguments
    X double {isnumeric}
    options.Color (1,3) double {mustBeNumeric, mustBeInRange(options.Color, 0, 1)} = [0.5 0.5 0.5]
    options.BandWidth double {mustBeScalarOrEmpty} = []
    options.DensityType {mustBeTextScalar, mustBeMember(options.DensityType, {'ks', 'rash'})} = 'ks'
    options.Box (1,1) double {mustBeNumeric} = 0
    options.BoxDodge (1,1) double {mustBeNumeric} = 0
    options.BoxDodgeAmount (1,1) double {mustBeNumeric} = 0
    options.Alpha (1,1) double {mustBeNumeric} = 1
    options.DotDodge (1,1) double {mustBeNumeric} = 0.4
    options.MatchBoxColor (1,1) double {mustBeNumeric} = 0
    options.MaxNesting (1,1) double {mustBeNumeric} = 5
    options.LineWidth (1,1) double {mustBeNumeric} = 2
    options.LowerBand (1,1) double {mustBeNumeric} = 1
    options.BoxColor (1,3) double {mustBeNumeric, mustBeInRange(options.BoxColor,0,1)} = [0 0 0];
    options.CloudEdgeColor (1,3) double {mustBeNumeric, mustBeInRange(options.CloudEdgeColor,0,1)} = [0 0 0];
    options.Parent = []
    options.Support {mustBeTextScalar, mustBeMember(options.Support, {'unbounded', 'positive'})} = 'positive';
end

%% check all the inputs and if they do not exist then revert to default settings
% then set/get all the inputs out of this structure
color               = options.Color;
density_type        = options.DensityType;
box_on              = options.Box;
box_dodge           = options.BoxDodge;
box_dodge_amount    = options.BoxDodgeAmount;
alpha               = options.Alpha;
dot_dodge_amount    = options.DotDodge;
box_col_match       = options.MatchBoxColor;
line_width          = options.LineWidth;
lwr_bnd             = options.LowerBand;
box_color           = options.BoxColor;
cloud_edge_col      = options.CloudEdgeColor;
band_width          = options.BandWidth;

% calculate kernel density
switch density_type
    case 'ks'
        [f, Xi, u] = ksdensity(X, 'Bandwidth', band_width, 'Support', options.Support);
        
    case 'rash'
        % must have https://github.com/CPernet/Robust_Statistical_Toolbox
        % for this to work
        
        % check for rst_RASH function (from Robust stats toolbox) in path, fail if not found
        assert(exist('rst_RASH', 'file') == 2, 'Could not compute density using RASH method. Do you have the Robust Stats toolbox on your path?');
        
        [Xi, f] = rst_RASH(X);
        u = NaN; % not sure how to handle this with RASH yet
end

if isempty(options.Parent)
    fig = figure('Name', 'Raincloud Plot v3', 'Color', 'w');
    ax = axes(fig,'NextPlot','add','FontName','Tahoma','XColor','k','YColor','k');
else
    ax = options.Parent;
    set(ax, 'NextPlot','add','FontName','Tahoma');
    fig = ax.Parent;
    ii = 0;
    while (~isa(fig, 'matlab.ui.Figure')) && (ii < options.MaxNesting)
        fig = fig.Parent;
        ii = ii + 1;
    end
    if ~isa(fig,'matlab.ui.Figure')
        error("Parent axes is too nested. Increase MaxNesting option or check that Parent is correct axes.");
    end
end

% density plot
h{1} = area(ax, Xi, f); 
set(h{1}, 'FaceColor', color);
set(h{1}, 'EdgeColor', cloud_edge_col);
set(h{1}, 'LineWidth', line_width);
set(h{1}, 'FaceAlpha', alpha);

% make some space under the density plot for the boxplot and raindrops
yl = get(ax, 'YLim');
set(ax, 'YLim', [-yl(2)*lwr_bnd yl(2)]);

% width of boxplot
wdth = yl(2) * 0.25;

% jitter for raindrops
jit = (rand(size(X)) - 0.5) * wdth;

% info for making boxplot
quartiles   = quantile(X, [0.25 0.75 0.5]);
iqr         = quartiles(2) - quartiles(1);
Xs          = sort(X);
i_right_whisk = Xs > (quartiles(1) - (1.5 * iqr));
if sum(i_right_whisk)==0
    i_right_whisk = round(numel(Xs)/2);
end
whiskers(1) = min(Xs(i_right_whisk));
i_left_whisk = Xs < (quartiles(2) + (1.5 * iqr));
if sum(i_left_whisk) == 0
    i_left_whisk = i_right_whisk;
end
whiskers(2) = max(Xs(i_left_whisk));
Y           = [quartiles whiskers];

% raindrops
if box_dodge
    drops_pos = (jit * 0.6) - yl(2) * dot_dodge_amount;
else
    drops_pos = jit - yl(2) / 2;
end

h{2} = scatter(ax, X, drops_pos);
h{2}.SizeData = 10;
h{2}.MarkerFaceColor = color;
h{2}.MarkerEdgeColor = 'none';

if box_on
    
    if box_col_match
        
        box_color = color;
        
    end
    
    if box_dodge
        box_pos = [Y(1) ((-yl(2) * box_dodge_amount) - (wdth * 0.3)) Y(2) - Y(1) (wdth * 0.6)];
        
        % median line
        h{4} = line(ax, [Y(3) Y(3)], [((-yl(2) * box_dodge_amount) - (wdth * 0.3)) ((-yl(2) * box_dodge_amount) + (wdth * 0.3))], 'col', box_color, 'LineWidth', line_width);
        
        % whiskers
        h{5} = line(ax, [Y(2) Y(5)], [(-yl(2) * box_dodge_amount) (-yl(2) * box_dodge_amount)], 'col', box_color, 'LineWidth', line_width);
        h{6} = line(ax, [Y(1) Y(4)], [(-yl(2) * box_dodge_amount) (-yl(2) * box_dodge_amount)], 'col', box_color, 'LineWidth', line_width);
    else
        box_pos = [Y(1) -yl(2)/2-(wdth * 0.5) Y(2)-Y(1) wdth];
        % median line
        h{4} = line(ax, [Y(3) Y(3)], [-yl(2)/2-(wdth * 0.5) -yl(2) / 2 + (wdth * 0.5)], 'col', box_color, 'LineWidth', line_width);
        
        % whiskers
        h{5} = line(ax, [Y(2) Y(5)], [-yl(2)/2 -yl(2)/2], 'col', box_color, 'LineWidth', line_width);
        h{6} = line(ax, [Y(1) Y(4)], [-yl(2)/2 -yl(2)/2], 'col', box_color, 'LineWidth', line_width);
    end
    % 'box' of 'boxplot'
    h{3} = rectangle(ax, 'Position', box_pos);
    set(h{3}, 'EdgeColor', box_color)
    set(h{3}, 'LineWidth', line_width);
    %set(h{3}, 'FaceColor', bxfacecl);
    
    ylim('manual') % set the ylim to capture all data
end
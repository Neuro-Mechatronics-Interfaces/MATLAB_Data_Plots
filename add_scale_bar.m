function h = add_scale_bar(ax, x0, y0, x1, y1, options)
%ADD_SCALE_BAR Adds a scale bar to an existing axis, with customizable units, color, and position.
%
% Syntax:
%   h = add_scale_bar(ax, x0, y0, x1, y1, options)
%
% Description:
%   This function adds a scale bar annotation to a specified axis, `ax`, allowing
%   customizable labels, colors, and units for both horizontal and vertical bars.
%   The function uses start and end points (`x0`, `y0`, `x1`, `y1`) to set the scale
%   bar dimensions, and labels the bars based on specified units and scale factors.
%
% Inputs:
%   ax  - Handle to the axis where the scale bar will be added.
%   x0  - Starting x-coordinate of the horizontal scale bar.
%   y0  - Starting y-coordinate of the vertical scale bar.
%   x1  - Ending x-coordinate of the horizontal scale bar.
%   y1  - Ending y-coordinate of the vertical scale bar.
%
% Options:
%   Color                - RGB color for the scale bar lines and text (default: [0.65 0.65 0.65]).
%   XBar                 - Logical to display the horizontal scale bar (default: true).
%   XUnits               - String label for horizontal units (default: 'sec').
%   XLabelScaleFactor    - Scaling factor for horizontal label values (default: 1).
%   XLabelRoundingLevel  - Integer rounding precision for horizontal label values (default: 0).
%   YBar                 - Logical to display the vertical scale bar (default: true).
%   YUnits               - String label for vertical units (default: 'μV').
%   YLabelScaleFactor    - Scaling factor for vertical label values (default: 1).
%   YLabelRoundingLevel  - Integer rounding precision for vertical label values (default: 0).
%   FontSize             - Font size for the scale bar labels (default: 12).
%   FontName             - Font name for the scale bar labels (default: 'Tahoma').
%
% Outputs:
%   h - Structure containing handles for each element of the scale bar:
%       * h.XBar  - Handle for the horizontal bar line (if displayed).
%       * h.XText - Handle for the horizontal bar label (if displayed).
%       * h.YBar  - Handle for the vertical bar line (if displayed).
%       * h.YText - Handle for the vertical bar label (if displayed).
%
% Example:
%   % Create a plot and add a scale bar with default settings
%   figure; ax = axes;
%   plot(rand(10, 1)); % Example plot
%   h = add_scale_bar(ax, 0.1, 0.1, 0.3, 0.5);
%
%   % Add a custom scale bar with specified colors, units, and labels
%   h = add_scale_bar(ax, 0.1, 0.1, 0.4, 0.5, 'Color', [0 0 1], ...
%                     'XUnits', 'ms', 'YUnits', 'mV', 'FontSize', 10, ...
%                     'XLabelScaleFactor', 1000, 'YLabelScaleFactor', 1000);
%
% Notes:
%   - The scale bar lengths are derived from the differences between `x0`, `x1` and `y0`, `y1`.
%   - Units, scaling, and label rounding are customizable to accommodate various contexts.
%   - Scale bar elements can be hidden by setting `XBar` or `YBar` options to false.
%
% See also: line, text, axes

arguments
    ax
    x0 (1,1) % Starting point of horizontal bar
    y0 (1,1) % Starting point of vertical bar
    x1 (1,1) % Ending point of horizontal bar
    y1 (1,1) % Ending point of vertical bar
    options.Color (1,3) double {mustBeInRange(options.Color,0,1)} = [0.65 0.65 0.65];
    options.XBar (1,1) logical = true;
    options.XUnits = 'sec';
    options.XLabelScaleFactor = 1;
    options.YBar (1,1) logical = true;
    options.YUnits = 'μV';
    options.YLabelScaleFactor = 1;
    options.FontSize = 12;
    options.FontName {mustBeTextScalar} = 'Tahoma';
    options.XLabelRoundingLevel (1,1) {mustBeInteger} = 0;
    options.YLabelRoundingLevel (1,1) {mustBeInteger} = 0;
end
% Add scalebar to plot
dx = abs(x1 - x0);
if ~isduration(dx)
    dx = round(dx.*options.XLabelScaleFactor,options.XLabelRoundingLevel)./options.XLabelScaleFactor;
    tx = sprintf('%g %s', dx * options.XLabelScaleFactor, options.XUnits); % time scalebar text
else
    dx = seconds(round(seconds(dx).*options.XLabelScaleFactor,options.XLabelRoundingLevel)./options.XLabelScaleFactor);
    tx = sprintf('%g %s', seconds(dx) * options.XLabelScaleFactor, options.XUnits); % time scalebar text
end
x1 = x0 + dx;
dy = abs(y1 - y0);
if ~isduration(dy)
    dy = round(dy.*options.YLabelScaleFactor, options.YLabelRoundingLevel)./options.YLabelScaleFactor;
    ty = sprintf('%g %s', dy * options.YLabelScaleFactor, options.YUnits); % vertical scalebar text
else
    dy = seconds(round(seconds(dy).*options.YLabelScaleFactor, options.YLabelRoundingLevel)./options.YLabelScaleFactor);
    ty = sprintf('%g %s', seconds(dy) * options.YLabelScaleFactor, options.YUnits); % vertical scalebar text
end
y1 = y0 + dy;

if options.XBar
    h.XBar = line(ax, [x0, x1], [y0, y0], 'Color',options.Color,'LineWidth',1.25,'LineStyle','-');
    h.XText = text(ax, x0+dx/2, y0-0.05*dy, tx, 'FontName',options.FontName,'HorizontalAlignment','center','VerticalAlignment','top','Color',options.Color,'FontSize',options.FontSize);
    h.YBar.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
if options.YBar
    h.YBar = line(ax, [x0, x0], [y0, y1], 'Color',options.Color,'LineWidth',1.25,'LineStyle','-');
    h.YText = text(ax, x0-0.05*dx, y0+dy/2, ty, 'FontName',options.FontName,'HorizontalAlignment','center','VerticalAlignment','bottom','Rotation',90,'Color',options.Color,'FontSize',options.FontSize);
    h.YBar.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
end
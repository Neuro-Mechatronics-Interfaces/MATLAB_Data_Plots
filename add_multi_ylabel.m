function h = add_multi_ylabel(ax, labels, options)
%ADD_MULTI_YLABEL Adds multiple y-axis labels to an axis with customized alignment, color, spacing, and font options.
%
% Syntax:
%   h = add_multi_ylabel(ax, labels, options)
%
% Description:
%   This function adds multiple y-axis labels to the specified axis, `ax`, with customizable options for 
%   alignment, color, font, and spacing. The labels can be aligned to either the left or right of the y-axis, 
%   spaced evenly within a specified y-limit range or distributed between the bottom and top axis limits.
%
% Inputs:
%   ax      - Handle to the axis to which labels will be added.
%   labels  - String array of label names to add along the y-axis.
%
% Options:
%   Color                - RGB color for each label, specified as either a single color (applied to all labels)
%                          or a unique color for each label (default: [0 0 0] for black).
%   Alignment            - Alignment side for labels: 'left' or 'right' (default: 'left').
%   OffsetXDistanceFraction - Distance to offset labels along the x-axis, relative to the x-limits of `ax`.
%                          Positive values place labels farther outside (default: 0.035).
%   FontName             - Font for the labels (default: "Consolas").
%   FontSize             - Font size for the labels (default: 8).
%   FontOptions          - Cell array of additional font options applied to each label 
%                          (e.g., {'FontWeight', 'bold'}, default: empty).
%   Spacing              - Controls vertical spacing of labels:
%                          'even' spaces labels from bottom to top with an even offset,
%                          'between' spaces labels directly between the bottom and top limits (default: 'even').
%   SpaceLim             - Limits for y-axis placement of the labels, specified as a two-element vector [min, max].
%                          Default values use the axis y-limits.
%
% Output:
%   h - Array of text object handles for each label.
%
% Example:
%   % Set up axis and add multiple y-axis labels with default options
%   ax = gca;
%   labels = ["Label 1", "Label 2", "Label 3"];
%   h = add_multi_ylabel(ax, labels);
%
%   % Customize y-axis labels with colors, right alignment, and larger font
%   colors = [0 0.5 0; 0 0 1; 1 0 0]; % Green, Blue, Red
%   h = add_multi_ylabel(ax, labels, 'Color', colors, 'Alignment', 'right', 'FontSize', 12);
%
%   % Add evenly spaced labels with custom y-axis limits and font options
%   h = add_multi_ylabel(ax, labels, 'Spacing', 'between', 'SpaceLim', [0, 10], 'FontOptions', {'FontWeight', 'bold'});
%
% Notes:
%   - Labels are rotated 90 degrees for vertical orientation.
%   - If `Color` is specified with multiple rows, the number of rows must match the number of labels.
%
% See also: text, axes, linspace
arguments
    ax (1,1)
    labels (1,:) string
    options.Color (:,3) double {mustBeInRange(options.Color,0,1)} = [0 0 0]; % Single color for each label, or can be given as unique color for each element in labels array
    options.Alignment {mustBeMember(options.Alignment, {'left', 'right'})} = 'left'; % Which y-axis side to align to
    options.OffsetXDistanceFraction = 0.035; % Ratio of XLim distance to offset (positive is further "outside" for given alignment)
    options.FontName {mustBeTextScalar} = "Consolas";
    options.FontSize (1,1) double {mustBePositive} = 8;
    options.FontOptions cell = {}; % Add any other text label properties you want here, they are applied to each label (i.e. {'FontWeight', 'bold', ...})
    options.Spacing {mustBeMember(options.Spacing, {'even', 'between'})} = 'even'; % 'even' - Offset from bottom/top of axes, with even spacing. 'between' - Even spacing, from exact bottom to exact top.
    options.SpaceLim (1,2) = [-inf, inf];
end

dx = diff(ax.XLim);
n = numel(labels);
spaceLim = options.SpaceLim;
if isinf(spaceLim(1))
    spaceLim(1) = ax.YLim(1);
end
if isinf(spaceLim(2))
    spaceLim(2) = ax.YLim(2);
end

switch options.Alignment
    case 'left'
        x0 = ax.XLim(1) - options.OffsetXDistanceFraction*dx;
    case 'right'
        x0 = ax.XLim(2) + options.OffsetXDistanceFraction*dx;
end

switch options.Spacing
    case 'even'
        y = linspace(spaceLim(1),spaceLim(2), n+1);
        y = y(1:(end-1)) + mean(diff(y))/2;
    case 'between'
        y = linspace(spaceLim(1), spaceLim(2), n);
end
if size(options.Color,1)==1
    cols = repmat(options.Color,n,1);
else
    cols = options.Color;
    if size(cols,1)~=n
        error("Must give only 1 Color or exactly matched colors for number of elements in labels (current: %d labels, but specified %d colors).", n, size(cols,1));
    end
end

h = gobjects(n,1);
for ii = 1:n
    h(ii) = text(ax, x0, y(ii), labels(ii), ...
        'FontName',options.FontName, ...
        'FontSize', options.FontSize, ...
        'Rotation', 90, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Color', cols(ii,:), ...
        'HitTest', 'off', ...
        'AffectAutoLimits', 'off', ...
        options.FontOptions{:});
end

end
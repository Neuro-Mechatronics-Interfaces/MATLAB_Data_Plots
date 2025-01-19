function h = add_multi_ylabel(ax, labels, options)
%ADD_MULTI_YLABEL Add multiple custom y-axis labels to an axis.
%
% Syntax:
%   h = plot.add_multi_ylabel(ax, labels);
%   h = plot.add_multi_ylabel(ax, labels, 'Name', value, ...);
%
% Description:
%   This function adds multiple text labels along the y-axis of a specified axis object (`ax`), either on the 
%   left or right side. The labels are evenly spaced or distributed between specified limits, with customizable 
%   text properties such as font, size, color, and alignment. It is useful for annotating plots with multiple 
%   descriptors along the y-axis.
%
% Inputs:
%   ax       - Handle to the axis where the y-axis labels will be added (scalar).
%   labels   - Array of label strings to display along the y-axis (1D string array).
%   options  - (Optional) Name-value arguments for customization:
%       Color                 - RGB color for the labels, either as a single color for all labels or a matrix 
%                               with one row per label (default: [0 0 0]).
%       Alignment             - Which side of the y-axis to align the labels ('left' or 'right', default: 'left').
%       OffsetXDistanceFraction - Offset of labels relative to the axis limits, as a fraction of the x-axis 
%                               range (default: 0.035).
%       FontName              - Font name for the labels (default: "Consolas").
%       FontSize              - Font size for the labels (default: 8).
%       FontOptions           - Cell array of additional text properties (default: {}, e.g., {'FontWeight', 'bold'}).
%       Spacing               - Spacing mode for the labels:
%                                 - 'even'    : Evenly spaced within `SpaceLim` (default).
%                                 - 'between' : Exactly at boundaries defined by `SpaceLim`.
%       SpaceLim              - [lower, upper] limits for label placement along the y-axis (default: [-inf, inf]).
%
% Outputs:
%   h        - Array of handles to the text objects created for the labels.
%
% Notes:
%   - The labels are rotated vertically (90 degrees) and placed outside the axis limits on the specified side.
%   - The `Color` input can be a single RGB triplet for all labels or a matrix with one row per label.
%   - The `SpaceLim` parameter can be used to restrict the placement of labels within a specific y-axis range.
%   - Labels do not affect the axis limits or respond to mouse clicks.
%
% Example:
%   % Example 1: Add labels to the left side of the y-axis
%   ax = axes;
%   plot(rand(10,1));
%   labels = ["Low", "Medium", "High"];
%   h = plot.add_multi_ylabel(ax, labels);
%
%   % Example 2: Add labels with custom spacing and colors
%   labels = ["A", "B", "C", "D"];
%   cols = [1 0 0; 0 1 0; 0 0 1; 0.5 0.5 0.5];
%   h = plot.add_multi_ylabel(ax, labels, ...
%               'Color', cols, ...
%               'Alignment','right',
%               'Spacing','between');
%
% See also: text, ylabel, axes

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
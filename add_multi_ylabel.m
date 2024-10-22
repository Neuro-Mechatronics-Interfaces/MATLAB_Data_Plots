function h = add_multi_ylabel(ax, labels, options)
%ADD_MULTI_YLABEL Adds multiple ylabels 
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
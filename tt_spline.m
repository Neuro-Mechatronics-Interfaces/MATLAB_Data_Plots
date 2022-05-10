function [ax, curve, goodness, output] = tt_spline(ax, TT, varName)
%TT_SPLINE  Plot smoothing spline fitted model from time table.
%
% Syntax:
%   ax = plot.tt_spline(TT, varName);
%   ax = plot.tt_spline(ax, TT, varName);
%
% Inputs:
%   ax  --  Optional, axes to plot on
%   TT  --  Timetable
%   varName -- Name of variable to plot against time from table
%
% Output:
%   ax  --  The axes that was plotted on
%   [curve, goodness, output] -- See `fit` builtin.
%
% See also: Contents

if nargin < 3
    varName = TT;
    TT = ax;
    ax = gca;
end

if numel(varName) > 1
    curve = cell(size(varName));
    goodness = cell(size(varName));
    output = cell(size(varName));
    for ii = 1:numel(varName)
        [ax, curve{ii}, goodness{ii}, output{ii}] = plot.tt_spline(ax, TT, varName(ii));
    end
    return;
end

dt = datenum(TT.Date);
data = TT.(varName);
mask = ~isnan(data);
[curve, goodness, output] = fit(dt(mask), data(mask), ...
    'smoothingspline');
axes(ax);
plot(curve, dt(mask), data(mask));
xlim(ax, [min(dt), max(dt)]);
labs = arrayfun(@(a)string(datestr(a)),ax.XTick);
set(ax, ...
    'NextPlot', 'add', ...
    'XTickLabel', labs, ...
    'XTickLabelRotation', 30, ...
    'XColor', 'k', 'YColor', 'k', 'LineWidth', 1.5, ...
    'FontName', 'Tahoma', 'FontSize', 11);
col = rand(1,3);
h = findobj(get(ax,'Children'), 'DisplayName', 'fitted curve');
set(h, 'LineWidth', 1.5, 'Tag', varName, 'DisplayName', sprintf('Fitted: %s', strrep(varName, '_', ' ')), 'Color', col);
h = findobj(get(ax,'Children'), 'DisplayName', 'data');
set(h, 'Tag', sprintf('data-%s', varName), 'DisplayName', sprintf('Data: %s', strrep(varName, '_', ' ')), 'Color', col);
xlabel(ax, 'Date', 'FontName', 'Tahoma', 'FontSize', 14, 'Color', 'k');
ylabel(ax, '');
legend(ax, 'Location','Best');
end
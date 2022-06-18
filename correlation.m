function fig = correlation(ax, tau, rho, varargin)
%CORRELATION  Plot correlation (optionally specifying an axes handle)
%
% Syntax:
%   plot.correlation(tau, rho);
%   fig = plot.correlation(ax, tau, rho, 'Name', value, ...);
%
% Inputs:
%   ax   - (Optional) target axes to plot on
%   tau  - Time lag values for correlation vector
%   rho  - Cross-correlation for each lag between two vectors of interest
%   varargin - (Optional) 'Name', value pairs. See default MATLAB plot
%               options for valid arguments.
%
% Output:
%   fig  - Figure handle
%
% See also: Contents

if isa(ax, 'matlab.graphics.axis.Axes')
    fig = get(ax, 'Parent');
    if ~isa(fig, 'matlab.ui.Figure')
        fig = get(fig, 'Parent');
    end
else
    if nargin >= 3
        varargin = [rho, varargin];
    end
    rho = tau;
    tau = ax;
    fig = figure('Name', 'Cross-Correlation', 'Color', 'w', 'Position', [250   150   375   372]);
    ax = axes(fig, 'NextPlot', 'add', 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1.5, 'FontName', 'Tahoma');
    xlabel(ax, '\tau (sec)', 'FontName', 'Tahoma', 'Color', 'k');
    ylabel(ax, '\rho_x_y(\tau)', 'FontName', 'Tahoma', 'Color', 'k');
    title(ax, 'Cross-Correlation Function', 'FontName', 'Tahoma', 'Color', 'k');
end

if isa(tau, 'duration')
    tau = seconds(tau);
end

if max(rho) > 1
    rho = rho ./ numel(tau); 
end

plot(ax, tau, rho, 'LineWidth', 1.5, 'Color', 'k', 'DisplayName', '\rho_x_y(\tau)', varargin{:});
[rho_max, i_max] = max(rho);
scatter(ax, tau(i_max), rho_max, 'SizeData', 14, 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', ...
    'DisplayName', sprintf("Max: <%4.2f sec, %g>", tau(i_max), rho_max));
legend(ax, 'FontName', 'Tahoma', "TextColor", 'black', 'Location', 'southoutside');
end
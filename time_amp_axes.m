function [fig, ax] = time_amp_axes(name, options)
%TIME_AMP_AXES  Return figure and axes handle for time-vs-amplitude axes
%
% Syntax:
%   [fig,ax] = plot.time_amp_axes(name);
%   [fig,ax] = plot.time_amp_axes(name, 'Name', value, ...);
%
% Inputs:
%   name - String, the title of the plot (and figure window name).
%
% Options:
%   XLim (1,2) double
%   YLim (1,2) double
%   FontName (1,1) string = "Tahoma"
%   ColorOrder (:,3) double = jet(100)
%   Position (1,4) double = [160   130   640   720]
%   
% Output:
%   fig - Figure handle
%   ax  - Axes handle
%
% See also: Contents, plot, plot.stacked_emg

arguments
    name (1,1) string = "Untitled";
    options.XLim (1,2) double = [nan nan];
    options.YLim (1,2) double = [nan nan];
    options.FontName (1,1) string = "Tahoma";
    options.ColorOrder (:,3) double = [0 0 0; 0 0 1; 1 0 0; 0.65 0.65 0.65];
    options.Position (1,4) double = [300   130   640   360];
    options.YUnits (1,1) string = "μV"
    options.XUnits (1,1) string = "ms"
end

fig = figure('Name', name, 'Color', 'w', 'Position', options.Position);
ax = axes(fig,...
    'NextPlot','add',...
    'FontName', options.FontName, ...
    'ColorOrder',options.ColorOrder);
xlabel(ax,sprintf('Time (%s)', options.XUnits),'FontName',options.FontName,'FontSize',14,'Color','k');
ylabel(ax,sprintf('Amplitude (%s)',options.YUnits),'FontName',options.FontName,'FontSize',14,'Color','k');
title(ax, strrep(name, '_', '-'), ...
    'FontName',options.FontName,'FontSize',20,'FontWeight','bold','Color','k');
if ~isnan(options.XLim(1))
    xlim(ax, options.XLim);
end
if ~isnan(options.YLim(1))
    ylim(ax, options.YLim);
end

end
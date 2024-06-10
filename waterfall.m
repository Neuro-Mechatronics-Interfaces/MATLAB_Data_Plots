function [fig, ax, h] = waterfall(t, Y, options)
%WATERFALL Make 2D waterfall plot.
%
% Syntax:
%   [fig, ax, h] = plot.waterfall(t, Y, 'Name', value, ...);
%
% Inputs:
%   t - Independent variable values (e.g. timesteps or spatial position). 
%           Should have same number of elements as rows of Y.
%   Y - nIndependentVarSamples x nRepetitions dependent variable traces to plot.
%
% Options:
%   'Parent' (def: []) Parent axes handle that waterfall plot goes onto.
%   'XLabel' (def: 'Time (ms)')
%   'YLabel' (def: 'Amplitude (µV)' Only added if AddScaleBar is true
%   'AddScaleBar' (def: true)
%   'CData'  (def: []) Color data either [1 x 3] or [nRepetitions x 3] to
%                       specify the axes ColorOrder property, setting the
%                       color of each line.
%   'YOffset' (def: 1000) The vertical offset between each trace.
%   'YScale'  (def: []) If empty, default is 0.1 * nRepetitions * YOffset.
%                       Only used if 'AddScaleBar' is true.
%   
% Output:
%   fig - Figure handle
%   ax  - Parent axes handle
%   h   - Array of line object handles. First element is "lowest" on y-axis.
%
% See also: Contents

arguments
    t (1,:) double
    Y (:,:) double
    options.Parent = [];
    options.XUnits {mustBeTextScalar} = 'ms';
    options.YUnits {mustBeTextScalar} = 'µV';
    options.YScale = [];
    options.AddScaleBar (1,1) logical = true;
    options.AddZeroMarker (1,1) logical = true;
    options.CData = [];
    options.YOffset (1,1) double = 100;
    options.Title {mustBeTextScalar} = "";
    options.Subtitle {mustBeTextScalar} = "";
    options.LineWidth (1,1) double = 1.25;
    options.LineOptions cell = {};
    options.FontName {mustBeTextScalar} = 'Tahoma';
    options.FontSize (1,1) double = 14;
    options.XLabelRoundingLevel (1,1) {mustBeInteger} = 0;
    options.YLabelRoundingLevel (1,1) {mustBeInteger} = 0;
    options.ZeroMarkerLabel {mustBeTextScalar} = 'Stim Onset';
    options.ZeroMarkerOptions cell = {};
end

if isempty(options.CData)
    cdata = spring(size(Y,2));
else
    cdata = options.CData;
    if size(cdata,1) == 1
        cdata = double(cm.umap(cdata,size(Y,2)))./255.0;
    elseif size(cdata,1)~=size(Y,2)
        error("Must have same number of rows in CData option as columns of Y.");
    end
end

if isempty(options.Parent)
    fig = figure('Name','Waterfall Plot', ...
        'Color','w', ...
        'Units','inches',...
        'Position',[0.5 0.5 3.25 5]);
    ax = axes(fig,...
        'NextPlot','add',...
        'FontName',options.FontName, ...
        'FontSize', options.FontSize, ...
        'YColor','none',...
        'Box','off',...
        'XColor','none');
else
    ax = options.Parent;
    fig = ax.Parent;
    iter = 0;
    while ~isa(fig,'matlab.ui.Figure') && (iter < 3) % Can only nest "3-deep"
        fig = fig.Parent;
        iter = iter + 1;
    end
end
ax.ColorOrder = cdata;

y0 = 0:options.YOffset:((size(Y,2)-1)*options.YOffset);
h = plot(ax,t,Y+y0, ...
    'LineWidth',options.LineWidth, ...
    options.LineOptions{:});

w0 = t(end)-t(1);
x0 = t(1)-0.025*w0;
h0 = round(0.1 * (size(Y,2)-1) * options.YOffset);

if options.AddScaleBar
    plot.add_scale_bar(ax,x0,0-options.YOffset/2,0,h0-options.YOffset/2, ...
        'XUnits',options.XUnits, ...
        'YUnits',options.YUnits, ...
        'FontSize', options.FontSize-4, ...
        'XLabelRoundingLevel', options.XLabelRoundingLevel, ...
        'YLabelRoundingLevel', options.YLabelRoundingLevel);
end

if options.AddZeroMarker
    xline(ax,0,...
        'Label',options.ZeroMarkerLabel,...
        'FontSize',options.FontSize,...
        'FontWeight','bold',...
        'Color','k', ...
        'LineStyle','--',...
        'LabelOrientation','aligned', ...
        'LabelHorizontalAlignment','left',...
        'LabelVerticalAlignment','middle',...
        options.ZeroMarkerOptions{:});
end

plot.add_titles(ax, options.Title, options.Subtitle, ...
    'FontName', options.FontName, ...
    'FontSize', options.FontSize);

end
function fig = raster(data, varargin)
%RASTER  Plot raster pulse trains from cell array of sample times
%
% Syntax:
%   fig = plot.raster(data, 'Name', value, ...);
%
% Inputs:
%   data - Cell array where each element is a vector of sample indices
%           denoting the time from the start of the record where there was
%           an event (such as MUAP or Spike) of interest.
%   varargin - (Optional) 'Name',value input argument pairs.
%         + 'CData' | [] | If empty then uses 'jet' colormap for pulse coloring
%               (or colormap specified in Colormap parameter; otherwise, give as
%                n x 3 array of colors where n is number of cells in `data`)
%         + 'Colormap' | @jet | The colormap function to use if no CData.
%         + 'Labels' | '' | Tick Label (1 for each cell in data)
%         + 'Position' | [200, 200, 400, 400] | Figure position
%         + 'PulseHeight' | 0.8 |  Fraction from 0-1 for height of each line. 
%           - (Can be set as scalar or vector for each cell element in data.)
%           - (Can also be set as a cell array, with individual elements
%                   matching to the "pulse height" of individual ticks in 
%                   the raster "pulse train"). 
%         + 'PulseLineWidth' | 0.9 | Line Width parameter for each pulse tick
%         + 'Reference' | [] | Reference signal(s)
%         + 'ReferenceColor' | [0.7 0.7 0.7] | Color of lightest reference (or nRef x 3 array of colors for each reference line).
%         + 'ReferenceLabels' | {} | Reference signal label(s)
%         + 'SampleRate' | 4000 |  Default from TMSi recordings
%         + 'Title' | '' | Title of figure and axes
%         + 'XLim' | [] | The x-limits set manually on figure axes.
%
% Output:
%   fig - Figure handle generated from data.
%
% See also: Contents, +ckc

p = inputParser;
p.addParameter('CData', []);
p.addParameter('Colormap', @jet);
p.addParameter('Labels', '');
p.addParameter('Position', [200, 200, 1000, 640]);
p.addParameter('PulseHeight', 0.8); % Fraction from 0-1 for height of each line. (Can be set as scalar or vector for each cell element in data.)
p.addParameter('PulseLineWidth', 0.9); % 'LineWidth' parameter for each line tick
p.addParameter('Reference', []); % Reference signal(s)
p.addParameter('ReferenceColor', [0.7, 0.7, 0.7]); % Color of lightest reference (or nRef x 3 array of colors for each reference line).
p.addParameter('ReferenceLabels', {});
p.addParameter('SampleRate', 4000); % Default from TMSi recordings
p.addParameter('Title', '');
p.addParameter('XLim', []); 
p.parse(varargin{:});
N = numel(data);

if isempty(p.Results.Title)
    namestr = 'Raster';
else
    namestr = sprintf('Raster: %s', p.Results.Title);
end

fig = figure(...
    'Name', namestr, ...
    'Color', 'w', ...
    'Position', p.Results.Position);
ax = axes(fig, ...
    'NextPlot', 'add', ...
    'FontName', 'Tahoma', ...
    'XColor', 'k', 'YColor', 'k', ...
    'YTick', 1:N, ...
    'YLim', [0.5, N+0.5], ...
    'YTickLabelRotation', 30, ...
    'Box', 'on');
xlabel(ax, 'Time (ms)', 'FontName', 'Tahoma', 'Color', 'k');
title(ax, namestr, 'FontName', 'Tahoma', 'Color', 'k');
if isempty(p.Results.Labels)
    ylabel(ax, '(Pulse Train)', ...
        'FontName', 'Tahoma', ...
        'FontSize', 10, ...
        'Color', [0.75 0.75 0.75]);
else
    ax.YTickLabels = p.Results.Labels;
end

fs = p.Results.SampleRate * 1e-3; % Sample rate, but samples per millisecond

nRef = size(p.Results.Reference, 1);
if nRef > 0
    tRef = (0:(size(p.Results.Reference,2)-1)) ./ fs;
    
    sRef = N / nRef; % Reference signals distributed along y-axis.
    lRef = string(p.Results.ReferenceLabels);
    if (nRef > 1)
        if numel(lRef) == 1
            lRef = repmat(lRef, nRef, 1);
        end
        if size(p.Results.ReferenceColor,1) == 1
            refColor = repmat(p.Results.ReferenceColor, nRef, 1) ./ ((1:nRef)');
        else
            refColor = p.Results.ReferenceColor;
        end
    else
        refColor = p.Results.ReferenceColor;
    end
    if isempty(lRef)
        XL = [tRef(1), tRef(end)]; % No need to "make room" for legend
    else
        XL = [tRef(1), tRef(end)+(tRef(end)-tRef(1))*0.2];
    end
    
    for ii = 1:nRef
        ref = p.Results.Reference(ii,:) - min(p.Results.Reference(ii,:));
        ref = ref ./ max(ref);
        yRef = sRef .* ref + sRef*(ii-1) + 0.5;
        cRef = refColor(ii,:);
        hRef = plot(ax, tRef, yRef, ...
            'LineWidth', max(5 - ii, 1), ...
            'Color', cRef);
        if ~isempty(lRef)
            hRef.DisplayName = lRef(ii);
        end
    end
else
    XL = [inf, -inf];
    lRef = string({});
end

% Use this to determine the x-bounds (if no ref. signal already).
xl = [inf, -inf];

% Set the pulse height for the individual pulse trains.
if ~iscell(p.Results.PulseHeight(1))
    hh = p.Results.PulseHeight ./ 2;
    if isscalar(hh)
        hh = repmat(hh, N, 1);
    end
else
    hh = p.Results.PulseHeight;
end

% Set cdata for each pulse train
if isempty(p.Results.CData)
    cdata = feval(p.Results.Colormap, N);
else
    cdata = p.Results.CData;
end

for ii = 1:N
    ts = sort(data{ii} ./ fs, 'ascend');
    n = numel(ts);
    ts = reshape(ts, 1, n);
    xl = [min(xl(1), ts(1)), max(xl(2), ts(end))];
    xts = [repmat(ts, 2, 1); nan(1, n)];
    if iscell(hh(ii))
        hhh = hh{ii} ./ 2;
    else
        hhh = hh(ii,:);
    end
    yts = [ones(1,n).*(ii-hhh); ...
           ones(1,n).*(ii+hhh); ...
           nan(1, n)];
    hl = line(ax, xts(:), yts(:), ...
        'Color', cdata(ii,:), ...
        'LineWidth', p.Results.PulseLineWidth);
    hl.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

if isempty(lRef)
    XL = [max(xl(1), XL(1)), min(xl(2), XL(2))];
else
    XL = [max(xl(1), XL(1)), min(xl(2)+diff(xl)*0.2, XL(2))];
    legend(ax, ...
        'FontName', 'Tahoma', ...
        'TextColor', 'black', ...
        'FontSize', 8, ...
        'Location', 'NorthEast');
end
if isempty(p.Results.XLim)
    ax.XLim = XL;
else
    ax.XLim = p.Results.XLim;
end

end
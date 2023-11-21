function [fig, L] = raster(data, options)
%RASTER  Plot raster pulse trains from cell array of sample times
%
% Syntax:
%   fig = plot.raster(data, 'Name', value, ...);
%   [fig, L] = plot.raster(__);
%
% Inputs:
%   data - Cell array where each element is a vector of sample indices
%           denoting the time from the start of the record where there was
%           an event (such as MUAP or Spike) of interest.
%
%  Options:
%         + 'BackgroundColor' | 'w' | Color in background of figure/axes.
%         + 'CData' | [] | If empty then uses 'jet' colormap for pulse coloring
%               (or colormap specified in Colormap parameter; otherwise, give as
%                n x 3 array of colors where n is number of cells in `data`)
%         + 'Colormap' | @jet | The colormap function to use if no CData.
%         + 'EventWaveform' | {} | Cell array, if included should have 1:1 correspondence of template event-waveforms to elements in data.
%         + 'Figure' | [] | Set this to figure handle to populate existing figure.
%         + 'Labels' | '' | Tick Label (1 for each cell in data)
%         + 'Layout' | [] | Set this to tiledlayout to populate existing layout.
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
%         + 'SmoothKernelWidth', | 0.25 | Width of smoothing kernel (seconds) 
%         + 'Title' | '' | Title of figure and axes
%         + 'UserData' | struct | UserData property of the returned figure handle.
%         + 'YLabel' | '(pulse train)' | YLabel text
%         + 'XLim' | [] | The x-limits set manually on figure axes.
%
% Output:
%   fig - Figure handle generated from data.
%   L   - Handle to tiledlayout graphics object at top-level of figure.
%
% See also: Contents, +ckc

arguments
    data cell
    options.BackgroundColor = 'w'; % Figure/axes background color
    options.CData = [];
    options.Colormap = @jet;
    options.EventWaveform cell = {};
    options.Figure = [];
    options.Labels = '';
    options.Layout = [];
    options.Position (1,4) double = [200 200 1000 640];
    options.PulseHeight double = 0.8 % Fraction from 0-1 for height of each line. (Can be set as scalar or vector for each cell element in data.)
    options.PulseLineWidth (1,1) double = 0.9 % 'LineWidth' parameter for each line tick
    options.Reference = [] % Reference signal(s)
    options.RatePercentileNormalizer (1,1) double {mustBeInRange(options.RatePercentileNormalizer,0,1)} = 0.8
    options.ReferenceColor (:,3) double {mustBeInRange(options.ReferenceColor, 0, 1)} = [0.7 0.7 0.7] % Color of lightest reference (or nRef x 3 array of colors for each reference line).
    options.ReferenceLabels = {}
    options.SampleRate (1,1) double = 4000 % Default from TMSi recordings
    options.SmoothKernelWidth (1,1) double = 0.25 % Smoothing kernel width (s)
    options.Title {mustBeTextScalar} = ''
    options.UserData (1,1) struct = struct
    options.YLabel {mustBeTextScalar} = '(pulse train)';
    options.XLim = [];
end

N = numel(data);

if isempty(options.Title)
    namestr = {'Pulse Trains'};
else
    if iscell(options.Title)
        namestr = options.Title;
    else
        if isa(options.Title, 'char')
            namestr = {options.Title};
        elseif isa(options.Title, 'string')
            namestr = options.Title;
        else
            error("Invalid type (%s) for parameter 'Title' (should be char, string, or cell).", class(options.Title));
        end
    end
end

if ~isempty(options.EventWaveform)
    if numel(options.EventWaveform) ~= numel(data)
        error("Must include a cell element (even an empty cell) 1:1 for each cell of raster data.");
    end
end

% Set cdata for each pulse train
if isempty(options.CData)
    cdata = feval(options.Colormap, N);
else
    cdata = options.CData;
end

if isempty(options.Layout) && isempty(options.Figure)
    fig = figure(...
        'Name', sprintf('Raster: %s', namestr{1}), ...
        'Color', options.BackgroundColor, ...
        'Position', options.Position, ...
        'UserData', options.UserData);
    L = get_layout(fig, N, options);
elseif isempty(options.Figure)
    L = options.Layout;
    fig = L.Parent;
elseif isempty(options.Layout)
    fig = options.Figure;
    L = get_layout(fig, N, options);
else
    fig = options.Figure;
    L = options.Layout;
end
raster_ax = nexttile(L, 1, [N 4]);
set(raster_ax, ...
    'NextPlot', 'add', ...
    'FontName', 'Tahoma', ...
    'Color', options.BackgroundColor, ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'YTick', 1:N, ...
    'YLim', [0.5, N+0.6], ...
    'YTickLabelRotation', 30, ...
    'Box', 'on');
xlabel(L, 'Time (s)', 'FontName', 'Tahoma', 'Color', 'k');
title(L, namestr{:}, 'FontName', 'Tahoma', 'Color', 'k');
ylabel(raster_ax, options.YLabel, ...
        'FontName', 'Tahoma', ...
        'FontSize', 10, ...
        'Color', [0.75 0.75 0.75]);
if ~isempty(options.Labels)
    raster_ax.YTickLabels = options.Labels;
end
if isnumeric(raster_ax.YTickLabels(1)) || ischar(raster_ax.YTickLabels(1))
    tmp = strings(size(raster_ax.YTickLabels));
    for ii = 1:numel(raster_ax.YTickLabels)
        if isnumeric(raster_ax.YTickLabels(ii))
            tmp(ii) = sprintf('\\color[rgb]{%f,%f,%f}%g', cdata(ii,1), cdata(ii,2), cdata(ii,3), raster_ax.YTickLabels(ii));
        else
            tmp(ii) = sprintf('\\color[rgb]{%f,%f,%f}%s', cdata(ii,1), cdata(ii,2), cdata(ii,3), raster_ax.YTickLabels(ii));
        end
    end
    raster_ax.YTickLabels = tmp;
else
    for ii = 1:numel(raster_ax.YTickLabels)
        if iscell(raster_ax.YTickLabels(ii))
            raster_ax.YTickLabels{ii} = sprintf('\\color[rgb]{%f,%f,%f}%s', cdata(ii,1), cdata(ii,2), cdata(ii,3), raster_ax.YTickLabels{ii});
        else
            raster_ax.YTickLabels(ii) = sprintf('\\color[rgb]{%f,%f,%f}%s', cdata(ii,1), cdata(ii,2), cdata(ii,3), raster_ax.YTickLabels(ii));
        end
    end
end

fs = options.SampleRate; % Sample rate

[nRef,nRefSamples] = size(options.Reference);
if nRef > nRefSamples
    ref_in = options.Reference';
    [nRef, nRefSamples] = size(ref_in);
else
    ref_in = options.Reference;
end
if nRef > 0
    tRef = (0:(nRefSamples-1)) ./ fs;
    
    sRef = N / nRef; % Reference signals distributed along y-axis.
    lRef = string(options.ReferenceLabels);
    if (nRef > 1)
        if numel(lRef) == 1
            lRef = repmat(lRef, nRef, 1);
        end
        if size(options.ReferenceColor,1) == 1
            refColor = repmat(options.ReferenceColor, nRef, 1) ./ ((1:nRef)');
        else
            refColor = options.ReferenceColor;
        end
    else
        refColor = options.ReferenceColor;
    end
    if isempty(lRef)
        XL = [tRef(1), tRef(end)]; % No need to "make room" for legend
    else
        XL = [tRef(1), tRef(end)+(tRef(end)-tRef(1))*0.2];
    end
    if isempty(options.EventWaveform)
        ref_ax = nexttile(L, N*4+1, [1 4]);
    else
        ref_ax = nexttile(L, N*5+1, [1 4]);
    end
    set(ref_ax, ...
        'NextPlot', 'add', ...
        'FontName', 'Tahoma', ...
        'Color', options.BackgroundColor, ...
        'XColor', 'k', ...
        'YColor', 'k', ...
        'YTick', 1:N, ...
        'YLim', [0.5, N+0.6], ...
        'YTickLabelRotation', 30, ...
        'Box', 'on');
    
    for ii = 1:nRef
        ref = ref_in(ii,:) - min(ref_in(ii,:));
        ref = ref ./ max(ref);
        yRef = sRef .* ref + sRef*(ii-1) + 0.5;
        cRef = refColor(ii,:);
        hRef = plot(ref_ax, tRef, yRef, ...
            'LineWidth', max(5 - ii, 1), ...
            'Color', cRef);
        if ~isempty(lRef)
            hRef.DisplayName = lRef(ii);
        end
    end
else
    XL = [-inf, inf];
    lRef = string({});
end

if ~isempty(options.EventWaveform)
    for ii = 1:numel(options.EventWaveform)
        event_ax = nexttile(L, N*5 - (ii-1)*5,[1 1]);
        set(event_ax, 'NextPlot','add','XColor','none','YColor','none');
        if isvector(options.EventWaveform{ii})
            plot(event_ax, options.EventWaveform{ii}, 'LineWidth', 3, 'Color', cdata(ii,:));
            title(event_ax, sprintf('Template-%d', ii), 'FontName','Tahoma','Color',cdata(ii,:));
        else
            mu = mean(options.EventWaveform{ii},2);
            t_mu = 1:numel(mu);
            sd = std(options.EventWaveform{ii},[],2);
            n_waves = size(options.EventWaveform{ii},2);
            errorbar(event_ax, t_mu, mu, mu+sd ./ (2*sqrt(n_waves)), mu-sd ./ (2*sqrt(n_waves)), 'Color', cdata(ii,:));
            title(event_ax, sprintf('Template-%d (N = %d)', ii, n_waves), 'FontName','Tahoma','Color',cdata(ii,:));
        end
    end
end

% Use this to determine the x-bounds (if no ref. signal already).
xl = [inf, -inf];

% Set the pulse height for the individual pulse trains.
if ~iscell(options.PulseHeight(1))
    hh = options.PulseHeight ./ 2;
    if isscalar(hh)
        hh = repmat(hh, N, 1);
    end
else
    hh = options.PulseHeight;
end

for ii = 1:N
    n = numel(data{ii});
    if n == 0
        xts = nan(3,1);
        n = 1; % "Pretend" so that we can put a hidden graphic object as placeholder.
    else
        ts = sort(data{ii} ./ fs, 'ascend');
        ts = reshape(ts, 1, n);
        xl = [min(xl(1), ts(1)), max(xl(2), ts(end))];
        xts = [repmat(ts, 2, 1); ...
               nan(1, n)]; % "Hack" to use NaN for separating the lines.
    end
    if iscell(hh(ii))
        hhh = hh{ii} ./ 2;
    else
        hhh = hh(ii,:);
    end
    yts = [ones(1,n).*(ii-hhh); ...
           ones(1,n).*(ii+hhh); ...
            nan(1,n)]; % "Hack" to use NaN for separating the lines.
    xt = 0:(1/fs):(max(ts)+0.1);
    xyt = zeros(size(xt));
    xyt(data{ii}) = 1;
    hhhy = mean(hhh(~isnan(hhh)));
    xyt = conv(xyt, hann(round(options.SmoothKernelWidth*fs)), 'same');
    xyt_s = sort(xyt, 'ascend');
    i_norm = round(numel(xyt_s)*options.RatePercentileNormalizer);
    xyt = (xyt./xyt_s(i_norm)).*(1-2*hhhy) + (ii+hhhy);

    hl = line(raster_ax, xts(:), yts(:), ...
        'Color', cdata(ii,:), ...
        'LineWidth', options.PulseLineWidth);
    hl.Annotation.LegendInformation.IconDisplayStyle = 'off';

    hl = line(raster_ax, xt, xyt, ...
        'Color', cdata(ii,:)*0.75, ...
        'LineStyle', '-', ...
        'LineWidth', options.PulseLineWidth*1.5);
    hl.Annotation.LegendInformation.IconDisplayStyle = 'off';
    if ii == N
        raster_ax.YLim(2) = max(xyt);
    end
end

if isempty(lRef)
    XL = [max(xl(1), XL(1)), min(xl(2), XL(2))];
else
    XL = [max(xl(1), XL(1)), min(xl(2)+diff(xl)*0.2, XL(2))];
    legend(raster_ax, ...
        'FontName', 'Tahoma', ...
        'TextColor', 'black', ...
        'FontSize', 8, ...
        'Location', 'NorthEast');
end
if exist('ref_ax','var')~=0
    linkaxes([raster_ax, ref_ax], 'x');
end
if isempty(options.XLim)
    raster_ax.XLim = XL;
else
    raster_ax.XLim = options.XLim;
end

    function L = get_layout(fig, N, options)
        if isempty(options.Reference) && isempty(options.EventWaveform)
            L = tiledlayout(fig, N, 4);
        elseif isempty(options.Reference)
            L = tiledlayout(fig, N, 5);
        elseif isempty(options.EventWaveform)
            L = tiledlayout(fig, N+1, 4);
        else
            L = tiledlayout(fig, N+1, 5);
        end
    end

end
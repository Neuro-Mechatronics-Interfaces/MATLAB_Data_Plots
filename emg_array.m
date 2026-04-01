function fig = emg_array(samples, options)
%EMG_ARRAY Plot emg array given samples as they come from TMSi SDK.
%
% Syntax:
%   fig = plot.emg_array(samples, 'Name', value, ...);
%
% Inputs:
%   samples - nChannels x nSamples array of sample data
%   options:
%
arguments
    samples (:, :) double
    options.ApplyFilter (1,1) logical = true;
    options.ApplyCommonReference (1,1) logical = true;
    options.ApplyZCAWhitening (1,1) logical = true; 
    options.CommonReferenceThreshold (1,1) double = 3.5;
    options.BadChannels (1,:) double = []; % Indices for bad channels (e.g. 1-64)
    options.CData = []; 
    options.FilterCutoff (1,:) double = 300; % Hz
    options.DifferentiatorWindow (1,1) double {mustBePositive, mustBeInteger} = 9;
    options.FilterType {mustBeMember(options.FilterType,{'IIR-HPF', 'FIR-Diff'})} = 'FIR-Diff';
    options.GridLayout (1,2) double = [8, 8]; % Number of rows, number of columns
    options.NumSubGrids (1,1) double {mustBePositive, mustBeInteger} = 2;
    options.RMSScalar (1,1) double = 6.5;
    options.RefTitle {mustBeTextScalar} = '';
    options.Reference (:,:) double = [];
    options.SampleRate = 2000; % Samples per second
    options.Subtitle {mustBeTextScalar} = '';
    options.ReferenceSamples = [];
    options.ReferenceTimes = [];
    options.EpochStartTime = nan;
    options.EpochStopTime = nan; 
    options.Title {mustBeTextScalar} = '';
    options.TLim = [-inf, inf]; % Time limits (seconds)
    options.XUnits = 'sec';
    options.XScaleGain (1,1) double = 1; 
    options.YUnits = 'μV'; 
    options.YScaleGain (1,1) double = 1; 
    options.YOffset = 20; % Offset between vertical traces
    options.TikhonovEpsilon (1,1) double = 5;
    options.LineOptions cell = {};
    options.SampleTimes = []; 
    options.ShowYScale (1,1) logical = true;
    options.ShowXScale (1,1) logical = true;
end

i_data = setdiff(1:size(samples,1),options.BadChannels);
samples(options.BadChannels,:) = 0;

if options.ApplyFilter
    switch options.FilterType
        case 'IIR-HPF'
            if isscalar(options.FilterCutoff)
                [b,a] = butter(4,options.FilterCutoff./(options.SampleRate/2),'high');
            else
                [b,a] = butter(4,options.FilterCutoff./(options.SampleRate/2),'bandpass');
            end
            fdata = filtfilt(b,a,samples')';
        case 'FIR-Diff'
            [~,b] = sgolay(2,options.DifferentiatorWindow);
            b = b(:,2); a = 1; 
            fdata = filter(b,a,samples,[],2);
            fdata(:,1:options.DifferentiatorWindow) = 0; % Zero out the filter transient from startup
            fdata = sgolayfilt(fdata,2,options.DifferentiatorWindow,[],2); 
    end

    if options.ApplyCommonReference
        fdata = doCMR(fdata, options.NumSubGrids, options.CommonReferenceThreshold);
    end
    
    if options.ApplyZCAWhitening
        fdata = zca_whiten_emg(fdata, ...
            'TikhonovEpsilon', options.TikhonovEpsilon); 
    end
else
    fdata = samples(i_data,:);
end

if isempty(options.SampleTimes)
    t_data = (0:(size(fdata,2)-1))./options.SampleRate;
else
    assert(numel(options.SampleTimes)==size(fdata,2),"Must have one time element for each column in sample data.");
    t_data = options.SampleTimes;
end
t_mask = (t_data >= options.TLim(1)) & (t_data < options.TLim(2));
t_data = t_data(1,t_mask);
t_offset = t_data(end)-t_data(1);

y_offset = repmat(options.YOffset.*(0:(options.GridLayout(1)-1))',1,options.GridLayout(2));
y_offset = y_offset(i_data)';

x_offset = repmat((0:(options.GridLayout(2)-1)).*t_offset.*1.05,options.GridLayout(1),1);
x_offset = x_offset(i_data)';
t_plot = repmat(t_data,numel(i_data),1) + x_offset;

y_data = fdata(:,t_mask);
y_plot = y_data+y_offset;
if isempty(options.CData)
    if options.NumSubGrids > 1
        grid_k = round(options.GridLayout(1)*options.GridLayout(2)/options.NumSubGrids);
        c_data = repmat([winter(grid_k); spring(grid_k)],ceil(options.NumSubGrids/2),1); 
    else
        c_data = turbo(size(samples,1));
        c_data = c_data(i_data,:);
    end
else
    c_data = options.CData;
end

fig = figure(...
    'Name','EMG Array Snippets', ...
    'Color','w', ...
    'Units', 'inches', ...
    'Position',[1 1 6.25 5.75]); 
if ~isempty(options.Reference)
    L = tiledlayout(fig, 4, 1);
else
    L = tiledlayout(fig, 3, 1);
end
ax = nexttile(L, 1, [3 1]);
x_max = max(t_plot(:));
x_min = min(t_plot(:));
x_min = x_min-.05.*(x_max-x_min);
% y_max = max(y_plot(:));
% y_min = min(y_plot(:));
y_min = -options.YOffset;
y_max = options.GridLayout(1)*options.YOffset;
% y_min = y_min -.05.*(y_max-y_min);
set(ax,'NextPlot','add',...
    'FontName','Tahoma',...
    'Clipping','off',...
    'ColorOrder',c_data, ...
    'XLim', [x_min, x_max], ...
    'XColor','none',...
    'YLim', [y_min, y_max], ...
    'YColor','none');
if ~isempty(options.ReferenceSamples) || ~isempty(options.ReferenceTimes)
    refSamples = [nan(numel(options.ReferenceTimes),1);options.ReferenceSamples];
    for ii = 1:numel(options.ReferenceTimes)
        [~,refSamples(ii)] = min(abs(t_data - options.ReferenceTimes(ii)));  
    end
    refSamples = unique(refSamples);
    refSamples = sort(refSamples,'ascend'); 
    plot(ax, t_plot', y_plot', ...
        'LineWidth', 2, 'MarkerSize', 24,  ...
        'Marker', '|', 'MarkerFaceColor', 'k',  ...
        'MarkerIndices', refSamples, ...
        'MarkerEdgeColor', 'k', options.LineOptions{:}); 
else
    plot(ax, t_plot',y_plot', options.LineOptions{:});
end
if ~isnan(options.EpochStartTime)
    [~,epochStartSample] = min(abs(t_data - options.EpochStartTime)); 
    plot(ax, t_plot', y_plot', ...
        'LineWidth', 1, 'LineStyle', 'none', 'MarkerSize', 18,  ...
        'Marker', '|', 'MarkerFaceColor', 'b',  ...
        'MarkerIndices', epochStartSample, ...
        'MarkerEdgeColor', 'b','DisplayName','START'); 
end
if ~isnan(options.EpochStopTime)
    [~,epochStopSample] = min(abs(t_data - options.EpochStopTime)); 
    plot(ax, t_plot', y_plot', ...
        'LineWidth', 1, 'LineStyle', 'none', 'MarkerSize', 18,  ...
        'Marker', '|', 'MarkerFaceColor', 'r',  ...
        'MarkerIndices', epochStopSample, ...
        'MarkerEdgeColor', 'r','DisplayName','STOP'); 
end

if options.ShowXScale
    t_scale = t_offset/2; 
    t_scale = round(t_scale*options.XScaleGain,1)/options.XScaleGain; 
    line(ax, [x_min, x_min + t_scale], [y_min, y_min], ...
        'Color','k','LineWidth',1.25);
    text(ax, x_min + t_scale, y_min, sprintf('%3.1f (%s)', t_scale*options.XScaleGain, options.XUnits), ...
        'FontName','Tahoma', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','top');
end

if options.ShowYScale
    if ~isnan(options.YOffset)
        r = round(options.YOffset/2,1);
    else
        r = rms(y_data(:)) * options.RMSScalar;
    end
    y_scale = round(r*options.YScaleGain)/options.YScaleGain;
    line(ax, [x_min, x_min], [y_min, y_min + y_scale], 'Color','k','LineWidth',1.25);
    text(ax, x_min, y_min+y_scale, sprintf('%.1f (%s)', y_scale*options.YScaleGain, options.YUnits), ...
        'FontName','Tahoma', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom');
end

if strlength(options.Title) > 0
    title(ax, options.Title, 'FontName','Tahoma','Color','k','FontWeight','bold');
end

if strlength(options.Subtitle) > 0
    subtitle(ax,options.Subtitle,'FontName','Tahoma','Color',[0.65 0.65 0.65]);
end

if ~isempty(options.Reference)
    ref_ax = nexttile(L, 4, [1 1]);
    c_data_ref = repmat(linspace(0,0.65,size(options.Reference,1))',1,3);
    set(ref_ax, 'NextPlot','add','FontName','Tahoma','ColorOrder', c_data_ref);
    plot(ref_ax, t_data, options.Reference(:,t_mask));
    if strlength(options.RefTitle) > 0
        title(ref_ax, options.RefTitle, 'FontName','Tahoma','Color',[0.35 0.35 0.35]);
    else
        title(ref_ax, '(Reference)', 'FontName','Tahoma','Color',[0.35 0.35 0.35]);
    end
    xlabel(ref_ax, 'Time (s)', 'FontName','Tahoma','Color','k'); 
end
fig.UserData = options;
end

function out = doCMR(in, numSubGrids, commonReferenceThreshold)
out = in; 
k = round(size(in,1)/numSubGrids);
vec = 1:k;
delta = std(in,[],2); 
sigma = std(in(:)); 
th = sigma*commonReferenceThreshold;
i_bad = find(delta>th);
i_good = setdiff(1:size(in,1),i_bad);
out(i_bad,:) = 0;
for ii = 1:k
    o = (ii-1)*k;
    sub = vec + o;
    smask = ismember(i_good,sub);
    out(smask,:) = in(smask,:) - mean(in(smask,:),1);
end
end


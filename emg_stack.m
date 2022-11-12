function fig = emg_stack(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, CHANNEL, varargin)
%EMG_STACK Create stack of EMG individual trial traces.
%
% Syntax:
%   fig = plot.emg_stack(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, CHANNEL);
%   fig = plot.emg_averages(___, 'Name', value,...);
%
% Inputs:
%   SUBJ  - Subject name (e.g. 'Frank' or "Frank")
%   YYYY  - Year (numeric or string, e.g. 2021 or "2021" or '2021')
%   MM    - Month (numeric or string, e.g. 11 or "11" or '11')
%   DD    - Day   (numeric or string, e.g. 18 or "18" or '18')
%   ARRAY - "A" or "B" (or 'A' or 'B') or ["A", "B"] or {'A', 'B'}
%   BLOCK - Recording parameter key (block; numeric or string, e.g. 0 or "0")
%   CHANNEL - Data channel
%   varargin - (Optional) 'Name', value parameter field/value pairs. See
%                       `pars` struct for specific parameter names. A few
%                       key parameters:
%       - 'Sync_Bit': (Default: 9) -- this determines which BIT is used
%           from the uint16 data on sync channels (TRIGGERS or STATUS).
%
%   If BLOCK is an array, then the returned output will be an array of
%   figure handles with a row for each element of BLOCK. If ARRAY is a
%   (MATLAB) array, then the returned output will have a column for each
%   element of ARRAY. For recording blocks that do not have any data or
%   block associated to a particular ARRAY, the corresponding figure handle
%   object will be an empty `gobjects` element.
%
% Output:
%   fig   - Figure handle. If no output is requested, figures are
%           automatically saved and deleted using `default.savefig`. The
%           output saved results are in `generated_data` and follow the
%           naming convention provided by the input arguments.
%
% See also: Contents, plot.emg_averages, plot.emg_waterfall

% Handle parsing of `pars`
if (numel(varargin) > 0) && isstruct(varargin{1})
    pars = varargin{1};
    varargin(1) = [];
else
    pars = plot.parameters('emg_stack');
end

pars = utils.parse_parameters(pars, varargin{:});
if ~isstruct(pars.Filtering)
    pars.Filtering = utils.get_default_filtering_pars(pars.Acquisition_Type, pars.EMG_Type, pars.Filtering);
end

if (numel(BLOCK) > 1) || (numel(ARRAY) > 1)
    % Can only reach here if char arguments were given instead of x directly
    fig = cell(numel(BLOCK), numel(ARRAY));
    for iB = 1:numel(BLOCK)
        for iA = 1:numel(ARRAY)
            if nargout > 0
                fig{iB, iA} = plot.emg_stack(SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), CHANNEL, pars);
            else
                plot.emg_stack(SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), CHANNEL, pars);
            end
        end
    end
    return;
end

% Load data in
if ~isempty(pars.Data)
    x = pars.Data;
else
    x = io.load_tmsi(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, pars.File_Type);
end

if numel(CHANNEL) > 1
    fig = gobjects(numel(CHANNEL),1);
    pars.Data = x;
    for iC = 1:numel(CHANNEL)
        if nargout > 0
            fig(iC) = plot.emg_stack(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, CHANNEL(iC), pars);
        else
            plot.emg_stack(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, CHANNEL(iC), pars);
        end
    end
    return;
end



tank = sprintf('%s_%04d_%02d_%02d', SUBJ, YYYY, MM, DD); % data "tank"
block = sprintf('%s_%s_%d', tank, ARRAY, BLOCK); % experimental "block" (recording within tank)
gen_data_folder = fullfile(pars.Output_Root, SUBJ, tank, num2str(BLOCK));

% Get trigger channel
channels = horzcat(x.channels{:});
if isnan(pars.Sync_Bit)
    sync_data_in_file = fullfile(gen_data_folder, sprintf('%s_sync.mat', x.name));
    if exist(sync_data_in_file, 'file')==0
        error('Plot:Sync', 'No sync data file (<strong>%s</strong>): must specify sync bit as non-NaN value!', sync_data_in_file);
    end
    in = load(sync_data_in_file, 'onset', 'offset', 'sync_data');
    stops = in.onset;
    trigs = in.offset;
else
    [stops, trigs, ~] = utils.parse_bit_sync(x, pars.Sync_Bit, gen_data_folder, pars.Inverted_Logic, pars.Trigger_Channel);
end

% Check if there are any trigger events.
if (numel(trigs) < 1) || (numel(stops) < 1)
    warning("Empty sync vector (trigs): check if TTL on TRIGGERS channel was present/parsed using correct bit.");
    fig = gobjects(1);
    return;
end

% Check that the first trigger onset is before the first "stop" onset.
if stops(1) < trigs(1)
    if numel(stops) > numel(trigs)
        stops(1) = [];
    else
        tmp = stops;
        stops = trigs;
        trigs = tmp;
    end
end

% Potentially select subset of trials to average
if ~isnan(pars.N_Trials)
    if pars.N_Trials(1) > numel(trigs)
        pars.N_Trials(1) = numel(trigs);
    end
    if numel(pars.N_Trials) == 1
        trials = 1:pars.N_Trials;
    else
        trials = reshape(pars.N_Trials,1,numel(pars.N_Trials));
    end 
    trigs = trigs(trials);
end

if strcmpi(pars.EMG_Type, 'Bipolar')
    iBip = contains({channels.alternative_name}, 'BIP')' & (sum(abs(x.samples-mean(x.samples, 2)),2) > eps);
    if sum(iBip) == 0
        fprintf(1, 'No BIPOLAR channels for recording: <strong>%s</strong>\n', block);
        return;
    else
        channels = channels(iBip);
    end
    data = x.samples(iBip, :)';
else
    iUnip = contains({channels.alternative_name}, 'UNI');
    channels = channels(iUnip);
    data = x.samples(iUnip, :)';
end

if pars.EMG_Filters_Applied==true
    z = data;
else
    if pars.Blank_Stim
        pars.Filtering.Apply_Stim_Blanking = true;
    end
    % Trigs is returned because the filtering function can exclude
    % out-of-bounds trigger sample indices based on stim-artifact-rejection
    % sample epoch width.
    [z, ~, pars.Filtering, trigs] = utils.apply_emg_filters(data, pars.Filtering, x.sample_rate, trigs, stops);
end

n_pre = -1 * round(pars.T(1) * 1e-3 * x.sample_rate); % Convert to seconds, then samples
n_post = round(pars.T(2) * 1e-3 * x.sample_rate);  % Convert to seconds, then samples


z = z(CHANNEL,:);
name = strrep(channels(CHANNEL).alternative_name, ' ', '');
try
    stim = utils.get_tmsi_stim_data(SUBJ, YYYY, MM, DD, ARRAY, BLOCK);
    chan_name = sprintf('%s (%s)', strrep(stim.map.Muscles.(name), '_', '-'), name);
catch
    chan_name = name;
end

% Trigs is returned because the averaging function can exclude
% out-of-bounds trigger sample indices based on snippet
% sampling epoch width (samples pre- and post-TTL marker).
[~, X, trigs] = math.triggered_average(trigs, z, n_pre, n_post, false, false, false);

% Align peaks (some trials may need assistance with stim event alignments
% due to manual searching of events)
if pars.Align_Peaks
    X = math.align_peaks(X, n_pre);
end

if pars.Subtract_Mean || pars.Filtering.Subtract_Cross_Trial_Mean
    pars.Filtering.Subtract_Cross_Trial_Mean = true;
    pars.Subtract_Mean = true;
    X = abs(X - mean(X,1)); 
end

N = size(X, 1);
t_sweep = (-n_pre:n_post)/x.sample_rate * 1e3; % Convert from samples to seconds, then milliseconds

% Generate figure
if isempty(pars.Axes)
    x = 60 + round(randn(1) * 5); % Screen position
    y = 60 + round(randn(1) * 5); % Screen position
    width = 700; % Width of figure
    height = 700; % Height of figure (by default in pixels)
    fig = figure('Position', [x y width height], 'Color', 'w');
    ax = axes(fig, 'NextPlot', 'add', ...
        'XColor', 'k', 'YColor', 'k', 'YDir', 'reverse', ...
        'LineWidth', 1.35, 'FontName', 'Tahoma', 'YLim', [-1 (N+2)], ...
        'XLim', [t_sweep(1), t_sweep(end)], 'Clipping', 'off');
else
    fig = pars.Axes.Parent;
    if ~isa(fig, 'matlab.ui.Figure')
        fig = fig.Parent; % In case axes is in a Layout or Panel or something
    end
    ax = pars.Axes;
    set(ax, 'NextPlot', 'add');
end

% Scale them so that most values should be < 1.
X = zscore(X, 0, 'all')./pars.Scale_Factor + (1:N)';

if isempty(pars.Series)
    pars.Series = 1:N;
    c = parula(N).*0.9;
else
    c = parula(numel(unique(pars.Series))).*0.9;
end
for ii = 1:N
    line(ax, t_sweep, X(ii,:), 'Color', c(pars.Series(ii),:), 'LineWidth', 1.5, 'DisplayName', sprintf('Trial-%d', ii));
end
if ~isempty(pars.XLim)
    xlim(ax, pars.XLim);
end
line(ax, ones(2,1) .* (ax.XLim(1) - 0.045 * diff(ax.XLim)), ...
    [ax.YLim(1) + 0.02 * diff(ax.YLim), ax.YLim(2) - 0.02 * diff(ax.YLim)], ...
    'MarkerIndices', 2, 'Marker', 'v', 'MarkerFaceColor', c(end,:), ...
    'Color', 'k', 'Tag', 'TrialOrder', 'LineWidth', 2.0, ...
    'MarkerSize', 16);
xline(ax, 0, 'r:', "Stim", 'LineWidth', 1.5);
% Generating figure title
str = utils.get_filtering_label_string(pars.Filtering);

if pars.Anonymize
    tmp = strsplit(block, '_');
    switch upper(string(SUBJ))
        case "MATS"
            block_a = strjoin(['NG', tmp(2:end)], '\\_');
        case "PULKIT"
            block_a = strjoin(['QH', tmp(2:end)], '\\_');
        case "CHAITANYA"
            block_a = strjoin(['DH', tmp(2:end)], '\\_');
        case "DOUG"
            block_a = strjoin(['EX', tmp(2:end)], '\\_');
        otherwise
            block_a = strjoin(tmp, '\\_');
    end
    title(ax, [char(block_a), ': ' str newline 'Chan:' chan_name ' (N = ' char(num2str(N))  ') Stacked'], ...
        'FontName', 'Tahoma', ...
        'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
else
    title(ax, [char(strrep(block, '_', '\_')), ': ' str newline 'Chan:' chan_name ' (N = ' char(num2str(N))  ') Stacked'], ...
        'FontName', 'Tahoma', ...
        'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
end

% Set figure view
xlabel('Time (ms)', pars.Font{:});
ylabel('Trial', pars.Font{:});

% If no figure then do not try to save regardless.
if isa(fig, 'matlab.graphics.GraphicsPlaceholder')
    return;
end

if nargout < 1
    tank = sprintf('%s_%04d_%02d_%02d', SUBJ, YYYY, MM, DD); % data "tank"
    block = sprintf('%s_%s_%d', tank, ARRAY, BLOCK); % experimental "block" (recording within tank)
    out_folder = fullfile(pars.Output_Root, SUBJ, tank, 'figures', 'Stack', pars.Filtering.Name);
    if exist(out_folder, 'dir') == 0
        try %#ok<TRYNC>
            mkdir(out_folder);
        end
    end
    out_name = fullfile(out_folder, sprintf('%s_%d_%d_%s-Ch%02d', block, round(pars.T(1)), round(pars.T(2)), pars.EMG_Type, CHANNEL));
    default.savefig(fig, out_name, sprintf("Stacked_EMG"), true);
    
    out_folder_2 = fullfile(pars.Output_Root, SUBJ, tank, num2str(BLOCK));
    if exist(out_folder_2, 'dir') == 0
        try %#ok<TRYNC>
            mkdir(out_folder_2);
        end
    end
    default.savefig(fig, fullfile(out_folder_2, block), sprintf('%d_%d_Stacked_%s-%s_%s-Ch%02d', round(pars.T(1)), round(pars.T(2)), pars.EMG_Type, pars.Filtering.Name, pars.EMG_Type, CHANNEL), false); 
end


end



function p = emg_waterfall(x, varargin)
%EMG_WATERFALL Create waterfall plots to show individual stimulus trial responses.
%
% Syntax:
%   p = plot.emg_waterfall(x, varargin);
%
% Inputs:
%   x        - double array containing the averaged segmented data (for the EMG grid)
%   chan     - specific channel to extract trigger aligned segmented data from
%   triggers - vector of trigger events to align data
%
%
% Output:
%   p - waterfall patch object handle
%
% Note: use this function with following steps: (Outdated steps 1-25-22)
%       x = load_tmsi_raw('Forrest',2022,1,25,'A',5);
%       output_path = create_raw_data_path('Forrest',2022,1,25,'A',5);
%       [onset, offset] = parse_bit_sync(x, 10, output_path, true, 'TRIGGERS');
%       [A,X, triggers] = triggered_average(onset, x.samples(1,:), 200, 800, true);
%
% Example:
%      p = plot_emg_waterfall('Forrest',2022,2,8,'A',68,'Data_Channel',14,'Sync_Bit',9,'N_Trials',99,...
%      'Trigger_Channel','TRIGGERS','Filtering',get_default_filtering_pars("TMSi","Array","Raw",...
%      'Apply_Virtual_Reference',false,"Apply_HPF",true,"HPF_Cutoff_Frequency",30))

% Define default parameters and handle `varargin`

if isa(x,'TMSiSAGA.Data') || isa(x, 'struct')
    str = split(x.name,'_');
    SUBJ = str{1};
    YYYY = str2double(str{2});
    MM = str2double(str{3});
    DD = str2double(str{4});
    ARRAY = str{5};
    BLOCK = str2double(str{6});
    pars = plot.parameters('emg_waterfall');
    pars.Sample_Rate = x.sample_rate;
    data_in = x.samples;
elseif isa(x, 'char') || isa(x, 'string')
    SUBJ = x;
    YYYY = varargin{1};
    MM = varargin{2};
    DD = varargin{3};
    ARRAY = varargin{4};
    BLOCK = varargin{5};
    varargin(1:5) = [];
    if (numel(varargin) > 0) && isstruct(varargin{1})
        pars = varargin{1};
        varargin(1) = [];
    else
        pars = plot.parameters('emg_waterfall'); 
    end
end

% Handle parsing of `pars`
pars = utils.parse_parameters(pars, varargin{:});
if ~isstruct(pars.Filtering)
    pars.Filtering = get_default_filtering_pars(pars.Acquisition_Type, pars.EMG_Type, pars.Filtering);
end
if pars.Verbose
    pars.Filtering.Verbose = pars.Verbose;
end

if (numel(BLOCK) > 1) || (numel(ARRAY) > 1)
    % Can only reach here if char arguments were given instead of x directly
    p = gobjects(numel(BLOCK), numel(ARRAY));
    for iB = 1:numel(BLOCK)
        for iA = 1:numel(ARRAY)
            if nargout > 0
                p(iB, iA) = plot.emg_waterfall(SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), pars);
            else
                plot.emg_waterfall(SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), pars);
            end
        end
    end
    return;
end

if ~(isa(x,'TMSiSAGA.Data') || isa(x, 'struct')) && isempty(pars.Trigger_Data)
    x = io.load_data(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, pars.File_Type, pars.Input_Root);
    data_in = x.samples;
    pars.Sample_Rate = x.sample_rate;
    if isempty(x)
        return;
    elseif numel(x) > 1
        warning('Should only return exactly one or no data matches. Check load_tmsi_raw. Block skipped.');
        return;
    end
end

if isempty(pars.Trigger_Data)
    tank = sprintf('%s_%04d_%02d_%02d', SUBJ, YYYY, MM, DD); % data "tank"
    block = sprintf('%s_%s_%d', tank, ARRAY, BLOCK); % experimental "block" (recording within tank)
    gen_data_folder = fullfile(pars.Output_Root, SUBJ, tank, num2str(BLOCK));
    sync_data_in_file = fullfile(gen_data_folder, sprintf('%s_sync.mat', x.name));
else
    block = x.name;
end

if exist(sync_data_in_file, 'file')==0
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
end

% Check if there are any trigger events.
if (numel(trigs) < 1) || (numel(stops) < 1)
    warning("Empty sync vector (trigs): check if TTL on TRIGGERS channel was present/parsed using correct bit.");
    p = gobjects(1);
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
    data = data_in(iBip, :)';
else
    iUnip = contains({channels.alternative_name}, 'UNI');
    data = data_in(iUnip, :)';
end

if pars.EMG_Filters_Applied==true
    if isempty(pars.Filtered_Data)
        z = data;
    else
        z = pars.Filtered_Data;
    end
else
    if isempty(pars.Filtered_Data)
        % Trigs is returned because the filtering function can exclude
        % out-of-bounds trigger sample indices based on stim-artifact-rejection
        % sample epoch width.
        [z, ~, pars.Filtering, trigs] = utils.apply_emg_filters(data, pars.Filtering, x.sample_rate, trigs, stops);
    else
        z = pars.Filtered_Data;
    end
end

[z, ~, pars.Filtering] = utils.apply_emg_filters(x, pars.Filtering, x.sample_rate, trigs, stops);

n_pre = -1 * round(pars.T(1) * 1e-3 * x.sample_rate); % Convert to seconds, then samples
n_post = round(pars.T(2) * 1e-3 * x.sample_rate);  % Convert to seconds, then samples

if isnan(pars.Data_Channel) || pars.Data_Channel==0
    % Select the second channel
    z = z(2,:);
    pars.Data_Channel = 2;
    chan_name = channels(64).alternative_name;
else
    z = z(pars.Data_Channel,:);
    chan_name = channels(pars.Data_Channel).alternative_name;
end

if isempty(pars.Trigger_Data)

    if isnan(pars.Data_Channel) || pars.Data_Channel==0
        % Trigs is returned because the averaging function can exclude
        % out-of-bounds trigger sample indices based on snippet
        % sampling epoch width (samples pre- and post-TTL marker).
        [~, X, trigs] = math.triggered_average(trigs, z, n_pre, n_post, false, false, false);
    else
        if strcmpi(pars.EMG_Type, 'Bipolar') && (isnan(pars.Data_Channel) || (pars.Data_Channel == 0))
            for iCh = 1:sum(iBip)
                [~, X, trigs] = math.triggered_average(trigs, z(iCh, :), n_pre, n_post, false, false, false);

            end
        else
            [~, X, trigs] = math.triggered_average(trigs, z, n_pre, n_post, false, false, false);
        end
    end
else
    X = pars.Trigger_Data;
end

N = size(X,1);

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

M = size(X,2);
t_sweep = (-n_pre:n_post)/pars.Sample_Rate * 1e3; % Convert from samples to seconds, then milliseconds
times = repmat(t_sweep,N,1);
trials = repmat(trials,M,1)';

% Generate figure
if isempty(pars.Axes)
    x = 10; % Screen position
    y = 30; % Screen position
    width = 700; % Width of figure
    height = 700; % Height of figure (by default in pixels)
    fig = figure( 'Name','Waterfall Plot',...
             'Units','Normalized', ...
             'Position',[0.1 0.1 0.8 0.8],...
             'Color', 'w');
    ax = axes(fig, 'NextPlot', 'add', ...
        'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', ...
        'LineWidth', 1.35, 'View', pars.View);

else
    fig = pars.Axes.Parent;
    if ~isa(fig, 'matlab.ui.Figure')
        fig = fig.Parent; % In case axes is in a Layout or Panel or something
    end
    ax = pars.Axes;
    set(ax, 'NextPlot', 'add', 'View', pars.View, 'XDir', 'reverse');
end
p = waterfall(ax, trials, times, X);
set(p, 'LineWidth', 1.5, 'FaceAlpha', 0.05, 'EdgeAlpha', 0.25);
colorbar(ax);
axis tight;
if ~isempty(pars.X_Lim)
    xlim(ax, pars.X_Lim);
end
if ~isempty(pars.Y_Lim)
    ylim(ax, pars.Y_Lim);
end
if ~isempty(pars.Z_Lim)
    zlim(ax, pars.Z_Lim); 
end
if ~isempty(pars.C_Lim)
    ax.CLim = pars.C_Lim; 
end
% Generate stim event line
% TO-DO


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
    title([char(block_a), ': ' str newline 'Chan:' chan_name ' (N = ' char(num2str(N))  ') Stacked'], ...
        'FontName', 'Tahoma', ...
        'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
else
    title([char(strrep(block, '_', '\_')), ': ' str newline 'Chan:' chan_name ' (N = ' char(num2str(N))  ') Stacked'], ...
        'FontName', 'Tahoma', ...
        'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
end


% Set figure view
% [caz,cel] = view(ax, [1 0.001 6]); % vertical stack
colormap(pars.Colormap);
%v = [5 4 4];
%[caz,cel] = view(v);
% fc = fig.Children;
ax.XTick = unique( round(ax.XTick) );
xlabel('Trial', pars.Font{:});
ylabel('Time (ms)', pars.Font{:});
zlabel('uV', pars.Font{:});


% If no figure then do not try to save regardless.
if isa(fig, 'matlab.graphics.GraphicsPlaceholder')
    return;
end

if nargout < 1
    % Second directory saves plots to folders separated by block, that way
    % the "../Figures/Waterfall/.." directory doesn't get cluttered
    tank = sprintf('%s_%04d_%02d_%02d', SUBJ, YYYY, MM, DD); % data "tank"
    block = sprintf('%s_%s_%d', tank, ARRAY, BLOCK); % experimental "block" (recording within tank)
    out_folder = fullfile(pars.Output_Root, SUBJ, tank, 'figures', 'Waterfall', pars.Filtering.Name);
    if exist(out_folder, 'dir') == 0
        try %#ok<TRYNC>
            mkdir(out_folder);
        end
    end
    out_name = fullfile(out_folder, sprintf('%s_%d_%d_', block, round(pars.T(1)), round(pars.T(2))));
    %default.savefig(fig, out_name, sprintf("Waterfall_EMG_%s_%s", string(chan_name), pars.Filtering.Name), true);
    
    out_folder_2 = fullfile(pars.Output_Root, SUBJ, tank, num2str(BLOCK));
    if exist(out_folder_2, 'dir') == 0
        try %#ok<TRYNC>
            mkdir(out_folder_2);
        end
    end
    default.savefig(fig, fullfile(out_folder_2, block), sprintf('%d_%d_Waterfall_%s-%s', round(pars.T(1)), round(pars.T(2)), pars.EMG_Type, pars.Filtering.Name), false); 
end

end



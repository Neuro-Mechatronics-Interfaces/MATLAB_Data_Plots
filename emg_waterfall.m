function p = emg_waterfall(x, varargin)
%PLOT_EMG_WATERFALL Create waterfall plots to show individual stimulus trial responses.
%
% Syntax:
%   p = emg_waterfall(x, varargin);
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

if (numel(varargin) == 1) && isstruct(varargin{1})
    pars = varargin{1};
else
    % Define default parameters and handle `varargin`
    pars = struct;
    pars.Acquisition_Type = "TMSi";
    pars.Align_Peaks = false;
    pars.Axes = [];
    pars.Colormap = 'spring';
    pars.Data_Channel = nan;
    pars.EMG_Filters_Applied = false;
    pars.EMG_Type = "Array"; % Can be: "Array" | "Bipolar"
    pars.Figure_Title = 'Waterfall';
    pars.File_Type = ".mat"; % Can be: ".mat" | ".poly5"
    pars.Filtering = get_default_filtering_pars("TMSi","Array","Raw", ...
        'Apply_Virtual_Reference',false,"Apply_HPF",true,"HPF_Cutoff_Frequency",30); % Return default filtering struct
    pars.Font = {'FontName', 'Tahoma', 'FontSize', 18, 'Color', 'k'};
    pars.Inverted_Logic = false;
    pars.N_Individual_Max = 10; % Max. number of individual traces to superimpose
    pars.N_Trials = 10;
    [pars.Output_Root, pars.Input_Root] = parameters('generated_data_folder', 'raw_data_folder');
    pars.Plot_Stim_Period = true; % Plot stim artifact with red stem lines?
    pars.Sample_Rate = 4000;    % Sample rate from acquisition
    pars.Sync_Bit = nan;        % The bit address for STIM sync TTL signal on TRIGGERS channel of TMSi.
    pars.T = [-15, 80];     % Time for epochs (milliseconds)
    pars.T_RMS = [30, 60];  % Time epoch for computing RMS
    pars.Trigger_Data = [];
    pars.Trigger_Channel = 'TRIGGER';
    pars.View = [95, 65]; 

    % % % Color limits as well as axes limits % % %
    pars.C_Lim = []; % If empty, use auto-scale, otherwise, fixed scale
    pars.X_Lim = []; % If empty, use auto-scale, otherwise, fixed scale
    pars.Y_Lim = []; % If empty, use auto-scale, otherwise, fixed scale
    pars.Z_Lim = []; % If empty, use auto-scale, otherwise, fixed scale.

    if isa(x,'TMSiSAGA.Data') || isa(x, 'struct')
        str = split(x.name,'_');
        SUBJ = str{1};
        YYYY = str2double(str{2});
        MM = str2double(str{3});
        DD = str2double(str{4});
        ARRAY = str{5};
        BLOCK = str2double(str{6});
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
    end
end

% Handle parsing of `pars`
pars = utils.parse_parameters(pars, varargin{:});
if ~isstruct(pars.Filtering)
    pars.Filtering = get_default_filtering_pars(pars.Acquisition_Type, pars.EMG_Type, pars.Filtering);
end

if (numel(BLOCK) > 1) || (numel(ARRAY) > 1)
    % Can only reach here if char arguments were given instead of x directly
    fig = gobjects(numel(BLOCK), numel(ARRAY));
    for iB = 1:numel(BLOCK)
        for iA = 1:numel(ARRAY)
            if nargout > 0
                fig(iB, iA) = plot.emg_waterfall(SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), pars);
            else
                plot.emg_waterfall(SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), pars);
            end
        end
    end
    return;
end

if ~(isa(x,'TMSiSAGA.Data') || isa(x, 'struct')) && isempty(pars.Trigger_Data)
    x = io.load_tmsi(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, pars.File_Type, pars.Input_Root);
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
else
    block = x.name;
end

% Get trigger channel
channels = horzcat(x.channels{:});
if isempty(pars.Trigger_Data)
    if isnan(pars.Sync_Bit)
        sync_data_in_file = fullfile(gen_data_folder, sprintf('%s_sync.mat', x.name));
        if exist(sync_data_in_file, 'file')==0
            error('No sync data file (<strong>%s</strong>): must specify sync bit as non-NaN value!', sync_data_in_file);
        end
        in = load(sync_data_in_file, 'onset', 'offset', 'sync_data');
        stops = in.onset;
        trigs = in.offset;
    else
        [stops, trigs, ~] = parse_bit_sync(x, pars.Sync_Bit, gen_data_folder, pars.Inverted_Logic, pars.Trigger_Channel);
    end
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
    data = data_in;
end

if pars.EMG_Filters_Applied==true
    z = data;
else
    [z, ~, pars.Filtering] = utils.apply_emg_filters(data, pars.Filtering, x.sample_rate, trigs, stops);
end

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
        [~, X, ~] = math.triggered_average(trigs, z, n_pre, n_post, false, false, false);
    else
        if strcmpi(pars.EMG_Type, 'Bipolar')
            for iCh = 1:sum(iBip)
                [~, X, ~] = math.triggered_average(trigs, z(iCh, :), n_pre, n_post, false, false, false);

            end
        else
            [~, X, ~] = math.triggered_average(trigs, z, n_pre, n_post, false, false, false);
        end
    end
else
    X = pars.Trigger_Data;
end

% Trim number of channels to show
if ~isnan(pars.N_Trials)
    if pars.N_Trials > size(X,1)
        pars.N_Trials = size(X,1);
    end
    X = X(1:pars.N_Trials,:);
end

% Align peaks (some trials may need assistance with stim event alignments
% due to manual searching of events)
if pars.Align_Peaks
    X = math.align_peaks(X, n_pre);
end

[N, M] = size(X);
t_sweep = (-n_pre:n_post)/pars.Sample_Rate * 1e3; % Convert from samples to seconds, then milliseconds
times = repmat(t_sweep,N,1);
trials = repmat(1:N,M,1)';

% Generate figure
if isempty(pars.Axes)
    x = 10; % Screen position
    y = 30; % Screen position
    width = 700; % Width of figure
    height = 700; % Height of figure (by default in pixels)
    fig = figure('Position', [x y width height], 'Color', 'w');
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
set(p, 'LineWidth', 1.5, 'FaceAlpha', 0.05, 'EdgeAlpha', 0.75);
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
str = get_filtering_label_string(pars.Filtering);
title([char(strrep(block, '_', '\_')), ': ' str newline 'Chan:' chan_name ' (N = ' char(num2str(pars.N_Trials))  ') Waterfall'], ...
    'FontName', 'Tahoma', ...
    'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');

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
    tank = sprintf('%s_%04d_%02d_%02d', SUBJ, YYYY, MM, DD); % data "tank"
    block = sprintf('%s_%s_%d', tank, ARRAY, BLOCK); % experimental "block" (recording within tank)
    out_folder = fullfile(pars.Output_Root, SUBJ, tank, 'figures', 'Waterfall', pars.Filtering.Name);
    if exist(out_folder, 'dir') == 0
        try %#ok<TRYNC>
            mkdir(out_folder);
        end
    end
    out_name = fullfile(out_folder, sprintf('%s_%d_%d_Ch%d', block, round(pars.T(1)), round(pars.T(2)), pars.Data_Channel));
    default.savefig(fig, out_name, sprintf("Waterfall_EMG"), true);
    
    out_folder_2 = fullfile(pars.Output_Root, SUBJ, tank, num2str(BLOCK));
    if exist(out_folder_2, 'dir') == 0
        try %#ok<TRYNC>
            mkdir(out_folder_2);
        end
    end
    default.savefig(fig, fullfile(out_folder_2, block), sprintf('%d_%d_Waterfall_%s', round(pars.T(1)), round(pars.T(2)), pars.EMG_Type, pars.Filtering.Name), false); 
end


end



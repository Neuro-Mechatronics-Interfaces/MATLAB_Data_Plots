function fig = emg_averages__unipolar_array(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, varargin)
%EMG_AVERAGES__UNIPOLAR_ARRAY  EMG processing for HD-EMG TMSi Array
%
% Syntax:
%   fig = plot.emg_averages__unipolar_array(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, varargin);
%
% This should be accessed via `"Array"` EMG_Type parameter in
% `plot.emg_averages`. 
%
% See also: Contents, plot.emg_averages

% % % % SEE PARS STRUCT BELOW % % % %

if (numel(varargin) > 0) && isstruct(varargin{1})
    pars = varargin{1};
    varargin(1) = [];
else
    pars = struct;
    pars.Acquisition_Type = "TMSi";
    pars.EMG_Type = "Array"; % Can be: "Array" | "Bipolar"
    pars.End_Linear_Fit = inf; % Set to some value in milliseconds to say when linear fit should end
    pars.File_Type = ".mat"; % Can be: ".mat" | ".poly5"
    pars.Filtering = utils.get_default_filtering_pars(); % Return default filtering struct
    pars.Inverted_Logic = true; % Set true to indicate that sync bit logic is inverted (default) or false if it is non-inverted.
    pars.Linear_Fit_Order = 1; % For Subtract_Linear_Fit polynomial detrend that is applied to individual trials (DIFFERENT FROM APPLY_POLYNOMIAL_DETREND IN APPLY_EMG_FILTERS)
    pars.Link_Axes = true;  % Set false to allow axes limits to adaptively set independently
    pars.N_SD_RMS = 3.5; % Number of times the RMS to multiply signal by when computing YLIM if it is not manually specified.
    pars.N_Individual_Max = 10; % Max. number of individual traces to superimpose
    [pars.Output_Root, pars.Input_Root] = parameters('generated_data_folder', 'raw_data_folder'); % Location where output figures are saved.
    pars.Plot_Stim_Period = true; % Plot stim artifact with red stem lines?
    pars.Style = "Shaded"; % Can be "Shaded" | "Individual"
    pars.Start_Linear_Fit = 4.75; % Start the linear fit subtraction (if Subtract_Linear_Fit is true) on or after the sample corresponding to this value (ms)
    pars.Subtract_Linear_Fit = true; % Set false to skip the linear-fit subtraction.
    pars.Sync_Bit = nan; % The bit address for STIM sync TTL signal on TRIGGERS channel of TMSi.
    pars.T = [-15, 80]; % Time for epochs (milliseconds)
    pars.T_RMS = [12, 60]; % Time epoch for computing RMS
    pars.Trigger_Channel = 'TRIGGER'; % Name of Trigger Channel
    pars.XLim = []; % If empty, use auto-scale, otherwise, fixed scale
    pars.YLim = []; % If empty, use auto-scale, otherwise, fixed scale

    % % % END DEFAULT PARS STRUCT FIELD DEFINITIONS % % %
end
% Handle parsing of `pars`
pars = utils.parse_parameters(pars, varargin{:});
if ~isstruct(pars.Filtering)
     pars.Filtering = utils.get_default_filtering_pars(pars.Acquisition_Type, pars.EMG_Type, pars.Filtering);
end

x = io.load_tmsi(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, pars.File_Type, pars.Input_Root);
if isempty(x)
    fig = gobjects(1);
    return;
elseif numel(x) > 1
    warning('Should only return exactly one or no data matches. Check load_tmsi_raw. Block skipped.');
    fig = gobjects(1);
    return;
end

tank = sprintf('%s_%04d_%02d_%02d', SUBJ, YYYY, MM, DD); % data "tank"
block = sprintf('%s_%s_%d', tank, ARRAY, BLOCK); % experimental "block" (recording within tank)
gen_data_folder = fullfile(pars.Output_Root, SUBJ, tank, num2str(BLOCK));

fig = default.figure(block, 'Position', [0.1 0.1 0.8 0.8]);

% v1.5: After adding `get_default_filtering_pars() etc.`
fig.UserData = struct('x', x, 'version', 2.0); % Associate thse data to the figure.
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
    triggers = in.sync_data;
else
    [stops, trigs, triggers] = utils.parse_bit_sync(x, pars.Sync_Bit, gen_data_folder, pars.Inverted_Logic, pars.Trigger_Channel);
end
[z, ~, pars.Filtering] = utils.apply_emg_filters(x, pars.Filtering, x.sample_rate, trigs, stops);

L = tiledlayout(fig, 8, 8);

% Set up axes labeling. Only label the bottom row and left column.
% lab = strcat(repmat("UNI-", 64, 1), num2str((1:64)', '%02d'));
lab = string({channels.alternative_name});
lab = lab(contains(lab, "UNI"));
xlab = grid.vec_to_grid(lab);
ylab = xlab;
xlab(setdiff(1:64, 8:8:64)) = "";
xlab = grid.grid_to_vec(xlab); % Put them back into order for indexing in loop.
ylab(9:end) = "";
ylab = grid.grid_to_vec(ylab);
iTile = flipud(reshape(1:64, 8, 8)');
iTile = iTile(:);

n_pre = -1 * round(pars.T(1) * 1e-3 * x.sample_rate); % Convert to seconds, then samples
n_post = round(pars.T(2) * 1e-3 * x.sample_rate);  % Convert to seconds, then samples
t_sweep = (-n_pre:n_post)/x.sample_rate * 1e3; % Convert from samples to seconds, then milliseconds
i_start_fit = find(t_sweep >= pars.Start_Linear_Fit, 1, 'first'); % Sample index to start linear fit
if isinf(pars.End_Linear_Fit)
    i_end_fit = numel(t_sweep);
else
    i_end_fit = find(t_sweep <= pars.End_Linear_Fit, 1, 'last'); % Sample index to end linear fit 
end
noise_bandwidth = 1e-6; % Compute noise bandwidth

for ich = 1:64
    ax = nexttile(L, iTile(ich));
    if ~pars.Filtering.Add_To_Plots(ich)
        continue;
    end
    [A, X, trigs] = math.triggered_average(trigs, z(ich, :), n_pre, n_post, false, false, false);
    T = math.triggered_average(trigs, triggers, n_pre, n_post, false, false, false);
    if pars.Subtract_Linear_Fit
        X = X'; % Transpose so columns are trials.
        X(i_start_fit:i_end_fit, :) = detrend(X(i_start_fit:i_end_fit, :), pars.Linear_Fit_Order);
        X = X'; % Return to original orientation.
    end
    if contains(lower(pars.Filtering.Name), "rectified")
        X = max(X - mean(X(:, t_sweep < 0), 2), zeros(size(X)));
        A = mean(X, 1);
    end
    if pars.Filtering.Apply_Max_Rescale
        X = X./max(max(abs(X))); 
        A = mean(X, 1);
    end
    
    if size(X, 2) ~= numel(t_sweep)
        warning('Timing mismatch - block skipped.');
        delete(fig); 
        fig = gobjects(1);
        return;
    end
    set(ax, 'NextPlot', 'add', 'FontName', 'Tahoma', 'FontSize', 8, ...
        'LineWidth', 1.25, 'ButtonDownFcn', @(src, evt)callback.handleAxesClick(src, evt), ...
        'UserData', struct('name', lab(ich), 'filtering', pars.Filtering, 'block', block, 'data', X));
    xlabel(ax, xlab(ich), 'FontName', 'Tahoma', 'Color', 'k', 'FontSize', 9);
    ylabel(ax, ylab(ich), 'FontName', 'Tahoma', 'Color', 'k', 'FontSize', 9);
    if ~isempty(pars.XLim)
        xlim(ax, pars.XLim);
    end
    if ~isempty(pars.YLim)
        ylim(ax, pars.YLim);
    else
        noise_bandwidth = max(noise_bandwidth, rms(A((t_sweep >= pars.T_RMS(1)) & (t_sweep <= pars.T_RMS(2)))) * pars.N_SD_RMS);
        if contains("rectified", lower(pars.Filtering.Name))
            ylim(ax, [0, noise_bandwidth]);
        else
            ylim(ax, [-noise_bandwidth, noise_bandwidth]);
        end
    end
    if pars.Style == "Individual"
        if size(X, 1) > pars.N_Individual_Max
             X = X(randsample(size(X, 1), pars.N_Individual_Max), :);
             if ich == 1 % Only mention this one time.
                fprintf(1, '\t->\tSuperimposing %d individual traces...\n', pars.N_Individual_Max);
             end
        end 
    end
    T(T > 0) = max(max(X));
    if pars.Plot_Stim_Period
        stem(ax, t_sweep, T, 'LineWidth', 2, 'Color', 'r', 'Marker', 'none', 'ButtonDownFcn', @(src, evt)callback.handleAxesClick(src.Parent, evt)); 
    end
    if pars.Style == "Individual"
        plot(ax, t_sweep, X, ...
            'LineWidth', 0.5, 'Color', [0.45 0.45 0.45], ...
            'Tag', 'Dispersion', ...
            'LineStyle', ':', 'ButtonDownFcn', @(src, evt)callback.handleAxesClick(src.Parent, evt));
        if isempty(pars.YLim)
            ylim(ax, [min(X(:, t_sweep > 12.5), [], 'all'), max(X(:, t_sweep > 12.5), [], 'all')]); 
        end
    else
        faceData = [1:(2*numel(A)), 1]; 
        xx = t_sweep(:);
        xx = [xx; flipud(xx)]; %#ok<AGROW>
        yy = A(:);
        sd = std(X, [], 1);
        sd = sd(:);
        yy = [yy-sd; flipud(yy+sd)];
        vertexData = horzcat(xx, yy);
        patch(ax, 'Faces', faceData, 'Vertices', vertexData, ...
            'FaceAlpha', 0.25, 'FaceColor', [0.25 0.25 0.25], ...
            'Tag', 'Dispersion', ...
            'DisplayName', '\pm1 SD');
    end
    plot(ax, t_sweep, A, 'LineWidth', 2, 'Color', 'k', ...
        'ButtonDownFcn', @(src, evt)callback.handleAxesClick(src.Parent, evt), ...
        'Tag', 'STA');  % Plot X for the individual ones
    
end
if pars.Link_Axes
    linkaxes(findobj(L.Children, 'type', 'axes'), 'xy'); % Share common limits.
else
    linkaxes(findobj(L.Children, 'type', 'axes'), 'off');
end
str = utils.get_filtering_label_string(pars.Filtering);
title(L, [char(strrep(block, '_', '\_')), ': ' str newline 'Stim Averages (solid line | N = ' char(num2str(numel(trigs)))  ')'], ...
    'FontName', 'Tahoma', ...
    'Color', 'k', 'FontSize', 14, 'FontWeight', 'bold');
finfo = strsplit(block, '_');
out_folder = fullfile(pars.Output_Root, SUBJ, tank, finfo{end}, '.averages');
out_name = fullfile(out_folder, sprintf('%s_%d_%d_Mean_Array_EMG__%s', block, round(pars.T(1)), round(pars.T(2)), pars.Filtering.Name));
set(fig, 'WindowKeyPressFcn', @(src, evt)callback.handleCommonWindowKeyPresses(src, evt, out_name));
end
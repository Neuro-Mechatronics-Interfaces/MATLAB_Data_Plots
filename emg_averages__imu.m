function fig = emg_averages__imu(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, varargin)
%EMG_AVERAGES__IMU  EMG processing for TMSi IMU channels
%
% Syntax:
%   fig = plot.emg_averages__imu(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, pars);
%
% This should be accessed via `"IMU"` EMG_Type parameter in
% `plot_emg_averages`. 
%
% See also: Contents, plot.emg_averages

if (numel(varargin) == 1) && isstruct(varargin{1})
    pars = varargin{1};
else
    pars = plot.parameters('emg_averages');
    pars.EMG_Type = "Accelerometer";
    % % % END DEFAULT PARS STRUCT FIELD DEFINITIONS % % %
    % Handle parsing of `pars`
    pars = utils.parse_parameters(pars, varargin{:});
end

if ~isstruct(pars.Filtering)
     pars.Filtering = get_default_filtering_pars(pars.Acquisition_Type, pars.EMG_Type, pars.Filtering);
end

[x, info] = io.load_data(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, pars.File_Type, pars.Input_Root);
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
fig = default.figure(block, 'Position', pars.Position);


% v1.4: After adding `get_default_filtering_pars() etc.`
fig.UserData = struct('x', x, 'version', pars.Version, 'pars', pars); % Associate thse data to the figure.
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
if (numel(trigs) < 1) || (numel(stops) < 1)
    warning("Empty sync vector (trigs): check if TTL on TRIGGERS channel was present/parsed using correct bit.");
    delete(fig);
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

iIMU = ( contains({channels.alternative_name}, 'X')' | contains({channels.alternative_name}, 'Y')' | ...
         contains({channels.alternative_name}, 'Z')' ) & (sum(abs(x.samples-mean(x.samples, 2)),2) > eps);
if sum(iIMU) == 0
    delete(fig);
    fig = gobjects(1);
    fprintf(1, 'No Accelerometer channels for recording: <strong>%s</strong>\n', block); 
    return;
else
    channels = channels(iIMU);
end
data = x.samples(iIMU, :)';
lab = string({channels.alternative_name});
if pars.Filtering.Apply_Virtual_Reference
    tmp = true;
    pars.Filtering.Apply_Virtual_Reference = false;
    mu = mean(x.samples(1:64, :), 1)';
    for iCh = 1:size(data, 2)
        data(:, iCh) = data(:, iCh) - (mu .* (rms(data(:, iCh)) ./ rms(mu)));
    end
else
    tmp = false;
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
if isempty(pars.Filtered_Data)
    % Trigs is returned because the filtering function can exclude
    % out-of-bounds trigger sample indices based on stim-artifact-rejection
    % sample epoch width.
    [z, ~, pars.Filtering, trigs] = utils.apply_emg_filters(data, pars.Filtering, x.sample_rate, trigs, stops);
else
    z = pars.Filtered_Data;
end


pars.Filtering.Apply_Virtual_Reference = tmp; % Revert
if (isnan(pars.N_Rows)) || (isnan(pars.N_Columns))
    L = tiledlayout(fig, 'flow');
else
    L = tiledlayout(fig, pars.N_Rows, pars.N_Columns); 
end


n_pre = -1 * round(pars.T(1) * 1e-3 * x.sample_rate); % Convert to seconds, then samples
n_post = round(pars.T(2) * 1e-3 * x.sample_rate);  % Convert to seconds, then samples
t_sweep = (-n_pre:n_post)/x.sample_rate * 1e3; % Convert from samples to seconds, then milliseconds
noise_bandwidth = 1e-6; % Compute noise bandwidth
i_start_fit = find(t_sweep >= pars.Start_Linear_Fit, 1, 'first'); % Sample index to start linear fit
if isinf(pars.End_Linear_Fit)
    i_end_fit = numel(t_sweep);
else
    i_end_fit = find(t_sweep <= pars.End_Linear_Fit, 1, 'last'); % Sample index to end linear fit 
end
stim = utils.get_tmsi_stim_data(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, pars.Input_Root);
ax = [];
if pars.Subtract_Mean || pars.Filtering.Subtract_Cross_Trial_Mean
    pars.Filtering.Subtract_Cross_Trial_Mean = true;
    pars.Subtract_Mean = true;
end
for iCh = 1:sum(iIMU)
    
    % Trigs is returned because the averaging function can exclude
    % out-of-bounds trigger sample indices based on snippet
    % sampling epoch width (samples pre- and post-TTL marker).
    [A, X, trigs] = math.triggered_average(trigs, z(iCh, :), n_pre, n_post, false, false, false);
    T = math.triggered_average(trigs, triggers, n_pre, n_post, false, false, false);
    if pars.Subtract_Mean
        X = abs(X - mean(X,1)); 
    end
    if all(abs(A) < eps) % Then this channel is "empty"
        continue;
    end
    if pars.Subtract_Linear_Fit
        X = X'; % Transpose so columns are trials.
        X(i_start_fit:i_end_fit, :) = detrend(X(i_start_fit:i_end_fit, :), pars.Linear_Fit_Order);
        X = X'; % Return to original orientation.
    end
    if contains(lower(pars.Filtering.Name), "rectified")
        X = max(X - mean(X(:, t_sweep < 0), 2), zeros(size(X)));
        A = mean(X, 1);
    end
    
    ax = nexttile(L);
    if size(X, 2) ~= numel(t_sweep)
        warning('Timing mismatch - block skipped.');
        delete(fig); 
        fig = gobjects(1);
        return;
    end
    set(ax, 'NextPlot', 'add', 'FontName', 'Tahoma', 'FontSize', 8, ...
        'LineWidth', 1.25, 'ButtonDownFcn', @(src, evt)callback.handleAxesClick(src, evt), ...
        'UserData', struct('name', lab(iCh), 'filtering', pars.Filtering, 'block', block));
    if ~isempty(pars.XLim)
        xlim(ax, pars.XLim);
    end
    if ~isempty(pars.YLim)
        ylim(ax, pars.YLim);
    else
        noise_bandwidth = max(noise_bandwidth, rms(A((t_sweep >= pars.T_RMS(1)) & (t_sweep <= pars.T_RMS(2)))) * pars.N_SD_RMS);
        if contains("rectified", lower(pars.Filtering.Name)) || pars.Filtering.Subtract_Cross_Trial_Mean
            ylim(ax, [0, noise_bandwidth]);
        else
            ylim(ax, [-noise_bandwidth, noise_bandwidth]);
        end
    end
    if pars.Style == "Individual"
        if size(X, 1) > pars.N_Individual_Max
             X = X(randsample(size(X, 1), pars.N_Individual_Max), :);
             if iCh == 1 % Only mention this one time.
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
            'Tag', 'Dispersion', ...
            'FaceAlpha', 0.25, 'FaceColor', [0.25 0.25 0.25], ...
            'DisplayName', '\pm1 SD');
    end
    plot(ax, t_sweep, A, 'LineWidth', 2, 'Color', 'k', ...
        'ButtonDownFcn', @(src, evt)callback.handleAxesClick(src.Parent, evt), ...
        'Tag', 'STA');  % Plot X for the individual ones
    ylabel(ax, channels(iCh).unit_name, 'FontName', 'Tahoma', 'Color', 'k');
    xlabel(ax, 'ms', 'FontName', 'Tahoma', 'Color', 'k');
    name = strrep(channels(iCh).alternative_name, ' ', '');
    try
        title(ax, strrep(stim.map.Muscles.(name), '_', '-'), 'FontName', 'Tahoma', 'Color', 'k', 'FontSize', 14, 'FontWeight', 'bold');
    catch
        title(ax, name, 'FontName', 'Tahoma', 'Color', 'k', 'FontSize', 14, 'FontWeight', 'bold');
    end
    if ax.XLim(1) >= 0
        utils.add_sd_threshold(ax, 6.5, t_sweep, A, 'LabelHorizontalAlignment', 'right');
    else
        utils.add_sd_threshold(ax, 6.5, t_sweep, A, 'LabelHorizontalAlignment', 'left');
    end
end
if isempty(ax)
    title(L, "No Bipolar EMG", 'FontName', 'Tahoma', ...
        'Color', 'r', 'FontSize', 14, 'FontWeight', 'bold');
else
    if pars.Link_Axes
        linkaxes(findobj(L.Children, 'type', 'axes'), 'xy'); % Share common limits.
    end
        str = utils.get_filtering_label_string(pars.Filtering);
        N = numel(trigs);
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
            title(L, [char(block_a), ': ' str newline ' (N = ' char(num2str(N))  ') IMU'], ...
                'FontName', 'Tahoma', ...
                'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
        else
            title(L, [char(strrep(block, '_', '\_')), ': ' str newline ' (N = ' char(num2str(N))  ') IMU'], ...
                'FontName', 'Tahoma', ...
                'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
        end

        finfo = strsplit(block, '_');
        out_folder = fullfile(pars.Output_Root, SUBJ, tank, finfo{end}, '.averages');
        out_name = fullfile(out_folder, sprintf('%s_%d_%d_Mean_Bipolar_EMG__%s', block, round(pars.T(1)), round(pars.T(2)), pars.Filtering.Name));
        set(fig, 'WindowKeyPressFcn', @(src, evt)callback.handleCommonWindowKeyPresses(src, evt, out_name));
    end
end
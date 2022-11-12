function fig = emg_rms(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, varargin)
%EMG_RMS  Plot RMS heatmap with filled contour lines.
%
% Syntax:
%   fig = plot.emg_rms(SUBJ, YYYY, MM, DD, ARRAY, BLOCK);
%   fig = plot.emg_rms(___, 'Name', value,...);
%
% Example 1:
%   % Run a single figure with default parameters.
%   fig = plot.emg_rms('Frank', 2021, 11, 18, "A", 43);
%   
% Example 2:
%   % Export a batch run of average figures.
%   plot.emg_rms('Frank', 2021, 11, 18, ["A", "B"], 0:105);
%
% Example 3:
%   % Modify specific parameters for a figure.
%   plot.emg_rms('Frank', ...
%       2021, 12, 9, "B", 156, ...
%       'T', [-22.5, 27.5], ... % Modifies epoch time (milliseconds)
%       'Output_Root', 'G:\Shared drives\NML_NHP\DARPA_N3\preliminary'); % Change where figures are saved 
%
% Inputs:
%   SUBJ  - Subject name (e.g. 'Frank' or "Frank")
%   YYYY  - Year (numeric or string, e.g. 2021 or "2021" or '2021')
%   MM    - Month (numeric or string, e.g. 11 or "11" or '11')
%   DD    - Day   (numeric or string, e.g. 18 or "18" or '18')
%   ARRAY - "A" or "B" (or 'A' or 'B') or ["A", "B"] or {'A', 'B'}
%   BLOCK - Recording parameter key (block; numeric or string, e.g. 0 or "0")
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
% See also: Contents, io.load_tmsi

% % % % SEE PARS STRUCT BELOW % % % %
if (numel(varargin) > 0) && isstruct(varargin{1})
    pars = varargin{1};
    varargin(1) = [];
else
    pars = plot.parameters('emg_rms'); 
end
% % % END DEFAULT PARS STRUCT FIELD DEFINITIONS % % %

% Handle parsing of `pars`
pars = utils.parse_parameters(pars, varargin{:});
if ~isstruct(pars.Filtering)
     pars.Filtering = utils.get_default_filtering_pars(pars.Acquisition_Type, pars.EMG_Type, pars.Filtering);
end


if (numel(ARRAY) > 1) || (numel(BLOCK) > 1)
    fig = gobjects(numel(BLOCK), numel(ARRAY));
    for iB = 1:numel(BLOCK)
        for iA = 1:numel(ARRAY)
            if nargout > 0
                fig(iB, iA) = plot.emg_rms(SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), pars);
            else
                plot.emg_rms(SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), pars);
            end
        end
    end
    return;
end

[YYYY, MM, DD] = utils.parse_date_args(YYYY, MM, DD);
% x = io.load_tmsi('Frank', 2021, 12, 9, "B", 155, ".poly5", "R:\NMLShare\raw_data\primate");
if isempty(pars.Data)
    x = io.load_tmsi(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, pars.File_Type, pars.Input_Root);
else
    x = pars.Data;
end

if isempty(x)
    fig = gobjects(1);
    return;
end

tank = sprintf('%s_%04d_%02d_%02d', SUBJ, YYYY, MM, DD); % data "tank"
block = sprintf('%s_%s_%d', tank, ARRAY, BLOCK); % experimental "block" (recording within tank)
gen_data_folder = fullfile(pars.Output_Root, SUBJ, tank, num2str(BLOCK));
sync_data_in_file = fullfile(gen_data_folder, sprintf('%s_sync.mat', x.name));

if exist(sync_data_in_file, 'file')==0
    if isnan(pars.Sync_Bit)
        error('No sync data file (<strong>%s</strong>): please run `parse_bit_sync` first! (or, set Sync_Bit to non-NaN value)', sync_data_in_file);
    end
    [stops, trigs] = utils.parse_bit_sync(x, pars.Sync_Bit, gen_data_folder);
else
    in = load(sync_data_in_file, 'offset', 'onset', 'sync_data');
    trigs = in.offset;
    stops = in.onset;
end
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
[Z, ~, pars.Filtering, trigs] = utils.apply_emg_filters(x, pars.Filtering, x.sample_rate, trigs, stops);
n_pre = -1 * round(pars.T(1) * 1e-3 * x.sample_rate); % Convert to seconds, then samples
n_post = round(pars.T(2) * 1e-3 * x.sample_rate);  % Convert to seconds, then samples
t_sweep = (-n_pre:n_post)/x.sample_rate * 1e3; % Convert from samples to seconds, then milliseconds

i_pre = t_sweep < pars.Pre_Stimulus_RMS_ms;
i_post = t_sweep > pars.Post_Stimulus_RMS_ms;
Zt = grid.triggered_array(Z', trigs, n_pre, n_post);
if pars.Subtract_Mean || pars.Filtering.Subtract_Cross_Trial_Mean
    Zt = Zt - mean(Zt,3); 
    pars.Filtering.Subtract_Cross_Trial_Mean = true;
    pars.Subtract_Mean = true;
end

Z_pre = squeeze(rms(Zt(:, i_pre, :), 2));
Z_post = squeeze(rms(Zt(:, i_post, :), 2));
Z_ratio = exp(mean(log(Z_post) - log(Z_pre), 2));
Z_ratio = reshape(Z_ratio, 8, 8);

if pars.Debug
    fprintf(1, '\t[DEBUG]\t->\tUsing <strong>%d</strong> PRE samples and <strong>%d</strong> POST samples.\n', sum(i_pre), sum(i_post));
    db_fig = figure('Name', 'RMS: Debugger', ...
        'Color', 'w', 'Units', 'Normalized', 'Position', [0.1 0.1 0.8 0.8]);
    db_L = tiledlayout(db_fig, 8, 8);
    iTile = flipud(reshape(1:64, 8, 8)');
    
    db_t_pre = repmat(t_sweep(i_pre), 1, size(Zt, 3));
    db_t_pre = db_t_pre(:);
    db_t_post = repmat(t_sweep(i_post), 1, size(Zt, 3));
    db_t_post = db_t_post(:);
    for ii = 1:64
        db_ax = nexttile(db_L, iTile(ii));
        lab = string(sprintf('UNI-%02d', ii));
        set(db_ax, 'NextPlot', 'add', 'XLim', pars.T, ...
            'FontName', 'Tahoma', 'FontSize', 6, ...
            'YScale', 'log', 'YLim', [1e-3, 1e3], ...
            'XTick', [pars.T(1) pars.Pre_Stimulus_RMS_ms pars.Post_Stimulus_RMS_ms pars.T(2)], 'XTickLabel', ["", "", "", ""], ...
            'YTick', [0.01, 1, 100], 'YTickLabel', ["", "", ""], ...
            'UserData', struct('name', lab, 'filtering', pars.Filtering, 'block', block, 'data', squeeze(Zt(ii, :, :))'), ...
            'ButtonDownFcn', @cb.handleAxesClick);
        
        db_tmp_pre = squeeze(Zt(ii, i_pre, :));
        scatter(db_ax, db_t_pre, db_tmp_pre(:), ...
            'DisplayName', 'Pre-Stim', ...
            'MarkerFaceColor', 'r', 'Marker', 'o', ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.015);
        db_tmp_pre_mu = mean(db_tmp_pre, 2);
        line(db_ax, t_sweep(i_pre), db_tmp_pre_mu, ...
            'DisplayName', 'E[Pre-Stim]', 'LineWidth', 2, 'Color', [0.2 0.2 0.2]);
        db_h = scatter(db_ax, mean(db_t_pre), mean(db_tmp_pre_mu), ...
            'MarkerFaceColor', 'r', ...
            'MarkerFaceAlpha', 0.5, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerEdgeAlpha', 0.75, ...
            'LineWidth', 1.5, ...
            'Marker', 's', ...
            'SizeData', round(mean(db_tmp_pre_mu.^2) * 10), ...
            'DisplayName', 'E[RMS_p_r_e]');
        db_h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        if ii == 1
            text(db_ax, mean(db_t_pre), 0.01, 'E[RMS_p_r_e]',...
                'BackgroundColor', 'k', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                'Color', 'w', 'FontWeight', 'bold', 'FontSize', 8, 'FontName', 'Tahoma'); 
        end
        
        db_tmp_post = squeeze(Zt(ii, i_post, :));
        scatter(db_ax, db_t_post, db_tmp_post(:), ...
            'DisplayName', 'Post-Stim', ...
            'MarkerFaceColor', 'b', 'Marker', 'o', ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.015);
        db_tmp_post_mu = mean(db_tmp_post, 2);
        line(db_ax, t_sweep(i_post), db_tmp_post_mu, ...
            'DisplayName', 'E[Post-Stim]', 'LineWidth', 2, 'Color', [0.2 0.2 0.2]);
        db_h = scatter(db_ax, mean(db_t_post), mean(db_tmp_post_mu), ...
            'MarkerFaceColor', 'b', ...
            'MarkerFaceAlpha', 0.5, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerEdgeAlpha', 0.75, ...
            'LineWidth', 1.5, ...
            'Marker', 's', ...
            'SizeData', round(mean(db_tmp_post_mu.^2) * 10), ...
            'DisplayName', 'E[RMS_p_o_s_t]');
        db_h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        if ii == 1
            text(db_ax, mean(db_t_post), 100, 'E[RMS_p_o_s_t]',...
                'BackgroundColor', 'k', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'Color', 'w', 'FontWeight', 'bold', 'FontSize', 8, 'FontName', 'Tahoma'); 
        end
        
        if rem(ii, 8) == 1
            xlabel(db_ax, 'ms', 'FontName', 'Tahoma', 'FontSize', 8, 'Color', 'k');
            set(db_ax, 'XTickLabel', db_ax.XTick);
        end
        if ceil(ii/8) == 1
            ylabel(db_ax, '\muV', 'FontName', 'Tahoma', 'FontSize', 8, 'Color', 'k'); 
            set(db_ax, 'YTickLabel', ["0.01", "1", "100"]);
        end
    end
    linkaxes(findobj(db_L.Children, 'Type', 'axes'), 'y');
    set(db_fig, 'WindowKeyPressFcn', @(src, evt)cb.handleCommonWindowKeyPresses(src, evt, fullfile(gen_data_folder, '.debug', sprintf('%s_RMS-Grid', block)), false));
    waitfor(db_fig);
    clear db_fig db_L db_ax db_t_pre db_t_post db_tmp_pre db_tmp_post db_tmp_pre_mu db_tmp_post_mu db_h;
end

% Generate figure
if isempty(pars.Axes)
    fig = figure('Name', 'RMS-Map',...
             'Units','Normalized', ...
             'Position',[0.1 0.1 0.8 0.8],...
             'Color', 'w');
    ax = axes(fig, ...
    'NextPlot', 'add', 'Color', 'r', ...
    'XColor', 'k', 'YColor', 'k', 'FontName', 'Arial', 'FontSize', 16, ...
    'YDir', 'normal', 'XTickLabel', 1:8:57, 'YTickLabel', 1:8);
else
    fig = pars.Axes.Parent;
    if ~isa(fig, 'matlab.ui.Figure')
        fig = fig.Parent; % In case axes is in a Layout or Panel or something
    end
    ax = pars.Axes;
    set(ax, 'NextPlot', 'add');
end


if ~isempty(pars.CLim)
    set(ax, 'CLim', pars.CLim);
    RMS_Max_Response_Ratio = pars.CLim(2);
else
    if ~isempty(pars.RMS_Max_Response_Ratio)
        pars.CLim = [1, pars.RMS_Max_Response_Ratio];
        RMS_Max_Response_Ratio = pars.RMS_Max_Response_Ratio;
        set(ax, 'CLim', pars.CLim);
    end
end
if ~isempty(pars.XLim)
    % old values : [0.5 8.5]
    set(ax, 'XLim', pars.XLim);
end
if ~isempty(pars.YLim)
    % old values: [0.5 8.5]
    set(ax, 'YLim', pars.YLim);
end

% Since we're now dealing in the image space where (0, 0) is at the bottom
% left, we need to flip our grid to reflect that.
Z_ratio = flipud(Z_ratio);
% Make any NaN values zero (artificially) to "clamp" the top and bottom
% edges of our patch if we use "Differential2" mode for filtering.
Z_ratio(isnan(Z_ratio)) = 0;
x0 = linspace(0.5, 8.5, 8);
y0 = fliplr(linspace(0.5, 8.5, 8));
[X, Y] = meshgrid(x0, y0');

xq = linspace(0.5, 8.5, 64);
yq = fliplr(linspace(0.5, 8.5, 64))';
[Xq, Yq] = meshgrid(xq, yq);
Z_ratio_q = interp2(X, Y, Z_ratio, Xq, Yq, 'spline');

contourf(ax, xq, yq, Z_ratio_q, 'LineWidth', 2);

if isempty(pars.CLim) && isempty(pars.RMS_Max_Response_Ratio)
    pars.CLim = get(ax, 'CLim'); 
    RMS_Max_Response_Ratio = pars.CLim(2); 
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
    title(ax, [char(block_a), ': ' str newline ' (N = ' char(num2str(N))  ') RMS'], ...
        'FontName', 'Tahoma', ...
        'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
else
    title(ax, [char(strrep(block, '_', '\_')), ': ' str newline ' (N = ' char(num2str(N))  ') RMS'], ...
        'FontName', 'Tahoma', ...
        'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
end

xlabel(ax, 'Channel', 'FontName', 'Tahoma', 'Color', 'k', 'FontSize', 11);
ylabel(ax, 'Channel', 'FontName', 'Tahoma', 'Color', 'k', 'FontSize', 11);
% text(ax, 0.25, 1.5,'Lateral, Distal', ...
%     'HorizontalAlignment','Right', ...
%     'FontName', 'Arial', 'Color', 'k',...
%     'FontSize', 14, 'Rotation', 90);
% text(ax, 9, 7.5,'Medial, Proximal', ...
%     'FontName', 'Arial', 'Color', 'k', ...
%     'FontSize', 14, 'Rotation', -90);

if ~isempty(pars.CLim)
    set(ax, 'CLim', pars.CLim);
    clim_min = pars.CLim(1);
    RMS_Max_Response_Ratio = pars.CLim(2);
else
    clim_min = 1;
    RMS_Max_Response_Ratio = max(max(Z_ratio));
end
if isempty(pars.RMS_Response_Ratio_Threshold)
    pars.RMS_Response_Ratio_Threshold = round(RMS_Max_Response_Ratio/2);
end

stim = utils.get_tmsi_stim_data(SUBJ, YYYY, MM, DD, ARRAY, BLOCK);
[xg, yg] = meshgrid(1:8, (8:-1:1)');
if isempty(pars.RMS_Response_Ratio_Threshold)
    pars.RMS_Response_Ratio_Threshold = 0.8 * ax.CLim(2);
end
idx = find(Z_ratio >= pars.RMS_Response_Ratio_Threshold);
if ~isempty(idx)
     channels = flipud(reshape(1:64, 8, 8));
     for ii = 1:numel(idx)
          ch = channels(idx(ii));
          chname = sprintf('UNI%02d', ch);
          str = strrep(stim.map.Muscles.(chname), '_', ' ');
          text(ax, xg(idx(ii)), yg(idx(ii)), str, ...
              'FontName', 'Tahoma', ...
              'FontSize', 10, ...
              'FontWeight', 'bold', ...
              'BackgroundColor', 'k', ...
              'Color', 'w', ...
              'HorizontalAlignment', 'center', ...
              'VerticalAlignment', 'middle');
     end
end
cbar = colorbar(ax);
cbar.Label.String = 'RMS_p_o_s_t/RMS_p_r_e';
if (RMS_Max_Response_Ratio > pars.RMS_Response_Ratio_Threshold) && (clim_min < pars.RMS_Response_Ratio_Threshold)
    cbar.Ticks = [clim_min, pars.RMS_Response_Ratio_Threshold, RMS_Max_Response_Ratio];
    cbar.TickLabels = [...
        string(sprintf('\\color{black}%3.2f', clim_min)), ...
        string(sprintf('\\color{red}%3.2f', pars.RMS_Response_Ratio_Threshold)), ...
        string(sprintf('\\color{black}%3.2f', RMS_Max_Response_Ratio))];
else
    cbar.Ticks = [clim_min, RMS_Max_Response_Ratio];
    cbar.TickLabels = [...
        string(sprintf('\\color{black}%3.2f', clim_min)), ...
        string(sprintf('\\color{black}%3.2f', RMS_Max_Response_Ratio))];
end
if nargout < 1
    tank = sprintf('%s_%04d_%02d_%02d', SUBJ, YYYY, MM, DD); % data "tank"
    block = sprintf('%s_%s_%d', tank, ARRAY, BLOCK); % experimental "block" (recording within tank)
    out_folder = fullfile(pars.Output_Root, SUBJ, tank, 'figures', 'RMS', pars.Filtering.Name);
    if exist(out_folder, 'dir') == 0
        try %#ok<TRYNC>
            mkdir(out_folder);
        end
    end
    out_name = fullfile(out_folder, sprintf('%s_%d_%d', block, round(pars.T(1)), round(pars.T(2))));
    default.savefig(fig, out_name, sprintf("%s_EMG", pars.EMG_Type), true);
    
    out_folder_2 = fullfile(pars.Output_Root, SUBJ, tank, num2str(BLOCK));
    if exist(out_folder_2, 'dir') == 0
        try %#ok<TRYNC>
            mkdir(out_folder_2);
        end
    end
    default.savefig(fig, fullfile(out_folder_2, block), sprintf('%d_%d_%s_EMG_%s', round(pars.T(1)), round(pars.T(2)), pars.EMG_Type, pars.Filtering.Name), false); 
end
%default.savefig(fig, fullfile(gen_data_folder, x.name, sprintf('%d_%d_', round(pars.T(1)), round(pars.T(2))), 'RMS-Map_', pars.Filtering.Name), nargout > 0);
end
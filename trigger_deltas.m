function fig = trigger_deltas(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, varargin)
%TRIGGER_DELTAS  Handle plotting stem display of times between each "trigger" (sync) sample.
%
% Syntax:
%   fig = plot.trigger_deltas(SUBJ, YYYY, MM, DD, ARRAY, BLOCK);
%   fig = plot.trigger_deltas(___, 'Name', value,...);
%
% Example 1:
%   % Run a single figure with default parameters.
%   fig = plot.trigger_deltas('Rascal', 2022, 2, 24, "A", 36);
%   
% Example 2:
%   % Export a batch run of average figures.
%   plot.trigger_deltas('Rascal', 2022, 2, 24, ["A", "B"], 0:47);
%
% Inputs:
%   SUBJ  - Subject name (e.g. 'Rascal' or "Rascal")
%   YYYY  - Year (numeric or string, e.g. 2022 or "2022" or '2022')
%   MM    - Month (numeric or string, e.g. 2 or "02" or '02')
%   DD    - Day   (numeric or string, e.g. 24 or "24" or '24')
%   ARRAY - "A" or "B" (or 'A' or 'B') or ["A", "B"] or {'A', 'B'}
%   BLOCK - Recording parameter key (block; numeric or string, e.g. 0 or "0")
%   varargin - (Optional) 'Name', value parameter field/value pairs. See
%                       `pars` struct for specific parameter names. A few
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
% See also: Contents, load_tmsi_raw, parse_bit_sync, parse_artifact_sync

pars = struct;
pars.Amplitude_Parameter = 'Current';
pars.Amplitude_Units = 'mA';
pars.Axes = []; % Specify Axes directly to superimpose multiple plots
pars.Axes_Properties = {}; % Cell array where each cell alternates as 'Name', value argument pairs for MATLAB `axes` constructor.
pars.Figure_Properties = {}; % Cell array where each cell alternates as 'Name', value inputs for MATLAB `figure` constructor.
pars.Font = {'FontName', 'Tahoma', 'FontSize', 18, 'Color', 'k'}; % Common "Font" argument for graphical text labels, titles, etc.
pars.Max_Stims_Shown = 50; % Maximum number to show
pars.Method = @parse_artifact_sync;  % Function handle to use for parsing sync signal. Must return [onset, offset] as first two arguments, which are vectors of sample indices.
pars.Output_Root = parameters('generated_data_folder');
pars.Sample_Rate = 4000;
pars.Stem_Properties = {}; % Cell array where each cell alternates as 'Name', value inputs for MATLAB `stem` constructor.
pars.Ymax = 100; % Maximum y-value (ms)
pars = utils.parse_parameters(pars, varargin{:});

if (numel(ARRAY) > 1) || (numel(BLOCK) > 1)
    if nargout > 0
        fig = gobjects(numel(BLOCK), numel(ARRAY));
    end
    for iB = 1:numel(BLOCK)
        for iA = 1:numel(ARRAY)
            if nargout > 0
                fig(iB, iA) = plot.trigger_deltas(SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), pars);
            else
                plot.trigger_deltas(SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), pars); 
            end
        end
    end
    return;
end
[YYYY, MM, DD] = utils.parse_date_args(YYYY, MM, DD);

[~, offset] = pars.Method(SUBJ, YYYY, MM, DD, ARRAY, BLOCK);
f = utils.get_block_name(SUBJ, YYYY, MM, DD, ARRAY, BLOCK);

dOffset = diff(offset) / (pars.Sample_Rate / 1e3);
if median(dOffset) > pars.Ymax
    pars.Ymax = 1.5 * median(dOffset);
end

if isempty(pars.Axes)
    fig = figure(...
        'Name', f.Block, 'Color', 'w', ...
        'Units', 'Normalized', 'Position', [0.15 0.45 0.6 0.3], ...
        pars.Figure_Properties{:});
    ax = axes(fig, 'NextPlot', 'add', 'LineWidth', 1.25, ...
        pars.Font{:}, ...
        'Color', 'w', ...
        'XColor', 'k', ...
        'YColor', 'k', ...
        'XLim', [1, min(numel(dOffset), pars.Max_Stims_Shown)], ...
        'YLim', [0, pars.Ymax], ...
        pars.Axes_Properties{:});
    title(ax, strrep(f.Block, "_", "\_"), ...
        pars.Font{:}, 'FontWeight', 'bold');
    xlabel(ax, 'Stim #', pars.Font{:});
    ylabel(ax, '\Deltat (ms)', pars.Font{:});
else
    ax = pars.Axes;
    fig = ax.Parent;
end

stem(ax, 1:numel(dOffset), dOffset, 'k-', ...
    'LineWidth', 1.25, ...
    'DisplayName', strrep(f.Block, '_', '-'), ...
    pars.Stem_Properties{:});

if nargout < 1
    out_folder = fullfile(pars.Output_Root, SUBJ, f.Tank, 'figures', 'Deltas');
    if exist(out_folder, 'dir') == 0
        try %#ok<*TRYNC>
            mkdir(out_folder);
        end
    end
    default.savefig(fig, fullfile(out_folder, f.Block), 'Deltas', true);
    out_folder_2 = fullfile(pars.Output_Root, SUBJ, f.Tank, num2str(BLOCK));
    if exist(out_folder_2, 'dir') == 0
        try
            mkdir(out_folder_2);
        end
    end
    default.savefig(fig, fullfile(out_folder_2, f.Block), 'Deltas', false); 
end

end
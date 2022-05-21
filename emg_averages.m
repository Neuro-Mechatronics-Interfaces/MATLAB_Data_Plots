function fig = emg_averages(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, varargin)
%EMG_AVERAGES  Handle plotting gridded averages (TMSi/HD-EMG array)
%
% Syntax:
%   fig = plot.emg_averages(SUBJ, YYYY, MM, DD, ARRAY, BLOCK);
%   fig = plot.emg_averages(___, 'Name', value,...);
%
% Example 1:
%   % Run a single figure with default parameters.
%   fig = plot.emg_averages('Frank', 2021, 11, 18, "A", 43);
%   
% Example 2:
%   % Export a batch run of average figures.
%   plot.emg_averages('Frank', 2021, 11, 18, ["A", "B"], 0:105);
%
% Example 3:
%   % Modify specific parameters for a figure.
%   plot.emg_averages('Frank', ...
%       2021, 11, 18, "A", 43, ...
%       'T', [-15, 80], ... % Modifies epoch time (milliseconds)
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
% See also: Contents, load_tmsi_raw

if ~isnumeric(YYYY)
    YYYY = str2double(YYYY);
end

if ~isnumeric(MM)
    MM = str2double(MM); 
end

if ~isnumeric(DD)
    DD = str2double(DD); 
end

if ~isnumeric(BLOCK)
    BLOCK = str2double(BLOCK); 
end

% % % % SEE PARS STRUCT BELOW % % % %
if (numel(varargin) == 1) && isstruct(varargin{1})
    pars = varargin{1};
else
    pars = struct;
    pars.Acquisition_Type = "TMSi";
    pars.EMG_Type = "Array"; % Can be: "Array" | "Bipolar"
    pars.End_Linear_Fit = inf; % Set to some value in milliseconds to say when linear fit should end
    pars.Filtering = get_default_filtering_pars(pars.Acquisition_Type, pars.EMG_Type, "Rectified", 'Polynomial_Detrend_Order', 7); % Return default filtering struct
    pars.File_Type = ".mat"; % Can be: ".mat" | ".poly5"
    pars.Inverted_Logic = true; % Set true to indicate that sync bit logic is inverted (default) or false if it is non-inverted.
    pars.Linear_Fit_Order = 1; % For Subtract_Linear_Fit polynomial detrend that is applied to individual trials (DIFFERENT FROM APPLY_POLYNOMIAL_DETREND IN APPLY_EMG_FILTERS)
    pars.Link_Axes = false;  % Set false to al axes limits to adaptively set independently
    pars.N_SD_RMS = 3.5; % Number of times the RMS to multiply signal by when computing YLIM if it is not manually specified.
    pars.N_Individual_Max = 10; % Max. number of individual traces to superimpose
    pars.N_Rows = nan; % Number of rows in grid layout
    pars.N_Columns = nan; % Number of columns in grid layout
    [pars.Output_Root, pars.Input_Root] = parameters('generated_data_folder', 'raw_data_folder'); % Location where output figures are saved.
    pars.Plot_Stim_Period = true; % Plot stim artifact with red stem lines?
    pars.Style = "Shaded"; % Can be: "Shaded" | "Individual"
    pars.Start_Linear_Fit = 4.75; % Start the linear fit subtraction (if Subtract_Linear_Fit is true) on or after the sample corresponding to this value (ms)
    pars.Subtract_Linear_Fit = false; % Set false to skip the linear-fit subtraction.
    pars.Sync_Bit = nan; % The bit address for STIM sync TTL signal on TRIGGERS channel of TMSi.
    pars.Use_Artifact_Sync = false;
    pars.Artifact_Sync_Rule = 'Mean'; % Synb_Bit alternative, choose the channel which has the greatest rms with likely stim artifact ( 'max_rms' | channel_number )
    pars.T = [-15, 80]; % Time for epochs (milliseconds)
    pars.T_RMS = [12, 60]; % Time epoch for computing RMS
    pars.Trigger_Channel = 'TRIGGER'; % Name of Trigger Channel
    pars.XLim = []; % If empty, use auto-scale, otherwise, fixed scale
    pars.YLim = []; % If empty, use auto-scale, otherwise, fixed scale
    % % % END DEFAULT PARS STRUCT FIELD DEFINITIONS % % %
    % Handle parsing of `pars`
    pars = utils.parse_parameters(pars, varargin{:});
    
end

if ~isstruct(pars.Filtering)
     pars.Filtering = get_default_filtering_pars(pars.Acquisition_Type, pars.EMG_Type, pars.Filtering);
end

if (numel(BLOCK) > 1) || (numel(ARRAY) > 1)
    if nargout < 1  % Then ensure we do not change "nargout" for actual implementation.
        for iB = 1:numel(BLOCK)
            for iA = 1:numel(ARRAY)
                plot.emg_averages(...
                    SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), pars);
            end
        end
    else
        fig = gobjects(numel(BLOCK), numel(ARRAY)); % Preallocate figure handle array.
        
        for iB = 1:numel(BLOCK)
            for iA = 1:numel(ARRAY)
                fig(iB, iA) = plot.emg_averages(...
                    SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), pars);
            end
        end
        
    end
    return;
end

pars.EMG_Type = fixCase(pars.EMG_Type);
switch pars.EMG_Type
    case 'Array'
        fig = plot.emg_averages__unipolar_array(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, pars);
    case 'Bipolar'
        fig = plot.emg_averages__bipolar(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, pars);
    case 'Accelerometer'
        fig = plot.emg_averages__imu(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, pars);
    otherwise
        error("I have not yet set up handling for EMG_Type == <strong>%s</strong>.\n", pars.EMG_Type);
end

% If no figure then do not try to save regardless.
if isa(fig, 'matlab.graphics.GraphicsPlaceholder')
    return;
end

if nargout < 1
    tank = sprintf('%s_%04d_%02d_%02d', SUBJ, YYYY, MM, DD); % data "tank"
    block = sprintf('%s_%s_%d', tank, ARRAY, BLOCK); % experimental "block" (recording within tank)
    out_folder = fullfile(pars.Output_Root, SUBJ, tank, 'figures', 'Averages', pars.Filtering.Name);
    if exist(out_folder, 'dir') == 0
        try %#ok<TRYNC>
            mkdir(out_folder);
        end
    end
    out_name = fullfile(out_folder, sprintf('%s_%d_%d', block, round(pars.T(1)), round(pars.T(2))));
    default.savefig(fig, out_name, sprintf("Mean_%s_EMG", pars.EMG_Type), true);
    
    out_folder_2 = fullfile(pars.Output_Root, SUBJ, tank, num2str(BLOCK));
    if exist(out_folder_2, 'dir') == 0
        try %#ok<TRYNC>
            mkdir(out_folder_2);
        end
    end
    default.savefig(fig, fullfile(out_folder_2, block), sprintf('%d_%d_Mean_%s_EMG_%s', round(pars.T(1)), round(pars.T(2)), pars.EMG_Type, pars.Filtering.Name), false); 
end

    
end
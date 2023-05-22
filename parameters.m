function varargout = parameters(varargin)
%PARAMETERS Return parameters struct, which sets default values for things like epoch durations etc.
%
% Example 1
%   pars = plot.parameters('emg_averages'); % Returns parameters for 'emg_averages'
%
% See also: Contents

pars = struct;
VERSION = "2.5.0";

%% For emg_averages
pars.emg_averages = struct;
pars.emg_averages.Acquisition_Type = "TMSi";
pars.emg_averages.Anonymize = true;
pars.emg_averages.Blank_Stim = false;
pars.emg_averages.Data = []; % If provided, skips re-loading the Data objects
pars.emg_averages.EMG_Type = "Array"; % Can be: "Array" | "Bipolar"
pars.emg_averages.End_Linear_Fit = inf; % Set to some value in milliseconds to say when linear fit should end
pars.emg_averages.Filtering = utils.get_default_filtering_pars(pars.emg_averages.Acquisition_Type, pars.emg_averages.EMG_Type, "Rectified"); % Return default filtering struct
pars.emg_averages.File_Type = ".mat"; % Can be: ".mat" | ".poly5"
pars.emg_averages.Filtered_Data = []; % If this is provided, skips the `utils.apply_emg_filters` step.
pars.emg_averages.Inverted_Logic = true; % Set true to indicate that sync bit logic is inverted (default) or false if it is non-inverted.
pars.emg_averages.Linear_Fit_Order = 1; % For Subtract_Linear_Fit polynomial detrend that is applied to individual trials (DIFFERENT FROM APPLY_POLYNOMIAL_DETREND IN APPLY_EMG_FILTERS)
pars.emg_averages.Link_Axes = false;  % Set false to al axes limits to adaptively set independently
pars.emg_averages.N_SD_RMS = 3.5; % Number of times the RMS to multiply signal by when computing YLIM if it is not manually specified.
pars.emg_averages.N_Individual_Max = 10; % Max. number of individual traces to superimpose
pars.emg_averages.N_Rows = nan; % Number of rows in grid layout
pars.emg_averages.N_Columns = nan; % Number of columns in grid layout
pars.emg_averages.N_Trials = nan;
[pars.emg_averages.Output_Root, pars.emg_averages.Input_Root] = parameters('generated_data_folder', 'raw_data_folder'); % Location where output figures are saved.
pars.emg_averages.Plot_Stim_Period = true; % Plot stim artifact with red stem lines?
pars.emg_averages.SNR_Sort = false; % rearrange the order of trials so responses with highest amplitudes are first
pars.emg_averages.Style = "Shaded"; % Can be: "Shaded" | "Individual"
pars.emg_averages.Start_Linear_Fit = 4.75; % Start the linear fit subtraction (if Subtract_Linear_Fit is true) on or after the sample corresponding to this value (ms)
pars.emg_averages.Subtract_Linear_Fit = false; % Set false to skip the linear-fit subtraction.
pars.emg_averages.Subtract_Mean = false; % Subtract mean?
pars.emg_averages.Sync_Bit = nan; % The bit address for STIM sync TTL signal on TRIGGERS channel of TMSi.
pars.emg_averages.T = [-60, 60]; % Time for epochs (milliseconds)
pars.emg_averages.T_RMS = [12, 60]; % Time epoch for computing RMS
pars.emg_averages.Trigger_Channel = 'TRIGGER'; % Name of Trigger Channel
pars.emg_averages.XLim = []; % If empty, use auto-scale, otherwise, fixed scale
pars.emg_averages.YLim = []; % If empty, use auto-scale, otherwise, fixed scale
pars.emg_averages.Version = VERSION;
pars.emg_averages.Verbose = false;
pars.emg_averages.Process_Steps = []; % If empty, use original filtering method, otherwise, use ordered filtering steps according to single char element or cell array

%% For emg_stack
pars.emg_stack = struct;
pars.emg_stack.Axes = [];
pars.emg_stack.Acquisition_Type = "TMSi";
pars.emg_stack.Align_Peaks = false;
pars.emg_stack.Anonymize = true;
pars.emg_stack.Blank_Stim = false;
pars.emg_stack.Data = [];
pars.emg_stack.EMG_Filters_Applied = false; 
pars.emg_stack.EMG_Type = "Array"; % Can be: "Array" | "Bipolar"
pars.emg_stack.Filtering = utils.get_default_filtering_pars(pars.emg_stack.Acquisition_Type, pars.emg_stack.EMG_Type, "Rectified"); % Return default filtering struct
pars.emg_stack.File_Type = ".mat"; % Can be: ".mat" | ".poly5"
pars.emg_stack.Filtered_Data = []; % If this is provided, skips the `utils.apply_emg_filters` step.
pars.emg_stack.Font = {'FontName', 'Tahoma', 'FontSize', 18, 'Color', 'k'};
pars.emg_stack.Inverted_Logic = true; % Set true to indicate that sync bit logic is inverted (default) or false if it is non-inverted.
pars.emg_stack.N_Rows = nan; % Number of rows in grid layout
pars.emg_stack.N_Columns = nan; % Number of columns in grid layout
pars.emg_stack.N_Trials = 30; % Number of trials (max)
[pars.emg_stack.Output_Root, pars.emg_stack.Input_Root] = parameters('generated_data_folder', 'raw_data_folder'); % Location where output figures are saved.
pars.emg_stack.Scale_Factor = 3.5;
pars.emg_stack.Series = [];
pars.emg_stack.Subtract_Linear_Fit = false; % Set false to skip the linear-fit subtraction.
pars.emg_stack.Subtract_Mean = false; % Subtract mean?
pars.emg_stack.Sync_Bit = nan; % The bit address for STIM sync TTL signal on TRIGGERS channel of TMSi.
pars.emg_stack.T = [-60, 60]; % Time for epochs (milliseconds)
pars.emg_stack.Trigger_Channel = 'TRIGGER'; % Name of Trigger Channel
pars.emg_stack.XLim = []; % If empty, use auto-scale, otherwise, fixed scale
pars.emg_stack.Version = VERSION;
pars.emg_stack.Verbose = false;
pars.emg_stack.Process_Steps = []; % If empty, use original filtering method, otherwise, use ordered filtering steps according to single char element or cell array

%% For emg_rms
pars.emg_rms = struct;
pars.emg_rms.Acquisition_Type = "TMSi";
pars.emg_rms.Anonymize = true;
pars.emg_rms.Blank_Stim = false;
pars.emg_rms.EMG_Type = "RMS"; % Can be: "Array" | "Bipolar" | "RMS"
pars.emg_rms.Data = [];
pars.emg_rms.Debug = false;
pars.emg_rms.File_Type = ".mat"; % Can be ".mat" | ".poly5"
pars.emg_rms.Filtering = utils.get_default_filtering_pars(pars.emg_rms.Acquisition_Type, pars.emg_rms.EMG_Type, "Rectified"); % Return default filtering struct
pars.emg_rms.Filtered_Data = []; % If this is provided, skips the `utils.apply_emg_filters` step.
pars.emg_rms.Axes = [];
[pars.emg_rms.Output_Root, pars.emg_rms.Input_Root] = parameters('generated_data_folder', 'raw_data_folder'); % Location where output figures are saved.
pars.emg_rms.N_Trials = nan;
pars.emg_rms.Pre_Stimulus_RMS_ms = -20; % Number of milliseconds relative to stimulus that is the latest the pre-stimulus window goes.
pars.emg_rms.Post_Stimulus_RMS_ms = 20; % Number of milliseconds after stimulus before computing RMS window
pars.emg_rms.RMS_Response_Ratio_Threshold = []; % Minimum amount to consider putting channel name on contour map
pars.emg_rms.RMS_Max_Response_Ratio = []; % Clip color scale above this value.
pars.emg_rms.Subtract_Mean = false; % Subtract cross-trial mean?
pars.emg_rms.Sync_Bit = nan; % The bit address for STIM sync TTL signal on TRIGGERS channel of TMSi.
pars.emg_rms.T = [-30, 30]; % Time for epochs (milliseconds)
pars.emg_rms.Trigger_Channel = 'TRIGGER'; % Name of Trigger Channel
pars.emg_rms.XLim = []; % If empty, use auto-scale, otherwise, fixed scale
pars.emg_rms.YLim = []; % If empty, use auto-scale, otherwise, fixed scale
pars.emg_rms.CLim = []; % If empty, use autoscale, otherwise, fixed scale
pars.emg_rms.Version = VERSION;
pars.emg_rms.Verbose = false;
pars.emg_rms.Process_Steps = []; % If empty, use original filtering method, otherwise, use ordered filtering steps according to single char element or cell array

%% For emg_waterfall
pars.emg_waterfall = struct;
pars.emg_waterfall.Acquisition_Type = "TMSi";
pars.emg_waterfall.Align_Peaks = false;
pars.emg_waterfall.Anonymize = true;
pars.emg_waterfall.Axes = [];
pars.emg_waterfall.Blank_Stim = false;
pars.emg_waterfall.Colormap = 'spring';
pars.emg_waterfall.Data = [];
pars.emg_waterfall.Data_Channel = nan;
pars.emg_waterfall.EMG_Filters_Applied = false;
pars.emg_waterfall.EMG_Type = "Array"; % Can be: "Array" | "Bipolar"
pars.emg_waterfall.Figure_Title = 'Waterfall';
pars.emg_waterfall.File_Type = ".mat"; % Can be: ".mat" | ".poly5"
pars.emg_waterfall.Filtering = utils.get_default_filtering_pars("TMSi","Array","Raw", ...
    'Apply_Virtual_Reference',true,"Apply_HPF",true,"HPF_Cutoff_Frequency",2,"HPF_Order",2); % Return default filtering struct
pars.emg_waterfall.Filtered_Data = []; % If this is provided, skips the `utils.apply_emg_filters` step.
pars.emg_waterfall.Font = {'FontName', 'Tahoma', 'FontSize', 18, 'Color', 'k'};
pars.emg_waterfall.Inverted_Logic = false;
% pars.emg_waterfall.N_Individual_Max = 10; % Max. number of individual traces to superimpose
pars.emg_waterfall.N_Trials = 20;
[pars.emg_waterfall.Output_Root, pars.emg_waterfall.Input_Root] = parameters('generated_data_folder', 'raw_data_folder');
pars.emg_waterfall.Plot_Stim_Period = true; % Plot stim artifact with red stem lines?
pars.emg_waterfall.Sample_Rate = 4000;    % Sample rate from acquisition
pars.emg_waterfall.Subtract_Mean = false; % Subtract cross-trial mean?
pars.emg_waterfall.Sync_Bit = nan;        % The bit address for STIM sync TTL signal on TRIGGERS channel of TMSi.
pars.emg_waterfall.T = [-60, 60];     % Time for epochs (milliseconds)
pars.emg_waterfall.T_RMS = [30, 60];  % Time epoch for computing RMS
pars.emg_waterfall.Trigger_Data = [];
pars.emg_waterfall.Trigger_Channel = 'TRIGGER';
pars.emg_waterfall.View = [95, 65]; 
pars.emg_waterfall.Version = VERSION;
pars.emg_waterfall.Process_Steps = []; % If empty, use original filtering method, otherwise, use ordered filtering steps according to single char element or cell array

% % % Color limits as well as axes limits % % %
pars.emg_waterfall.C_Lim = []; % If empty, use auto-scale, otherwise, fixed scale
pars.emg_waterfall.X_Lim = []; % If empty, use auto-scale, otherwise, fixed scale
pars.emg_waterfall.Y_Lim = []; % If empty, use auto-scale, otherwise, fixed scale
pars.emg_waterfall.Z_Lim = []; % If empty, use auto-scale, otherwise, fixed scale.
pars.emg_waterfall.Verbose = false;

%% For impedances
% Impedance parameters
pars.impedance.FileTag = '';
pars.impedance.CLim = [0 150]; % kOhms
pars.impedance.AxLim = [0.5 8.5];
pars.impedance.Colormap = cm.map('greenred'); % Like TMSi one
pars.impedance.Colorscale = 'log';
pars.impedance.Version = VERSION;
pars.impedance.Verbose = false;
pars.impedance.Figure_Process_Title = nan; % The substring in the figure title that displays the filtering processes used (accepts manual overwrite)

%% Handle parsing
N = numel(varargin);
if nargout == 1
    if rem(N, 2) == 1
        varargout = {pars.(varargin{end})};
        return;
    else
        f = fieldnames(pars);
        for iV = 1:2:N
            idx = strcmpi(f, varargin{iV});
            if sum(idx) == 1
               pars.(f{idx}) = varargin{iV+1}; 
            end
        end
        varargout = {pars};
        return;
    end
else
    f = fieldnames(pars);
    varargout = cell(1, nargout);
    for iV = 1:numel(varargout)
        idx = strcmpi(f, varargin{iV});
        if sum(idx) == 1
            varargout{iV} = pars.(f{idx}); 
        else
            error('Could not find parameter: %s', varargin{iV}); 
        end
    end
end

end

function fig = raw_tmsi_channel(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, CHANNEL, varargin)
%RAW_TMSI_CHANNEL Handle plotting single-channel raw TMSi data.
%
% Syntax:
%   fig = plot.raw_tmsi_channel(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, CHANNEL);
%   fig = plot.raw_tmsi_channel(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, CHANNEL, 'Name', value, ...);
%
% Example 1:
%   % Run channel 24 from block 43, using default parameters.
%   fig = plot.raw_tmsi_channel('Frank', 2021, 11, 18, "A", 43, 24);
%
% Inputs:
%   SUBJ  - Subject name (e.g. 'Frank' or "Frank")
%   YYYY  - Year (numeric or string, e.g. 2021 or "2021" or '2021')
%   MM    - Month (numeric or string, e.g. 11 or "11" or '11')
%   DD    - Day   (numeric or string, e.g. 18 or "18" or '18')
%   ARRAY - "A" or "B" (or 'A' or 'B') or ["A", "B"] or {'A', 'B'}
%   BLOCK - Recording parameter key (block; numeric or string, e.g. 0 or "0")
%   CHANNEL - The 1-indexed channel to record from (e.g. 24).
%   varargin - (Optional) 'Name', value parameter field/value pairs. See
%                       `pars` struct for specific parameter names.
%
%   If BLOCK is an array, then the returned output will be an array of
%   figure handles with a row for each element of BLOCK. If ARRAY is a
%   (MATLAB) array, then the returned output will have a column for each
%   element of ARRAY. For recording blocks that do not have any data or
%   block associated to a particular ARRAY, the corresponding figure handle
%   object will be an empty `gobjects` element. The third dimension of the
%   returned output figure handle array refers to matched elements of
%   CHANNEL, if CHANNEL is given as a vector.
%
% Output:
%   fig   - Figure handle. If no output is requested, figures are
%           automatically saved and deleted using `default.savefig`. The
%           output saved results are in `generated_data` and follow the
%           naming convention provided by the input arguments.

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

pars = struct;
pars.Acquisition_Type = "TMSi";
pars.EMG_Type = "Array"; % Can be: "Array" | "Bipolar"
pars.File_Type = ".mat"; % Can be: ".mat" | ".poly5"
pars.Filtering = utils.get_default_filtering_pars(); % Return default filtering struct
pars.N_SD_RMS = 2.5; % Number of times the RMS to multiply signal by when computing YLIM if it is not manually specified.
pars.N_Individual_Max = 10; % Max. number of individual traces to superimpose
[pars.Output_Root, pars.Input_Root] = parameters('generated_data_folder', 'raw_data_folder'); % Location where output figures are saved.
pars.Sync_Bit = 9; % The bit address for STIM sync TTL signal on TRIGGERS channel of TMSi.
pars.T = [-15, 80]; % Time for epochs (milliseconds)
pars.T_RMS = [30, 60]; % Time epoch for computing RMS
pars.Trigger_Channel = 'TRIGGER'; % Name of Trigger Channel
pars.XLim = []; % If empty, use auto-scale, otherwise, fixed scale
pars.YLim = []; % If empty, use auto-scale, otherwise, fixed scale

% % % END DEFAULT PARS STRUCT FIELD DEFINITIONS % % %

% Handle parsing of `pars`
pars = utils.parse_parameters(pars, varargin{:});
if ~isstruct(pars.Filtering)
     pars.Filtering = utils.get_default_filtering_pars(pars.Acquisition_Type, pars.EMG_Type, pars.Filtering);
end

if (numel(BLOCK) > 1) || (numel(ARRAY) > 1) || (numel(CHANNEL) > 1)
    if nargout < 1  % Then ensure we do not change "nargout" for actual implementation.
        for iB = 1:numel(BLOCK)
            for iA = 1:numel(ARRAY)
                for iC = 1:numel(CHANNEL)
                    plot.raw_tmsi_channel(...
                        SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), CHANNEL(iC), pars);
                end
            end
        end
    else
        fig = gobjects(numel(BLOCK), numel(ARRAY), numel(CHANNEL)); % Preallocate figure handle array.
        
        for iB = 1:numel(BLOCK)
            for iA = 1:numel(ARRAY)
                for iC = 1:numel(CHANNEL)
                    fig(iB, iA, iC) = plot.raw_tmsi_channel(...
                        SUBJ, YYYY, MM, DD, ARRAY(iA),BLOCK(iB), CHANNEL(iC), pars);
                end
            end
        end
        
    end
    return;
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
fig = default.figure(block, 'Position', [0.1 0.1 0.8 0.8]);
fig.UserData = struct('x', x, 'version', 2.0); % Associate thse data to the figure.


end

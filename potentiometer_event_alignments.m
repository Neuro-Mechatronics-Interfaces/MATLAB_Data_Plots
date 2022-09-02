function fig = potentiometer_event_alignments(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, ALIGNMENT, varargin)
%POTENTIOMETER_EVENT_ALIGNMENTS Handle plotting potentiometers from wrist task against parsed events for a given exported alignment.
%
% Syntax:
%   fig = plot.potentiometer_event_alignments(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, ALIGNMENT);
%   plot.potentiometer_event_alignments(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, ALIGNMENT, 'Name', value, ...);
%
% Example 1:
%   % Run channel 24 from block 43, using default parameters.
%   fig = plot.potentiometer_event_alignments('Spencer', 2021, 8, 15, "A", 1, "MOVE_THRESHOLD");
%
% Inputs:
%   SUBJ  - Subject name (e.g. 'Frank' or "Frank")
%   YYYY  - Year (numeric or string, e.g. 2021 or "2021" or '2021')
%   MM    - Month (numeric or string, e.g. 11 or "11" or '11')
%   DD    - Day   (numeric or string, e.g. 18 or "18" or '18')
%   ARRAY - "A" or "B" (or 'A' or 'B') or ["A", "B"] or {'A', 'B'}
%   BLOCK - Recording parameter key (block; numeric or string, e.g. 0 or "0")
%   ALIGNMENT - "MOVE_ONSET" | "MOVE_THRESHOLD" | "OVERSHOOT_ONSET" | "REWARD_ONSET"
%   varargin - (Optional) 'Name', value parameter field/value pairs. See
%                       `pars` struct for specific parameter names.
%
% Output:
%   fig   - Figure handle. If no output is requested, figures are
%           automatically saved and deleted using `default.savefig`. The
%           output saved results are in `generated_data` and follow the
%           naming convention provided by the input arguments.
%
% See also: Contents, io.load_tmsi_potentiometers

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
    [pars.Output_Root, pars.Input_Root, pars.Version] = parameters('generated_data_folder', 'raw_data_folder', 'version'); % Location where output figures are saved.
    pars.XLim = []; % If empty, use auto-scale, otherwise, fixed scale
    pars.YLim = []; % If empty, use auto-scale, otherwise, fixed scale
    % Handle parsing of `pars`
    pars = utils.parse_parameters(pars, varargin{:});
end
% % % END DEFAULT PARS STRUCT FIELD DEFINITIONS % % %

ALIGNMENT = string(ALIGNMENT);

if (numel(BLOCK) > 1) || (numel(ARRAY) > 1) || (numel(ALIGNMENT) > 1)
    if nargout < 1  % Then ensure we do not change "nargout" for actual implementation.
        for iB = 1:numel(BLOCK)
            for iA = 1:numel(ARRAY)
                for iC = 1:numel(ALIGNMENT)
                    plot.potentiometer_event_alignments(...
                        SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), ALIGNMENT(iC), pars);
                end
            end
        end
    else
        fig = gobjects(numel(BLOCK), numel(ARRAY), numel(CHANNEL)); % Preallocate figure handle array.
        
        for iB = 1:numel(BLOCK)
            for iA = 1:numel(ARRAY)
                for iC = 1:numel(CHANNEL)
                    fig(iB, iA, iC) = plot.potentiometer_event_alignments(...
                        SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), ALIGNMENT(iC), pars);
                end
            end
        end
        
    end
    return;
end

% Load in potentiometer data
data = io.load_tmsi_potentiometers(SUBJ, YYYY, MM, DD, ARRAY, BLOCK);

% % % TODO: Add the rest of this! % % %


fig = default.figure(block, 'Position', [0.1 0.1 0.8 0.8]);
ax = default.axes(fig, 'XLabel', "Time (s)", 'YLabel', "Angle (\circ)");

fig.UserData = struct('function', 'plot.potentiometer_event_alignments', ...
                      'version', pars.Version); % Associate these data to the figure.


end

function fig = target_classification_tracking(T, subj, m)
%TARGET_CLASSIFICATION_TRACKING Plot tracking through time for F1-score
%
% Syntax:
%   fig = plot.trial_durations(T, subj, m);
%
% Inputs:
%   T - Table with trial data.
%   SUBJ - Name of experimental subject.
%   m - string month (e.g. "October")
%
% Output:
%   fig - Figure handle
%
% See also: Contents

if nargin < 2
    subj = "Spencer"; 
end

if nargin < 3
    m = "October"; 
end

if subj == "Spencer"
    t_subj = "Monkey-S";
else
    t_subj = "Monkey-R";
end

fig = figure('Name', 'Target Classification', 'Color', 'w'); 
ax = axes(fig, 'NextPlot', 'add', 'XColor', 'k', 'YColor', 'k', 'FontName', 'Tahoma');
plot(ax, T.Date(T.Subject==subj), T.Macro_F1(T.Subject==subj), ...
    'DisplayName', 'Macro-F1', ...
    'LineWidth', 1.5);
plot(ax, T.Date(T.Subject==subj), T.Macro_Recall(T.Subject==subj), ...
    'DisplayName', 'Macro-Recall', ...
    'LineWidth', 1.5);
plot(ax, T.Date(T.Subject==subj), T.Macro_Precision(T.Subject==subj), ...
    'DisplayName', 'Macro-Precision', ...
    'LineWidth', 1.5);
legend(ax, 'Location', 'southwest', 'TextColor', 'black', ...
    'FontName', 'Tahoma');
title(ax, sprintf('%s %s Training', t_subj, m), 'Color', 'k', ...
    'FontName', 'Tahoma');
ylim(ax, [0 1]);
end
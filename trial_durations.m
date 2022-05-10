function fig = trial_durations(T, SUBJ, m)
%TRIAL_DURATIONS Plot bar graphs showing distribution of trial durations.
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
    SUBJ = "Spencer"; 
end

if nargin < 3
    m = "October"; 
end

if SUBJ == "Spencer"
    t_subj = "Monkey-S";
else
    t_subj = "Monkey-R";
end

T = T(T.Subject == SUBJ, :);
fig = figure('Name', 'Wrist Task Trial Durations', 'Color', 'w'); 
ax = axes(fig, 'NextPlot', 'add', 'XColor', 'k', 'YColor', 'k', 'FontName', 'Tahoma');
yyaxis left;
ylabel(ax, 'Duration (s)');
ylim([0 3]);
plot(ax, T.Date, T.Median_Success_Duration, 'DisplayName', 'Trial Duration', 'LineWidth', 1.5);
yyaxis right;
ylabel(ax, 'E[Mistakes]');
ylim([0 1]);
plot(ax, T.Date, T.Mean_Success_Overshoots, 'DisplayName', 'Mean # Mistakes', 'LineWidth', 1.5);

title(ax, sprintf('%s %s Training', t_subj, m), 'Color', 'k', ...
    'FontName', 'Tahoma');

end

function fig = response_curves(S)
%RESPONSE_CURVES Creates set of response curves for parsed "response" metric data.
%
% Syntax:
%   fig = response_curves(S);
%
% Inputs:
%   S - Table with columns: ["Power_Ratio", "Frequency", "Amplitude", "Width", "Channel"]
% Output:
%   fig - Figure handle
%
% See also: Contents

[uFreq, matched, iU] = unique(S.Frequency);
uWidth = S.Width(matched);
C = rand(numel(uFreq), 3);
s = strcat(num2str(uFreq), " Hz | ", num2str(uWidth), " ms");

fig = figure(...
    'Name', 'Response Curves', ...
    'Color', 'w', ...
    'Units','Normalized',...
    'Position', [0.1 0.1 0.4 0.8]);
L = tiledlayout(fig, 3, 1);
ax = nexttile(L);
set(ax, 'NextPlot', 'add', 'XColor', 'k', 'YColor', 'k', ...
    'FontName', 'Tahoma', 'FontSize', 10);
colormap(ax, C);
for ii = 1:numel(uFreq)
    scatter(ax, S.Frequency(iU == ii), S.Power_Ratio(iU == ii), 'MarkerFaceAlpha', 0.25, ...
        'CData', iU(iU == ii), 'MarkerFaceColor', 'flat', 'DisplayName', s(ii));
end
xlabel(ax, 'Frequency (Hz)', 'FontName', 'Tahoma', 'Color', 'k');
ylabel(ax, 'Ratio', 'FontName', 'Tahoma', 'Color', 'k');
legend(ax, 'Location', 'northeast');

ax = nexttile(L);
set(ax, 'NextPlot', 'add', 'XColor', 'k', 'YColor', 'k', ...
    'FontName', 'Tahoma', 'FontSize', 10);
colormap(ax, C);
scatter(ax, S.Width, S.Power_Ratio, 'MarkerFaceAlpha', 0.25,...
    'CData', iU, 'MarkerFaceColor', 'flat');
xlabel(ax, 'Pulse Width (ms)', 'FontName', 'Tahoma', 'Color', 'k');
ylabel(ax, 'Ratio', 'FontName', 'Tahoma', 'Color', 'k');

ax = nexttile(L);
set(ax, 'NextPlot', 'add', 'XColor', 'k', 'YColor', 'k', ...
    'FontName', 'Tahoma', 'FontSize', 10);

uAmp = unique(S.Amplitude);
n = 20;
test = S(1:n, :);
test.Amplitude = linspace(min(uAmp), max(uAmp), n)';
for ii = 1:numel(uFreq)
    idx = S.Frequency == uFreq(ii);
    test.Frequency = ones(n, 1) * uFreq(ii);
    test.Width = ones(n, 1) * uWidth(ii);
    h = scatter(ax, S.Amplitude(idx), S.Power_Ratio(idx), ...
        'MarkerEdgeColor', C(ii, :), ...
        'MarkerFaceColor', C(ii, :), ...
        'MarkerFaceAlpha', 0.25, ...
        'DisplayName', s(ii));
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

xlabel(ax, 'Current Amplitude (mA)', 'FontName', 'Tahoma', 'Color', 'k');
ylabel(ax, 'Ratio', 'FontName', 'Tahoma', 'Color', 'k');
title(L, 'Per-Channel Post-Stim to Pre-Stim RMS Ratios', ...
    'FontName', 'Tahoma', 'Color', 'k');
end

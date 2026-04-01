function h = add_groupwise_sig_brackets(ax, pairIdx, labels, xCenters, groupTopY, opts)
%ADD_GROUPWISE_SIG_BRACKETS Draw non-overlapping significance brackets.
%
% h = plot.add_groupwise_sig_brackets(ax, pairIdx, labels, xCenters, groupTopY)
%
% Inputs
%   ax        : target axes
%   pairIdx   : Nx2 matrix of group indices, e.g. [1 2; 1 3; 2 4]
%   labels    : Nx1 string/cellstr with text to show above each bracket
%   xCenters  : 1xG or Gx1 vector of x-center positions for each group
%   groupTopY : 1xG or Gx1 vector giving current top y per group
%
% Name-value options
%   opts.BaseGapFrac     : gap above tallest data in axis-fraction units
%   opts.LevelGapFrac    : vertical spacing between bracket levels
%   opts.TickHeightFrac  : bracket vertical tick height
%   opts.FontSize        : text font size
%   opts.LineWidth       : bracket line width
%   opts.Color           : bracket/text color
%   opts.FontWeight      : text font weight
%
% Notes
%   - Brackets are automatically stacked to avoid overlap when x-ranges
%     intersect.
%   - The function expands the y-limits if needed.

arguments
    ax (1,1) matlab.graphics.axis.Axes
    pairIdx (:,2) double
    labels
    xCenters (:,1) double
    groupTopY (:,1) double
    opts.BaseGapFrac (1,1) double = 0.04
    opts.LevelGapFrac (1,1) double = 0.06
    opts.TickHeightFrac (1,1) double = 0.015
    opts.FontSize (1,1) double = 11
    opts.LineWidth (1,1) double = 1.2
    opts.Color = [0 0 0]
    opts.FontWeight char {mustBeMember(opts.FontWeight,{'normal','bold'})} = 'bold'
end

if isempty(pairIdx)
    h = gobjects(0);
    return
end

labels = string(labels);
xCenters  = xCenters(:);
groupTopY = groupTopY(:);

holdState = ishold(ax);
hold(ax,'on');

yl = ylim(ax);
yRange = yl(2) - yl(1);
if yRange <= 0
    yRange = 1;
end

baseGap    = opts.BaseGapFrac    * yRange;
levelGap   = opts.LevelGapFrac   * yRange;
tickHeight = opts.TickHeightFrac * yRange;

nPairs = size(pairIdx,1);

% Convert group indices to x-ranges and sort by span, then by left edge
leftX  = nan(nPairs,1);
rightX = nan(nPairs,1);
spanX  = nan(nPairs,1);
baseY  = nan(nPairs,1);

for i = 1:nPairs
    i1 = pairIdx(i,1);
    i2 = pairIdx(i,2);

    x1 = xCenters(i1);
    x2 = xCenters(i2);

    leftX(i)  = min(x1,x2);
    rightX(i) = max(x1,x2);
    spanX(i)  = rightX(i) - leftX(i);

    baseY(i) = max(groupTopY(min(i1,i2):max(i1,i2))) + baseGap;
end

[~, order] = sortrows([spanX, leftX], [1 2]);

% Assign stacking levels so overlapping x-ranges do not share a level
levelAssigned = nan(nPairs,1);
levelRanges = {};  % each cell contains [left right] rows already in that level

for k = 1:nPairs
    i = order(k);
    placed = false;

    for level = 1:numel(levelRanges)
        existing = levelRanges{level};
        overlaps = any(leftX(i) <= existing(:,2) & rightX(i) >= existing(:,1));
        if ~overlaps
            levelAssigned(i) = level;
            levelRanges{level}(end+1,:) = [leftX(i), rightX(i)];
            placed = true;
            break
        end
    end

    if ~placed
        levelAssigned(i) = numel(levelRanges) + 1;
        levelRanges{end+1} = [leftX(i), rightX(i)];
    end
end

% Draw brackets
h = gobjects(nPairs,2); % line + text
maxNeededY = yl(2);

for i = 1:nPairs
    i1 = pairIdx(i,1);
    i2 = pairIdx(i,2);

    x1 = xCenters(i1);
    x2 = xCenters(i2);

    y = baseY(i) + (levelAssigned(i)-1)*levelGap;

    % bracket
    h(i,1) = plot(ax, ...
        [x1 x1 x2 x2], ...
        [y y+tickHeight y+tickHeight y], ...
        'Color', opts.Color, ...
        'LineWidth', opts.LineWidth);

    % label
    h(i,2) = text(ax, mean([x1 x2]), y + tickHeight, labels(i), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', opts.FontSize, ...
        'FontWeight', opts.FontWeight, ...
        'Color', opts.Color);

    maxNeededY = max(maxNeededY, y + tickHeight + 0.04*yRange);
end

if maxNeededY > yl(2)
    ylim(ax, [yl(1), maxNeededY]);
end

if ~holdState
    hold(ax,'off');
end
end
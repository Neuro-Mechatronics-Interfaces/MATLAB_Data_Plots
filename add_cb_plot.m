function [fig,ax,hg,hLeg] = add_cb_plot(x, y, opts)
%ADD_CB_PLOT Adds confidence-band with centered-line and corresponding legend entries.
%
% Syntax:
%   [fig,ax,hg,hLeg] = plot.add_cb_plot(x, y, 'Name', value, ...);
%
% Inputs:
%   x (:,1) double = linspace(0,2*pi,201)';
%   y (:,1) double = sin(x);
%
% Options:
%   'Axes' = [];
%   'YPatchOffset' (1,2) double = [-0.5, 0.5];
%   'YConfidenceBounds' (:,2) double = [nan, nan];
%   'PatchFaceColor' = [0.5 0.5 0.5];
%   'PatchFaceAlpha' (1,1) double = 0.75;
%   'PatchEdgeColor' = 'none';
%   'PatchEdgeAlpha' = 0.5;
%   'LineColor' (1,3) double = [0 0 0];
%   'LineWidth' (1,1) double = 2;
%   'LineStyle' {mustBeTextScalar} = '-';
%   'FontName' {mustBeTextScalar} = 'Arial';
%   'GroupName' {mustBeTextScalar} = 'My Group';
%   'AddLegend' (1,1) logical = true;
%
% Output:
%   fig  - Figure handle
%   ax   - Axes handle
%   hg   - Group object handle
%   hLeg - Custom Composite Legend Group handle

arguments
    x (:,1) double = linspace(0,2*pi,201)';
    y (:,1) double = sin(x);
    opts.Axes = [];
    opts.YPatchOffset (1,2) double = [-0.5, 0.5];
    opts.YConfidenceBounds (:,2) double = [nan, nan];
    opts.PatchFaceColor = [0.5 0.5 0.5];
    opts.PatchFaceAlpha (1,1) double = 0.25;
    opts.PatchEdgeColor = 'none';
    opts.PatchEdgeAlpha (1,1) double = 0.25;
    opts.LineColor (1,3) double = [0 0 0];
    opts.LineWidth (1,1) double = 2;
    opts.LineStyle {mustBeTextScalar} = '-';
    opts.FontName {mustBeTextScalar} = 'Arial';
    opts.GroupName {mustBeTextScalar} = 'My Group';
    opts.AddLegend (1,1) logical = true;
    opts.ResizeAxes (1,1) logical = true;
    opts.PatchWidthFrac (1,1) double = 0.15;  % width of patch icon
    opts.PatchHeightFrac (1,1) double  = 0.06;
    opts.TextOffsetFrac (1,1) double = 0.02; % gap between patch and text in x
    opts.LeftMarginFrac (1,1) double = 0.02;
    opts.RightMarginFrac (1,1) double = 0.02; % margin to right edge of axes
    opts.TopMarginFrac (1,1) double = 0.02; % margin to top edge of axes
    opts.BottomMarginFrac (1,1) double = 0.02;
    opts.EntryMarginFrac (1,1) double = 0.01;
    opts.EntryHeightFrac (1,1) double = 0.08;
    opts.HorizontalExpansionFrac (1,1) double = 0.5;
    opts.VerticalExpansionFrac (1,1) double = 0.5;
    opts.LegendLocation {mustBeMember(opts.LegendLocation,["northeast","north","northwest","west","southwest","south","southeast","east"])} = "northeast";
end

% -------------------------------------------------------------------------
% Axes / figure setup
% -------------------------------------------------------------------------
if isempty(opts.Axes)
    fig = figure('Color','w', ...
        'Name','Test HG-Legend', ...
        'Units','inches', ...
        'Position',[1 1 5 5]);
    ax = axes(fig, ...
        'NextPlot','add', ...
        'FontName',opts.FontName, ...
        'XColor','k', ...
        'YColor','k', ...
        'XLim',[min(x), max(x)]);
else
    ax = opts.Axes;
    fig = ancestor(ax,'figure');
end

% Make sure we don't accidentally clear things
hold(ax, 'on');

% -------------------------------------------------------------------------
% Main group + patch + line
% -------------------------------------------------------------------------
hg = hggroup(ax, 'DisplayName', opts.GroupName);
% Turn off built-in legend handling – we'll draw our own composite legend
hg.Annotation.LegendInformation.IconDisplayStyle = 'off';

faces = 1:(2*numel(x));
faces(end+1) = 1;

xx = [x; flipud(x)];

if size(opts.YConfidenceBounds,1) == numel(x)
    yy = [opts.YConfidenceBounds(:,1); flipud(opts.YConfidenceBounds(:,2))];
else
    if isempty(y)
        error("Cannot provide empty y with missing 'YConfidenceBounds'!");
    end
    yy = [y + opts.YPatchOffset(1); flipud(y) + opts.YPatchOffset(2)];
end

verts = [xx, yy];

% Patch
ph = patch('Faces',faces, ...
    'Vertices',verts, ...
    'EdgeColor',opts.PatchEdgeColor, ...
    'FaceColor',opts.PatchFaceColor, ...
    'FaceAlpha',opts.PatchFaceAlpha, ...
    'EdgeAlpha',opts.PatchEdgeAlpha, ...
    'Parent',hg);
% Ensure these don't get their own built-in legend entries
ph.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Line
if ~isempty(y)
    lh = plot(ax, x, y, ...
        'Color',opts.LineColor, ...
        'LineWidth',opts.LineWidth, ...
        'LineStyle',opts.LineStyle, ...
        'Parent',hg);
    lh.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

% -------------------------------------------------------------------------
% Custom "legend" inside the same axes, using a composite entry per group
% -------------------------------------------------------------------------
hLeg = findobj(ax.Children,'Tag','CustomCompositeLegendGroup');

if opts.AddLegend

    % ============================================
    % Determine legend location + expansion rules
    % ============================================
    location = opts.LegendLocation;

    % Defaults
    verticalDirection = -1;      % stack downward
    needYExpansion    = true;    % most locations require Y-expansion
    needXExpansion    = false;
    horizontalMode    = false;   % true only for north/south

    switch location
        case "northeast"
            anchorXFrac = 1 - opts.RightMarginFrac;
            anchorYFrac = 1 - opts.TopMarginFrac;
            verticalDirection = -1;
            needYExpansion = true;

        case "northwest"
            anchorXFrac = opts.LeftMarginFrac;
            anchorYFrac = 1 - opts.TopMarginFrac;
            verticalDirection = -1;
            needYExpansion = true;

        case "southeast"
            anchorXFrac = 1 - opts.RightMarginFrac;
            anchorYFrac = opts.BottomMarginFrac;
            verticalDirection = +1;
            needYExpansion = true;

        case "southwest"
            anchorXFrac = opts.LeftMarginFrac;
            anchorYFrac = opts.BottomMarginFrac;
            verticalDirection = +1;
            needYExpansion = true;

        case "east"
            anchorXFrac = 1 - opts.RightMarginFrac;
            anchorYFrac =  0.5;
            verticalDirection = -1;
            needYExpansion = false;
            needXExpansion = true;

        case "west"
            anchorXFrac = opts.LeftMarginFrac;
            anchorYFrac = 0.5;
            verticalDirection = -1;
            needYExpansion = false;
            needXExpansion = true;

        case "north"
            horizontalMode = true;       % horizontal centered layout
            anchorXFrac = 0.50;          % CENTERED
            anchorYFrac = 1 - opts.TopMarginFrac;
            verticalDirection = -1;
            needYExpansion = true;

        case "south"
            horizontalMode = true;       % horizontal centered layout
            anchorXFrac = 0.50;
            anchorYFrac = opts.BottomMarginFrac;
            verticalDirection = +1;
            needYExpansion = true;
    end

    % ======================================================
    % Retrieve/create legend group
    % ======================================================
    if isempty(hLeg)
        hLeg = hggroup(ax, ...
            'Tag','CustomCompositeLegendGroup', ...
            'HandleVisibility','on', ...
            'HitTest','off');
        uistack(hLeg, 'top');
    else
        hLeg = hLeg(1);
    end

    % Count existing entries
    existingEntries = findobj(hLeg.Children,'Type','hggroup',...
        'Tag','CustomCompositeLegendEntry');
    nEntries = numel(existingEntries);
    nTot     = nEntries + 1;

    % ======================================================
    % Axis expansion (X or Y depending on location)
    % ======================================================
    xlim = get(ax,'XLim');
    ylim = get(ax,'YLim');

    origYLim = ylim;
    origXLim = xlim;

    dx = diff(xlim);
    dy = diff(ylim);

    % Base metrics in data units
    marginH       = opts.EntryMarginFrac * dy;
    entryHeightH  = opts.EntryHeightFrac * dy;
    patchH        = opts.PatchHeightFrac * dy;

    % Total needed for vertical stacks
    totalLegendHeight = nTot * entryHeightH + 2*marginH;

    % Expand Y when needed
    if needYExpansion && opts.ResizeAxes
        if verticalDirection < 0   % expand upward
            ylim(2) = ylim(2) + totalLegendHeight;
        else                       % expand downward
            ylim(1) = ylim(1) - totalLegendHeight;
        end
        set(ax,'YLim',ylim);
    end

    % Expand X when needed (east/west)
    if needXExpansion && opts.ResizeAxes
        extraX = opts.HorizontalExpansionFrac * dx;
        if location == "east"
            xlim(2) = xlim(2) + extraX;
        else
            xlim(1) = xlim(1) - extraX;
        end
        set(ax,'XLim',xlim);
    end

    % Recompute after expansion
    xlim = get(ax,'XLim'); dx = diff(xlim);
    ylim = get(ax,'YLim'); dy = diff(ylim);

    % ======================================================
    % Compute anchor position in data coordinates
    % ======================================================
    anchorX = xlim(1) + anchorXFrac * dx;
    anchorY = ylim(1) + anchorYFrac * dy;

    % ======================================================
    % Create new entry group
    % ======================================================
    entryGroup = hggroup('Parent',hLeg, ...
        'Tag','CustomCompositeLegendEntry', ...
        'HitTest','off');

    % ======================================================
    % Vertical stack mode (corners + east/west)
    % ======================================================
    if ~horizontalMode

        % compute yCenter for this entry
        yCenter = anchorY + verticalDirection * (nEntries * entryHeightH);

        % horizontal geometry
        patchW = opts.PatchWidthFrac * diff(origXLim);
        textOffsetX = opts.TextOffsetFrac * dx;

        x1 = anchorX;
        x2 = anchorX + patchW;
        y1 = yCenter - patchH/2;
        y2 = yCenter + patchH/2;

        % patch
        patch(entryGroup,[x1 x2 x2 x1],[y1 y1 y2 y2],opts.PatchFaceColor,...
            'EdgeColor',opts.PatchEdgeColor,...
            'FaceAlpha',opts.PatchFaceAlpha,...
            'EdgeAlpha',opts.PatchEdgeAlpha,...
            'Clipping','off','HitTest','off');

        % line
        plot(entryGroup,[x1 x2],[yCenter yCenter],...
            'Color',opts.LineColor,...
            'LineWidth',opts.LineWidth,...
            'LineStyle',opts.LineStyle,...
            'Clipping','off','HitTest','off');

        % text
        t = text(entryGroup,x2+textOffsetX,yCenter,opts.GroupName,...
            'FontName',opts.FontName,...
            'HorizontalAlignment','left',...
            'VerticalAlignment','middle',...
            'Interpreter','none','Clipping','off','HitTest','off');

        % Spilloff correction (right side)
        ext = get(t,'Extent');
        xRight = ext(1) + ext(3);
        rightMarginX = xlim(2) - opts.RightMarginFrac*dx;
        if xRight > rightMarginX
            deltaX = rightMarginX - xRight;
            % shift entire entry
            kids = entryGroup.Children;
            for k=1:numel(kids)
                if isprop(kids(k),'XData')
                    kids(k).XData = kids(k).XData + deltaX;
                elseif isa(kids(k),'matlab.graphics.primitive.Text')
                    pos = kids(k).Position; pos(1)=pos(1)+deltaX; kids(k).Position=pos;
                end
            end
        end
        % ------------------------------------------------------------------
        % Enforce common left alignment for ALL vertical-stack text entries
        % ------------------------------------------------------------------
        allVertEntries = findobj(hLeg.Children, 'Tag','CustomCompositeLegendEntry');

        % Collect all text objects
        allText = [];
        for e = allVertEntries'
            txt = findobj(e.Children, 'Type','text');
            if ~isempty(txt)
                allText = [allText; txt];
            end
        end

        if ~isempty(allText)
            % Extract all left X positions
            allExt = arrayfun(@(tt)get(tt,'Extent'), allText, 'UniformOutput', false);
            allExt = cat(1, allExt{:});
            allLeftX = allExt(:,1);

            % The minimum left X becomes the aligned left boundary
            commonLeftX = min(allLeftX);

            % Shift each entry so its text starts at commonLeftX
            for e = allVertEntries'
                txt = findobj(e.Children, 'Type','text');
                if isempty(txt), continue; end

                ext = get(txt,'Extent');
                shiftX = commonLeftX - ext(1);

                kids = e.Children;
                for k = 1:numel(kids)
                    if isa(kids(k),'matlab.graphics.primitive.Text')
                        pos = kids(k).Position;
                        pos(1) = pos(1) + shiftX;
                        kids(k).Position = pos;
                    elseif isprop(kids(k),'XData')
                        kids(k).XData = kids(k).XData + shiftX;
                    end
                end
            end
        end

        uistack(hLeg,'top');
        return
    end

    % ======================================================
    % Horizontal centered mode (north/south)
    % ======================================================
    % Compute the icon widths for this entry
    patchW      = opts.PatchWidthFrac * origXLim;
    textOffsetX = opts.TextOffsetFrac * dx;

    % First create text invisibly to measure width
    tTmp = text(ax,0,0,opts.GroupName,'FontName',opts.FontName,...
        'Visible','off','Interpreter','none');
    extTxt = get(tTmp,'Extent'); textW = extTxt(3); delete(tTmp);

    entryTotalW = patchW + textOffsetX + textW;

    % Determine total width across ALL entries
    existingWidths = [];
    if ~isempty(existingEntries)
        for e = existingEntries'
            txt = findobj(e.Children,'Type','text');
            if isempty(txt), continue; end
            ext = get(txt,'Extent');
            txtW = ext(3);
            pw = opts.PatchWidthFrac * origXLim;
            to = opts.TextOffsetFrac * dx;
            existingWidths(end+1) = pw + to + txtW;
        end
    end
    allWidths = [existingWidths, entryTotalW];
    legendWide = sum(allWidths) + (numel(allWidths)-1)*opts.HorizontalSpacingFrac*dx;

    % Compute left boundary for centered layout
    leftX = anchorX - legendWide/2;

    % Starting X for this entry
    startX = leftX + sum(allWidths(1:end-1)) + ...
        (numel(allWidths)-1)*opts.HorizontalSpacingFrac*dx;

    % Vertical position
    yCenter = anchorY + verticalDirection*(entryHeightH/2);

    % Draw patch
    x1 = startX;
    x2 = startX + patchW;
    y1 = yCenter - patchH/2;
    y2 = yCenter + patchH/2;

    patch(entryGroup,[x1 x2 x2 x1],[y1 y1 y2 y2],opts.PatchFaceColor,...
        'EdgeColor',opts.PatchEdgeColor,...
        'FaceAlpha',opts.PatchFaceAlpha,...
        'EdgeAlpha',opts.PatchEdgeAlpha,...
        'Clipping','off','HitTest','off');

    % Draw line
    plot(entryGroup,[x1 x2],[yCenter yCenter],...
        'Color',opts.LineColor,...
        'LineWidth',opts.LineWidth,...
        'LineStyle',opts.LineStyle,...
        'Clipping','off','HitTest','off');

    % Draw text (left aligned beside this icon)
    t = text(entryGroup, x2+textOffsetX, yCenter, opts.GroupName,...
        'FontName',opts.FontName,...
        'HorizontalAlignment','left',...
        'VerticalAlignment','middle',...
        'Interpreter','none','Clipping','off','HitTest','off');

    % ======= ensure everything stays inside XLim =======
    ext = get(t,'Extent');
    xRight = ext(1) + ext(3);
    if xRight > xlim(2)
        deltaX = xlim(2) - xRight - opts.RightMarginFrac*dx;
        kids = entryGroup.Children;
        for k=1:numel(kids)
            if isprop(kids(k),'XData')
                kids(k).XData = kids(k).XData + deltaX;
            elseif isa(kids(k),'matlab.graphics.primitive.Text')
                pos = kids(k).Position; pos(1)=pos(1)+deltaX; kids(k).Position=pos;
            end
        end
    end

    uistack(hLeg,'top');
end


end
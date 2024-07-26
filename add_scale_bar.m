function h = add_scale_bar(ax, x0, y0, x1, y1, options)
%ADD_SCALE_BAR  Add scale bar to an existing axes
%
% Syntax:
%   h = plot.add_scale_bar(ax, x0, y0, x1, y1, 'Name', value, ...);

arguments
    ax
    x0 (1,1) % Starting point of horizontal bar
    y0 (1,1) % Starting point of vertical bar
    x1 (1,1) % Ending point of horizontal bar
    y1 (1,1) % Ending point of vertical bar
    options.Color (1,3) double {mustBeInRange(options.Color,0,1)} = [0.65 0.65 0.65];
    options.XUnits = 'sec';
    options.XLabelScaleFactor = 1;
    options.YUnits = 'μV';
    options.YLabelScaleFactor = 1;
    options.FontSize = 12;
    options.FontName {mustBeTextScalar} = 'Tahoma';
    options.XLabelRoundingLevel (1,1) {mustBeInteger} = 0;
    options.YLabelRoundingLevel (1,1) {mustBeInteger} = 0;
end
% Add scalebar to plot
dx = abs(x1 - x0);
if ~isduration(dx)
    dx = round(dx.*options.XLabelScaleFactor,options.XLabelRoundingLevel)./options.XLabelScaleFactor;
    tx = sprintf('%g %s', dx * options.XLabelScaleFactor, options.XUnits); % time scalebar text
else
    dx = seconds(round(seconds(dx).*options.XLabelScaleFactor,options.XLabelRoundingLevel)./options.XLabelScaleFactor);
    tx = sprintf('%g %s', seconds(dx) * options.XLabelScaleFactor, options.XUnits); % time scalebar text
end
x1 = x0 + dx;
dy = abs(y1 - y0);
if ~isduration(dy)
    dy = round(dy.*options.YLabelScaleFactor, options.YLabelRoundingLevel)./options.YLabelScaleFactor;
    ty = sprintf('%g %s', dy * options.YLabelScaleFactor, options.YUnits); % vertical scalebar text
else
    dy = seconds(round(seconds(dy).*options.YLabelScaleFactor, options.YLabelRoundingLevel)./options.YLabelScaleFactor);
    ty = sprintf('%g %s', seconds(dy) * options.YLabelScaleFactor, options.YUnits); % vertical scalebar text
end
y1 = y0 + dy;

h.XBar = line(ax, [x0, x0], [y0, y1], 'Color',options.Color,'LineWidth',1.25,'LineStyle','-');
h.YBar = line(ax, [x0, x1], [y0, y0], 'Color',options.Color,'LineWidth',1.25,'LineStyle','-');
h.YText = text(ax, x0-0.05*dx, y0+dy/2, ty, 'FontName',options.FontName,'HorizontalAlignment','center','VerticalAlignment','bottom','Rotation',90,'Color',options.Color,'FontSize',options.FontSize);
h.XText = text(ax, x0+dx/2, y0-0.05*dy, tx, 'FontName',options.FontName,'HorizontalAlignment','center','VerticalAlignment','top','Color',options.Color,'FontSize',options.FontSize);
h.XBar.Annotation.LegendInformation.IconDisplayStyle = 'off';
h.YBar.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
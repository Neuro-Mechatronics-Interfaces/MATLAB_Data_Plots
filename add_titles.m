function add_titles(ax, Title, Subtitle, options)
%ADD_TITLES  Adds title and subtitle to axes
%
% Syntax:
%   plot.add_titles(ax, Title, Subtitle, 'Name', value, ...);
arguments
    ax
    Title {mustBeText}
    Subtitle {mustBeText}
    options.FontName {mustBeTextScalar} = "Tahoma";
    options.FontSize (1,1) double = 14;
    options.TitleColor (1,3) double {mustBeInRange(options.TitleColor,0,1)} = [0 0 0];
    options.SubtitleColor (1,3) double {mustBeInRange(options.SubtitleColor,0,1)} = [0.65 0.65 0.65];
end

if strlength(Title) > 0
    title(ax, Title, ...
        'FontName',options.FontName, ...
        'Color',options.TitleColor, ...
        'FontWeight','bold', ...
        'FontSize', options.FontSize+2);
end

if strlength(Subtitle) > 0
    subtitle(ax, Subtitle, ...
        'FontName', options.FontName, ...
        'Color', options.SubtitleColor, ...
        'FontSize', options.FontSize);
end

end
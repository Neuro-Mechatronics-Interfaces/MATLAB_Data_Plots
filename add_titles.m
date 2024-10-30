function add_titles(ax, Title, Subtitle, options)
%ADD_TITLES Adds a title and subtitle to an axis with customizable font and color options.
%
% Syntax:
%   add_titles(ax, Title, Subtitle, options)
%
% Description:
%   This utility function adds a title and subtitle to a specified axis (`ax`) with customizable
%   font properties, including font name, size, and color for both the title and subtitle.
%   The title is displayed in bold, with the subtitle shown in a lighter color below it.
%
% Inputs:
%   ax       - Handle to the axis where the title and subtitle will be added.
%   Title    - Text for the title (must be a character vector or string scalar).
%   Subtitle - Text for the subtitle (must be a character vector or string scalar).
%
% Options:
%   FontName      - Font name for the title and subtitle text (default: "Tahoma").
%   FontSize      - Base font size for the title and subtitle text. The title font size will be
%                   slightly larger (default: 14).
%   TitleColor    - RGB color for the title text (default: [0 0 0] for black).
%   SubtitleColor - RGB color for the subtitle text (default: [0.65 0.65 0.65] for gray).
%
% Example:
%   % Set up an axis and add a title and subtitle
%   ax = gca;
%   Title = "Main Plot Title";
%   Subtitle = "This is a descriptive subtitle.";
%   add_titles(ax, Title, Subtitle);
%
%   % Customize font and color for title and subtitle
%   add_titles(ax, Title, Subtitle, 'FontName', "Arial", 'FontSize', 16, ...
%              'TitleColor', [0.2 0.4 0.6], 'SubtitleColor', [0.5 0.5 0.5]);
%
% Notes:
%   - The function automatically adjusts the title font size by increasing it by 2 points relative 
%     to `FontSize` for better distinction between title and subtitle.
%   - The title is shown in bold, while the subtitle is displayed in regular weight.
%
% See also: title, subtitle, axes
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
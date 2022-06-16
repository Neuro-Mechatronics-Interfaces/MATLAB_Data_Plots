function fig = impedances_hdemg(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, varargin)
%IMPEDANCES_HDEMG  Generate impedance plot in grid orientation for HD-EMG
%
% Syntax:
%   fig = plot.impedances_hdemg(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, 'Name', value, ...);
%
% Inputs:
%   SUBJ  - Subject name (e.g. 'Frank' or "Frank")
%   YYYY  - Year (numeric or string, e.g. 2021 or "2021" or '2021')
%   MM    - Month (numeric or string, e.g. 11 or "11" or '11')
%   DD    - Day   (numeric or string, e.g. 18 or "18" or '18')
%   ARRAY - "A" or "B" (or 'A' or 'B') or ["A", "B"] or {'A', 'B'}
%   BLOCK - Recording parameter key (block; numeric or string, e.g. 0 or "0")
%   varargin - (Optional) 'Name', value parameter field/value pairs.
%
%       -> 'DataFolder' - Folder where raw data lives.
%       -> 'FigureFolder' - Folder where figures are saved if no output is
%               requested as an argument.
%       -> 'FileTag' (default is '') - "Tag" to append to "Impedances"
%               figure filename.
%       -> 'CLim' (default is [0 150]) - Sets colormap scaling.
%       -> 'AxLim' (default is [0.5 8.5]) - Sets limits on axes
%       -> 'Colormap' (default is cm.map('redgreen'))
%       -> 'Colorscale' (default is 'log', can be 'linear' or 'log').
%
% Output:
%   fig - Figure handle. If not requested, auto-saves and deletes figure.
%
% See also: Contents

pars = parameters('Impedance');
pars = utils.parse_parameters(pars, varargin{:});

if (numel(ARRAY) > 1) || (numel(BLOCK) > 1)
    fig = gobjects(numel(BLOCK), numel(ARRAY));
    for iB = 1:numel(BLOCK)
        for iA = 1:numel(ARRAY)
            if nargout > 0
                fig(iB, iA) =  plot.impedances_hdemg(SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), pars);
            else
                plot.impedances_hdemg(SUBJ, YYYY, MM, DD, ARRAY(iA), BLOCK(iB), pars);
            end
        end
    end
    if nargout == 0
        clear fig;
    end
    return;
end
[YYYY, MM, DD] = utils.parse_date_args(YYYY, MM, DD);
map = io.load_muscle_map(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, pars.DataFolder);
f = utils.get_block_name(SUBJ, YYYY, MM, DD, ARRAY, BLOCK, ...
        'rootdir_raw', pars.DataFolder, ...
        'rootdir_gen', pars.FigureFolder);
fig = figure('Name', 'Impedance Map', 'Color', 'w', ...
    'Position', [116         207        1022         640]);
ax = axes(fig, 'NextPlot', 'add', ...
    'XColor', [0.75 0.75 0.75], ...
    'YColor', [0.75 0.75 0.75], ...
    'FontName', 'Tahoma', 'FontSize', 16, ...
    'XLim', pars.AxLim, 'YLim', pars.AxLim, 'CLim', pars.CLim,...
    'Colormap', pars.Colormap, 'Colorscale', pars.Colorscale, ...
    'YDir', 'normal', 'XTickLabel', 1:8:57, 'YTickLabel', 1:8);
colormap(ax, 'summer');
cbar = colorbar(ax, 'FontName', 'Tahoma', ...
    'Color', 'k', 'FontSize', 14, 'Box', 'off');
cbar.Label.String = "kΩ";
img = image(ax, 1:8, 1:8, nan(8, 8), 'CDataMapping','scaled');
t_str = sprintf("%s: %s | Array %s", f.Animal, strrep(f.Date, "_", "-"), ARRAY);
title(ax,  t_str, 'FontName', 'Tahoma', 'Color', 'k');
subtitle(ax, 'HD-EMG Electrode Impedance', 'FontName', 'Tahoma', ...
    'Color', [0.65 0.65 0.65]);

channels = flipud(reshape((1:64)', 8, 8));
Z = nan(8, 8);
[xg, yg] = meshgrid(1:8, 8:-1:1);
for ii = 1:numel(Z)
    chname = sprintf('UNI%02d', channels(ii));
    mname = strrep(map.Muscles.(chname), "_", "-");
    Z(ii) = map.Impedances.(chname);
    text(ax, xg(ii), yg(ii), sprintf('%s: %g', mname, Z(ii)), ...
          'FontName', 'Tahoma', ...
          'FontSize', 8, ...
          'FontWeight', 'bold', ...
          'Color', 'k', ...
          'HorizontalAlignment', 'center', ...
          'VerticalAlignment', 'middle');
end
xlabel(ax, 'Channel Index', 'FontName', 'Tahoma', 'Color', [0.75 0.75 0.75]);
ylabel(ax, 'Channel Index', 'FontName', 'Tahoma', 'Color', [0.75 0.75 0.75]);
img.CData = Z;

if nargout < 1
    
    outfolder = fullfile(f.Generated.Tank, 'figures', 'Impedances');
    if exist(outfolder, 'dir')==0
        mkdir(outfolder);
    end
    default.savefig(fig, fullfile(outfolder, sprintf('Impedances-%s', ARRAY)), pars.FileTag); 
end

end
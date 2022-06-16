function fig = synergy(s, tag)
%SYNERGY Plot muscle synergies from nnmf results
%
% Syntax:
%   fig = plot.synergy(s, tag);
%
% Inputs:
%   s   - Struct with fields:
%           'W'  -  left-hand matrix from NNMF
%           'H'  -  right-hand matrix from NNMF
%           'ts' -  timesteps corresponding to rows of 'W'
%
%   tag - (Optional; default is "EMG Synergies") -- Tag for saving figure
%
% Output:
%   fig - Figure handle. If not requested, then default behavior is to use
%           "tag" to generate a file name and save/delete the figure.
%           -> This latter behavior only works if you have a 'parameters'
%               function in your workspace folder that handles the
%               parameter 'preliminary_output_folder'.
%
% See also: Contents

if numel(s) > 1
    fig = gobjects(size(s));
    if size(tag, 1) == 1
        tag = repmat(string(tag), numel(s));
    end
    for ii = 1:numel(s)
        if nargout > 1
            if nargin > 1
                fig(ii) = plot.synergy(s(ii), tag(ii));
            else
                fig(ii) = plot.synergy(s(ii));
            end
        else
            if nargin > 1
                plot.synergy(s(ii), tag(ii));
            else
                plot.synergy(s(ii));
            end
        end        
    end
    return;
end

fig = figure('Name', 'NNMF EMG Synergies', 'Color', 'w', ...
    'Units','Normalized', 'Position', [0.2 0.2 0.5 0.5]);
k = sum(sum(s.H,2) > 0);
L = tiledlayout(fig, k, 3);
for iSyn = 1:k
    ax = nexttile(L, (iSyn-1)*3 + 2, [1 2]);
    set(ax, ... 'YLim', [0 250], ...
        'NextPlot', 'add', 'FontName', 'Arial', ...
        'XColor', 'k', 'YColor', 'k', 'LineWidth', 1.5); 
    plot(ax, s.ts * 1e3, s.W(:, iSyn), ...
        'Color', ax.ColorOrder(iSyn,:), 'LineWidth', 2); 
    xlabel(ax, 'Time (ms)'); 
    ylabel(ax, '(\muV)'); 
    title(ax, sprintf('Synergy-%g', iSyn)); 
    
    ax = nexttile(L, (iSyn-1)*3 + 1, [1 1]);
    format_axes_EMG(ax);
    imagesc(ax, flipud(reshape(s.H(iSyn, :), 8, 8))); 
    colorbar(ax);
end

if nargout < 1
    if nargin < 2
        tag = "EMG_Synergies";
    end
    out_folder = fullfile(parameters('preliminary_output_folder'), 'Synergies');
    if exist(out_folder, 'dir') == 0
        mkdir(out_folder);
    end
    default.savefig(fig, out_folder, tag);
end
end
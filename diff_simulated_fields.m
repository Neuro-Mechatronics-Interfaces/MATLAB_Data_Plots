function fig = diff_simulated_fields(SUBJ, YYYY, MM, DD, BLOCK_A, BLOCK_B, varargin)
%DIFF_SIMULATED_FIELDS  Plots difference between simulated fields from Block_A and Block_B
%
% fig = plot.diff_simulated_fields(SUBJ, YYYY, MM, DD, BLOCK_A, BLOCK_B);
%
% Inputs:
%   SUBJ  - Subject name (e.g. 'Frank' or "Frank")
%   YYYY  - Year (numeric or string, e.g. 2021 or "2021" or '2021')
%   MM    - Month (numeric or string, e.g. 11 or "11" or '11')
%   DD    - Day   (numeric or string, e.g. 18 or "18" or '18')
%   BLOCK_A - Recording parameter key for first block. "Red" diff values
%               (positive) correspond to spatial regions in which this
%               block field is strong compared to BLOCK_B.
%   BLOCK_B - Recording parameter key for the second block. "Blue" diff
%               values (negative) correspond to spatial regions in which
%               this block field is strong compared to BLOCK_A.
%   varargin - (Optional) 'Name', value parameter field/value pairs. 
%
% Output:
%   fig   - Figure handle. If no output is requested, figures are
%           automatically saved and deleted using `default.savefig`. The
%           output saved results are in `generated_data` and follow the
%           naming convention provided by the input arguments.
%
% See also: Contents, load_block_simulation

pars = struct;
pars.CLim_Delta = [];
pars.CMap_Delta = 'bone';
pars.CLim_Rho = [0 9];
pars.CMap_Rho = 'winter';
[pars.Output_Root, pars.Version] = parameters('generated_data_folder', 'version');

% Handle parsing of `pars`
pars = utils.parse_parameters(pars, varargin{:});

if (numel(BLOCK_A) > 1) || (numel(BLOCK_B) > 1)
    if numel(BLOCK_A) ~= numel(BLOCK_B)
        error('If BLOCK_A or BLOCK_B is an array, then both BLOCK_A and BLOCK_B must have an equal number of elements.');
    end
    if nargout < 1
        for iB = 1:numel(BLOCK_A)
            plot.diff_simulated_fields(SUBJ, YYYY, MM, DD, BLOCK_A(iB), BLOCK_B(iB), pars);
        end
    else
        fig = gobjects(size(BLOCK_A));
        for iB = 1:numel(BLOCK_A)
            fig(iB) = plot.diff_simulated_fields(SUBJ, YYYY, MM, DD, BLOCK_A(iB), BLOCK_B(iB), pars);
        end
    end
    return;
end
    
[YYYY, MM, DD] = utils.parse_date_args(YYYY, MM, DD);

% Get the two simulated datasets
F_a = io.load_block_simulation(SUBJ, YYYY, MM, DD, BLOCK_A);
F_b = io.load_block_simulation(SUBJ, YYYY, MM, DD, BLOCK_B);

% Compute the deltas over full mesh
delta_Z = F_a.CData - F_b.CData;
[deltaGlobalMax, ~] = max(abs(delta_Z(:)));
[deltaMax, iDeltaMax] = max(delta_Z(:));
[deltaMin, iDeltaMin] = min(delta_Z(:));

% Compute products over full mesh using normalized values
z_a = (F_a.CData-mean(F_a.CData(:))) ./ std(F_a.CData(:));
z_b = (F_b.CData-mean(F_b.CData(:)))./std(F_b.CData(:));
rho_Z = z_a .* z_b;
[rhoMax, iRhoMax] = max(rho_Z(:));

% Convert mesh to milliseconds
X = F_a.x*1e3;
Y = F_a.y*1e3; 

% Get pattern names
[~, fa, ~] = fileparts(F_a.pattern_file);
fa = strrep(fa, '_', '\_');
[~, fb, ~] = fileparts(F_b.pattern_file);
fb = strrep(fb, '_', '\_');

% Generate figure
fig = figure('Name', 'Simulated Field Differences',...
         'Units','Normalized', ...
         'Position',[0.1 0.1 0.8 0.8],...
         'Color', 'w', ...
         'UserData', struct('Version', pars.Version, 'pars', pars));
L = tiledlayout(fig, 1, 2);
% First, show the delta grid
ax = nexttile(L);
set(ax, 'XLim', [min(X) max(X)], 'YLim', [min(Y), max(Y)], ...
    'FontName', 'Tahoma', 'FontSize', 14, 'XColor', 'k', 'YColor', 'k', ...
    'LineWidth', 1.5, 'NextPlot', 'add', 'XDir', 'reverse', ...
    'CLim', [-deltaGlobalMax, deltaGlobalMax]);
colormap(ax, pars.CMap_Delta);
imagesc(ax, X, Y, delta_Z);
Y = F_a.X*1e3;
X = F_a.Y*1e3;
scatter(ax, X(iDeltaMax), Y(iDeltaMax), 16, 'b', 'Marker', 'o', 'LineWidth', 2.5);
txt = sprintf('   \\leftarrow <%3.1f, %3.1f>: %3.1fmA/cm^2', round(X(iDeltaMax), 1), round(Y(iDeltaMax), 1), round(deltaMax, 1));
text(ax, X(iDeltaMax), Y(iDeltaMax), txt, 'Color', 'b', ...
    'FontWeight', 'bold', 'FontSize', 14, ...
    'HorizontalAlignment', 'left', 'FontName', 'Tahoma');
scatter(ax, X(iDeltaMin), Y(iDeltaMin), 16, 'r', 'Marker', 'o', 'LineWidth', 2.5);
txt = sprintf('   \\leftarrow <%3.1f, %3.1f>: %3.1fmA/cm^2', round(X(iDeltaMin), 1), round(Y(iDeltaMin), 1), round(deltaMin, 1));
text(ax, X(iDeltaMin), Y(iDeltaMin), txt, 'Color', 'r', ...
    'FontWeight', 'bold', 'FontSize', 14, ...
    'HorizontalAlignment', 'left', 'FontName', 'Tahoma');

xlabel(ax, 'ML (mm)', 'FontName', 'Tahoma', 'FontSize', 18, 'FontWeight', 'bold');
ylabel(ax, 'AP (mm)', 'FontName', 'Tahoma', 'FontSize', 18, 'FontWeight', 'bold');
title(ax, '\Delta(A, B)', ...
    'FontName', 'Tahoma', 'FontSize', 20, 'FontWeight', 'bold');
c = colorbar(ax);
set(c.Label, ...
    'FontName', 'Tahoma', ...
    'FontSize', 14, ...
    'Color', 'k', ...
    'String', 'mA/cm^2');
if ~isempty(pars.CLim_Delta)
    set(ax, 'CLim', pars.CLim_Delta); 
end

% Next, show the cross-correlation grid
ax = nexttile(L);
% Convert mesh to milliseconds
X = F_b.x*1e3;
Y = F_b.y*1e3; 
set(ax, 'XLim', [min(X) max(X)], 'YLim', [min(Y), max(Y)], ...
    'FontName', 'Tahoma', 'FontSize', 14, 'XColor', 'k', 'YColor', 'k', ...
    'LineWidth', 1.5, 'NextPlot', 'add', 'XDir', 'reverse');
colormap(ax, pars.CMap_Rho);
imagesc(ax, X, Y, rho_Z);
Y = F_b.X*1e3;
X = F_b.Y*1e3;
scatter(ax, X(iRhoMax), Y(iRhoMax), 16, 'w', 'Marker', 'o', 'LineWidth', 1.5);
txt = sprintf('   \\leftarrow <%3.1f, %3.1f>: %3.1f', round(X(iRhoMax), 1), round(Y(iRhoMax), 1), round(rhoMax, 1));
text(ax, X(iRhoMax), Y(iRhoMax), txt, 'Color', 'w', ...
    'FontWeight', 'bold', 'FontSize', 14, ...
    'HorizontalAlignment', 'left', 'FontName', 'Tahoma');

xlabel(ax, 'ML (mm)', ...
    'FontName', 'Tahoma', ...
    'FontSize', 18,       ...
    'FontWeight', 'bold');
title(ax, 'z(A) .* z(B)', ...
    'FontName', 'Tahoma', ...
    'FontSize', 20, ...
    'FontWeight', 'bold');
title(L, ['||J|| Comparison' newline sprintf('(A) %s vs (B) %s', fa, fb)], ...
    'FontName', 'Tahoma', 'FontSize', 24, 'FontWeight', 'bold');
c = colorbar(ax);
set(c.Label, ...
    'FontName', 'Tahoma', ...
    'FontSize', 14, ...
    'Color', 'k', ...
    'String', '(z-score^2)');
if ~isempty(pars.CLim_Rho)
    set(ax, 'CLim', pars.CLim_Rho); 
end

if nargout < 1
    tank = sprintf('%s_%04d_%02d_%02d', SUBJ, YYYY, MM, DD); % data "tank"
    out_folder = fullfile(pars.Output_Root, SUBJ, tank, 'figures', 'Stim_Differences');
    if exist(out_folder, 'dir') == 0
        try %#ok<TRYNC>
            mkdir(out_folder);
        end
    end
    out_name = fullfile(out_folder, sprintf('%s_Stim_Difference', tank));
    default.savefig(fig, out_name, sprintf("%s-%s", num2str(BLOCK_A), num2str(BLOCK_B)), false);
end

end
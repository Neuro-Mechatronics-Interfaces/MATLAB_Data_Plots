function fig = tiled_correlogram(bin_centers, count_data, options)
%TILED_CORRELOGRAM Plot single axes, tiled, or batch of tiled correlograms.
%
% Inputs:
%   bin_centers (1,:) double
%   count_data - Can be vector same size as bin_centers; 
%                   cell array of such vectors (makes tiled layout), or 
%                   cell array of such cell arrays (returns multiple 
%                                                   figures in array).
arguments
    bin_centers (1,:) double
    count_data
    options.Title {mustBeTextScalar} = "";
    options.Subtitle {mustBeTextScalar} = "";
    options.BarOptions (1,:) cell = {};
    options.XLabel {mustBeTextScalar} = "Lag (ms)";
    options.YLabel {mustBeTextScalar} = "# per Spike";
    options.LabelTrigger string = "Spike";
    options.LabelTarget string = "MUAP";
    options.MinNonzeroBins (1,1) double = 10;
end

if iscell(count_data)
    if iscell(count_data{1})
        fig = gobjects(numel(count_data),1);
        for ii = 1:numel(count_data)
            fig(ii) = figure('Name', 'Correlogram', 'Color', 'w');
            L = tiledlayout(fig(ii),'flow');
            for ik = 1:numel(count_data{ii})
                if nnz(count_data{ii}{ik}) < options.MinNonzeroBins
                    continue;
                end
                ax = nexttile(L);
                set(ax,'NextPlot','add','FontName','Tahoma','FontSize',14);
                bar(ax, bin_centers, count_data{ii}{ik}, ...
                    'FaceColor', 'k', 'EdgeColor', 'none', ...
                    options.BarOptions{:});
                if strlength(options.Title) < 1
                    title(ax, sprintf("%s-%02d", options.LabelTarget, ik), ...
                        'FontName','Tahoma','Color','k','FontSize',6);
                else
                    title(ax, sprintf("%s-%02d -> %s-%02d", options.LabelTrigger, ii, options.LabelTarget, ik), ...
                        'FontName','Tahoma','Color','k','FontSize',6);
                end
                
            end
            xlabel(L, options.XLabel, 'FontName','Tahoma','Color','k');
            ylabel(L, options.YLabel, 'FontName','Tahoma','Color','k');
            if strlength(options.Title) < 1
                titleStr = sprintf("Trigger: %s-%02d", options.LabelTrigger, ii);
            else
                titleStr = options.Title;
            end
            plot.add_titles(L, titleStr, options.Subtitle);
        end
    else
        fig = figure('Name', 'Correlogram', 'Color', 'w');
        L = tiledlayout(fig,'flow');
        for ik = 1:numel(count_data)
            if nnz(count_data{ik})<options.MinNonzeroBins
                continue;
            end
            ax = nexttile(L);
            set(ax,'NextPlot','add','FontName','Tahoma','FontSize',14);
            bar(ax, bin_centers, count_data{ik}, 'FaceColor', 'k', 'EdgeColor', 'none', ...
                options.BarOptions{:});
            title(ax, sprintf("%s-%02d", options.LabelTarget, ik), ...
                'FontName','Tahoma','Color','k','FontSize',6);
        end
        xlabel(L, options.XLabel, 'FontName','Tahoma','Color','k');
        ylabel(L, options.YLabel, 'FontName','Tahoma','Color','k');
        plot.add_titles(L, options.Title, options.Subtitle);
    end
else
    fig = figure('Name', 'Correlogram', 'Color', 'w');
    ax = axes(fig,'NextPlot','add','FontName','Tahoma','FontSize',14);
    bar(ax, bin_centers, count_data, 'FaceColor', 'k', 'EdgeColor', 'none', ...
                 options.BarOptions{:});
    xlabel(ax, options.XLabel, 'FontName','Tahoma','Color','k');
    ylabel(ax, options.YLabel, 'FontName','Tahoma','Color','k');
    plot.add_titles(ax, options.Title, options.Subtitle);
end
end
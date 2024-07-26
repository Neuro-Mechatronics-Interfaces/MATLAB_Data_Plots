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
    options.CDataFaces (:,3) = [nan nan nan];
    options.CDataEdges (:,3) = [nan nan nan];
    options.XUnits {mustBeTextScalar} = "ms";
    options.YLabel {mustBeTextScalar} = "# per Spike";
    options.LabelTrigger string = "Spike";
    options.LabelTarget string = "MUAP";
    options.MinNonzeroBins (1,1) double = 10;
    options.AddPeakLatency (1,1) logical = true;
    options.Threshold (1,1) double = 2;
end


if strlength(options.Subtitle) < 1
    subtitle_str = sprintf('Mean + %g SD', options.Threshold);
else
    subtitle_str = options.Subtitle;
end

if iscell(count_data)
    if iscell(count_data{1})
        if isnan(options.CDataFaces(1,1))
            cdata_faces = zeros(numel(count_data{1}),3);
        else
            cdata_faces = options.CDataFaces;
        end
        if isnan(options.CDataEdges(1,1)) 
            cdata_edges = zeros(numel(count_data),3);
        else
            cdata_edges = options.CDataEdges;
        end
        fig = gobjects(numel(count_data),1);
        for ii = 1:numel(count_data)
            fig(ii) = figure('Name', 'Correlogram', 'Color', 'w', ...
                'Units', 'inches', 'Position', [2 2 8 5], ...
                'UserData', struct('name',string(sprintf("%s-%02d",options.LabelTrigger,ii))));
            L = tiledlayout(fig(ii),'flow');
            for ik = 1:numel(count_data{ii})
                if nnz(count_data{ii}{ik}) < options.MinNonzeroBins
                    continue;
                end
                ax = nexttile(L);
                set(ax,'NextPlot','add','FontName','Tahoma','FontSize',12);
                bar(ax, bin_centers, count_data{ii}{ik}, 'BarWidth', 1, ...
                    'FaceColor', cdata_faces(ik,:), 'EdgeColor', 'none', ...
                    options.BarOptions{:});
                if options.AddPeakLatency
                    mu = mean(count_data{ii}{ik});
                    sd = std(count_data{ii}{ik});
                    thresh = mu+options.Threshold*sd;
                    yline(ax, thresh, 'LineStyle', ':', ...
                        'LineWidth', 1.25, ...
                        'Color', 'k');
                    [pk,loc] = max(count_data{ii}{ik});
                    if pk > thresh
                        line(ax,bin_centers(loc),pk,'Color','b','LineStyle','none','Marker','*','MarkerIndices',1);
                        % text(ax,bin_centers(loc),pk+0.1*thresh,...
                        %     sprintf("%.2f%ss/%s",round(pk,2),options.LabelTarget,options.LabelTrigger), ...
                        %     'FontSize', 6, 'Color', 'b', 'BackgroundColor', 'w');
                        text(ax,bin_centers(loc),pk-0.05*thresh,...
                            sprintf("%.1f %s", bin_centers(loc), options.XUnits), ...
                            'FontSize', 6, 'Color', 'b', 'BackgroundColor', 'w');
                    end
                end
                if strlength(options.Title) < 1
                    title(ax, sprintf("%s-%02d", options.LabelTarget, ik), ...
                        'FontName','Tahoma','Color',cdata_faces(ik,:),'FontSize',6);
                else
                    title(ax, sprintf("%s-%02d -> %s-%02d", options.LabelTrigger, ii, options.LabelTarget, ik), ...
                        'FontName','Tahoma','Color','k','FontSize',6);
                end
                
            end
            xlabel(L, sprintf("Lag (%s)",options.XUnits), 'FontName','Tahoma','Color','k');
            ylabel(L, options.YLabel, 'FontName','Tahoma','Color','k');
            if strlength(options.Title) < 1
                titleStr = sprintf("Trigger: %s-%02d", options.LabelTrigger, ii);
                cdata_Title = cdata_edges(ii,:);
            else
                titleStr = options.Title;
                cdata_Title = [0 0 0];
            end
            plot.add_titles(L, titleStr, subtitle_str, "TitleColor", cdata_Title);
        end
    else
        if isnan(options.CDataFaces(1,1))
            cdata_faces = zeros(numel(count_data),3);
        else
            cdata_faces = options.CDataFaces;
        end
        if isnan(options.CDataEdges(1,1)) 
            cdata_edges = [0 0 0];
        else
            cdata_edges = options.CDataEdges;
        end
        fig = figure('Name', 'Correlogram', 'Color', 'w', ...
                'Units', 'inches', 'Position', [2 2 8 5], ...
                'UserData', struct('name',string(sprintf("%s-Correlograms",options.LabelTrigger))));
        L = tiledlayout(fig,'flow');
        for ik = 1:numel(count_data)
            if nnz(count_data{ik})<options.MinNonzeroBins
                continue;
            end
            ax = nexttile(L);
            set(ax,'NextPlot','add','FontName','Tahoma','FontSize',12);
            bar(ax, bin_centers, count_data{ik}, ...
                'FaceColor', cdata_faces(ik,:), ...
                'EdgeColor', cdata_edges, ...
                options.BarOptions{:});
            title(ax, sprintf("%s-%02d", options.LabelTarget, ik), ...
                'FontName','Tahoma','Color',cdata_faces(ik,:),'FontSize',6);
            if options.AddPeakLatency
                mu = mean(count_data{ik});
                sd = std(count_data{ik});
                thresh = mu+options.Threshold*sd;
                yline(ax, thresh, 'LineStyle', ':', ...
                    'LineWidth', 1.25, ...
                    'Color', 'k');
                [pk,loc] = max(count_data{ik});
                if pk > thresh
                        line(ax,bin_centers(loc),pk,'Color','b','LineStyle','none','Marker','*','MarkerIndices',1);
                    % text(ax,bin_centers(loc),pk+0.1*thresh,...
                    %     sprintf("%.2f%ss/%s",round(pk,2),options.LabelTarget,options.LabelTrigger), ...
                    %     'FontSize', 6, 'Color', 'b', 'BackgroundColor', 'w');
                    text(ax,bin_centers(loc),pk-0.05*thresh,...
                        sprintf("%.1f %s", bin_centers(loc), options.XUnits), ...
                        'FontSize', 6, 'Color', 'b', 'BackgroundColor', 'w');
                end
            end
        end
        xlabel(L, sprintf("Lag (%s)",options.XUnits), 'FontName','Tahoma','Color','k');
        ylabel(L, options.YLabel, 'FontName','Tahoma','Color','k');
        plot.add_titles(L, options.Title, subtitle_str);
    end
else
    if isnan(options.CDataFaces(1,1))
        cdata_faces = [0 0 0];
    else
        cdata_faces = options.CDataFaces;
    end
    if isnan(options.CDataEdges(1,1)) 
        cdata_edges = [0 0 0];
    else
        cdata_edges = options.CDataEdges;
    end
    fig = figure('Name', 'Correlogram', 'Color', 'w', ...
                'Units', 'inches', 'Position', [2 2 8 5], ...
                'UserData', struct('name',string(sprintf("%s-Correlograms",options.LabelTrigger))));
    ax = axes(fig,'NextPlot','add','FontName','Tahoma','FontSize',14);
    bar(ax, bin_centers, count_data, ...
        'FaceColor', cdata_faces, 'EdgeColor',cdata_edges, ...
                 options.BarOptions{:});
    xlabel(ax, sprintf("Lag (%s)",options.XUnits), 'FontName','Tahoma','Color','k');
    ylabel(ax, options.YLabel, 'FontName','Tahoma','Color','k');
    if options.AddPeakLatency
        mu = mean(count_data);
        sd = std(count_data);
        thresh = mu+options.Threshold*sd;
        yline(ax, thresh, 'LineStyle', ':', ...
            'LineWidth', 1.25, ...
            'Color', 'k');
        [pk,loc] = max(count_data);
        if pk > thresh
            line(ax,bin_centers(loc),pk,'Color','b','LineStyle','none','Marker','*','MarkerIndices',1);
            % text(ax,bin_centers(loc),pk+0.1*thresh,...
            %     sprintf("%.2f%ss/%s",round(pk,2),options.LabelTarget,options.LabelTrigger), ...
            %     'FontSize', 6, 'Color', 'b', 'BackgroundColor', 'w');
            text(ax,bin_centers(loc),pk-0.05*thresh,...
                sprintf("%.1f %s", bin_centers(loc), options.XUnits), ...
                'FontSize', 6, 'Color', 'b', 'BackgroundColor', 'w');
        end
    end
    plot.add_titles(ax, options.Title, subtitle_str);
end
end
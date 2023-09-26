function [fig, handles] = direct_form_I(b, a, options)
%DIRECT_FORM_I Plot Direct Form-I implementation block diagram and pole/zero plot.
%
% Syntax:
%   [fig, handles] = plot.direct_form_I(b,a,'Name',value,...);
%
% Example:
%   [b,a] = butter(2,0.1,'high');
%   fig = plot.direct_form_I(b,a);
%
% This would plot the block diagram and poles/zeroes for a 2nd order
% butterworth highpass filter with a normalized cutoff frequency of 0.1.
% The Direct-Form-I implementation is less-susceptible to
% numerical-instability issues related to the precision of pole values, but
% requires twice as many taps as the Direct-Form-II realization of the same
% system. 
%
% Inputs:
%   b - System transfer function numerator polynomial coefficients 
%           (i.e. row vector as returned by [b,a] = butter(...))
%   a - System transfer function denominator polynomial coefficients
%           (i.e. row vector as returned by [b,a] = butter(...))
%
% Options:
%   'CoefficientThreshold' (1,1) double = 0.01 -- Minimum polynomial coefficient to plot.
%   'NFrequencies' (1,1) double = 512 -- Number of frequencies in frequency-response plot.
%   'Parent' {matlab.graphics.layout.TiledChartLayout'} -- Can specify the tiled layout container for these axes.
%   'SampleRate' (1,1) double = nan -- Set to numeric value to make frequency axes use Hz on x-axis instead of normalized frequencies.
%   'Threshold' (1,1) double = -6 -- Determine where to draw fiducial threshold on frequency-response axes. Set this to NaN to turn off threshold fiducials.
%   'Transposed' (1,1) logical = false -- Set true to plot transposed Direct-Form-I implementation. 
%   'Type' {mustBeTextScalar, mustBeMember(options.Type, {'Cut', 'Pass'})} = 'Cut' -- 'Cut' for LPF or HPF, 'Pass' for Bandpass
%
% Output:
%   fig - Figure handle containing the block diagram and pole/zero plot.
%   handles - Struct containing other graphics handles like the layout and individual axes.
%
% See also: Contents, plot.direct_form_II

arguments
    b (1,:) double
    a (1,:) double
    options.CoefficientThreshold (1,1) double = 0.001
    options.NFrequencies (1,1) double = 512;
    options.Parent {mustBeA(options.Parent, 'matlab.graphics.layout.TiledChartLayout')} = matlab.graphics.layout.TiledChartLayout.empty();
    options.SampleRate (1,1) double = nan;
    options.Threshold (1,1) double = -6;
    options.Transposed (1,1) logical = false;
    options.Type {mustBeTextScalar, mustBeMember(options.Type, {'Cut', 'Pass'})} = 'Cut';
    options.YLimFreqz (1,2) double = [-30, 10];
end

% % Get the poles, zeros, and filter order % %
b = round(b,3);
a = round(a,3);
[z,p,~] = tf2zp(b,a);
ord = numel(b) - 1;
if strcmpi(options.Type, 'Pass')
    u_samp = -5:150;
    u = zeros(size(u_samp));
    leg_loc = 'southoutside';
else
    u_samp = -5:25;
    u = zeros(size(u_samp));
    leg_loc = 'southoutside';
end
u(u_samp >= 0) = 1;

% %

handles = struct;
if isempty(options.Parent)
    fig = figure('Name', 'Direct Form-I Filter Implementation', ...
        'Color', 'w', 'Position', [600 270 800 720]);
    handles.layout = tiledlayout(fig, 2, 3);
        if options.Transposed
            title(handles.layout,'Transposed Direct Form-I Filter Characteristics', ...
                'FontName','Tahoma','Color',[0.65 0.65 0.65], 'FontWeight','bold');
        else
            title(handles.layout,'Direct Form-I Filter Characteristics', ...
                'FontName','Tahoma','Color',[0.65 0.65 0.65], 'FontWeight','bold');
        end
else
    fig = get(options.Parent, 'Parent');
    while ~isa(fig, 'matlab.ui.Figure')
        fig = get(fig,'Parent');
    end
    handles.layout = options.Parent;
    handles.layout.GridSize = [2 3];
end

% % % Plot the frequency response of this filter % % %
handles.freqz_ax = nexttile(handles.layout, 1, [1 1]);
set(handles.freqz_ax, ...
    'FontName','Tahoma', ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'NextPlot', 'add');
if isnan(options.SampleRate)
    [h,f] = freqz(b,a,options.NFrequencies);
    f = f ./ pi;
    set(handles.freqz_ax, 'XLim', [0, 1]);
    xlabel(handles.freqz_ax, '\omega (Normalized Frequency)', ...
        'FontName','Tahoma',"Color",'k');
    f_unit = '\omega';
else
    [h,f] = freqz(b,a,options.NFrequencies,options.SampleRate);
    set(handles.freqz_ax, 'XLim', [0, options.SampleRate/2]);
    xlabel(handles.freqz_ax, 'f (Hz)', ...
        'FontName','Tahoma',"Color",'k');
    f_unit = 'Hz';
end
P = 20*log10(abs(h));
handles.freqz = plot(handles.freqz_ax, f, P, 'LineWidth', 1.5, 'Color', 'k');
if strcmpi(options.Type, 'Cut')
    if ~isnan(options.Threshold)
        idx = find(P > options.Threshold, 1, 'last');
        handles.fc_vert = xline(handles.freqz_ax, f(idx), 'Color','m','LineWidth',1.2, 'Label', sprintf('%3.1f %s',round(f(idx),1), f_unit),'LabelVerticalAlignment','bottom');
        handles.fc_horz = yline(handles.freqz_ax, P(idx), 'Color','m','LineWidth',1.2, 'Label', sprintf('%3.1f dB', round(P(idx),1)));
    else
        handles.fc_vert = [];
        handles.fc_horz = [];
    end
elseif strcmpi(options.Type, 'Pass')
    [~,idx] = max(P);
    handles.fc_vert = xline(handles.freqz_ax, f(idx), 'Color','m','LineWidth',1.2, 'Label', sprintf('%3.1f %s',round(f(idx),1), f_unit),'LabelVerticalAlignment','top');
    handles.fc_horz = yline(handles.freqz_ax, P(idx), 'Color','m','LineWidth',1.2, 'Label', sprintf('%3.1f dB', round(P(idx),1)));
else
    handles.fc_vert = [];
    handles.fc_horz = [];
end
yl = [max(options.YLimFreqz(1), handles.freqz_ax.YLim(1)), max(handles.freqz_ax.YLim(2), options.YLimFreqz(2))];
set(handles.freqz_ax, 'YLim', yl);
title(handles.freqz_ax, 'Frequency Response', 'FontName','Tahoma','Color','k');
ylabel(handles.freqz_ax, 'Magnitude (dB)', 'FontName','Tahoma','Color','k');

% % % Plot the pole-zero diagram for this filter % % %
handles.pz_ax = nexttile(handles.layout, 2, [1 1]);
set(handles.pz_ax, ...
    'YLim', [-1.75, 1.75], ...
    'YTick', [-1, -0.5, 0, 0.5, 1], ...
    'YAxisLocation', 'right', ...
    'XLim', [-1.75, 1.75], ...
    'XTick', [-1, -0.5, 0, 0.5, 1], ...
    'FontName','Tahoma', ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'NextPlot', 'add');
xlabel(handles.pz_ax, 'Real Part', ...
        'FontName','Tahoma',"Color",'k');
ylabel(handles.pz_ax, 'Imaginary Part', ...
        'FontName','Tahoma',"Color",'k');
hx = xline(handles.pz_ax, 0, 'LineStyle', ':', 'Color', [0.5 0.5 0.5]);
hx.Annotation.LegendInformation.IconDisplayStyle = 'off';
hy = yline(handles.pz_ax, 0, 'LineStyle', ':', 'Color', [0.5 0.5 0.5]);
hy.Annotation.LegendInformation.IconDisplayStyle = 'off';
theta_c = linspace(-pi,pi,61);
plot(handles.pz_ax, cos(theta_c), sin(theta_c), 'Color', 'k', 'LineWidth', 1.5, 'Displayname','Unit Circle');
for ii = 1:numel(z)
    handles.zero(ii) = plot(handles.pz_ax, real(z(ii)), imag(z(ii)),'Marker','o','Color',[0.35 0.35 0.85],'LineWidth',1.25,'DisplayName','Zero', 'LineStyle','none');
    handles.pole(ii) = plot(handles.pz_ax, real(p(ii)), imag(p(ii)),'Marker','x','Color',[0.85 0.35 0.35],'LineWidth',1.25,'DisplayName','Pole', 'LineStyle','none');
    if ii > 1
        handles.zero(ii).Annotation.LegendInformation.IconDisplayStyle = 'off';
        handles.pole(ii).Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end
legend(handles.pz_ax, 'Location','southoutside', ...
    'FontName','Tahoma','TextColor','black', ...
    'Orientation','horizontal');
title(handles.pz_ax, 'Pole-Zero Plot', 'FontName','Tahoma','Color','k');

% % % Plot the step-response for this filter % % %
handles.step_ax = nexttile(handles.layout, 3, [1 1]);
set(handles.step_ax, ...
    'YLim', [-0.25, 1.25], ...
    'YTick',[0, 0.5, 1.0], ...
    'XLim', [u_samp(1), u_samp(end)], ...
    'FontName','Tahoma', ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'NextPlot', 'add');
plot(handles.step_ax, u_samp, u, 'Color', [0.65 0.65 0.65], ...
    'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Ideal Step');
u_act = filter(b,a,u);
plot(handles.step_ax, u_samp, u_act, 'Color', [0.35 0.35 0.35], ...
    'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'Actual Step');
title(handles.step_ax, 'Step Response', 'FontName','Tahoma','Color','k','FontWeight','bold');
xlabel(handles.step_ax, 'n (Sample)', 'FontName','Tahoma','Color','k');
ylabel(handles.step_ax, 'Amplitude (a.u.)', 'FontName','Tahoma','Color','k');
legend(handles.step_ax, 'FontName',"Tahoma", 'TextColor','black','Location',leg_loc);

% % % Plot the block-diagram for this filter % % %
handles.bd_ax = nexttile(handles.layout, 4, [1 3]);
set(handles.bd_ax, ...
    'YLim', [-ord-0.5, 0.5], ...
    'XLim', [-2.5, 2.5], ...
    'FontName','Tahoma', ...
    'XColor', 'none', ...
    'YColor', 'none', ...
    'NextPlot', 'add');
title(handles.bd_ax, 'Block Diagram', 'FontName','Tahoma','Color','k');

if ~options.Transposed
    % Put the blocks on the axes %
    block_scl = 0.3;
    triangle_scl = 0.5;
    circ_scl = 0.4;
    handles.bd_text.a = gobjects(numel(a)-1,1);
    handles.bd_text.b = gobjects(numel(b), 1);
    for ii = 1:numel(b)
        if abs(b(ii)) > options.CoefficientThreshold
            handles.bd_text.b(ii) = add_forward_gain_block(handles.bd_ax, ii-1, b(ii), triangle_scl);
        end
        if ii > 1
            add_delay_block(handles.bd_ax, ii-1, block_scl);
            if abs(a(ii)) > options.CoefficientThreshold
                handles.bd_text.a(ii-1) = add_feedback_gain_block(handles.bd_ax, ii-1, a(ii), triangle_scl);
            end
        end
        if any([abs(b(ii:end)) > options.CoefficientThreshold, abs(a(ii:end)) > options.CoefficientThreshold])
            add_sum_block(handles.bd_ax, ii-1, circ_scl);
        end
    end
    
    % Connect the blocks with arrows %
    % % Input % %
    text(handles.bd_ax,-2.5,0.1,'x[n]','Color','k','FontName','Tahoma','FontAngle','italic','Rotation',45);
    line(handles.bd_ax,[-2.5,-1-triangle_scl*0.5],[0 0],'Color','k','LineStyle','-','LineWidth',2);
    line(handles.bd_ax,[-1+triangle_scl*0.65,-circ_scl*0.65],[0 0],'Color','k','LineStyle','-','Marker','>','MarkerFaceColor','k','MarkerIndices',2,'MarkerSize',10,'LineWidth',2.5,'MarkerEdgeColor','none');
    % % Output % %
    text(handles.bd_ax,2+circ_scl*0.75,0.1,'y[n]','Color','k','FontName','Tahoma','FontAngle','italic','Rotation',45);
    line(handles.bd_ax,[circ_scl*0.5,2.5],[0 0],'Color','k','LineStyle','-','Marker','>','MarkerFaceColor','k','MarkerIndices',2,'MarkerSize',10,'LineWidth',2.5,'MarkerEdgeColor','none');
    for ii = 1:ord
        % Middle column (sums) - vertical %
        if any([abs(b((ii+1):end)) > options.CoefficientThreshold, abs(a((ii+1):end)) > options.CoefficientThreshold])
            line(handles.bd_ax,[0, 0], [-ii+0.5*circ_scl, -(ii-1)-circ_scl*0.75],'Color','k','LineStyle','-','LineWidth',2,'Marker','^','MarkerFaceColor','k','MarkerIndices',2,'MarkerSize',10,'MarkerEdgeColor','none');
        end
        % Forward delay - horizontal %
        if abs(b(ii+1)) > options.CoefficientThreshold
            line(handles.bd_ax,[-1+0.65*triangle_scl, -0.65*circ_scl], [-ii,-ii],'Color','k','LineStyle','-','LineWidth',2,'Marker','>','MarkerFaceColor','k','MarkerIndices',2,'MarkerSize',10,'MarkerEdgeColor','none');
            line(handles.bd_ax,[-2, -1-0.5*triangle_scl],[-ii, -ii],'Color','k','LineStyle','-','LineWidth',2.5);
        end
        % Forward delay - vertical %
        if any(abs(b((ii+1):end)) > options.CoefficientThreshold)
            line(handles.bd_ax,[-2, -2],[-(ii-1), -ii+0.5+0.5*block_scl],'Color','k','LineStyle','-','LineWidth',2.5);
            line(handles.bd_ax,[-2, -2],[-ii+0.5-0.5*block_scl, -ii],'Color','k','LineStyle','-','LineWidth',2.5);
        end
        % Feedback delay - horizontal %
        if abs(a(ii+1)) > options.CoefficientThreshold
            line(handles.bd_ax,[1-0.65*triangle_scl, 0.65*circ_scl], [-ii,-ii],'Color','k','LineStyle','-','LineWidth',2,'Marker','<','MarkerFaceColor','k','MarkerIndices',2,'MarkerSize',10,'MarkerEdgeColor','none');
            line(handles.bd_ax,[2, 1+0.5*triangle_scl],[-ii, -ii],'Color','k','LineStyle','-','LineWidth',2.5);
        end
        % Feedback delay - vertical %
        if any(abs(a((ii+1):end)) > options.CoefficientThreshold)
            line(handles.bd_ax,[2, 2],[-(ii-1), -ii+0.5+0.5*block_scl],'Color','k','LineStyle','-','LineWidth',2.5);
            line(handles.bd_ax,[2, 2],[-ii+0.5-0.5*block_scl, -ii],'Color','k','LineStyle','-','LineWidth',2.5);
        end
    end
end

    function add_delay_block(ax, index, scl)
        yc = -1*index+0.5;
        rectangle(ax, 'Position', [-2-0.5*scl, yc-0.5*scl, scl, scl], ...
            'LineWidth', 1.25, 'LineStyle', '-', 'EdgeColor', 'k');
        text(ax, -2, yc, 'z^{-1}', 'FontName','Tahoma','Color','k',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
        rectangle(ax, 'Position', [2-0.5*scl, yc-0.5*scl, scl, scl], ...
            'LineWidth', 1.25, 'LineStyle', '-', 'EdgeColor', 'k');
        text(ax, 2, yc, 'z^{-1}', 'FontName','Tahoma','Color','k',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end

    function h_a_txt = add_feedback_gain_block(ax, index, value, scl)
        xc = 1;
        yc = -1*index;
        draw_triangle(ax, xc, yc, -1, scl);
        h_a_txt = text(ax, xc+scl/8, yc, sprintf('-a_{%d}=%5.2f',index,-value), ...
            'FontName','Tahoma','FontSize',6,'Color','r',...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle');
    end

    function add_sum_block(ax, index, scl)
        xc = 0;
        yc = -1*index;
        rectangle(ax, ...
            'Position', [xc-0.5*scl, yc-0.5*scl, scl, scl], ...
            'Curvature', [1 1], ...
            'LineWidth', 1.25, ...
            'LineStyle', '-', ...
            'EdgeColor', 'k');
        text(ax, xc, yc, '+', 'FontName','Tahoma','Color','k',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end

    function h_b_txt = add_forward_gain_block(ax, index, value, scl)
        xc = -1;
        yc = -1*index;
        draw_triangle(ax, xc, yc, 1, scl);
        h_b_txt = text(ax, xc-scl/6, yc, sprintf('b_{%d}=%5.3f',index,value), ...
            'FontName','Tahoma','FontSize',6,'Color','b',...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle');
    end

    function draw_triangle(ax, xc, yc, d, scl)
        yt = [yc-0.5*scl, yc+0.5*scl, yc, yc-0.5*scl];
        xt = [xc-d*0.5*scl, xc-d*0.5*scl, xc+d*0.65*scl, xc-d*0.5*scl];
        if d > 0
            c = 'b';
        else
            c = 'r';
        end
        h_tri = line(ax, xt, yt, 'Color', c, 'LineStyle', '-', 'LineWidth', 1.25);
        h_tri.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end

end
function fig = emg_array(samples, options)
%EMG_ARRAY Plot emg array given samples as they come from TMSi SDK.
%
% Syntax:
%   fig = plot.emg_array(samples, 'Name', value, ...);
%
% Inputs:
%   samples - nChannels x nSamples array of sample data
%   options:
%       
arguments
    samples (:, :) double
%     options.ApplyNotch (1,1) double = true; % Set true to apply notch filter (and harmonics)
%     options.NotchFundamental (1,1) double = 60; % Fundamental frequency of notch filter (Hz)
%     options.NotchWidth (1,1) double = 0.5; % Width of notch filter stops (Hz)
    options.BadChannels (1,:) double = []; % Indices for bad channels (e.g. 1-64)
    options.FilterCutoff (1,:) double = 300; % Hz
    options.SampleRate = 4000; % Samples per second
    options.YOffset = 30; % Offset between vertical traces
    options.TLim = [-inf, inf]; % Time limits (seconds)
    options.GridLayout (1,2) double = [8, 8]; % Number of rows, number of columns
    options.RMSScalar (1,1) double = 6.5;
end

i_data = setdiff(2:65,options.BadChannels+1);

% if options.ApplyNotch
%     Y = fft(samples(i_data,:)');
%     freqs = linspace(0,options.SampleRate,size(Y,1));
%     freq_suppress = options.NotchFundamental:options.NotchFundamental:options.SampleRate;
%     i_suppress = false(size(freqs));
%     for ii = 1:numel(freq_suppress)
%         i_suppress = i_suppress | (abs(freqs-freq_suppress(ii)) < options.NotchWidth);
%     end
%     Y(i_suppress,:) = 0;
%     samples(i_data,:) = abs(ifft(Y)');
% end

if numel(options.FilterCutoff) == 1
    [b,a] = butter(4,options.FilterCutoff./(options.SampleRate/2),'high');
else
    [b,a] = butter(4,options.FilterCutoff./(options.SampleRate/2),'bandpass');
end

fdata = filtfilt(b,a,samples(i_data,:)')';
cref = mean(fdata,1);

t_data = (0:(size(fdata,2)-1))./4000;
t_mask = (t_data >= options.TLim(1)) & (t_data < options.TLim(2));
t_data = t_data(1,t_mask);
t_offset = t_data(end)-t_data(1);

y_offset = repmat(options.YOffset.*(0:(options.GridLayout(1)-1))',1,options.GridLayout(1));
y_offset = y_offset(i_data-1)';

x_offset = repmat((0:(options.GridLayout(2)-1)).*t_offset.*1.05,options.GridLayout(2),1);
x_offset = x_offset(i_data-1)';
t_plot = repmat(t_data,numel(i_data),1) + x_offset;

y_data = fdata(:,t_mask)-cref(t_mask);
y_plot = y_data+y_offset;
c_data = turbo(64);
c_data = c_data(i_data-1,:);



fig = figure(...
    'Name','EMG Array Snippets', ...
    'Color','w', ...
    'Position',[299   348   564   554]); 
x_max = max(t_plot(:));
x_min = min(t_plot(:));
x_min = x_min-.05.*(x_max-x_min);
y_max = max(y_plot(:));
y_min = min(y_plot(:));
y_min = y_min -.05.*(y_max-y_min);
ax = axes(fig,'NextPlot','add',...
    'FontName','Tahoma',...
    'ColorOrder',c_data, ...
    'XLim', [x_min, x_max], ...
    'XColor','none',...
    'YLim', [y_min, y_max], ...
    'YColor','none'); 
plot(ax, t_plot',y_plot'); 
t_scale = round(t_offset/2,1);
line(ax, [x_min, x_min + t_scale], [y_min, y_min], ...
    'Color','k','LineWidth',1.25);
text(ax, x_min + t_scale, y_min, sprintf('%3.1f (s)', t_scale), ...
    'FontName','Tahoma', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom');

y_scale = round(rms(y_data(:)).*options.RMSScalar);
line(ax, [x_min, x_min], [y_min, y_min + y_scale], 'Color','k','LineWidth',1.25);
text(ax, x_min, y_min+y_scale, sprintf('%d (\\muV)', y_scale), ...
    'FontName','Tahoma', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom');
end
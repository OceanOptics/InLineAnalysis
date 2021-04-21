function DiagnosticPlot(data, instrument, level)
% Plot all level of processing 3D spectrums of BB and AC sensors to
% quality check along processing
%
%% Author: Guillaume BOurdin
% Date: 26 Nov. 2019
%
% INPUT:
%   data: <NxM table> data containing:
%     - <1xM datenum> time vector
%     - <N-1xM double> data
%   instrument: <char> instrument name
%   level: <1xL cellstr> processing levels to plot
%%
if contains(instrument,'BB')
  wl = data.lambda;
elseif contains(instrument,'AC')
  wla = data.lambda_a;
  wlc = data.lambda_c;
else
  error('Intrument not supported')
end

for j = 1:length(level)
  switch level{j}
    case 'raw'
      if contains(instrument,'AC')
%         day_to_plot = [min(data.data.dt) min(data.data.dt)+0.4];
        day_to_plot = [max([min(data.raw.tsw.dt) min(data.raw.fsw.dt)]) ...
          max([min(data.raw.tsw.dt) min(data.raw.fsw.dt)]) + 0.4];
        tabletoplot = {'tsw', 'fsw'};
        toplot = {'a','c'};
      else
%         day_to_plot = [min(data.data.dt) min(data.data.dt)+8];
        day_to_plot = [max([min(data.raw.tsw.dt) min(data.raw.fsw.dt)]) ...
          max([min(data.raw.tsw.dt) min(data.raw.fsw.dt)]) + 8];
        tabletoplot = {'tsw', 'fsw'};
        toplot = {'beta'};
      end
    case 'bin'
      day_to_plot = [min([data.bin.tsw.dt; data.bin.fsw.dt]) max([data.bin.tsw.dt; data.bin.fsw.dt])];
      tabletoplot = {'tsw', 'fsw'};
      if contains(instrument,'AC')
        toplot = {'a','c'};
      else
        toplot = {'beta'};
      end
    case 'qc'
      day_to_plot = [min([data.qc.tsw.dt; data.qc.fsw.dt]) max([data.qc.tsw.dt; data.qc.fsw.dt])];
      tabletoplot = {'tsw', 'fsw'};
      if contains(instrument,'AC')
        toplot = {'a','c'};
      else
        toplot = {'beta'};
      end
    case 'prod'
      fieldna = fieldnames(data.prod);
      if isempty(fieldna)
        error('No product loaded')
      end
      day_to_plot = [min(data.prod.(fieldna{1}).dt) max(data.prod.(fieldna{1}).dt)];      
      tabletoplot = fieldna(~structfun(@isempty, data.prod))'; % contains(fieldna, {'p','g'}) & 
      if contains(instrument,'AC')
        toplot = {['a' tabletoplot{~contains(tabletoplot, 'QCfailed')}] ...
          ['c' tabletoplot{~contains(tabletoplot, 'QCfailed')}]};
      else
        toplot = {['beta' tabletoplot{1}], ['bb' tabletoplot{1}]};
      end
    otherwise
    error('Level not supported')
  end
  
  for i = 1:size(tabletoplot,2)
    if ~isempty(data.(level{j}).(tabletoplot{i}))
      for k = 1:size(toplot,2)
        sel = data.(level{j}).(tabletoplot{i}).dt >= day_to_plot(1) & ...
          data.(level{j}).(tabletoplot{i}).dt <= day_to_plot(2);
        if contains(instrument,'AC') && contains(toplot{k}, 'a')
          wl = wla;
        elseif contains(instrument,'AC') && contains(toplot{k}, 'c')
          wl = wlc;
        end
        visProd3D(wl, data.(level{j}).(tabletoplot{i}).dt(sel), ...
            data.(level{j}).(tabletoplot{i}).(toplot{k})(sel,:), ...
            false, 'Wavelength', false, j*i+k*10);
        zlabel([(level{j}) ' ' toplot{k} ' ' tabletoplot{i} ' (m^{-1})']);
        xlabel('lambda (nm)');
        ylabel('time');
        % plot plan at 676 nm to check shift in chl a peak wavelength
        if contains(instrument,'AC') && contains(toplot{k}, 'a')
          hold on
          zsc = zlim;
          ysc = ylim;
          surf([676 676; 676 676], ...
            [max(ysc) min(ysc); max(ysc) min(ysc)], ...
            [zsc(2) zsc(2); zsc(1) zsc(1)], 'FaceColor', [0.2 0.8 0.2], ...
            'FaceAlpha',0.3, 'EdgeAlpha',0.2)
          text(676, max(data.(level{j}).(tabletoplot{i}).dt(sel)), zsc(2)-zsc(2)/5, '676 nm')
          hold off
        end
      end
    end
  end
end

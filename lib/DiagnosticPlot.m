function DiagnosticPlot(data, instrument, level)
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
        day_to_plot = [min(data.data.dt) min(data.data.dt)+0.4];
        tabletoplot = {'tsw', 'fsw'};
        toplot = {'a','c'};
      else
        day_to_plot = [min(data.data.dt) min(data.data.dt)+8];
        tabletoplot = {'tsw', 'fsw'};
        toplot = {'beta'};
      end
    case 'bin'
      day_to_plot = [min(data.data.dt) max(data.data.dt)];
      tabletoplot = {'tsw', 'fsw'};
      if contains(instrument,'AC')
        toplot = {'a','c'};
      else
        toplot = {'beta'};
      end
    case 'qc'
      day_to_plot = [min(data.data.dt) max(data.data.dt)];
      tabletoplot = {'tsw', 'fsw'};
      if contains(instrument,'AC')
        toplot = {'a','c'};
      else
        toplot = {'beta'};
      end
    case 'prod'
      day_to_plot = [min(data.data.dt) max(data.data.dt)];
      tabletoplot = {'p'};
      if contains(instrument,'AC')
        toplot = {'ap','cp'};
      else
        toplot = {'betap','bbp'};
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
      end
    end
  end
end

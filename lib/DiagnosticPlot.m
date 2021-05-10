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
  instrument = 'BB';
elseif contains(instrument,'AC')
  wla = data.lambda_a;
  wlc = data.lambda_c;
  instrument = 'AC';
else
  error('Intrument not supported')
end

for j = 1:length(level)
  fieldna = fieldnames(data.(level{j}));
  if strcmp(fieldna, 'bad')
    data.(level{j}).bad = [];
  end
  tabletoplot = fieldna(~structfun(@isempty, data.(level{j})));
  sztoplot = table(0.4, 8, 'VariableNames', {'AC', 'BB'});
  switch level{j}
    case 'raw'
      szdt = NaN(size(tabletoplot, 1), 2);
      for i = 1:size(tabletoplot, 1)
        szdt(i, 1) = min(data.(level{j}).(tabletoplot{i}).dt);
        szdt(i, 2) = min(data.(level{j}).(tabletoplot{i}).dt);
      end
      day_to_plot = [max(szdt(:,1)) max(szdt(:,2)) + sztoplot.(instrument)];
      if contains(instrument,'AC')
        toplot = {'a','c'};
      else
        toplot = {'beta'};
      end
    case {'bin', 'qc'}
      szdt = NaN(size(tabletoplot, 1), 2);
      for i = 1:size(tabletoplot, 1)
        szdt(i, 1) = min(data.(level{j}).(tabletoplot{i}).dt);
        szdt(i, 2) = max(data.(level{j}).(tabletoplot{i}).dt);
      end
      day_to_plot = [min(szdt(:,1)) max(szdt(:,2))];
      if contains(instrument,'AC')
        toplot = {'a','c'};
      else
        toplot = {'beta'};
      end
    case 'prod'
      szdt = NaN(size(tabletoplot, 1), 2);
      for i = 1:size(tabletoplot, 1)
        szdt(i, 1) = min(data.(level{j}).(tabletoplot{i}).dt);
        szdt(i, 2) = max(data.(level{j}).(tabletoplot{i}).dt);
      end
      day_to_plot = [min(szdt(:,1)) max(szdt(:,2))];
      if contains(instrument,'AC')
        toplot = [cellfun(@(x) ['a' x], tabletoplot, 'un', 0) ...
          cellfun(@(x) ['c' x], tabletoplot, 'un', 0)];
        toplot(contains(toplot, 'QCfailed')) = {'ap', 'cp'};
%         toplot = {['a' tabletoplot{~contains(tabletoplot, 'QCfailed')}] ...
%           ['c' tabletoplot{~contains(tabletoplot, 'QCfailed')}]};
      else
        toplot = {['beta' tabletoplot{1}], ['bb' tabletoplot{1}]};
      end
  end
  
  for i = 1:size(tabletoplot,1)
    if ~isempty(data.(level{j}).(tabletoplot{i}))
      for k = 1:size(toplot,2)
        if strcmp(tabletoplot{i}, 'diw') % if DI plot entire dataset
          sel = true(size(data.(level{j}).(tabletoplot{i}).dt));
        else
          sel = data.(level{j}).(tabletoplot{i}).dt >= day_to_plot(1) & ...
            data.(level{j}).(tabletoplot{i}).dt <= day_to_plot(2);
        end
        if contains(instrument,'AC') && contains(toplot{k}, 'a')
          wl = wla;
        elseif contains(instrument,'AC') && contains(toplot{k}, 'c')
          wl = wlc;
        end
        if size(data.(level{j}).(tabletoplot{i}).dt(sel),1) > 1
          visProd3D(wl, data.(level{j}).(tabletoplot{i}).dt(sel), ...
            data.(level{j}).(tabletoplot{i}).(toplot{i, k})(sel,:), ...
            false, 'Wavelength', false, j*i+k*10);
          zlabel([(level{j}) ' ' toplot{i, k} ' (' tabletoplot{i} ') (m^{-1})']);
          xlabel('lambda (nm)');
          ylabel('time');
        else
          visProd2D(wl, data.(level{j}).(tabletoplot{i}).dt(sel), ...
            data.(level{j}).(tabletoplot{i}).(toplot{i, k})(sel,:), ...
            false, j*i+k*10);
          ylabel([(level{j}) ' ' toplot{i, k} ' (' tabletoplot{i} ') (m^{-1})']);
          xlabel('lambda (nm)');
          title(datestr(data.(level{j}).(tabletoplot{i}).dt(sel)));
        end
        % plot plan at 676 nm to check shift in chl a peak wavelength
        if contains(instrument,'AC') && contains(toplot{k}, 'a') && ...
            ~strcmp(tabletoplot{i}, 'diw')
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

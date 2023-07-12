function user_selection = SpectralQC(data, instrument, level, save_figure, prefix, toClean)
% Plot all level of processing 3D spectrums of BB and AC sensors to
% quality check while processing
% Option to save the plots (save_figure)
% Option to hand pick bad spectra (toClean)
%
%% Author: Guillaume Bourdin
% Date: 26 Nov. 2019
%
% INPUT:
%   data: <NxM table> data containing:
%     - <1xM datenum> time vector
%     - <N-1xM double> data
%   instrument: <char> instrument name
%   level: <1xL cellstr> processing levels to plot
%
% OPTIONAL INPUT:
%   save_figure: <1x1 boolean> to save plot or not
%   prefix: <char> prefix to add in the filename of the saved plot
%   toClean: <1x2 cellstr> name of the table and variable to clean hand-picking bad spectra 
%
% OUTPUT:
%   user_selection: <Px1 cellstr> hand-picked time of bad spectra (if variable 'toClean' is input)
%%
if nargin < 3
  error('Not enough input argument')
elseif nargin == 3
  save_figure = false;
  prefix = 'plot';
  toClean = {'',''};
elseif nargin == 4
  prefix = 'plot';
  toClean = {'',''};
elseif nargin == 5
  toClean = {'',''};
elseif nargin > 6
  error('Too many input argument')
end
  
if contains(instrument,'BB')
  wl = data.lambda;
  instrument = 'BB';
  ylab = 'Lambda (nm)';
elseif contains(instrument,'AC')
  wla = data.lambda_a;
  wlc = data.lambda_c;
  instrument = 'AC';
  ylab = 'Lambda (nm)';
elseif contains(instrument,'LISST')
  theta = data.theta;
  instrument = 'LISST';
  ylab = 'theta (°)';
else
  error('Intrument not supported')
end

ACsize_to_plot = 50000;
BBsize_to_plot = 120000;
LISSTsz_to_plot = 120000;

user_selection = [];
for j = 1:length(level)
  if ~isempty(toClean{1})
    tabletoplot = toClean(1);
    % close figure not to clean
    fig_toclose = [(1:4) .* (1:2)' + 1 * 10; (1:4) .* (1:2)' + 2 * 10];
    fig_toclose = fig_toclose(:);
    for c = 1:size(fig_toclose, 1)
      if ishghandle(fig_toclose(c))
        close(fig_toclose(c))
      end
    end
  else
    % get fieldname of data.level structure
    tabletoplot = fieldnames(data.(level{j}));
    % remove fieldname of empty table
    tabletoplot = tabletoplot(~structfun(@isempty, data.(level{j})));
    tabletoplot(strcmp(tabletoplot, 'bad')) = [];
    tabletoplot(strcmp(tabletoplot, 'FiltStat')) = [];
  end
  sztoplot = table(0.4, 8, 8, 'VariableNames', {'AC', 'BB', 'LISST'});
  switch level{j}
    case 'raw'
      szdt = NaN(size(tabletoplot, 1), 2);
      for i = 1:size(tabletoplot, 1)
        szdt(i, 1) = min(data.(level{j}).(tabletoplot{i}).dt);
        szdt(i, 2) = min(data.(level{j}).(tabletoplot{i}).dt);
      end
      day_to_plot = [max(szdt(:,1)) max(szdt(:,2)) + sztoplot.(instrument)];
      if contains(instrument,'AC')
        if ~isempty(toClean{1})
          toplot = toClean(2);
        else
          toplot = repmat({'a','c'}, size(tabletoplot));
        end
        if strcmp(toClean{end}, 'all')
          toClean{end} = 'a';
        end
      elseif contains(instrument,'BB') || contains(instrument,'LISST')
        if ~isempty(toClean{1})
          toplot = toClean(2);
        else
          toplot = repmat({'beta'}, size(tabletoplot));
        end
        if strcmp(toClean{end}, 'all')
          toClean{end} = 'beta';
        end
      end
    case {'bin', 'qc'}
      szdt = NaN(size(tabletoplot, 1), 2);
      for i = 1:size(tabletoplot, 1)
        szdt(i, 1) = min(data.(level{j}).(tabletoplot{i}).dt);
        szdt(i, 2) = max(data.(level{j}).(tabletoplot{i}).dt);
      end
      day_to_plot = [min(szdt(:,1)) max(szdt(:,2))];
      if contains(instrument,'AC')
        if ~isempty(toClean{1})
          toplot = toClean(2);
        else
          toplot = repmat({'a','c'}, size(tabletoplot));
        end
        if strcmp(toClean{end}, 'all')
          toClean{end} = 'a';
        end
      elseif contains(instrument,'BB') || contains(instrument,'LISST')
        if ~isempty(toClean{1})
          toplot = toClean(2);
        else
          toplot = repmat({'beta'}, size(tabletoplot));
        end
        if strcmp(toClean{end}, 'all')
          toClean{end} = 'beta';
        end
      end
    case 'prod'
      szdt = NaN(size(tabletoplot, 1), 2);
      for i = 1:size(tabletoplot, 1)
        szdt(i, 1) = min(data.(level{j}).(tabletoplot{i}).dt);
        szdt(i, 2) = max(data.(level{j}).(tabletoplot{i}).dt);
      end
      day_to_plot = [min(szdt(:,1)) max(szdt(:,2))];
      if contains(instrument,'AC')
        if ~isempty(toClean{1})
          toplot = toClean(2);
        else
          toplot = [cellfun(@(x) ['a' x], tabletoplot, 'un', 0) ...
            cellfun(@(x) ['c' x], tabletoplot, 'un', 0)];
        end
        if strcmp(toClean{end}, 'all')
          toClean{end} = cellfun(@(x) ['a' x], tabletoplot, 'un', 0);
        end
        if any(contains(toplot, 'QCfailed'))
          toplot(contains(toplot, 'QCfailed')) = {'ap', 'cp'};
        end
      elseif contains(instrument,'LISST')
        if ~isempty(toClean{1})
          toplot = toClean(2);
        else
          toplot = [cellfun(@(x) ['beta' x], tabletoplot, 'un', 0) 'PSD' 'VSD'];
        end
        if strcmp(toClean{end}, 'all')
          toClean{end} = cellfun(@(x) ['beta' x], tabletoplot, 'un', 0);
        end
        if any(contains(toplot, 'QCfailed'))
          toplot(contains(toplot, 'QCfailed')) = {'betap', 'bbp'};
        end
      elseif contains(instrument,'BB')
        if ~isempty(toClean{1})
          toplot = toClean(2);
        else
          toplot = [cellfun(@(x) ['beta' x], tabletoplot, 'un', 0) ...
            cellfun(@(x) ['bb' x], tabletoplot, 'un', 0)];
        end
        if strcmp(toClean{end}, 'all')
          toClean{end} = cellfun(@(x) ['bb' x], tabletoplot, 'un', 0);
        end
        if any(contains(toplot, 'QCfailed'))
          toplot(contains(toplot, 'QCfailed')) = {'betap', 'bbp'};
        end
      end
  end
  
  for i = 1:size(tabletoplot,1)
    if ~isempty(data.(level{j}).(tabletoplot{i}))
      for k = 1:size(toplot,2)
        if strcmp(tabletoplot{i}, 'diw') % if DI plot entire dataset
          sel = true(size(data.(level{j}).(tabletoplot{i}).dt));
        else
          sel = false(size(data.(level{j}).(tabletoplot{i}).dt));
          if contains(instrument,'AC') && size(sel,1) < ACsize_to_plot
            sel(1:end) = true;
          elseif contains(instrument,'AC')
            warning('Large dataset, only the first %i %s %s %s spectrum were plotted to save computer memory', ...
              ACsize_to_plot, level{j}, tabletoplot{i}, toplot{i, k})
            sel(1:ACsize_to_plot) = true;
          elseif contains(instrument,'LISST') && size(sel,1) < LISSTsz_to_plot
            sel(1:end) = true;
          elseif contains(instrument,'LISST')
            warning('Large dataset, only the first %i %s %s %s spectrum were plotted to save computer memory', ...
              LISSTsz_to_plot, level{j}, tabletoplot{i}, toplot{i, k})
            sel(1:LISSTsz_to_plot) = true;
          elseif contains(instrument,'BB') && size(sel,1) < BBsize_to_plot
            sel(1:end) = true;
          elseif contains(instrument,'BB')
            warning('Large dataset, only the first %i %s %s %s spectrum were plotted to save computer memory', ...
              BBsize_to_plot, level{j}, tabletoplot{i}, toplot{i, k})
            sel(1:BBsize_to_plot) = true;
          end
%           sel = data.(level{j}).(tabletoplot{i}).dt >= day_to_plot(1) & ...
%             data.(level{j}).(tabletoplot{i}).dt <= day_to_plot(2);
        end
        sel = sel & any(~isnan(data.(level{j}).(tabletoplot{i}).(toplot{i, k})),2);
%         if strcmp(tabletoplot{i}, 'diw') % if DI plot entire dataset
%           sel = true(size(data.(level{j}).(tabletoplot{i}).dt));
%         else
%           sel = data.(level{j}).(tabletoplot{i}).dt >= day_to_plot(1) & ...
%             data.(level{j}).(tabletoplot{i}).dt <= day_to_plot(2);
%         end
        if contains(instrument,'AC') && contains(toplot{i, k}, 'a')
          wl = wla;
        elseif contains(instrument,'AC') && contains(toplot{i, k}, 'c')
          wl = wlc;
        elseif contains(instrument,'LISST')
          wl = theta;
        end
        if sum(sel) > 1
          fh = visProd3D(wl, data.(level{j}).(tabletoplot{i}).dt(sel), ...
            data.(level{j}).(tabletoplot{i}).(toplot{i, k})(sel,:), ...
            false, 'Wavelength', false, j*i+k*10);
          zlabel([(level{j}) ' ' toplot{i, k} ' (' tabletoplot{i} ') (m^{-1})']);
          xlabel(ylab);
          ylabel('time');
          title([instrument ' ' level{j} ' ' tabletoplot{i}], 'FontSize', 16);
        elseif sum(sel) < 2
          fh = visProd2D(wl, data.(level{j}).(tabletoplot{i}).dt(sel), ...
            data.(level{j}).(tabletoplot{i}).(toplot{i, k})(sel,:), ...
            false, j*i+k*10);
          ylabel([(level{j}) ' ' toplot{i, k} ' (' tabletoplot{i} ') (m^{-1})']);
          xlabel(ylab);
          title([instrument ' ' level{j} ' ' tabletoplot{i} ' ' datestr(data.(level{j}).(tabletoplot{i}).dt(sel))], ...
            'FontSize', 16);
        end
        % plot plan at 676 nm to check shift in chl a peak wavelength
        if contains(instrument,'AC') && contains(toplot{i, k}, 'a') && ...
            ~strcmp(tabletoplot{i}, 'diw') && ~isempty(data.(level{j}).(tabletoplot{i}).dt(sel))
          hold on
          ap676 = find(abs(data.lambda_a - 676) == min(abs(data.lambda_a - 676)), 1, 'first');
          zsc = prctile(data.(level{j}).(tabletoplot{i}).(toplot{i, k})(sel,ap676), [0.01 99.99]);
          ysc = ylim;
          surf([676 676; 676 676], [max(ysc) min(ysc); max(ysc) min(ysc)], ...
            [zsc(2)+0.05*zsc(2) zsc(2)+0.05*zsc(2); zsc(1)-abs(0.05*zsc(2)) zsc(1)-abs(0.05*zsc(2))], ...
            'FaceColor', [0.2 0.8 0.2], 'FaceAlpha',0.3, 'EdgeAlpha',0.2)
          text(676, max(data.(level{j}).(tabletoplot{i}).dt(sel)), zsc(2)+0.05*zsc(2)+0.005, '676 nm')
          hold off
        end
        if any(strcmp(toClean{1}, tabletoplot{i})) && any(strcmp(toClean{2}, toplot{i, k}))
          title(['\fontsize{20}\color{red}' instrument ' ' level{j} ' ' tabletoplot{i} ' spectral QC:' newline '\fontsize{16}\color{black}Select single spectrum using datatip and press "d" to delete' newline '(press q to save and quit)'], 'interpreter', 'tex');
          [ ~, ~, ~, user_sel ] = guiSelectOnTimeSeries(fh);
          user_selection = [user_selection; user_sel];
        end
        if save_figure
          pause(0.01)
          if ishghandle(fh) && get(fh,'Number') == j*i+k*10
            if ~isfolder(fullfile(data.path.prod, 'plots'))
              mkdir(fullfile(data.path.prod, 'plots'))
            end
            filename = fullfile(data.path.prod, 'plots', [prefix '_' ...
              datestr(day_to_plot(1), 'yyyymmdd') '_' datestr(day_to_plot(2), 'yyyymmdd') '_' ...
              data.model data.sn '_' level{j} '_' toplot{i, k} '_' tabletoplot{i}]);
            savefig(fh, filename, 'compact')
          end
        end
      end
    end
  end
end

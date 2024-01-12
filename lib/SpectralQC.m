function user_selection = SpectralQC(data, dt_toclean, instrument, level, save_figure, prefix, toClean)
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
%   dt_toclean: <1x2 datenum> days to clean
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
elseif nargin == 4
  save_figure = false;
  prefix = 'plot';
  toClean = {'',''};
elseif nargin == 5
  prefix = 'plot';
  toClean = {'',''};
elseif nargin == 6
  toClean = {'',''};
elseif nargin > 7
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

ACsize_to_plot = 100000; % 50000
BBsize_to_plot = 220000; % 220000
LISSTsz_to_plot = 220000; % 220000

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
  if isempty(tabletoplot)
      warning('No data to plot in %s level %s: ignored', instrument, level{j})
    return
  end
  for i = 1:size(tabletoplot, 1)
    if isempty(data.(level{j}).(tabletoplot{i}))
      warning('Instrument %s level %s table %s empty: ignored', instrument, level{j}, tabletoplot{i})
    end
  end
  switch level{j}
    case 'raw'
      if contains(instrument,'AC')
        if ~isempty(toClean{1})
          toplot = toClean(2);
        else
          toplot = repmat({'a','c'}, size(tabletoplot));
        end
        if strcmp(toClean{end}, 'all')
          toplot{end} = 'a';
          toClean{end} = 'a';
        end
      elseif contains(instrument,'BB') || contains(instrument,'LISST')
        if ~isempty(toClean{1})
          toplot = toClean(2);
        else
          toplot = repmat({'beta'}, size(tabletoplot));
        end
        if strcmp(toClean{end}, 'all')
          toplot{end} = 'beta';
          toClean{end} = 'beta';
        end
      end
    case {'bin', 'qc'}
      if contains(instrument,'AC')
        if ~isempty(toClean{1})
          toplot = toClean(2);
        else
          toplot = repmat({'a','c'}, size(tabletoplot));
        end
        if strcmp(toClean{end}, 'all')
          toClean{end} = 'a';
          toplot{end} = 'a';
        end
      elseif contains(instrument,'BB') || contains(instrument,'LISST')
        if ~isempty(toClean{1})
          toplot = toClean(2);
        else
          toplot = repmat({'beta'}, size(tabletoplot));
        end
        if strcmp(toClean{end}, 'all')
          toplot{end} = 'beta';
          toClean{end} = 'beta';
        end
      end
    case 'prod'
      if contains(instrument,'AC')
        if ~isempty(toClean{1})
          toplot = toClean(2);
        else
          toplot = [cellfun(@(x) ['a' x], tabletoplot, 'un', 0) ...
            cellfun(@(x) ['c' x], tabletoplot, 'un', 0)];
        end
        if strcmp(toClean{end}, 'all')
          toplot(end) = cellfun(@(x) ['a' x], tabletoplot, 'un', 0);
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
          toplot(end) = cellfun(@(x) ['beta' x], tabletoplot, 'un', 0);
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
          toplot(end) = cellfun(@(x) ['bb' x], tabletoplot, 'un', 0);
          toClean{end} = cellfun(@(x) ['bb' x], tabletoplot, 'un', 0);
        end
        if any(contains(toplot, 'QCfailed'))
          toplot(contains(toplot, 'QCfailed')) = {'betap', 'bbp'};
        end
      end
  end
  
  for i = 1:size(tabletoplot,1)
    if ~isempty(data.(level{j}).(tabletoplot{i}))
      % keep only day2clean without margin
      day2clean = data.(level{j}).(tabletoplot{i}).dt >= min(dt_toclean) & data.(level{j}).(tabletoplot{i}).dt < max(dt_toclean)+1;
      sel = false(size(data.(level{j}).(tabletoplot{i}).dt));
      for k = 1:size(toplot,2)
        % select 0.5-99.5th percentiles spectra to plot
        grp = datenum(dateshift(datetime(data.(level{j}).(tabletoplot{i}).dt, 'ConvertFrom', 'datenum'), 'start', 'hour'));
        ugrp = unique(grp);
        outliers = false(size(data.(level{j}).(tabletoplot{i}).dt));
        for g = 1:size(ugrp, 1)
          prct_1_99 = prctile(data.(level{j}).(tabletoplot{i}).(toplot{i, k})(grp == ugrp(g),:), [0.5 99.5], 1);
          outliers(grp == ugrp(g) & any(data.(level{j}).(tabletoplot{i}).(toplot{i, k}) < prct_1_99(1,:) | data.(level{j}).(tabletoplot{i}).(toplot{i, k}) > prct_1_99(2,:), 2)) = true;
        end
        % add total 5-95th percentiles to the spectra to plot
        prct_5_95 = prctile(data.(level{j}).(tabletoplot{i}).(toplot{i, k}), [5 95], 1);
        outliers(any(data.(level{j}).(tabletoplot{i}).(toplot{i, k}) < prct_5_95(1,:) | data.(level{j}).(tabletoplot{i}).(toplot{i, k}) > prct_5_95(2,:), 2)) = true;
        % reduce number of spectra to plot as a function of spectra non-outliers
        foo = false(size(outliers));
        keep_pts = ceil(log10(sum(~outliers))*log10(sum(~outliers))/5);
        foo(1:keep_pts:end) = true;
        non_outliers = ~outliers & foo;
        pts_toclean = find(day2clean & (outliers | non_outliers));
        if contains(instrument,'AC') && size(pts_toclean,1) > ACsize_to_plot
          warning('Large dataset, only the first %i %s %s %s spectrum were plotted to save computer memory', ...
            ACsize_to_plot, level{j}, tabletoplot{i}, toplot{i, k})
          pts_toclean = pts_toclean(1:ACsize_to_plot);
        elseif contains(instrument,'LISST') && size(pts_toclean,1) > LISSTsz_to_plot
          warning('Large dataset, only the first %i %s %s %s spectrum were plotted to save computer memory', ...
            LISSTsz_to_plot, level{j}, tabletoplot{i}, toplot{i, k})
          pts_toclean = pts_toclean(1:LISSTsz_to_plot);
        elseif contains(instrument,'BB') && size(pts_toclean,1) > BBsize_to_plot
          warning('Large dataset, only the first %i %s %s %s spectrum were plotted to save computer memory', ...
            BBsize_to_plot, level{j}, tabletoplot{i}, toplot{i, k})
          pts_toclean = pts_toclean(1:BBsize_to_plot);
        end
        sel(pts_toclean) = true;
        sel = sel & any(~isnan(data.(level{j}).(tabletoplot{i}).(toplot{i, k})),2);
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
        elseif sum(sel) == 1
          fh = visProd2D(wl, data.(level{j}).(tabletoplot{i}).dt(sel), ...
            data.(level{j}).(tabletoplot{i}).(toplot{i, k})(sel,:), ...
            false, j*i+k*10);
          ylabel([(level{j}) ' ' toplot{i, k} ' (' tabletoplot{i} ') (m^{-1})']);
          xlabel(ylab);
          title([instrument ' ' level{j} ' ' tabletoplot{i} ' ' datestr(data.(level{j}).(tabletoplot{i}).dt(sel))], ...
            'FontSize', 16);
        else
          continue
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
              datestr(min(data.(level{j}).(tabletoplot{i}).dt(sel)), 'yyyymmdd') '_' ...
              datestr(max(data.(level{j}).(tabletoplot{i}).dt(sel)), 'yyyymmdd') '_' ...
              data.model data.sn '_' level{j} '_' toplot{i, k} '_' tabletoplot{i}]);
            savefig(fh, filename, 'compact')
          end
        end
      end
    end
  end
end

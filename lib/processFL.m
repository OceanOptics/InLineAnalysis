function [p, g, FiltStat] = processFL(param, tot, filt_qc, filt_raw, filt_bad, di, di_method, filt_method, fth, fth_constants, days2run)
% Note DI is not interpolated as it's assume to be stable in time
% BB3 parameters is a structure
%   param.lambda <1x3 double> wavelength (nm)
%   param.theta <1x1 double> scattering angle (degree)
%   param.dark  <1x3 double> dark
%   param.slope <1x3 double> slope

if exist('fth', 'var')
  % check FTH data
  if ~exist('fth_constants', 'var')
    % Assume most recent FlowControl software
    SWITCH_FILTERED = 1;
    SWITCH_TOTAL = 0;
  else
    SWITCH_FILTERED = fth_constants.SWITCH_FILTERED;
    SWITCH_TOTAL = fth_constants.SWITCH_TOTAL;
  end
  if strcmp(filt_method, '25percentil')
    if isempty(fth.qc.tsw)
      error('No FLOW qc data loaded, when "filt_method == 25percentil", FLOW qc data level should be loaded')
    end
    if isempty(filt_qc)
      error('No BB3 qc data loaded, when "filt_method == 25percentil", BB3 qc data level should be loaded')
    end
    fth_temp = fth.qc.tsw;
  elseif strcmp(filt_method, 'exponential_fit')
    if isempty(fth.raw.tsw)
      error('No FLOW raw data loaded, when "filt_method == exponential_fit", FLOW raw data level should be loaded')
    end
    if isempty(filt_raw)
      error('No BB3 raw data loaded, when "filt_method == exponential_fit", BB3 raw data level should be loaded')
    end
    fth_temp = fth.raw.tsw;
  end
  
  % sort filt_qc data
  filt_qc = sortrows(filt_qc, 'dt');
  fth_temp.swt(fth_temp.swt > 0) = 1;
  % delete duplicats
  [~, L, ~] = unique(fth_temp.dt,'first');
  indexToDump = not(ismember(1:numel(fth_temp.dt),L));
  fth_temp(indexToDump, :) = [];
  
  % interpolate fth_temp.swt onto binned data to fill missing flow data
  fth_interp = table([tot.dt; fth_temp.dt; filt_qc.dt], 'VariableNames', {'dt'});
  % delete duplicats
  [~, L, ~] = unique(fth_interp.dt,'first');
  indexToDump = not(ismember(1:numel(fth_interp.dt),L));
  fth_interp(indexToDump, :) = [];
  % sort dates
  [~,b] = sort(fth_interp.dt);
  fth_interp.dt = fth_interp.dt(b,:);
  fth_interp.swt = interp1(fth_temp.dt, fth_temp.swt, fth_interp.dt, 'previous');
  fth_interp.swt = fth_interp.swt > 0;
  % Find switch events from total to filtered
  sel_start = find(fth_interp.swt(1:end-1) == SWITCH_TOTAL & fth_interp.swt(2:end) == SWITCH_FILTERED);
  % Find switch events from filtered to total
  sel_end = find(fth_interp.swt(1:end-1) == SWITCH_FILTERED & fth_interp.swt(2:end) == SWITCH_TOTAL);
  % Verify selections of filtered period
  if sel_start(1) > sel_end(1); sel_end(1) = []; end
  if sel_start(end) > sel_end(end); sel_end(end+1) = size(fth_interp.swt,1); end
  if size(sel_start,1) ~= size(sel_end,1); error('Inconsistent fth data'); end

  % Compute filtered period median
  filt_avg = table(NaN(size(sel_start)), 'VariableNames', {'dt'});
  filt_avg.fchl = NaN(size(filt_avg,1), 1);
  filt_avg.fchl_avg_sd = NaN(size(filt_avg,1), 1);
  filt_avg.fchl_avg_n = NaN(size(filt_avg,1), 1);
  switch filt_method
    case '25percentil'
      for i=1:size(sel_start, 1)
        sel_filt = fth_interp.dt(sel_start(i)) <= filt_qc.dt & filt_qc.dt <= fth_interp.dt(sel_end(i));
        if sum(sel_filt) > 0
          foo = filt_qc(sel_filt,:);
          if sum(sel_filt) == 1
            filt_avg.dt(i) = foo.dt;
            filt_avg.fchl(i,:) = foo.fchl;
            filt_avg.fchl_avg_sd(i,:) = foo.fchl_avg_sd;
            filt_avg.fchl_avg_n(i,:) = foo.fchl_avg_n;
          else
            perc25 = foo.fchl > prctile(foo.fchl, 25, 1);
            foo.fchl_avg_sd(perc25) = NaN;
            foo.fchl(perc25) = NaN;
            % compute average of all values smaller than 25th percentile for each filter event
            filt_avg.dt(i) = mean(foo.dt(any(~perc25, 2)), 'omitnan');
            filt_avg.fchl(i,:) = mean(foo.fchl, 1, 'omitnan');
            filt_avg.fchl_avg_sd(i,:) = mean(foo.fchl_avg_sd, 1, 'omitnan');
            filt_avg.fchl_avg_n(i) = sum(foo.fchl_avg_n(any(~isnan(foo.fchl), 2)), 'omitnan');
          end
        end
      end
    case 'exponential_fit'
      % Based on method in: Dall’Olmo, G., Westberry, T.K., Behrenfeld, M.J., Boss, 
      %       E., Slade, W.H., 2009. Direct contribution of phytoplankton-sized particles 
      %       to optical backscattering in the open ocean. Biogeosciences Discuss 6, 291–340. 
      %       https://doi.org/10.5194/bgd-6-291-2009
      fprintf('Fitting exponential to filter events ... \n')
      filt_avg.dt = (fth_interp.dt(sel_start) + fth_interp.dt(sel_end)) ./ 2;
      [filt_avg, FiltStat] = FiltExpFit('fchl', filt_avg, filt_raw, filt_bad, fth_interp.dt(sel_start), fth_interp.dt(sel_end));
      fprintf('Done\n')
      % run 25 percentile method on failed exponential fits
      for i=1:size(sel_start, 1)
        sel_filt = fth_interp.dt(sel_start(i)) <= filt_qc.dt & filt_qc.dt <= fth_interp.dt(sel_end(i));
        if sum(sel_filt) > 0
          if any(~FiltStat.exitflag(i,:))
            foo = filt_qc(sel_filt,:);
            if sum(sel_filt) == 1
              filt_avg.dt(i) = foo.dt;
              filt_avg.fchl(i,~FiltStat.exitflag(i,:)) = foo.fchl(:, ~FiltStat.exitflag(i,:));
            else
              perc25 = foo.beta > prctile(foo.beta, 25, 1);
              foo.fchl_avg_sd(perc25) = NaN;
              foo.fchl(perc25) = NaN;
              % compute average of all values smaller than 25th percentile for each filter event
              filt_avg.dt(i) = mean(foo.dt(any(~perc25, 2)), 'omitnan');
              filt_avg.fchl(i, ~FiltStat.exitflag(i,:)) = mean(foo.fchl(:,~FiltStat.exitflag(i,:)), 1, 'omitnan');
              filt_avg.fchl_avg_sd(i, ~FiltStat.exitflag(i,:)) = mean(foo.fchl_avg_sd(:,~FiltStat.exitflag(i,:)), 1, 'omitnan');
              filt_avg.fchl_avg_n(i) = sum(foo.fchl_avg_n(any(~isnan(foo.fchl(:,~FiltStat.exitflag(i,:))), 2)), 'omitnan');
            end
          end
        else
          filt_avg.fchl(i,:) = NaN;
        end
      end
    otherwise
      error('filter event method "filt_method" not supported')
  end
  filt_avg(all(isnan(filt_avg.fchl), 2), :) = [];
else
  filt_avg = filt_qc;
end

% Interpolate filtered on total linearly
filt_interp = table(tot.dt, 'VariableNames', {'dt'});
filt_interp.fchl = interp1(filt_avg.dt, filt_avg.fchl, filt_interp.dt);%, 'linear', 'extrap');
filt_interp.fchl_avg_sd = interp1(filt_avg.dt, filt_avg.fchl_avg_sd, filt_interp.dt);%, 'linear', 'extrap');

% id only day to run in all tables to plot
filt_interp_id = filt_interp.dt >= min(days2run) & filt_interp.dt < max(days2run)+1;
tot_id = tot.dt >= min(days2run) & tot.dt < max(days2run)+1;
filt_avg_id = filt_avg.dt >= min(days2run) & filt_avg.dt < max(days2run)+1;
% plot
if exist('visFlag', 'file') && exist('fth', 'var')
  fh = visFlag([], filt_interp(filt_interp_id, :), tot(tot_id, :), [], filt_avg(filt_avg_id, :), [], 'fchl', round(size(tot.fchl, 2)/2), ...
    [], fth_temp, fth.view.spd_variable);
  title('Check filter event interpolation, press q to continue', 'FontSize', 14)
  legend('Filtered interpolated', 'Total', 'Filtered median', 'Flow rate',...
    'AutoUpdate','off', 'FontSize', 12)
  guiSelectOnTimeSeries(fh);
elseif exist('visFlag', 'file')
  fh = visFlag([], filt_interp(filt_interp_id, :), tot(tot_id, :), [], filt_avg(filt_avg_id, :), [], 'fchl', round(size(tot.fchl, 2)/2), [], []);
  title('Check filter event interpolation, press q to continue', 'FontSize', 14)
  legend('Filtered interpolated', 'Total', 'Filtered median', 'Flow rate',...
    'AutoUpdate','off', 'FontSize', 12)
  guiSelectOnTimeSeries(fh);
end

% Compute fchl particulate
p = table(tot.dt, 'VariableNames', {'dt'});
p.fchlp = tot.fchl - filt_interp.fchl;
p.fchlp_avg_sd = tot.fchl_avg_sd;

% Calibrate beta_p (counts to scientific units)
p.chl = param.slope .* p.fchlp; % Dark independent

% Propagate error
%   Note: Error is not propagated through Scattering & Residual temperature
%         correction as required by SeaBASS
p.chl_sd = param.slope .* sqrt(tot.fchl_avg_sd.^2 + filt_interp.fchl_avg_sd.^2);
p.chl_n = tot.chl_avg_n;

% remove empty lines
p(isnan(p.chl), :)=[];

% remove negative values
p(any(p.chl < 0, 2),:) = [];


%% ag & cg
if ~isempty(di)
  if strcmp(di_method, 'best_di')
    % select DIW with lowest a or c values between 550-650nm
    di_orig = di;
    di_dt = datetime(di_orig.dt, 'ConvertFrom', 'datenum');
    best_di = NaN(size(di_orig,1), 1);
    for i = 1:size(di_orig,1)
      if i == 1 || i == size(di_orig,1)
        iddi = abs(di_dt(i) - di_dt) < hours(72);
      else
        iddi = abs(di_dt(i) - di_dt) < hours(36);
      end
      lowest_di = di_orig.fchl == min(di_orig.fchl(iddi, :), [], 1);
      foo = find(sum(lowest_di, 2) == max(sum(lowest_di, 2)));
      di.fchl(i, :) = di_orig.fchl(foo, :);
      di.fchl_avg_sd(i, :) = di_orig.fchl_avg_sd(foo, :);
      
      best_di(i) = foo(1);
    end
  end
  
  % remove when a and c are full of NaNs
  filt_avg(all(isnan(filt_avg.fchl), 2) & all(isnan(filt_avg.fchl), 2),:) = [];
  
  % Interpolate filtered on Total
  di_interp = table(filt_avg.dt, 'VariableNames', {'dt'});
  di_interp.fchl = interp1(di.dt, di.fchl, di_interp.dt, 'linear', 'extrap');
  di_interp.fchl_avg_sd = interp1(di.dt, di.fchl_avg_sd, di_interp.dt, 'linear', 'extrap');

  % Dissolved = Filtered - DI
  g = table(filt_avg.dt, 'VariableNames', {'dt'});
  g.ag = filt_avg.fchl - di_interp.fchl;
  fprintf('Done\n')

  % Propagate error
  %   Note: Error is not propagated through Scattering & Residual temperature
  %         correction as required by SeaBASS
  g.fchlg_sd = sqrt(filt_avg.fchl_avg_sd.^2 + di_interp.fchl_avg_sd.^2);
  g.fchlg_n = filt_avg.fchl_avg_n;
  fprintf('Done\n')
  
else
  g = table();
end
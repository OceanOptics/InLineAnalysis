function [p, g] = processTAU(tot, filt, di, cdom_base, fth, fth_constants, interpolation_method, di_method)
%% cp
% check FTH data
if ~exist('fth_constants', 'var')
  % Assume most recent FlowControl software
  SWITCH_FILTERED = 1;
  SWITCH_TOTAL = 0;
else
  SWITCH_FILTERED = fth_constants.SWITCH_FILTERED;
  SWITCH_TOTAL = fth_constants.SWITCH_TOTAL;
end
% remove duplicates
[~, L, ~] = unique(fth.qc.tsw.dt,'first');
indexToDump = not(ismember(1:numel(fth.qc.tsw.dt), L));
if sum(indexToDump) > 0
  fprintf('Warning: %i identical dates in FLOW data => deleted\n', sum(indexToDump))
  fth.qc.tsw(indexToDump, :) = [];
end
% remove duplicates
[~, L, ~] = unique(tot.dt,'first');
indexToDump = not(ismember(1:numel(tot.dt), L));
if sum(indexToDump) > 0
  fprintf('Warning: %i identical dates in total data => deleted\n', sum(indexToDump))
  tot(indexToDump, :) = [];
end
% remove duplicates
[~, L, ~] = unique(filt.dt,'first');
indexToDump = not(ismember(1:numel(filt.dt), L));
if sum(indexToDump) > 0
  fprintf('Warning: %i identical dates in filtered data => deleted\n', sum(indexToDump))
  filt(indexToDump, :) = [];
end
% interpolate fth.qc.tsw.swt onto binned data to fill missing flow data
fth_interp = table([tot.dt; fth.qc.tsw.dt; filt.dt], 'VariableNames', {'dt'});
[~,b] = sort(fth_interp.dt); % sort dates
fth_interp.dt = fth_interp.dt(b,:);
fth_interp.swt = interp1(fth.qc.tsw.dt, fth.qc.tsw.swt, fth_interp.dt, 'previous');%, 'linear', 'extrap');
fth_interp.swt = fth_interp.swt > 0;
% Find switch events from total to filtered
sel_start = find(fth_interp.swt(1:end-1) == SWITCH_TOTAL & fth_interp.swt(2:end) == SWITCH_FILTERED);
% Find switch events from filtered to total
sel_end = find(fth_interp.swt(1:end-1) == SWITCH_FILTERED & fth_interp.swt(2:end) == SWITCH_TOTAL);
% Verify selections of filtered period
if sel_start(1) > sel_end(1); sel_end(1) = []; end
if sel_start(end) > sel_end(end); sel_end(end+1) = size(fth_interp.swt,1); end
if size(sel_start,1) ~= size(sel_end,1); error('Inconsistent fth data'); end
% interpolate filtered event

% Compute filtered period median
filt_avg = table((fth_interp.dt(sel_start) + fth_interp.dt(sel_end)) ./ 2, 'VariableNames', {'dt'});
filt_avg.Tau = NaN(size(filt_avg,1), 1);
filt_avg.Beamc = NaN(size(filt_avg,1), 1);
filt_avg.Tau_avg_sd = NaN(size(filt_avg,1), 1);
filt_avg.Beamc_avg_sd = NaN(size(filt_avg,1), 1);
filt_avg.Tau_avg_n = NaN(size(filt_avg,1), 1);
filt_avg.Beamc_avg_n = NaN(size(filt_avg,1), 1);
for i=1:size(sel_start, 1)
  sel_filt = fth_interp.dt(sel_start(i)) <= filt.dt & filt.dt <= fth_interp.dt(sel_end(i));
  foo = filt(sel_filt,:);
  if sum(sel_filt) == 1
    filt_avg.Tau(i,:) = foo.Tau;
    filt_avg.Beamc(i,:) = foo.Beamc;
    filt_avg.Tau_avg_sd(i,:) = foo.Tau_avg_sd;
    filt_avg.Beamc_avg_sd(i,:) = foo.Beamc_avg_sd;
    filt_avg.Tau_avg_n(i) = foo.Tau_avg_n;
    filt_avg.Beamc_avg_n(i) = foo.Beamc_avg_n;
  else
    foo.Tau_avg_sd(foo.Tau < prctile(foo.Tau, 75, 1)) = NaN;
    foo.Beamc_avg_sd(foo.Beamc > prctile(foo.Beamc, 25, 1)) = NaN;
    foo.Tau(foo.Tau < prctile(foo.Tau, 75, 1)) = NaN;
    foo.Beamc(foo.Beamc > prctile(foo.Beamc, 25, 1)) = NaN;
    % compute average of all values smaller than 25th percentile for each filter event
    filt_avg.Tau(i,:) = mean(foo.Tau, 1, 'omitnan');
    filt_avg.Beamc(i,:) = mean(foo.Beamc, 1, 'omitnan');
    filt_avg.Tau_avg_sd(i,:) = mean(foo.Tau_avg_sd, 1, 'omitnan');
    filt_avg.Beamc_avg_sd(i,:) = mean(foo.Beamc_avg_sd, 1, 'omitnan');
    filt_avg.Tau_avg_n(i) = sum(foo.Tau_avg_n(any(~isnan(foo.Tau), 2)), 'omitnan');
    filt_avg.Beamc_avg_n(i) = sum(foo.Beamc_avg_n(any(~isnan(foo.Beamc), 2)), 'omitnan');
  end
end
filt_avg(all(isnan(filt_avg.Tau), 2) | all(isnan(filt_avg.Beamc), 2), :) = [];

% check if cdom data loaded
if strcmp(interpolation_method, 'CDOM')
  if ~isempty(cdom_base)
    if ~any(cdom_base.dt >= min([tot.dt; filt.dt]) & cdom_base.dt <= max([tot.dt; filt.dt]))
      fprintf('Warning: CDOM dates do not correspond to ACS dates: interpolation switched to "linear"\n')
      interpolation_method = 'linear';
    end
  else
    fprintf('Warning: CDOM data not loaded: interpolation switched to "linear"\n')
    interpolation_method = 'linear';
  end
end

switch interpolation_method
  case 'CDOM'
    % Require both CDOM & Switch position
    % remove duplicates
    [~, L, ~] = unique(cdom_base.dt,'first');
    indexToDump = not(ismember(1:numel(cdom_base.dt), L));
    if sum(indexToDump) > 0
      fprintf('Warning: %i identical dates in CDOM data => deleted\n', sum(indexToDump))
      cdom_base(indexToDump, :) = [];
    end
    cdom = table();
    foo_dt = datetime(cdom_base.dt, 'ConvertFrom', 'datenum');
    cdom.dt = (datenum(foo_dt(1):median(diff(foo_dt)):foo_dt(end)))';
    cdom.fdom = interp1(cdom_base.dt, cdom_base.fdom, cdom.dt, 'linear');
    % keep only leg data
    cdom_base = cdom_base(cdom_base.dt >= min([tot.dt; filt.dt]) & cdom_base.dt <= max([tot.dt; filt.dt]), :);
    cdom = cdom(cdom.dt >= min([tot.dt; filt.dt]) & cdom.dt <= max([tot.dt; filt.dt]), :);
    
    % smooth FDOM data
    cdom_base.mv_fdom = movmean(cdom_base.fdom, 30);
    cdom.mv_fdom = movmean(cdom.fdom, 10);
%     cdom.fdom(~ismember(cdom.dt, cdom_base.dt),:) = NaN;
    % interpolate FDOM onto filt_avg
    filt_avg.cdom_b = interp1(cdom_base.dt, cdom_base.mv_fdom, filt_avg.dt, 'linear');
    filt_avg.cdom_b(isnan(filt_avg.cdom_b)) = interp1(cdom_base.dt, cdom_base.mv_fdom, filt_avg.dt(isnan(filt_avg.cdom_b)), ...
      'nearest', 'extrap');
    filt_avg.cdom = interp1(cdom.dt, cdom.mv_fdom, filt_avg.dt, 'linear');
    filt_avg.cdom(isnan(filt_avg.cdom)) = interp1(cdom.dt, cdom.mv_fdom, filt_avg.dt(isnan(filt_avg.cdom)), ...
      'nearest', 'extrap');
    
    % Use simple mathematical function to interpolate based on CDOM
    n_periods = size(filt_avg,1)-1;
    filt_interp = table(tot.dt, 'VariableNames', {'dt'});
    filt_interp.cdom = interp1(cdom.dt, cdom.mv_fdom, tot.dt, 'linear');
    filt_interp.cdom(isnan(filt_interp.cdom)) = interp1(cdom.dt, ...
      cdom.mv_fdom, tot.dt(isnan(filt_interp.cdom)), 'nearest', 'extrap'); % as independent from tot|filt period
    filt_interp.Tau = NaN(size(filt_interp,1), 1);
    filt_interp.Beamc = NaN(size(filt_interp,1), 1);
    % For each period going from t0 to t1, starting and finishing by a filtered time
    for i=1:n_periods
      it0 = i; it1 = i + 1;
      it = filt_avg.dt(it0) <= filt_interp.dt & filt_interp.dt <= filt_avg.dt(it1);
      if any(it)
        it_filt_interp = filt_interp(it,:);

        % linearly interpolate cdom between filter average
        dt_filtavg_it = filt_avg.dt(it0:it1,:);
        cdom_filtavg_it = filt_avg.cdom(it0:it1,:);
        if sum(isnan(cdom_filtavg_it)) == 1 && ~all(isnan(filt_interp.cdom(it,:)))
          foo = it_filt_interp(~isnan(it_filt_interp.cdom),:);
          cdom_filtavg_it(isnan(cdom_filtavg_it)) = foo.cdom(find(abs(foo.dt - dt_filtavg_it(isnan(cdom_filtavg_it))) == ...
            min(abs(foo.dt - dt_filtavg_it(isnan(cdom_filtavg_it)))), 1, 'first'));
        end
        lin_cdom = interp1(dt_filtavg_it, cdom_filtavg_it, it_filt_interp.dt, 'linear');
        % get cdom variance to the linear interpolation
        Xt = lin_cdom ./ it_filt_interp.cdom;
        % fill missing values in cdom variance by interpolating linearly
        if any(isnan(Xt))
          foo_Xt = [1; Xt; 1];
          foo_Xt_dt = [dt_filtavg_it(1); it_filt_interp.dt; dt_filtavg_it(2)];
          % remove duplicates
          if foo_Xt_dt(1) == foo_Xt_dt(2)
            foo_Xt_dt(1) = foo_Xt_dt(2) - median(diff(foo_Xt_dt), 'omitnan');
          end
          if foo_Xt_dt(end) == foo_Xt_dt(end-1)
            foo_Xt_dt(end) = foo_Xt_dt(end-1) + median(diff(foo_Xt_dt), 'omitnan');
          end
          foo_Xt = fillmissing(foo_Xt, 'linear','SamplePoints',foo_Xt_dt);
          Xt = foo_Xt(2:end-1);
        end
        % linearly interpolate a and c
        lin_Tau = interp1(dt_filtavg_it, filt_avg.Tau(it0:it1,:), it_filt_interp.dt, 'linear');
        lin_Beamc = interp1(dt_filtavg_it, filt_avg.Beamc(it0:it1,:), it_filt_interp.dt, 'linear');
        foo_interp_Tau = NaN(size(lin_Tau));
        foo_interp_Beamc = NaN(size(lin_Beamc));
        
        pos_interp_Tau = lin_Tau ./ Xt;
        pos_interp_Beamc = lin_Beamc ./ Xt;
        neg_interp_Tau = lin_Tau .* Xt;
        neg_interp_Beamc = lin_Beamc .* Xt;
        
        foo_interp_Tau(lin_Tau < 0) = pos_interp_Tau(lin_Tau < 0);
        foo_interp_Beamc(lin_Beamc > 0) = pos_interp_Beamc(lin_Beamc > 0);
        
        foo_interp_Tau(lin_Tau > 0) = neg_interp_Tau(lin_Tau > 0);
        foo_interp_Beamc(lin_Beamc < 0) = neg_interp_Beamc(lin_Beamc < 0);

        filt_interp.Tau(it, :) = foo_interp_Tau;
        filt_interp.Beamc(it, :) = foo_interp_Beamc;
      end
    end
    
    % Note std is interpolated linearly (not using the cdom function)
    filt_interp.Tau_avg_sd = interp1(filt.dt, filt.Tau_avg_sd, filt_interp.dt, 'linear');
    filt_interp.Beamc_avg_sd = interp1(filt.dt, filt.Beamc_avg_sd, filt_interp.dt, 'linear');
    
    % remove interpolation when there is no data
    filt_interp.Tau_avg_sd(all(isnan(tot.Tau),2), :) = NaN;
    filt_interp.Tau(all(isnan(tot.Tau),2), :) = NaN;
    filt_interp.Beamc_avg_sd(all(isnan(tot.Beamc),2), :) = NaN;
    filt_interp.Beamc(all(isnan(tot.Beamc),2), :) = NaN;
    
  case 'linear'
    % Interpolate filtered on total linearly
    filt_interp = table(tot.dt, 'VariableNames', {'dt'});
    filt_interp.Tau = interp1(filt_avg.dt, filt_avg.Tau, filt_interp.dt, 'linear');
    filt_interp.Beamc = interp1(filt_avg.dt, filt_avg.Beamc, filt_interp.dt, 'linear');
    filt_interp.Tau_avg_sd = interp1(filt_avg.dt, filt_avg.Tau_avg_sd, filt_interp.dt, 'linear');
    filt_interp.Beamc_avg_sd = interp1(filt_avg.dt, filt_avg.Beamc_avg_sd, filt_interp.dt, 'linear');
  otherwise
  error('Method not supported.');
end

% Remove lines full of NaNs or with inf data
sel2rm = any(~isfinite(tot.Tau),2) | any(~isfinite(tot.Beamc),2)| all(isnan(tot.Tau),2) | ...
         all(isnan(tot.Beamc),2) | all(isnan(filt_interp.Tau),2) | all(isnan(filt_interp.Beamc),2);
tot(sel2rm,:) = [];
filt_interp(sel2rm,:) = [];

if strcmp(interpolation_method, 'CDOM')
  % remove cdom interpolation when there is no a and c data
  cdom(~ismember(cdom.dt, filt_interp.dt), :) = [];
end

if exist('visFlag', 'file')
  fh = visFlag([], filt_interp, tot, [], filt_avg, [], 'Beamc', round(size(tot.Tau, 2)/2), [], [], [], true);
  title('Check filter event interpolation, press q to continue', 'FontSize', 14)
  if strcmp(interpolation_method, 'CDOM')
    ax1 = gca;
    ax1.YColor = [0	205	205]/255;
    scatter(cdom.dt, cdom.mv_fdom, 15, [0	205	205]/255, 'filled')
    ylabel('FDOM (volts)')
    legend('Filtered interpolated', 'Total', 'Filtered median', 'smoothed FDOM',...
      'AutoUpdate','off', 'FontSize', 12)
  else
    legend('Filtered interpolated', 'Total', 'Filtered median',...
      'AutoUpdate','off', 'FontSize', 12)
  end
  guiSelectOnTimeSeries(fh);
end

% Particulate = Total - FSW
p = table(tot.dt, 'VariableNames', {'dt'});
p.Tau = tot.Tau + filt_interp.Tau;
p.Beamcp = tot.Beamc - filt_interp.Beamc;

% Propagate error (using geometric mean of measurement errors)
%   Note: Error is not propagated through Scattering & Residual temperature
%         correction as required by SeaBASS
p.Tau_sd = sqrt(tot.Tau_avg_sd.^2 + filt_interp.Tau_avg_sd.^2);
p.Beamcp_sd = sqrt(tot.Beamc_avg_sd.^2 + filt_interp.Beamc_avg_sd.^2);
p.Tau_n = tot.Tau_avg_n;
p.Beamcp_n = tot.Beamc_avg_n;

% % delete ap spectrum full of NaNs
% p(all(isnan(p.Taup),2),:) = [];
% p(all(isnan(p.Beamcp),2),:) = [];

%% ag & cg
if ~isempty(di)
  if strcmp(di_method, 'best_di')
    % select DIW with lowest a or c values between 550-650nm
    di_orig = di;
    di_dt = datetime(di_orig.dt, 'ConvertFrom', 'datenum');
    best_di_Tau = NaN(size(di_orig,1), 1);
    best_di_Beamc = NaN(size(di_orig,1), 1);
    for i = 1:size(di_orig,1)
      if i == 1 || i == size(di_orig,1)
        iddi = abs(di_dt(i) - di_dt) < hours(72);
      else
        iddi = abs(di_dt(i) - di_dt) < hours(36);
      end
      highest_di_Tau = di_orig.Tau == max(di_orig.Tau(iddi, :), [], 1);
      lowest_di_Beamc = di_orig.Beamc == min(di_orig.Beamc(iddi, :), [], 1);
      foo_Tau = find(sum(highest_di_Tau, 2) == max(sum(highest_di_Tau, 2)));
      foo_Beamc = find(sum(lowest_di_Beamc, 2) == max(sum(lowest_di_Beamc, 2)));
      di.Tau(i, :) = di_orig.Tau(foo_Tau, :);
      di.Beamc(i, :) = di_orig.Beamc(foo_Beamc, :);
      di.Tau_avg_sd(i, :) = di_orig.Tau_avg_sd(foo_Tau, :);
      di.Beamc_avg_sd(i, :) = di_orig.Beamc_avg_sd(foo_Beamc, :);
      
      best_di_Tau(i) = foo_Tau(1);
      best_di_Beamc(i) = foo_Beamc(1);
    end
  end
  
  % remove when a and c are full of NaNs
  filt_avg(all(isnan(filt_avg.Tau), 2) & all(isnan(filt_avg.Beamc), 2),:) = [];
  
  % Interpolate filtered on Total
  di_interp = table(filt_avg.dt, 'VariableNames', {'dt'});
  di_interp.Tau = interp1(di.dt, di.Tau, di_interp.dt, 'linear', 'extrap');
  di_interp.Beamc = interp1(di.dt, di.Beamc, di_interp.dt, 'linear', 'extrap');
  di_interp.Tau_avg_sd = interp1(di.dt, di.Tau_avg_sd, di_interp.dt, 'linear', 'extrap');
  di_interp.Beamc_avg_sd = interp1(di.dt, di.Beamc_avg_sd, di_interp.dt, 'linear', 'extrap');

  % Dissolved = Filtered - DI
  g = table(filt_avg.dt, 'VariableNames', {'dt'});
  g.Taug = filt_avg.Tau - di_interp.Tau;
  g.Beamcg = filt_avg.Beamc - di_interp.Beamc;

  % Propagate error
  %   Note: Error is not propagated through Scattering & Residual temperature
  %         correction as required by SeaBASS
  g.Taug_sd = sqrt(filt_avg.Tau_avg_sd.^2 + di_interp.Tau_avg_sd.^2);
  g.Beamcg_sd = sqrt(filt_avg.Beamc_avg_sd.^2 + di_interp.Beamc_avg_sd.^2);
  g.Taug_n = filt_avg.Tau_avg_n;
  g.Beamcg_n = filt_avg.Beamc_avg_n;
  
else
  g = table();
end
end




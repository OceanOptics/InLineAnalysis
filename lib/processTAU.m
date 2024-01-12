function [p, g] = processTAU(tot, filt, di, cdom_base, fth, fth_constants, interpolation_method, di_method, days2run)
  %% default wavelength (TODO add lambda as dynamic variable)
  lambda.c = 650;

  % check FTH data
  if ~exist('fth_constants', 'var')
    % Assume most recent FlowControl software
    SWITCH_FILTERED = 1;
    SWITCH_TOTAL = 0;
  else
    SWITCH_FILTERED = fth_constants.SWITCH_FILTERED;
    SWITCH_TOTAL = fth_constants.SWITCH_TOTAL;
  end
  flow_data = fth.qc.tsw;
  
  % round time stamp and remove time duplicates
  tot = round_timestamp(tot);
  filt = round_timestamp(filt);
  flow_data = round_timestamp(flow_data);
  
  % check if fCDOM data loaded if fCDOM interpolation
  if strcmp(interpolation_method, 'CDOM')
    cdom_tbl_name = fieldnames(CDOM.prod);
    if isempty(cdom_tbl_name)
      error('No fDOM prod data loaded')
    end
    % Require both CDOM & Switch position
    if ~isempty(CDOM.prod.(cdom_tbl_name{1}))
      if ~isempty(CDOM.prod.(cdom_tbl_name{1}))
        cdom_base = CDOM.prod.(cdom_tbl_name{1});
      end
      if ~any(cdom_base.dt >= min([tot.dt; filt.dt]) & cdom_base.dt <= max([tot.dt; filt.dt]))
        fprintf('Warning: fCDOM dates do not correspond to ACS dates: interpolation switched to "linear"\n')
        interpolation_method = 'linear';
      else
        % round time stamp and remove time duplicates
        cdom_base = round_timestamp(cdom_base);
      end
    else
      fprintf('Warning: fCDOM data not loaded: interpolation switched to "linear"\n')
      interpolation_method = 'linear';
    end
  end
  % check if TSG data loaded
  if ~isempty(tsg)
    if ~isempty(tsg.prod.a)
      tsg_data = tsg.prod.a;
    elseif ~isempty(tsg.qc.tsw)
      tsg_data = tsg.qc.tsw;
    else
      error('No TSG qc or prod data loaded')
    end
    if ~any(tsg_data.dt >= min([tot.dt; filt.dt]) & tsg_data.dt <= max([tot.dt; filt.dt]))
      fprintf('Warning: TSG dates do not correspond to ACS dates: salinity correction not performed\n')
      tsg_loaded = false;
    else
      tsg_loaded = true;
      % round time stamp and remove time duplicates
      tsg_data = round_timestamp(tsg_data);
    end
  else
    tsg_loaded = false;
  end
  % interpolate flow_data.swt onto binned data to fill missing flow data
  fth_interp = table([tot.dt; flow_data.dt; filt.dt], 'VariableNames', {'dt'});
  [~,b] = sort(fth_interp.dt); % sort dates
  fth_interp.dt = fth_interp.dt(b,:);
  fth_interp.swt = interp1(flow_data.dt, flow_data.swt, fth_interp.dt, 'previous');%, 'linear', 'extrap');
  fth_interp.swt = fth_interp.swt > 0;
  % Find switch events from total to filtered
  sel_start = find(fth_interp.swt(1:end-1) == SWITCH_TOTAL & fth_interp.swt(2:end) == SWITCH_FILTERED);
  % Find switch events from filtered to total
  sel_end = find(fth_interp.swt(1:end-1) == SWITCH_FILTERED & fth_interp.swt(2:end) == SWITCH_TOTAL);
  % Verify selections of filtered period
  if sel_start(1) > sel_end(1); sel_end(1) = []; end
  if sel_start(end) > sel_end(end); sel_end(end+1) = size(fth_interp.swt,1); end
  if size(sel_start,1) ~= size(sel_end,1); error('Inconsistent fth data'); end

  %%%%%%%%%%%%%%%%%%%%%%% Compute filtered average %%%%%%%%%%%%%%%%%%%%%%%
  % TODO maybe: write a filter_average function that can be use with every instrument

  % create filter average data table
  filt_avg = table(NaN(size(sel_start)), 'VariableNames', {'dt'});
  filt_avg.start = fth_interp.dt(sel_start);
  filt_avg.end = fth_interp.dt(sel_end);
  filt_avg.Tau = NaN(size(filt_avg,1), 1);
  filt_avg.Tau_avg_sd = NaN(size(filt_avg,1), 1);
  filt_avg.Tau_avg_n = NaN(size(filt_avg,1), 1);
  filt_avg.Beamc = NaN(size(filt_avg,1), 1);
  filt_avg.Beamc_avg_sd = NaN(size(filt_avg,1), 1);
  filt_avg.Beamc_avg_n = NaN(size(filt_avg,1), 1);

  % create filter interpolation table
  filt_interp = table(tot.dt, 'VariableNames', {'dt'});
  filt_interp.Tau = NaN(size(tot.dt));
  filt_interp.Beamc = NaN(size(tot.dt));

  % add T/S to filt and filt_interp if TSG loaded
  if tsg_loaded
    % interpolate linear T and S on filt_interp table
    filt_interp = addvars(filt_interp, interp_extrap(tsg_data, filt_interp.dt, tsg.temperature_variable, 30, true, 'linear', 'nearest'), ...
      interp_extrap(tsg_data, filt_interp.dt, 's', 30, true, 'linear', 'nearest'), 'NewVariableNames', {'t','s'}, 'After', 'dt');
    % interpolate T and S on filtered data
    filt = addvars(filt, interp_extrap(tsg_data, filt.dt, tsg.temperature_variable, 30, true, 'linear', 'nearest'), ...
      interp_extrap(tsg_data, filt.dt, 's', 30, true, 'linear', 'nearest'), 'NewVariableNames', {'t','s'}, 'After', 'dt');
    % add T and S variables on filt_avg table
    filt_avg = addvars(filt_avg, NaN(size(filt_avg,1),1), NaN(size(filt_avg,1),1), 'NewVariableNames', {'t','s'}, 'After', 'dt');
    % add T and S variables on filt_avg table
    filt_avg = addvars(filt_avg, NaN(size(filt_avg,1),1), 'NewVariableNames', 't', 'After', 'dt');
    filt_avg = addvars(filt_avg, NaN(size(filt_avg,1),1), 'NewVariableNames', 's', 'After', 'dt');
    % remove psiS and psiT from tot and filt with Tref = median(T) of whole processed section
    % (NO CORRECTION FOR TAU)
    Tref = median([filt.t; filt_interp.t], 'omitnan');
    c_psiS = interp1(psi.wl, psi.c_psiS, lambda.c, 'spline');
    c_psiT = interp1(psi.wl, psi.psiT, lambda.c, 'spline');
    % Correct independently only when t and s are available to avoid getting NaN when t or s are NaN
    %%%%%%%%%%%% TODO add flag when filt or tot are not corrected (t or s not available) %%%%%%%%%%%%
    filt_t_nan = isnan(filt.t);
    filt_s_nan = isnan(filt.s);
    tot_t_nan = isnan(filt_interp.t);
    tot_s_nan = isnan(filt_interp.s);
    filt.Beamc(~filt_t_nan, :) = filt.Beamc(~filt_t_nan, :) - (c_psiT .* (filt.t(~filt_t_nan) - Tref));
    filt.Beamc(~filt_s_nan, :) = filt.Beamc(~filt_s_nan, :) - (c_psiS .* filt.s(~filt_s_nan));
    tot.Beamc(~tot_t_nan, :) = tot.Beamc(~tot_t_nan, :) - (c_psiT .* (filt_interp.t(~tot_t_nan) - Tref));
    tot.Beamc(~tot_s_nan, :) = tot.Beamc(~tot_s_nan, :) - (c_psiS .* filt_interp.s(~tot_s_nan));
    % filt.Beamc = filt.Beamc - (c_psiS .* filt.s + c_psiT .* (filt.t - Tref)); % old code
    % tot.Beamc = tot.Beamc - (c_psiS .* filt_interp.s + c_psiT .* (filt_interp.t - Tref)); % old code
  end

  % prepare fCDOM data if CDOM interpolation
  if strcmp(interpolation_method, 'CDOM')
    % interpolate spline and extrapolate nearest fdom on filt_interp table
    filt_interp = addvars(filt_interp, interpspline_extrapnearest(cdom_base, filt_interp.dt, 'fdom', 30, false), ...
      'NewVariableNames', 'fdom', 'After', 'dt');
    % interpolate fdom on filtered data
    filt = addvars(filt, interpspline_extrapnearest(cdom_base, filt.dt, 'fdom', 30, false), ...
      'NewVariableNames', 'fdom', 'After', 'dt');
    % add fdom variables on filt_avg table
    filt_avg = addvars(filt_avg, NaN(size(filt_avg,1),1), 'NewVariableNames', 'fdom', 'After', 'dt');
  end

  for i=1:size(sel_start, 1)
    sel_filt = fth_interp.dt(sel_start(i)) <= filt.dt & filt.dt <= fth_interp.dt(sel_end(i));
    foo = filt(sel_filt,:);
    if sum(sel_filt) == 1
      filt_avg.dt(i) = foo.dt;
      filt_avg.Tau(i,:) = foo.Tau;
      filt_avg.Beamc(i,:) = foo.Beamc;
      filt_avg.Tau_avg_sd(i,:) = foo.Tau_avg_sd;
      filt_avg.Beamc_avg_sd(i,:) = foo.Beamc_avg_sd;
      filt_avg.Tau_avg_n(i) = foo.Tau_avg_n;
      filt_avg.Beamc_avg_n(i) = foo.Beamc_avg_n;
      % repeat for T/S if tsg is loaded
      if tsg_loaded
        filt_avg.t(i,:) = foo.t;
        filt_avg.s(i,:) = foo.s;
      end
      % repeat for fdom if interpolation method == CDOM
      if strcmp(interpolation_method, 'CDOM')
        filt_avg.fdom(i,:) = foo.fdom;
      end
    else
      Tau_perc75 = foo.Tau < prctile(foo.Tau, 75, 1);
      Beamc_perc25 = foo.Beamc > prctile(foo.Beamc, 25, 1);
      foo.Tau_avg_sd(Tau_perc75) = NaN;
      foo.Beamc_avg_sd(Beamc_perc25) = NaN;
      foo.Tau(Tau_perc75) = NaN;
      foo.Beamc(Beamc_perc25) = NaN;
      % compute average of all values smaller than 25th percentile for each filter event
      filt_avg.dt(i) = mean(foo.dt(any(~a_perc25, 2) | any(~c_perc25, 2)), 'omitnan');
      if any(any(any(~a_perc25, 2) | any(~c_perc25, 2), 2))
        filt_avg.start(i) = min(foo.dt(any(~a_perc25, 2) | any(~c_perc25, 2)));
        filt_avg.end(i) = max(foo.dt(any(~a_perc25, 2) | any(~c_perc25, 2)));
      else
        filt_avg.start(i) = NaN;
        filt_avg.end(i) = NaN;
      end
      filt_avg.dt(i) = mean(foo.dt(any(~Tau_perc75, 2)), 'omitnan');
      filt_avg.Tau(i,:) = mean(foo.Tau, 1, 'omitnan');
      filt_avg.Beamc(i,:) = mean(foo.Beamc, 1, 'omitnan');
      filt_avg.Tau_avg_sd(i,:) = mean(foo.Tau_avg_sd, 1, 'omitnan');
      filt_avg.Beamc_avg_sd(i,:) = mean(foo.Beamc_avg_sd, 1, 'omitnan');
      filt_avg.Tau_avg_n(i) = sum(foo.Tau_avg_n(any(~isnan(foo.Tau), 2)), 'omitnan');
      filt_avg.Beamc_avg_n(i) = sum(foo.Beamc_avg_n(any(~isnan(foo.Beamc), 2)), 'omitnan');
      % repeat for T/S if tsg is loaded
      if tsg_loaded
        foo.t(foo.t > prctile(foo.t, 25, 1)) = NaN;
        foo.s(foo.s > prctile(foo.s, 25, 1)) = NaN;
        % compute average of all values smaller than 25th percentile for each filter event
        filt_avg.t(i) = mean(foo.t, 1, 'omitnan');
        filt_avg.s(i) = mean(foo.s, 1, 'omitnan');
      end
      % repeat for fdom if interpolation method == CDOM
      if strcmp(interpolation_method, 'CDOM')
        foo.fdom(foo.fdom > prctile(foo.fdom, 25, 1)) = NaN;
        % compute average of all values smaller than 25th percentile for each filter event
        filt_avg.fdom(i) = mean(foo.fdom, 1, 'omitnan');
      end
    end
  end
  % remove empty filter events
  rm_filter_event = all(isnan(filt_avg.Tau), 2) | all(isnan(filt_avg.Beamc), 2);
  % sel_start(rm_filter_event) = [];
  % sel_end(rm_filter_event) = [];
  filt_avg(rm_filter_event, :) = [];
    
  
  switch interpolation_method
    case 'CDOM'
      % interpolate ag and cg using fcdom
      filt_interp = taugbeamcg_fdom_interpolation(tot, filt_interp, filt, CDOM.dark, filt_avg);
  
      % TODO: propagate error from regression and filter event STD
      filt_interp.Tau_avg_sd = interp1(filt.dt, filt.Tau_avg_sd, filt_interp.dt, 'linear');
      filt_interp.Beamc_avg_sd = interp1(filt.dt, filt.Beamc_avg_sd, filt_interp.dt, 'linear');
      
      % remove interpolation when there is no data
      filt_interp.Tau_avg_sd(all(isnan(tot.Tau),2), :) = NaN;
      filt_interp.Tau(all(isnan(tot.Tau),2), :) = NaN;
      filt_interp.Beamc_avg_sd(all(isnan(tot.Beamc),2), :) = NaN;
      filt_interp.Beamc(all(isnan(tot.Beamc),2), :) = NaN;
      
    case 'linear'
      % Interpolate filtered on total linearly
      filt_interp.Tau = interp1(filt_avg.dt, filt_avg.Tau, filt_interp.dt, 'linear');
      filt_interp.Tau = fillmissing(filt_interp.Tau, 'nearest');
      filt_interp.Beamc = interp1(filt_avg.dt, filt_avg.Beamc, filt_interp.dt, 'linear');
      filt_interp.Beamc = fillmissing(filt_interp.Beamc, 'nearest');
      filt_interp.Tau_avg_sd = interp1(filt_avg.dt, filt_avg.Tau_avg_sd, filt_interp.dt, 'linear');
      filt_interp.Tau_avg_sd = fillmissing(filt_interp.Tau_avg_sd, 'nearest');
      filt_interp.Beamc_avg_sd = interp1(filt_avg.dt, filt_avg.Beamc_avg_sd, filt_interp.dt, 'linear');
      filt_interp.Beamc_avg_sd = fillmissing(filt_interp.Beamc_avg_sd, 'nearest');
      regression_stats = struct();
    otherwise
    error('Method not supported.');
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
            cdom_filtavg_it(isnan(cdom_filtavg_it)) = foo.Beamcdom(find(abs(foo.dt - dt_filtavg_it(isnan(cdom_filtavg_it))) == ...
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
  
  if exist('visFlag', 'file')
    % id only day to run in all tables to plot
    filt_interp_id = filt_interp.dt >= min(days2run) & filt_interp.dt < max(days2run)+1;
    tot_id = tot.dt >= min(days2run) & tot.dt < max(days2run)+1;
    filt_avg_id = filt_avg.dt >= min(days2run) & filt_avg.dt < max(days2run)+1;
    % plot
    fh = visFlag([], filt_interp(filt_interp_id, :), tot(tot_id, :), [], filt_avg(filt_avg_id, :), [], 'Beamc', round(size(tot.Tau, 2)/2), [], [], [], true);
    title('Check filter event interpolation, press q to continue', 'FontSize', 14)
    if strcmp(interpolation_method, 'CDOM')
      ax1 = gca;
      ax1.YColor = [0	205	205]/255;
      scatter(filt_interp.dt(filt_interp_id), filt_interp.fdom(filt_interp_id), '.', 'MarkerEdgeColor', [0	205	205]/255, 'MarkerEdgeAlpha', 0.5)
      plot(filt_avg.dt(filt_avg_id), filt_avg.fdom(filt_avg_id), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0	205	205]/255)
      ylabel('FDOM (volts)')
      legend('Filtered interpolated', 'Total', 'Filtered percentile average', 'fCDOM', ...
        'fCDOM filtered percentile average', 'AutoUpdate','off', 'FontSize', 12)
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


function [filt_interp] = taugbeamcg_fdom_interpolation(tot, filt_interp, filt, cdom_dark, filt_avg)
  % Reconstruct ag and cg from fCDOM data
  [filt_avg, filt_interp, ~] = fdom_taugbeamcg_model(filt, filt_interp, filt_avg); % regression_stats
  n_periods = size(filt_avg,1)-1;

  % interpolate ag and cg based on fdom: Guillaume's method (global)
  % allocate new variables
  filt_interp.Taug_fact = NaN(size(filt_interp.Tau));
  filt_interp.Beamcg_fact = NaN(size(filt_interp.Beamc));
  filt_interp.fCDOM_mix_Tau_cluster = false(size(filt_interp.dt));
  filt_interp.fCDOM_mix_Beamc_cluster = false(size(filt_interp.dt));
  filt_interp.fCDOM_Tau_cluster_chg = false(size(filt_interp.dt));
  filt_interp.fCDOM_Beamc_cluster_chg = false(size(filt_interp.dt));
  filt_interp.flag_linear_interp = false(size(filt_interp.dt));

  %%%%% REWRITTEN ON 2023-12-29
  % linearly interpolate a, c, and fdom
  filt_interp.lin_Tau = interp1(filt_avg.dt, filt_avg.a, filt_interp.dt, 'linear');
  filt_interp.lin_Beamc = interp1(filt_avg.dt, filt_avg.c, filt_interp.dt, 'linear');
  filt_interp.lin_fdom = interp1(filt_avg.dt, filt_avg.fdom, filt_interp.dt, 'linear');
  filt_interp.lin_Tau = fillmissing(filt_interp.lin_Tau, 'nearest');
  filt_interp.lin_Beamc = fillmissing(filt_interp.lin_Beamc, 'nearest');
  filt_interp.lin_fdom = fillmissing(filt_interp.lin_fdom, 'nearest');
  for i=1:n_periods
    it0 = i; it1 = i + 1;
    it = filt_avg.dt(it0) <= filt_interp.dt & filt_interp.dt <= filt_avg.dt(it1);
    if any(it)
      it_filt_interp = filt_interp(it,:);
      % interpolation based on fCDOM if enough dynamic range in fCDOM
      if any((filt_interp.fdom(it) - cdom_dark) / cdom_dark > -0.05*2000)
        % compute interpolated ag factor
        Taug_fact0 = (filt_avg.a(it0,:) - filt_avg.intercept_Tau(it0,:)) ./ filt_avg.fdom(it0);
        Taug_fact1 = (filt_avg.a(it1,:) - filt_avg.intercept_Tau(it1,:)) ./ filt_avg.fdom(it1);
        filt_interp.Taug_fact(it,:) = interp1(filt_avg.dt(it0:it1,:), [abs(Taug_fact0); abs(Taug_fact1)], it_filt_interp.dt, 'linear');
        % compute interpolated cg factor
        Beamcg_fact0 = (filt_avg.c(it0,:) - filt_avg.intercept_Beamc(it0,:)) ./ filt_avg.fdom(it0);
        Beamcg_fact1 = (filt_avg.c(it1,:) - filt_avg.intercept_Beamc(it1,:)) ./ filt_avg.fdom(it1);
        filt_interp.Beamcg_fact(it,:) = interp1(filt_avg.dt(it0:it1,:), [Beamcg_fact0; Beamcg_fact1], it_filt_interp.dt, 'linear');

        % flag when mixed a clusters within filter event 0 or filter event 1
        if any(filt_avg.cluster_Tau_flag(it0) | filt_avg.cluster_Tau_flag(it1))
          fprintf('Tau filter interpolation #%i flagged: mixed cluster within start and/or end filter event(s)\n', i)
          filt_interp.fCDOM_mix_Tau_cluster(i) = true;
        end
        % flag when mixed c clusters within filter event 0 or filter event 1
        if any(filt_avg.cluster_Beamc_flag(it0) | filt_avg.cluster_Beamc_flag(it1))
          fprintf('Beamc filter interpolation #%i flagged: mixed cluster within start and/or end filter event(s)\n', i)
          filt_interp.fCDOM_mix_Beamc_cluster(i) = true;
        end
        % flag when a clusters change between filter event 0 and filter event 1
        if ~all(isnan(filt_avg.cluster_weighted_Tau(it0,:))) && ~all(isnan(filt_avg.cluster_weighted_Tau(it1,:))) && ...
            any(filt_avg.cluster_weighted_Tau(it0,:) ~= filt_avg.cluster_weighted_Tau(it1, :)) && any(any(~isnan(tot.Tau(it, :)), 2))
          fprintf('Tau filter interpolation #%i flagged: change of clusters between start and end filter events\n', i)
          filt_interp.fCDOM_Tau_cluster_chg(it) = true;
        end
        % flag when c clusters change between filter event 0 and filter event 1
        if ~all(isnan(filt_avg.cluster_weighted_Beamc(it0,:))) && ~all(isnan(filt_avg.cluster_weighted_Beamc(it1,:))) && ...
            any(filt_avg.cluster_weighted_Beamc(it0,:) ~= filt_avg.cluster_weighted_Beamc(it1, :)) && any(any(~isnan(tot.Beamc(it, :)), 2))
          fprintf('Beamc filter interpolation #%i flagged: change of clusters between start and end filter events\n', i)
          filt_interp.fCDOM_Beamc_cluster_chg(it) = true;
        end
        % set linear interpolation boolean
        filt_interp.flag_linear_interp(it) = false;
      else
        % linearly interpolate ag and cg if fdom is constant between it0 and it1
        fprintf('Filter interpolation #%i flagged: linearly interpolated\n', i)
        filt_interp.flag_linear_interp(it) = true;
      end
    end
  end
  % fill missing Taug_fact and Beamcg_fact
  filt_interp.Taug_fact = fillmissing(filt_interp.Taug_fact, 'nearest');
  filt_interp.Beamcg_fact = fillmissing(filt_interp.Beamcg_fact, 'nearest');
  % compute Taug and Beamcg
  filt_interp.Tau(~filt_interp.flag_linear_interp, :) = filt_interp.lin_Tau(~filt_interp.flag_linear_interp, :) + ...
    filt_interp.Taug_fact(~filt_interp.flag_linear_interp, :) .* ...
    (filt_interp.fdom(~filt_interp.flag_linear_interp) - filt_interp.lin_fdom(~filt_interp.flag_linear_interp)); 
  filt_interp.Beamc(~filt_interp.flag_linear_interp, :) = filt_interp.lin_Beamc(~filt_interp.flag_linear_interp, :) + ...
    filt_interp.Beamcg_fact(~filt_interp.flag_linear_interp, :) .* ...
    (filt_interp.fdom(~filt_interp.flag_linear_interp) - filt_interp.lin_fdom(~filt_interp.flag_linear_interp));
  % linear interpolation when fdom interpolation was not used
  filt_interp.Tau(filt_interp.flag_linear_interp, :) = filt_interp.lin_Tau(filt_interp.flag_linear_interp, :);
  filt_interp.Beamc(filt_interp.flag_linear_interp, :) = filt_interp.lin_Beamc(filt_interp.flag_linear_interp, :);
end

%% 
function [filt_avg, filt_interp, regress_stats] = fdom_taugbeamcg_model(filt, filt_interp, filt_avg)
  
  clusters = table();
  % iterate dbscan with multiple epsilon until nRMSE of all cluster < 10%
  fprintf('Finding best epsilon for dbscan Kernel Density clustering of a(filt)/fCDOM(filt) ... ')
  k = 1;
  epsilons = 1:-0.01:0.01;
  keep_going = true;
  max_nb_clusters = day(datetime(max(filt_interp.dt) - min(filt_interp.dt), 'ConvertFrom', 'datenum'));
  eval_cluster = table();
  eval_cluster.nRMSE_Tau = NaN(length(epsilons), max_nb_clusters);
  eval_cluster.nRMSE_Beamc = NaN(length(epsilons), max_nb_clusters);
  eval_cluster.slope_Tau = NaN(length(epsilons), max_nb_clusters);
  eval_cluster.slope_Beamc = NaN(length(epsilons), max_nb_clusters);
  epsilon_Tau = epsilons(1);
  epsilon_Beamc = epsilons(1);
  % normalize the time variable between 0 to 0.5
  time_var = (filt.dt-min(filt.dt))/max(filt.dt-min(filt.dt)) * 0.5;
  data_tocluster_Tau = [filt.Tau./filt.fdom time_var]; % time_var filt.t filt.s
  data_tocluster_Beamc = [filt.Beamc./filt.fdom time_var]; % time_var filt.t filt.s
  % find minimum cluster size depending on the number of days loaded
  % round(size(time_var, 1) / (max(day(datetime(filt.dt-min(filt.dt)+1, 'ConvertFrom', 'datenum'))) * 24))
  minimum_nb_point_per_cluster = 60;
  while k <= length(epsilons) && keep_going
    clusters.Tau = dbscan(data_tocluster_Tau, epsilons(k), minimum_nb_point_per_cluster, 'Distance', 'euclidean');
    % clusters.Tau(clusters.Tau == -1) = NaN;
    clusters.Beamc = dbscan(data_tocluster_Beamc, epsilons(k), minimum_nb_point_per_cluster, 'Distance', 'euclidean');
    clusters.Beamc(clusters.Beamc == -1) = NaN;
    % regress a&c with fdom
    [regress_stats_Tau, regress_stats_Beamc] = regress_acfilt(filt.Tau, filt.Beamc, filt.fdom, clusters);
    eval_cluster.nRMSE_Tau(k, 1:size(regress_stats_Tau.nRMSE, 2)) = regress_stats_Tau.nRMSE;
    eval_cluster.nRMSE_Beamc(k, 1:size(regress_stats_Beamc.nRMSE, 2)) = regress_stats_Beamc.nRMSE;
    eval_cluster.slope_Tau(k, 1:size(regress_stats_Tau.slope, 2)) = regress_stats_Tau.slope;
    eval_cluster.slope_Beamc(k, 1:size(regress_stats_Beamc.slope, 2)) = regress_stats_Beamc.slope;
    if k >= 2
      % get epsilon_a
      if all(eval_cluster.nRMSE_Tau(k,1:size(regress_stats_Tau.nRMSE, 2)) < 10, 2)
        epsilon_Tau = epsilons(k);
      end
      % get epsilon_c
      if all(eval_cluster.nRMSE_Beamc(k,1:size(regress_stats_Beamc.nRMSE, 2)) < 10, 2)
        epsilon_Beamc = epsilons(k);
      end
      % if both a and c found stop iterations
      if all(eval_cluster.nRMSE_Tau(k,1:size(regress_stats_Tau.nRMSE, 2)) < 10, 2) && ...
           all(eval_cluster.nRMSE_Beamc(k,1:size(regress_stats_Beamc.nRMSE, 2)) < 10, 2)
        keep_going = false;
      elseif size(regress_stats_Tau.nRMSE, 2) > max_nb_clusters || size(regress_stats_Beamc.nRMSE, 2) > max_nb_clusters
        keep_going = false;
        if epsilon_Tau == epsilons(1)
          epsilon_Tau = epsilons(k);
        end
        if epsilon_Beamc == epsilons(1)
          epsilon_Beamc = epsilons(k);
        end
      end
    end
    k = k + 1;
  end
  fprintf('done\n')

  % run clustering with the best epsilon
  fprintf('Kernel Density clustering of a(filt)/fCDOM(filt) ... ')
  clusters = table();
  clusters.Tau = dbscan(data_tocluster_Tau, epsilon_Tau, minimum_nb_point_per_cluster, 'Distance', 'euclidean');
  clusters.Tau(clusters.Tau == -1) = NaN;
  clusters.Beamc = dbscan(data_tocluster_Beamc, epsilon_Beamc, minimum_nb_point_per_cluster, 'Distance', 'euclidean');
  clusters.Beamc(clusters.Beamc == -1) = NaN;
  fprintf('done\n')

  clusters.Tau = fillmissing(clusters.Tau, 'nearest');
  clusters.Beamc = fillmissing(clusters.Beamc, 'nearest');

  % regress a&c with fdom
  regress_stats = struct();
  [regress_stats.Tau, regress_stats.Beamc] = regress_taubeamcfilt(filt.Tau, filt.Beamc, filt.fdom, clusters);

  % plot clusters
  fprintf('Plotting clustering results ... ')
  figure(33); subplot(2, 2, 1);
  gscatter(datetime(filt.dt, 'ConvertFrom', 'datenum'), filt.Tau ./ filt.fdom, ...
    clusters.Tau, [], 'vo<s>pdh+*x', 8, 'on', 'Time', 'a/fdom');
  title('Tau cluster time series'); set(gca, 'Fontsize', 12)
  subplot(2, 2, 2); gsc = gscatter(filt.fdom, filt.Tau, clusters.Tau, [], 'vo<s>pdh+*x');
  xlimit = get(gca, 'XLim');
  ylimit = get(gca, 'YLim');
  for c = 1:size(regress_stats.Tau.slope, 2)
    rl = refline(regress_stats.Tau.slope(:, c), regress_stats.Tau.intercept(:, c));
    set(rl, 'Color', gsc(c).Color, 'LineWidth', 1)
  end
  set(gca, 'XLim', xlimit)
  set(gca, 'YLim', ylimit)
  xlabel('fdom'); ylabel('\tau'); title('Tau (\tau) cluster: a/fdom'); set(gca, 'Fontsize', 12)
  leg = findobj(gcf, 'Type', 'Legend');
  title(leg,'Clusters')

  subplot(2, 2, 3);
  gscatter(datetime(filt.dt, 'ConvertFrom', 'datenum'), filt.Beamc ./ filt.fdom, ...
    clusters.Beamc, [], 'vo<s>pdh+*x', 8, 'on', 'Time', 'c/fdom');
  title('Beam-c cluster time series'); set(gca, 'Fontsize', 12)
  subplot(2, 2, 4); gsc = gscatter(filt.fdom, filt.Beamc, clusters.Beamc, [], 'vo<s>pdh+*x');
  xlimit = get(gca, 'XLim');
  ylimit = get(gca, 'YLim');
  for c = 1:size(regress_stats.Beamc.slope, 2)
    rl = refline(regress_stats.Beamc.slope(:, c), regress_stats.Beamc.intercept(:, c));
    set(rl, 'Color', gsc(c).Color, 'LineWidth', 1)
  end
  set(gca, 'XLim', xlimit)
  set(gca, 'YLim', ylimit)
  xlabel('fdom'); ylabel('Beam-c'); title('Beam-c cluster: Beam-c/fdom'); set(gca, 'Fontsize', 12)
  leg = findobj(gcf, 'Type', 'Legend');
  leg(1).String = strrep(leg(1).String, 'data', 'regression cluster ');
  leg(3).String = strrep(leg(3).String, 'data', 'regression cluster ');
  title(leg,'Clusters')
  fprintf('done\n')

  drawnow

  %%% apply the regression to get Tau and Beamc dissolved from fdom
  groups_Tau = unique(clusters.Tau(~isnan(clusters.Tau)));
  groups_Beamc = unique(clusters.Beamc(~isnan(clusters.Beamc)));

  % Weight cluster for each filter event in case clusters are changing
  % during a filter event (unlikely with the new dbscan method but let's keep it in case)
  filt_avg = addvars(filt_avg, NaN(size(filt_avg,1),size(groups_Tau,1)), NaN(size(filt_avg,1),size(groups_Beamc,1)), ...
    'NewVariableNames', {'cluster_weighted_Tau', 'cluster_weighted_Beamc'}, 'after', 'end');
  filt_avg = addvars(filt_avg, false(size(filt_avg.dt)), false(size(filt_avg.dt)), ...
    'NewVariableNames', {'cluster_Beamc_flag','cluster_Tau_flag'}, 'after', 'end');
  for i=1:size(filt_avg, 1)
    sel_filt = filt_avg.start(i) <= filt.dt & filt.dt <= filt_avg.end(i);
    foo = clusters(sel_filt,:);
    % flag if clusters are not all the same within a single filter event
    if size(unique(foo.Tau(~isnan(foo.Tau))), 1) > 1
      filt_avg.cluster_Tau_flag(i) = true;
    end
    if size(unique(foo.Beamc(~isnan(foo.Beamc))), 1) > 1
      filt_avg.cluster_Beamc_flag(i) = true;
    end
    % weight clusters proportionally for each filter events
    for j = 1:size(groups_Tau,1)
      filt_avg.cluster_weighted_Tau(i,j) = sum(foo.Tau == groups_Tau(j)) ./ sum(~isnan(foo.Tau));
    end
    for j = 1:size(groups_Beamc,1)
      filt_avg.cluster_weighted_Beamc(i,j) = sum(foo.Beamc == groups_Beamc(j)) ./ sum(~isnan(foo.Beamc));
    end
  end
  % weight a slope and intercept
  filt_avg.slope_Tau = zeros(size(filt_avg, 1), size(regress_stats.Tau, 1));
  filt_avg.intercept_Tau = zeros(size(filt_avg, 1), size(regress_stats.Tau, 1));
  for j = 1:size(groups_Tau, 1)
    filt_avg.slope_Tau = filt_avg.slope_Tau + regress_stats.Tau.slope(:,j)'.*filt_avg.cluster_weighted_Tau(:,j);
    filt_avg.intercept_Tau = filt_avg.intercept_Tau + regress_stats.Tau.intercept(:,j)'.*filt_avg.cluster_weighted_Tau(:,j);
  end
  % weight c slope and intercept
  filt_avg.slope_Beamc = zeros(size(filt_avg, 1), size(regress_stats.Beamc, 1));
  filt_avg.intercept_Beamc = zeros(size(filt_avg, 1), size(regress_stats.Beamc, 1));
  for j = 1:size(groups_Beamc, 1)
    filt_avg.slope_Beamc = filt_avg.slope_Beamc + regress_stats.Beamc.slope(:,j)'.*filt_avg.cluster_weighted_Beamc(:,j);
    filt_avg.intercept_Beamc = filt_avg.intercept_Beamc + regress_stats.Beamc.intercept(:,j)'.*filt_avg.cluster_weighted_Beamc(:,j);
  end
  % interpolate slope and intercept cluster specific onto filt_interp.dt
  filt_interp.slope_interp_Tau = interpspline_extrapnearest(filt_avg, filt_interp.dt, 'slope_Tau');
  filt_interp.intercept_interp_a = interpspline_extrapnearest(filt_avg, filt_interp.dt, 'intercept_Tau');
  filt_interp.slope_interp_Beamc = interpspline_extrapnearest(filt_avg, filt_interp.dt, 'slope_Beamc');
  filt_interp.intercept_interp_c = interpspline_extrapnearest(filt_avg, filt_interp.dt, 'intercept_Beamc');
end


%% regression between a/c and other variable in filter events
function [tau_reg, c_reg] = regress_acfilt(Tau, Beamc, ancillary, clusters)
  if nargin < 4
    clusters = table();
    clusters.Tau = ones(size(Tau, 1), 1);
    clusters.Beamc = ones(size(Beamc, 1), 1);
  end
  % robust linear regression between filt a&c and ancillary variable
  tau_reg = array2table(NaN(size(Tau, 2), 5), 'VariableNames', {'slope', 'intercept', 'RMSE', 'nRMSE', 'R2'});
  c_reg = array2table(NaN(size(Beamc, 2), 5), 'VariableNames', {'slope', 'intercept', 'RMSE', 'nRMSE', 'R2'});
  groups_Tau = unique(clusters.Tau(~isnan(clusters.Tau)));
  groups_Beamc = unique(clusters.Beamc(~isnan(clusters.Beamc)));
  for i = 1:size(Tau, 2)
    for j = 1:size(groups_Tau, 1)
      id_grp_tau = clusters.Tau == groups_Tau(j);
      % regress variable with filtered water absorption
      [stats.b, stats.stats] = robustfit(ancillary(id_grp_tau), Tau(id_grp_tau, i));
      [~, MSGID] = lastwarn();
      warning('off', MSGID)
      tau_reg.slope(i, j) = stats.b(2);
      tau_reg.intercept(i, j) = stats.b(1);
      tau_reg.RMSE(i, j) = stats.stats.robust_s;
      % calculated normalize RMSE by average and R2
      id_nonan = all(~isnan([Tau(:, i) ancillary]), 2);
      tau_reg.nRMSE(i, j) = tau_reg.RMSE(i, j)/abs(mean(Tau(id_grp_tau & id_nonan, i), 'omitnan'))*100;
      tau_reg.R2(i, j) = corr(Tau(id_grp_tau & id_nonan, i), stats.b(1) + stats.b(2) * ancillary(id_grp_tau & id_nonan))^2;
    end
    for j = 1:size(groups_Beamc, 1)
      id_grp_c = clusters.Beamc == groups_Beamc(j);
      % regress variable with filtered water attenuation
      [stats.b, stats.stats] = robustfit(ancillary(id_grp_c), Beamc(id_grp_c, i));
      [~, MSGID] = lastwarn();
      warning('off', MSGID)
      c_reg.slope(i, j) = stats.b(2);
      c_reg.intercept(i, j) = stats.b(1);
      c_reg.RMSE(i, j) = stats.stats.robust_s;
      % calculated normalize RMSE by average and R2
      id_nonan = all(~isnan([Beamc(:, i) ancillary]), 2);
      c_reg.nRMSE(i, j) = c_reg.RMSE(i, j)/abs(mean(Beamc(id_grp_c & id_nonan, i), 'omitnan'))*100;
      c_reg.R2(i, j) = corr(Beamc(id_grp_c & id_nonan, i), stats.b(1) + stats.b(2) * ancillary(id_grp_c & id_nonan))^2;
    end
  end
end


%%
function data_out = interpspline_extrapnearest(data_in, dt_vector, var_tointerp, max_missing_length, extrap_bool)
  if nargin < 5
    extrap_bool = true;
  end
  % remove row full of NaNs
  data_in(all(isnan(data_in.(var_tointerp)), 2), :) = [];
  % convert dt_vector to datetime
  datetime_vector = datetime(dt_vector, 'ConvertFrom', 'datenum');
  dt = (datenum(min(datetime_vector):median(diff(datetime_vector)):max(datetime_vector)))';
  % id extrapolation
  extrapolated_id = isnan(interp1(data_in.dt, data_in.(var_tointerp), dt, 'nearest'));
  % interpolate spline
  data_out = interp1(data_in.dt, data_in.(var_tointerp), dt, 'spline');
  % replace extrapolated data by NaN
  data_out(extrapolated_id) = NaN;
  % fill missing data with nearest interpolation
  if extrap_bool
    data_out = fillmissing(data_out, 'nearest');
  end
  % replace interpolated values over gaps > max_missing_length by NaN
  if nargin >= 4
    if any(isnan(data_out))
      if ~isempty(max_missing_length)
        missing_data = ~ismember(dt, data_in.dt);
        t = [true; diff(missing_data) ~= 0];
        k = diff(find([t; true])) .* missing_data(t);
        long_nan = k(cumsum(t)) > max_missing_length;
        data_out(long_nan, :) = NaN;
      end
    end
  end
  data_out = interp1(dt, data_out, dt_vector, 'nearest');
end


%%
function data_out = round_timestamp(data_in)
  data_out = table();
  % make sure data_in.dt in rounded to the time binning frequency
  datetime_data_in_dt = datetime(data_in.dt, 'ConvertFrom', 'datenum');
  % get time binning frequency
  Tbin_data_in = median(diff(datetime_data_in_dt));
  if Tbin_data_in >= hours(1)
    data_out.dt = datenum(dateshift(datetime_data_in_dt, 'start', 'hours'));
  elseif Tbin_data_in >= minutes(1)
    % round start/end time to minute
    data_out.dt = datenum(dateshift(datetime_data_in_dt, 'start', 'minutes'));
  elseif Tbin_data_in >= seconds(1)
    % round start/end time to seconds
    data_out.dt = datenum(dateshift(datetime_data_in_dt, 'start', 'seconds'));
  else
    error('automatic detection of sampling rate detected a frequency not supported: check interpspline_extrapnearest function in processACS.m')
  end
  % remove duplicates
  [~, L, ~] = unique(data_in.dt,'first');
  indexToDump = not(ismember(1:numel(data_in.dt), L));
  if sum(indexToDump) > 0
    data_in(indexToDump, :) = [];
  end
  % remove duplicates
  [~, L, ~] = unique(data_out.dt,'first');
  indexToDump = not(ismember(1:numel(data_out.dt), L));
  if sum(indexToDump) > 0
    data_out(indexToDump, :) = [];
  end
  % interpolate data on rounded datetime
  vars = data_in.Properties.VariableNames;
  vars(strcmp(vars, 'dt')) = [];
  for v = vars
    data_out.(v{:}) = interp1(data_in.dt, data_in.(v{:}), data_out.dt, 'linear');
  end
end


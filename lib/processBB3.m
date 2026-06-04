function [p, g, FiltStat] = processBB3(param, tot, filt_qc, filt_raw, filt_bad, di, ...
  TSG, di_method, filt_method, FTH, fth_constants, AC, CDOM, days2run)
  % Note DI is not interpolated as it's assume to be stable in time
  % BB3 parameters is a structure
  %   param.lambda <1x3 double> wavelength (nm)
  %   param.theta <1x1 double> scattering angle (degree)
  %   param.dark  <1x3 double> dark
  %   param.slope <1x3 double> slope
  
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
    if isempty(FTH.qc.tsw)
      error('No FLOW qc data loaded, when "filt_method == 25percentil", FLOW qc data level should be loaded')
    end
    if isempty(filt_qc)
      error('No BB3 qc data loaded, when "filt_method == 25percentil", BB3 qc data level should be loaded')
    end
    flow_data = FTH.qc.tsw;
  elseif strcmp(filt_method, 'exponential_fit')
    if isempty(FTH.raw.tsw)
      error('No FLOW raw data loaded, when "filt_method == exponential_fit", FLOW raw data level should be loaded')
    end
    if isempty(filt_raw)
      error('No BB3 raw data loaded, when "filt_method == exponential_fit", BB3 raw data level should be loaded')
    end
    flow_data = FTH.raw.tsw;
  end
  if islogical(flow_data.(FTH.view.swt_variable)); flow_data.(FTH.view.swt_variable) = double(flow_data.(FTH.view.swt_variable)); end
  flow_data.(FTH.view.swt_variable)(flow_data.(FTH.view.swt_variable) > 0) = 1;

  % round time stamp and remove time duplicates
  tot = round_timestamp(tot, minutes(1));
  filt_qc = round_timestamp(filt_qc, minutes(1));
  if ~isempty(filt_raw)
    filt_raw = round_timestamp(filt_raw);
  end
  if ~isempty(filt_bad)
    filt_bad = round_timestamp(filt_bad);
  end
  flow_data = round_timestamp(flow_data);
  
  % check if TSG data loaded
  tsg_data = [];
  if ~isempty(TSG)
    tsg_tbl_name = fieldnames(TSG.prod);
    if ~isempty(tsg_tbl_name)
      if ~isempty(TSG.prod.(tsg_tbl_name{1}))
        tsg_data = TSG.prod.(tsg_tbl_name{1});
      end
    end
    if isempty(tsg_data)
      if ~isempty(TSG.qc.tsw)
        tsg_data = TSG.qc.tsw;
      end
    end
    if ~isempty(tsg_data)
      if ~any(tsg_data.dt >= min([tot.dt; filt_qc.dt]) & tsg_data.dt <= max([tot.dt; filt_qc.dt]))
        warning(['TSG dates do not correspond to BB dates: T/S correction not applied before removing dissolved part from total.' newline ...
          'T/S assumed to be constant between filter and total events.n'])
        tsg_data = [];
      else
        % round time stamp and remove time duplicates
        tsg_data = round_timestamp(tsg_data, minutes(1));
      end
    end
  end
  if isempty(tsg_data)
    warning('No TSG qc or prod data loaded')
  end

  % check if AC data loaded
  if ~isempty(AC)
    ac_tbl_name = fieldnames(AC.prod);
    if isempty(ac_tbl_name)
      warning('No AC prod data loaded')
      ac_p = [];
      ac_g = [];
    else
      if any(strcmp(ac_tbl_name, 'p'))
        if isempty(AC.prod.p)
          warning('ap prod data not loaded')
          ac_p = [];
        else
          ac_p = AC.prod.p;
          if ~any(ac_p.dt >= min(tot.dt) & ac_p.dt <= max(tot.dt))
            warning('ap dates do not correspond to HBB dates')
            ac_p = [];
          else
            % round time stamp and remove time duplicates
            ac_p = round_timestamp(ac_p, minutes(1));
          end
        end
      end
      if any(strcmp(ac_tbl_name, 'g'))
        if isempty(AC.prod.g)
          warning('ag prod data not loaded')
          ac_g = [];
        else
          ac_g = AC.prod.g;
          if ~any(ac_g.dt >= min(tot.dt) & ac_g.dt <= max(tot.dt))
            warning('ag dates do not correspond to HBB dates')
            ac_g = [];
          else
            % round time stamp and remove time duplicates
            ac_g = round_timestamp(ac_g, minutes(1));
          end
        end
      end
    end
  else
    warning('No AC prod data loaded')
    ac_p = [];
    ac_g = [];
  end

  % if fdom_ag_parameters is input, check if FDOM data is loaded and interpolate on fth_interp
  fdom_data = [];
  if ~isempty(param.fdom_ag_parameters)
    if ~isempty(CDOM)
      fdom_tbl_name = fieldnames(CDOM.prod);
      if ~isempty(fdom_tbl_name)
        if ~isempty(CDOM.prod.(fdom_tbl_name{1}))
          fdom_data = CDOM.prod.(fdom_tbl_name{1});
        end
      end
      if ~isempty(fdom_data)
        if ~any(fdom_data.dt >= min([tot.dt; filt_qc.dt]) & fdom_data.dt <= max([tot.dt; filt_qc.dt]))
          warning('FDOM dates do not correspond to BB dates')
        else
          % round time stamp and remove time duplicates
          fdom_data = round_timestamp(fdom_data, minutes(1));
          % smooth fDOM if fDOM sensor is WSCD
          if strcmp(CDOM.model, 'WSCD') || strcmp(CDOM.prefix, 'WSCD')
            fdom_temp = table();
            foo_dt = [min([tot.dt; filt_qc.dt]) max([tot.dt; filt_qc.dt])];
            fdom_temp.dt = dateshift((foo_dt(1):minutes(1):foo_dt(2))','start','minute');
            id_nan = isnan(fdom_data.fdom);
            % interpolate fdom on fdom_temp
            fdom_temp.fdom = interp1(fdom_data.dt, fdom_data.fdom, fdom_temp.dt, 'linear');
            fdom_temp.fdom = fillmissing(fdom_temp.fdom, 'nearest');
            % smooth WSCD cdom data removing frequencies higher than 15min
            d = designfilt('lowpassiir', 'FilterOrder', 1, 'PassbandFrequency', 1/(60*15), 'SampleRate', 1/60);
            fdom_temp.fdom = filtfilt(d, fdom_temp.fdom);
            fdom_data.fdom = interp1(fdom_temp.dt, fdom_temp.fdom, fdom_data.dt, 'linear');
            fdom_data.fdom(id_nan) = NaN;
          end
        end
      end
      if ~isempty(param.fdom_ag_parameters) && ~isempty(fdom_data)
        replace_consecutive_nan = 30;
        fth_interp = merge_timeseries(fth_interp, fdom_data, 'fdom', '', replace_consecutive_nan);
      end
    else
      warning('No FDOM prod data loaded')
    end
  else
    warning('No FDOM prod data loaded')
  end

  % interpolate flow_data.swt onto binned data to fill missing flow data
  fth_interp = table([tot.dt; flow_data.dt; filt_qc.dt], 'VariableNames', {'dt'});
  fth_interp = round_timestamp(fth_interp, minutes(1));
  % sort dates
  fth_interp = sortrows(fth_interp, 'dt'); % sort dates
  fth_interp.swt = interp1(flow_data.dt, flow_data.(FTH.view.swt_variable), fth_interp.dt, 'previous');%, 'linear', 'extrap');
  fth_interp.swt = fth_interp.swt > 0;
  % Find switch events from total to filtered
  sel_start = find(fth_interp.swt(1:end-1) == SWITCH_TOTAL & fth_interp.swt(2:end) == SWITCH_FILTERED);
  % Find switch events from filtered to total
  sel_end = find(fth_interp.swt(1:end-1) == SWITCH_FILTERED & fth_interp.swt(2:end) == SWITCH_TOTAL);
  % Verify selections of filtered period
  if sel_start(1) > sel_end(1); sel_end(1) = []; end
  if sel_start(end) > sel_end(end); sel_end(end+1) = size(fth_interp.swt,1); end
  if size(sel_start,1) ~= size(sel_end,1); error('Inconsistent switch start/end data'); end

  % create filter average data table
  filt_avg = table(NaN(size(sel_start)), 'VariableNames', {'dt'});
  filt_avg.start = fth_interp.dt(sel_start);
  filt_avg.end = fth_interp.dt(sel_end);
  filt_avg.beta = NaN(size(filt_avg,1), size(param.lambda, 2));
  filt_avg.beta_avg_sd = NaN(size(filt_avg,1), size(param.lambda, 2));
  filt_avg.beta_avg_n = NaN(size(filt_avg,1), size(param.lambda, 2));
  % add flag for T/S correction (default true: data NOT T/S corrected)
  filt_avg.flag_deltaTSnotcorrected_filt = true(size(filt_avg,1), 1);

  % create filter interpolation table
  filt_interp = table(tot.dt, 'VariableNames', {'dt'});
  filt_interp.beta = NaN(size(tot.beta));
  filt_interp.beta_avg_sd = NaN(size(tot.beta));
  filt_interp.flag_deltaTSnotcorrected_total = true(size(tot.dt)); % by default set flag to true until data are corrected
  filt_interp.flag_deltaTSnotcorrected_filt0 = true(size(tot.dt)); % by default set flag to true until data are corrected
  filt_interp.flag_deltaTSnotcorrected_filt1 = true(size(tot.dt)); % by default set flag to true until data are corrected

  % add flag for T/S correction (default true: data NOT T/S corrected)
  filt_qc.flag_deltaTSnotcorrected_filt = true(size(filt_qc.dt));
  % add T/S to filt and filt_interp if TSG loaded
  if ~isempty(tsg_data)
    % interpolate linear T and S on filt_interp table
    filt_interp = addvars(filt_interp, interp_extrap(tsg_data, filt_interp.dt, TSG.temperature_variable, 30, true, 'linear', 'nearest'), ...
      interp_extrap(tsg_data, filt_interp.dt, TSG.salinity_variable, 30, true, 'linear', 'nearest'), 'NewVariableNames', {'t','s'}, 'After', 'dt');
    % interpolate T and S on filtered data
    filt_qc = addvars(filt_qc, interp_extrap(tsg_data, filt_qc.dt, TSG.temperature_variable, 30, true, 'linear', 'nearest'), ...
      interp_extrap(tsg_data, filt_qc.dt, TSG.salinity_variable, 30, true, 'linear', 'nearest'), 'NewVariableNames', {'t','s'}, 'After', 'dt');
    % add T and S variables on filt_avg table
    filt_avg = addvars(filt_avg, NaN(size(filt_avg,1),1), NaN(size(filt_avg,1),1), 'NewVariableNames', {'t','s'}, 'After', 'dt');
    % correct T/S in tot
    tot_ts_nan = isnan(filt_interp.t) | isnan(filt_interp.s);
    % Get beta seawater from Zhang et al. 2009
    tot_beta_sw = NaN(size(filt_interp.beta));
    for j = 1:size(filt_interp,1)
      tot_beta_sw(j, :) = betasw_ZHH2009(param.lambda, filt_interp.t(j), param.theta, filt_interp.s(j));
    end
    % substract beta seawater from total event
    tot.beta(~tot_ts_nan, :) = tot.beta(~tot_ts_nan, :) - tot_beta_sw(~tot_ts_nan, :);
    % flag when T/S not available in total events
    filt_interp.flag_deltaTSnotcorrected_total(~tot_ts_nan) = false; % change flag if corrected
  end

  % prepare fCDOM data if available
  if ~isempty(param.fdom_ag_parameters) && ~isempty(fdom_data)
    % interpolate spline and extrapolate nearest fdom on filt_interp table
    filt_interp = addvars(filt_interp, interp_extrap(fdom_data, filt_interp.dt, 'fdom', 30, false, 'linear'), ...
      'NewVariableNames', 'fdom', 'After', 'dt');
    % interpolate fdom on filtered data
    filt_qc = addvars(filt_qc, interp_extrap(fdom_data, filt_qc.dt, 'fdom', 30, false, 'linear'), ...
      'NewVariableNames', 'fdom', 'After', 'dt');
    % add fdom variables on filt_avg table
    filt_avg = addvars(filt_avg, NaN(size(filt_avg,1),1), 'NewVariableNames', 'fdom', 'After', 'dt');
  end

  switch filt_method
    case '25percentil'
      for i=1:size(sel_start, 1)
        sel_filt = fth_interp.dt(sel_start(i)) <= filt_qc.dt & filt_qc.dt <= fth_interp.dt(sel_end(i));
        if sum(sel_filt) > 0
          foo = filt_qc(sel_filt,:);
          if sum(sel_filt) == 1
            filt_avg.dt(i) = foo.dt;
            filt_avg.beta(i,:) = foo.beta;
            filt_avg.beta_avg_sd(i,:) = foo.beta_avg_sd;
            filt_avg.beta_avg_n(i,:) = foo.beta_avg_n;
            % repeat for T/S if tsg is loaded
            if ~isempty(tsg_data)
              filt_avg.t(i) = foo.t;
              filt_avg.s(i) = foo.s;
            end
          else
            perc25 = foo.beta > prctile(foo.beta, 25, 1);
            foo.beta_avg_sd(perc25) = NaN;
            foo.beta(perc25) = NaN;
            % compute average of all values smaller than 25th percentile for each filter event
            filt_avg.dt(i) = mean(foo.dt(any(~perc25, 2)), 'omitnan');
            filt_avg.beta(i,:) = mean(foo.beta, 1, 'omitnan');
            filt_avg.beta_avg_sd(i,:) = mean(foo.beta_avg_sd, 1, 'omitnan');
            filt_avg.beta_avg_n(i) = sum(foo.beta_avg_n(any(~isnan(foo.beta), 2)), 'omitnan');
            % repeat for T/S if tsg is loaded
            if ~isempty(tsg_data)
              % compute average of all values for each filter event
              filt_avg.t(i) = mean(foo.t, 1, 'omitnan');
              filt_avg.s(i) = mean(foo.s, 1, 'omitnan');
            end
          end
        end
      end
    case 'exponential_fit'
      % Based on method in: Dall’Olmo, G., Westberry, T.K., Behrenfeld, M.J., Boss, 
      %       E., Slade, W.H., 2009. Direct contribution of phytoplankton-sized particles 
      %       to optical backscattering in the open ocean. Biogeosciences Discuss 6, 291–340. 
      %       https://doi.org/10.5194/bgd-6-291-2009
      fprintf('Fitting exponential to filter events ... ')
      filt_avg.dt = mean([fth_interp.dt(sel_start) fth_interp.dt(sel_end)], 2);
      [filt_avg, FiltStat] = FiltExpFit('beta', filt_avg, filt_raw, filt_bad, fth_interp.dt(sel_start), fth_interp.dt(sel_end));
      fprintf('Done\n')      
      % run 25 percentile method on failed exponential fits
      for i=1:size(sel_start, 1)
        sel_filt = fth_interp.dt(sel_start(i)) <= filt_qc.dt & filt_qc.dt <= fth_interp.dt(sel_end(i));
        if sum(sel_filt) > 0
          foo = filt_qc(sel_filt,:);
          if sum(sel_filt) == 1
            filt_avg.dt(i) = foo.dt;
            % repeat for T/S if tsg is loaded
            if ~isempty(tsg_data)
              filt_avg.t(i) = foo.t;
              filt_avg.s(i) = foo.s;
            end
            if any(~FiltStat.exitflag(i,:))
              % filt_avg.beta(i,~FiltStat.exitflag(i,:)) = foo.beta(:,~FiltStat.exitflag(i,:));
              filt_avg.beta(i,:) = foo.beta;
            end
          else
            perc25 = foo.beta > prctile(foo.beta, 25, 1);
            foo.beta_avg_sd(perc25) = NaN;
            foo.beta(perc25) = NaN;
            % compute average of all values smaller than 25th percentile for each filter event
            filt_avg.dt(i) = mean(foo.dt(any(~perc25, 2)), 'omitnan');
            % repeat for T/S if tsg is loaded
            if ~isempty(tsg_data)
              % compute average of all values for each filter event
              filt_avg.t(i) = mean(foo.t, 1, 'omitnan');
              filt_avg.s(i) = mean(foo.s, 1, 'omitnan');
            end
            if any(~FiltStat.exitflag(i,:))
              % filt_avg.beta(i,~FiltStat.exitflag(i,:)) = mean(foo.beta(:,~FiltStat.exitflag(i,:)), 1, 'omitnan');
              % filt_avg.beta_avg_sd(i,~FiltStat.exitflag(i,:)) = mean(foo.beta_avg_sd(:,~FiltStat.exitflag(i,:)), 1, 'omitnan');
              filt_avg.beta(i,:) = mean(foo.beta, 1, 'omitnan');
              filt_avg.beta_avg_sd(i,:) = mean(foo.beta_avg_sd, 1, 'omitnan');
            end
            % filt_avg.beta_avg_n(i) = sum(foo.beta_avg_n(any(~isnan(foo.beta(:,~FiltStat.exitflag(i,:))), 2)), 'omitnan');
            filt_avg.beta_avg_n(i) = sum(foo.beta_avg_n(any(~isnan(foo.beta), 2)), 'omitnan');
          end
        else
          filt_avg.beta(i,:) = NaN;
        end
      end
    otherwise
      error('filter event method "filt_method" not supported')
  end
  filt_avg(all(isnan(filt_avg.beta), 2), :) = [];
  
  if ~isempty(tsg_data)
    % correct T/S in tot
    filt_ts_nan = isnan(filt_avg.t) | isnan(filt_avg.s);
    % Get beta seawater from Zhang et al. 2009
    filt_beta_sw = NaN(size(filt_avg.beta));
    for j = 1:size(filt_avg,1)
      filt_beta_sw(j, :) = betasw_ZHH2009(param.lambda, filt_avg.t(j), param.theta, filt_avg.s(j));
    end
    % substract beta seawater from total event
    filt_avg.beta(~filt_ts_nan, :) = filt_avg.beta(~filt_ts_nan, :) - filt_beta_sw(~filt_ts_nan, :);
    % flag when T/S not available in total events
    filt_avg.flag_deltaTSnotcorrected_filt(~filt_ts_nan) = false; % change flag if corrected
  end

  % Find switch events from filtered to total in fth_interp table
  total_events = table();
  total_events.start = NaN(size(filt_avg, 1), 1);
  total_events.end = NaN(size(filt_avg, 1), 1);
  for i = 1:size(total_events, 1)
    % find start/end of filter events in filt_interp
    if i == 1 && ~isempty(find(filt_interp.dt <= filt_avg.dt(i), 1, 'last'))
      total_events.start(i) = 1;
      total_events.end(i) = find(filt_interp.dt <= filt_avg.dt(i), 1, 'last');
    elseif i > 1
      total_events.start(i) = find(filt_interp.dt >= filt_avg.dt(i-1), 1, 'first');
      total_events.end(i) = find(filt_interp.dt <= filt_avg.dt(i), 1, 'last');
    end
  end
  if total_events.end(i) < size(filt_interp, 1)
    total_events = [total_events; table(NaN, NaN, 'VariableNames',{'start', 'end'})];
    total_events.start(i+1) = find(filt_interp.dt >= filt_avg.dt(i), 1, 'first');
    total_events.end(i+1) = size(filt_interp, 1);
  end
  total_events(all(isnan([total_events.start total_events.end]), 2), :) = [];

  % adjust T/S correction flags only for remaining filter events
  for i=1:size(total_events, 1) % sel_start
    % adjust flag for T/S correction
    if i == 1 && total_events.start(i) > 1 && ~filt_avg.flag_deltaTSnotcorrected_filt(i)
      filt_interp.flag_deltaTSnotcorrected_filt1(total_events.start(i):total_events.end(i)) = false;
    end
    if i > 1 && ~filt_avg.flag_deltaTSnotcorrected_filt(i-1)
      filt_interp.flag_deltaTSnotcorrected_filt0(total_events.start(i):total_events.end(i)) = false;
    end
    if i < size(total_events, 1) && ~filt_avg.flag_deltaTSnotcorrected_filt(i)
      filt_interp.flag_deltaTSnotcorrected_filt1(total_events.start(i):total_events.end(i)) = false;
    end
    if i == size(total_events, 1) && total_events.end(i) < size(filt_interp, 1)
      filt_interp.flag_deltaTSnotcorrected_filt0(total_events.start(i):total_events.end(i)) = false;
    end
  end

  % Interpolate filtered on total linearly
  filt_interp.beta = interp1(filt_avg.dt, filt_avg.beta, filt_interp.dt, 'linear');
  filt_interp.beta = fillmissing(filt_interp.beta, 'nearest', 'SamplePoints', filt_interp.dt);
  filt_interp.beta_avg_sd = interp1(filt_avg.dt, filt_avg.beta_avg_sd, filt_interp.dt, 'linear');
  filt_interp.beta_avg_sd = fillmissing(filt_interp.beta_avg_sd, 'nearest', 'SamplePoints', filt_interp.dt);
  
  % id only day to run in all tables to plot
  filt_interp_id = filt_interp.dt >= min(days2run) & filt_interp.dt < max(days2run)+days(1);
  tot_id = tot.dt >= min(days2run) & tot.dt < max(days2run)+days(1);
  filt_avg_id = filt_avg.dt >= min(days2run) & filt_avg.dt < max(days2run)+days(1);
  
  % plot
  if exist('visFlag', 'file') && exist('FTH', 'var')
    fh = visFlag([], filt_interp(filt_interp_id, :), tot(tot_id, :), [], filt_avg(filt_avg_id, :), [], 'beta', round(size(tot.beta, 2)/2), ...
      [], flow_data, FTH.view.spd_variable);
    title('Check filter event interpolation, press q to continue', 'FontSize', 14)
    legend('Filtered interpolated', 'Total', 'Filtered median', 'Flow rate',...
      'AutoUpdate','off', 'FontSize', 12)
    guiSelectOnTimeSeries(fh);
  elseif exist('visFlag', 'file')
    fh = visFlag([], filt_interp(filt_interp_id, :), tot(tot_id, :), [], filt_avg(filt_avg_id, :), [], 'beta', round(size(tot.beta, 2)/2), [], []);
    title('Check filter event interpolation, press q to continue', 'FontSize', 14)
    legend('Filtered interpolated', 'Total', 'Filtered median', 'Flow rate',...
      'AutoUpdate','off', 'FontSize', 12)
    guiSelectOnTimeSeries(fh);
  end
  
  % Compute beta particulate
  p = table(tot.dt, 'VariableNames', {'dt'});
  
  % Interpolate X_p with values from Sullivan et al. 2013
  theta_ref = 90:10:170;
  % From Zhang et al. (2021) taking into account angular dispersion:
  % Xiaodong Zhang, Edouard Leymarie, Emmanuel Boss, and Lianbo Hu, "Deriving the angular response function for backscattering sensors," Appl. Opt. 60, 8676-8687 (2021)
  X_p_ref = [0.709 0.884 1.024 1.121 1.174 1.171 1.138 1.039 0.931];
  % % From Sullivan et al. 2013:
  % X_p_ref = [0.684 0.858 1.000 1.097 1.153 1.167 1.156 1.131 1.093];
  X_p = interp1(theta_ref, X_p_ref, param.theta, 'spline');
  
  % attenuation correction and backscattering computation
  [p.betap, p.bbp, flags] = betap_Kcorrection(tot, filt_interp, param.lambda, param.k_exp, ...
    X_p, ac_p, ac_g, AC.lambda_a, AC.lambda_a, param.fdom_ag_parameters);
  
  % Calibrate beta_p (counts to m-1)
  p.betap = param.slope .* p.betap; % Dark independent
  p.bbp = param.slope .* p.bbp; % Dark independent

  % pass flags into flags table
  flags.deltaTSnotcorrected_filt0 = filt_interp.flag_deltaTSnotcorrected_filt0;
  flags.deltaTSnotcorrected_filt1 = filt_interp.flag_deltaTSnotcorrected_filt1;
  flags.deltaTSnotcorrected_total = filt_interp.flag_deltaTSnotcorrected_total;

  % Propagate error
  %   Note: Error is not propagated through Scattering & Residual temperature
  %         correction as required by SeaBASS
  p.betap_sd = param.slope .* sqrt(tot.beta_avg_sd.^2 + filt_interp.beta_avg_sd.^2);
  p.betap_n = tot.beta_avg_n;
  
  % Derive Gamma_bbp (does not support NaN values)
  % Correct bu on March 5, 2018, FitSpectra_HM2 does not accept NaNs in spectra
  flags(all(isnan(p.bbp),2), :) = [];
  p(all(isnan(p.bbp),2), :) = [];
  sel = all(~isnan(p.bbp), 2);
  fprintf('Computing Gamma_bbp ... ')
  [~, p.gamma_bbp(sel)] = FitSpectra_HM2(param.lambda, p.bbp(sel, :));
  p.gamma_bbp(p.gamma_bbp < -0.5) = NaN;
  fprintf('Done\n')
  
  % Estimate POC and Cphyto from bbp
  [p.poc, ~, ~, p.cphyto, ~, ~,] = estimatePOC_Cphyto(p.bbp, [470 532 650], 'soccom');
  
  % remove negative values
  flags(any(p.bbp < 0, 2), :) = [];
  p(any(p.bbp < 0, 2), :) = [];
  % set flag column
  p.flag_bit = set_flagbit(flags);
  % flag_info = read_flagbit(p.flag_bit, 'BB');

  if nargout > 1 && any(~isempty(di) | strcmp(di_method,'SW_scattering')) && ~isempty(tsg_data)
    % correct T/S in filt_qc
    filt_ts_nan = isnan(filt_qc.t) | isnan(filt_qc.s);
    % Get beta seawater from Zhang et al. 2009
    filt_beta_sw = NaN(size(filt_qc.beta));
    for j = 1:size(filt_qc,1)
      filt_beta_sw(j, :) = betasw_ZHH2009(param.lambda, filt_qc.t(j), param.theta, filt_qc.s(j));
    end
    % substract beta seawater from filter events
    filt_qc.beta(~filt_ts_nan, :) = filt_qc.beta(~filt_ts_nan, :) - filt_beta_sw(~filt_ts_nan, :);
    % flag when T/S not available in filt_qc events
    filt_qc.flag_deltaTSnotcorrected_filt(~filt_ts_nan) = false; % change flag if corrected
    
    % Correct T in DIW
    di.flag_deltaTnotcorrected_diw = true(size(di.dt));
    if ~isempty(di)
      diw_t = interp1(filt_qc.dt, filt_qc.t, di.dt, "linear");
      diw_t_nan = isnan(diw_t);
      % Get beta temperature during DIW from Zhang et al. 2009 and using water temperature reading from HyperBB
      beta_tdiw = NaN(size(filt_qc.beta));
      for j = 1:size(filt_qc,1)
        beta_tdiw(j, :) = betasw_ZHH2009(param.lambda, diw_t(j), param.theta, 0);
      end
      di.beta(~diw_t_nan, :) = di.beta(~diw_t_nan, :) - beta_tdiw(~diw_t_nan, :);
      % flag when T not available in DIW events
      di.flag_deltaTnotcorrected_diw(~diw_t_nan) = false;
    end

    % create beta dissolved table
    g = table(filt_qc.dt, 'VariableNames', {'dt'});
    switch di_method
      case 'interpolate'
        % Interpolate DI on Filtered
        %     + recommended if sensor drift with time
        di_pp = table(filt_avg.dt, 'VariableNames', {'dt'});
        di_pp.beta = interp1(di.dt, di.beta, di_pp.dt, "linear");
        di_pp.beta_avg_sd = interp1(di.dt, di.beta_avg_sd, di_pp.dt, "linear");
        % Compute beta dissolved
        g.betag = filt_qc.beta - di_pp.beta;
        % Propagate error
        %   Note: Error is not propagated through Zhang betasw_ZHH2009
        g.betag_sd = sqrt(filt_qc.beta_avg_sd.^2 + di_pp.beta_avg_sd.^2);
      case 'constant'
        % Average all given DI samples
        %     + recommended if no drift are observed with sensor
        di_pp = table(NaN, 'VariableNames', {'dt'});
        % Get not NaN DI values
        di_beta_sel = di.beta(all(~isnan(di.beta),2));
        di_beta_avg_sd_sel = di.beta_avg_sd(all(~isnan(di.beta),2));
        % Select only DI values within 5th and 75th percentile
        foo = prctile(di_beta_sel,[5, 75],1);
        avg_pl(1,:) = foo(1,:); % low percentile
        avg_ph(1,:) = foo(2,:); % high percentile
        avg_sel = any(avg_pl(1,:) <= di_beta_sel & di_beta_sel <= avg_ph(1,:),2);
        % Average values
        di_pp.beta = mean(di_beta_sel(avg_sel,:));
        di_pp.beta_avg_sd = mean(di_beta_avg_sd_sel(avg_sel,:));
        % Compute beta dissolved
        g.betag = filt_qc.beta;
        g.flag_deltaTSnotcorrected_filt = filt_qc.flag_deltaTSnotcorrected_filt;
        % Propagate error
        %   Note: Error is not propagated through Zhang betasw_ZHH2009
        g.betag_sd = filt_qc.beta_avg_sd;
      case 'SW_scattering'
      otherwise
        error('Method not supported.');
    end
    g.betag_n = filt_avg.beta_avg_n;
    
    % Compute bbg and gamma_bbg
    g.bbg = 2 * pi * X_p .* g.betag;
    
    % Derive Gamma_bbp (does not support NaN values)
    % Correct bu on March 5, 2018, FitSpectra_HM2 does not accept NaNs in spectra
    fprintf('Computing beta filtered slope... ')
    g(all(isnan(g.bbg),2),:) = [];
    sel = ~any(isnan(g.bbg), 2); % select wavelenght with no NaN
    [~, g.gamma_bbg(sel, :)] = FitSpectra_HM2(param.lambda, g.bbg(sel, :));

    % % set flag column TODO
    % g.flag_bit = set_flagbit(flag_g);
    % % flag_info = read_flagbit(g.flag_bit, 'ACS');

    fprintf('Done\n')
  else
    g = table();
  end
end

function [poc, poc_lower, poc_upper, cphyto, cphyto_lower, cphyto_upper] = estimatePOC_Cphyto(bbp, lambda, method)
  %ESTIMATE_POC_CPHYTO Particulate Organic Carbon (POC) and Cphyto are leanearly
  %   proportional to particulate backscattering bbp, various empirical relationship exist,
  %   few of them are implemented in this function.
  %
  % /!\ The calculations used are applicable only in the top layer
  %     with a maximum depth defined by max(MLD, Zeu).
  %
  %Syntax:  [ poc ] = estimate_poc( bbp, lambda, method, true )
  %Syntax:  [ poc, cphyto ] = estimate_poc( bbp, lambda, method, true )
  %
  %Inputs:
  %    Required:
  %        bbp NxM double corresponding to the values of the VSF at one angle in m^{-1}
  %    Optional:
  %        lambda 1x1 or 1xM double corresponding to the wavelength in nm
  %           default: 700
  %        method string of the name of the method to use for poc estimation
  %           default: 'soccom'
  %           soccom: POC = 3.23e4 x bbp(700) + 2.76
  %           an emprirical relationship built for the SOCCOM floats based on
  %           the relationship between the first profile of the floats and
  %           in-situ measurements taken during deployement
  %           (cruises: PS89, P16S and IN2015v1)
  %           NAB08_down or NAB08_up: Specific to North Atlantic in Spring
  %           based on empirical relationship (n=321), with data points
  %           ranging between 0-600 m, recommend downast
  %
  %Outputs:
  %     - poc NxM double corresponding to poc in mg.m^{-3}
  %     - poc_lower NxM double with corresponding to the lower poc estimation
  %     - poc_upper NxM double with corresponding to the upper poc estimation
  %     - cphyto NxM double corresponding to cphyto in mg.m^{-3}
  %     - cphyto_lower NxM double with corresponding to the lower cphyto estimation
  %     - cphyto_upper NxM double with corresponding to the upper cphyto estimation
  %
  %Examples:
  % [poc] = estimate_poc(bbp);
  % [poc] = estimate_poc(bbp, 700,'soccom', true);
  % [poc, cphyto] = estimate_poc(bbp);
  % [poc, cphyto] = estimate_poc(bbp, 700,'soccom', false);
  %
  %References:
  %   - Graff, J.R. et al., 2015. Analytical phytoplankton carbon measurements 
  %   spanning diverse ecosystems. Deep Sea Research Part I: Oceanographic 
  %   Research Papers 102, 16–25. https://doi.org/10.1016/j.dsr.2015.04.006
  
  %   - Cetinic I. et al., 2012. Particulate organic carbon and inherent optical
  %   properties during 2008 North Atlantic bloom experiment.
  %   J. Geophys. Res. Ocean. 117, doi:10.1029/2011JC007771.
  
  %   - Boss E. et al., 2013. The characteristics of particulate absorption,
  %   scattering and attenuation coefficients in the surface ocean;
  %   Contribution of the Tara Oceans expedition. Methods in Oceanography, 7:52?62
  %   ISSN 22111220. doi: 10.1016/j.mio.2013.11.002.
  %   URL http://dx.doi. org/10.1016/j.mio.2013.11.002.
  %
  % Tested with: Matlab R2015b & R2020b
  %
  % Author: Nils Haentjens, Ms, University of Maine
  % modified Guillaume Bourdin
  % Email: nils.haentjens@maine.edu / guillaume.bourdin@maine.edu
  % Created: February 5th 2016
  % modified: March 2020
  
  % Check Nargin
  if nargin > 4
     error('Too many input arguments.')
  elseif nargin < 1
     error('Not enough input arguments.')
  end
  % Set default param
  if ~exist('lambda','var')
    lambda = 700;
  end
  if ~exist('method','var') || isempty(method)
    method = 'soccom';
  end
  
  % Check size of input/content of input
  if size(bbp,2) ~= size(lambda,2) || size(lambda,1) ~= 1
    error('bbp should be NxM and lambda should be 1xM');
  end
  
  % Resize lambda
  lambda = bsxfun(@times, ones(size(bbp)), lambda);
  
  % Estimate poc
  switch method
    case 'soccom'
      % switch to bbp(700)
      bbp_700 = bbp .* (700 ./ lambda) .^ (-0.78);
      % estimate poc from bbp(700)
      poc = 3.23 * 10^4 * bbp_700 + 2.76;
      poc_lower = poc * 0.95;
      poc_upper = poc * 1.05;
    case 'NAB08_up'
      % upcast
      bbp_700 = bbp .* (700 ./ lambda) .^ (-0.78);
      poc = 43317 * bbp_700 - 18.4;
      poc_lower = (43317-2092) * bbp_700 - (18.4+5.8);
      poc_upper = (43317+2092) * bbp_700 - (18.4-5.8);
    case 'NAB08_down'
      % downcast
      bbp_700 = bbp .* (700 ./ lambda) .^ (-0.78);
      poc = 35422 * bbp_700 - 14.4; 
      poc_lower = (35422-1754) * bbp_700 - (14.4+5.8);
      poc_upper = (35422+1754) * bbp_700 - (14.4-5.8);
    otherwise
      error('Unknown method %s', method);
  end
  % Estimate cphyto
  % switch to bbp(470)
  bbp_470 = bbp .* (470 ./ lambda) .^ (-0.78);
  % estimate cphyto from bbp(470)
  cphyto = 12128 * bbp_470 + 0.59;
  cphyto_lower = cphyto * 0.95;
  cphyto_upper = cphyto * 1.05;
end


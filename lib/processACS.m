function [p, g, bad, regression_stats] = processACS(lambda, tot, filt, param, modelG50, modelmphi, di, ...
  CDOM, fth, fth_constants, interpolation_method, di_method, scattering_correction, compute_ad_aphi, tsg, days2run)
  % NOTE: wavelength of c are interpolated to wavelength of a
  % load psi_s & psi_t from Sullivan et al. 2006 for salinity & temprature correction
  psi = table();
  psi.wl = [400;402;404;406;408;410;412;414;416;418;420;422;424;426;428;430;432;434;436;438;440;442;444;446;448;450;452;454;456;458;460;462;464;466;468;470;472;474;476;478;480;482;484;486;488;490;492;494;496;498;500;502;504;506;508;510;512;514;516;518;520;522;524;526;528;530;532;534;536;538;540;542;544;546;548;550;552;554;556;558;560;562;564;566;568;570;572;574;576;578;580;582;584;586;588;590;592;594;596;598;600;602;604;606;608;610;612;614;616;618;620;622;624;626;628;630;632;634;636;638;640;642;644;646;648;650;652;654;656;658;660;662;664;666;668;670;672;674;676;678;680;682;684;686;688;690;692;694;696;698;700;702;704;706;708;710;712;714;716;718;720;722;724;726;728;730;732;734;736;738;740;742;744;746;748;750];
  psi.psiT = [0.0001;0.0001;0.0001;0.0001;0;0;0;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0002;0.0003;0.0003;0.0004;0.0005;0.0006;0.0006;0.0007;0.0008;0.0009;0.001;0.001;0.001;0.001;0.001;0.0009;0.0009;0.0008;0.0007;0.0006;0.0006;0.0005;0.0004;0.0003;0.0003;0.0002;0.0001;0.0001;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0002;0.0002;0.0002;0.0001;0.0001;0.0001;0;0;-0.0001;-0.0001;-0.0001;-0.0001;-0.0001;-0.0001;0;0;0.0001;0.0002;0.0003;0.0005;0.0007;0.0009;0.0013;0.0017;0.0021;0.0026;0.0032;0.0038;0.0045;0.0054;0.0063;0.0073;0.0083;0.0094;0.0104;0.0113;0.0121;0.0128;0.0133;0.0136;0.0136;0.0133;0.0129;0.0124;0.0116;0.0107];
  psi.a_psiS = [3.0e-05;3.0e-05;3.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;0;0;0;0;0;0;0;0;0;0;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;0;1.0e-05;2.0e-05;3.0e-05;3.0e-05;4.0e-05;5.0e-05;5.0e-05;6.0e-05;6.0e-05;6.0e-05;6.0e-05;6.0e-05;6.0e-05;6.0e-05;5.0e-05;5.0e-05;5.0e-05;5.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;3.0e-05;3.0e-05;3.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;1.0e-05;1.0e-05;1.0e-05;0;0;-1.0e-05;-2.0e-05;-3.0e-05;-4.0e-05;-6.0e-05;-7.0e-05;-8.0e-05;-9.0e-05;-0.00011;-0.00012;-0.00014;-0.00015;-0.00016;-0.00017;-0.00018;-0.00019;-0.00020;-0.00020;-0.00020;-0.00021;-0.00021;-0.00021;-0.00021;-0.00021;-0.00020;-0.00017;-0.00013;-8.0e-05;-1.0e-05;7.0e-05;0.00016;0.00026;0.00037;0.00046;0.00054;0.00061;0.00067];
  psi.c_psiS = [-1.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-4.0e-05;-3.0e-05;-3.0e-05;-2.0e-05;-1.0e-05;0;1.0e-05;1.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;0;0;0;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-3.0e-05;-4.0e-05;-5.0e-05;-6.0e-05;-6.0e-05;-8.0e-05;-9.0e-05;-0.00010;-0.00011;-0.00013;-0.00014;-0.00016;-0.00017;-0.00018;-0.00019;-0.00020;-0.00021;-0.00022;-0.00022;-0.00023;-0.00023;-0.00023;-0.00024;-0.00024;-0.00024;-0.00024;-0.00022;-0.00021;-0.00017;-0.00012;-6.0e-05;2.0e-05;0.00012;0.00022;0.00031;0.00041;0.00049;0.00056;0.00062];
  psi.sigma_psiT = [0.0002;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0;0;0;0;0.0001;0.0001;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0;0;0;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0002;0.0002;0.0003;0.0003;0.0004;0.0004;0.0004;0.0004;0.0005;0.0005;0.0006;0.0006;0.0007;0.0007;0.0007;0.0006;0.0005;0.0004;0.0003;0.0003;0.0004;0.0005;0.0006;0.0007;0.0008;0.0009];
  psi.c_sigma_psiS = [4e-005;4e-005;4e-005;4e-005;4e-005;4e-005;4e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;0;-1e-005;-2e-005;-3e-005;-4e-005;-6e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;2e-005;2e-005;2e-005;3e-005;3e-005;3e-005];
  psi.a_sigma_psiS = [3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;NaN;NaN;NaN;NaN;NaN;NaN;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;2e-005;2e-005;2e-005;3e-005];
  
  %% ap & cp
  % check FTH data
  if ~exist('fth_constants', 'var')
    % Assume most recent FlowControl software
    SWITCH_FILTERED = 1;
    SWITCH_TOTAL = 0;
  else
    SWITCH_FILTERED = fth_constants.SWITCH_FILTERED;
    SWITCH_TOTAL = fth_constants.SWITCH_TOTAL;
  end
  % check scattering correction method
  if ~any(strcmp(scattering_correction, {'Zaneveld1994_proportional', 'Rottgers2013_semiempirical', 'ZaneveldRottgers_blended'}))
    warning('only scattering correction supported: %s', strjoin({'Zaneveld1994_proportional', 'Rottgers2013_semiempirical', 'ZaneveldRottgers_blended'}, ', '))
    error('%s residual temperature and scattering correction not supported', scattering_correction)
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
        % smooth fDOM if fDOM sensor is WSCD
        if strcmp(CDOM.model, 'WSCD') || strcmp(CDOM.prefix, 'WSCD')
          ancillary = table();
          foo_dt = [datetime(min([tot.dt; filt.dt]), 'ConvertFrom', 'datenum') datetime(max([tot.dt; filt.dt]), 'ConvertFrom', 'datenum')];
          ancillary.dt = datenum(dateshift((foo_dt(1):minutes(1):foo_dt(2))','start','minute'));
          id_nan = isnan(cdom_base.fdom);
          % create new fdom table interpolating for gaps < 60 min
          ancillary.fdom = interp1(cdom_base.dt, cdom_base.fdom, ancillary.dt, 'linear');
          ancillary.fdom = fillmissing(ancillary.fdom, 'nearest');
          % smooth WSCD cdom data removing frequencies higher than 15min
          d = designfilt('lowpassiir', 'FilterOrder', 1, 'PassbandFrequency', 1/(60*15), 'SampleRate', 1/60);
          ancillary.fdom = filtfilt(d, ancillary.fdom);
          cdom_base.fdom = interp1(ancillary.dt, ancillary.fdom, cdom_base.dt, 'linear');
          cdom_base.fdom(id_nan) = NaN;
        end
      end
      if ~any(cdom_base.dt >= min([tot.dt; filt.dt]) & cdom_base.dt <= max([tot.dt; filt.dt]))
        warning('fCDOM dates do not correspond to ACS dates: interpolation switched to "linear"')
        interpolation_method = 'linear';
      else
        % round time stamp and remove time duplicates
        cdom_base = round_timestamp(cdom_base);
      end
    else
      warning('fCDOM data not loaded: interpolation switched to "linear"')
      interpolation_method = 'linear';
    end
    if ~isfield(param, 'min_nb_pts_per_cluster')
      warning('Field "min_nb_pts_per_cluster" missing: set to default 100')
      param.min_nb_pts_per_cluster = 100;
    end
    if ~isfield(param, 'time_weight_for_cluster')
      warning('Field "time_weight_for_cluster" missing: set to default 0.9')
      param.time_weight_for_cluster = 0.9;
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
  fth_interp = round_timestamp(fth_interp);
  fth_interp = sortrows(fth_interp, 'dt'); % sort dates
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
  filt_avg.a = NaN(size(filt_avg,1), size(lambda.a, 2));
  filt_avg.c = NaN(size(filt_avg,1), size(lambda.c, 2));
  filt_avg.a_avg_sd = NaN(size(filt_avg,1), size(lambda.a, 2));
  filt_avg.c_avg_sd = NaN(size(filt_avg,1), size(lambda.c, 2));
  filt_avg.a_avg_n = NaN(size(filt_avg,1), 1);
  filt_avg.c_avg_n = NaN(size(filt_avg,1), 1);

  % create filter interpolation table
  filt_interp = table(tot.dt, 'VariableNames', {'dt'});
  filt_interp.a = NaN(size(tot.a));
  filt_interp.c = NaN(size(tot.c));

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
    % remove psiS and psiT from tot and filt with Tref = median(T) of whole processed section
    Tref = median([filt.t; filt_interp.t], 'omitnan');
    a_psiS = interp1(psi.wl, psi.a_psiS, lambda.a, 'spline');
    c_psiS = interp1(psi.wl, psi.c_psiS, lambda.c, 'spline');
    a_psiT = interp1(psi.wl, psi.psiT, lambda.a, 'spline');
    c_psiT = interp1(psi.wl, psi.psiT, lambda.c, 'spline');
    % Correct independently only when t and s are available to avoid getting NaN when t or s are NaN
    %%%%%%%%%%%% TODO add flag when filt or tot are not corrected (t or s not available) %%%%%%%%%%%%
    filt_t_nan = isnan(filt.t);
    filt_s_nan = isnan(filt.s);
    tot_t_nan = isnan(filt_interp.t);
    tot_s_nan = isnan(filt_interp.s);
    filt.a(~filt_t_nan, :) = filt.a(~filt_t_nan, :) - (a_psiT .* (filt.t(~filt_t_nan) - Tref));
    filt.c(~filt_t_nan, :) = filt.c(~filt_t_nan, :) - (c_psiT .* (filt.t(~filt_t_nan) - Tref));
    filt.a(~filt_s_nan, :) = filt.a(~filt_s_nan, :) - (a_psiS .* filt.s(~filt_s_nan));
    filt.c(~filt_s_nan, :) = filt.c(~filt_s_nan, :) - (c_psiS .* filt.s(~filt_s_nan));
    % filt.a = filt.a - (a_psiS .* filt.s + a_psiT .* (filt.t - Tref)); % old code
    % filt.c = filt.c - (c_psiS .* filt.s + c_psiT .* (filt.t - Tref)); % old code
    tot.a(~tot_t_nan, :) = tot.a(~tot_t_nan, :) - (a_psiT .* (filt_interp.t(~tot_t_nan) - Tref));
    tot.c(~tot_t_nan, :) = tot.c(~tot_t_nan, :) - (c_psiT .* (filt_interp.t(~tot_t_nan) - Tref));
    tot.a(~tot_s_nan, :) = tot.a(~tot_s_nan, :) - (a_psiS .* filt_interp.s(~tot_s_nan));
    tot.c(~tot_s_nan, :) = tot.c(~tot_s_nan, :) - (c_psiS .* filt_interp.s(~tot_s_nan));
    % tot.a = tot.a - (a_psiS .* filt_interp.s + a_psiT .* (filt_interp.t - Tref)); % old code
    % tot.c = tot.c - (c_psiS .* filt_interp.s + c_psiT .* (filt_interp.t - Tref)); % old code
  end

  % prepare fCDOM data if CDOM interpolation
  if strcmp(interpolation_method, 'CDOM')
    % interpolate spline and extrapolate nearest fdom on filt_interp table
    filt_interp = addvars(filt_interp, interp_extrap(cdom_base, filt_interp.dt, 'fdom', 30, false, 'linear'), ...
      'NewVariableNames', 'fdom', 'After', 'dt');
    % interpolate fdom on filtered data
    filt = addvars(filt, interp_extrap(cdom_base, filt.dt, 'fdom', 30, false, 'linear'), ...
      'NewVariableNames', 'fdom', 'After', 'dt');
    % add fdom variables on filt_avg table
    filt_avg = addvars(filt_avg, NaN(size(filt_avg,1),1), 'NewVariableNames', 'fdom', 'After', 'dt');
  end
  
  for i=1:size(sel_start, 1)
    sel_filt = fth_interp.dt(sel_start(i)) <= filt.dt & filt.dt <= fth_interp.dt(sel_end(i));
    foo = filt(sel_filt,:);
    % filter average of a and c
    if sum(sel_filt) == 1
      filt_avg.dt(i) = foo.dt;
      filt_avg.start(i) = foo.dt;
      filt_avg.end(i) = foo.dt;
      filt_avg.a(i,:) = foo.a;
      filt_avg.c(i,:) = foo.c;
      filt_avg.a_avg_sd(i,:) = foo.a_avg_sd;
      filt_avg.c_avg_sd(i,:) = foo.c_avg_sd;
      filt_avg.a_avg_n(i) = foo.a_avg_n;
      filt_avg.c_avg_n(i) = foo.c_avg_n;
      % repeat for T/S if tsg is loaded
      if tsg_loaded
        filt_avg.t(i,:) = foo.t;
        filt_avg.s(i,:) = foo.s;
      end
      % repeat for fdom if interpolation method == CDOM
      if strcmp(interpolation_method, 'CDOM')
        filt_avg.fdom(i,:) = foo.fdom;
      end
    elseif sum(sel_filt) > 1
      a_perc25 = foo.a > prctile(foo.a, 25, 1);
      c_perc25 = foo.c > prctile(foo.c, 25, 1);
      foo.a_avg_sd(a_perc25) = NaN;
      foo.c_avg_sd(c_perc25) = NaN;
      foo.a(a_perc25) = NaN;
      foo.c(c_perc25) = NaN;
      % compute average of all values smaller than 25th percentile for each filter event
      filt_avg.dt(i) = mean(foo.dt(any(~a_perc25, 2) | any(~c_perc25, 2)), 'omitnan');
      if any(any(any(~a_perc25, 2) | any(~c_perc25, 2), 2))
        filt_avg.start(i) = min(foo.dt(any(~a_perc25, 2) | any(~c_perc25, 2)));
        filt_avg.end(i) = max(foo.dt(any(~a_perc25, 2) | any(~c_perc25, 2)));
      else
        filt_avg.start(i) = NaN;
        filt_avg.end(i) = NaN;
      end
      filt_avg.a(i,:) = mean(foo.a, 1, 'omitnan');
      filt_avg.c(i,:) = mean(foo.c, 1, 'omitnan');
      filt_avg.a_avg_sd(i,:) = mean(foo.a_avg_sd, 1, 'omitnan');
      filt_avg.c_avg_sd(i,:) = mean(foo.c_avg_sd, 1, 'omitnan');
      filt_avg.a_avg_n(i) = sum(foo.a_avg_n(any(~isnan(foo.a), 2)), 'omitnan');
      filt_avg.c_avg_n(i) = sum(foo.c_avg_n(any(~isnan(foo.c), 2)), 'omitnan');
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
  rm_filter_event = all(isnan(filt_avg.a), 2) | all(isnan(filt_avg.c), 2);
  % sel_start(rm_filter_event) = [];
  % sel_end(rm_filter_event) = [];
  filt_avg(rm_filter_event, :) = [];
  
  switch interpolation_method
    case 'CDOM'
      % interpolate ag and cg using fcdom
      filt_interp = agcg_fdom_interpolation(tot, filt_interp, filt, lambda, CDOM.dark, filt_avg, param.min_nb_pts_per_cluster, param.time_weight_for_cluster);
  
      % TODO: propagate error from regression and filter event STD
      filt_interp.a_avg_sd = interp1(filt.dt, filt.a_avg_sd, filt_interp.dt, 'linear');
      filt_interp.c_avg_sd = interp1(filt.dt, filt.c_avg_sd, filt_interp.dt, 'linear');
      
      % remove interpolation when there is no data
      filt_interp.a_avg_sd(all(isnan(tot.a),2), :) = NaN;
      filt_interp.a(all(isnan(tot.a),2), :) = NaN;
      filt_interp.c_avg_sd(all(isnan(tot.c),2), :) = NaN;
      filt_interp.c(all(isnan(tot.c),2), :) = NaN;
      
      % % check interpolation if bad replace by linear interpolation TODO add
      % % this in agcg_fdom_interpolation and apply it by section instead of by point
      % test_interp = table();
      % test_interp.ap = tot.a - filt_interp.a;
      % test_interp.cp = tot.c - filt_interp.c;
      % id_bad_interp_a = sum(test_interp.ap < 0, 2) > size(lambda.a, 2)/3;
      % id_bad_interp_c = sum(test_interp.cp < 0, 2) > size(lambda.c, 2)/3;
      % filt_interp.a(id_bad_interp_a, :) = filt_interp.lin_a(id_bad_interp_a, :);
      % filt_interp.c(id_bad_interp_c, :) = filt_interp.lin_c(id_bad_interp_c, :);
      
    case 'linear'
      % Interpolate filtered on total linearly
      filt_interp.a = interp1(filt_avg.dt, filt_avg.a, filt_interp.dt, 'linear');
      filt_interp.a = fillmissing(filt_interp.a, 'nearest', 'SamplePoints', filt_interp.dt);
      filt_interp.c = interp1(filt_avg.dt, filt_avg.c, filt_interp.dt, 'linear');
      filt_interp.c = fillmissing(filt_interp.c, 'nearest', 'SamplePoints', filt_interp.dt);
      filt_interp.a_avg_sd = interp1(filt_avg.dt, filt_avg.a_avg_sd, filt_interp.dt, 'linear');
      filt_interp.a_avg_sd = fillmissing(filt_interp.a_avg_sd, 'nearest', 'SamplePoints', filt_interp.dt);
      filt_interp.c_avg_sd = interp1(filt_avg.dt, filt_avg.c_avg_sd, filt_interp.dt, 'linear');
      filt_interp.c_avg_sd = fillmissing(filt_interp.c_avg_sd, 'nearest', 'SamplePoints', filt_interp.dt);
      regression_stats = struct();
    otherwise
    error('Method not supported.');
  end

  % Remove lines full of NaNs or with inf data
  sel2rm = any(~isfinite(tot.a),2) | any(~isfinite(tot.c),2)| all(isnan(tot.a),2) | ...
           all(isnan(tot.c),2) | all(isnan(filt_interp.a),2) | all(isnan(filt_interp.c),2);
  tot(sel2rm,:) = [];
  filt_interp(sel2rm,:) = [];
  
  if exist('visFlag', 'file')
    % id only day to run in all tables to plot
    filt_interp_id = filt_interp.dt >= min(days2run) & filt_interp.dt < max(days2run)+1;
    tot_id = tot.dt >= min(days2run) & tot.dt < max(days2run)+1;
    filt_avg_id = filt_avg.dt >= min(days2run) & filt_avg.dt < max(days2run)+1;
    % plot
    fh = visFlag([], filt_interp(filt_interp_id, :), tot(tot_id, :), [], filt_avg(filt_avg_id, :), [], 'a', round(size(tot.a, 2)/2), [], []);
    title('Check filter event interpolation, press q to continue', 'FontSize', 14)
    if strcmp(interpolation_method, 'CDOM')
      ax1 = gca;
      ax1.YColor = [0	205	205]/255;
      scatter(filt_interp.dt(filt_interp_id), filt_interp.fdom(filt_interp_id), '.', 'MarkerEdgeColor', [0	205	205]/255, 'MarkerEdgeAlpha', 0.5)
      plot(filt_avg.dt(filt_avg_id), filt_avg.fdom(filt_avg_id), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0	205	205]/255)
      ylabel('FDOM (volts)')
      legend('Filtered interpolated', 'Total', 'Filtered (avg filter event)', 'fCDOM (raw)', ...
        'fCDOM (avg filter event)', 'AutoUpdate','off', 'FontSize', 12)
    else
      legend('Filtered interpolated', 'Total', 'Filtered median',...
        'AutoUpdate','off', 'FontSize', 12)
    end
    guiSelectOnTimeSeries(fh);
  end
  
  % figure();
  % subplot(1, 3, 1)
  % scatter(filt_interp.a(filt_interp.flag_linear_interp,40), filt_interp.s(filt_interp.flag_linear_interp), 20, 'k+')
  % scatter(filt_interp.a(~filt_interp.flag_linear_interp,40), filt_interp.s(~filt_interp.flag_linear_interp), 20, 'bo', 'filled')
  % ylabel('Salinity')
  % xlabel('ag')
  % subplot(1, 3, 2)
  % scatter(filt_interp.c(filt_interp.flag_linear_interp,40), filt_interp.s(filt_interp.flag_linear_interp), 20, 'k+')
  % scatter(filt_interp.c(~filt_interp.flag_linear_interp,40), filt_interp.s(~filt_interp.flag_linear_interp), 20, 'bo', 'filled')
  % ylabel('Salinity')
  % xlabel('cg')
  % subplot(1, 3, 3)
  % scatter(filt_interp.fdom(filt_interp.flag_linear_interp), filt_interp.s(filt_interp.flag_linear_interp), 20, 'k+')
  % scatter(filt_interp.fdom(~filt_interp.flag_linear_interp), filt_interp.s(~filt_interp.flag_linear_interp), 20, 'bo', 'filled')
  % ylabel('Salinity')
  % xlabel('fdom')

  % Particulate = Total - FSW
  p = table(tot.dt, 'VariableNames', {'dt'});
  p.ap = tot.a - filt_interp.a;
  p.cp = tot.c - filt_interp.c;
  
  if size(lambda.a, 2) > 50 % perform two separate corrections for ap and cp only for ACS data, not AC9
    % Interpolate wavelengths for Scattering & Residual temperature correction
    ap_for_cpresiduals_corr = interp1(lambda.a', p.ap', lambda.c', 'linear', 'extrap')';
    cp_for_apresiduals_corr = interp1(lambda.c', p.cp', lambda.a', 'linear', 'extrap')';
    % ap Scattering & Residual temperature correction
    % cp Residual correction (for efficiency use the one computed from ap as it should be the same)
    fprintf('ap %s residual temperature and scattering correction ... ', strrep(scattering_correction, '_', ' '))
    if strcmp(scattering_correction, 'Rottgers2013_semiempirical')
      [p.ap, ~] = ResidualTempScatterCorrRottgers_semiempirical(p.ap, cp_for_apresiduals_corr, lambda.a, psi);
    elseif strcmp(scattering_correction, 'Zaneveld1994_proportional')
      [p.ap, ~] = ResidualTempScatterCorrZaneveld_proportional(p.ap, cp_for_apresiduals_corr, lambda.a, psi);
    elseif strcmp(scattering_correction, 'ZaneveldRottgers_blended')
      [p.ap, ~] = ResidualTempScatterCorrZaneveldRottgers_blended(p.ap, cp_for_apresiduals_corr, lambda.a, psi);
    end
    fprintf('Done\n')
    fprintf('cp %s residual temperature and scattering correction ... ', strrep(scattering_correction, '_', ' '))
    if strcmp(scattering_correction, 'Rottgers2013_semiempirical')
      [~, p.cp] = ResidualTempScatterCorrRottgers_semiempirical(ap_for_cpresiduals_corr, p.cp, lambda.c, psi);
    elseif strcmp(scattering_correction, 'Zaneveld1994_proportional')
      [~, p.cp] = ResidualTempScatterCorrZaneveld_proportional(ap_for_cpresiduals_corr, p.cp, lambda.c, psi);
    elseif strcmp(scattering_correction, 'ZaneveldRottgers_blended')
      [~, p.cp] = ResidualTempScatterCorrZaneveldRottgers_blended(ap_for_cpresiduals_corr, p.cp, lambda.c, psi);
    end
    fprintf('Done\n')
  else
    fprintf('ap and cp %s residual temperature and scattering correction ... ', strrep(scattering_correction, '_', ' '))
    if strcmp(scattering_correction, 'Rottgers2013_semiempirical')
      [p.ap, p.cp] = ResidualTempScatterCorrRottgers_semiempirical(p.ap, p.cp, lambda.ref, psi);
    elseif strcmp(scattering_correction, 'Zaneveld1994_proportional')
      [p.ap, p.cp] = ResidualTempScatterCorrZaneveld_proportional(p.ap, p.cp, lambda.ref, psi);
    elseif strcmp(scattering_correction, 'ZaneveldRottgers_blended')
      [p.ap, p.cp] = ResidualTempScatterCorrZaneveldRottgers_blended(p.ap, p.cp, lambda.ref, psi);
    end
  end
  
  % Remove lines full of NaNs (Rottgers2013_semiempirical potentially fail when temperature correct is too large)
  sel2rm = all(isnan(p.ap), 2) | all(isnan(p.cp), 2);
  p(sel2rm, :) = [];
  tot(sel2rm, :) = [];
  filt_interp(sel2rm, :) = [];

  % Propagate error (using geometric mean of measurement errors)
  %   Note: Error is not propagated through Scattering & Residual temperature
  %         correction as required by SeaBASS
  p.ap_sd = sqrt(tot.a_avg_sd.^2 + filt_interp.a_avg_sd.^2);
  p.cp_sd = sqrt(tot.c_avg_sd.^2 + filt_interp.c_avg_sd.^2);
  p.ap_n = tot.a_avg_n;
  p.cp_n = tot.c_avg_n;
  
  % % delete ap spectra full of NaNs
  % p(all(isnan(p.ap),2),:) = [];
  % p(all(isnan(p.cp),2),:) = [];
  
  % Unsmoothing ACS spectra
  % Ron Zaneveld, WET Labs, Inc., 2005
  if size(lambda.a, 2) > 50 % unsmooth only ACS, not AC9
    p = unsmoothACS(p, lambda);
  end
  
  % Auto QC on product spectra
  wla_430 = lambda.a(find(lambda.a <= 430, 1,'last')); % find lower and closest to 430nm wavelength
  wla_700 = lambda.a(find(lambda.a >= 700, 1,'first')); % find higher and closest to 700nm wavelength
  
  % replace values at each end of the spectra when < -0.0015
  p.ap_sd(p.ap < -0.0015 & lambda.a < wla_430) = NaN;
  % p.ap_sd(p.ap < -0.0015 & lambda.a >= wla_700) = NaN;
  p.ap(p.ap < -0.0015 & lambda.a < wla_430) = NaN;
  % p.ap(p.ap < -0.0015 & lambda.a >= wla_700) = NaN;
  
  % set flag matrix
  flag_varname = {...
    'fCDOM_mix_a_cluster', 'fCDOM_mix_c_cluster', 'fCDOM_a_cluster_chg', ...
    'fCDOM_c_cluster_chg', 'flag_linear_interp', 'cp_neg','ap_neg', ...
    'ap_shape','ap_bubbles','cp_bubbles','ap430_700_neg','cp_over10','noisy600_650',...
    'ap460_640_04_450','positive_ap450_570','poc_flag','chl_ap676lh_flag',...
    'gamma_flag','chl_Halh_flag','HH_mphi_flag','HH_G50_flag','chlratio_flag',...
    'gamma_suspicious','poc_suspicious','chl_ap676lh_suspicious',...
    'chl_Halh_suspicious','HH_G50_mphi_suspicious'};
  flag = array2table(false(size(p,1), length(flag_varname)), 'VariableNames', flag_varname);
  if strcmp(interpolation_method, 'CDOM')
    flag.fCDOM_mix_a_cluster = filt_interp.fCDOM_mix_a_cluster;
    flag.fCDOM_mix_c_cluster = filt_interp.fCDOM_mix_c_cluster;
    flag.fCDOM_a_cluster_chg = filt_interp.fCDOM_a_cluster_chg;
    flag.fCDOM_c_cluster_chg = filt_interp.fCDOM_c_cluster_chg;
    flag.flag_linear_interp = filt_interp.flag_linear_interp;
  else
    flag.flag_linear_interp = true(size(p,1), 1);
  end

  % delete cp spectra when any cp < -0.0015
  todelete = any(p.cp < -0.0015, 2);
  if sum(todelete) > 0
    fprintf('%.2f%% (%i) spectra failed auto-QC: cp < -0.0015\n', ...
      sum(todelete) / size(p, 1) * 100, sum(todelete))
  end
  flag.cp_neg(todelete) = true;
  bad = [p(todelete, :) table(repmat({'cp < -0.0015'}, ...
    sum(todelete), 1), 'VariableNames', {'QC_failed'})];
  p.cp(todelete, :) = NaN;

  % flag ap spectra when ap430-700 < -0.0015
  toflag = any(p.ap < -0.0015 & lambda.a >= wla_430 & lambda.a <= wla_700, 2);
  if sum(toflag) > 0
    fprintf('%.2f%% (%i) spectra flagged: ap 430-700 < -0.0015\n', ...
      sum(toflag) / size(p, 1) * 100, sum(toflag))
  end
  flag.ap430_700_neg(toflag) = true;
  % bad = [bad; p(toflag, :) table(repmat({'ap 430-700 < -0.0015'}, ...
  %   sum(toflag), 1), 'VariableNames', {'QC_failed'})];
  % p.ap(toflag, :) = NaN;
  
  % flag ap spectra when any ap < -0.01
  toflag = any(p.ap < -0.01 & lambda.a >= wla_430 & lambda.a <= wla_700, 2);
  if sum(toflag) > 0
    fprintf('%.2f%% (%i) spectra flagged: ap 430-700 < -0.01\n', ...
      sum(toflag) / size(p, 1) * 100, sum(toflag))
  end
  flag.ap_neg(toflag) = true;
  % bad = [bad; p(todelete, :) table(repmat({'ap 430-700 < -0.01'}, ...
  %   sum(todelete), 1), 'VariableNames', {'QC_failed'})];
  % p.ap(todelete, :) = NaN;
  
  % flag when ap650 > ap676: high NAP
  toflag = any(mean(p.ap(:, lambda.a >= 640 & lambda.a <= 655), 2, 'omitnan') > ...
    mean(p.ap(:, lambda.a >= 670 & lambda.a <= 680), 2, 'omitnan'), 2);
  if sum(toflag) > 0
    fprintf('%.2f%% (%i) spectra flagged: ap650 > ap676 high NAP\n', ...
      sum(toflag) / size(p, 1) * 100, sum(toflag))
  end
  flag.ap_shape(toflag) = true;
  % bad = [bad; p(toflag, :) table(repmat({'ap650 > ap676: high NAP'}, ...
  %   sum(toflag), 1), 'VariableNames', {'QC_failed'})];
  % p.cp(toflag, :) = NaN;
  
  % flag attenuation spectra when cp > 10
  toflag = any(p.cp > 10, 2);
  if sum(toflag) > 0
    fprintf('%.2f%% (%i) spectra flagged: p.cp > 10\n', ...
      sum(toflag) / size(p, 1) * 100, sum(toflag))
  end
  flag.cp_over10(toflag) = true;
  % bad = [bad; p(toflag, :) table(repmat({'p.cp > 10'}, ...
  %   sum(toflag), 1), 'VariableNames', {'QC_failed'})];
  % p.cp(toflag, :) = NaN;
  
  % find wavelength below and above which ap and cp are unrealistic and replace by NaNs
  if size(lambda.a, 2) > 50 % clean only ACS data, not AC9
    % replace unrealistic ap and cp red wavelenght by NaN
    [~, foo] = min(p.cp(:, lambda.c > 715),[],2);
    minred_wlc = sum(lambda.c <= 715) + foo;
    [~, foo] = min(abs(p.ap(:, lambda.a > 715)),[],2);
    minred_wla = sum(lambda.a <= 715) + foo;
    p.cp(lambda.c > lambda.c(minred_wlc)' & lambda.c > 720) = NaN;
    p.cp_sd(lambda.c > lambda.c(minred_wlc)' & lambda.c > 720) = NaN;
    p.ap(lambda.a > lambda.a(minred_wla)' & lambda.a > 715) = NaN;
    p.ap_sd(lambda.a > lambda.a(minred_wla)' & lambda.a > 715) = NaN;
    
    % replace unrealistic ap in blue wavelength by NaN
    blue_wl = p.ap;
    blue_wl(:,lambda.a > 550) = NaN;
    blue_wl_var = [zeros(size(blue_wl, 1), 1) abs(diff(blue_wl,[],2))]; % get absolute derivative over wavelengths  
    cutoffblue_wla = NaN(size(blue_wl_var, 1), 1);
    for i = 1:size(cutoffblue_wla,1)
      foo = lambda.a(find(blue_wl_var(i,:) > 6 * mean(blue_wl_var(i, lambda.a > 450), 2, 'omitnan'), 1, 'last'));
      if ~isempty(foo)
        cutoffblue_wla(i) = foo;
      else
        cutoffblue_wla(i) = min(lambda.a);
      end
    end
    if ~all(isnan(cutoffblue_wla))
      p.ap_sd(lambda.a < cutoffblue_wla & lambda.a <= 450) = NaN;
      p.ap(lambda.a < cutoffblue_wla & lambda.a <= 450) = NaN;
    end
    
    % Auto QC from ratio std 600-650 / 640-676
    ap_450 = p.ap(:, find(lambda.a >= 450, 1,'first'));
    ratiod600_ap450 = sum(abs(diff(p.ap(:, lambda.a > 600 & lambda.a <= 650), 2, 2)),2) ./ ap_450;
    
    % get automatic threshold for ratio std 600-650 / 640-676
    fudge_list = (0.1:0.1:10)';
    ndel_spec = NaN(size(fudge_list));
    for i=1:size(fudge_list,1)
      toflag = ratiod600_ap450 > fudge_list(i) * median(ratiod600_ap450, 'omitnan');
      ndel_spec(i) = sum(toflag);
    end
    fudge_factor = fudge_list(find(abs(diff(ndel_spec)) == min(abs(diff(ndel_spec))) & ...
      ndel_spec(2:end) < 0.001 * max(ndel_spec), 1,  'first')); % 0.05 threshold on first derivative of number of spectra deleted
    if isempty(fudge_factor)
      fudge_factor = fudge_list(find(abs(diff(ndel_spec)) == min(abs(diff(ndel_spec))) & ...
        ndel_spec(2:end) < 0.01 * max(ndel_spec), 1,  'first')); % 0.05 threshold on first derivative of number of spectra deleted
    end
    if ~isempty(fudge_factor)
      toflag = ratiod600_ap450 > fudge_factor * median(ratiod600_ap450);
      if sum(toflag) > 0
        fprintf('%.2f%% (%i) spectra flagged: sum(abs(d(ap)/d(lambda(600-650)))) / ap450nm\n', ...
          sum(toflag) / size(p, 1) * 100, sum(toflag))
      end
  
    %   delete bad spectra
      flag.noisy600_650(toflag) = true;
      bad = [bad; p(toflag, :) ...
        table(repmat({'sum(abs(d(ap)/d(lambda(600-650)))) / ap_{450nm}'}, sum(toflag), 1), ...
        'VariableNames', {'QC_failed'})];
  %     p.ap(toflag, :) = NaN;
    end
    % flag potential bubbles using the step in the center of the ap spectra
    foo_difstep = diff(p.ap(:, lambda.a > 560 & lambda.a <= 600),1, 2) ./ diff(lambda.a(lambda.a > 560 & lambda.a <= 600));
    foo_difref = diff(p.ap(:, lambda.a > 500 & lambda.a <= 560),1, 2) ./ diff(lambda.a(lambda.a > 500 & lambda.a <= 560));
    toflag = any(abs(foo_difstep) > 3 * median(abs(foo_difref), 2, 'omitnan'), 2);
    if sum(toflag)
      fprintf('Signal of bubbles detected in %.2f%% (%i) of the ap spectra, flagged: step in ap between 550 and 600 nm\n', ...
        sum(toflag) / size(p, 1) * 100, sum(toflag))
      flag.ap_bubbles(toflag) = true;
    end
    % flag potential bubbles using the step in the center of the cp spectra
    foo_difstep = diff(p.cp(:, lambda.c > 560 & lambda.c <= 600),1, 2) ./ diff(lambda.c(lambda.c > 560 & lambda.c <= 600));
    foo_difref = diff(p.cp(:, lambda.c > 500 & lambda.c <= 560),1, 2) ./ diff(lambda.c(lambda.c > 500 & lambda.c <= 560));
    toflag = any(abs(foo_difstep) > 3 * median(abs(foo_difref), 2, 'omitnan'), 2);
    if sum(toflag)
      fprintf('Signal of bubbles detected in %.2f%% (%i) of the cp spectra, flagged: step in cp between 550 and 600 nm\n', ...
        sum(toflag) / size(p, 1) * 100, sum(toflag))
      flag.cp_bubbles(toflag) = true;
    end
  end
  
  % Auto QC when a positive first derivatives of ap over
  % wavelentght between 460 and 640 nm is larger than 0.4 times ap at 450nm
  % or second derivatives of ap over wavelentght between 460 and 640 nm is
  % larger than 0.006
  ap_450 = p.ap(:, find(lambda.a >= 450, 1,'first'));
  d460_640 = diff(p.ap(:, lambda.a > 460 & lambda.a <= 640),[],2);
  % delete bad spectra
  toflag = any(d460_640 > 0.4 * ap_450,2) | any(abs(diff(d460_640,[],2)) > 0.05, 2);
  if sum(toflag)
    fprintf('%.2f%% (%i) spectra flagged: d(ap)/d(lambda460-640) > 0.4 * ap_{450nm} | abs(d"(ap)/d(lambda460-640)) > 0.05)\n', ...
      sum(toflag) / size(p, 1) * 100, sum(toflag))
  end
  flag.ap460_640_04_450(toflag) = true;
  bad = [bad; p(toflag, :) table(repmat({'d(ap)/d(lambda460-640) > 0.4 * ap_{450nm} | abs(d"(ap)/d(lambda460-640)) > 0.05)'}, ...
    sum(toflag), 1), 'VariableNames', {'QC_failed'})];
  % p.ap(toflag, :) = NaN;
  
  % Auto QC when ap spectra contains 4 consecutive positive first derivatives of ap over
  % wavelentght between 485 and 570 nm
  d485_570 = diff(p.ap(:, lambda.a > 485 & lambda.a <= 570),[],2);
  pos_d485_570 = d485_570 > 0;
  toflag = false(size(pos_d485_570, 1), 1);
  N = 4; % Required number of consecutive numbers following a first one
  for i = 1:size(toflag,1)
    t = [false pos_d485_570(i,:) false];
    if any(find(diff(t)==-1)-find(diff(t)==1)>=N) % First t followed by >=N consecutive numbers
      toflag(i) = true;
    end
  end
  if sum(toflag)
    fprintf('%.2f%% (%i) spectra flagged: 4 consecutive d(ap)/d(lambda485-570) > 0\n', ...
      sum(toflag) / size(p, 1) * 100, sum(toflag))
  end
  flag.positive_ap450_570(toflag) = true;
  bad = [bad; p(toflag, :) table(repmat({'4 consecutive d(ap)/d(lambda485-570) > 0'}, ...
    sum(toflag), 1), 'VariableNames', {'QC_failed'})];
  bad = sortrows(bad, 'dt');
  % p.ap(toflag, :) = NaN;
  
  % % Auto QC when ap spectra contains 5 consecutive positive first derivatives of cp over
  % % wavelentght between 485 and 570 nm
  % d485_570 = diff(p.cp(:, lambda.c > 485 & lambda.c <= 570),[],2);
  % pos_d485_570 = d485_570 > 0;
  % todelete = false(size(pos_d485_570, 1), 1);
  % N = 5; % Required number of consecutive numbers following a first one
  % for i = 1:size(todelete,1)
  %   t = [false pos_d485_570(i,:) false];
  %   if any(find(diff(t)==-1)-find(diff(t)==1)>=N) % First t followed by >=N consecutive numbers
  %     todelete(i) = true;
  %   end
  % end
  % if sum(todelete)
  %   fprintf('%.2f%% (%i) spectra failed auto-QC: 5 consecutive d(cp)/d(lambda485-570) > 0\n', ...
  %     sum(todelete) / size(p, 1) * 100, sum(todelete))
  % end
  % bad = [bad; p(todelete, :) table(repmat({'5 consecutive d(cp)/d(lambda485-570) > 0'}, ...
  %   sum(todelete), 1), 'VariableNames', {'QC_failed'})];
  % bad = sortrows(bad, 'dt');
  % p.cp(todelete, :) = NaN;
  
  % run gaussian decomposition
  agaus = GaussDecomp(p, lambda.a, compute_ad_aphi);
  p = [p agaus];
  
  fprintf('Calculating Chl line height, POC & gamma ... ')
  % Derive standard products from ap and cp
  % Derive POC (Specific to region)
  % 	The particulate organic carbon (POC) is computed using the particulate attenuation at 660 nm Using the global relationship from Gardner et al. (2006)
  % Gardner, W.D., Mishonov, A., Richardson, M.J., 2006. Global POC concentrations from in-situ and satellite data. Deep Sea Res. II 53, 718?740.
  cp660 = interp1(lambda.c, p.cp',660,'linear')';
  p.poc = cp660.*380;
  flag.poc_flag(p.poc < 0) = true;
  % p.poc(p.poc < 0) = NaN;
  
  % Derive Chl (Line heigh at 676 compared to 650 and 715)
  % 	Chlorophyll a (chl) is computed using the particulate absorption line height at 676 nm and the global relationship from Tara Ocean (Boss et al. 2013)
  % REFERENCES:
  % Emmanuel Boss, Marc Picheral, Thomas Leeuw, Alison Chase, Eric Karsenti, Gabriel Gorsky, Lisa Taylor, Wayne Slade, Josephine Ras, Herve Claustre, 2013.The characteristics of particulate absorption, scattering and attenuation coefficients in the surface ocean; Contribution of the Tara Oceans expedition, Methods in Oceanography.
  ap_a = interp1(lambda.a, p.ap',[650 676 715],'linear')';
  p.ap(ap_a(:,1) > ap_a(:,2), :) = NaN; % deleted unrealistic spectra
  ap_a(ap_a(:,1) > ap_a(:,2), :) = NaN; % deleted unrealistic spectra
  p.ap676_lh = ap_a(:,2)-(39/65*ap_a(:,1)+26/65*ap_a(:,3));  % ap_a650-(39/65*ap_a(:,1)+26/65*ap_a(:,3));
  flag.chl_ap676lh_flag(real(p.ap676_lh) ~= p.ap676_lh) = true;
  p.ap676_lh(real(p.ap676_lh) ~= p.ap676_lh) = NaN;
  flag.ap676_lh(p.ap676_lh < 0) = true;
  p.ap676_lh(p.ap676_lh < 0) = NaN;
  p.chl_ap676lh = 157*p.ap676_lh.^1.22;
  
  % 3.3 Derive Gamma (does not support NaN values) (Boss et al. 2001)
  % REFERENCES:
  % Boss, E., W.S. Pegau, W.D. Gardner, J.R.V. Zaneveld, A.H. Barnard., M.S. Twardowski, G.C. Chang, and T.D. Dickey, 2001. Spectral particulate attenuation and particle size distribution in the bottom boundary layer of a continental shelf. Journal of Geophysical Research, 106, 9509-9516.
  % [~,p.gamma] = FitSpectra_HM2(lambda.ref(:,1:end-2),p.cp(:,1:end-2));
  % Correct bu on March 5, 2018, FitSpectra_HM2 does not accept NaNs in cp
  p.gamma = NaN(size(p, 1), 1);
  cp_temp_gam = p.cp(~all(isnan(p.cp),2), :);
  sel = ~any(isnan(cp_temp_gam));
  [~,temp_gam] = FitSpectra_HM2(lambda.c(:,sel), cp_temp_gam(:,sel));
  p.gamma(~all(isnan(p.cp),2)) = temp_gam;
  flag.gamma_flag(p.gamma < 0) = true;
  % p.gamma(p.gamma < 0) = NaN;
  fprintf('Done\n')
  
  fprintf('Calculating chlorophyll from cp (H_alh) ... ')
  % Chlorophyll absorption (alh) and phytoplankton size eigenvectors (P) inferred from cp
  % Houskeeper, H.F., Draper, D., Kudela, R.M., Boss, E., 2020. Chlorophyll absorption and phytoplankton size information inferred from hyperspectral particulate beam attenuation. Appl. Opt. 59, 6765. https://doi.org/10.1364/AO.396832
  % First compute hskpr P parameters (put link to github of hkpr do not include in github).
  [p.Halh, Pr] = houskeeperetal2020(lambda.c, p.cp, false);
  flag.chl_Halh_flag(p.Halh < 0) = true;
  p.Halh(p.Halh < 0) = NaN;
  p.chl_Halh = 157*p.Halh.^1.22;
  fprintf('Done\n')
  
  fprintf('Estimating G50 and mphi (slope of PSD) ... ')
  % Use fit from Haëntjens et al. 2021v22 to get median average
  % cross-sectional area (G50) and slope of phytoplankton size distribution
  % (in abundance) (mphi):
  p.HH_G50 = modelG50.predictFcn(Pr');
  p.HH_mphi = modelmphi.predictFcn(Pr');
  flag.HH_mphi_flag(p.HH_mphi > 0 | p.HH_mphi < -8) = true;
  % p.HH_mphi(p.HH_mphi > 0 | p.HH_mphi < -8) = NaN;
  flag.HH_G50_flag(p.HH_G50 < 0 | p.HH_G50 > 500) = true;
  % p.HH_G50(p.HH_G50 < 0 | p.HH_G50 > 50) = NaN;
  fprintf('Done\n')
  
  % Extra flags for suspicious data
  % gamma_suspicious
  flag.gamma_suspicious(p.gamma < 0.4) = true;
  % poc_suspicious
  flag.poc_suspicious(p.poc > 1000) = true;
  % chl_ap676lh_suspicious
  flag.chl_ap676lh_suspicious(p.chl_ap676lh > 30) = true;
  % chl_Halh_suspicious if ap676_lh is one order of magnitude different from Halh
  flag.chlratio_flag(p.ap676_lh./p.Halh <= 0.1 | p.ap676_lh./p.Halh >= 10 | p.chl_Halh > 30) = true;
  % HH_G50_mphi_suspicious
  flag.HH_G50_mphi_suspicious(p.HH_G50 < 0.3 | p.HH_mphi < -5 | p.HH_G50 > 50) = true;
  
  % set flag column
  p.flag_bit = set_flagbit(flag);
  % flag_info = read_flagbit(p.flag_bit, 'ACS');
  
  %% ag & cg
  if ~isempty(di)
    if strcmp(di_method, 'best_di')
      % select DIW with lowest a or c values between 550-650nm
      di_orig = di;
      di_dt = datetime(di_orig.dt, 'ConvertFrom', 'datenum');
      best_di_a = NaN(size(di_orig,1), 1);
      best_di_c = NaN(size(di_orig,1), 1);
      for i = 1:size(di_orig,1)
        if i == 1 || i == size(di_orig,1)
          iddi = abs(di_dt(i) - di_dt) < hours(72);
        else
          iddi = abs(di_dt(i) - di_dt) < hours(36);
        end
        lowest_di_a = di_orig.a(:, lambda.a >= 550 & lambda.a <= 650) == ...
          min(di_orig.a(iddi, lambda.a >= 550 & lambda.a <= 650), [], 1);
        lowest_di_c = di_orig.c(:, lambda.c >= 550 & lambda.c <= 650) == ...
          min(di_orig.c(iddi, lambda.c >= 550 & lambda.c <= 650), [], 1);
        foo_a = find(sum(lowest_di_a, 2) == max(sum(lowest_di_a, 2)));
        foo_c = find(sum(lowest_di_c, 2) == max(sum(lowest_di_c, 2)));
        di.a(i, :) = di_orig.a(foo_a, :);
        di.c(i, :) = di_orig.c(foo_c, :);
        di.a_avg_sd(i, :) = di_orig.a_avg_sd(foo_a, :);
        di.c_avg_sd(i, :) = di_orig.c_avg_sd(foo_c, :);
        
        best_di_a(i) = foo_a(1);
        best_di_c(i) = foo_c(1);
      end
    end
    
    % remove when a and c are full of NaNs
    filt_avg(all(isnan(filt_avg.a), 2) & all(isnan(filt_avg.c), 2),:) = [];
    
    % Interpolate filtered on Total
    di_interp = table(filt_avg.dt, 'VariableNames', {'dt'});
    di_interp.a = interp1(di.dt, di.a, di_interp.dt, 'linear', 'extrap');
    di_interp.c = interp1(di.dt, di.c, di_interp.dt, 'linear', 'extrap');
    di_interp.a_avg_sd = interp1(di.dt, di.a_avg_sd, di_interp.dt, 'linear', 'extrap');
    di_interp.c_avg_sd = interp1(di.dt, di.c_avg_sd, di_interp.dt, 'linear', 'extrap');
  
    % Dissolved = Filtered - DI
    g = table(filt_avg.dt, 'VariableNames', {'dt'});
    g.ag = filt_avg.a - di_interp.a;
    g.cg = filt_avg.c - di_interp.c;
    % Interpolate wavelength of c on a 
    % g.cg = cell2mat(arrayfun(@(i) interp1(c_wl, g.cg(i,:), a_wl, 'linear', 'extrap'), 1:size(g,1), 'UniformOutput', false)');
    g.ag = interp1(lambda.a', g.ag', lambda.ref', 'linear', 'extrap')';
    g.cg = interp1(lambda.c', g.cg', lambda.ref', 'linear', 'extrap')';
    fprintf('Correcting for temperature & salinity dependence ... ')
    % Temperature & Salinity Correction (No Scattering correction needed)
    [g.ag, g.cg] = TemperatureAndSalinityDependence(g.ag, g.cg, lambda.a, lambda.c, psi);
    fprintf('Done\n')
  
    % Propagate error
    %   Note: Error is not propagated through Scattering & Residual temperature
    %         correction as required by SeaBASS
    g.ag_sd = sqrt(filt_avg.a_avg_sd.^2 + di_interp.a_avg_sd.^2);
    g.cg_sd = sqrt(filt_avg.c_avg_sd.^2 + di_interp.c_avg_sd.^2);
    g.ag_n = filt_avg.a_avg_n;
    g.cg_n = filt_avg.c_avg_n;
    
  %   visProd3D(lambda.a, g.dt, g.ag, false, 'Wavelength', false, 70);
  %   title('ag with auto best DI 72h')
  %   saveGraph('ag_with_auto_best_DI_72h', 'fig')
  %   
  %   visProd3D(lambda.c, g.dt, g.cg, false, 'Wavelength', false, 71);
  %   title('cg with auto best DI 72h')
  %   saveGraph('cg_with_auto_best_DI_72h', 'fig')
  
  %   visProd3D(lambda.a, di_interp.dt, di_interp.a, false, 'Wavelength', false, 72);
  %   visProd3D(lambda.c, di_interp.dt, di_interp.c, false, 'Wavelength', false, 73);
  
    fprintf('Exponential fit to ag and cg ... ')
    sel_a = lambda.a < 660 & ~any(isnan(g.ag(~all(isnan(g.ag),2), :)));
    sel_c = lambda.c < 660 & ~any(isnan(g.cg(~all(isnan(g.cg),2), :)));
    [g.y_intercp_fit_ag, g.base_fit_ag, ~, ~, g.RMSE_fit_ag] = FitExp(lambda.a(sel_a), ...
      g.ag(:, sel_a), g.ag_sd(:, sel_a));
    % add fit flag
    g.ag_fitflag = false(size(g, 1), 1);
    g.ag_fitflag(g.RMSE_fit_ag > 0.0025) = true;
    [g.y_intercp_fit_cg, g.base_fit_cg, ~, ~, g.RMSE_fit_cg] = FitExp(lambda.c(sel_c), ...
      g.cg(:, sel_c), g.cg_sd(:, sel_c));
    % add fit flag
    g.cg_fitflag = false(size(g, 1), 1);
    g.cg_fitflag(g.RMSE_fit_cg > 0.0025) = true;
    fprintf('Done\n')
    
  %   % QC with ag and cg spectra (limited testing on the QC)
  %   g(g.ag(:,1) < 0 & g.cg(:,end-3) < -0.005, :) = [];
  else
    g = table();
  end
end


%% Temperature And Salinity Dependence correction
function [a_corr, c_corr, a_slope, c_slope] = TemperatureAndSalinityDependence(a, c, wl_a, wl_c, psi)
  % Note that a and c does not have to be on the same wavelength, however the function would need to be edited for that

  % Interpolate Sullivan values on the current ACS
  a_psiT=interp1(psi.wl, psi.psiT, wl_a,'spline'); % PCHIP or SPLINE -> better than linear
  c_psiT=interp1(psi.wl, psi.psiT, wl_c,'spline'); % PCHIP or SPLINE -> better than linear
  a_psiS=interp1(psi.wl, psi.a_psiS, wl_a,'spline');
  c_psiS=interp1(psi.wl, psi.c_psiS, wl_c,'spline');
  % Center psiS on 0 instead of +/- 0.001
  a_psiS = a_psiS - median(a_psiS(wl_a <= 590));
  c_psiS = c_psiS - median(c_psiS(wl_c <= 590));
  
  % Parameters of minization routine
  opts = optimset('fminsearch');      
  opts = optimset(opts,'MaxIter',20000000); 
  opts = optimset(opts,'MaxFunEvals',20000); % 20000
  opts = optimset(opts,'TolX',1e-8);
  opts = optimset(opts,'TolFun',1e-8);
  
  % Wavelength selection
  iwla = wl_a >= 450;
  iwlc = wl_c >= 450;
  % iwl = wl >= 700;
  
  % Init minimization parameters
  deltaT = 10;
  deltaS = 20;
  amp = a(find(~any(isnan(a),2), 1, 'first'),20);
  slope = 0.014;
  
  % Init loop
  n = size(a,1);
  a_deltaT = NaN(n,1);
  a_deltaS = NaN(n,1);
  a_amp = NaN(n,1);
  c_deltaT = NaN(n,1);
  c_deltaS = NaN(n,1);
  c_amp = NaN(n,1);
  a_slope = NaN(n,1);
  c_slope = NaN(n,1);
  
  % Run minimization
  for k=1:n
    % Force Temperature, Salinity, and Slope
    x = fminsearch(@costfun_TSD, [deltaT, deltaS, amp, slope], opts, a(k,iwla), a_psiT(iwla), a_psiS(iwla), wl_a(iwla));
    a_deltaT(k) = x(1); a_deltaS(k) = x(2); a_amp(k) = x(3); a_slope(k) = x(4);
    x = fminsearch(@costfun_TSD, [deltaT, deltaS, amp, slope], opts, c(k,iwlc), c_psiT(iwlc), c_psiS(iwlc), wl_c(iwlc));
    c_deltaT(k) = x(1); c_deltaS(k) = x(2); c_amp(k) = x(3); c_slope(k) = x(4);
    % Without Salinity forcing
  %   x = fminsearch(@costfun_TSD, [deltaT, amp, slope], opts, a(k,iwla), a_psiT(iwla), a_psiS(iwla), wl_a(iwla));
  %   a_deltaT(k) = x(1); a_amp(k) = x(2); a_slope(k) = x(3);
  %   x = fminsearch(@costfun_TSD, [deltaT, amp, slope], opts, c(k,iwlc), c_psiT(iwlc), c_psiS(iwlc), wl_c(iwlc));
  %   c_deltaT(k) = x(1); c_amp(k) = x(2); c_slope(k) = x(3);
  end
  
  % Apply correction
  % a_corr = a - a_deltaT.*a_psiT - deltaS.*a_psiS;
  % c_corr = c - c_deltaT.*c_psiT - deltaS.*c_psiS;
  a_corr = a - a_deltaT.*a_psiT - a_deltaS.*a_psiS;
  c_corr = c - c_deltaT.*c_psiT - c_deltaS.*c_psiS;
  
  % Display the forcing parameters
  % disp([a_deltaT, c_deltaT, a_deltaS, c_deltaS, a_slope, c_slope])
end


%%
function cost = costfun_TSD(x, a, psiT, psiS, wl)
  % Force Temperature, Salinity, and Slope
  cost = sum((a - psiT.*x(1) - psiS.*x(2) - x(3).*exp(-x(4)*(wl-450))).^2);
  % Without Salinity forcing
  % cost = sum((a - psiT.*x(1) - x(2).*exp(-x(3)*(wl-450))).^2);
end


%% Residual Temperature And Scattering Correction (Zaneveld 1994 proportional)
function [ap_corr, cp_corr] = ResidualTempScatterCorrZaneveld_proportional(ap, cp, wl, psi)
  % Function from Emmanuel Boss improved by Nils Haëntjens
  % Assumes negligible ap in NIR after Zaneveld 1994 method 3
  psiT = interp1(psi.wl, psi.psiT, wl);
  
  % Parameters of minization routine
  opts = optimset('fminsearch');
  % opts = optimset(opts,'NonlEqnAlgorithm', 'gn'); % Does not work on R2017a
  opts = optimset(opts,'MaxIter',20000000);
  opts = optimset(opts,'MaxFunEvals',20000);
  opts = optimset(opts,'TolX',1e-8);
  opts = optimset(opts,'TolFun',1e-8);
  
  % Find Near Infrared & references
  iNIR = 710 <= wl &  wl <= 750;  % spectral srange for optimization (710 to 750nm)
  if isempty(iNIR); error('Unable to perform correction as no wavelength available in NIR.'); end
  % Find nearest wavelength to greater than 730 nm to use as reference for correction
  iref = abs(wl - 715) == min(abs(wl - 715)); % 715 730
  % If ACS spectra do not go up to 730 nm take the closest wavelength to 730 nm
  if isempty(iref); [~, iref] = max(wl); end % works as there is data in iNIR so lowest wavelength is 710
  
  % Initialize output arrays
  deltaT = NaN(size(ap,1),1);
  
 % Init routine parameter scattering correction
  bp = cp - ap;
  
  % Run minimization routine on good spectra only
  sel = find(all(isfinite(ap),2));
  for k = sel'
    deltaT(k) = fminsearch(@costFun_RTSC, 0, opts, ap(k,:), bp(k,:), psiT, iNIR, iref);         
  end
  ap_corr = ap - psiT.*deltaT - ((ap(:,iref) - psiT(iref).*deltaT) ./ bp(:,iref)) .* bp; 
  cp_corr = cp - psiT.*deltaT;
end


%% Residual Temperature And Scattering Correction (Rottgers2013 semiempirical)
function [ap_corr, cp_corr] = ResidualTempScatterCorrRottgers_semiempirical(ap, cp, wl, psi)
  % Function from Emmanuel Boss after Rottgers et al. 2013
  % Allows absorption in NIR in case of high NAP
  psiT = interp1(psi.wl, psi.psiT, wl);
  
  % Parameters of minization routine
  opts = optimset('fminsearch');
  % opts = optimset(opts,'NonlEqnAlgorithm', 'gn'); % Does not work on R2017a
  opts = optimset(opts,'MaxIter',20000000); 
  opts = optimset(opts,'MaxFunEvals',20000);
  opts = optimset(opts,'TolX',1e-8);
  opts = optimset(opts,'TolFun',1e-8);
  
  % Find Near Infrared & references
  iNIR = 710 <= wl &  wl <= 750;  % spectral srange for optimization (710 to 750nm)
  if isempty(iNIR); error('Unable to perform correction as no wavelength available in NIR.'); end
  % Find nearest wavelength to greater than 730 nm to use as reference for correction
  iref = abs(wl - 715) == min(abs(wl - 715)); % 715 730
  % If ACS spectra do not go up to 730 nm take the closest wavelength to 730 nm
  if isempty(iref); [~, iref] = max(wl); end % works as there is data in iNIR so lowest wavelength is 710
  
  % Initialize output arrays
  deltaT = NaN(size(ap,1),1);
  
  % Init routine parameter scattering correction
  bp = cp - ap;
  
  % Run minimization routine on good spectra only
  sel = find(all(isfinite(ap),2));
  for k = sel'
    deltaT(k) = fminsearch(@costFun_RTSC, 0, opts, ap(k,:), bp(k,:), psiT, iNIR, iref);         
  end
  % apply temperature correction and replace negative values at iref by NaN to avoid complex solution leakage
  ap_Tcorr = ap - psiT.*deltaT;
  cp_corr = cp - psiT.*deltaT;
  ap_Tcorr_iref = ap_Tcorr(:,iref);
  ap_Tcorr_iref(ap_Tcorr_iref < 0) = NaN;
  % apply flat scattering correction
  ap_corr = ap_Tcorr - (ap_Tcorr_iref - 0.212*ap_Tcorr_iref.^1.135);
end

%% Residual Temperature And Scattering Correction (Zaneveld 1994 proportional + NIR offset from Rottgers2013 semiempirical)
function [ap_corr, cp_corr] = ResidualTempScatterCorrZaneveldRottgers_blended(ap, cp, wl, psi)
  % Function from Emmanuel Boss after Zaneveld 1994 method 3 and Rottgers et al. 2013
  % Scattering correction + allows absorption in NIR in case of high NAP
  psiT = interp1(psi.wl, psi.psiT, wl);
  
  % Parameters of minization routine
  opts = optimset('fminsearch');
  % opts = optimset(opts,'NonlEqnAlgorithm', 'gn'); % Does not work on R2017a
  opts = optimset(opts,'MaxIter',20000000); 
  opts = optimset(opts,'MaxFunEvals',20000);
  opts = optimset(opts,'TolX',1e-8);
  opts = optimset(opts,'TolFun',1e-8);
  
  % Find Near Infrared & references
  iNIR = 710 <= wl &  wl <= 750;  % spectral srange for optimization (710 to 750nm)
  if isempty(iNIR); error('Unable to perform correction as no wavelength available in NIR.'); end
  % Find nearest wavelength to greater than 730 nm to use as reference for correction
  iref = abs(wl - 715) == min(abs(wl - 715)); % 715 730
  % If ACS spectra do not go up to 730 nm take the closest wavelength to 730 nm
  if isempty(iref); [~, iref] = max(wl); end % works as there is data in iNIR so lowest wavelength is 710
  
  % Initialize output arrays
  deltaT = NaN(size(ap,1),1);
  
  % Init routine parameter scattering correction
  bp = cp - ap;
  
  % Run minimization routine on good spectra only
  sel = find(all(isfinite(ap),2));
  for k = sel'
    % deltaT(k) = fminsearch(@costFun_RTSC_RottgersBoss, 0, opts, ap(k,:), bp(k,:), psiT, iNIR, iref);         
    deltaT(k) = fminsearch(@costFun_RTSC, 0, opts, ap(k,:), bp(k,:), psiT, iNIR, iref);         
  end
  % apply temperature correction and replace negative values at iref by NaN to avoid complex solution leakage
  ap_Tcorr = ap - psiT.*deltaT;
  cp_corr = cp - psiT.*deltaT;
  ap_Tcorr_iref = ap_Tcorr(:,iref);
  ap_Tcorr_iref(ap_Tcorr_iref < 0) = NaN;
  bp_Tcorr = cp_corr - ap_Tcorr;
  % apply scattering correction
  ap_corr = ap_Tcorr - ap_Tcorr_iref ./ bp_Tcorr(:,iref) .* bp_Tcorr + 0.212*ap_Tcorr_iref.^1.135;
end

%% Residual Temperature And Scattering Correction cost function 
function cost = costFun_RTSC(deltaT, ap, bp, psiT, iNIR, iref)
  cost = sum(abs(ap(iNIR) - psiT(iNIR).*deltaT - ((ap(iref)-psiT(iref).*deltaT)./bp(iref)).*bp(iNIR)));
end

%% Temperature And Salinity Correction
function [a_ts, c_ts] = TemperatureAndSalinityCorrection(a, c, a_wl, c_wl, delta_t, delta_s, psi)

  %interpolate literature psiT values to acs wavelengths
  a_psi_t = interp1(psi.wl, psi.psiT, a_wl, 'linear', 'extrap');
  c_psi_t = interp1(psi.wl, psi.psiT, c_wl, 'linear', 'extrap');
  a_psi_s = interp1(psi.wl, psi.a_psiS, a_wl, 'linear', 'extrap');
  c_psi_s = interp1(psi.wl, psi.c_psiS, c_wl, 'linear', 'extrap');
  %a_sigma_psi_t = interp1(psi.wl, psi.a_sigma_psiS, wl_a, 'linear', 'extrap');
  %c_sigma_psi_t = interp1(psi.wl, psi.c_sigma_psiS, wl_c, 'linear', 'extrap');
    
  %correct acdata_raw for temp-dependent water absorbance, propagate error
  %into del_raw.  Output acdata_t and del_t.
   
  a_ts = a - (a_psi_t(ones(size(delta_t,1),1),:) .* delta_t(:,ones(size(a_psi_t,2),1)))...
    + (a_psi_s(ones(size(delta_s,1),1),:) .* delta_s(:,ones(size(a_psi_s,2),1)));
  c_ts = c - (c_psi_t(ones(size(delta_t,1),1),:) .* delta_t(:,ones(size(c_psi_t,2),1)))...
    + (c_psi_s(ones(size(delta_s,1),1),:) .* delta_s(:,ones(size(c_psi_s,2),1)));
  
  %c_ts = c - (c_psi_t .* delta_t + c_psi_s .* delta_s);
  %a_del_t(:,1) = ((del_rawa(:,1)).^2 + (a_sigma_psi_t(:,1).*delta_t).^2).^(1/2);
  %c_del_t(:,1) = ((del_rawc(:,1)).^2 + (c_sigma_psi_t(:,1).*delta_t).^2).^(1/2);
end


%%
function a_corr = ScatteringCorrection(a_wl, c_wl, a, c, method)
  % Interpolate wavelength of c on a
  c = interp1(c_wl', c, a_wl', 'linear', 'extrap');
  
  switch method
    case 'flat'
      % spectrally flat correction
      a_730 = interp1(wl, a, 730, 'linear', 'extrap');
      a_corr = a - a_730;
    case 'varying'
      % spectrally varying scattering correction
      b = c - a;
      a_730 = interp1(wl, a, 730, 'linear', 'extrap');
      b_730 = interp1(wl, b, 730, 'linear', 'extrap');
      a_corr = a - b .* a_730 ./ b_730;
    case 'rottgers'
      % Rottgers et al., 2013
      a_715 = interp1(wl, a, 715, 'linear', 'extrap');
      c_715 = interp1(wl, c, 715, 'linear', 'extrap');
      a_corr = a - a_715 .* (1/0.56 .* c - a)./(1/0.56 .* c_715 - a_715);
    otherwise
      error('Method not supported');
  end
end


%%
function acs_unsmoothed = unsmoothACS(acs_data, lambda)
  % AC-S "un-smoothing" and spectral decomposition method
  % Ron Zaneveld, WET Labs, Inc., 2005
  % Ali Chase, University of Maine, 2014
  %
  % Unsmoothing method from spectral decomposition in:
  % Chase, A., et al., Decomposition of in situ particulate absorption
  % spectra. Methods in Oceanography (2014), http://dx.doi.org/10.1016/j.mio.2014.02.022
  %%
  todo = acs_data.Properties.VariableNames(contains(acs_data.Properties.VariableNames, ...
    {'ap', 'cp'}) & ~contains(acs_data.Properties.VariableNames, {'_sd','_n'}));
  
  acs_unsmoothed = table();
  acs_unsmoothed.dt = acs_data.dt;
  for i = todo
    fprintf([i{:} ' unsmoothing ... '])
    % Set up filter factors at every 0.1 nm from 1 to 799 nm, with center
    % wavelength at centwavel (i.e. at the data wavelengths)
    wavelength = .1:.1:799; % Thus index of 1 nm = 10; 356 nm= 3560;
    SIG1 = (-9.845*10^-8.*lambda.(i{:}(1)).^3 + 1.639*10^-4*lambda.(i{:}(1)).^2 - 7.849*10^-2*lambda.(i{:}(1)) + 25.24)/2.3547 ;
    for j = 1:size(lambda.(i{:}(1)),2)
      for jkl = 1:size(wavelength,2)
        filtfunc(jkl,j) = (1/(sqrt(2*pi)*SIG1(j)))*exp(-0.5*((wavelength(jkl)-lambda.(i{:}(1))(j))/SIG1(j)).^2); % First term normalizes area under the curve to 1.
      end
    end
  
    % Convolve the measurement with the fiter factors add the difference to
    % the measured spectrum to get the first corrected spectrum.
    % This is the corrected absorption spectrum "ap".
    minwavel = min(lambda.(i{:}(1)));
    maxwavel = max(lambda.(i{:}(1)));
  
    centwavel = minwavel:.1:maxwavel;% The range of centwavel is 0.1 nm.
    splinap = spline(lambda.(i{:}(1)), acs_data.(i{:}), centwavel); % Spline the measured data to every 0.1 nm.
    % We need data from 0 to 799 nm to multiply by filtfac.
    absspec = zeros(size(acs_data.(i{:}),1), size(wavelength,2));
    absspec(:, minwavel*10:maxwavel*10) = splinap;
    absspec(:, 1:minwavel*10-1) = ones(1, size(1:minwavel*10-1,2)) .* absspec(:, minwavel*10);
    aspecprime = absspec';
  
    meassignal6 = NaN(size(aspecprime, 2), size(lambda.(i{:}(1)), 2));
    parfor j = 1:size(aspecprime, 2)        
      measur2 = aspecprime(:,j) .* filtfunc; % the measured signal for every filter factor.
      meassignal6(j,:) = 0.1 * sum(measur2); % The measured spectrum at a wavelength i is the sum of what a filter measured at
    end
    acs_unsmoothed.(i{:}) = acs_data.(i{:}) - meassignal6 + acs_data.(i{:});
    fprintf('Done\n')
  end
  
  for i = 1:size(acs_data,2)
    if ~any(strcmp(acs_unsmoothed.Properties.VariableNames, acs_data.Properties.VariableNames{i}))
      acs_unsmoothed = [acs_unsmoothed acs_data(:,i)];
    end
  end
end


%%
function agaus = GaussDecomp(p, lambda, compute_ad_aphi)
  % Gaussian decomposition
  % Ron Zaneveld, WET Labs, Inc., 2005
  % Ali Chase, University of Maine, 2014
  % Adaptation to InLineAnalysis: Guillaume Bourdin, March 2021
  %
  % Reference:
  % Chase, A., et al., Decomposition of in situ particulate absorption
  % spectra. Methods in Oceanography (2014), http://dx.doi.org/10.1016/j.mio.2014.02.022
  
  %% identify lines full of NaNs to reconstruct table of same size
  idnan = all(isnan(p.ap),2);
  % delete row full of NaN
  p(idnan, :) = [];
  
  % delete row and columns full of nans
  lambda(all(isnan(p.ap),1)) = [];
  p.ap_sd(:, all(isnan(p.ap),1)) = [];
  p.ap(:, all(isnan(p.ap),1)) = [];
  
  % Peak center values ("peak_loc") determined using a interative approach that allows the location to vary (uses the matlab
  % function LSQNONLIN), and are rounded to nearest integer. Sigma values ("lsqsig") are determined similarly. FWHM = sigma*2.355
  peak_loc = [406,434,453,470,492,523,550,584,617,638,660,675];
  lsqsig = [16,12,12,13,16,14,14,16,13,11,11,10];
  onenm = 400:1:720;
  
  fprintf('Gaussian decomposition ... ')
  
  % extrapolate NaN tails with constant nearest non-NaN value
  ap_filled = fillmissing(p.ap, 'nearest', 2, 'SamplePoints', lambda);
  ap_sd_filled = fillmissing(p.ap_sd, 'nearest', 2, 'SamplePoints', lambda);
  
  % interpolate the un-smoothed ap spectra to one nm resolution
  acorr2onenm = interp1(lambda, ap_filled', onenm, 'spline');
  
  % define the matrix of component Gaussian functions using the peaks
  % and widths (sigma) above
  coef2 = exp(-0.5 .* (((onenm .* ones(size(peak_loc,2),1))' - peak_loc .* ...
    ones(size(onenm,2),1)) ./ lsqsig) .^ 2);
  
  % define a function for non-algal particles and concatenate this to the Gaussian matrix
  coef2nap = exp(-0.01 * (onenm - 400));
  coef2 = [coef2nap', coef2];
  
  % normalize both the component functions and the measured ap
  % spectra by the uncertainty (std dev) in the ap spectra
  ap_sd_int = interp1(lambda, ap_sd_filled', onenm, 'linear', 'extrap');
  acorr2onenm_new = (acorr2onenm ./ ap_sd_int);
  
  amps = NaN(size(p,1), size(coef2,2));
  sumspec_temp = NaN(size(coef2,1), size(p,1));
  % compspec_temp = NaN(size(coef2,1), size(coef2,2), size(p,1));
  
  parfor i = 1:size(acorr2onenm_new,2)
    coef2_new = coef2 ./ ap_sd_int(:,i);
    % Inversion analysis
    amps(i, :) = lsqnonneg(coef2_new, acorr2onenm_new(:, i));
    % Using the inverted amplitudes, build a matrix of new component
    % spectra (compspec) and the sum of the Gaussian and nap functions (sumspec):
    sumspec_temp(:, i) = sum(amps(i, :) .* coef2, 2);
  %   compspec_temp(:, :, i) = amps(i, :) .* coef2;
  end
  
  % interpolate back to the original resolution
  % compspec = interp1(onenm', compspec_temp, lambda, 'spline');
  sumspec = interp1(onenm, sumspec_temp, lambda, 'spline')';
  
  uncertainty = sum(abs(p.ap(:, lambda > 440 & lambda < 705) - ...
    sumspec(:, lambda > 440 & lambda < 705)), 2, 'omitnan') / ...
      size(lambda(lambda > 440 & lambda < 705),2);
    
  % reconstruct matrix with NaN in the same place
  agaus = NaN(size(idnan,1), 14);
  agaus(~idnan,:) = [amps uncertainty];
  agaus = array2table(agaus, 'VariableNames', ...
      [{'ad_model400'} cellfun(@(x) ['agaus' x], cellstr(num2str(peak_loc'))', 'un', 0) {'agaus-mae'}]);
  
  if compute_ad_aphi
    % Compute a non-algal and a phytoplankton from Zheng, G., and D. Stramski (2013) model
    % Reference: 
    % - Zheng, G., and D. Stramski (2013), A model based on stacked-constraints approach for partitioning the light absorption coefficient of seawater 
    % into phytoplankton and non-phytoplankton components, J. Geophys. Res. Oceans, 118, 2155?2174, doi:10.1002/jgrc.20115.
    % - Zheng, G., and D. Stramski (2013), A model for partitioning the light absorption coefficient of suspended marine particles into phytoplankton 
    % and nonalgal components, J. Geophys. Res. Oceans, 118, 2977?2991, doi:10.1002/jgrc.20206.
    wl_ZS13 = [400 412 420 430 443 450 467 490 500 510 550 555 630 650 670 700];
    qwl = unique([wl_ZS13(:)' 442 676]);
    ap_ZS13 = interp1(lambda, ap_filled', qwl, 'linear', 'extrap'); % need to extrapolate for 400
  
  %   % vectorized partition_ap (DOESN'T WORK YET)
  %   agaus.ad_ZS13 = NaN(size(agaus, 1), size(qwl, 2));
  %   agaus.aphi_ZS13 = NaN(size(agaus, 1), size(qwl, 2));
  %   [agaus.ad_ZS13, agaus.aphi_ZS13] = partition_ap_vec(ap_ZS13, qwl', 50);  % Must be in colum direction
    
    % Run partition_ap
    ad_ZS13 = NaN(size(ap_ZS13, 2), size(qwl, 2));
    aphi_ZS13 = NaN(size(ap_ZS13, 2), size(qwl, 2));
  %   parfor i = 1:size(ap_ZS13, 2)
    for i = progress(1:size(ap_ZS13, 2))
      [ad_ZS13(i,:), aphi_ZS13(i,:)] = partition_ap(ap_ZS13(:, i), qwl', 50);
    end
  %   figure; hold on
  %   plot(qwl, ad)
  %   plot(qwl, aphi_ZS13)
  
    agaus.ad_ZS13 = NaN(size(agaus, 1), size(qwl, 2));
    agaus.aphi_ZS13 = NaN(size(agaus, 1), size(qwl, 2));
    agaus.ad_ZS13(~idnan,:) = ad_ZS13;
    agaus.aphi_ZS13(~idnan,:) = aphi_ZS13;
    agaus.Properties.VariableUnits = repmat({'1/m'}, 1, size(agaus, 2));
  end
  fprintf('Done\n')
end

function [filt_interp] = agcg_fdom_interpolation(tot, filt_interp, filt, lambda, ...
  cdom_dark, filt_avg, min_nb_pts_per_cluster, time_weight_for_cluster)
  % Reconstruct ag and cg from fCDOM data
  [filt_avg, filt_interp, ~] = fdom_agcg_model(filt, filt_interp, lambda, filt_avg, ...
    min_nb_pts_per_cluster, time_weight_for_cluster);
  
  n_periods = size(filt_avg,1)-1;

  % visProd3D(lambda.a, filt_interp.dt, filt_interp.intercept_interp_a, false, 'Wavelength', false, 27);
  % title('intercept a')
  % visProd3D(lambda.a, filt_interp.dt, filt_interp.intercept_interp_c, false, 'Wavelength', false, 28);
  % title('intercept c')
  % 
  % visProd3D(lambda.a, filt_interp.dt, filt_interp.slope_interp_a, false, 'Wavelength', false, 27);
  % title('slope a')
  % visProd3D(lambda.a, filt_interp.dt, filt_interp.slope_interp_c, false, 'Wavelength', false, 28);
  % title('slope c')

  % % interpolate based on fCDOM: Emmanuel's method (local)
  % % allocate fdom0, a0, and c0 varialbes
  % filt_interp.fdom0 = NaN(size(filt_interp.fdom));
  % filt_interp.a0 = NaN(size(filt_interp.a));
  % filt_interp.c0 = NaN(size(filt_interp.c));
  % 
  % % interpolate based on fCDOM: Emmanuel's method
  % delta_acf = table(NaN(n_periods, 1), 'VariableNames', {'dt'});
  % delta_acf.delta_af = NaN(n_periods, size(filt_avg.a, 2));
  % delta_acf.delta_cf = NaN(n_periods, size(filt_avg.a, 2));
  % % Prepare coefficients for each total period going from t0 to t1, starting and finishing by a filtered time
  % for i=1:n_periods
  %   it0 = i; it1 = i + 1;
  %   it = filt_avg.dt(it0) <= filt_interp.dt & filt_interp.dt <= filt_avg.dt(it1);
  %   if any(it)
  %     % interpolation based on fCDOM if enough dynamic range in fCDOM
  %     if any((filt_interp.fdom(it) - cdom_dark) / cdom_dark > 0.05)
  %       % prepare delta a / delta f and delta c / delta f
  %       delta_acf.dt(i) = mean([filt_avg.dt(it0) filt_avg.dt(it1)], 'omitnan');
  %       delta_acf.delta_af(i, :) = (filt_avg.a(it1,:) - filt_avg.a(it0,:)) / (filt_avg.fdom(it1) - filt_avg.fdom(it0));
  %       delta_acf.delta_cf(i, :) = (filt_avg.c(it1,:) - filt_avg.c(it0,:)) / (filt_avg.fdom(it1) - filt_avg.fdom(it0));
  %       filt_interp.fdom0(it) = repmat(filt_avg.fdom(it0), sum(it), 1);
  %       filt_interp.a0(it,:) = repmat(filt_avg.a(it0,:), sum(it), 1);
  %       filt_interp.c0(it,:) = repmat(filt_avg.c(it0,:), sum(it), 1);
  %       filt_interp.flag_linear_interp(it) = false;
  %     else
  %       filt_interp.flag_linear_interp(it) = true;
  %     end
  %   end
  % end
  % delta_acf(all(isnan(delta_acf.delta_af), 2) | all(isnan(delta_acf.delta_cf), 2), :) = [];
  % % interpolated delta a / delta fdom over total events
  % filt_interp.delta_af = interp_extrap(delta_acf, filt_interp.dt, 'delta_af', [], true, 'linear', 'nearest');
  % filt_interp.delta_cf = interp_extrap(delta_acf, filt_interp.dt, 'delta_cf', [], true, 'linear', 'nearest');
  % % extrapolate a0, c0, fdom0
  % filt_interp.fdom0 = fillmissing(filt_interp.fdom0, 'nearest', 'SamplePoints', filt_interp.dt);
  % filt_interp.a0 = fillmissing(filt_interp.a0, 'nearest', 'SamplePoints', filt_interp.dt);
  % filt_interp.c0 = fillmissing(filt_interp.c0, 'nearest', 'SamplePoints', filt_interp.dt);
  % % compute ag and cg
  % filt_interp.a = filt_interp.al + filt_interp.delta_af .* (filt_interp.fdom - filt_interp.fdoml);
  % filt_interp.c = filt_interp.cl + filt_interp.delta_cf .* (filt_interp.fdom - filt_interp.fdoml);
  % % linearly interpolate ag and cg
  % filt_interp.lin_a = interp1(filt_avg.dt, filt_avg.a, filt_interp.dt, 'linear');
  % filt_interp.lin_c = interp1(filt_avg.dt, filt_avg.c, filt_interp.dt, 'linear');
  % % replace fdom interpolation by linear for filter events where fdom is constant
  % filt_interp.a(filt_interp.flag_linear_interp) = filt_interp.lin_a(filt_interp.flag_linear_interp);
  % filt_interp.c(filt_interp.flag_linear_interp) = filt_interp.lin_c(filt_interp.flag_linear_interp);

  % visProd3D(lambda.a, filt_interp.dt, filt_interp.a, false, 'Wavelength', false, 27);
  % title('ag fdom interpolated')
  % visProd3D(lambda.a, filt_interp.dt, filt_interp.lin_a, false, 'Wavelength', false, 28);
  % title('ag linearly interpolated')
  % 
  % visProd3D(lambda.c, filt_interp.dt, filt_interp.c, false, 'Wavelength', false, 29);
  % title('cg fdom interpolated')
  % visProd3D(lambda.c, filt_interp.dt, filt_interp.lin_c, false, 'Wavelength', false, 30);
  % title('cg linearly interpolated')

  % interpolate ag and cg based on fdom: Guillaume's method (global)
  % allocate new variables
  filt_interp.ag_fact = NaN(size(filt_interp.a));
  filt_interp.cg_fact = NaN(size(filt_interp.c));
  filt_interp.fCDOM_mix_a_cluster = false(size(filt_interp.dt));
  filt_interp.fCDOM_mix_c_cluster = false(size(filt_interp.dt));
  filt_interp.fCDOM_a_cluster_chg = false(size(filt_interp.dt));
  filt_interp.fCDOM_c_cluster_chg = false(size(filt_interp.dt));
  filt_interp.flag_linear_interp = false(size(filt_interp.dt));

  %%%%% REWRITTEN ON 2023-12-29
  % linearly interpolate a, c, and fdom
  % filt_interp.lin_fdom = fillmissing(filt_avg.fdom, 'linear', 'SamplePoints', filt_avg.dt);
  filt_interp.lin_a = interp1(filt_avg.dt, filt_avg.a, filt_interp.dt, 'linear');
  filt_interp.lin_c = interp1(filt_avg.dt, filt_avg.c, filt_interp.dt, 'linear');
  filt_interp.lin_fdom = interp1(filt_avg.dt, filt_avg.fdom, filt_interp.dt, 'linear');
  filt_interp.lin_a = fillmissing(filt_interp.lin_a, 'nearest', 'SamplePoints', filt_interp.dt);
  filt_interp.lin_c = fillmissing(filt_interp.lin_c, 'nearest', 'SamplePoints', filt_interp.dt);
  filt_interp.lin_fdom = fillmissing(filt_interp.lin_fdom, 'nearest', 'SamplePoints', filt_interp.dt);
  for i=1:n_periods
    it0 = i; it1 = i + 1;
    it = filt_avg.dt(it0) <= filt_interp.dt & filt_interp.dt <= filt_avg.dt(it1);
    if any(it)
      it_filt_interp = filt_interp(it,:);
      % interpolation based on fCDOM if enough dynamic range in fCDOM
      if any((filt_interp.fdom(it) - cdom_dark) / cdom_dark > -0.05*2000) && ...
          (any(~isnan(filt_interp.fdom(it))) || (any(~isnan(filt_interp.a(it))) && any(~isnan(filt_interp.c(it)))))
        % compute interpolated ag factor
        ag_fact0 = (filt_avg.a(it0,:) - filt_avg.intercept_a(it0,:)) ./ filt_avg.fdom(it0);
        ag_fact1 = (filt_avg.a(it1,:) - filt_avg.intercept_a(it1,:)) ./ filt_avg.fdom(it1);
        filt_interp.ag_fact(it,:) = interp1(filt_avg.dt(it0:it1,:), [abs(ag_fact0); abs(ag_fact1)], it_filt_interp.dt, 'linear');
        % compute interpolated cg factor
        cg_fact0 = (filt_avg.c(it0,:) - filt_avg.intercept_c(it0,:)) ./ filt_avg.fdom(it0);
        cg_fact1 = (filt_avg.c(it1,:) - filt_avg.intercept_c(it1,:)) ./ filt_avg.fdom(it1);
        filt_interp.cg_fact(it,:) = interp1(filt_avg.dt(it0:it1,:), [cg_fact0; cg_fact1], it_filt_interp.dt, 'linear');
        % % compute ag and cg (method using only a0 constant accross total event
        % filt_interp.a(it, :) = filt_avg.a(it0,:) + ag_fact .* (filt_interp.fdom(it,:) - filt_avg.fdom(it0,:));
        % filt_interp.c(it, :) = filt_avg.c(it0,:) + cg_fact .* (filt_interp.fdom(it,:) - filt_avg.fdom(it0,:));

        %%%%% AS WRITTEN BEFORE 2023-12-29
        % % compute interpolated fdom0, a0, and c0 variables
        % filt_interp.fdoml(it) = interp1(filt_avg.dt(it0:it1,:), [filt_avg.fdom(it0); filt_avg.fdom(it1)], it_filt_interp.dt, 'linear');
        % filt_interp.al(it,:) = interp1(filt_avg.dt(it0:it1,:), [filt_avg.a(it0,:); filt_avg.a(it1,:)], it_filt_interp.dt, 'linear');
        % filt_interp.cl(it,:) = interp1(filt_avg.dt(it0:it1,:), [filt_avg.c(it0,:); filt_avg.c(it1,:)], it_filt_interp.dt, 'linear');

        % flag when mixed a clusters within filter event 0 or filter event 1
        if any(filt_avg.cluster_a_flag(it0) | filt_avg.cluster_a_flag(it1))
          fprintf('a filter interpolation #%i flagged: mixed cluster within start and/or end filter event(s)\n', i)
          filt_interp.fCDOM_mix_a_cluster(i) = true;
        end
        % flag when mixed c clusters within filter event 0 or filter event 1
        if any(filt_avg.cluster_c_flag(it0) | filt_avg.cluster_c_flag(it1))
          fprintf('c filter interpolation #%i flagged: mixed cluster within start and/or end filter event(s)\n', i)
          filt_interp.fCDOM_mix_c_cluster(i) = true;
        end
        % flag when a clusters change between filter event 0 and filter event 1
        if ~all(isnan(filt_avg.cluster_weighted_a(it0,:))) && ~all(isnan(filt_avg.cluster_weighted_a(it1,:))) && ...
            any(filt_avg.cluster_weighted_a(it0,:) ~= filt_avg.cluster_weighted_a(it1, :)) && any(any(~isnan(tot.a(it, :)), 2))
          fprintf('a filter interpolation #%i flagged: change of clusters between start and end filter events\n', i)
          filt_interp.fCDOM_a_cluster_chg(it) = true;
        end
        % flag when c clusters change between filter event 0 and filter event 1
        if ~all(isnan(filt_avg.cluster_weighted_c(it0,:))) && ~all(isnan(filt_avg.cluster_weighted_c(it1,:))) && ...
            any(filt_avg.cluster_weighted_c(it0,:) ~= filt_avg.cluster_weighted_c(it1, :)) && any(any(~isnan(tot.c(it, :)), 2))
          fprintf('c filter interpolation #%i flagged: change of clusters between start and end filter events\n', i)
          filt_interp.fCDOM_c_cluster_chg(it) = true;
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
  % fill missing ag_fact and cg_fact
  filt_interp.ag_fact = fillmissing(filt_interp.ag_fact, 'nearest', 'SamplePoints', filt_interp.dt);
  filt_interp.cg_fact = fillmissing(filt_interp.cg_fact, 'nearest', 'SamplePoints', filt_interp.dt);
  % compute ag and cg
  filt_interp.a(~filt_interp.flag_linear_interp, :) = filt_interp.lin_a(~filt_interp.flag_linear_interp, :) + ...
    filt_interp.ag_fact(~filt_interp.flag_linear_interp, :) .* ...
    (filt_interp.fdom(~filt_interp.flag_linear_interp) - filt_interp.lin_fdom(~filt_interp.flag_linear_interp)); 
  filt_interp.c(~filt_interp.flag_linear_interp, :) = filt_interp.lin_c(~filt_interp.flag_linear_interp, :) + ...
    filt_interp.cg_fact(~filt_interp.flag_linear_interp, :) .* ...
    (filt_interp.fdom(~filt_interp.flag_linear_interp) - filt_interp.lin_fdom(~filt_interp.flag_linear_interp));
  % linear interpolation when fdom interpolation was not used
  filt_interp.a(filt_interp.flag_linear_interp, :) = filt_interp.lin_a(filt_interp.flag_linear_interp, :);
  filt_interp.c(filt_interp.flag_linear_interp, :) = filt_interp.lin_c(filt_interp.flag_linear_interp, :);

  % figure; hold on
  % scatter(datetime(filt_interp.dt, 'convertfrom', 'datenum'), ag_fact_tot(:,70), 20, 'ko')
  % scatter(datetime(filt_interp.dt(~filt_interp.flag_linear_interp), 'convertfrom', 'datenum'), filt_interp.delta_af(~filt_interp.flag_linear_interp,70), 20, 'r+')
  % 
  % figure; hold on
  % scatter(filt_interp.dt, cg_fact_tot(:,70), 20, 'ko')
  % scatter(filt_interp.dt(~filt_interp.flag_linear_interp), filt_interp.delta_cf(~filt_interp.flag_linear_interp,70), 20, 'r+')

  % visProd3D(lambda.a, filt_interp.dt(it), filt_interp.a(it,:), false, 'Wavelength', false, i+1);
  % title('ag fdom interpolated')
  % visProd3D(lambda.a, filt_interp.dt(it), filt_interp.lin_a(it,:), false, 'Wavelength', false, i+2);
  % title('ag linearly interpolated')
  % 
  % visProd3D(lambda.c, filt_interp.dt(it), filt_interp.c(it,:), false, 'Wavelength', false, i+3);
  % title('cg fdom interpolated')
  % visProd3D(lambda.c, filt_interp.dt(it), filt_interp.lin_c(it,:), false, 'Wavelength', false, i+4);
  % title('cg linearly interpolated')
  % 
  % visProd3D(lambda.a, filt_interp.dt(it), filt_interp.a(it,:) - filt_interp.lin_a(it,:), false, 'Wavelength', false, i+5);
  % title('ag difference')
  % visProd3D(lambda.a, filt_interp.dt(it), filt_interp.c(it,:) - filt_interp.lin_c(it,:), false, 'Wavelength', false, i+6);
  % title('cg difference')
  %
  % visProd3D(lambda.a, filt_interp.dt, filt_interp.a, false, 'Wavelength', false, i+1);
  % title('ag fdom interpolated')
  % visProd3D(lambda.c, filt_interp.dt, filt_interp.c, false, 'Wavelength', false, i+2);
  % title('cg fdom interpolated')

  % yaxr = [it_filt_interp.fdom; filt_avg.fdom(it0:it1)];
  % 
  % wl_toplot = 40;
  % figure(1)
  % clf
  % subplot(1,2,1)
  % yyaxis('left')
  % hold on
  % plot_dt = datetime(it_filt_interp.dt, 'ConvertFrom', 'datenum');
  % scatter(plot_dt, ...
  %   filt_interp.a(it,wl_toplot), 40, 'filled', 'MarkerFaceColor', 'b','MarkerFaceAlpha', 0.5)
  % scatter(datetime(filt_avg.dt, 'ConvertFrom', 'datenum'), ...
  %   filt_avg.a(:,wl_toplot), 100, 'filled', 'MarkerFaceColor', [255	193	37]/255,'MarkerFaceAlpha', 1)
  % xlim([plot_dt(1)-hours(0.5) plot_dt(end)+hours(0.5)])
  % ylabel('ag')
  % yyaxis('right')
  % hold on
  % scatter(datetime(filt_interp.dt, 'ConvertFrom', 'datenum'), ...
  %   filt_interp.fdom, 40, 'filled', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5, 'Marker', 'd')
  % scatter(datetime(filt_avg.dt, 'ConvertFrom', 'datenum'), ...
  %   filt_avg.fdom, 100, 'filled', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5, 'Marker', 'd')
  % xlim([plot_dt(1)-hours(0.5) plot_dt(end)+hours(0.5)])
  % ylim([min(yaxr)-0.001 max(yaxr)+0.001])
  % ylabel('fdom')
  % legend('ag fdom interpolated', 'ag filter average', 'fdom smoothed', 'Location', 'Best')
  % subplot(1,2,2)
  % yyaxis('left')
  % hold on
  % scatter(plot_dt, ...
  %   filt_interp.c(it,wl_toplot), 40, 'filled', 'MarkerFaceColor', 'r','MarkerFaceAlpha', 0.5)
  % scatter(datetime(filt_avg.dt, 'ConvertFrom', 'datenum'), ...
  %   filt_avg.c(:,wl_toplot), 100, 'filled', 'MarkerFaceColor', [255	193	37]/255,'MarkerFaceAlpha', 1)
  % xlim([plot_dt(1)-hours(0.5) plot_dt(end)+hours(0.5)])
  % ylabel('cg')
  % yyaxis('right')
  % hold on
  % scatter(datetime(filt_interp.dt, 'ConvertFrom', 'datenum'), ...
  %   filt_interp.fdom, 40, 'filled', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5, 'Marker', 'd')
  % scatter(datetime(filt_avg.dt, 'ConvertFrom', 'datenum'), ...
  %   filt_avg.fdom, 100, 'filled', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5, 'Marker', 'd')
  % xlim([plot_dt(1)-hours(0.5) plot_dt(end)+hours(0.5)])
  % ylim([min(yaxr)-0.001 max(yaxr)+0.001])
  % ylabel('fdom')
  % legend('cg fdom interpolated', 'cg filter average', 'fdom smoothed', 'Location', 'Best')

  % visProd3D(lambda.a, filt_interp.dt, filt_interp.a, false, 'Wavelength', false, 23);
  % title('ag fdom interpolated')
  % visProd3D(lambda.a, filt_interp.dt, filt_interp.lin_a, false, 'Wavelength', false, 24);
  % title('ag linearly interpolated')
  % 
  % visProd3D(lambda.c, filt_interp.dt, filt_interp.c, false, 'Wavelength', false, 25);
  % title('cg fdom interpolated')
  % visProd3D(lambda.c, filt_interp.dt, filt_interp.lin_c, false, 'Wavelength', false, 26);
  % title('cg linearly interpolated')

end

%% 
function [filt_avg, filt_interp, regress_stats] = fdom_agcg_model(filt, filt_interp, ...
  lambda, filt_avg, min_nb_pts_per_cluster, time_weight_for_cluster)
  
  % min_nb_pts_per_cluster = 100;
  % time_weight_for_cluster = 0.9;

  clusters = table();
  % define reference wavelength to do regression (450nm)
  id450_a = abs(lambda.a - 450) == min(abs(lambda.a - 450));
  id450_c = abs(lambda.c - 450) == min(abs(lambda.c - 450));

  % iterate dbscan with multiple epsilon until nRMSE of all cluster < 10%
  fprintf('Finding best epsilon for dbscan Kernel Density clustering of a(filt)/fCDOM(filt) ... ')
  k = 1;
  epsilons = 1:-0.01:0.01;
  keep_going = true;
  max_nb_clusters = day(datetime(max(filt_interp.dt) - min(filt_interp.dt), 'ConvertFrom', 'datenum'));
  eval_cluster = table();
  eval_cluster.nRMSE_a450 = NaN(length(epsilons), max_nb_clusters);
  eval_cluster.nRMSE_c450 = NaN(length(epsilons), max_nb_clusters);
  eval_cluster.slope_a450 = NaN(length(epsilons), max_nb_clusters);
  eval_cluster.slope_c450 = NaN(length(epsilons), max_nb_clusters);
  epsilon_a = epsilons(1);
  epsilon_c = epsilons(1);
  % normalize the time variable between 0 to 0.5
  time_var = (filt.dt-min(filt.dt))/max(filt.dt-min(filt.dt)) * time_weight_for_cluster;
  data_tocluster_a = [filt.a(:, id450_a)./filt.fdom time_var]; % time_var filt.t filt.s
  data_tocluster_c = [filt.c(:, id450_c)./filt.fdom time_var]; % time_var filt.t filt.s
  if any(~isnan(data_tocluster_a(:, 1))) || any(~isnan(data_tocluster_c(:, 1)))
    % find minimum cluster size depending on the number of days loaded
    % round(size(time_var, 1) / (max(day(datetime(filt.dt-min(filt.dt)+1, 'ConvertFrom', 'datenum'))) * 24))
    while k <= length(epsilons) && keep_going
      clusters.a = dbscan(data_tocluster_a, epsilons(k), min_nb_pts_per_cluster, 'Distance', 'euclidean');
      clusters.a(clusters.a == -1) = NaN;
      clusters.c = dbscan(data_tocluster_c, epsilons(k), min_nb_pts_per_cluster, 'Distance', 'euclidean');
      clusters.c(clusters.c == -1) = NaN;
      % regress a&c with fdom
      [regress_stats_a, regress_stats_c] = regress_acfilt(filt.a(:, id450_a), filt.c(:, id450_c), filt.fdom, clusters);
      eval_cluster.nRMSE_a450(k, 1:size(regress_stats_a.nRMSE, 2)) = regress_stats_a.nRMSE;
      eval_cluster.nRMSE_c450(k, 1:size(regress_stats_c.nRMSE, 2)) = regress_stats_c.nRMSE;
      eval_cluster.slope_a450(k, 1:size(regress_stats_a.slope, 2)) = regress_stats_a.slope;
      eval_cluster.slope_c450(k, 1:size(regress_stats_c.slope, 2)) = regress_stats_c.slope;
      if k >= 2
        % get epsilon_a
        if all(eval_cluster.nRMSE_a450(k,1:size(regress_stats_a.nRMSE, 2)) < 10, 2)
          epsilon_a = epsilons(k);
        end
        % get epsilon_c
        if all(eval_cluster.nRMSE_c450(k,1:size(regress_stats_c.nRMSE, 2)) < 10, 2)
          epsilon_c = epsilons(k);
        end
        % if both a and c found stop iterations
        if all(eval_cluster.nRMSE_a450(k,1:size(regress_stats_a.nRMSE, 2)) < 10, 2) && ...
             all(eval_cluster.nRMSE_c450(k,1:size(regress_stats_c.nRMSE, 2)) < 10, 2)
          keep_going = false;
        elseif size(regress_stats_a.nRMSE, 2) > max_nb_clusters || size(regress_stats_c.nRMSE, 2) > max_nb_clusters
          keep_going = false;
          if epsilon_a == epsilons(1)
            epsilon_a = epsilons(k);
          end
          if epsilon_c == epsilons(1)
            epsilon_c = epsilons(k);
          end
        end
      end
      k = k + 1;
    end
    fprintf('done\n')
    % eval_cluster.nRMSE_a450(:,k:end) = [];
    % eval_cluster.nRMSE_c450(:,k:end) = [];
    % eval_cluster.slope_a450(:,k:end) = [];
    % eval_cluster.slope_c450(:,k:end) = [];
    % eval_cluster(k:end, :) = [];

    % run clustering with the best epsilon
    fprintf('Kernel Density clustering of a(filt)/fCDOM(filt) ... ')
    clusters.a = dbscan(data_tocluster_a, epsilon_a, min_nb_pts_per_cluster, 'Distance', 'euclidean');
    clusters.a(clusters.a == -1) = NaN;
    clusters.c = dbscan(data_tocluster_c, epsilon_c, min_nb_pts_per_cluster, 'Distance', 'euclidean');
    clusters.c(clusters.c == -1) = NaN;
    fprintf('done\n')
  
    clusters.a = fillmissing(clusters.a, 'nearest', 'SamplePoints', filt.dt);
    clusters.c = fillmissing(clusters.c, 'nearest', 'SamplePoints', filt.dt);

    % regress a&c with fdom
    regress_stats = struct();
    [regress_stats.a, regress_stats.c] = regress_acfilt(filt.a, filt.c, filt.fdom, clusters);

    % plot clusters
    fprintf('Plotting clustering results ... ')
    figure(33); subplot(2, 2, 1);
    gscatter(datetime(filt.dt, 'ConvertFrom', 'datenum'), filt.a(:, id450_a) ./ filt.fdom, ...
      clusters.a, [], 'vo<s>pdh+*x', 8, 'on', 'Time', 'a/fdom');
    title('a cluster time series'); set(gca, 'Fontsize', 12)
    subplot(2, 2, 2); gsc = gscatter(filt.fdom, filt.a(:,id450_a), clusters.a, [], 'vo<s>pdh+*x');
    xlimit = get(gca, 'XLim');
    ylimit = get(gca, 'YLim');
    for c = 1:size(regress_stats.a.slope, 2)
      rl = refline(regress_stats.a.slope(id450_a, c), regress_stats.a.intercept(id450_a, c));
      set(rl, 'Color', gsc(c).Color, 'LineWidth', 1)
    end
    set(gca, 'XLim', xlimit)
    set(gca, 'YLim', ylimit)
    xlabel('fdom'); ylabel('a'); title('a cluster: a/fdom'); set(gca, 'Fontsize', 12)
    leg = findobj(gcf, 'Type', 'Legend');
    title(leg,'Clusters')
  
    subplot(2, 2, 3);
    gscatter(datetime(filt.dt, 'ConvertFrom', 'datenum'), filt.c(:, id450_c) ./ filt.fdom, ...
      clusters.c, [], 'vo<s>pdh+*x', 8, 'on', 'Time', 'c/fdom');
    title('c cluster time series'); set(gca, 'Fontsize', 12)
    subplot(2, 2, 4); gsc = gscatter(filt.fdom, filt.c(:,id450_c), clusters.c, [], 'vo<s>pdh+*x');
    xlimit = get(gca, 'XLim');
    ylimit = get(gca, 'YLim');
    for c = 1:size(regress_stats.c.slope, 2)
      rl = refline(regress_stats.c.slope(id450_c, c), regress_stats.c.intercept(id450_c, c));
      set(rl, 'Color', gsc(c).Color, 'LineWidth', 1)
    end
    set(gca, 'XLim', xlimit)
    set(gca, 'YLim', ylimit)
    xlabel('fdom'); ylabel('c'); title('c cluster: c/fdom'); set(gca, 'Fontsize', 12)
    leg = findobj(gcf, 'Type', 'Legend');
    leg(1).String = strrep(leg(1).String, 'data', 'regression cluster ');
    leg(3).String = strrep(leg(3).String, 'data', 'regression cluster ');
    title(leg,'Clusters')
    fprintf('done\n')
  
    drawnow
  
    % % find wavelength where fdom fit and temperature fit work poorly
    % bad_fit_a = all(regress_stats.a.R2 < prctile(regress_stats.a.R2, 5, 1), 2);
    % bad_fit_c = all(regress_stats.c.R2 < prctile(regress_stats.c.R2, 5, 1), 2);
    
    %%% apply the regression to get a and c dissolved from fdom
    groups_a = unique(clusters.a(~isnan(clusters.a)));
    groups_c = unique(clusters.c(~isnan(clusters.c)));
  
    % Weight cluster for each filter event in case clusters are changing
    % during a filter event (unlikely with the new dbscan method but let's keep it in case)
    filt_avg = addvars(filt_avg, NaN(size(filt_avg,1),size(groups_a,1)), NaN(size(filt_avg,1),size(groups_c,1)), ...
      'NewVariableNames', {'cluster_weighted_a', 'cluster_weighted_c'}, 'after', 'end');
    filt_avg = addvars(filt_avg, false(size(filt_avg.dt)), false(size(filt_avg.dt)), ...
      'NewVariableNames', {'cluster_a_flag', 'cluster_c_flag'}, 'after', 'end');
    for i=1:size(filt_avg, 1)
      sel_filt = filt_avg.start(i) <= filt.dt & filt.dt <= filt_avg.end(i);
      foo = clusters(sel_filt,:);
      % flag if clusters are not all the same within a single filter event
      if size(unique(foo.a(~isnan(foo.a))), 1) > 1
        filt_avg.cluster_a_flag(i) = true;
      end
      if size(unique(foo.c(~isnan(foo.c))), 1) > 1
        filt_avg.cluster_c_flag(i) = true;
      end
      % weight clusters proportionally for each filter events
      for j = 1:size(groups_a,1)
        filt_avg.cluster_weighted_a(i,j) = sum(foo.a == groups_a(j)) ./ sum(~isnan(foo.a));
      end
      for j = 1:size(groups_c,1)
        filt_avg.cluster_weighted_c(i,j) = sum(foo.c == groups_c(j)) ./ sum(~isnan(foo.c));
      end
    end
    % weight a slope and intercept
    filt_avg.slope_a = zeros(size(filt_avg, 1), size(regress_stats.a, 1));
    filt_avg.intercept_a = zeros(size(filt_avg, 1), size(regress_stats.a, 1));
    for j = 1:size(groups_a, 1)
      filt_avg.slope_a = filt_avg.slope_a + regress_stats.a.slope(:,j)'.*filt_avg.cluster_weighted_a(:,j);
      filt_avg.intercept_a = filt_avg.intercept_a + regress_stats.a.intercept(:,j)'.*filt_avg.cluster_weighted_a(:,j);
    end
    % weight c slope and intercept
    filt_avg.slope_c = zeros(size(filt_avg, 1), size(regress_stats.c, 1));
    filt_avg.intercept_c = zeros(size(filt_avg, 1), size(regress_stats.c, 1));
    for j = 1:size(groups_c, 1)
      filt_avg.slope_c = filt_avg.slope_c + regress_stats.c.slope(:,j)'.*filt_avg.cluster_weighted_c(:,j);
      filt_avg.intercept_c = filt_avg.intercept_c + regress_stats.c.intercept(:,j)'.*filt_avg.cluster_weighted_c(:,j);
    end
    % interpolate slope and intercept cluster specific onto filt_interp.dt
    filt_interp.slope_interp_a = interp_extrap(filt_avg, filt_interp.dt, 'slope_a', [], true, 'linear', 'nearest');
    filt_interp.intercept_interp_a = interp_extrap(filt_avg, filt_interp.dt, 'intercept_a', [],  true, 'linear', 'nearest');
    filt_interp.slope_interp_c = interp_extrap(filt_avg, filt_interp.dt, 'slope_c', [],  true, 'linear', 'nearest');
    filt_interp.intercept_interp_c = interp_extrap(filt_avg, filt_interp.dt, 'intercept_c', [],  true, 'linear', 'nearest');
  else
    % clusters.a = ones(size(data_tocluster_a));
    % clusters.c = ones(size(data_tocluster_c));
    % regress_stats = [];
    if all(isnan(filt.fdom))
      error('All fDOM filter event data loaded are NaNs')
    end
    if all(isnan(data_tocluster_a))
      error('All "a" filter event data loaded are NaNs')
    end
    if all(isnan(data_tocluster_c))
      error('All "c" filter event data loaded are NaNs')
    end
  end

end


%% regression between a/c and other variable in filter events
function [a_reg, c_reg] = regress_acfilt(a, c, ancillary, clusters)
  if nargin < 4
    clusters = table();
    clusters.a = ones(size(a, 1), 1);
    clusters.c = ones(size(c, 1), 1);
  end
  % robust linear regression between filt a&c and ancillary variable
  a_reg = array2table(NaN(size(a, 2), 5), 'VariableNames', {'slope', 'intercept', 'RMSE', 'nRMSE', 'R2'});
  c_reg = array2table(NaN(size(c, 2), 5), 'VariableNames', {'slope', 'intercept', 'RMSE', 'nRMSE', 'R2'});
  groups_a = unique(clusters.a(~isnan(clusters.a)));
  groups_c = unique(clusters.c(~isnan(clusters.c)));
  for i = 1:size(a, 2)
    for j = 1:size(groups_a, 1)
      id_grp_a = clusters.a == groups_a(j);
      % regress variable with filtered water absorption
      [stats.b, stats.stats] = robustfit(ancillary(id_grp_a), a(id_grp_a, i));
      [~, MSGID] = lastwarn();
      warning('off', MSGID)
      a_reg.slope(i, j) = stats.b(2);
      a_reg.intercept(i, j) = stats.b(1);
      a_reg.RMSE(i, j) = stats.stats.robust_s;
      % calculated normalize RMSE by average and R2
      id_nonan = all(~isnan([a(:, i) ancillary]), 2);
      a_reg.nRMSE(i, j) = a_reg.RMSE(i, j)/abs(mean(a(id_grp_a & id_nonan, i), 'omitnan'))*100;
      a_reg.R2(i, j) = corr(a(id_grp_a & id_nonan, i), stats.b(1) + stats.b(2) * ancillary(id_grp_a & id_nonan))^2;
    end
    for j = 1:size(groups_c, 1)
      id_grp_c = clusters.c == groups_c(j);
      % regress variable with filtered water attenuation
      [stats.b, stats.stats] = robustfit(ancillary(id_grp_c), c(id_grp_c, i));
      [~, MSGID] = lastwarn();
      warning('off', MSGID)
      c_reg.slope(i, j) = stats.b(2);
      c_reg.intercept(i, j) = stats.b(1);
      c_reg.RMSE(i, j) = stats.stats.robust_s;
      % calculated normalize RMSE by average and R2
      id_nonan = all(~isnan([c(:, i) ancillary]), 2);
      c_reg.nRMSE(i, j) = c_reg.RMSE(i, j)/abs(mean(c(id_grp_c & id_nonan, i), 'omitnan'))*100;
      c_reg.R2(i, j) = corr(c(id_grp_c & id_nonan, i), stats.b(1) + stats.b(2) * ancillary(id_grp_c & id_nonan))^2;
    end
  end
end

%%
function data_out = interp_extrap(data_in, dt_vector, var_tointerp, max_missing_length, extrap_bool, interp_method, extrap_method)
  if nargin < 4
    max_missing_length = 30;
    extrap_bool = true;
    interp_method = 'linear';
    extrap_method = 'nearest';
  elseif nargin < 5
    extrap_bool = true;
    interp_method = 'linear';
    extrap_method = 'nearest';
  elseif nargin < 6
    interp_method = 'linear';
    extrap_method = 'nearest';
  elseif nargin < 7
    extrap_method = 'nearest';
  end
  % remove row full of NaNs
  data_in(all(isnan(data_in.(var_tointerp)), 2), :) = [];
  % convert dt_vector to datetime
  datetime_vector = datetime(dt_vector, 'ConvertFrom', 'datenum');
  dt = (datenum(min(datetime_vector):median(diff(datetime_vector)):max(datetime_vector)))';
  % id extrapolation
  extrapolated_id = isnan(interp1(data_in.dt, data_in.(var_tointerp), dt, 'nearest'));
  % find consecutive NaN longer than max_missing_length
  if ~isempty(max_missing_length) && any(any(isnan(data_in.(var_tointerp)), 2))
    missing_data = ~ismember(dt, data_in.dt);
    t = [true; diff(missing_data) ~= 0];
    k = diff(find([t; true])) .* missing_data(t);
    long_nan = k(cumsum(t)) > max_missing_length;
  end
  % interpolate
  data_out = interp1(data_in.dt, data_in.(var_tointerp), dt, interp_method);
  % replace extrapolated data by NaN
  data_out(extrapolated_id) = NaN;
  % fill missing data with extrapolation method
  if extrap_bool
    data_out = fillmissing(data_out, extrap_method, 'SamplePoints', dt);
  end
  % replace interpolated values over gaps > max_missing_length by NaN
  if ~isempty(max_missing_length) && any(any(isnan(data_in.(var_tointerp)), 2))
    if any(long_nan)
      data_out(long_nan, :) = NaN;
    end
  end
  data_out = interp1(dt, data_out, dt_vector, 'nearest');
end

%% moved to separate function to use it elsewhere
% function data_out = round_timestamp(data_in)
%   data_out = table();
%   % make sure data_in.dt in rounded to the time binning frequency
%   datetime_data_in_dt = datetime(data_in.dt, 'ConvertFrom', 'datenum');
%   % get time binning frequency
%   Tbin_data_in = median(diff(datetime_data_in_dt));
%   if Tbin_data_in >= hours(1)
%     data_out.dt = datenum(dateshift(datetime_data_in_dt, 'start', 'hours'));
%   elseif Tbin_data_in >= minutes(1)
%     % round start/end time to minute
%     data_out.dt = datenum(dateshift(datetime_data_in_dt, 'start', 'minutes'));
%   elseif Tbin_data_in >= seconds(1)
%     % round start/end time to seconds
%     data_out.dt = datenum(dateshift(datetime_data_in_dt, 'start', 'seconds'));
%   else
%     error('automatic detection of sampling rate detected a frequency not supported: check round_timestamp function in processACS.m')
%   end
%   % remove duplicates
%   [~, L, ~] = unique(data_in.dt,'first');
%   indexToDump = not(ismember(1:numel(data_in.dt), L));
%   if sum(indexToDump) > 0
%     data_in(indexToDump, :) = [];
%   end
%   % remove duplicates
%   [~, L, ~] = unique(data_out.dt,'first');
%   indexToDump = not(ismember(1:numel(data_out.dt), L));
%   if sum(indexToDump) > 0
%     data_out(indexToDump, :) = [];
%   end
%   % interpolate data on rounded datetime
%   vars = data_in.Properties.VariableNames;
%   vars(strcmp(vars, 'dt')) = [];
%   for v = vars
%     data_out.(v{:}) = interp1(data_in.dt, data_in.(v{:}), data_out.dt, 'linear');
%   end
% end

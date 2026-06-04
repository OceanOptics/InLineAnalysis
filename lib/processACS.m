function [p, g, bad, DIW_biofouling_correction] = processACS(lambda, tot, filt, param, modelG50, modelmphi, di, ...
  CDOM, FTH, fth_constants, interpolation_method, di_method, scattering_correction, compute_ad_aphi, TSG, days2run)
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
  if ~any(strcmp(scattering_correction, {'Zaneveld1994_proportional','Rottgers2013_semiempirical','Semiempirical_blended1','Semiempirical_blended2','Semiempirical_blended3'}))
    warning('only scattering correction supported: %s', strjoin({'Zaneveld1994_proportional','Rottgers2013_semiempirical','Semiempirical_blended1','Semiempirical_blended2','Semiempirical_blended3'}, ', '))
    error('%s residual temperature and scattering correction not supported', scattering_correction)
  end
  flow_data = FTH.qc.tsw;
  
  if ~isempty(tot)
    if ~isdatetime(tot.dt)
      tot.dt = datetime(tot.dt, 'ConvertFrom', 'datenum');
    end
  else
    error('No AC total data loaded')
  end
  if ~isempty(filt)
    if ~isdatetime(filt.dt)
      filt.dt = datetime(filt.dt, 'ConvertFrom', 'datenum');
    end
  else
    error('No AC filtered data loaded')
  end
  if ~isempty(di)
    if ~isdatetime(di.dt)
      di.dt = datetime(di.dt, 'ConvertFrom', 'datenum');
    end
  end

  % round time stamp and remove time duplicates
  tot = round_timestamp(tot, minutes(1));
  filt = round_timestamp(filt, minutes(1));
  flow_data = round_timestamp(flow_data, minutes(1));

  % check if fCDOM data loaded if fCDOM interpolation
  cdom_base = [];
  if ~isempty(CDOM)
    cdom_tbl_name = fieldnames(CDOM.prod);
    if isempty(cdom_tbl_name)
      error('No fDOM prod data loaded')
    end
    % Require both CDOM & Switch position
    if ~isempty(CDOM.prod.(cdom_tbl_name{1}))
      cdom_base = CDOM.prod.(cdom_tbl_name{1});
      if ~isdatetime(cdom_base.dt)
        cdom_base.dt = datetime(cdom_base.dt, 'ConvertFrom', 'datenum');
      end
      if ~any(cdom_base.dt >= min([tot.dt; filt.dt]) & cdom_base.dt <= max([tot.dt; filt.dt]))
        warning('fCDOM dates do not correspond to ACS dates: interpolation switched to "linear"')
        interpolation_method = 'linear';
      else
        % round time stamp and remove time duplicates
        cdom_base = round_timestamp(cdom_base, minutes(1));
        % smooth fDOM if fDOM sensor is WSCD
        if strcmp(CDOM.model, 'WSCD') || strcmp(CDOM.prefix, 'WSCD')
          fdom_temp = table();
          foo_dt = [min([tot.dt; filt.dt]) max([tot.dt; filt.dt])];
          fdom_temp.dt = dateshift((foo_dt(1):minutes(1):foo_dt(2))','start','minute');
          id_nan = isnan(cdom_base.fdom);
          % interpolate fdom on fdom_temp
          fdom_temp.fdom = interp1(cdom_base.dt, cdom_base.fdom, fdom_temp.dt, 'linear');
          fdom_temp.fdom = fillmissing(fdom_temp.fdom, 'nearest');
          % smooth WSCD cdom data removing frequencies higher than 15min
          d = designfilt('lowpassiir', 'FilterOrder', 1, 'PassbandFrequency', 1/(60*15), 'SampleRate', 1/60);
          fdom_temp.fdom = filtfilt(d, fdom_temp.fdom);
          cdom_base.fdom = interp1(fdom_temp.dt, fdom_temp.fdom, cdom_base.dt, 'linear');
          cdom_base.fdom(id_nan) = NaN;
        end
      end
    else
      warning('fCDOM data not loaded: interpolation switched to "linear"')
      interpolation_method = 'linear';
    end
    % if ~isfield(param, 'min_nb_pts_per_cluster')
    %   fprintf('Minimum number of points per cluster set to 200\n')
    %   param.min_nb_pts_per_cluster = 200;
    % end
  end

  % check if TSG data loaded
  if ~isempty(TSG)
    tsg_tbl_name = fieldnames(TSG.prod);
    if ~isempty(tsg_tbl_name)
      if ~isempty(TSG.prod.(tsg_tbl_name{1}))
        tsg_loaded = true;
        tsg_data = TSG.prod.(tsg_tbl_name{1});
      else
        tsg_loaded = false;
      end
    else
      tsg_loaded = false;
    end
    if ~tsg_loaded
      if ~isempty(TSG.qc.tsw)
        tsg_loaded = true;
        tsg_data = TSG.qc.tsw;
      else
        tsg_loaded = false;
      end
    end
    if tsg_loaded
      if ~isdatetime(tsg_data.dt)
        tsg_data.dt = datetime(tsg_data.dt, 'ConvertFrom', 'datenum');
      end
      if ~any(tsg_data.dt >= min([tot.dt; filt.dt]) & tsg_data.dt <= max([tot.dt; filt.dt]))
        warning(['TSG dates do not correspond to ACS dates: T/S correction not applied before removing dissolved part from total.' newline ...
          'T/S assumed to be constant between filter and total events.' newline ...
          'Remaining T variation will be corrected using the residual temperature correction\n'])
        tsg_loaded = false;
      else
        tsg_loaded = true;
        % round time stamp and remove time duplicates
        tsg_data = round_timestamp(tsg_data, minutes(1));
      end
    end
  else
    tsg_loaded = false;
  end
  if ~tsg_loaded
    warning('No TSG qc or prod data loaded')
  end

  % interpolate flow_data.swt onto binned data to fill missing flow data
  fth_interp_dt = (dateshift(min([tot.dt; flow_data.dt; filt.dt]), 'Start', 'minute'):minutes(1):...
    dateshift(max([tot.dt; flow_data.dt; filt.dt]), 'Start', 'minute')+minutes(1))';
  fth_interp = table(fth_interp_dt, 'VariableNames', {'dt'});
  % fth_interp = table([tot.dt; flow_data.dt; filt.dt], 'VariableNames', {'dt'});
  fth_interp = round_timestamp(fth_interp, minutes(1));
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

  %%%%%%%%%%%%%%%%%%%%%%% Compute filtered average %%%%%%%%%%%%%%%%%%%%%%%
  % TODO maybe: write a filter_average function that can be use with every instrument

  % create filter average data table
  filt_avg = table(NaT(size(sel_start)), 'VariableNames', {'dt'});
  filt_avg.start = fth_interp.dt(sel_start);
  filt_avg.end = fth_interp.dt(sel_end);
  filt_avg.a = NaN(size(filt_avg,1), size(lambda.a, 2));
  filt_avg.c = NaN(size(filt_avg,1), size(lambda.c, 2));
  filt_avg.a_avg_sd = NaN(size(filt_avg,1), size(lambda.a, 2));
  filt_avg.c_avg_sd = NaN(size(filt_avg,1), size(lambda.c, 2));
  filt_avg.a_avg_n = NaN(size(filt_avg,1), 1);
  filt_avg.c_avg_n = NaN(size(filt_avg,1), 1);
  % add flag for T/S correction (default true: data NOT T/S corrected)
  filt_avg.flag_deltaTSnotcorrected_filt = true(size(filt_avg,1), 1);

  % create filter interpolation table
  filt_interp = table(tot.dt, 'VariableNames', {'dt'});
  filt_interp = merge_timeseries(filt_interp, flow_data, FTH.view.spd_variable, '', 5, minutes(1));
  filt_interp = renamevars(filt_interp, FTH.view.spd_variable, 'flow_rate');
  % filt_interp.flow_rate = interp1(flow_data.dt, flow_data.(FTH.view.swt_variable), filt_interp.dt, 'previous');%, 'linear', 'extrap');
  filt_interp.a = NaN(size(tot.a));
  filt_interp.c = NaN(size(tot.c));
  filt_interp.flag_deltaTSnotcorrected_total = true(size(tot.dt)); % by default set flag to true until data are corrected
  filt_interp.flag_deltaTSnotcorrected_filt0 = true(size(tot.dt)); % by default set flag to true until data are corrected
  filt_interp.flag_deltaTSnotcorrected_filt1 = true(size(tot.dt)); % by default set flag to true until data are corrected

  % add flag for T/S correction (default true: data NOT T/S corrected)
  filt.flag_deltaTSnotcorrected_filt = true(size(filt.dt));
  % add T/S to filt and filt_interp if TSG loaded
  if tsg_loaded
    % interpolate linear T and S on filt_interp table
    filt_interp = addvars(filt_interp, interp_extrap(tsg_data, filt_interp.dt, TSG.temperature_variable, 30, true, 'linear', 'nearest'), ...
      interp_extrap(tsg_data, filt_interp.dt, TSG.salinity_variable, 30, true, 'linear', 'nearest'), 'NewVariableNames', {'t','s'}, 'After', 'dt');
    % interpolate T and S on filtered data
    filt = addvars(filt, interp_extrap(tsg_data, filt.dt, TSG.temperature_variable, 30, true, 'linear', 'nearest'), ...
      interp_extrap(tsg_data, filt.dt, TSG.salinity_variable, 30, true, 'linear', 'nearest'), 'NewVariableNames', {'t','s'}, 'After', 'dt');
    % add T and S variables on filt_avg table
    filt_avg = addvars(filt_avg, NaN(size(filt_avg,1),1), NaN(size(filt_avg,1),1), 'NewVariableNames', {'t','s'}, 'After', 'dt');
    % remove psiS and psiT from tot and filt with Tref = median(T) of whole processed section
    Tref = median([filt.t; filt_interp.t], 'omitnan');
    a_psiS = interp1(psi.wl, psi.a_psiS, lambda.a, 'spline');
    c_psiS = interp1(psi.wl, psi.c_psiS, lambda.c, 'spline');
    a_psiT = interp1(psi.wl, psi.psiT, lambda.a, 'spline');
    c_psiT = interp1(psi.wl, psi.psiT, lambda.c, 'spline');
    % Correct independently only when t and s are available to avoid getting NaN when t or s are NaN
    % correct T/S in filt
    filt_t_nan = isnan(filt.t);
    filt_s_nan = isnan(filt.s);
    filt_ts_nan = isnan(filt.t) | isnan(filt.s);
    filt.a(~filt_t_nan, :) = filt.a(~filt_t_nan, :) - (a_psiT .* (filt.t(~filt_t_nan) - Tref));
    filt.a(~filt_s_nan, :) = filt.a(~filt_s_nan, :) - (a_psiS .* filt.s(~filt_ts_nan));
    % filt.a(~filt_ts_nan, :) = filt.a(~filt_ts_nan, :) - (a_psiT .* (filt.t(~filt_ts_nan) - Tref)) - (a_psiS .* filt.s(~filt_ts_nan));
    filt.c(~filt_t_nan, :) = filt.c(~filt_t_nan, :) - (c_psiT .* (filt.t(~filt_t_nan) - Tref));
    filt.c(~filt_s_nan, :) = filt.c(~filt_s_nan, :) - (c_psiS .* filt.s(~filt_s_nan));
    % filt.c(~filt_ts_nan, :) = filt.c(~filt_ts_nan, :) - (c_psiT .* (filt.t(~filt_ts_nan) - Tref)) - (c_psiS .* filt.s(~filt_ts_nan));
    % correct T/S in tot
    tot_t_nan = isnan(filt_interp.t);
    tot_s_nan = isnan(filt_interp.s);
    tot_ts_nan = isnan(filt_interp.t) | isnan(filt_interp.s);
    tot.a(~tot_t_nan, :) = tot.a(~tot_t_nan, :) - (a_psiT .* (filt_interp.t(~tot_t_nan) - Tref));
    tot.a(~tot_s_nan, :) = tot.a(~tot_s_nan, :) - (a_psiS .* filt_interp.s(~tot_s_nan));
    % tot.a(~tot_ts_nan, :) = tot.a(~tot_ts_nan, :) - (a_psiT .* (filt_interp.t(~tot_ts_nan) - Tref)) - (a_psiS .* filt_interp.s(~tot_ts_nan));
    tot.c(~tot_t_nan, :) = tot.c(~tot_t_nan, :) - (c_psiT .* (filt_interp.t(~tot_t_nan) - Tref));
    tot.c(~tot_s_nan, :) = tot.c(~tot_s_nan, :) - (c_psiS .* filt_interp.s(~tot_s_nan));
    % tot.c(~tot_ts_nan, :) = tot.c(~tot_ts_nan, :) - (c_psiT .* (filt_interp.t(~tot_ts_nan) - Tref)) - (c_psiS .* filt_interp.s(~tot_ts_nan));
    % adjust T/S correction flags
    % filt.flag_deltaTnotcorrected_filt(~filt_t_nan) = false; % change flag if corrected
    % filt.flag_deltaSnotcorrected_filt(~filt_s_nan) = false; % change flag if corrected
    filt.flag_deltaTSnotcorrected_filt(~filt_ts_nan) = false; % change flag if corrected
    % filt_interp.flag_deltaTnotcorrected_total(~tot_t_nan) = false; % change flag if corrected
    % filt_interp.flag_deltaSnotcorrected_total(~tot_s_nan) = false; % change flag if corrected
    filt_interp.flag_deltaTSnotcorrected_total(~tot_ts_nan) = false; % change flag if corrected
  end

  % prepare fCDOM data if available
  if ~isempty(cdom_base)
    % interpolate spline and extrapolate nearest fdom on filt_interp table
    filt_interp = addvars(filt_interp, interp_extrap(cdom_base, filt_interp.dt, 'fdom', 30, false, 'linear'), ...
      'NewVariableNames', 'fdom', 'After', 'dt');
    filt_interp = addvars(filt_interp, interp_extrap(cdom_base, filt_interp.dt, 'fdom_sd', 30, false, 'linear'), ...
      'NewVariableNames', 'fdom_sd', 'After', 'fdom');
    filt_interp = addvars(filt_interp, interp_extrap(cdom_base, filt_interp.dt, 'fdom_n', 30, false, 'nearest'), ...
      'NewVariableNames', 'fdom_n', 'After', 'fdom_sd');
    % interpolate fdom on filtered data
    filt = addvars(filt, interp_extrap(cdom_base, filt.dt, 'fdom', 30, false, 'linear'), ...
      'NewVariableNames', 'fdom', 'After', 'dt');
    filt = addvars(filt, interp_extrap(cdom_base, filt.dt, 'fdom_sd', 30, false, 'linear'), ...
      'NewVariableNames', 'fdom_sd', 'After', 'fdom');
    filt = addvars(filt, interp_extrap(cdom_base, filt.dt, 'fdom_n', 30, false, 'linear'), ...
      'NewVariableNames', 'fdom_n', 'After', 'fdom_sd');
    % add fdom variables on filt_avg table
    filt_avg = addvars(filt_avg, NaN(size(filt_avg,1),1), 'NewVariableNames', 'fdom', 'After', 'dt');
  end

  % compute filter averages
  for i=1:size(sel_start, 1) % 37
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
      filt_avg.flag_deltaTSnotcorrected_filt(i) = foo.flag_deltaTSnotcorrected_filt;
      % repeat for T/S if tsg is loaded
      if tsg_loaded
        filt_avg.t(i) = foo.t;
        filt_avg.s(i) = foo.s;
      end
      % repeat for fdom available
      if any(strcmp(foo.Properties.VariableNames, 'fdom'))
        filt_avg.fdom(i) = foo.fdom;
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
      % change flag only if all filtered data are T/S corrected
      if ~any(foo.flag_deltaTSnotcorrected_filt) 
        filt_avg.flag_deltaTSnotcorrected_filt(i) = false;
      end
      % repeat for T/S if tsg is loaded
      if tsg_loaded
        foo.t(all(a_perc25, 2) & all(c_perc25, 2)) = NaN;
        foo.s(all(a_perc25, 2) & all(c_perc25, 2)) = NaN;
        % compute average of all values smaller than 25th percentile for each filter event
        filt_avg.t(i) = mean(foo.t, 1, 'omitnan');
        filt_avg.s(i) = mean(foo.s, 1, 'omitnan');
      end
      % repeat for fdom if available
      if any(strcmp(foo.Properties.VariableNames, 'fdom'))
        foo.fdom(foo.fdom > prctile(foo.fdom, 25, 1)) = NaN;
        % compute average of all values smaller than 25th percentile for each filter event
        filt_avg.fdom(i) = mean(foo.fdom, 1, 'omitnan');
      end
    end
  end
  % remove empty filter events
  rm_filter_event = all(isnan(filt_avg.a), 2) & all(isnan(filt_avg.c), 2);
  filt_avg(rm_filter_event, :) = [];

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
  
  switch interpolation_method
    case 'CDOM'
      % % interpolate ag and cg using fcdom
      % cluster_param = table();
      % cluster_param.time_weight = [0; 0.001; 0.01; 0.1; 1; 5; 10];
      % cluster_param.perc_negative_slope = NaN(size(cluster_param.time_weight));
      % cluster_param.perc_not_clustered = NaN(size(cluster_param.time_weight));
      % cluster_param.filt_interp = cell(size(cluster_param.time_weight));
      % cluster_param.filt = cell(size(cluster_param.time_weight));
      % cluster_param.regress_stats = cell(size(cluster_param.time_weight));
      % % get correlation between ag/cg and fdom using one cluster
      % [cluster_param.filt_interp{1}, cluster_param.filt{1}, cluster_param.regress_stats{1}] = agcg_fdom_interpolation(...
      %   tot, filt_interp, filt, lambda, filt_avg, param.min_nb_pts_per_cluster, cluster_param.time_weight(1), false);
      % cluster_param.perc_negative_slope(1) = (sum(cluster_param.filt_interp{1}.flag_a_negative_slope(:)) + ...
      %   sum(cluster_param.filt_interp{1}.flag_c_negative_slope(:))) / (size(filt_interp, 1)*2);
      % cluster_param.perc_not_clustered(1) = (sum(cluster_param.filt_interp{1}.flag_a_filt0_not_clustered) + ...
      %   sum(cluster_param.filt_interp{1}.flag_c_filt0_not_clustered)) / (size(filt_interp, 1)*2);
      % % Try multiple weights for time variable used for clustering and keep the one with the least perc_negative_slope,
      % % the least perc_not_clustered, and the smallest time_weight
      % if cluster_param.perc_negative_slope(1) > 0
      %   for i = 2:size(cluster_param,1)
      %     [cluster_param.filt_interp{i}, cluster_param.filt{i}, cluster_param.regress_stats{i}] = agcg_fdom_interpolation(...
      %       tot, filt_interp, filt, lambda, filt_avg, param.min_nb_pts_per_cluster, cluster_param.time_weight(i), true);
      %     cluster_param.perc_negative_slope(i) = (sum(cluster_param.filt_interp{i}.flag_a_negative_slope(:)) + ...
      %       sum(cluster_param.filt_interp{i}.flag_c_negative_slope(:))) / (size(filt_interp, 1)*2);
      %     cluster_param.perc_not_clustered(i) = (sum(cluster_param.filt_interp{i}.flag_a_filt0_not_clustered) + ...
      %       sum(cluster_param.filt_interp{i}.flag_c_filt0_not_clustered)) / (size(filt_interp, 1)*2);
      %   end
      % end
      % 
      % cluster_param = sortrows(cluster_param, {'perc_negative_slope','perc_not_clustered','time_weight'});
      % fprintf('%.3f time weight selected for clustering\n', cluster_param.time_weight(1))
      % filt_interp = cluster_param.filt_interp{1};
      % filt = cluster_param.filt{1};
      % regression_stats = cluster_param.regress_stats{1};
      % 
      % % id 450 nm wl to plot
      % id450_a = abs(lambda.a - 450) == min(abs(lambda.a - 450));
      % id450_c = abs(lambda.c - 450) == min(abs(lambda.c - 450));
      % 
      % % replace nan by 1 when all clusters are NaN
      % if all(isnan(filt.a_clusters))
      %   filt.a_clusters = ones(size(filt.a_clusters));
      % end
      % if all(isnan(filt.c_clusters))
      %   filt.c_clusters = ones(size(filt.c_clusters));
      % end
      % 
      % % plot clusters
      % fprintf('Plotting clustering results ... ')
      % figure(33); subplot(2, 2, 1);
      % gscatter(filt.dt, filt.a(:, id450_a) ./ filt.fdom, filt.a_clusters, [], ...
      %   'vo<s>pdh+*x', 8, 'on', 'Time', 'a/fdom');
      % title('a cluster time series'); set(gca, 'Fontsize', 12)
      % subplot(2, 2, 2); gsc = gscatter(filt.fdom, filt.a(:,id450_a), filt.a_clusters, [], 'vo<s>pdh+*x');
      % xlimit = get(gca, 'XLim');
      % ylimit = get(gca, 'YLim');
      % for c = 1:size(regression_stats.a.slope, 2)
      %   rl = refline(regression_stats.a.slope(id450_a, c), regression_stats.a.intercept(id450_a, c));
      %   set(rl, 'Color', gsc(c).Color, 'LineWidth', 1)
      % end
      % set(gca, 'XLim', xlimit)
      % set(gca, 'YLim', ylimit)
      % xlabel('fdom'); ylabel('a'); title('a cluster: a/fdom'); set(gca, 'Fontsize', 12)
      % leg = findobj(gcf, 'Type', 'Legend');
      % title(leg,'Clusters')
      % 
      % subplot(2, 2, 3);
      % gscatter(filt.dt, filt.c(:, id450_c) ./ filt.fdom, filt.c_clusters, [], ...
      %   'vo<s>pdh+*x', 8, 'on', 'Time', 'c/fdom');
      % title('c cluster time series'); set(gca, 'Fontsize', 12)
      % subplot(2, 2, 4); gsc = gscatter(filt.fdom, filt.c(:,id450_c), filt.c_clusters, [], 'vo<s>pdh+*x');
      % xlimit = get(gca, 'XLim');
      % ylimit = get(gca, 'YLim');
      % for c = 1:size(regression_stats.c.slope, 2)
      %   rl = refline(regression_stats.c.slope(id450_c, c), regression_stats.c.intercept(id450_c, c));
      %   set(rl, 'Color', gsc(c).Color, 'LineWidth', 1)
      % end
      % set(gca, 'XLim', xlimit)
      % set(gca, 'YLim', ylimit)
      % xlabel('fdom'); ylabel('c'); title('c cluster: c/fdom'); set(gca, 'Fontsize', 12)
      % leg = findobj(gcf, 'Type', 'Legend');
      % leg(1).String = strrep(leg(1).String, 'data', 'regression cluster ');
      % leg(3).String = strrep(leg(3).String, 'data', 'regression cluster ');
      % title(leg,'Clusters')
      % fprintf('done\n')
      % drawnow
      
      % get correlation between dag/dcg and dfdom
      [filt_interp, filt, fit_unc] = agcg_fdom_interpolation(filt_interp, filt, lambda, filt_avg);

      % compute filter event standard deviation to be used in final error propagation
      filt_interp.a_avg_sd = interp1(filt.dt, filt.a_avg_sd, filt_interp.dt, 'linear');
      filt_interp.c_avg_sd = interp1(filt.dt, filt.c_avg_sd, filt_interp.dt, 'linear');
      
      % remove interpolation when there is no data
      filt_interp.a_avg_sd(all(isnan(tot.a),2), :) = NaN;
      filt_interp.a(all(isnan(tot.a),2), :) = NaN;
      filt_interp.c_avg_sd(all(isnan(tot.c),2), :) = NaN;
      filt_interp.c(all(isnan(tot.c),2), :) = NaN;
      
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
      % regression_stats = struct();
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
    filt_interp_id = filt_interp.dt >= min(days2run) & filt_interp.dt < max(days2run)+days(1);
    tot_id = tot.dt >= min(days2run) & tot.dt < max(days2run)+days(1);
    filt_avg_id = filt_avg.dt >= min(days2run) & filt_avg.dt < max(days2run)+days(1);
    % plot
    fh = figure(52);
    tiledlayout(2, 1, "TileSpacing", "compact", "Padding", "compact")
    nexttile % nxt = 
    if any(strcmp(filt_avg.Properties.VariableNames, 'fdom'))
      visFlag([], filt_interp(filt_interp_id, :), tot(tot_id, :), [], filt_avg(filt_avg_id, :), [], 'a', ...
        round(size(tot.a, 2)/2), [], [], 'spd', false, fh, filt_interp(filt_interp_id,:), filt_avg(filt_avg_id,:));
    else
      visFlag([], filt_interp(filt_interp_id, :), tot(tot_id, :), [], filt_avg(filt_avg_id, :), [], 'a', round(size(tot.a, 2)/2), [], [], 'spd', false, fh);
    end

    title('Check filter event interpolation, press q to continue', 'FontSize', 14)
    nexttile % nxt = 
    if any(strcmp(filt_avg.Properties.VariableNames, 'fdom'))
      visFlag([], filt_interp(filt_interp_id, :), tot(tot_id, :), [], filt_avg(filt_avg_id, :), [], 'c', ...
        round(size(tot.a, 2)/2), [], [], 'spd', false, fh, filt_interp(filt_interp_id,:), filt_avg(filt_avg_id,:));
    else
      visFlag([], filt_interp(filt_interp_id, :), tot(tot_id, :), [], filt_avg(filt_avg_id, :), [], 'c', round(size(tot.a, 2)/2), [], [], 'spd', false, fh);
    end
    drawnow
    % guiSelectOnTimeSeries(fh);
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

  % % keep apm for intermediate stats
  % pathf = '/Volumes/Extreme SSD/TaraEuropa/prod';
  % filename = sprintf('TaraEuropa_InLine_ACS_%s_%s_apm_%s_interpolation_v20260302.mat', ...
  %   datetime(min(p.dt), 'Format','yyyyMMdd'),datetime(max(p.dt), 'Format','yyyyMMdd'),interpolation_method);
  % save(fullfile(pathf, filename), 'p')





  
  if size(lambda.a, 2) > 50 % perform two separate corrections for ap and cp only for ACS data, not AC9
    % Interpolate wavelengths for Scattering & Residual temperature correction
    ap_for_cpresiduals_corr = interp1(lambda.a', p.ap', lambda.c', 'linear', 'extrap')';
    cp_for_apresiduals_corr = interp1(lambda.c', p.cp', lambda.a', 'linear', 'extrap')';
    % ap Scattering & Residual temperature correction
    % cp Residual correction (for efficiency use the one computed from ap as it should be the same)
    fprintf('ap %s residual temperature and scattering correction ', strrep(scattering_correction, '_', ' '))
    switch scattering_correction
      case 'Rottgers2013_semiempirical'
        [p.ap, ~, filt_interp.flag_Tresidual] = ResidualTempScatterCorrRottgers_semiempirical(p.ap, cp_for_apresiduals_corr, lambda.a, psi, p.dt);
        fprintf(' Done\n')
        fprintf('cp %s residual temperature and scattering correction ', strrep(scattering_correction, '_', ' '))
        [~, p.cp, ~] = ResidualTempScatterCorrRottgers_semiempirical(ap_for_cpresiduals_corr, p.cp, lambda.c, psi, p.dt);
        nap_offset = true;
        fprintf(' Done\n')
      case 'Zaneveld1994_proportional'
        [p.ap, ~, filt_interp.flag_Tresidual] = ResidualTempScatterCorrZaneveld_proportional(p.ap, cp_for_apresiduals_corr, lambda.a, psi, p.dt);
        fprintf(' Done\n')
        fprintf('cp %s residual temperature and scattering correction ', strrep(scattering_correction, '_', ' '))
        [~, p.cp, ~] = ResidualTempScatterCorrZaneveld_proportional(ap_for_cpresiduals_corr, p.cp, lambda.c, psi, p.dt);
        nap_offset = false;
        fprintf(' Done\n')
      case 'Semiempirical_blended1'
        [p.ap, ~, filt_interp.flag_Tresidual] = ResidualTempScatterCorrSemiempirical_blended1(p.ap, cp_for_apresiduals_corr, lambda.a, psi, p.dt);
        fprintf(' Done\n')
        fprintf('cp %s residual temperature and scattering correction ', strrep(scattering_correction, '_', ' '))
        [~, p.cp, ~] = ResidualTempScatterCorrSemiempirical_blended1(ap_for_cpresiduals_corr, p.cp, lambda.c, psi, p.dt);
        nap_offset = true;
        fprintf(' Done\n')
      case 'Semiempirical_blended2'
        [p.ap, ~, filt_interp.flag_Tresidual] = ResidualTempScatterCorrSemiempirical_blended2(p.ap, cp_for_apresiduals_corr, lambda.a, psi, p.dt);
        fprintf(' Done\n')
        fprintf('cp %s residual temperature and scattering correction ', strrep(scattering_correction, '_', ' '))
        [~, p.cp, ~] = ResidualTempScatterCorrSemiempirical_blended2(ap_for_cpresiduals_corr, p.cp, lambda.c, psi, p.dt);
        nap_offset = true;
        fprintf(' Done\n')
      case 'Semiempirical_blended3'
        [p.ap, ~, filt_interp.flag_Tresidual] = ResidualTempScatterCorrSemiempirical_blended3(p.ap, cp_for_apresiduals_corr, lambda.a, psi, p.dt);
        fprintf(' Done\n')
        fprintf('cp %s residual temperature and scattering correction ', strrep(scattering_correction, '_', ' '))
        [~, p.cp, ~] = ResidualTempScatterCorrSemiempirical_blended3(ap_for_cpresiduals_corr, p.cp, lambda.c, psi, p.dt);
        nap_offset = true;
        fprintf(' Done\n')
      otherwise
        error('Residual temperature and scattering correction "%s" not supported', scattering_correction)
    end
  else
    fprintf('ap and cp %s residual temperature and scattering correction ', strrep(scattering_correction, '_', ' '))
    switch scattering_correction
      case 'Rottgers2013_semiempirical'
        [p.ap, p.cp, filt_interp.flag_Tresidual] = ResidualTempScatterCorrRottgers_semiempirical(p.ap, p.cp, lambda.ref, psi, p.dt);
        nap_offset = true;
        fprintf(' Done\n')
      case 'Zaneveld1994_proportional'
        [p.ap, p.cp, filt_interp.flag_Tresidual] = ResidualTempScatterCorrZaneveld_proportional(p.ap, p.cp, lambda.ref, psi, p.dt);
        nap_offset = false;
        fprintf(' Done\n')
      case 'Semiempirical_blended1'
        [p.ap, p.cp, filt_interp.flag_Tresidual] = ResidualTempScatterCorrSemiempirical_blended1(p.ap, p.cp, lambda.ref, psi, p.dt);
        nap_offset = true;
      case 'Semiempirical_blended2'
        [p.ap, p.cp, filt_interp.flag_Tresidual] = ResidualTempScatterCorrSemiempirical_blended2(p.ap, p.cp, lambda.ref, psi, p.dt);
        nap_offset = true;
        fprintf(' Done\n')
      case 'Semiempirical_blended3'
        [p.ap, p.cp, filt_interp.flag_Tresidual] = ResidualTempScatterCorrSemiempirical_blended3(p.ap, p.cp, lambda.ref, psi, p.dt);
        nap_offset = true;
        fprintf(' Done\n')
      otherwise
        error('Residual temperature and scattering correction "%s" not supported', scattering_correction)
    end
    fprintf(' Done\n')
  end
  
  % Remove lines full of NaNs (Rottgers2013_semiempirical potentially fail when temperature correct is too large)
  sel2rm = all(isnan(p.ap), 2) | all(isnan(p.cp), 2);
  p(sel2rm, :) = [];
  tot(sel2rm, :) = [];
  filt_interp(sel2rm, :) = [];

  % % Propagate error (using geometric mean of measurement errors) - OLD METHOD
  % %   Note: Error is not propagated through Scattering & Residual temperature
  % %         correction as required by SeaBASS
  % p.ap_sd = sqrt(tot.a_avg_sd.^2 + filt_interp.a_avg_sd.^2);
  % p.cp_sd = sqrt(tot.c_avg_sd.^2 + filt_interp.c_avg_sd.^2);
  
  % Propagate error - new analytic method (Bourdin et al. 2026)
  p.ap_unc = NaN(size(filt_interp.a));
  p.cp_unc = NaN(size(filt_interp.c));
  dth = hours(filt_avg.dt-min(filt_avg.dt));
  % compute da/dt and dc/dt with dt == delta between filter events in hours
  da_dt = diff(filt_avg.a, [], 1) ./ diff(dth);
  dc_dt = diff(filt_avg.c, [], 1) ./ diff(dth);
  % compute delta time in hours with the closest filter event (and fdom if interpolation_method == CDOM)
  if strcmp(interpolation_method, 'CDOM')
    foo_dt = [filt_avg.dt; filt_interp.dt(~isnan(filt_interp.fdom))];
    % remove duplicates
    [~, L, ~] = unique(foo_dt,'first');
    indexToDump = not(ismember(1:numel(foo_dt), L));
    if any(indexToDump); foo_dt(indexToDump, :) = []; end
    delta_time_filt = hours(abs(filt_interp.dt - interp1(foo_dt, foo_dt, filt_interp.dt, 'nearest')));
  else
    delta_time_filt = hours(abs(filt_interp.dt - interp1(filt_avg.dt, filt_avg.dt, filt_interp.dt, 'nearest')));
  end
  % compute uncertainty of linear interpolation 
  sigma_lin_interp_a = prctile(abs(da_dt), 50, 1) .* delta_time_filt;
  sigma_lin_interp_c = prctile(abs(dc_dt), 50, 1) .* delta_time_filt;

  % Populate uncertainties
  sigma_interp_a = sigma_lin_interp_a;
  sigma_interp_c = sigma_lin_interp_c;
  if strcmp(interpolation_method, 'CDOM')
    sigma_interp_a(~filt_interp.flag_linear_interp_a,:) = filt_interp.ag_interp_unc(~filt_interp.flag_linear_interp_a,:);
    sigma_interp_c(~filt_interp.flag_linear_interp_c,:) = filt_interp.cg_interp_unc(~filt_interp.flag_linear_interp_c,:);
  end

  % divide total std and filtered std by number of independant spectra based on flow rate if available, or assuming a flow rate of 2LPM if no flow is available
  % add up total standard error, filtered standard error, and interpolation standard error
  filt_interp.flow_rate(~isnan(filt_interp.flow_rate)) = 2;
  acs_tube_volume = 0.05; % 0.05 L
  p.ap_unc = sqrt((tot.a_avg_sd./sqrt(filt_interp.flow_rate./acs_tube_volume)).^2 + ...
    (filt_interp.a_avg_sd./sqrt(filt_interp.flow_rate./acs_tube_volume)).^2 + sigma_interp_a.^2);
  p.cp_unc = sqrt((tot.c_avg_sd./sqrt(filt_interp.flow_rate./acs_tube_volume)).^2 + ...
    (filt_interp.c_avg_sd./sqrt(filt_interp.flow_rate./acs_tube_volume)).^2 + sigma_interp_c.^2);


  % sela = lambda.a < 550;
  % % selc = lambda.c < 550;
  % 
  % % compute delta time in hours with the closest filter event (and fdom if interpolation_method == CDOM)
  % delta_time_filt = hours(abs(filt_interp.dt - interp1(filt_avg.dt, filt_avg.dt, filt_interp.dt, 'nearest')));
  % % compute uncertainty of linear interpolation 
  % sigma_lin_interp_a = abs(prctile(da_dt, 95, 1)) .* delta_time_filt;
  % sigma_lin_interp_c = abs(prctile(dc_dt, 95, 1)) .* delta_time_filt;
  % % Populate uncertainties
  % sigma_interp_a = sigma_lin_interp_a;
  % sigma_interp_c = sigma_lin_interp_c;
  % id_flow = ~isnan(filt_interp.flow_rate);
  % % divide total std and filtered std by number of independant spectra based on flow rate if available, or assuming a flow rate of 2LPM if no flow is available
  % % add up total standard error, filtered standard error, and linear interpolation standard error
  % p.ap_unc(id_flow, :) = sqrt((tot.a_avg_sd(id_flow, :)./sqrt(filt_interp.flow_rate(id_flow)./0.05)).^2 + ...
  %   (filt_interp.a_avg_sd(id_flow, :)./sqrt(filt_interp.flow_rate(id_flow)./0.05)).^2 + sigma_interp_a(id_flow, :).^2); % 0.05 L == volume of ACS tubes
  % p.ap_unc(~id_flow, :) = sqrt((tot.a_avg_sd(~id_flow, :)./sqrt(2./0.05)).^2 + ...
  %   (filt_interp.a_avg_sd(~id_flow, :)./sqrt(2./0.05)).^2 + sigma_interp_a(~id_flow, :).^2); % 0.05 L == volume of ACS tubes
  % p.cp_unc(id_flow, :) = sqrt((tot.c_avg_sd(id_flow, :)./sqrt(filt_interp.flow_rate(id_flow)./0.05)).^2 + ...
  %   (filt_interp.c_avg_sd(id_flow, :)./sqrt(filt_interp.flow_rate(id_flow)./0.05)).^2 + sigma_interp_c(id_flow, :).^2); % 0.05 L == volume of ACS tubes
  % p.cp_unc(~id_flow, :) = sqrt((tot.a_avg_sd(~id_flow, :)./sqrt(2./0.05)).^2 + ...
  %   (filt_interp.a_avg_sd(~id_flow, :)./sqrt(2./0.05)).^2 + sigma_interp_a(~id_flow, :).^2); % 0.05 L == volume of ACS tubes
  % 
  % sc_r = [];
  % figure(12); clf; subplot(2,2,1); hold on
  % % perr = prctile(p.ap_unc./abs(p.ap).*100, 75) - median(p.ap_unc./abs(p.ap), 'omitmissing')*100;
  % % nerr = median(p.ap_unc./abs(p.ap), 'omitmissing')*100 - prctile(p.ap_unc./abs(p.ap).*100, 25);
  % % sc_r(1) = errshaded(lambda.a, median(p.ap_unc./abs(p.ap), 'omitmissing')*100, [perr; nerr], 'k', 0.1, '--', 1.5);
  % sc_r(1) = plot(lambda.a, median(p.ap_unc./abs(p.ap), 'omitmissing')*100,  '--k', 'LineWidth', 1.5);
  % 
  % sc_abs = [];
  % subplot(2,2,2); hold on
  % % perr = prctile(p.ap_unc, 97.5) - median(p.ap_unc, 'omitmissing');
  % % nerr = median(p.ap_unc, 'omitmissing') - prctile(p.ap_unc, 2.5);
  % % sc_abs(1) = errshaded(lambda.a, median(p.ap_unc, 'omitmissing'), [perr; nerr], 'k', 0.1, '--', 1.5);
  % sc_abs(1) = plot(lambda.a, median(p.ap_unc, 'omitmissing'), '--k', 'LineWidth', 1.5);
  % 
  % err_r = p.ap_unc(:,sela)./abs(p.ap(:,sela))*100;
  % subplot(2,2,3); hold on
  % h = histogram(err_r(err_r<200));
  % h.EdgeColor = h.FaceColor;
  % h.EdgeAlpha = 0.5;
  % h.FaceAlpha = 0.5;
  % 
  % err_a = p.ap_unc(:,sela);
  % subplot(2,2,4); hold on
  % h = histogram(err_a(err_a<0.015));
  % h.EdgeColor = h.FaceColor;
  % h.EdgeAlpha = 0.5;
  % h.FaceAlpha = 0.5;
  % 
  % delta_time_filt = hours(abs(filt_interp.dt - interp1([filt_avg.dt; filt_interp.dt(~isnan(filt_interp.fdom))], [filt_avg.dt; filt_interp.dt(~isnan(filt_interp.fdom))], filt_interp.dt, 'nearest')));
  % % compute uncertainty of linear interpolation 
  % sigma_lin_interp_a = abs(prctile(da_dt, 95, 1)) .* delta_time_filt;
  % sigma_lin_interp_c = abs(prctile(dc_dt, 95, 1)) .* delta_time_filt;
  % % Populate uncertainties
  % sigma_interp_a = sigma_lin_interp_a;
  % sigma_interp_c = sigma_lin_interp_c;
  % sigma_interp_a(~filt_interp.flag_linear_interp_a,:) = filt_interp.ag_interp_unc(~filt_interp.flag_linear_interp_a,:);
  % sigma_interp_c(~filt_interp.flag_linear_interp_c,:) = filt_interp.cg_interp_unc(~filt_interp.flag_linear_interp_c,:);
  % id_flow = ~isnan(filt_interp.flow_rate);
  % % divide total std and filtered std by number of independant spectra based on flow rate if available, or assuming a flow rate of 2LPM if no flow is available
  % % add up total standard error, filtered standard error, and linear interpolation standard error
  % p.ap_unc(id_flow, :) = sqrt((tot.a_avg_sd(id_flow, :)./sqrt(filt_interp.flow_rate(id_flow)./0.05)).^2 + ...
  %   (filt_interp.a_avg_sd(id_flow, :)./sqrt(filt_interp.flow_rate(id_flow)./0.05)).^2 + sigma_interp_a(id_flow, :).^2); % 0.05 L == volume of ACS tubes
  % p.ap_unc(~id_flow, :) = sqrt((tot.a_avg_sd(~id_flow, :)./sqrt(2./0.05)).^2 + ...
  %   (filt_interp.a_avg_sd(~id_flow, :)./sqrt(2./0.05)).^2 + sigma_interp_a(~id_flow, :).^2); % 0.05 L == volume of ACS tubes
  % p.cp_unc(id_flow, :) = sqrt((tot.c_avg_sd(id_flow, :)./sqrt(filt_interp.flow_rate(id_flow)./0.05)).^2 + ...
  %   (filt_interp.c_avg_sd(id_flow, :)./sqrt(filt_interp.flow_rate(id_flow)./0.05)).^2 + sigma_interp_c(id_flow, :).^2); % 0.05 L == volume of ACS tubes
  % p.cp_unc(~id_flow, :) = sqrt((tot.a_avg_sd(~id_flow, :)./sqrt(2./0.05)).^2 + ...
  %   (filt_interp.a_avg_sd(~id_flow, :)./sqrt(2./0.05)).^2 + sigma_interp_a(~id_flow, :).^2); % 0.05 L == volume of ACS tubes
  % 
  % subplot(2,2,1); 
  % % perr = prctile(p.ap_unc./abs(p.ap).*100, 75) - median(p.ap_unc./abs(p.ap), 'omitmissing')*100;
  % % nerr = median(p.ap_unc./abs(p.ap), 'omitmissing')*100 - prctile(p.ap_unc./abs(p.ap).*100, 25);
  % % sc_r(2) = errshaded(lambda.a, median(p.ap_unc./abs(p.ap), 'omitmissing')*100, [perr; nerr], 'b', 0.1, '-', 1.5);
  % sc_r(2) = plot(lambda.a, median(p.ap_unc./abs(p.ap), 'omitmissing')*100, '-b', 'LineWidth', 1.5);
  % ylabel('a_p % standard error')
  % xlabel('\lambda')
  % legend(sc_r, 'Linear interpolation', 'fDOM interpolation')
  % saveGraph('/Users/gui/Desktop/percentage_error_lin_fdom_interp', 'fig')
  % 
  % subplot(2,2,2);
  % % perr = prctile(p.ap_unc, 97.5) - median(p.ap_unc, 'omitmissing');
  % % nerr = median(p.ap_unc, 'omitmissing') - prctile(p.ap_unc, 2.5);
  % % sc_abs(2) = errshaded(lambda.a, median(p.ap_unc, 'omitmissing'), [perr; nerr], 'b', 0.1, '-', 1.5);
  % sc_abs(2) = plot(lambda.a, median(p.ap_unc, 'omitmissing'), '-b', 'LineWidth', 1.5);
  % ylabel('a_p standard error [m^{-1}]')
  % xlabel('\lambda')
  % legend(sc_abs, 'Linear interpolation', 'fDOM interpolation')
  % saveGraph('/Users/gui/Desktop/absolute_error_lin_fdom_interp', 'fig')
  % 
  % err_r = p.ap_unc(:,1:40)./abs(p.ap(:,1:40))*100;
  % subplot(2,2,3);
  % h = histogram(err_r(err_r<200));
  % h.EdgeColor = h.FaceColor;
  % h.EdgeAlpha = 0.5;
  % h.FaceAlpha = 0.5;
  % legend('Linear interpolation', 'fDOM interpolation')
  % xlabel('a_p % error')
  % 
  % err_a = p.ap_unc(:,sela);
  % subplot(2,2,4);
  % h = histogram(err_a(err_a<0.015));
  % h.EdgeColor = h.FaceColor;
  % h.EdgeAlpha = 0.5;
  % h.FaceAlpha = 0.5;
  % legend('Linear interpolation', 'fDOM interpolation')
  % xlabel('a_p absolute error')






  % populate number of spectra
  p.ap_n = tot.a_avg_n;
  p.cp_n = tot.c_avg_n;
  if any(strcmp(filt_interp.Properties.VariableNames, 'fdom'))
    p.fdom = filt_interp.fdom;
  end
  
  % % delete ap spectra full of NaNs
  % p(all(isnan(p.ap),2),:) = [];
  % p(all(isnan(p.cp),2),:) = [];
  
  % Model ag and cg based on FCDOM and Ra_FCDOM(lambda) uncertainties
  if strcmp(interpolation_method, 'CDOM') & any(strcmp(filt_interp.Properties.VariableNames, 'slope_interp_a'))
    p.ag_modelled = filt_interp.slope_interp_a .* p.fdom;
    p.cg_modelled = filt_interp.slope_interp_c .* p.fdom;
    % save Ra_FCDOM(lambda)
    da_dfdom_tbl = table();
    da_dfdom_tbl.wla = lambda.a(:);
    da_dfdom_tbl.da_dfdom = median(filt_interp.slope_interp_a)';
    da_dfdom_tbl.wlc = lambda.c(:);
    da_dfdom_tbl.dc_dfdom = median(filt_interp.slope_interp_c)';
    cruise_name = strsplit(FTH.path.prod, filesep);
    save(fullfile(FTH.path.prod, sprintf('%s_da_dfdom_%s_%s.mat', cruise_name{end-1}, datetime(min(p.dt),'Format','yyyyMMdd'), datetime(max(p.dt),'Format','yyyyMMdd'))), 'da_dfdom_tbl','fit_unc')
  end

  % Unsmoothing ACS spectra
  % Ron Zaneveld, WET Labs, Inc., 2005
  if size(lambda.a, 2) > 50 % unsmooth only ACS, not AC9
    p = unsmoothACS(p, lambda);
  end
  
  % Auto QC on product spectra
  wla_430 = lambda.a(find(lambda.a <= 430, 1,'last')); % find lower and closest to 430nm wavelength
  wla_700 = lambda.a(find(lambda.a >= 700, 1,'first')); % find higher and closest to 700nm wavelength
  
  % replace values at each end of the spectra when < -0.0015
  p.ap_unc(p.ap < -0.0015 & lambda.a < wla_430) = NaN;
  % p.ap_unc(p.ap < -0.0015 & lambda.a >= wla_700) = NaN;
  p.ap(p.ap < -0.0015 & lambda.a < wla_430) = NaN;
  % p.ap(p.ap < -0.0015 & lambda.a >= wla_700) = NaN;

  % set flag matrix
  % flag_varname = {...
  %   'deltaTS_not_corrected_total','deltaTS_not_corrected_filt0','deltaTS_not_corrected_filt1',...
  %   'fCDOM_mix_a_cluster0', 'fCDOM_mix_a_cluster1','fCDOM_mix_c_cluster0', 'fCDOM_mix_c_cluster1',...
  %   'fCDOM_a_cluster_chg', 'fCDOM_c_cluster_chg', 'flag_linear_interp_a', 'flag_linear_interp_c',...
  %   'a_filt0_not_clustered', 'a_filt1_not_clustered', 'c_filt0_not_clustered', 'c_filt1_not_clustered',...
  %   'flag_a_negative_slope', 'flag_c_negative_slope',...
  %   'cp_neg','ap_neg','ap_shape','ap_bubbles','cp_bubbles','ap430_700_neg','cp_over10',...
  %   'noisy600_650','ap460_640_04_450','positive_ap450_570','poc_flag','chl_ap676lh_flag',...
  %   'gamma_flag','chl_Halh_flag','HH_mphi_flag','HH_G50_flag','chlratio_flag',...
  %   'gamma_suspicious','poc_suspicious','chl_ap676lh_suspicious',...
  %   'chl_Halh_suspicious','HH_G50_mphi_suspicious'};
  flag_varname = {...
    'deltaTS_not_corrected_total','deltaTS_not_corrected_filt0','deltaTS_not_corrected_filt1',...
    'flag_linear_interp_a', 'flag_linear_interp_c','flag_Tresidual',...
    'cp_neg','ap_neg','ap_shape','ap_bubbles','cp_bubbles','ap430_700_neg','cp_over10',...
    'noisy600_650','ap460_640_04_450','positive_ap450_570','poc_flag','chl_ap676lh_flag',...
    'gamma_flag','chl_Halh_flag','HH_mphi_flag','HH_G50_flag','chlratio_flag',...
    'gamma_suspicious','poc_suspicious','chl_ap676lh_suspicious',...
    'chl_Halh_suspicious','HH_G50_mphi_suspicious'};
  flag = table('Size', [size(p, 1) size(flag_varname,2)], ...
    'VariableTypes', repmat({'logical'},1,size(flag_varname,2)), 'VariableNames', flag_varname);
  if strcmp(interpolation_method, 'CDOM')
    % flag.fCDOM_mix_a_cluster0 = filt_interp.fCDOM_mix_a_cluster0;
    % flag.fCDOM_mix_a_cluster1 = filt_interp.fCDOM_mix_a_cluster1;
    % flag.fCDOM_mix_c_cluster0 = filt_interp.fCDOM_mix_c_cluster0;
    % flag.fCDOM_mix_c_cluster1 = filt_interp.fCDOM_mix_c_cluster1;
    % flag.fCDOM_a_cluster_chg = filt_interp.fCDOM_a_cluster_chg;
    % flag.fCDOM_c_cluster_chg = filt_interp.fCDOM_c_cluster_chg;
    flag.flag_linear_interp_a = filt_interp.flag_linear_interp_a;
    flag.flag_linear_interp_c = filt_interp.flag_linear_interp_c;
    flag.flag_Tresidual = filt_interp.flag_Tresidual;
    % flag.a_filt0_not_clustered = filt_interp.flag_a_filt0_not_clustered;
    % flag.a_filt1_not_clustered = filt_interp.flag_a_filt1_not_clustered;
    % flag.c_filt0_not_clustered = filt_interp.flag_c_filt0_not_clustered;
    % flag.c_filt1_not_clustered = filt_interp.flag_c_filt1_not_clustered;
    % flag.cluster_a_negative_slope = filt_interp.flag_a_negative_slope;
    % flag.cluster_c_negative_slope = filt_interp.flag_c_negative_slope;
  else
    flag.flag_linear_interp_a = true(size(p,1), 1);
    flag.flag_linear_interp_c = true(size(p,1), 1);
  end
  flag.deltaTS_not_corrected_total = filt_interp.flag_deltaTSnotcorrected_total;
  flag.deltaTS_not_corrected_filt0 = filt_interp.flag_deltaTSnotcorrected_filt0;
  flag.deltaTS_not_corrected_filt1 = filt_interp.flag_deltaTSnotcorrected_filt1;

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
    p.cp_unc(lambda.c > lambda.c(minred_wlc)' & lambda.c > 720) = NaN;
    p.ap(lambda.a > lambda.a(minred_wla)' & lambda.a > 715) = NaN;
    p.ap_unc(lambda.a > lambda.a(minred_wla)' & lambda.a > 715) = NaN;
    
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
      p.ap_unc(lambda.a < cutoffblue_wla & lambda.a <= 450) = NaN;
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
    % flag for the step in the center of the ap spectra
    % idstep = lambda.a > 560 & lambda.a <= 600;
    % idref = lambda.a > 500 & lambda.a <= 560;
    idstep = lambda.a > 550 & lambda.a <= 600;
    idref = lambda.a > 450 & lambda.a <= 550;
    id440 = abs(lambda.a - 440) == min(abs(lambda.a - 440));
    foo_difstep = diff(p.ap(:,idstep),1, 2) ./ diff(lambda.a(idstep));
    foo_difref = diff(p.ap(:,idref),1, 2) ./ diff(lambda.a(idref));
    % toflag = any(abs(foo_difstep) > 3 * prctile(abs(foo_difref), 95, 2) & abs(foo_difstep) > 0.005, 2); % old version before 2025-07-21
    toflag = any(max(abs(foo_difstep),[],2)  > 3 * prctile(abs(foo_difref),50,2) & max(abs(foo_difstep),[],2) > 0.005.*p.ap(:, id440), 2);
    if sum(toflag)
      fprintf('Signal of bubbles/large particles detected in %.2f%% (%i) of the ap spectra, flagged: step in ap between 550 and 600 nm\n', ...
        sum(toflag) / size(p, 1) * 100, sum(toflag))
      flag.ap_bubbles(toflag) = true;
    end
    % flag for the step in the center of the cp spectra
    % idstep = lambda.c > 560 & lambda.c <= 600;
    % idref = lambda.c > 500 & lambda.c <= 560;
    idstep = lambda.c > 550 & lambda.c <= 600;
    idref = lambda.c > 450 & lambda.c <= 550;
    id440 = abs(lambda.c - 440) == min(abs(lambda.c - 440));
    foo_difstep = diff(p.cp(:,idstep),1, 2) ./ diff(lambda.c(idstep));
    foo_difref = diff(p.cp(:,idref),1, 2) ./ diff(lambda.c(idref));
    % toflag = any(abs(foo_difstep) > 3 * prctile(abs(foo_difref), 95, 2) & abs(foo_difstep) > 0.005, 2); % old version before 2025-07-21
    toflag = any(max(abs(foo_difstep),[],2)  > 3 * prctile(abs(foo_difref),50,2) & max(abs(foo_difstep),[],2) > 0.005.*p.cp(:, id440), 2);
    if sum(toflag)
      fprintf('Signal of bubbles/large particles detected in %.2f%% (%i) of the cp spectra, flagged: step in cp between 550 and 600 nm\n', ...
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
  agaus = GaussDecomp(p, lambda.a, compute_ad_aphi, nap_offset);
  p = [p agaus];
  
  fprintf('Computing Chl line height, POC & gamma ... ')
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
  ap_a = interp1(lambda.a, p.ap',[650 676 715],'spline')';
  p.ap676_lh = ap_a(:,2)-(39/65*ap_a(:,1)+26/65*ap_a(:,3));  % ap_a650-(39/65*ap_a(:,1)+26/65*ap_a(:,3));
  flag.chl_ap676lh_flag(real(p.ap676_lh) ~= p.ap676_lh) = true;
  p.ap676_lh(real(p.ap676_lh) ~= p.ap676_lh) = NaN;
  flag.chl_ap676lh_flag(p.ap676_lh < 0) = true;
  p.ap676_lh(p.ap676_lh < 0) = NaN;
  p.chl_ap676lh = 157*p.ap676_lh.^1.22;
  
  % 3.3 Derive Gamma (does not support NaN values) (Boss et al. 2001)
  % REFERENCES:
  % Boss, E., W.S. Pegau, W.D. Gardner, J.R.V. Zaneveld, A.H. Barnard., M.S. Twardowski, G.C. Chang, and T.D. Dickey, 2001. Spectral particulate attenuation and particle size distribution in the bottom boundary layer of a continental shelf. Journal of Geophysical Research, 106, 9509-9516.
  p.gamma = NaN(size(p, 1), 1);
  cp_temp_gam = p.cp(~all(isnan(p.cp),2), :);
  sel = ~any(isnan(cp_temp_gam));
  [~,temp_gam] = FitSpectra_HM2(lambda.c(:,sel), cp_temp_gam(:,sel));
  p.gamma(~all(isnan(p.cp),2)) = temp_gam;
  flag.gamma_flag(p.gamma < 0) = true;
  % p.gamma(p.gamma < 0) = NaN;
  fprintf('Done\n')
  
  fprintf('Computing chlorophyll from cp (H_alh) ... ')
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
      best_di = table();
      best_di.a = NaN(size(di_orig,1), 1);
      best_di.c = NaN(size(di_orig,1), 1);
      best_di.assign_loc = cell(size(di_orig,1), 1);
      for i = 1:size(di_orig,1)
        if i == 1 || i == size(di_orig,1)
          best_di.assign_loc{i} = abs(di.dt(i) - di.dt) < hours(240); % 120 96 72 48
        else
          best_di.assign_loc{i} = abs(di.dt(i) - di.dt) < hours(220); % 72 48 36 24
        end
        % lowest_di_a = di_orig.a(:, lambda.a >= 550 & lambda.a <= 650) == ...
        %   min(di_orig.a(iddi, lambda.a >= 550 & lambda.a <= 650), [], 1);
        % lowest_di_c = di_orig.c(:, lambda.c >= 550 & lambda.c <= 650) == ...
        %   min(di_orig.c(iddi, lambda.c >= 550 & lambda.c <= 650), [], 1);
        % foo_a = find(sum(lowest_di_a, 2) == max(sum(lowest_di_a, 2)));
        % foo_c = find(sum(lowest_di_c, 2) == max(sum(lowest_di_c, 2)));
        best_di.a(i) = find(mean(di_orig.a(:, lambda.a >= 550 & lambda.a <= 650), 2) == ...
          min(mean(di_orig.a(best_di.assign_loc{i}, lambda.a >= 550 & lambda.a <= 650), 2), [], 1), 1, 'first');
        best_di.c(i) = find(mean(di_orig.c(:, lambda.c >= 550 & lambda.c <= 650), 2) == ...
          min(mean(di_orig.c(best_di.assign_loc{i}, lambda.c >= 550 & lambda.c <= 650), 2), [], 1), 1, 'first');
      end
      for i = 1:size(di_orig,1)
        di.a(best_di.assign_loc{i}, :) = repmat(di_orig.a(best_di.a(i), :), sum(best_di.assign_loc{i}), 1);
        di.c(best_di.assign_loc{i}, :) = repmat(di_orig.c(best_di.c(i), :), sum(best_di.assign_loc{i}), 1);
        di.a_avg_sd(best_di.assign_loc{i}, :) = repmat(di_orig.a_avg_sd(best_di.a(i), :), sum(best_di.assign_loc{i}), 1);
        di.c_avg_sd(best_di.assign_loc{i}, :) = repmat(di_orig.c_avg_sd(best_di.c(i), :), sum(best_di.assign_loc{i}), 1);
      end
    
      % visProd3D(lambda.a, di.dt, di.a, false, 'Wavelength', false, 5); zlabel('a_diw (m^{-1})');
      % visProd3D(lambda.c, di.dt, di.c, false, 'Wavelength', false, 6); zlabel('c_diw (m^{-1})');

    end

    % duplicate DI from the start of each section to the end of each DI section
    di_end = di(1:end-1, :);
    di_end.dt = di.dt(2:end) - minutes(1);
    di = [di; di_end];
    di = sortrows(di, 'dt');
    if max(di.dt) < max(filt_avg.dt)
      di.dt(end) = max(filt_avg.dt) + minutes(1);
    end

% visProd3D(lambda.a, di.dt, di.a, false); zlabel('a_f (m^{-1})');
% visProd3D(lambda.a, di_interp.dt, di_interp.a, false); zlabel('a_f (m^{-1})');

    % remove when a and c are full of NaNs
    filt_avg(all(isnan(filt_avg.a), 2) & all(isnan(filt_avg.c), 2),:) = [];
    
    % Interpolate DI on filtered
    di_interp = table(filt_avg.dt, 'VariableNames', {'dt'});
    di_interp.a = interp1(di.dt, di.a, di_interp.dt, 'linear', 'extrap');
    di_interp.c = interp1(di.dt, di.c, di_interp.dt, 'linear', 'extrap');
    di_interp.a_avg_sd = interp1(di.dt, di.a_avg_sd, di_interp.dt, 'linear', 'extrap');
    di_interp.c_avg_sd = interp1(di.dt, di.c_avg_sd, di_interp.dt, 'linear', 'extrap');
  
    % Dissolved = Filtered - DI
    g = table(filt_avg.dt, 'VariableNames', {'dt'});
    g.ag = filt_avg.a - di_interp.a;
    g.cg = filt_avg.c - di_interp.c;
    % set flag matrix
    flag_g_varname = {'deltaTSnotcorrected_filt'};
    flag_g = table('Size', [size(g,1) size(flag_g_varname,2)], ...
      'VariableTypes', repmat({'logical'},1,size(flag_g_varname,2)), 'VariableNames', flag_g_varname);
    flag_g.deltaTSnotcorrected_filt = filt_avg.flag_deltaTSnotcorrected_filt;

    % Interpolate wavelength of c on a 
    % g.cg = cell2mat(arrayfun(@(i) interp1(c_wl, g.cg(i,:), a_wl, 'linear', 'extrap'), 1:size(g,1), 'UniformOutput', false)');
    g.ag = interp1(lambda.a', g.ag', lambda.ref', 'linear', 'extrap')';
    g.cg = interp1(lambda.c', g.cg', lambda.ref', 'linear', 'extrap')';
    fprintf('Correcting for temperature & salinity dependence ... ')
    % Temperature & Salinity Correction (No Scattering correction needed)
    [g.ag, g.cg] = TemperatureAndSalinityDependence(g.ag, g.cg, lambda.a, lambda.c, psi, flag_g.deltaTSnotcorrected_filt);
    fprintf('Done\n')
  
    % correct dissolved data for biofouling
    ida450 = abs(lambda.a - 450) == min(abs(lambda.a - 450));
    idc450 = abs(lambda.c - 450) == min(abs(lambda.c - 450));
    [data_corrected, DIW_biofouling_correction] = BiofoulingCorrection(g, lambda, ...
      {'ag', 'cg'}, find(ida450), flow_data, flow_data.(FTH.view.spd_variable), cdom_base);
    
% visProd3D(lambda.a, data_corrected.dt, data_corrected.ag, false, 'Wavelength', false, 2); zlabel('a_g (m^{-1})');
% visProd3D(lambda.c, data_corrected.dt, data_corrected.cg, false, 'Wavelength', false, 3); zlabel('c_g (m^{-1})');
% visProd3D(lambda.c, data_corrected.dt, data_corrected.ag - data_corrected.cg, false, 'Wavelength', false, 4); zlabel('a_g - c_g (m^{-1})');
% 
% idwl = abs(lambda.a - 450) == min(abs(lambda.a - 450));
% figure; subplot(1,3,1)
% plot_linreg(filt_avg.fdom, data_corrected.ag(:, idwl), 'robust', 'linear', false)
% xlabel('fCDOM (v)')
% ylabel('a_g 450 nm [m^-1]')
% subplot(1,3,2)
% plot_linreg(filt_avg.fdom, data_corrected.cg(:, idwl), 'robust', 'linear', false)
% xlabel('fCDOM (v)')
% ylabel('c_g 450 nm [m^-1]')
% subplot(1,3,3)
% scatter(filt_avg.dt, data_corrected.ag(:, idwl) - data_corrected.cg(:, idwl), 50, 'filled')
% ylabel('a_g450 - c_g450 [m^-1]')
% set(gca, 'FontSize', 14)

    % Propagate error
    %   Note: Error is not propagated through temprature and salinity dependance correction as required by SeaBASS
    data_corrected.ag_sd = sqrt(filt_avg.a_avg_sd.^2 + di_interp.a_avg_sd.^2);
    data_corrected.cg_sd = sqrt(filt_avg.c_avg_sd.^2 + di_interp.c_avg_sd.^2);
    data_corrected.ag_n = filt_avg.a_avg_n;
    data_corrected.cg_n = filt_avg.c_avg_n;
    if any(strcmp(filt_avg.Properties.VariableNames, 'fdom'))
      data_corrected.fdom = filt_avg.fdom;
    end

% visProd3D(lambda.a, data_corrected.dt, data_corrected.ag, false, 'Wavelength', false, 5); zlabel('a_g raw (m^{-1})');
% visProd3D(lambda.a, data_corrected.dt, data_corrected.cg, false, 'Wavelength', false, 6); zlabel('c_g raw (m^{-1})');

    % force all ag > 0 to fit exponential
    forcedto0 = data_corrected;
    forcedto0.ag = data_corrected.ag - min(data_corrected.ag, [], 2);
    % force all cg > 0 to fit exponential
    forcedto0.cg = data_corrected.cg - min(data_corrected.cg, [], 2);

% visProd3D(lambda.a, forcedto0.dt, forcedto0.ag, false, 'Wavelength', false, 5); zlabel('a_g forced 0 (m^{-1})');
% visProd3D(lambda.a, forcedto0.dt, forcedto0.cg, false, 'Wavelength', false, 6); zlabel('c_g forced 0 (m^{-1})');
    
    % First fit exponential to estimate the offset to correct for in NIR
    fprintf('Finding ag and cg NIR offset ... \n')
    sel_a = lambda.a >= 440 & lambda.a < 580 & ~any(isnan(forcedto0.ag(~all(isnan(forcedto0.ag),2), :)));
    sel_c = lambda.c >= 440 & lambda.c < 580 & ~any(isnan(forcedto0.cg(~all(isnan(forcedto0.cg),2), :)));
    % fit exponential on ag
    [forcedto0.y_intercp_fit_ag, forcedto0.base_fit_ag, ~, ~, forcedto0.RMSE_fit_ag] = FitExp(lambda.a(sel_a), ...
      forcedto0.ag(:, sel_a), forcedto0.ag_sd(:, sel_a));
    % fit exponential on cg
    [forcedto0.y_intercp_fit_cg, forcedto0.base_fit_cg, ~, ~, forcedto0.RMSE_fit_cg] = FitExp(lambda.c(sel_c), ...
      forcedto0.cg(:, sel_c), forcedto0.cg_sd(:, sel_c));

    % rebuild ag from exponential fit parameters
    ag_rebuilt = NaN(size(forcedto0.ag));
    cg_rebuilt = NaN(size(forcedto0.cg));
    % define exponential function
    expfun = @(p, xd) p(1) * exp(p(2) * (xd - 440));
    for i = 1:size(forcedto0, 1)
      ag_rebuilt(i, :) = expfun([forcedto0.y_intercp_fit_ag(i) forcedto0.base_fit_ag(i)], lambda.a);
      cg_rebuilt(i, :) = expfun([forcedto0.y_intercp_fit_cg(i) forcedto0.base_fit_cg(i)], lambda.c);
    end
    % add back minimum offset removed before exponential fit
    ag_rebuilt = ag_rebuilt + min(data_corrected.ag, [], 2);
    cg_rebuilt = cg_rebuilt + min(data_corrected.cg, [], 2);

% visProd3D(lambda.a, forcedto0.dt, ag_rebuilt, false, 'Wavelength', false, 5); zlabel('a_g forced 0 rebuilt (m^{-1})');
% visProd3D(lambda.c, forcedto0.dt, cg_rebuilt, false, 'Wavelength', false, 6); zlabel('c_g forced 0 rebuilt (m^{-1})');

    % remove offset in NIR to the entire ag/cg spectra to force to 0
    forced_NiRto0 = data_corrected;
    NIR_offset_ag = ag_rebuilt(:, end);
    NIR_offset_cg = cg_rebuilt(:, end);
    forced_NiRto0.ag = data_corrected.ag - NIR_offset_ag;
    forced_NiRto0.cg = data_corrected.cg - NIR_offset_cg;

% visProd3D(lambda.a, forced_NiRto0.dt, forced_NiRto0.ag, false, 'Wavelength', false, 5); zlabel('a_g forced 0 at 650 nm (m^{-1})');
% visProd3D(lambda.c, forced_NiRto0.dt, forced_NiRto0.cg, false, 'Wavelength', false, 6); zlabel('c_g forced 0 at 650 nm (m^{-1})');


% % rebuild ag from exponential fit parameters
% ag_rebuilt_f = NaN(size(forcedto0.ag));
% cg_rebuilt_f = NaN(size(forcedto0.cg));
% % define exponential function
% expfun = @(p, xd) p(1) * exp(p(2) * (xd - 440));
% for i = 1:size(forced_NiRto0, 1)
%   ag_rebuilt_f(i, :) = expfun([y_intercp_fit_ag(i) base_fit_ag(i)], lambda.a);
%   cg_rebuilt_f(i, :) = expfun([y_intercp_fit_cg(i) base_fit_cg(i)], lambda.c);
% end
% 
% visProd3D(lambda.a, forcedto0.dt, ag_rebuilt_f, false, 'Wavelength', false, 7); zlabel('a_g forced 0 rebuilt_f (m^{-1})');
% visProd3D(lambda.c, forcedto0.dt, cg_rebuilt_f, false, 'Wavelength', false, 8); zlabel('c_g forced 0 rebuilt_f (m^{-1})');

% 
% idwl = abs(lambda.a - 450) == min(abs(lambda.a - 450));
% figure; subplot(1,3,1)
% plot_linreg(filt_avg.fdom, forced_NiRto0.ag(:, idwl), 'robust', 'linear', false)
% xlabel('fCDOM (v)')
% ylabel('a_g450 forced 0 at 650 nm[m^-1]')
% subplot(1,3,2)
% plot_linreg(filt_avg.fdom, forced_NiRto0.cg(:, idwl), 'robust', 'linear', false)
% xlabel('fCDOM (v)')
% ylabel('c_g450 forced 0 at 650 nm [m^-1]')
% subplot(1,3,3)
% scatter(filt_avg.dt, forced_NiRto0.ag(:, idwl) - forced_NiRto0.cg(:, idwl), 50, 'filled')
% ylabel('a_g450 - c_g450 forced 0 at 650 nm [m^-1]')
% set(gca, 'FontSize', 14)
% 
% figure; hold on
% scatter(g.dt, g.base_fit_ag, 50, 'filled')
% scatter(g.dt, g.base_fit_cg, 50, 'filled')
% figure; hold on
% scatter(g.dt, g.y_intercp_fit_ag, 50, 'filled')
% scatter(g.dt, g.y_intercp_fit_cg, 50, 'filled')
    
    % if absolute mean difference between ag and cg is lower after forcing use forced data otherwise original data
    if mean(abs(forced_NiRto0.ag(:, ida450) - forced_NiRto0.cg(:, idc450)), 'omitmissing') < ...
        mean(abs(data_corrected.ag(:, ida450) - data_corrected.cg(:, idc450)), 'omitmissing')
      g = forced_NiRto0;
      g.ag_NIRoffset = NIR_offset_ag;
      g.cg_NIRoffset = NIR_offset_cg;
    else
      g = data_corrected;
      g.ag_NIRoffset = zeros(size(g.dt));
      g.cg_NIRoffset = zeros(size(g.dt));
    end
    
  %   visProd3D(lambda.a, g.dt, g.ag, false, 'Wavelength', false, 70);
  %   title('ag with auto best DI 72h')
  %   saveGraph('ag_with_auto_best_DI_72h', 'fig')
  %   
  %   visProd3D(lambda.c, g.dt, g.cg, false, 'Wavelength', false, 71);
  %   title('cg with auto best DI 72h')
  %   saveGraph('cg_with_auto_best_DI_72h', 'fig')
  
  %   visProd3D(lambda.a, di_interp.dt, di_interp.a, false, 'Wavelength', false, 72);
  %   visProd3D(lambda.c, di_interp.dt, di_interp.c, false, 'Wavelength', false, 73);
  
    % First fit final exponential
    fprintf('Exponential fit to ag and cg ... \n')
    % fit exponential on ag
    [g.y_intercp_fit_ag, g.base_fit_ag, ~, ~, g.RMSE_fit_ag] = FitExp(lambda.a(sel_a), ...
      g.ag(:, sel_a), g.ag_sd(:, sel_a));
    % add fit flag
    g.ag_fitflag = false(size(g, 1), 1);
    g.ag_fitflag(g.RMSE_fit_ag > 0.0025) = true;
    % fit exponential on cg
    [g.y_intercp_fit_cg, g.base_fit_cg, ~, ~, g.RMSE_fit_cg] = FitExp(lambda.c(sel_c), ...
      g.cg(:, sel_c), g.cg_sd(:, sel_c));
    % add fit flag
    g.cg_fitflag = false(size(g, 1), 1);
    g.cg_fitflag(g.RMSE_fit_cg > 0.0025) = true;
    fprintf('Done\n')

    % keep only data 440 < data < 580 nm (blue wl too noisy and salinity & temperature impact > 580 nm)
    g.ag(:, ~sel_a) = NaN;
    g.ag_sd(:, ~sel_a) = NaN;
    g.cg(:, ~sel_c) = NaN;
    g.cg_sd(:, ~sel_c) = NaN;

  %   % QC with ag and cg spectra (limited testing on the QC)
  %   g(g.ag(:,1) < 0 & g.cg(:,end-3) < -0.005, :) = [];

    % set flag column
    g.flag_bit = set_flagbit(flag_g);
    % flag_info = read_flagbit(g.flag_bit, 'ACS');
    fprintf('Done\n')
  else
    g = table();
  end
end


%% Temperature And Salinity Dependence correction
function [a_corr, c_corr, a_slope, c_slope] = TemperatureAndSalinityDependence(a, c, wl_a, wl_c, psi, salinity_correction)
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
    if salinity_correction(k)
      % Force Temperature, Salinity, and Slope
      x = fminsearch(@costfun_TSD, [deltaT, deltaS, amp, slope], opts, a(k,iwla), a_psiT(iwla), a_psiS(iwla), wl_a(iwla), salinity_correction(k));
      a_deltaT(k) = x(1); a_deltaS(k) = x(2); a_amp(k) = x(3); a_slope(k) = x(4);
      x = fminsearch(@costfun_TSD, [deltaT, deltaS, amp, slope], opts, c(k,iwlc), c_psiT(iwlc), c_psiS(iwlc), wl_c(iwlc), salinity_correction(k));
      c_deltaT(k) = x(1); c_deltaS(k) = x(2); c_amp(k) = x(3); c_slope(k) = x(4);
    else
      % Without Salinity forcing
      x = fminsearch(@costfun_TSD, [deltaT, amp, slope], opts, a(k,iwla), a_psiT(iwla), a_psiS(iwla), wl_a(iwla), salinity_correction(k));
      a_deltaT(k) = x(1); a_amp(k) = x(2); a_slope(k) = x(3);
      x = fminsearch(@costfun_TSD, [deltaT, amp, slope], opts, c(k,iwlc), c_psiT(iwlc), c_psiS(iwlc), wl_c(iwlc), salinity_correction(k));
      c_deltaT(k) = x(1); c_amp(k) = x(2); c_slope(k) = x(3);
    end
  end
  
  % Apply correction
  if salinity_correction
    a_corr = a - a_deltaT.*a_psiT - a_deltaS.*a_psiS;
    c_corr = c - c_deltaT.*c_psiT - c_deltaS.*c_psiS;
  else
    a_corr = a - a_deltaT.*a_psiT;
    c_corr = c - c_deltaT.*c_psiT;
  end
  
  % Display the forcing parameters
  % disp([a_deltaT, c_deltaT, a_deltaS, c_deltaS, a_slope, c_slope])
end

% cost function for Temperature And Salinity Dependence correction
function cost = costfun_TSD(x, spectra, psiT, psiS, wl, salinity_correction)
  if salinity_correction
    % Force Temperature, Salinity, and Slope
    cost = sum((spectra - psiT.*x(1) - psiS .* x(2) - x(3) .* exp(-x(4)*(wl-450))).^2);
  else
    % Without Salinity forcing
    cost = sum((spectra - psiT.*x(1) - x(2).*exp(-x(3)*(wl-450))).^2);
  end
end

%% Residual Temperature And Scattering Correction (Zaneveld 1994 proportional)
function [ap_final, cp_final, flag_Tresidual] = ResidualTempScatterCorrZaneveld_proportional(ap_approx, cp_approx, wl, psi, dt)
  fprintf('...')
  % Function from Emmanuel Boss, Nils Haëntjens, and Guillaume Bourdin after Zaneveld 1994 method 3: assumes negligible ap in NIR
  % Find Near Infrared & references
  iNIR = 710 <= wl &  wl <= 750;  % spectral srange for optimization (710 to 750nm)
  if isempty(iNIR); error('Unable to perform correction as no wavelength available in NIR.'); end
  % Find nearest wavelength to greater than 715 nm to use as reference for correction
  iref = find(abs(wl - 715) == min(abs(wl - 715)), 1, 'first'); % 715 730
  % If ACS spectra do not go up to 715 nm take the closest wavelength to 715 nm
  if isempty(iref); [~, iref] = max(wl); end % works as there is data in iNIR so lowest wavelength is 710
  
  % approximate bp
  bp_approx = cp_approx - ap_approx;
  
  % correct for residual temperature signal
  [ap_Tcorr, cp_final, flag_Tresidual] = ResidualTemperatureCorrection(wl, ap_approx, cp_approx, psi, iNIR, iref, dt, NaN(size(dt)), 'Zaneveld1994_proportional');
  
  % apply flat scattering correction
  ap_final = ap_Tcorr - ap_Tcorr(:,iref) ./ bp_approx(:,iref) .* bp_approx; 
end

%% Residual Temperature And Scattering Correction (Rottgers2013 semiempirical)
function [ap_final, cp_final, flag_Tresidual] = ResidualTempScatterCorrRottgers_semiempirical(ap_approx, cp_approx, wl, psi, dt)
  fprintf('...')
  % Function from Emmanuel Boss and Guillaume Bourdin after Rottgers et al. 2013: allows absorption in NIR in case of high NAP
  % Find Near Infrared & references
  iNIR = 710 <= wl &  wl <= 750;  % spectral srange for optimization (710 to 750nm)
  if isempty(iNIR); error('Unable to perform correction as no wavelength available in NIR.'); end
  % Find nearest wavelength to greater than 715 nm to use as reference for correction
  iref = find(abs(wl - 715) == min(abs(wl - 715)), 1, 'first'); % 715 730
  % If ACS spectra do not go up to 715 nm take the closest wavelength to 715 nm
  if isempty(iref); [~, iref] = max(wl); end % works as there is data in iNIR so lowest wavelength is 710
  
  % Define empirical function of scattering correction at lambda reference
  ap715 = @(ap) 0.212 * ap .^ 1.135;

  % correct for residual temperature signal
  [ap_Tcorr, cp_final, flag_Tresidual] = ResidualTemperatureCorrection(wl, ap_approx, cp_approx, psi, iNIR, iref, dt, ap715(ap_approx(:,iref)), 'Rottgers2013_semiempirical');
  
  % select ap reference and replace negative values by NaN to prevent leakage
  ap_Tcorr_iref = ap_Tcorr(:,iref);
  ap_Tcorr_iref(ap_Tcorr_iref < 0) = NaN;
  % apply flat scattering correction
  % ap_final = ap_Tcorr - (ap_Tcorr_iref - 0.212*ap_Tcorr_iref.^1.135);
  ap_final = ap_Tcorr - (ap_Tcorr_iref - ap715(ap_Tcorr(:,iref)));
end

%% Blended 1 residual Temperature And Scattering Correction using cp shape (Slade 2015 + NIR offset from Rottgers2013 semiempirical and Bourdin et al. in prep)
% function [ap_final, cp_final, flag_Tresidual] = ResidualTempScatterCorrSlade_blended(ap_approx, cp_approx, wl, psi, dt)
function [ap_final, cp_final, flag_Tresidual] = ResidualTempScatterCorrSemiempirical_blended1(ap_approx, cp_approx, wl, psi, dt)
  % Function from Emmanuel Boss and Guillaume Bourdin after Slade and Boss 2015 and Rottgers et al. 2013: 
  % Proportional scattering correction + allows absorption in NIR in case of high NAP
  fprintf('.')
  % Find Near Infrared & references
  iNIR = 710 <= wl &  wl <= 750;  % spectral srange for optimization (710 to 750nm)
  if isempty(iNIR); error('Unable to perform correction as no wavelength available in NIR.'); end
  % Find nearest wavelength to greater than 730 nm to use as reference for correction
  iref = find(abs(wl - 715) == min(abs(wl - 715)), 1, 'first'); % 715 730
  % If ACS spectra do not go up to 715 nm take the closest wavelength to 715 nm
  if isempty(iref); [~, iref] = max(wl); end % works as there is data in iNIR so lowest wavelength is 710
  
  % Define empirical function of scattering correction at lambda reference
  ap715 = @(ap) 0.0831 * ap;

  % correct for residual temperature signal
  [ap_Tcorr, cp_Tcorr, flag_Tresidual] = ResidualTemperatureCorrection(wl, ap_approx, cp_approx, psi, iNIR, iref, dt, ap715(ap_approx(:,iref)), 'Semiempirical_blended1');
  fprintf('.')

  % Apply scattering correction adjusting for acceptance angle in c tube after Rottgers et al. 2013: 
  % Evaluation of scatter corrections for ac-9 absorption measurements in coastal waters, DOI: 10.1016/j.mio.2013.11.001

  % % Out of Slade and Boss, 2015
  % ec = 1; % acceptance angle correcton
  % ec = 1/0.56; % acceptance angle correcton

  % % succesful attempt with bp
  % ap_corr = ap_Tcorr - (ap_Tcorr(:,iref) - ap715(ap_Tcorr(:,iref))) .* (ec.*cp_Tcorr - ap_Tcorr) ./ (ec.*cp_Tcorr(:,iref) - ap_Tcorr(:,iref));
  % succesful attempt with cp instead of bp
  ap_corr = ap_Tcorr - (ap_Tcorr(:,iref) - ap715(ap_Tcorr(:,iref))) .* cp_Tcorr ./ cp_Tcorr(:,iref);
  fprintf('.')

  % correct again for residual temperature signal after scattering correction
  [ap_final, cp_final, flag_Tresidual2] = ResidualTemperatureCorrection(wl, ap_corr, cp_Tcorr, psi, iNIR, iref, dt, ap715(ap_approx(:,iref)), 'Semiempirical_blended1');
  flag_Tresidual(flag_Tresidual2) = true;


  % figure; hold on
  % plot(wl, ap_Tcorr(100, :))
  % plot(wl, ap_corr(100, :))
  % plot(wl, ap_final_Z(100, :))
  % legend('ap\_Tcorr','ap\_final','ap\_final_{Zaneveld}')


  % zoom_in = wl > 550;
  % visProd3D(wl(zoom_in), dt, ap_approx(:,zoom_in), false, 'Wavelength', false, 12);
  % title('ap')
  % visProd3D(wl(zoom_in), dt, cp_approx(:,zoom_in), false, 'Wavelength', false, 22);
  % title('cp')
  % 
  % zoom_in = wl > 550;
  % visProd3D(wl(zoom_in), dt, ap_corr(:,zoom_in), false, 'Wavelength', false, 13);
  % title('ap')
  % visProd3D(wl(zoom_in), dt, cp_Tcorr(:,zoom_in), false, 'Wavelength', false, 23);
  % title('cp')
  % 
  % zoom_in = wl > 550;
  % visProd3D(wl(zoom_in), dt, ap_final(:,zoom_in), false, 'Wavelength', false, 14);
  % title('ap')
  % visProd3D(wl(zoom_in), dt, cp_final(:,zoom_in), false, 'Wavelength', false, 24);
  % title('cp')
  % 
  % 
  % 
  % visProd3D(wl, dt, ap_approx, false, 'Wavelength', false, 12);
  % title('ap')
  % visProd3D(wl, dt, cp_approx, false, 'Wavelength', false, 22);
  % title('cp')
  % 
  % visProd3D(wl, dt, ap_corr, false, 'Wavelength', false, 13);
  % title('ap')
  % visProd3D(wl, dt, cp_Tcorr, false, 'Wavelength', false, 23);
  % title('cp')
  % 
  % visProd3D(wl, dt, ap_final, false, 'Wavelength', false, 14);
  % title('ap')
  % visProd3D(wl, dt, cp_final, false, 'Wavelength', false, 24);
  % title('cp')
  % 
  % 
  % 
  % visProd3D(wl, dt, ap_corr-ap_final, false, 'Wavelength', false, 14);
  % title('ap')
  % visProd3D(wl, dt, cp_Tcorr-cp_final, false, 'Wavelength', false, 24);
  % title('cp')
  % 
  % % ap_iref_function = @(ap) 0.0854 * ap;
  % % ap_iref_function = @(ap) 0.0854 * ap;
  % 
  % 
  % 
  % 
  % 
  % 
  % % Select good spectra to run minimization routine on
  % ap_approx2 = NaN(size(ap_approx));
  % cp_approx2 = NaN(size(cp_approx));
  % % ap_it_diff = NaN(size(ap_approx));
  % 
  % max_iteration = 10;
  % i = 1;
  % 
  % 
  % 
  % 
  % 
  % 
  % 
  % figure(11); hold on
  % while i <= max_iteration && ((all(isnan(ap_approx2),"all") && all(isnan(cp_approx2),"all")) || any(abs(ap_it_diff) > 0.001, "all"))
  %   % Init routine parameter scattering correction
  %   bp_approx = cp_approx - ap_approx;
  %   sel = find(all(isfinite(ap_approx),2));
  %   for k = sel'
  %     deltaT(k) = fminsearch(@costFun_RTSC, 0, opts, ap_approx(k,:), bp_approx(k,:), psiT, iNIR, iref);         
  %   end
  %   % Interpolate linearly deltaT in case of missing few wavelength data in spectra (up to 10 min gaps)
  %   deltaT = fillmissing(deltaT,'linear',1,'MaxGap',minutes(10),'SamplePoints',dt);
  %   flag_deltaT(sel) = false;
  %   % Apply temperature correction and replace negative values at iref by NaN to avoid complex solution leakage
  %   ap_Tcorr = ap_approx - psiT.*deltaT;
  %   cp_approx2 = cp_approx - psiT.*deltaT;
  %   ap_Tcorr_iref = ap_Tcorr(:,iref);
  %   ap_Tcorr_iref(ap_Tcorr_iref < 0) = NaN;
  %   bp_Tcorr = cp_approx2 - ap_Tcorr;
  %   % Apply scattering correction
  %   ap_approx2 = ap_Tcorr - ap_Tcorr_iref ./ bp_Tcorr(:,iref) .* bp_Tcorr + 0.0854*ap_Tcorr_iref; % 0.212*ap_Tcorr_iref.^1.135;
  %   % Check difference with previous iteration
  %   ap_it_diff = ap_approx - ap_approx2;
  %   cp_it_diff = cp_approx - cp_approx2;
  % 
  %   ap_approx = ap_approx2;
  %   cp_approx = cp_approx2;
  %   i = i + 1;
  % 
  %   % figure(15)
  %   % plot(wl, mean(ap_it_diff,1,'omitmissing'))
  %   % ylabel('ap difference between iteration')
  %   % 
  %   % figure;
  %   % histogram(ap_it_diff)
  %   % xlabel('ap difference between iteration')
  % 
  %   figure(11); subplot(1,2,1); hold on
  %   plot(wl, mean(ap_it_diff,1,'omitmissing'))
  %   ylabel('ap difference between iteration')
  %   subplot(1,2,2); hold on
  %   plot(wl, mean(cp_it_diff,1,'omitmissing'))
  %   ylabel('ap difference between iteration')
  % 
  %   figure; subplot(1,2,1)
  %   histogram(ap_it_diff)
  %   xlabel('ap difference between iteration')
  %   subplot(1,2,2)
  %   histogram(cp_it_diff)
  %   xlabel('cp difference between iteration')
  % 
  %   fprintf('.')
  % end
  % 
  % zoom_in = wl > 550;
  % visProd3D(wl(zoom_in), dt, ap_approx(:,zoom_in), false, 'Wavelength', false, 12+i);
  % title('ap')
  % visProd3D(wl(zoom_in), dt, cp_approx(:,zoom_in), false, 'Wavelength', false, 22+i);
  % title('cp')
  % 
  % 
  % cp_corr = cp_approx2;
  % ap_corr = ap_approx2;
end

%% Blended 3 residual Temperature And Scattering Correction using cp shape (Zaneveld 1994 proportional + NIR offset from Rottgers2013 semiempirical and Bourdin et al. in prep)
% function [ap_final, cp_final, flag_Tresidual] = ResidualTempScatterCorrZaneveld_cp_blended(ap_approx, cp_approx, wl, psi, dt)
function [ap_final, cp_final, flag_Tresidual] = ResidualTempScatterCorrSemiempirical_blended2(ap_approx, cp_approx, wl, psi, dt)
  % Function from Emmanuel Boss and Guillaume Bourdin after Zaneveld 1994 method 3 modified to use cp instead of bp and Rottgers et al. 2013: 
  % Proportional scattering correction + allows absorption in NIR in case of high NAP
  fprintf('.')
  % Find Near Infrared & references
  iNIR = 710 <= wl &  wl <= 750;  % spectral srange for optimization (710 to 750nm)
  if isempty(iNIR); error('Unable to perform correction as no wavelength available in NIR.'); end
  % Find nearest wavelength to greater than 730 nm to use as reference for correction
  iref = find(abs(wl - 715) == min(abs(wl - 715)), 1, 'first'); % 715 730
  % If ACS spectra do not go up to 715 nm take the closest wavelength to 715 nm
  if isempty(iref); [~, iref] = max(wl); end % works as there is data in iNIR so lowest wavelength is 710
  
  % Define empirical function of scattering correction at lambda reference (From Bourdin et al. 2025)
  ap715 = @(ap) 0.0831 * ap;

  % correct for residual temperature signal
  [ap_Tcorr, cp_Tcorr, flag_Tresidual] = ResidualTemperatureCorrection(wl, ap_approx, cp_approx, psi, iNIR, iref, dt, ap715(ap_approx(:,iref)), 'Semiempirical_blended2');
  fprintf('.')

  % Apply scattering correction adjusting for acceptance angle in c tube after Rottgers et al. 2013:
  % Evaluation of scatter corrections for ac-9 absorption measurements in coastal waters, DOI: 10.1016/j.mio.2013.11.001
  ap_corr = ap_Tcorr - ap_Tcorr(:,iref) .* cp_Tcorr ./ cp_Tcorr(:,iref) + ap715(ap_Tcorr(:,iref));
  fprintf('.')

  % correct again for residual temperature signal after scattering correction
  [ap_final, cp_final, flag_Tresidual2] = ResidualTemperatureCorrection(wl, ap_corr, cp_Tcorr, psi, iNIR, iref, dt, ap715(ap_approx(:,iref)), 'Semiempirical_blended2');
  flag_Tresidual(flag_Tresidual2) = true;
end

%% Blended 3 residual Temperature And Scattering Correction using bp shape (Zaneveld 1994 proportional + NIR offset from Rottgers2013 semiempirical and Bourdin et al. in prep)
% function [ap_final, cp_final, flag_Tresidual] = ResidualTempScatterCorrZaneveld_bp_blended(ap_approx, cp_approx, wl, psi, dt)
function [ap_final, cp_final, flag_Tresidual] = ResidualTempScatterCorrSemiempirical_blended3(ap_approx, cp_approx, wl, psi, dt)
  % Function from Emmanuel Boss and Guillaume Bourdin after Zaneveld 1994 method 3 and Rottgers et al. 2013: 
  % Proportional scattering correction + allows absorption in NIR in case of high NAP
  fprintf('.')
  % Find Near Infrared & references
  iNIR = 710 <= wl &  wl <= 750;  % spectral srange for optimization (710 to 750nm)
  if isempty(iNIR); error('Unable to perform correction as no wavelength available in NIR.'); end
  % Find nearest wavelength to greater than 730 nm to use as reference for correction
  iref = find(abs(wl - 715) == min(abs(wl - 715)), 1, 'first'); % 715 730
  % If ACS spectra do not go up to 715 nm take the closest wavelength to 715 nm
  if isempty(iref); [~, iref] = max(wl); end % works as there is data in iNIR so lowest wavelength is 710
  
  % Define empirical function of scattering correction at lambda reference (From Bourdin et al. 2025)
  ap715 = @(ap) 0.0831 * ap;

  % correct for residual temperature signal
  [ap_Tcorr, cp_Tcorr, flag_Tresidual] = ResidualTemperatureCorrection(wl, ap_approx, cp_approx, psi, iNIR, iref, dt, ap715(ap_approx(:,iref)), 'Semiempirical_blended3');
  fprintf('.')

  % Apply scattering correction adjusting for acceptance angle in c tube after Rottgers et al. 2013:
  % Evaluation of scatter corrections for ac-9 absorption measurements in coastal waters, DOI: 10.1016/j.mio.2013.11.001
  ec = 1/0.56; % acceptance angle correcton
  ap_corr = (1 - ap_Tcorr(:,iref) ./ (cp_Tcorr(:,iref).*ec - ap715(ap_Tcorr(:,iref)))).^(-1) .* ...
    (ap_Tcorr - ap_Tcorr(:,iref) .* cp_Tcorr.*ec ./ (cp_Tcorr(:,iref).*ec - ap715(ap_Tcorr(:,iref))) + ap715(ap_Tcorr(:,iref)));
  fprintf('.')
  
  % correct again for residual temperature signal after scattering correction
  [ap_final, cp_final, flag_Tresidual2] = ResidualTemperatureCorrection(wl, ap_corr, cp_Tcorr, psi, iNIR, iref, dt, ap715(ap_approx(:,iref)), 'Semiempirical_blended3');
  flag_Tresidual(flag_Tresidual2) = true;
end

%% Residual Temperature Correction
function [ap_Tcorr, cp_Tcorr, flag_Tresidual] = ResidualTemperatureCorrection(wl, ap, cp, psi, iNIR, iref, dt, NiRoffset, scattering_correction)
  % interpolate psiT on ACs lambda
  psiT = interp1(psi.wl, psi.psiT, wl);
  % fillmissing data when ACS wavelength < 400 nm
  psiT = fillmissing(psiT,'linear',2, 'SamplePoints',wl);
  % Parameters of minization routine
  opts = optimset('fminsearch');
  % opts = optimset(opts,'NonlEqnAlgorithm', 'gn'); % Does not work on R2017a
  opts = optimset(opts,'MaxIter',20000000);
  opts = optimset(opts,'MaxFunEvals',20000);
  opts = optimset(opts,'TolX',1e-8);
  opts = optimset(opts,'TolFun',1e-8);
  % % compute bp
  % bp = cp - ap;
  % % Select only spectra with no NaN
  % sel = find(all(isfinite(bp),2));
  % initialize flag for spectra with NaNs
  flag_Tresidual = false(size(dt));
  deltaTresiduals = NaN(size(dt));

  switch scattering_correction
    case 'Zaneveld1994_proportional'
      % Select only spectra with no NaN
      sel = find(all(isfinite(ap),2) & all(isfinite(cp),2));
      for k = sel'
        deltaTresiduals(k) = fminsearch(@costFun_RTSC_Zaneveld, 0, opts, ap(k,:), cp(k,:), psiT, iNIR, iref);
      end
    case 'Rottgers2013_semiempirical'
      % Select only spectra with no NaN
      sel = find(all(isfinite(ap),2));
      for k = sel'
        foo = fminsearch(@costFun_RTSC_Rottgers, [0 NiRoffset(k)], opts, ap(k,:), psiT, iNIR, iref);
        deltaTresiduals(k) = foo(1);
      end
    case 'Semiempirical_blended1'
      % Select only spectra with no NaN
      sel = find(all(isfinite(ap),2) & all(isfinite(cp),2));
      for k = sel'
        foo = fminsearch(@costFun_RTSC_blended1, [0 NiRoffset(k)], opts, ap(k,:), cp(k,:), psiT, iNIR, iref);
        deltaTresiduals(k) = foo(1);
      end
    case 'Semiempirical_blended2'
      % Select only spectra with no NaN
      sel = find(all(isfinite(ap),2) & all(isfinite(cp),2));
      for k = sel'
        foo = fminsearch(@costFun_RTSC_blended2, [0 NiRoffset(k)], opts, ap(k,:), cp(k,:), psiT, iNIR, iref);
        deltaTresiduals(k) = foo(1);
      end
    case 'Semiempirical_blended3'
      % compute bp
      bp = cp - ap;
      % Select only spectra with no NaN
      sel = find(all(isfinite(bp),2));
      for k = sel'
        foo = fminsearch(@costFun_RTSC_blended3, [0 NiRoffset(k)], opts, ap(k,:), cp(k,:), 1/0.56, psiT, iNIR, iref);
        deltaTresiduals(k) = foo(1);
      end
    otherwise
      error('Residual temperature and scattering correction "%s" not supported', scattering_correction)
  end

  % for k = sel'
  %   if isnan(NiRoffset(k))
  %     deltaTresiduals(k) = fminsearch(@costFun_RTSC, 0, opts, ap(k,:), bp(k,:), psiT, iNIR, iref, NiRoffset(k));
  %   else
  %     foo = fminsearch(@costFun_RTSC, [0 NiRoffset(k)], opts, ap(k,:), bp(k,:), psiT, iNIR, iref, NiRoffset(k));
  %     deltaTresiduals(k) = foo(1);
  %   end
  % end

  % Interpolate linearly deltaT in case of missing few wavelength data in spectra (up to 10 min gaps)
  deltaTresiduals = fillmissing(deltaTresiduals,'linear',1,'MaxGap',minutes(10),'SamplePoints',dt);
  % flag spectra not corrected for residual T
  idflag = true(size(flag_Tresidual));
  idflag(sel) = false;
  flag_Tresidual(idflag) = true;
  % flag_Tresidual(~all(isfinite(bp),2)) = true;
  % Apply temperature correction and replace negative values at iref by NaN to avoid complex solution leakage
  ap_Tcorr = ap - psiT .* deltaTresiduals;
  cp_Tcorr = cp - psiT .* deltaTresiduals;
end

% % Residual Temperature Correction cost function
% function cost = costFun_RTSC(x0, ap, bp, psiT, iNIR, iref, NiRoffset)
%   if isnan(NiRoffset)
%     cost = sum(abs(ap(iNIR) - psiT(iNIR).*x0(1) - ((ap(iref)-psiT(iref).*x0(1))./bp(iref)).*bp(iNIR)));
%   else
%     cost = sum(abs(ap(iNIR) - psiT(iNIR).*x0(1) - ((ap(iref)-psiT(iref).*x0(1))./bp(iref)).*bp(iNIR) + x0(2)));
%   end
% end
% % Residual Temperature Correction cost function
% function cost = costFun_RTSC(deltaT, ap, bp, psiT, iNIR, iref)
%   cost = sum(abs(ap(iNIR) - psiT(iNIR).*deltaT - ((ap(iref)-psiT(iref).*deltaT)./bp(iref)).*bp(iNIR)));
% end

  % % blended 1
  % ap_corr = ap_Tcorr - (ap_Tcorr(:,iref) - ap715(ap_Tcorr(:,iref))) .* cp_Tcorr ./ cp_Tcorr(:,iref);
  % 
  % % blended 2
  % ap_corr = ap_Tcorr - ap_Tcorr(:,iref) .* cp_Tcorr ./ cp_Tcorr(:,iref) + ap715(ap_Tcorr(:,iref));
  % 
  % % blended 3
  % ap_corr = ((1 - ap_Tcorr(:,iref) ./ (cp_Tcorr(:,iref).*ec - ap715(ap_Tcorr(:,iref))))).^(-1) .* ...
  %   (ap_Tcorr - ap_Tcorr(:,iref) .* cp_Tcorr.*ec ./ (cp_Tcorr(:,iref).*ec - ap715(ap_Tcorr(:,iref))) + ap715(ap_Tcorr(:,iref)));
  

% Residual Temperature Correction cost function (Zaneveld)
function cost = costFun_RTSC_Zaneveld(x0, ap, cp, psiT, iNIR, iref)
  % ap_final = ap_Tcorr - ap_Tcorr(:,iref) ./ bp_approx(:,iref) .* bp_approx;
  % cost = sum(abs(ap(iNIR) - psiT(iNIR).*x0 - ((ap(iref)-psiT(iref).*x0)./bp(iref)).*bp(iNIR))); old zaneveld cost function
  cost = sum(abs(ap(iNIR)-psiT(iNIR).*x0 - ((ap(iref)-psiT(iref).*x0)./(cp(iref)-(ap(iref)-psiT(iref).*x0(1))).*(cp(iNIR)-(ap(iNIR)-psiT(iNIR).*x0(1))))));
end
% Residual Temperature Correction cost function (Rottgers)
function cost = costFun_RTSC_Rottgers(x0, ap, psiT, iNIR, iref)
  % ap_final = ap_Tcorr - (ap_Tcorr_iref - ap715(ap_Tcorr(:,iref)));
  % cost = sum(abs(ap(iNIR)-psiT(iNIR).*x0(1) - ((ap(iref)-psiT(iref).*x0(1)) - (x0(2)-psiT(iref).*x0(1)))));
  % cost = sum(abs(ap(iNIR)-psiT(iNIR).*x0(1) - ((ap(iref)-psiT(iref).*x0(1)) - x0(2))));
  cost = sum(abs(ap(iNIR)-psiT(iNIR).*x0(1) - (ap(iref)-psiT(iref).*x0(1))));
end
% Residual Temperature Correction cost function (blended1)
function cost = costFun_RTSC_blended1(x0, ap, cp, psiT, iNIR, iref)
  % ap_corr = ap_Tcorr - (ap_Tcorr(:,iref) - ap715(ap_Tcorr(:,iref))) .* cp_Tcorr ./ cp_Tcorr(:,iref);
  % cost = sum(abs((ap(iNIR)-psiT(iNIR).*x0(1)) - ((ap(iref)-psiT(iref).*x0(1)) - (x0(2)-psiT(iref).*x0(1))) .* ...
  %   (cp(iNIR)-psiT(iNIR).*x0(1)) ./ (cp(iref)-psiT(iref).*x0(1))));
  % cost = sum(abs((ap(iNIR)-psiT(iNIR).*x0(1)) - ((ap(iref)-psiT(iref).*x0(1)) - x0(2)) .* ...
  %   (cp(iNIR)-psiT(iNIR).*x0(1)) ./ (cp(iref)-psiT(iref).*x0(1))));
  cost = sum(abs((ap(iNIR)-psiT(iNIR).*x0(1)) - (ap(iref)-psiT(iref).*x0(1)) .* (cp(iNIR)-psiT(iNIR).*x0(1)) ./ (cp(iref)-psiT(iref).*x0(1))));
end
% Residual Temperature Correction cost function (blended2)
function cost = costFun_RTSC_blended2(x0, ap, cp, psiT, iNIR, iref)
  % ap_corr = ap_Tcorr - ap_Tcorr(:,iref) .* cp_Tcorr ./ cp_Tcorr(:,iref) + ap715(ap_Tcorr(:,iref));
  % cost = sum(abs((ap(iNIR)-psiT(iNIR).*x0(1)) - (ap(iref)-psiT(iref).*x0(1)) .* (cp(iNIR)-psiT(iNIR).*x0(1)) ./ ...
  %   (cp(iref)-psiT(iref).*x0(1)) + (x0(2)-psiT(iref).*x0(1))));
  % cost = sum(abs((ap(iNIR)-psiT(iNIR).*x0(1)) - (ap(iref)-psiT(iref).*x0(1)) .* (cp(iNIR)-psiT(iNIR).*x0(1)) ./ ...
  %   (cp(iref)-psiT(iref).*x0(1)) + x0(2)));
  cost = sum(abs((ap(iNIR)-psiT(iNIR).*x0(1)) - (ap(iref)-psiT(iref).*x0(1)) .* (cp(iNIR)-psiT(iNIR).*x0(1)) ./ (cp(iref)-psiT(iref).*x0(1))));
end
% Residual Temperature Correction cost function (blended3)
function cost = costFun_RTSC_blended3(x0, ap, cp, ec, psiT, iNIR, iref)
  % ap_corr = ((1 - ap_Tcorr(:,iref) ./ (cp_Tcorr(:,iref).*ec - ap715(ap_Tcorr(:,iref))))).^(-1) .* ...
  %   (ap_Tcorr - ap_Tcorr(:,iref) .* cp_Tcorr.*ec ./ (cp_Tcorr(:,iref).*ec - ap715(ap_Tcorr(:,iref))) + ap715(ap_Tcorr(:,iref)));
  % cost = sum(abs(((1 - (ap(iNIR)-psiT(iNIR).*x0(1)) ./ ((cp(iNIR)-psiT(iNIR).*x0(1)).*ec - (x0(2)-psiT(iref).*x0(1))))).^(-1) .* ...
  %   ((ap(iNIR)-psiT(iNIR).*x0(1)) - (ap(iref)-psiT(iref).*x0(1)) .* (cp(iNIR)-psiT(iNIR).*x0(1)).*ec ./ ((cp(iref)-psiT(iref).*x0(1)).*ec - ...
  %   (x0(2)-psiT(iref).*x0(1))) + (x0(2)-psiT(iref).*x0(1)))));
  % cost = sum(abs(((1 - (ap(iNIR)-psiT(iNIR).*x0(1)) ./ ((cp(iNIR)-psiT(iNIR).*x0(1)).*ec - x0(2)))).^(-1) .* ...
  %   ((ap(iNIR)-psiT(iNIR).*x0(1)) - (ap(iref)-psiT(iref).*x0(1)) .* (cp(iNIR)-psiT(iNIR).*x0(1)).*ec ./ ((cp(iref)-psiT(iref).*x0(1)).*ec - ...
  %   x0(2)) + x0(2))));
  cost = sum(abs((1 - (ap(iNIR)-psiT(iNIR).*x0(1)) ./ ((cp(iNIR)-psiT(iNIR).*x0(1)).*ec)).^(-1) .* ...
    ((ap(iNIR)-psiT(iNIR).*x0(1)) - (ap(iref)-psiT(iref).*x0(1)) .* (cp(iNIR)-psiT(iNIR).*x0(1)).*ec ./ ((cp(iref)-psiT(iref).*x0(1)).*ec))));
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
    {'ap', 'cp'}) & ~contains(acs_data.Properties.VariableNames, {'_sd','_se','_unc','_n'}));
  
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
function agaus = GaussDecomp(p, lambda, compute_ad_aphi, nap_offset)
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
  p.ap_unc(:, all(isnan(p.ap),1)) = [];
  p.ap(:, all(isnan(p.ap),1)) = [];
  if nap_offset
    ap_715 = interp1(lambda, p.ap', 715, 'linear')';
    p.ap = p.ap - ap_715;
  end
  
  % Peak center values ("peak_loc") determined using a interative approach that allows the location to vary (uses the matlab
  % function LSQNONLIN), and are rounded to nearest integer. Sigma values ("lsqsig") are determined similarly. FWHM = sigma*2.355
  peak_loc = [406,434,453,470,492,523,550,584,617,638,660,675];
  lsqsig = [16,12,12,13,16,14,14,16,13,11,11,10];
  onenm = 400:1:720;
  
  fprintf('Gaussian decomposition ... ')
  
  % extrapolate NaN tails with constant nearest non-NaN value
  ap_filled = fillmissing(p.ap, 'nearest', 2, 'SamplePoints', lambda);
  ap_unc_filled = fillmissing(p.ap_unc, 'nearest', 2, 'SamplePoints', lambda);
  
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
  ap_unc_int = interp1(lambda, ap_unc_filled', onenm, 'linear', 'extrap');
  acorr2onenm_new = (acorr2onenm ./ ap_unc_int);
  
  amps = NaN(size(p,1), size(coef2,2));
  sumspec_temp = NaN(size(coef2,1), size(p,1));
  % compspec_temp = NaN(size(coef2,1), size(coef2,2), size(p,1));
  
  parfor i = 1:size(acorr2onenm_new,2)
    coef2_new = coef2 ./ ap_unc_int(:,i);
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

function [filt_interp, filt, fit_unc] = agcg_fdom_interpolation(filt_interp, filt, lambda, filt_avg)
  % Compute da, dc, and dfdom
  % da = diff(filt_avg.a,[],1);
  % dc = diff(filt_avg.c,[],1);
  % dfdom = diff(filt_avg.fdom,[],1);
  da = diff(filt.a,[],1);
  dc = diff(filt.c,[],1);
  dfdom = diff(filt.fdom,[],1);
  % keep only data with high a, c, and fdom dynamic range to minimize impact of drift and biofouling on slopes
  ida = abs(lambda.a - 420) == min(abs(lambda.a - 420));
  idc = abs(lambda.c - 420) == min(abs(lambda.c - 420));
  % keep da and dc only when da|dc > 0.001 m^(-1) and dfdom > 0.001
  da(da(:, ida) < 0.001 | dfdom < 0.001, :) = NaN;
  dc(dc(:, idc) < 0.001 | dfdom < 0.001, :) = NaN;
  % keep dfdom when dfdom > 0.001 volts (instrument resolution after binning)
  dfdom(dfdom < 0.001) = NaN;

  % figure; scatter(dfdom, da(:,ida), 20, 'filled'); set(gca, 'XScale', 'log', 'YScale', 'log')
  % figure; histogram(da(:,ida))
  % figure; histogram(dfdom)
  
  % Two methods to compute the slopes giving comparable results
  % % METHOD 1: Robust linear regression between da and dfdom and dc and dfdom
  % slope_da_dfdom = NaN(size(lambda.a));
  % slope_dc_dfdom = NaN(size(lambda.a));
  % stats = struct();
  % for i = 1:size(lambda.a,2)
  %   [stats.b, stats.stats] = robustfit(dfdom, da(:,i));
  %   sslope_da_dfdom(i) = stats.b(2);
  %   [stats.b, stats.stats] = robustfit(dfdom, dc(:,i));
  %   slope_dc_dfdom(i) = stats.b(2);
  % end

  % METHOD 2: median or mean of percentiles of da/dfdom and dc/dfdom ratios
  da_dfdom = da ./ dfdom;
  dc_dfdom = dc ./ dfdom;
  % % Using 30th and 25th percentiles
  % slope_da_dfdom = prctile(da_dfdom, 30, 1);
  % slope_dc_dfdom = prctile(dc_dfdom, 25, 1);
  % Using the mean of 0 to 50th percentiles of da_dfdom and 0 to 40th percentiles of dc_dfdom
  ida_avg = da_dfdom > prctile(da_dfdom, 0, 1) & da_dfdom < prctile(da_dfdom, 50, 1);
  idc_avg = dc_dfdom > prctile(dc_dfdom, 0, 1) & dc_dfdom < prctile(dc_dfdom, 40, 1);
  % Test including everything
  % ida_avg = da_dfdom >= prctile(da_dfdom, 0, 1) & da_dfdom <= prctile(da_dfdom, 100, 1);
  % idc_avg = dc_dfdom >= prctile(dc_dfdom, 0, 1) & dc_dfdom <= prctile(dc_dfdom, 100, 1);
  da_dfdom(~ida_avg) = NaN;
  dc_dfdom(~idc_avg) = NaN;
  slope_da_dfdom = mean(da_dfdom, 1, 'omitmissing');
  slope_dc_dfdom = mean(dc_dfdom, 1, 'omitmissing');


  % Fit exponential between 410 and 550 nm
  fprintf('Exponential fit on da/dfdom and dc/dfdom ... \n')
  sel_a = lambda.a > 405 & lambda.a < 550;
  sel_c = lambda.c > 405 & lambda.c < 550;
  [da_dfdom_intercp_fit, da_dfdom_base_fit, fit_unc.SSE_a, fit_unc.MSE_a, fit_unc.RMSE_a, fit_unc.nRMSE_a] = FitExp(lambda.a(sel_a), slope_da_dfdom(sel_a), []);
  [dc_dfdom_intercp_fit, dc_dfdom_base_fit, fit_unc.SSE_c, fit_unc.MSE_c, fit_unc.RMSE_c, fit_unc.nRMSE_c] = FitExp(lambda.c(sel_c), slope_dc_dfdom(sel_c), []);
  fprintf('Done\n')
  
  % rebuild da/dfdom and dc/dfdom on full spectra
  expfun = @(p, xd) p(1) * exp(p(2) * (xd - 440));
  da_dfdom_fit = expfun([da_dfdom_intercp_fit da_dfdom_base_fit], lambda.a);
  dc_dfdom_fit = expfun([dc_dfdom_intercp_fit dc_dfdom_base_fit], lambda.c);
  
  % populate filt_interp database
  filt_interp.slope_interp_a = repmat(da_dfdom_fit, size(filt_interp,1), 1);
  filt_interp.slope_interp_c = repmat(dc_dfdom_fit, size(filt_interp,1), 1);
  
  % plot da/dfdom and dc/dfdom and Tara Europa's bentchop spectrophotometer slopes
  % TE_param = readtable('/Volumes/Extreme SSD/TaraEuropa/DeviceFiles/parameter_to_convert_FDOM_to_ag_TaraEuropa_2023_2024.csv');
  TE_param = table();
  TE_param.wl = [700	699	698	697	696	695	694	693	692	691	690	689	688	687	686	685	684	683	682	681	680	679	678	677	676	675	674	673	672	671	670	669	668	667	666	665	664	663	662	661	660	659	658	657	656	655	654	653	652	651	650	649	648	647	646	645	644	643	642	641	640	639	638	637	636	635	634	633	632	631	630	629	628	627	626	625	624	623	622	621	620	619	618	617	616	615	614	613	612	611	610	609	608	607	606	605	604	603	602	601	600	599	598	597	596	595	594	593	592	591	590	589	588	587	586	585	584	583	582	581	580	579	578	577	576	575	574	573	572	571	570	569	568	567	566	565	564	563	562	561	560	559	558	557	556	555	554	553	552	551	550	549	548	547	546	545	544	543	542	541	540	539	538	537	536	535	534	533	532	531	530	529	528	527	526	525	524	523	522	521	520	519	518	517	516	515	514	513	512	511	510	509	508	507	506	505	504	503	502	501	500	499	498	497	496	495	494	493	492	491	490	489	488	487	486	485	484	483	482	481	480	479	478	477	476	475	474	473	472	471	470	469	468	467	466	465	464	463	462	461	460	459	458	457	456	455	454	453	452	451	450	449	448	447	446	445	444	443	442	441	440	439	438	437	436	435	434	433	432	431	430	429	428	427	426	425	424	423	422	421	420	419	418	417	416	415	414	413	412	411	410	409	408	407	406	405	404	403	402	401	400	399	398	397	396	395	394	393	392	391	390	389	388	387	386	385	384	383	382	381	380	379	378	377	376	375	374	373	372	371	370	369	368	367	366	365	364	363	362	361	360	359	358	357	356	355	354	353	352	351	350	349	348	347	346	345	344	343	342	341	340	339	338	337	336	335	334	333	332	331	330	329	328	327	326	325	324	323	322	321	320	319	318	317	316	315	314	313	312	311	310	309	308	307	306	305	304	303	302	301	300	299	298	297	296	295	294	293	292	291	290	289	288	287	286	285	284	283	282	281	280	279	278	277	276	275	274	273	272	271	270	269	268	267	266	265	264	263	262	261	260	259	258	257	256	255	254	253	252	251	250	249	248	247	246	245	244	243	242	241	240	239	238	237	236	235	234	233	232	231	230];
  TE_param.intercept = [-0.00062712	-0.000168307	0.000440015	0.001265337	0.001781034	0.001904382	0.002541238	0.002749404	0.003667174	0.003390618	0.004065411	0.005169505	0.005265082	0.005492812	0.005466388	0.006032192	0.005945831	0.006499932	0.006550342	0.007022579	0.007389955	0.00762924	0.008648955	0.008010777	0.009130219	0.00947055	0.009273344	0.009114811	0.009373714	0.00977936	0.008610542	0.010250619	0.009696448	0.009557338	0.010174578	0.009961212	0.00928291	0.009950914	0.009477283	0.009066386	0.008845121	0.008691017	0.008863094	0.00854934	0.009526717	0.009744255	0.01012761	0.010196946	0.010572502	0.01017496	0.010983867	0.010619283	0.01160285	0.011757633	0.011163521	0.011559563	0.011618443	0.011349325	0.012424583	0.012124843	0.012222371	0.012021504	0.011872705	0.012365508	0.012056551	0.012666259	0.012394772	0.012128165	0.01343131	0.012256362	0.012778959	0.012968785	0.013512474	0.013219101	0.013057791	0.013132307	0.013332956	0.013510574	0.013483958	0.013383176	0.012906849	0.013647882	0.013422855	0.013365735	0.01355293	0.013579312	0.01341752	0.01391735	0.013357953	0.013084092	0.013903311	0.013664796	0.013077641	0.012804929	0.011993884	0.012079912	0.011005249	0.010551916	0.010482686	0.010135502	0.010243765	0.010094909	0.011394044	0.011688503	0.012121652	0.012065581	0.012088836	0.012506222	0.012486244	0.01280454	0.013462966	0.012883767	0.012962375	0.013379794	0.013534126	0.013043624	0.013876131	0.013142417	0.013917719	0.014253176	0.013946659	0.014662849	0.014613984	0.014952453	0.014820813	0.015073096	0.015404107	0.015401313	0.015226158	0.015439163	0.015956998	0.016069413	0.016323638	0.016160639	0.01660261	0.016479347	0.016299284	0.016836126	0.017213198	0.017023241	0.016838795	0.017482468	0.017082995	0.017450211	0.017682585	0.017403935	0.017924902	0.018001611	0.017947614	0.017456228	0.01780269	0.018027378	0.01829827	0.018017295	0.018268537	0.018446308	0.018410103	0.018671339	0.018729986	0.018497864	0.019152654	0.019359972	0.019074855	0.019309673	0.019504451	0.01942247	0.019878786	0.019571953	0.019655931	0.020210295	0.020090269	0.020017801	0.019230525	0.020221958	0.020456805	0.020484032	0.020833639	0.020862527	0.020932246	0.020989409	0.020835551	0.021110026	0.021221561	0.020875341	0.021043114	0.021481421	0.021530241	0.021568662	0.021874617	0.021906158	0.021959576	0.022188692	0.022240884	0.021994048	0.022670711	0.022976957	0.022989854	0.022878742	0.023523613	0.023397784	0.023602925	0.023600232	0.023910274	0.023574739	0.023769148	0.023638144	0.024416777	0.024238039	0.024675273	0.025039616	0.025108509	0.025412944	0.025610739	0.025520092	0.02586117	0.026250689	0.025736981	0.025879902	0.026437935	0.026491274	0.026908571	0.026341293	0.026985855	0.027203561	0.027403335	0.027531668	0.027497542	0.027937995	0.028144751	0.028252665	0.028279002	0.028496017	0.028387903	0.028912802	0.028726303	0.029900889	0.029602425	0.029901441	0.029984891	0.030291919	0.030499672	0.030827532	0.030301554	0.030583337	0.031154395	0.031476151	0.031700286	0.031891874	0.032116568	0.032162349	0.032897762	0.033161288	0.032793596	0.033149607	0.033710346	0.034044092	0.033973224	0.034279939	0.034172989	0.034972436	0.035053755	0.035533228	0.035948434	0.036564765	0.036019826	0.036662648	0.037354278	0.037158153	0.037181576	0.037726603	0.038521029	0.03825438	0.039132886	0.03938187	0.039096152	0.040253772	0.04029499	0.039786767	0.040587436	0.041190321	0.041314472	0.041421212	0.042305245	0.042955864	0.043099716	0.042856019	0.042987514	0.043601079	0.044183466	0.044007919	0.044083966	0.042895616	0.045025979	0.044492889	0.044011646	0.043638426	0.041906282	0.041264298	0.041553778	0.042528853	0.042441972	0.042677435	0.04279519	0.042268425	0.04237839	0.042321008	0.042076479	0.04277613	0.043162114	0.042730409	0.043391878	0.043253598	0.044377595	0.043243488	0.04470216	0.045267275	0.044773245	0.045002259	0.044558273	0.04507876	0.046634984	0.046534815	0.048093783	0.047410288	0.048719589	0.048504105	0.049907362	0.049988885	0.051113235	0.051086588	0.051565192	0.052067603	0.052249338	0.053660257	0.054415937	0.055043461	0.055529795	0.055231897	0.057242044	0.056163875	0.056271329	0.05652474	0.057967872	0.056923599	0.058481769	0.076827366	0.078738433	0.081615104	0.083227248	0.087304381	0.086734674	0.091293691	0.094100948	0.094430326	0.098079411	0.100394298	0.103044563	0.106364963	0.108737286	0.112965999	0.117275047	0.119918652	0.122158204	0.125382749	0.130470433	0.133592206	0.135983203	0.14114225	0.142560269	0.146693883	0.150567828	0.071450201	0.075510385	0.076685355	0.079670803	0.083289946	0.086998241	0.088651933	0.090694415	0.094357534	0.095293659	0.101207192	0.104442744	0.10955638	0.114625824	0.118973372	0.11888609	0.124694978	0.127320404	0.133701542	0.137405794	0.143932443	0.147368478	0.155602966	0.1617951	0.168755232	0.171146676	0.141152682	0.150091431	0.158947522	0.166393207	0.183042034	0.194317161	0.207443395	0.218858235	0.230311822	0.241286717	0.258013581	0.274162347	0.289437797	0.304403192	0.323899054	0.340362579	0.382613771	0.408817793	0.433319236	0.455145065	0.4768694	0.49833128	0.517912829	0.539166345	0.560823142	0.582547386	0.601762456	0.622397145	0.642917517	0.663340435	0.683988665	0.704601166	0.725528801	0.744132679	0.808043619	0.83336855	0.857915143	0.887314594	1.186381322	1.209490261	1.233123512	1.258818592	1.281501095	1.305294443	1.331191046	1.349295298	1.372678714	1.398284889	1.42510988	1.451967831	1.480450893	1.510288414	1.544778959	1.579203947	1.620037166	1.660605212	1.705616337	1.758136553	1.812372805	1.874312661	1.943272876	2.009907912	2.086396044	2.171070023	2.265519977	2.36843692	2.469393762	2.579110698	2.642322941	2.779903648	2.880163205	3.078941949	3.325339336	3.631728653];
  TE_param.slope = [0.013963641	0.014342945	0.013906888	0.01469178	0.013978119	0.013821395	0.01439463	0.014989059	0.013726671	0.015471164	0.015842676	0.013068498	0.014190455	0.014433903	0.015228531	0.015359335	0.017745874	0.016570844	0.017367728	0.01679041	0.016596225	0.016858289	0.015011722	0.016542247	0.014226081	0.013961759	0.015619655	0.01755128	0.016951385	0.016397343	0.01943629	0.017823232	0.018068392	0.018713674	0.016213778	0.01850826	0.021147478	0.018322398	0.01910029	0.019219676	0.0222295	0.020662593	0.022193354	0.024422261	0.022497305	0.025569537	0.025686953	0.027136429	0.027333598	0.0277238	0.027269141	0.030090043	0.028627744	0.028861175	0.030386317	0.029756702	0.029212247	0.031721604	0.029290137	0.031156692	0.032109128	0.033415941	0.033150544	0.034292896	0.034907512	0.034793796	0.034403338	0.037112291	0.032399309	0.036857207	0.035935233	0.037799234	0.035946779	0.038443823	0.039437381	0.039364821	0.039187085	0.038520478	0.042127595	0.040884744	0.044314124	0.043347715	0.043031069	0.044803268	0.044409438	0.046204674	0.047303261	0.047331302	0.047859599	0.049769951	0.047112722	0.049910795	0.051691958	0.051945912	0.054144644	0.054597708	0.057378383	0.057857168	0.058338268	0.059813294	0.060941249	0.062594288	0.062077757	0.062815866	0.06416546	0.064996872	0.066131435	0.067154401	0.068147175	0.067908118	0.065905578	0.070723063	0.072088145	0.071025129	0.071588059	0.074555601	0.07319133	0.076504361	0.074797157	0.075634235	0.078592095	0.077740168	0.079222851	0.079786176	0.08132405	0.082377992	0.082469195	0.085310618	0.085331707	0.088251294	0.086857496	0.088662002	0.089470404	0.090185564	0.091339681	0.0935923	0.094678382	0.094433263	0.095368786	0.098482134	0.099364051	0.0995197	0.101350791	0.101984963	0.103604828	0.105541843	0.106036466	0.107513412	0.108757887	0.111605572	0.111557759	0.113608292	0.114668332	0.117023502	0.116982764	0.118118715	0.121165806	0.122014895	0.123338455	0.125094212	0.126911359	0.127084299	0.130321388	0.132220768	0.1335815	0.135554355	0.13492174	0.139976501	0.141698688	0.142393717	0.14475249	0.146389506	0.150997731	0.148244308	0.152107721	0.154325776	0.154919641	0.157686416	0.161676699	0.162303689	0.164709503	0.166479228	0.169343754	0.170967593	0.173664778	0.174293267	0.177304512	0.177880099	0.180581734	0.185343545	0.185703588	0.187918433	0.190957074	0.195527165	0.196642437	0.197297316	0.199822445	0.203408007	0.206296256	0.208930445	0.212052126	0.215293194	0.217199228	0.219311191	0.22329917	0.226906248	0.22872212	0.233339946	0.23500989	0.237545811	0.240866854	0.244838302	0.24884399	0.250232284	0.25559344	0.257682634	0.262929783	0.266288268	0.270591032	0.272739429	0.278153342	0.284780958	0.288093266	0.291936844	0.29549467	0.300501926	0.305254297	0.309608007	0.314886255	0.319743295	0.324985545	0.330088726	0.335407915	0.341059571	0.345771414	0.350613018	0.358250691	0.363742573	0.369954517	0.375460596	0.381495439	0.387741575	0.394485977	0.401247859	0.407615284	0.413499334	0.421132901	0.428766606	0.435818058	0.444687427	0.450579116	0.458322139	0.467613292	0.474855747	0.482748057	0.491255304	0.501179309	0.508933606	0.519679269	0.528743979	0.538077369	0.549105311	0.558657604	0.568696015	0.580467522	0.591315658	0.602101609	0.614377234	0.626005749	0.63676594	0.648919302	0.661838657	0.674425089	0.685536782	0.700162062	0.710961373	0.724703967	0.738720613	0.753096627	0.766392838	0.780969464	0.795945125	0.810252037	0.824960643	0.842947498	0.858013005	0.879202467	0.891218743	0.908322817	0.927720673	0.945993529	0.978561822	0.980292675	1.007886082	1.034766914	1.061212528	1.114327243	1.138565094	1.163083689	1.182665877	1.204009969	1.233849171	1.257476592	1.282125907	1.308202254	1.339335656	1.366447014	1.395542772	1.420979028	1.450367158	1.480141059	1.510658355	1.539109871	1.57011345	1.603408919	1.632472373	1.666832669	1.695510649	1.729985107	1.763628659	1.796870706	1.833124723	1.866357894	1.902086152	1.937718349	1.975290499	2.012237061	2.050961449	2.089614192	2.132971724	2.169679524	2.211239918	2.258357827	2.302843677	2.339072455	2.387946761	2.436771476	2.489389547	2.530753716	2.584064714	2.632939628	2.692291111	2.742616143	2.802352713	2.865422232	2.70241995	2.737733419	2.771257102	2.820036298	2.84474157	2.93745649	2.967665361	3.01223245	3.083825853	3.140555295	3.179408208	3.244088428	3.29639754	3.351924466	3.412744727	3.467799503	3.532092262	3.600195467	3.668925061	3.732406862	3.792603037	3.875942197	3.941013595	4.0176498	4.09281159	4.164334767	5.113363565	5.219787612	5.315087868	5.418303304	5.524444511	5.632108556	5.740339165	5.84895858	5.960095609	6.074270423	6.190511905	6.308340044	6.426516823	6.551596293	6.673383374	6.799267055	6.93049489	7.067193819	7.198660744	7.339961501	7.478150094	7.624556247	7.764982202	7.924549071	8.065557529	8.247112158	8.657700835	8.847383397	9.035759688	9.219699083	9.382904237	9.566291179	9.742034058	9.933674472	10.1288433	10.34658365	10.52498149	10.7317399	10.94081729	11.16694334	11.36085956	11.59230093	11.47208169	11.61902904	11.7461494	11.86988482	12.04275433	12.24141351	12.40413694	12.58228576	12.73535768	12.9427355	13.08116149	13.22986386	13.42408394	13.56213483	13.77577653	13.96647327	14.16315178	14.37402422	13.75094368	13.80628363	13.83990255	13.81417076	10.32054043	10.41958508	10.54194863	10.65193904	10.76951607	10.88040343	10.99758304	11.09919226	11.20274984	11.31381259	11.4258159	11.54103986	11.65201808	11.76760519	11.88566535	12.00559499	12.13330516	12.26991393	12.40774103	12.54513213	12.70569124	12.86364812	13.04926016	13.24373758	13.46486094	13.70860484	13.98291038	14.31802582	14.72472411	15.26240453	16.48046381	17.23783388	18.45567704	19.12545123	19.92350948	20.6215564];
  poo = [];
  figure(36); clf; hold on
  poo(1) = plot(TE_param.wl(TE_param.wl > 400), TE_param.slope(TE_param.wl > 400),'-g','LineWidth',1.5);
  poo(2) = errshaded(lambda.a, slope_da_dfdom, std(da_dfdom, 'omitmissing')./sum(~isnan(da_dfdom),1), 'b', 0.2, '--', 1);
  poo(3) = errshaded(lambda.a, slope_dc_dfdom, std(dc_dfdom, 'omitmissing')./sum(~isnan(dc_dfdom),1), 'r', 0.2, '--', 1);
  poo(4) = plot(lambda.a, da_dfdom_fit, '-b','LineWidth',1.5);
  poo(5) = plot(lambda.c, dc_dfdom_fit, '-r','LineWidth',1.5);
  xlabel('\lambda')
  ylabel('\Deltaa_g/\DeltaF_{DOM} and \Deltac_g/\DeltaF_{DOM}')
  legend(poo, 'TaraEuropa benchtop \Deltaa_g/\DeltaF_{DOM}','ACS mean(prctile(\Deltaa_g/\DeltaF_{DOM}))','ACS mean(prctile(\Deltac_g/\DeltaF_{DOM}))','fit ACS \Deltaa_g/\DeltaF_{DOM}','fit ACS \Deltac_g/\DeltaF_{DOM}')

  % linear interpolation of a, c, and fdom
  filt_interp.lin_a = interp1(filt_avg.dt, filt_avg.a, filt_interp.dt, 'linear');
  filt_interp.lin_c = interp1(filt_avg.dt, filt_avg.c, filt_interp.dt, 'linear');
  filt_interp.lin_fdom = interp1(filt_avg.dt, filt_avg.fdom, filt_interp.dt, 'linear');
  filt_interp.lin_a = fillmissing(filt_interp.lin_a, 'nearest', 'SamplePoints', filt_interp.dt);
  filt_interp.lin_c = fillmissing(filt_interp.lin_c, 'nearest', 'SamplePoints', filt_interp.dt);
  filt_interp.lin_fdom = fillmissing(filt_interp.lin_fdom, 'nearest', 'SamplePoints', filt_interp.dt);
  
  % create linear interpolation flag true as default
  filt_interp.flag_linear_interp_a = true(size(filt_interp.dt));
  filt_interp.flag_linear_interp_c = true(size(filt_interp.dt));
  % set linear interpolation flag false when fdom is available
  filt_interp.flag_linear_interp_a(~isnan(filt_interp.fdom)) = false;
  filt_interp.flag_linear_interp_c(~isnan(filt_interp.fdom)) = false;
  % set linear interpolation flag true when median(da/dfdom) < 0 or median(dc/dfdom) < 0
  filt_interp.flag_linear_interp_a(filt_interp.slope_interp_a < 0) = true;
  filt_interp.flag_linear_interp_c(filt_interp.slope_interp_c < 0) = true;

  % compute ag and cg with fCDOM interpolation
  filt_interp.a(~filt_interp.flag_linear_interp_a, :) = filt_interp.lin_a(~filt_interp.flag_linear_interp_a, :) + ...
    (filt_interp.fdom(~filt_interp.flag_linear_interp_a) - filt_interp.lin_fdom(~filt_interp.flag_linear_interp_a)) .* ...
    filt_interp.slope_interp_a(~filt_interp.flag_linear_interp_a, :);
  filt_interp.c(~filt_interp.flag_linear_interp_c, :) = filt_interp.lin_c(~filt_interp.flag_linear_interp_c, :) + ...
    (filt_interp.fdom(~filt_interp.flag_linear_interp_c) - filt_interp.lin_fdom(~filt_interp.flag_linear_interp_c)) .* ...
    filt_interp.slope_interp_c(~filt_interp.flag_linear_interp_c, :);

  % compute uncertainty associated with FDOM interpolation
  filt_interp.ag_interp_unc = std(da_dfdom, 'omitmissing')./sqrt(sum(~isnan(da_dfdom),1)) .* abs(filt_interp.fdom - filt_interp.lin_fdom) + ...
    da_dfdom_fit .* filt_interp.fdom_sd ./ sqrt(filt_interp.fdom_n) .* sqrt(2);
  filt_interp.cg_interp_unc = std(dc_dfdom, 'omitmissing')./sqrt(sum(~isnan(dc_dfdom),1)) .* abs(filt_interp.fdom - filt_interp.lin_fdom) + ...
    da_dfdom_fit .* filt_interp.fdom_sd ./ sqrt(filt_interp.fdom_n) .* sqrt(2);

  % compute ag and cg with linear interpolation when fCDOM does not work merge filt_avg rows into 
  % filt_interp to linearly interpolate ag and cg only for missing data (in case fCDOM interpolation 
  % was available for half of filter event)
  append_filt_interp = filt_interp(end-size(filt_avg, 1)+1:end, :);
  var = filt_interp.Properties.VariableNames';
  var2mv = {'dt','fdom','t','s','a','c'};
  for i = 1:size(var, 1)
    if any(strcmp(var{i}, var2mv))
      append_filt_interp.(var{i}) = filt_avg.(var{i});
    else
      if isdatetime(append_filt_interp.(var{i}))
        append_filt_interp.(var{i}) = NaT(size(append_filt_interp.(var{i})));
      elseif islogical(append_filt_interp.(var{i}))
        append_filt_interp.(var{i}) = false(size(append_filt_interp.(var{i})));
      elseif isnumeric(append_filt_interp.(var{i}))
        append_filt_interp.(var{i}) = NaN(size(append_filt_interp.(var{i})));
      end
    end
  end
  % remove 1 second from appended filt_avg.dt to prevent duplicates
  append_filt_interp.dt = append_filt_interp.dt - seconds(1);
  % check if any duplicates remaining, if yes, remove 1 second
  iddup = ismember(append_filt_interp.dt, filt_interp.dt);
  while any(iddup)
    append_filt_interp.dt(iddup) = append_filt_interp.dt(iddup) - seconds(1);
    iddup = ismember(append_filt_interp.dt, filt_interp.dt);
  end
  filt_interp = [filt_interp; append_filt_interp];

  filt_interp = sortrows(filt_interp, 'dt');
  % linear interpolation when fdom interpolation was not used
  filt_interp.a = fillmissing(filt_interp.a, 'linear', 'SamplePoints', filt_interp.dt, 'EndValues', 'nearest');
  filt_interp.c = fillmissing(filt_interp.c, 'linear', 'SamplePoints', filt_interp.dt, 'EndValues', 'nearest');
  % remove filt_avg rows added for filling missing linear interpolation
  filt_interp(ismember(filt_interp.dt, append_filt_interp.dt), :) = [];
end

%%
% function [filt_interp, filt, regress_stats] = agcg_fdom_interpolation(tot, filt_interp, filt, lambda, ...
%   filt_avg, min_nb_pts_per_cluster, time_weight_for_cluster, cluster_bool)
%   % Reconstruct ag and cg from fCDOM data
%   [filt_avg, filt_interp, filt, regress_stats] = fdom_agcg_model(filt, filt_interp, lambda, filt_avg, ...
%     min_nb_pts_per_cluster, time_weight_for_cluster, cluster_bool);
% 
%   n_periods = size(filt_avg,1)-1;
% 
%   % visProd3D(lambda.a, filt_interp.dt, filt_interp.intercept_interp_a, false, 'Wavelength', false, 27);
%   % title('intercept a')
%   % visProd3D(lambda.a, filt_interp.dt, filt_interp.intercept_interp_c, false, 'Wavelength', false, 28);
%   % title('intercept c')
%   % 
%   % visProd3D(lambda.a, filt_interp.dt, filt_interp.slope_interp_a, false, 'Wavelength', false, 27);
%   % title('slope a')
%   % visProd3D(lambda.a, filt_interp.dt, filt_interp.slope_interp_c, false, 'Wavelength', false, 28);
%   % title('slope c')
% 
%   % % interpolate based on fCDOM: Emmanuel's method (local)
%   % % allocate fdom0, a0, and c0 varialbes
%   % filt_interp.fdom0 = NaN(size(filt_interp.fdom));
%   % filt_interp.a0 = NaN(size(filt_interp.a));
%   % filt_interp.c0 = NaN(size(filt_interp.c));
%   % 
%   % % interpolate based on fCDOM: Emmanuel's method
%   % delta_acf = table(NaN(n_periods, 1), 'VariableNames', {'dt'});
%   % delta_acf.delta_af = NaN(n_periods, size(filt_avg.a, 2));
%   % delta_acf.delta_cf = NaN(n_periods, size(filt_avg.a, 2));
%   % % Prepare coefficients for each total period going from t0 to t1, starting and finishing by a filtered time
%   % for i=1:n_periods
%   %   it0 = i; it1 = i + 1;
%   %   it = filt_avg.dt(it0) <= filt_interp.dt & filt_interp.dt <= filt_avg.dt(it1);
%   %   if any(it)
%   %     % interpolation based on fCDOM if enough dynamic range in fCDOM
%   %     if any((filt_interp.fdom(it) - cdom_dark) / cdom_dark > 0.05)
%   %       % prepare delta a / delta f and delta c / delta f
%   %       delta_acf.dt(i) = mean([filt_avg.dt(it0) filt_avg.dt(it1)], 'omitnan');
%   %       delta_acf.delta_af(i, :) = (filt_avg.a(it1,:) - filt_avg.a(it0,:)) / (filt_avg.fdom(it1) - filt_avg.fdom(it0));
%   %       delta_acf.delta_cf(i, :) = (filt_avg.c(it1,:) - filt_avg.c(it0,:)) / (filt_avg.fdom(it1) - filt_avg.fdom(it0));
%   %       filt_interp.fdom0(it) = repmat(filt_avg.fdom(it0), sum(it), 1);
%   %       filt_interp.a0(it,:) = repmat(filt_avg.a(it0,:), sum(it), 1);
%   %       filt_interp.c0(it,:) = repmat(filt_avg.c(it0,:), sum(it), 1);
%   %       filt_interp.flag_linear_interp(it) = false;
%   %     else
%   %       filt_interp.flag_linear_interp(it) = true;
%   %     end
%   %   end
%   % end
%   % delta_acf(all(isnan(delta_acf.delta_af), 2) | all(isnan(delta_acf.delta_cf), 2), :) = [];
%   % % interpolated delta a / delta fdom over total events
%   % filt_interp.delta_af = interp_extrap(delta_acf, filt_interp.dt, 'delta_af', [], true, 'linear', 'nearest');
%   % filt_interp.delta_cf = interp_extrap(delta_acf, filt_interp.dt, 'delta_cf', [], true, 'linear', 'nearest');
%   % % extrapolate a0, c0, fdom0
%   % filt_interp.fdom0 = fillmissing(filt_interp.fdom0, 'nearest', 'SamplePoints', filt_interp.dt);
%   % filt_interp.a0 = fillmissing(filt_interp.a0, 'nearest', 'SamplePoints', filt_interp.dt);
%   % filt_interp.c0 = fillmissing(filt_interp.c0, 'nearest', 'SamplePoints', filt_interp.dt);
%   % % compute ag and cg
%   % filt_interp.a = filt_interp.al + filt_interp.delta_af .* (filt_interp.fdom - filt_interp.fdoml);
%   % filt_interp.c = filt_interp.cl + filt_interp.delta_cf .* (filt_interp.fdom - filt_interp.fdoml);
%   % % linearly interpolate ag and cg
%   % filt_interp.lin_a = interp1(filt_avg.dt, filt_avg.a, filt_interp.dt, 'linear');
%   % filt_interp.lin_c = interp1(filt_avg.dt, filt_avg.c, filt_interp.dt, 'linear');
%   % % replace fdom interpolation by linear for filter events where fdom is constant
%   % filt_interp.a(filt_interp.flag_linear_interp) = filt_interp.lin_a(filt_interp.flag_linear_interp);
%   % filt_interp.c(filt_interp.flag_linear_interp) = filt_interp.lin_c(filt_interp.flag_linear_interp);
% 
%   % visProd3D(lambda.a, filt_interp.dt, filt_interp.a, false, 'Wavelength', false, 27);
%   % title('ag fdom interpolated')
%   % visProd3D(lambda.a, filt_interp.dt, filt_interp.lin_a, false, 'Wavelength', false, 28);
%   % title('ag linearly interpolated')
%   % 
%   % visProd3D(lambda.c, filt_interp.dt, filt_interp.c, false, 'Wavelength', false, 29);
%   % title('cg fdom interpolated')
%   % visProd3D(lambda.c, filt_interp.dt, filt_interp.lin_c, false, 'Wavelength', false, 30);
%   % title('cg linearly interpolated')
% 
%   % interpolate ag and cg based on fdom: Guillaume's method (global)
%   % allocate new variables
%   filt_interp.fCDOM_mix_a_cluster0 = false(size(filt_interp.dt));
%   filt_interp.fCDOM_mix_a_cluster1 = false(size(filt_interp.dt));
%   filt_interp.fCDOM_mix_c_cluster0 = false(size(filt_interp.dt));
%   filt_interp.fCDOM_mix_c_cluster1 = false(size(filt_interp.dt));
%   filt_interp.fCDOM_a_cluster_chg = false(size(filt_interp.dt));
%   filt_interp.fCDOM_c_cluster_chg = false(size(filt_interp.dt));
%   filt_interp.flag_linear_interp_a = false(size(filt_interp.dt));
%   filt_interp.flag_linear_interp_c = false(size(filt_interp.dt));
%   filt_interp.flag_a_filt0_not_clustered = false(size(filt_interp.dt));
%   filt_interp.flag_a_filt1_not_clustered = false(size(filt_interp.dt));
%   filt_interp.flag_c_filt0_not_clustered = false(size(filt_interp.dt));
%   filt_interp.flag_c_filt1_not_clustered = false(size(filt_interp.dt));
% 
%   %%%%% RE-WRITTEN ON 2023-12-29 and RE-WRITTEN AGAIN ON 2024-11-23
%   % linear interpolation of a, c, and fdom
%   filt_interp.lin_a = interp1(filt_avg.dt, filt_avg.a, filt_interp.dt, 'linear');
%   filt_interp.lin_c = interp1(filt_avg.dt, filt_avg.c, filt_interp.dt, 'linear');
%   filt_interp.lin_fdom = interp1(filt_avg.dt, filt_avg.fdom, filt_interp.dt, 'linear');
%   filt_interp.lin_a = fillmissing(filt_interp.lin_a, 'nearest', 'SamplePoints', filt_interp.dt);
%   filt_interp.lin_c = fillmissing(filt_interp.lin_c, 'nearest', 'SamplePoints', filt_interp.dt);
%   filt_interp.lin_fdom = fillmissing(filt_interp.lin_fdom, 'nearest', 'SamplePoints', filt_interp.dt);
% 
%   for i=1:n_periods
%     it0 = i; it1 = i + 1;
%     it = filt_avg.dt(it0) <= filt_interp.dt & filt_interp.dt <= filt_avg.dt(it1);
%     if any(it)
%       % set flags for fCDOM interpolation
%       if any(~isnan(filt_interp.fdom(it)))
%         % flag when mixed a clusters within filter event 0 or filter event 1
%         idnonnan = any(~isnan(tot.a), 2);
%         if any(it & idnonnan)
%           if filt_avg.mix_cluster_a_flag(it0)
%             % fprintf('a filter interpolation #%i flagged: mixed cluster within start filter event\n', i)
%             filt_interp.fCDOM_mix_a_cluster0(it & idnonnan) = true;
%           end
%           if filt_avg.mix_cluster_a_flag(it1)
%             % fprintf('a filter interpolation #%i flagged: mixed cluster within end filter event\n', i)
%             filt_interp.fCDOM_mix_a_cluster1(it & idnonnan) = true;
%           end
%           % flag when a clusters change between filter event 0 and filter event 1
%           if any(~isnan(filt_avg.cluster_weighted_a(it0,:))) && any(~isnan(filt_avg.cluster_weighted_a(it1,:))) && ...
%               any(filt_avg.cluster_weighted_a(it0,:) ~= filt_avg.cluster_weighted_a(it1, :))
%             % fprintf('a filter interpolation #%i flagged: change of clusters between start and end filter events\n', i)
%             filt_interp.fCDOM_a_cluster_chg(it & any(idnonnan)) = true;
%           end
%           % flag when previous filter event was not clustered
%           if filt_avg.a_not_clustered_flag(it0)
%             filt_interp.flag_a_filt0_not_clustered(it & idnonnan) = filt_avg.a_not_clustered_flag(it0);
%           end
%           % flag when next filter event was not clustered
%           if filt_avg.a_not_clustered_flag(it1)
%             filt_interp.flag_a_filt1_not_clustered(it & idnonnan) = filt_avg.a_not_clustered_flag(it1);
%           end
%           % set linear interpolation boolean for missing fdom data
%           filt_interp.flag_linear_interp_a(it & isnan(filt_interp.fdom)) = true;
%         end
%         % flag when mixed c clusters within filter event 0 or filter event 1
%         idnonnan = any(~isnan(tot.c), 2);
%         if any(it & idnonnan)
%           if filt_avg.mix_cluster_c_flag(it0)
%             % fprintf('c filter interpolation #%i flagged: mixed cluster within start filter event\n', i)
%             filt_interp.fCDOM_mix_c_cluster0(it & idnonnan) = true;
%           end
%           if filt_avg.mix_cluster_c_flag(it1)
%             % fprintf('c filter interpolation #%i flagged: mixed cluster within end filter event\n', i)
%             filt_interp.fCDOM_mix_c_cluster1(it & idnonnan) = true;
%           end
%           % flag when c clusters change between filter event 0 and filter event 1
%           if ~all(isnan(filt_avg.cluster_weighted_c(it0,:))) && ~all(isnan(filt_avg.cluster_weighted_c(it1,:))) && ...
%               any(filt_avg.cluster_weighted_c(it0,:) ~= filt_avg.cluster_weighted_c(it1, :))
%             % fprintf('c filter interpolation #%i flagged: change of clusters between start and end filter events\n', i)
%             filt_interp.fCDOM_c_cluster_chg(it & any(idnonnan)) = true;
%           end
%           % flag when previous filter event was not clustered
%           if filt_avg.c_not_clustered_flag(it0)
%             filt_interp.flag_c_filt0_not_clustered(it & idnonnan) = filt_avg.c_not_clustered_flag(it0);
%           end
%           % flag when next filter event was not clustered
%           if filt_avg.c_not_clustered_flag(it1)
%             filt_interp.flag_c_filt1_not_clustered(it & idnonnan) = filt_avg.c_not_clustered_flag(it1);
%           end
%           % set linear interpolation boolean for missing fdom data
%           filt_interp.flag_linear_interp_c(it & isnan(filt_interp.fdom)) = true;
%         end
%       else
%         % linearly interpolate ag and cg if no fdom available
%         filt_interp.flag_linear_interp_a(it) = true;
%         filt_interp.flag_linear_interp_c(it) = true;
%       end
%     end
%   end
% 
%   % set linear interpolation flag true when slope is negative
%   filt_interp.flag_linear_interp_a(filt_interp.flag_a_negative_slope) = true;
%   filt_interp.flag_linear_interp_c(filt_interp.flag_c_negative_slope) = true;
% 
%   % compute ag and cg with fCDOM
%   filt_interp.a(~filt_interp.flag_linear_interp_a, :) = filt_interp.lin_a(~filt_interp.flag_linear_interp_a, :) + ...
%     (filt_interp.fdom(~filt_interp.flag_linear_interp_a) - filt_interp.lin_fdom(~filt_interp.flag_linear_interp_a)) .* ...
%     filt_interp.slope_interp_a(~filt_interp.flag_linear_interp_a, :);
%   filt_interp.c(~filt_interp.flag_linear_interp_c, :) = filt_interp.lin_c(~filt_interp.flag_linear_interp_c, :) + ...
%     (filt_interp.fdom(~filt_interp.flag_linear_interp_c) - filt_interp.lin_fdom(~filt_interp.flag_linear_interp_c)) .* ...
%     filt_interp.slope_interp_c(~filt_interp.flag_linear_interp_c, :);
% 
%   % compute ag and cg with linear interpolation when fCDOM does not work merge filt_avg rows into 
%   % filt_interp to linearly interpolate ag and cg only for missing data (in case fCDOM interpolation 
%   % was available for half of filter event)
%   append_filt_interp = filt_interp(end-size(filt_avg, 1)+1:end, :);
%   var = filt_interp.Properties.VariableNames';
%   var2mv = {'dt','fdom','t','s','a','c'};
%   for i = 1:size(var, 1)
%     if any(strcmp(var{i}, var2mv))
%       append_filt_interp.(var{i}) = filt_avg.(var{i});
%     else
%       if isdatetime(append_filt_interp.(var{i}))
%         append_filt_interp.(var{i}) = NaT(size(append_filt_interp.(var{i})));
%       elseif islogical(append_filt_interp.(var{i}))
%         append_filt_interp.(var{i}) = false(size(append_filt_interp.(var{i})));
%       elseif isnumeric(append_filt_interp.(var{i}))
%         append_filt_interp.(var{i}) = NaN(size(append_filt_interp.(var{i})));
%       end
%     end
%   end
%   % remove 1 second from appended filt_avg.dt to prevent duplicates
%   append_filt_interp.dt = append_filt_interp.dt - seconds(1);
%   % check if any duplicates remaining, if yes, remove 1 second
%   iddup = ismember(append_filt_interp.dt, filt_interp.dt);
%   while any(iddup)
%     append_filt_interp.dt(iddup) = append_filt_interp.dt(iddup) - seconds(1);
%     iddup = ismember(append_filt_interp.dt, filt_interp.dt);
%   end
%   filt_interp = [filt_interp; append_filt_interp];
% 
%   filt_interp = sortrows(filt_interp, 'dt');
%   % linear interpolation when fdom interpolation was not used
%   filt_interp.a = fillmissing(filt_interp.a, 'linear', 'SamplePoints', filt_interp.dt, 'EndValues', 'nearest');
%   filt_interp.c = fillmissing(filt_interp.c, 'linear', 'SamplePoints', filt_interp.dt, 'EndValues', 'nearest');
%   % remove filt_avg rows added for filling missing linear interpolation
%   filt_interp(ismember(filt_interp.dt, append_filt_interp.dt), :) = [];
% 
%   % visProd3D(lambda.a, filt_interp.dt, filt_interp.a, false, 'Wavelength', false, 2);
%   % title('ag fdom interpolated')
%   % visProd3D(lambda.c, filt_interp.dt, filt_interp.c, false, 'Wavelength', false, 3);
%   % title('cg fdom interpolated')
%   % 
%   % visProd3D(lambda.a, filt_interp.dt, tot.a - filt_interp.a, false, 'Wavelength', false, 2);
%   % title('ap fdom interp')
%   % visProd3D(lambda.c, filt_interp.dt, tot.c - filt_interp.c, false, 'Wavelength', false, 3);
%   % title('cp fdom interp')
%   % 
%   % visProd3D(lambda.a, filt_interp.dt, tot.a - filt_interp.lin_a, false, 'Wavelength', false, 4);
%   % title('ap linear interp')
%   % visProd3D(lambda.c, filt_interp.dt, tot.c - filt_interp.lin_c, false, 'Wavelength', false, 5);
%   % title('cp linear interp')
%   % 
%   % visProd3D(lambda.a, filt_interp.dt, (tot.a - filt_interp.lin_a) - (tot.a - filt_interp.a), false, 'Wavelength', false, 4);
%   % title('ap diff')
%   % visProd3D(lambda.c, filt_interp.dt, (tot.c - filt_interp.lin_c) - (tot.c - filt_interp.c), false, 'Wavelength', false, 5);
%   % title('cp diff')
% 
%   % % Fill missing values (start/end of analyzed period) with nearest filt_interp data
%   % filt_interp.a = fillmissing(filt_interp.a, 'nearest', 'SamplePoints', filt_interp.dt);
%   % filt_interp.c = fillmissing(filt_interp.c, 'nearest', 'SamplePoints', filt_interp.dt);
% 
%   % visProd3D(lambda.a, filt_interp.dt, filt_interp.a, false, 'Wavelength', false, 4);
%   % title('ag fdom interpolated and nearest')
%   % visProd3D(lambda.c, filt_interp.dt, filt_interp.c, false, 'Wavelength', false, 5);
%   % title('cg fdom interpolated and nearest')
% 
% 
% 
%   % figure; hold on
%   % scatter(filt_interp.dt, ag_fact_tot(:,70), 20, 'ko')
%   % scatter(filt_interp.dt(~filt_interp.flag_linear_interp), filt_interp.delta_af(~filt_interp.flag_linear_interp,70), 20, 'r+')
%   % 
%   % figure; hold on
%   % scatter(filt_interp.dt, cg_fact_tot(:,70), 20, 'ko')
%   % scatter(filt_interp.dt(~filt_interp.flag_linear_interp), filt_interp.delta_cf(~filt_interp.flag_linear_interp,70), 20, 'r+')
% 
%   % visProd3D(lambda.a, filt_interp.dt(it), filt_interp.a(it,:), false, 'Wavelength', false, i+1);
%   % title('ag fdom interpolated')
%   % visProd3D(lambda.a, filt_interp.dt(it), filt_interp.lin_a(it,:), false, 'Wavelength', false, i+2);
%   % title('ag linearly interpolated')
%   % 
%   % visProd3D(lambda.c, filt_interp.dt(it), filt_interp.c(it,:), false, 'Wavelength', false, i+3);
%   % title('cg fdom interpolated')
%   % visProd3D(lambda.c, filt_interp.dt(it), filt_interp.lin_c(it,:), false, 'Wavelength', false, i+4);
%   % title('cg linearly interpolated')
%   % 
%   % visProd3D(lambda.a, filt_interp.dt(it), filt_interp.a(it,:) - filt_interp.lin_a(it,:), false, 'Wavelength', false, i+5);
%   % title('ag difference')
%   % visProd3D(lambda.a, filt_interp.dt(it), filt_interp.c(it,:) - filt_interp.lin_c(it,:), false, 'Wavelength', false, i+6);
%   % title('cg difference')
%   %
%   % visProd3D(lambda.a, filt_interp.dt, filt_interp.a, false, 'Wavelength', false, i+1);
%   % title('ag fdom interpolated')
%   % visProd3D(lambda.c, filt_interp.dt, filt_interp.c, false, 'Wavelength', false, i+2);
%   % title('cg fdom interpolated')
% 
%   % yaxr = [it_filt_interp.fdom; filt_avg.fdom(it0:it1)];
%   % 
%   % wl_toplot = 40;
%   % figure(1)
%   % clf
%   % subplot(1,2,1)
%   % yyaxis('left')
%   % hold on
%   % scatter(it_filt_interp.dt, filt_interp.a(it,wl_toplot), 40, 'filled', 'MarkerFaceColor', 'b','MarkerFaceAlpha', 0.5)
%   % scatter(filt_avg.dt, filt_avg.a(:,wl_toplot), 100, 'filled', 'MarkerFaceColor', [255	193	37]/255,'MarkerFaceAlpha', 1)
%   % xlim([it_filt_interp.dt(1)-hours(0.5) it_filt_interp.dt(end)+hours(0.5)])
%   % ylabel('ag')
%   % yyaxis('right')
%   % hold on
%   % scatter(filt_interp.dt, filt_interp.fdom, 40, 'filled', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5, 'Marker', 'd')
%   % scatter(filt_avg.dt, filt_avg.fdom, 100, 'filled', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5, 'Marker', 'd')
%   % xlim([it_filt_interp.dt(1)-hours(0.5) it_filt_interp.dt(end)+hours(0.5)])
%   % ylim([min(yaxr)-0.001 max(yaxr)+0.001])
%   % ylabel('fdom')
%   % legend('ag fdom interpolated', 'ag filter average', 'fdom smoothed', 'Location', 'Best')
%   % subplot(1,2,2)
%   % yyaxis('left')
%   % hold on
%   % scatter(it_filt_interp.dt, filt_interp.c(it,wl_toplot), 40, 'filled', 'MarkerFaceColor', 'r','MarkerFaceAlpha', 0.5)
%   % scatter(filt_avg.dt, filt_avg.c(:,wl_toplot), 100, 'filled', 'MarkerFaceColor', [255	193	37]/255,'MarkerFaceAlpha', 1)
%   % xlim([it_filt_interp.dt(1)-hours(0.5) it_filt_interp.dt(end)+hours(0.5)])
%   % ylabel('cg')
%   % yyaxis('right')
%   % hold on
%   % scatter(filt_interp.dt, filt_interp.fdom, 40, 'filled', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5, 'Marker', 'd')
%   % scatter(filt_avg.dt, filt_avg.fdom, 100, 'filled', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5, 'Marker', 'd')
%   % xlim([it_filt_interp.dt(1)-hours(0.5) it_filt_interp.dt(end)+hours(0.5)])
%   % ylim([min(yaxr)-0.001 max(yaxr)+0.001])
%   % ylabel('fdom')
%   % legend('cg fdom interpolated', 'cg filter average', 'fdom smoothed', 'Location', 'Best')
% 
%   % visProd3D(lambda.a, filt_interp.dt, filt_interp.a, false, 'Wavelength', false, 23);
%   % title('ag fdom interpolated')
%   % visProd3D(lambda.a, filt_interp.dt, filt_interp.lin_a, false, 'Wavelength', false, 24);
%   % title('ag linearly interpolated')
%   % 
%   % visProd3D(lambda.c, filt_interp.dt, filt_interp.c, false, 'Wavelength', false, 25);
%   % title('cg fdom interpolated')
%   % visProd3D(lambda.c, filt_interp.dt, filt_interp.lin_c, false, 'Wavelength', false, 26);
%   % title('cg linearly interpolated')
% 
% end
% 
% %% 
% function [filt_avg, filt_interp, filt, regress_stats] = fdom_agcg_model(filt, filt_interp, ...
%   lambda, filt_avg, min_nb_pts_per_cluster, time_weight_for_cluster, cluster_bool)
% 
% %   % regress diff(a) & diff(c) with diff(fdom)
% % 
% %   regress_stats = struct();
% %   [regress_stats.a, regress_stats.c] = regress_acfilt(diff(filt.a), diff(filt.c), diff(filt.fdom));
% % 
% % figure(3)
% % plot_linreg(diff(filt.a(:, 1)), diff(filt.fdom), 'robust', 'linear', false, false)
% % plot_linreg(diff(filt.c(:, 1)), diff(filt.fdom), 'robust', 'linear', false, false)
% % 
% % figure(1); hold on
% % plot((405:4:725)', regress_stats.a.slope)
% % plot((405:4:725)', regress_stats.c.slope)
% % 
% % figure(2); hold on
% % plot((405:4:725)', regress_stats.a.intercept)
% % plot((405:4:725)', regress_stats.c.intercept)
% 
%   % cluster method
%   if cluster_bool
%     % define reference wavelength to do regression (450nm)
%     id450_a = abs(lambda.a - 450) == min(abs(lambda.a - 450));
%     id450_c = abs(lambda.c - 450) == min(abs(lambda.c - 450));
% 
%     % iterate dbscan with multiple epsilon until nRMSE of all cluster < 10%
%     fprintf('Time variable weight %.3f: Finding best epsilon for dbscan Kernel Density clustering of a(filt)/fCDOM(filt) ... ', time_weight_for_cluster)
%     k = 1;
%     epsilons = 1:-0.01:0.01;
%     keep_going_a = true;
%     keep_going_c = true;
%     max_nb_clusters = round(days(max(filt_interp.dt) - min(filt_interp.dt))*2);
%     eval_cluster = table();
%     eval_cluster.nRMSE_a450 = NaN(length(epsilons), max_nb_clusters);
%     eval_cluster.nRMSE_c450 = NaN(length(epsilons), max_nb_clusters);
%     eval_cluster.slope_a450 = NaN(length(epsilons), max_nb_clusters);
%     eval_cluster.slope_c450 = NaN(length(epsilons), max_nb_clusters);
%     eval_cluster.a_perc_not_clustered = NaN(length(epsilons), 1);
%     eval_cluster.c_perc_not_clustered = NaN(length(epsilons), 1);
%     eval_cluster.a_perc_negative_slope = NaN(length(epsilons), 1);
%     eval_cluster.c_perc_negative_slope = NaN(length(epsilons), 1);
%     % normalize the time variable between 0 to 0.5
%     time_var = (filt.dt-min(filt.dt))/max(filt.dt-min(filt.dt)) * time_weight_for_cluster;
%     data_tocluster_a = [filt.a(:, id450_a)./filt.fdom time_var]; % time_var filt.t filt.s
%     data_tocluster_c = [filt.c(:, id450_c)./filt.fdom time_var]; % time_var filt.t filt.s
% 
%     if any(~isnan(data_tocluster_a(:, 1))) || any(~isnan(data_tocluster_c(:, 1)))
%       % find minimum cluster size depending on the number of days loaded
%       % round(size(time_var, 1) / (max(day(filt.dt-min(filt.dt)+1)) * 24))
%       while k <= length(epsilons) && keep_going_a && keep_going_c
%         if keep_going_a
%           filt.a_clusters = dbscan(data_tocluster_a, epsilons(k), min_nb_pts_per_cluster, 'Distance', 'euclidean');
%           filt.a_clusters(filt.a_clusters == -1) = NaN;
%         end
%         if keep_going_c
%           filt.c_clusters = dbscan(data_tocluster_c, epsilons(k), min_nb_pts_per_cluster, 'Distance', 'euclidean');
%           filt.c_clusters(filt.c_clusters == -1) = NaN;
%         end
%         % regress a&c with fdom
%         [regress_stats_a, regress_stats_c] = regress_acfilt(filt.a(:, id450_a), filt.c(:, id450_c), filt.fdom, filt.a_clusters, filt.c_clusters);
%         if keep_going_a
%           eval_cluster.nRMSE_a450(k, 1:size(regress_stats_a.nRMSE, 2)) = regress_stats_a.nRMSE;
%           eval_cluster.slope_a450(k, 1:size(regress_stats_a.slope, 2)) = regress_stats_a.slope;
%           eval_cluster.a_perc_not_clustered(k) = sum(isnan(filt.a_clusters))/size(filt,1)*100;
%           if any(regress_stats_a.slope <= 0)
%             eval_cluster.a_perc_negative_slope(k) = sum(sum(filt.a_clusters==find(regress_stats_a.slope<=0))) / size(filt,1)*100;
%           else
%             eval_cluster.a_perc_negative_slope(k) = 0;
%           end
%         end
%         if keep_going_c
%           eval_cluster.slope_c450(k, 1:size(regress_stats_c.slope, 2)) = regress_stats_c.slope;
%           eval_cluster.nRMSE_c450(k, 1:size(regress_stats_c.nRMSE, 2)) = regress_stats_c.nRMSE;
%           eval_cluster.c_perc_not_clustered(k) = sum(isnan(filt.c_clusters))/size(filt,1)*100;
%           if any(regress_stats_c.slope <= 0)
%             eval_cluster.c_perc_negative_slope(k) = sum(sum(filt.c_clusters==find(regress_stats_c.slope<=0))) / size(filt,1)*100;
%           else
%             eval_cluster.c_perc_negative_slope(k) = 0;
%           end
%         end
%         if k >= 2
%           % get epsilon_a
%           if all(eval_cluster.slope_a450(k,1:size(regress_stats_a.slope, 2)) > 0) && ...
%               all(eval_cluster.nRMSE_a450(k,1:size(regress_stats_a.nRMSE, 2)) < 25, 2)
%             % epsilon_a = epsilons(k);
%             keep_going_a = false;
%           elseif size(regress_stats_a.slope, 2) > max_nb_clusters
%             keep_going_a = false;
%           end
%           % get epsilon_c
%           if all(eval_cluster.slope_c450(k,1:size(regress_stats_c.slope, 2)) > 0) && ...
%               all(eval_cluster.nRMSE_c450(k,1:size(regress_stats_c.nRMSE, 2)) < 25, 2)
%             % epsilon_c = epsilons(k);
%             keep_going_c = false;
%           elseif size(regress_stats_c.slope, 2) > max_nb_clusters
%             keep_going_c = false;
%           end
%         end
%         k = k + 1;
%       end
%       % find best epsilon_a based on minimum % data clustered with slope < 0 and minimum % data not clustered
%       [~, ord_a] = sortrows([eval_cluster.a_perc_negative_slope eval_cluster.a_perc_not_clustered], [1 2]);
%       epsilon_a = epsilons(ord_a(1));
%       % find best epsilon_c based on minimum % data clustered with slope < 0 and minimum % data not clustered
%       [~, ord_c] = sortrows([eval_cluster.c_perc_negative_slope eval_cluster.c_perc_not_clustered], [1 2]);
%       epsilon_c = epsilons(ord_c(1));
%       fprintf('done\n')
% 
%       % run clustering with the best epsilon
%       fprintf('Kernel Density clustering of a(filt)/fCDOM(filt) ... ')
%       % cluster absorption
%       filt.a_not_clustered_flag = false(size(filt, 1), 1);
%       filt.a_clusters = dbscan(data_tocluster_a, epsilon_a, min_nb_pts_per_cluster, 'Distance', 'euclidean');
%       filt.a_not_clustered_flag(filt.a_clusters == -1) = true;
%       if ~all(filt.a_clusters == -1)
%         filt.a_clusters(filt.a_clusters == -1) = NaN;
%         filt.a_clusters = fillmissing(filt.a_clusters, 'nearest', 'SamplePoints', filt.dt);
%       end
%       % cluster attenuation
%       filt.c_not_clustered_flag = false(size(filt, 1), 1);
%       filt.c_clusters = dbscan(data_tocluster_c, epsilon_c, min_nb_pts_per_cluster, 'Distance', 'euclidean');
%       filt.c_not_clustered_flag(filt.c_clusters == -1) = true;
%       if ~all(filt.a_clusters == -1)
%         filt.c_clusters(filt.c_clusters == -1) = NaN;
%         filt.c_clusters = fillmissing(filt.c_clusters, 'nearest', 'SamplePoints', filt.dt);
%       end
%       fprintf('done\n')
%     else
%       % regress_stats = [];
%       if all(isnan(filt.fdom))
%         error('All fDOM filter event data loaded are NaNs')
%       end
%       if all(isnan(data_tocluster_a))
%         error('All "a" filter event data loaded are NaNs')
%       end
%       if all(isnan(data_tocluster_c))
%         error('All "c" filter event data loaded are NaNs')
%       end
%     end
%   else
%     filt.a_clusters = ones(size(filt.dt));
%     filt.a_not_clustered_flag = false(size(filt, 1), 1);
%     filt.c_clusters = ones(size(filt.dt));
%     filt.c_not_clustered_flag = false(size(filt, 1), 1);
%   end
%   % regress a&c with fdom
%   regress_stats = struct();
%   [regress_stats.a, regress_stats.c] = regress_acfilt(filt.a, filt.c, filt.fdom, filt.a_clusters, filt.c_clusters);
% 
%   % % flag data part of cluster with negative slope
%   % clusters.a_negative_slope = false(size(filt, 1), 1);
%   % clusters.c_negative_slope = false(size(filt, 1), 1);
%   % if any(regress_stats_a.slope <= 0)
%   %   clusters.a_negative_slope = filt.a_clusters == find(regress_stats_a.slope<=0);
%   % end
%   % if any(regress_stats_c.slope <= 0)
%   %   clusters.c_negative_slope = filt.c_clusters == find(regress_stats_c.slope<=0);
%   % end
% 
%   % % find wavelength where fdom fit and temperature fit work poorly
%   % bad_fit_a = all(regress_stats.a.R2 < prctile(regress_stats.a.R2, 5, 1), 2);
%   % bad_fit_c = all(regress_stats.c.R2 < prctile(regress_stats.c.R2, 5, 1), 2);
% 
%   %%% apply the regression to get a and c dissolved from fdom
%   groups_a = unique(filt.a_clusters(~isnan(filt.a_clusters)));
%   groups_c = unique(filt.c_clusters(~isnan(filt.c_clusters)));
% 
%   % Weight cluster for each filter event in case clusters are changing
%   % during a filter event (unlikely with dbscan method but let's keep it in case)
%   filt_avg = addvars(filt_avg, NaN(size(filt_avg,1),size(groups_a,1)), NaN(size(filt_avg,1),size(groups_c,1)), ...
%     false(size(filt_avg.dt)), false(size(filt_avg.dt)),...
%     false(size(filt_avg.dt)), false(size(filt_avg.dt)),...
%     false(size(filt_avg.dt)), false(size(filt_avg.dt)),...
%     'NewVariableNames', {'cluster_weighted_a', 'cluster_weighted_c',...
%     'mix_cluster_a_flag', 'mix_cluster_c_flag',...
%     'cluster_a_slope_flag', 'cluster_c_slope_flag',...
%     'a_not_clustered_flag', 'c_not_clustered_flag'}, 'after', 'end');
%   for i=1:size(filt_avg, 1)
%     sel_filt = filt_avg.start(i) <= filt.dt & filt.dt <= filt_avg.end(i);
%     foo = filt(sel_filt,:);
%     % flag if clusters are not all the same within a single filter event
%     if size(unique(foo.a_clusters(~isnan(foo.a_clusters))), 1) > 1
%       filt_avg.mix_cluster_a_flag(i) = true;
%     end
%     if size(unique(foo.c_clusters(~isnan(foo.c_clusters))), 1) > 1
%       filt_avg.mix_cluster_c_flag(i) = true;
%     end
%     % flag if any a_filt data was not clustered in dbscan
%     if any(foo.a_not_clustered_flag)
%       filt_avg.a_not_clustered_flag(i) = true;
%     end
%     % flag if any a_filt data was not clustered in dbscan
%     if any(foo.c_not_clustered_flag)
%       filt_avg.c_not_clustered_flag(i) = true;
%     end
%     % weight clusters proportionally for each filter events
%     for j = 1:size(groups_a,1)
%       filt_avg.cluster_weighted_a(i,j) = sum(foo.a_clusters == groups_a(j)) ./ sum(~isnan(foo.a_clusters));
%     end
%     for j = 1:size(groups_c,1)
%       filt_avg.cluster_weighted_c(i,j) = sum(foo.c_clusters == groups_c(j)) ./ sum(~isnan(foo.c_clusters));
%     end
%   end
%   % weight a slope and intercept
%   filt_avg.slope_a = zeros(size(filt_avg, 1), size(regress_stats.a, 1));
%   filt_avg.intercept_a = zeros(size(filt_avg, 1), size(regress_stats.a, 1));
%   % Define same exponential function as in FitExp to get the slopes computed from the base and intercept parameters
%   expfun = @(p, xd) p(1) * exp(p(2) * (xd - 440));
%   for j = 1:size(groups_a, 1)
%     % Fit exponential function to slopes vs lambda in the blue to reconstruct the slope values in the red where S/N is very low
%     [y_intercp_a_slopefit, base_a_slopefit, ~, ~, ~] = FitExp(lambda.a(lambda.a < 500), ...
%       regress_stats.a.slope(lambda.a < 500, j)', regress_stats.a.RMSE(lambda.a < 500, j)');
%     fitted_a_slope = expfun([y_intercp_a_slopefit base_a_slopefit], lambda.a);
% 
%     % using slopes directly
%     filt_avg.slope_a = filt_avg.slope_a + fitted_a_slope.*filt_avg.cluster_weighted_a(:,j);
%     filt_avg.intercept_a = filt_avg.intercept_a + regress_stats.a.intercept(:,j)'.*filt_avg.cluster_weighted_a(:,j);
% 
%     % figure; 
%     % subplot(2,1,1); plot(lambda.a, regress_stats.a.slope)
%     % subplot(2,1,2); plot(lambda.a, regress_stats.a.intercept)
% 
%     % figure(3)
%     % subplot(1, 2, 1)
%     % plot(lambda.a, regress_stats.a.slope, 'bo');  hold on;   %plot your raw data
%     % plot(lambda.a, fitted_a_slope, 'r-');  %plot the fit data
%     % ylabel('a slope of a_{filt}-fdom correlation')
%     % xlabel('\lambda')
% 
%   end
%   % weight c slope and intercept
%   filt_avg.slope_c = zeros(size(filt_avg, 1), size(regress_stats.c, 1));
%   filt_avg.intercept_c = zeros(size(filt_avg, 1), size(regress_stats.c, 1));
%   for j = 1:size(groups_c, 1)
%     % Fit exponential function to slopes vs lambda in the blue to reconstruct the slope values in the red where S/N is very low
%     [y_intercp_c_slopefit, base_c_slopefit, ~, ~, ~] = FitExp(lambda.c(lambda.c < 500), ...
%       regress_stats.c.slope(lambda.c < 500, j)', regress_stats.c.RMSE(lambda.c < 500, j)');
%     fitted_c_slope = expfun([y_intercp_c_slopefit base_c_slopefit], lambda.c);
% 
%     % using slopes directly
%     filt_avg.slope_c = filt_avg.slope_c + fitted_c_slope.*filt_avg.cluster_weighted_c(:,j);
%     filt_avg.intercept_c = filt_avg.intercept_c + regress_stats.c.intercept(:,j)'.*filt_avg.cluster_weighted_c(:,j);
% 
%     % figure; 
%     % subplot(4,2,2); plot(lambda.c, regress_stats.c.slope)
%     % subplot(4,2,4); plot(lambda.c, regress_stats.c.intercept)
% 
%     % subplot(1, 2, 2)
%     % plot(lambda.c, regress_stats.c.slope, 'bo');  hold on;   %plot your raw data
%     % plot(lambda.c, fitted_c_slope, 'r-');  %plot the fit data
%     % ylabel('c slope of c_{filt}-fdom correlation')
%     % xlabel('\lambda')
% 
%   end
%   % interpolate slope and intercept cluster specific onto filt_interp.dt
%   filt_interp.slope_interp_a = interp_extrap(filt_avg, filt_interp.dt, 'slope_a', [], true, 'linear', 'nearest');
%   filt_interp.intercept_interp_a = interp_extrap(filt_avg, filt_interp.dt, 'intercept_a', [],  true, 'linear', 'nearest');
%   filt_interp.slope_interp_c = interp_extrap(filt_avg, filt_interp.dt, 'slope_c', [],  true, 'linear', 'nearest');
%   filt_interp.intercept_interp_c = interp_extrap(filt_avg, filt_interp.dt, 'intercept_c', [],  true, 'linear', 'nearest');
%   % flag data for which a/fdom and/or c/fdom correlation has a negative slope
%   filt_interp.flag_a_negative_slope = false(size(filt_interp.dt));
%   filt_interp.flag_c_negative_slope = false(size(filt_interp.dt));
%   filt_interp.flag_a_negative_slope(any(filt_interp.slope_interp_a < 0, 2)) = true;
%   filt_interp.flag_c_negative_slope(any(filt_interp.slope_interp_c < 0, 2)) = true;
% end
% 
% 
% %% regression between a/c and other variable in filter events (including clusters)
% function [a_reg, c_reg] = regress_acfilt(a, c, ancillary, a_clusters, c_clusters)
%   if nargin < 4
%     a_clusters = ones(size(a, 1), 1);
%     c_clusters = ones(size(c, 1), 1);
%   end
%   % robust linear regression between filt a&c and ancillary variable
%   a_reg = array2table(NaN(size(a, 2), 5), 'VariableNames', {'slope', 'intercept', 'RMSE', 'nRMSE', 'R2'});
%   c_reg = array2table(NaN(size(c, 2), 5), 'VariableNames', {'slope', 'intercept', 'RMSE', 'nRMSE', 'R2'});
%   groups_a = unique(a_clusters(~isnan(a_clusters)));
%   groups_c = unique(c_clusters(~isnan(c_clusters)));
%   for i = 1:size(a, 2)
%     for j = 1:size(groups_a, 1)
%       id_grp_a = a_clusters == groups_a(j);
%       % regress variable with filtered water absorption
%       [stats.b, stats.stats] = robustfit(ancillary(id_grp_a), a(id_grp_a, i));
%       [~, MSGID] = lastwarn();
%       warning('off', MSGID)
%       a_reg.slope(i, j) = stats.b(2);
%       a_reg.intercept(i, j) = stats.b(1);
%       a_reg.RMSE(i, j) = stats.stats.robust_s;
%       % calculated normalize RMSE by average and R2
%       id_nonan = all(~isnan([a(:, i) ancillary]), 2);
%       a_reg.nRMSE(i, j) = a_reg.RMSE(i, j)/abs(mean(a(id_grp_a & id_nonan, i), 'omitnan'))*100;
%       a_reg.R2(i, j) = corr(a(id_grp_a & id_nonan, i), stats.b(1) + stats.b(2) * ancillary(id_grp_a & id_nonan))^2;
%     end
%     for j = 1:size(groups_c, 1)
%       id_grp_c = c_clusters == groups_c(j);
%       % regress variable with filtered water attenuation
%       [stats.b, stats.stats] = robustfit(ancillary(id_grp_c), c(id_grp_c, i));
%       [~, MSGID] = lastwarn();
%       warning('off', MSGID)
%       c_reg.slope(i, j) = stats.b(2);
%       c_reg.intercept(i, j) = stats.b(1);
%       c_reg.RMSE(i, j) = stats.stats.robust_s;
%       % calculated normalize RMSE by average and R2
%       id_nonan = all(~isnan([c(:, i) ancillary]), 2);
%       c_reg.nRMSE(i, j) = c_reg.RMSE(i, j)/abs(mean(c(id_grp_c & id_nonan, i), 'omitnan'))*100;
%       c_reg.R2(i, j) = corr(c(id_grp_c & id_nonan, i), stats.b(1) + stats.b(2) * ancillary(id_grp_c & id_nonan))^2;
%     end
%   end
% end
% 
% 
% % %% regression between a/c and other variable in filter events
% % function [a_reg, c_reg] = regress_acfilt(a, c, ancillary)
% %   % robust linear regression between filt a&c and ancillary variable
% %   a_reg = array2table(NaN(size(a, 2), 5), 'VariableNames', {'slope', 'intercept', 'RMSE', 'nRMSE', 'R2'});
% %   c_reg = array2table(NaN(size(c, 2), 5), 'VariableNames', {'slope', 'intercept', 'RMSE', 'nRMSE', 'R2'});
% %   for i = 1:size(a, 2)
% %     % regress variable with filtered water absorption
% %     [stats.b, stats.stats] = robustfit(ancillary, a(:, i));
% %     [~, MSGID] = lastwarn();
% %     warning('off', MSGID)
% %     a_reg.slope(i, :) = stats.b(2);
% %     a_reg.intercept(i, :) = stats.b(1);
% %     a_reg.RMSE(i, :) = stats.stats.robust_s;
% %     % calculated normalize RMSE by average and R2
% %     id_nonan = all(~isnan([a(:, i) ancillary]), 2);
% %     a_reg.nRMSE(i, :) = a_reg.RMSE(i, :)/abs(mean(a(id_nonan, i), 'omitnan'))*100;
% %     a_reg.R2(i, :) = corr(a(id_nonan, i), stats.b(1) + stats.b(2) * ancillary(id_nonan))^2;
% % 
% %     % regress variable with filtered water attenuation
% %     [stats.b, stats.stats] = robustfit(ancillary, c(:, i));
% %     [~, MSGID] = lastwarn();
% %     warning('off', MSGID)
% %     c_reg.slope(i, :) = stats.b(2);
% %     c_reg.intercept(i, :) = stats.b(1);
% %     c_reg.RMSE(i, :) = stats.stats.robust_s;
% %     % calculated normalize RMSE by average and R2
% %     id_nonan = all(~isnan([c(:, i) ancillary]), 2);
% %     c_reg.nRMSE(i, :) = c_reg.RMSE(i, :)/abs(mean(c(id_nonan, i), 'omitnan'))*100;
% %     c_reg.R2(i, :) = corr(c(id_nonan, i), stats.b(1) + stats.b(2) * ancillary(id_nonan))^2;
% %   end
% % end



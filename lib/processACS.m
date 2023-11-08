function [p, g, bad, regression_stats] = processACS(lambda, tot, filt, param, modelG50, modelmphi, di, ...
  CDOM, fth, fth_constants, interpolation_method, di_method, scattering_correction, compute_ad_aphi, tsg)
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
  if ~any(strcmp(scattering_correction, {'Zaneveld1994', 'Kostakis2022'}))
    error('%s residual temperature and scattering correction not supported', scattering_correction)
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
  % check if fCDOM data loaded if fCDOM interpolation
  if strcmp(interpolation_method, 'CDOM')
    % Require both CDOM & Switch position
    if ~isempty(CDOM.prod.a) || ~isempty(CDOM.qc.tsw)
      if ~isempty(CDOM.prod.a)
        cdom_base = CDOM.prod.a;
      else
        cdom_base = CDOM.qc.tsw;
      end
      if ~any(cdom_base.dt >= min([tot.dt; filt.dt]) & cdom_base.dt <= max([tot.dt; filt.dt]))
        fprintf('Warning: fCDOM dates do not correspond to ACS dates: interpolation switched to "linear"\n')
        interpolation_method = 'linear';
      else
        % remove duplicates
        [~, L, ~] = unique(cdom_base.dt,'first');
        indexToDump = not(ismember(1:numel(cdom_base.dt), L));
        if sum(indexToDump) > 0
          fprintf('Warning: %i identical dates in fCDOM data => deleted\n', sum(indexToDump))
          cdom_base(indexToDump, :) = [];
        end
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
    else
      tsg_data = tsg.qc.tsw;
    end
    if ~any(tsg_data.dt >= min([tot.dt; filt.dt]) & tsg_data.dt <= max([tot.dt; filt.dt]))
      fprintf('Warning: TSG dates do not correspond to ACS dates: salinity correction not performed\n')
      tsg_loaded = false;
    else
      tsg_loaded = true;
      % remove duplicates
      [~, L, ~] = unique(tsg_data.dt,'first');
      indexToDump = not(ismember(1:numel(tsg_data.dt), L));
      if sum(indexToDump) > 0
        fprintf('Warning: %i identical dates in TSG data => deleted\n', sum(indexToDump))
        tsg_data(indexToDump, :) = [];
      end
    end
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

  % create ancillary table and prepare T/S if loaded
  ancillary = table();
  foo_dt = [datetime(min([tot.dt; filt.dt]), 'ConvertFrom', 'datenum') datetime(max([tot.dt; filt.dt]), 'ConvertFrom', 'datenum')];
  ancillary.dt = datenum(dateshift((foo_dt(1):minutes(1):foo_dt(2))','start','minute'));
  if tsg_loaded
    % interpolate T and S on ancillary table
    ancillary.t = interp1(tsg_data.dt, tsg_data.(tsg.temperature_variable), ancillary.dt, 'spline');
    ancillary.t = fillmissing(ancillary.t, 'nearest');
    ancillary.s = interp1(tsg_data.dt, tsg_data.s, ancillary.dt, 'spline');
    ancillary.s = fillmissing(ancillary.s, 'nearest');

    % replace interpolated values over gaps > 60 min by NaN
    missing_tsg= ~ismember(ancillary.dt, tsg_data.dt);
    t = [true; diff(missing_tsg) ~= 0];
    k = diff(find([t; true])) .* missing_tsg(t);
    long_nan = k(cumsum(t)) > 120;
    ancillary.t(long_nan) = NaN;
    ancillary.s(long_nan) = NaN;
    % % fill long gaps with nearest values
    % ancillary.t = fillmissing(ancillary.t, 'nearest');
    % ancillary.s = fillmissing(ancillary.s, 'nearest');

    % interpolate spline T on filt_interp table 
    filt_interp = addvars(filt_interp, interp1(ancillary.dt, ancillary.t, filt_interp.dt, 'spline'), ...
      'NewVariableNames', 't', 'After', 'dt');
    % interpolate spline S on filt_interp table
    filt_interp = addvars(filt_interp, interp1(ancillary.dt, ancillary.s, filt_interp.dt, 'spline'), ...
      'NewVariableNames', 's', 'After', 'dt');

    % interpolate T and S on filtered data 
    filt = addvars(filt, interp1(ancillary.dt, ancillary.t, filt.dt, 'spline'), 'NewVariableNames', 't', 'After', 'dt');
    filt = addvars(filt, interp1(tsg_data.dt, tsg_data.s, filt.dt, 'spline'), 'NewVariableNames', 's', 'After', 'dt');
    
    % add T and S variables on filt_avg table
    filt_avg = addvars(filt_avg, NaN(size(filt_avg,1),1), 'NewVariableNames', 's', 'After', 'dt');
    filt_avg = addvars(filt_avg, NaN(size(filt_avg,1),1), 'NewVariableNames', 't', 'After', 'dt');

    % remove psiS and psiT with Tref = median(T) of whole processed section
    Tref = median(ancillary.t, 'omitnan');
    a_psiS = interp1(psi.wl, psi.a_psiS, lambda.a, 'spline');
    c_psiS = interp1(psi.wl, psi.c_psiS, lambda.c, 'spline');
    a_psiT = interp1(psi.wl, psi.psiT, lambda.a, 'spline');
    c_psiT = interp1(psi.wl, psi.psiT, lambda.c, 'spline');
    filt.a = filt.a - (a_psiS .* filt.s + a_psiT .* (filt.t - Tref));
    filt.c = filt.c - (c_psiS .* filt.s + c_psiT .* (filt.t - Tref));
    tot.a = tot.a - (a_psiS .* filt_interp.s + a_psiT .* (filt_interp.t - Tref));
    tot.c = tot.c - (c_psiS .* filt_interp.s + c_psiT .* (filt_interp.t - Tref));
  end

  % prepare fCDOM data if CDOM interpolation
  if strcmp(interpolation_method, 'CDOM')
    % create new fdom table interpolating for gaps < 60 min
    ancillary.fdom = interp1(cdom_base.dt, cdom_base.fdom, ancillary.dt, 'spline');
    ancillary.fdom = fillmissing(ancillary.fdom, 'nearest');
    % smooth cdom data removing frequencies higher than 10min
    d = designfilt('lowpassiir','FilterOrder',1, 'PassbandFrequency',1/(60*10),'SampleRate', 1/60);
    ancillary.fdom_10smooth = filtfilt(d, ancillary.fdom);

    % replace interpolated values over gaps > 60 min by NaN
    missing_cdom = ~ismember(ancillary.dt, cdom_base.dt);
    t = [true; diff(missing_cdom) ~= 0];
    k = diff(find([t; true])) .* missing_cdom(t);
    long_nan = k(cumsum(t)) > 120;
    ancillary.fdom(long_nan) = NaN;
    ancillary.fdom_10smooth(long_nan) = NaN;
    % fill long gaps with nearest values to force agcg_fdom_interpolation function to do linear interpolation
    ancillary.fdom = fillmissing(ancillary.fdom, 'nearest');
    ancillary.fdom_10smooth = fillmissing(ancillary.fdom_10smooth, 'nearest');

    % interpolate spline and extrapolate nearest fdom on filt_interp table
    filt_interp = addvars(filt_interp, interp1(ancillary.dt, ancillary.fdom, filt_interp.dt, 'spline'), ...
      'NewVariableNames', 'fdom', 'After', 'dt');
    filt_interp.fdom = fillmissing(filt_interp.fdom, 'nearest');

    % interpolate fdom on filtered data
    filt = addvars(filt, interp1(ancillary.dt, ancillary.fdom, filt.dt, 'spline'), ...
      'NewVariableNames', 'fdom', 'After', 'dt');
    filt.fdom = fillmissing(filt.fdom, 'nearest');

    % add fdom variables on filt_avg table
    filt_avg = addvars(filt_avg, NaN(size(filt_avg,1),1), 'NewVariableNames', 'fdom_10smooth', 'After', 'dt');
    filt_avg = addvars(filt_avg, NaN(size(filt_avg,1),1), 'NewVariableNames', 'fdom', 'After', 'dt');
  end

  % figure; hold on
  % yyaxis('left')
  % scatter(datetime(filt.dt, 'ConvertFrom', 'datenum'), filt.t, 20, 'ro', 'filled')
  % scatter(datetime(filt_interp.dt, 'ConvertFrom', 'datenum'), filt_interp.t, 20, 'b+')
  % ylabel('Temperature')
  % yyaxis('right')
  % scatter(datetime(filt.dt, 'ConvertFrom', 'datenum'), filt.fdom, 20, 'co', 'filled')
  % scatter(datetime(filt_interp.dt, 'ConvertFrom', 'datenum'), filt_interp.fdom, 20, 'co', 'filled')
  % ylabel('fCDOM')
  
  for i=1:size(sel_start, 1)
    sel_filt = fth_interp.dt(sel_start(i)) <= filt.dt & filt.dt <= fth_interp.dt(sel_end(i));
    foo = filt(sel_filt,:);
    if strcmp(interpolation_method, 'CDOM') || tsg_loaded
      sel_anc = fth_interp.dt(sel_start(i)) <= ancillary.dt & ancillary.dt <= fth_interp.dt(sel_end(i));
      foo_anc = ancillary(sel_anc,:);
    end
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
      if tsg_loaded
        filt_avg.t(i,:) = foo_anc.t;
        filt_avg.s(i,:) = foo_anc.s;
      end
      if strcmp(interpolation_method, 'CDOM')
        filt_avg.fdom(i,:) = foo_anc.fdom;
        filt_avg.fdom_10smooth(i,:) = foo_anc.fdom_10smooth;
      end
    elseif sum(sel_filt) > 1
      a_perc25 = foo.a > prctile(foo.a, 25, 1);
      c_perc25 = foo.c > prctile(foo.c, 25, 1);
      foo.a_avg_sd(a_perc25) = NaN;
      foo.c_avg_sd(c_perc25) = NaN;
      foo.a(a_perc25) = NaN;
      foo.c(c_perc25) = NaN;
      % compute average of all values smaller than 25th percentile for each filter event
      filt_avg.a(i,:) = mean(foo.a, 1, 'omitnan');
      filt_avg.c(i,:) = mean(foo.c, 1, 'omitnan');
      filt_avg.dt(i) = mean(foo.dt(any(a_perc25, 2) | any(c_perc25, 2)), 'omitnan');
      filt_avg.start(i) = min(foo.dt(any(a_perc25, 2) | any(c_perc25, 2)));
      filt_avg.end(i) = max(foo.dt(any(a_perc25, 2) | any(c_perc25, 2)));
      filt_avg.a_avg_sd(i,:) = mean(foo.a_avg_sd, 1, 'omitnan');
      filt_avg.c_avg_sd(i,:) = mean(foo.c_avg_sd, 1, 'omitnan');
      filt_avg.a_avg_n(i) = sum(foo.a_avg_n(any(~isnan(foo.a), 2)), 'omitnan');
      filt_avg.c_avg_n(i) = sum(foo.c_avg_n(any(~isnan(foo.c), 2)), 'omitnan');
      if tsg_loaded
        foo_anc.t(foo_anc.t > prctile(foo_anc.t, 25, 1)) = NaN;
        foo_anc.s(foo_anc.s > prctile(foo_anc.s, 25, 1)) = NaN;
        % compute average of all values smaller than 25th percentile for each filter event
        filt_avg.t(i,:) = mean(foo_anc.t, 1, 'omitnan');
        filt_avg.s(i,:) = mean(foo_anc.s, 1, 'omitnan');
      end
      if strcmp(interpolation_method, 'CDOM')
        foo_anc.fdom(foo_anc.fdom > prctile(foo_anc.fdom, 25, 1)) = NaN;
        foo_anc.fdom_10smooth(foo_anc.fdom_10smooth > prctile(foo_anc.fdom_10smooth, 25, 1)) = NaN;
        % compute average of all values smaller than 25th percentile for each filter event
        filt_avg.fdom(i,:) = mean(foo_anc.fdom, 1, 'omitnan');
        filt_avg.fdom_10smooth(i,:) = mean(foo_anc.fdom_10smooth, 1, 'omitnan');
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
      [filt_interp] = agcg_fdom_interpolation(filt_interp, filt, lambda, CDOM.dark, filt_avg);
  
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
      filt_interp.c = interp1(filt_avg.dt, filt_avg.c, filt_interp.dt, 'linear');
      filt_interp.a_avg_sd = interp1(filt_avg.dt, filt_avg.a_avg_sd, filt_interp.dt, 'linear');
      filt_interp.c_avg_sd = interp1(filt_avg.dt, filt_avg.c_avg_sd, filt_interp.dt, 'linear');
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
    fh = visFlag([], filt_interp, tot, [], filt_avg, [], 'a', round(size(tot.a, 2)/2), [], []);
    title('Check filter event interpolation, press q to continue', 'FontSize', 14)
    if strcmp(interpolation_method, 'CDOM')
      ax1 = gca;
      ax1.YColor = [0	205	205]/255;
      scatter(filt_interp.dt, filt_interp.fdom, '.', 'MarkerEdgeColor', [0	205	205]/255, 'MarkerEdgeAlpha', 0.5)
      plot(filt_avg.dt, filt_avg.fdom, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0	205	205]/255)
      ylabel('FDOM (volts)')
      legend('Filtered interpolated', 'Total', 'Filtered percentile average', 'fCDOM', ...
        'fCDOM filtered percentile average', 'AutoUpdate','off', 'FontSize', 12)
    else
      legend('Filtered interpolated', 'Total', 'Filtered median',...
        'AutoUpdate','off', 'FontSize', 12)
    end
    guiSelectOnTimeSeries(fh);
  end
  
% figure(1);
% subplot(1, 3, 1)
% scatter(filt_interp.a(filt_interp.lin_bool,40), filt_interp.s(filt_interp.lin_bool), 20, 'k+')
% scatter(filt_interp.a(~filt_interp.lin_bool,40), filt_interp.s(~filt_interp.lin_bool), 20, 'bo', 'filled')
% ylabel('Salinity')
% xlabel('ag')
% subplot(1, 3, 2)
% scatter(filt_interp.c(filt_interp.lin_bool,40), filt_interp.s(filt_interp.lin_bool), 20, 'k+')
% scatter(filt_interp.c(~filt_interp.lin_bool,40), filt_interp.s(~filt_interp.lin_bool), 20, 'bo', 'filled')
% ylabel('Salinity')
% xlabel('cg')
% subplot(1, 3, 3)
% scatter(filt_interp.fdom(filt_interp.lin_bool), filt_interp.s(filt_interp.lin_bool), 20, 'k+')
% scatter(filt_interp.fdom(~filt_interp.lin_bool), filt_interp.s(~filt_interp.lin_bool), 20, 'bo', 'filled')
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
    fprintf(['ap ' scattering_correction ' residual temperature and scattering correction ... '])
    if strcmp(scattering_correction, 'Kostakis2022')
      [p.ap, ~] = ResidualTemperatureAndScatteringCorrectionKostakis(p.ap, cp_for_apresiduals_corr, lambda.a, psi);
    elseif strcmp(scattering_correction, 'Zaneveld1994')
      [p.ap, ~] = ResidualTemperatureAndScatteringCorrectionZaneveld(p.ap, cp_for_apresiduals_corr, lambda.a, psi);
    end
    fprintf('Done\n')
    fprintf(['cp ' scattering_correction ' residual temperature and scattering correction ... '])
    if strcmp(scattering_correction, 'Kostakis2022')
      [~, p.cp] = ResidualTemperatureAndScatteringCorrectionKostakis(ap_for_cpresiduals_corr, p.cp, lambda.c, psi);
    elseif strcmp(scattering_correction, 'Zaneveld1994')
      [~, p.cp] = ResidualTemperatureAndScatteringCorrectionZaneveld(ap_for_cpresiduals_corr, p.cp, lambda.c, psi);
    end
    fprintf('Done\n')
  else
    fprintf(['ap and cp ' scattering_correction ' residual temperature and scattering correction ... '])
    if strcmp(scattering_correction, 'Kostakis2022')
      [p.ap, p.cp] = ResidualTemperatureAndScatteringCorrectionKostakis(p.ap, p.cp, lambda.ref, psi);
    elseif strcmp(scattering_correction, 'Zaneveld1994')
      [p.ap, p.cp] = ResidualTemperatureAndScatteringCorrectionZaneveld(p.ap, p.cp, lambda.ref, psi);
    end
  end
  
  % Remove lines full of NaNs (Kostakis2022 potentially fail when temperature correct is too large
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
  flag = array2table(false(size(p,1), 21), 'VariableNames', {'cp_neg','ap_neg',...
    'ap_shape','ap_bubbles','ap430_700_neg','cp_over10','noisy600_650',...
    'ap460_640_04_450','positive_ap450_570','poc_flag','chl_ap676lh_flag',...
    'gamma_flag','chl_Halh_flag','HH_mphi_flag','HH_G50_flag','chlratio_flag',...
    'gamma_suspicious','poc_suspicious','chl_ap676lh_suspicious',...
    'chl_Halh_suspicious','HH_G50_mphi_suspicious'});
  
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
    % flag potential bubbles using the step in the center of the spectra
    foo_difstep = diff(p.ap(:, lambda.a > 560 & lambda.a <= 600),1, 2) ./ diff(lambda.a(lambda.a > 560 & lambda.a <= 600));
    foo_difref = diff(p.ap(:, lambda.a > 500 & lambda.a <= 560),1, 2) ./ diff(lambda.a(lambda.a > 500 & lambda.a <= 560));
    toflag = any(abs(foo_difstep) > 3 * median(abs(foo_difref), 2, 'omitnan'), 2);
    if sum(toflag)
      fprintf('Signal of bubbles detected in %.2f%% (%i) of the spectra, flagged: step in ap between 550 and 600 nm\n', ...
        sum(toflag) / size(p, 1) * 100, sum(toflag))
      flag.ap_bubbles(toflag) = true;
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
  flag.HH_G50_flag(p.HH_G50 < 0 | p.HH_G50 > 50) = true;
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
  flag.HH_G50_mphi_suspicious(p.HH_G50 < 0.3 | p.HH_mphi < -5 | p.HH_G50 > 20) = true;
  
  % set flag column
  p.flag_bit = set_flagbit(flag);
  
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
  %   [g.ag, g.cg] = ResidualTemperatureAndScatteringCorrection(g.ag, g.cg, lambda.ref);
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


%% Residual Temperature And Scattering Correction (Zaneveld)
function [ap_corr, cp_corr] = ResidualTemperatureAndScatteringCorrectionZaneveld(ap, cp, wl, psi)
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
  iref = find(730 <= wl, 1,'first'); % 715 730
  % If ACS spectra do not go up to 730 nm take the closest wavelength to 730 nm
  if isempty(iref); [~, iref] = max(wl); end % works as there is data in iNIR so lowest wavelength is 710
  
  % Initialize output arrays
  deltaT = NaN(size(ap,1),1);
  
  % Init routine parameters
  bp = cp - ap;
  
  % Run minimization routine on good spectra only
  sel = find(all(isfinite(ap),2));
  for k = sel'
    deltaT(k) = fminsearch(@costFun_RTSC, 0, opts, ap(k,:), bp(k,:), psiT, iNIR, iref);         
  end
  ap_corr = ap - psiT.*deltaT - ((ap(:,iref) - psiT(iref).*deltaT) ./ bp(:,iref)) .* bp; 
  cp_corr = cp - psiT.*deltaT;
end


%% Residual Temperature And Scattering Correction (Kostakis)
function [ap_corr, cp_corr] = ResidualTemperatureAndScatteringCorrectionKostakis(ap, cp, wl, psi)
  % Function from Emmanuel Boss after Kostakis et al. 2022
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
  iref = find(730 <= wl, 1,'first'); % 715 730
  % If ACS spectra do not go up to 730 nm take the closest wavelength to 730 nm
  if isempty(iref); [~, iref] = max(wl); end % works as there is data in iNIR so lowest wavelength is 710
  
  % Initialize output arrays
  deltaT = NaN(size(ap,1),1);
  
  % Init routine parameters
  bp = cp - ap;
  
  % Run minimization routine on good spectra only
  sel = find(all(isfinite(ap),2));
  for k = sel'
    deltaT(k) = fminsearch(@costFun_RTSC, 0, opts, ap(k,:), bp(k,:), psiT, iNIR, iref);         
  end
  % split up temperature correction base and exponent and replace negative base by NaN to avoid complex solution leakage
  Tb = 0.212*(ap(:,iref) - psiT(iref).*deltaT);
  Tb(Tb < 0) = NaN;
  % apply correction
  ap_corr = ap - psiT.*deltaT - (ap(:,iref) - psiT(iref).*deltaT) + Tb.^1.135;
  cp_corr = cp - psiT.*deltaT;
end

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
  ap_filled = fillmissing(p.ap, 'nearest', 2);
  ap_sd_filled = fillmissing(p.ap_sd, 'nearest', 2);
  
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
  
  %   % vectorized partition_ap (DOESN'T WORK)
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

function [filt_interp] = agcg_fdom_interpolation(filt_interp, filt, lambda, cdom_dark, filt_avg)
  % Reconstruct ag and cg from fCDOM data
  [filt_avg, filt_interp, ~] = fdom_agcg_model(filt, filt_interp, lambda, filt_avg); % regression_stats
  n_periods = size(filt_avg,1)-1;
  filt_interp.lin_bool = false(size(filt_interp.fdom));

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
  %     % interpolation based on fCDOM - fCDOM_dark > 20% of fCDOM_dark
  %     if any((filt_interp.fdom(it) - cdom_dark) / cdom_dark > 0.2)
  %       % prepare delta a / delta f and delta c / delta f
  %       delta_acf.dt(i) = mean([filt_avg.dt(it0) filt_avg.dt(it1)], 'omitnan');
  %       delta_acf.delta_af(i, :) = (filt_avg.a(it1,:) - filt_avg.a(it0,:)) / (filt_avg.fdom_10smooth(it1) - filt_avg.fdom_10smooth(it0));
  %       delta_acf.delta_cf(i, :) = (filt_avg.c(it1,:) - filt_avg.c(it0,:)) / (filt_avg.fdom_10smooth(it1) - filt_avg.fdom_10smooth(it0));
  %       % prepare fdom0, a0 and c0
  %       filt_interp.fdom0(it) = repmat(filt_avg.fdom_10smooth(it0), sum(it), 1);
  %       filt_interp.a0(it,:) = repmat(filt_avg.a(it0,:), sum(it), 1);
  %       filt_interp.c0(it,:) = repmat(filt_avg.c(it0,:), sum(it), 1);
  %       filt_interp.lin_bool(it) = false;
  %     else
  %       filt_interp.lin_bool(it) = true;
  %     end
  %   end
  % end
  % delta_acf(all(isnan(delta_acf.delta_af), 2) | all(isnan(delta_acf.delta_cf), 2), :) = [];
  % % interpolated delta a / delta fdom over total events
  % filt_interp.delta_af = interp1(delta_acf.dt, delta_acf.delta_af, filt_interp.dt, 'linear');
  % filt_interp.delta_af = fillmissing(filt_interp.delta_af, 'nearest');
  % filt_interp.delta_cf = interp1(delta_acf.dt, delta_acf.delta_cf, filt_interp.dt, 'linear');
  % filt_interp.delta_cf = fillmissing(filt_interp.delta_cf, 'nearest');
  % % extrapolate a0, c0, fdom0
  % filt_interp.fdom0 = fillmissing(filt_interp.fdom0, 'nearest');
  % filt_interp.a0 = fillmissing(filt_interp.a0, 'nearest');
  % filt_interp.c0 = fillmissing(filt_interp.c0, 'nearest');
  % % compute ag and cg
  % filt_interp.a = filt_interp.a0 + filt_interp.delta_af .* (filt_interp.fdom - filt_interp.fdom0);
  % filt_interp.c = filt_interp.c0 + filt_interp.delta_cf .* (filt_interp.fdom - filt_interp.fdom0);
  % % linearly interpolate ag and cg
  % filt_interp.lin_a = interp1(filt_avg.dt, filt_avg.a, filt_interp.dt, 'linear');
  % filt_interp.lin_c = interp1(filt_avg.dt, filt_avg.c, filt_interp.dt, 'linear');
  % % replace fdom interpolation by linear for filter events where fdom is constant
  % filt_interp.a(filt_interp.lin_bool) = filt_interp.lin_a(filt_interp.lin_bool);
  % filt_interp.c(filt_interp.lin_bool) = filt_interp.lin_c(filt_interp.lin_bool);

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
  for i=1:n_periods
    it0 = i; it1 = i + 1;
    it = filt_avg.dt(it0) <= filt_interp.dt & filt_interp.dt <= filt_avg.dt(it1);
    if any(it)
      it_filt_interp = filt_interp(it,:);
      % interpolation based on fCDOM - fCDOM_dark > 20% of fCDOM_dark
      if any((filt_interp.fdom(it) - cdom_dark) / cdom_dark > 0.5)
        % compute interpolated ag factor
        ag_fact0 = (filt_avg.a(it0,:) - filt_avg.intercept_a(it0,:)) ./ filt_avg.fdom_10smooth(it0);
        ag_fact1 = (filt_avg.a(it1,:) - filt_avg.intercept_a(it1,:)) ./ filt_avg.fdom_10smooth(it1);
        ag_fact = interp1(filt_avg.dt(it0:it1,:), [ag_fact0; ag_fact1], it_filt_interp.dt, 'linear');
        % % compute interpolated cg factor
        % cg_fact0 = (filt_avg.c(it0,:) - filt_avg.intercept_c(it0,:)) ./ filt_avg.fdom_10smooth(it0);
        % cg_fact1 = (filt_avg.c(it1,:) - filt_avg.intercept_c(it1,:)) ./ filt_avg.fdom_10smooth(it1);
        % cg_fact = interp1(filt_avg.dt(it0:it1,:), [cg_fact0; cg_fact1], it_filt_interp.dt, 'linear');
        % compute ag and cg
        filt_interp.a(it, :) = filt_avg.a(it0,:) + ag_fact .* (filt_interp.fdom(it,:) - filt_avg.fdom_10smooth(it0,:));
        filt_interp.c(it, :) = filt_avg.c(it0,:) + ag_fact .* (filt_interp.fdom(it,:) - filt_avg.fdom_10smooth(it0,:));
        filt_interp.lin_bool(it) = false;
      else
        % linearly interpolate ag and cg if fdom is constant between it0 and it1
        filt_interp.a(it, :) = interp1(filt_avg.dt(it0:it1,:), filt_avg.a(it0:it1,:), it_filt_interp.dt, 'linear');
        filt_interp.c(it, :) = interp1(filt_avg.dt(it0:it1,:), filt_avg.c(it0:it1,:), it_filt_interp.dt, 'linear');
        filt_interp.lin_bool(it) = true;
      end
    end
  end

  % figure; hold on
  % scatter(datetime(filt_interp.dt, 'convertfrom', 'datenum'), ag_fact_tot(:,70), 20, 'ko')
  % scatter(datetime(filt_interp.dt(~filt_interp.lin_bool), 'convertfrom', 'datenum'), filt_interp.delta_af(~filt_interp.lin_bool,70), 20, 'r+')
  % 
  % figure; hold on
  % scatter(filt_interp.dt, cg_fact_tot(:,70), 20, 'ko')
  % scatter(filt_interp.dt(~filt_interp.lin_bool), filt_interp.delta_cf(~filt_interp.lin_bool,70), 20, 'r+')

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
  % cdom = table();
  % cdom.dt = datetime([filt.dt; filt_interp.dt], 'ConvertFrom', 'datenum');
  % cdom.fdom = [filt.fdom; filt_interp.fdom];
  % 
  % figure;
  % scatter(cdom.dt, cdom.fdom, 20, 'filled')



  % idp = filt_interp.dt >= min(it_filt_interp.dt) & filt_interp.dt <= max(it_filt_interp.dt);
  % figure(i+7)
  % scatter(datetime(filt_interp.dt(idp), 'ConvertFrom', 'datenum'), filt_interp.fdom(idp), 20, 'filled')
  % ylabel('fdom smoothed')
  % 
  % wl_toplot = 40;
  % % figure(1)
  % figure(i+10)
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
  % scatter(datetime(cdom.dt, 'ConvertFrom', 'datenum'), ...
  %   cdom.fdom, 40, 'filled', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5, 'Marker', 'd')
  % xlim([plot_dt(1)-hours(0.5) plot_dt(end)+hours(0.5)])
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
  % scatter(datetime(cdom.dt, 'ConvertFrom', 'datenum'), ...
  %   cdom.fdom, 40, 'filled', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5, 'Marker', 'd')
  % xlim([plot_dt(1)-hours(0.5) plot_dt(end)+hours(0.5)])
  % ylabel('fdom')
  % legend('cg fdom interpolated', 'cg filter average', 'fdom smoothed', 'Location', 'Best')

  % linearly interpolate a and c
  filt_interp.lin_a = interp1(filt_avg.dt, filt_avg.a, filt_interp.dt, 'linear');
  filt_interp.lin_c = interp1(filt_avg.dt, filt_avg.c, filt_interp.dt, 'linear');

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
function [filt_avg, filt_interp, regress_stats] = fdom_agcg_model(filt, filt_interp, lambda, filt_avg)
  


  % filt.dt = datetime(filt.dt, 'ConvertFrom', 'datenum');
  % cluster
  clusters = table();
  % clusters.a = NaN(size(filt.a));
  % clusters.c = NaN(size(filt.c));
  % for i = id450_a
  %   clusters.a(:,i) = kmeans(filt.a(:, i) ./ filt.fdom, 3);
  %   clusters.c(:,i) = kmeans(filt.c(:, i) ./ filt.fdom, 3);
  %   % figure; hold on
  %   % gscatter(filt.a(:,i), filt.fdom, clusters.a(:,i), [], 'xod')
  % end


  % id450_a = find(lambda.a > 450 & lambda.a < 550);
  id450_a = abs(lambda.a - 450) == min(abs(lambda.a - 450));
  id450_c = abs(lambda.c - 450) == min(abs(lambda.c - 450));
  
  % number of clusters to try
  klist = 1:5;
  kmeanfunc = @(X,K)(kmeans(X, K));
  % evaluate optimum number of cluster for a
  eva = evalclusters(filt.a(:, id450_a) ./ filt.fdom, kmeanfunc, 'silhouette','klist', klist);
  clusters.a = kmeans(filt.a(:, id450_a) ./ filt.fdom, eva.OptimalK);
  % evaluate optimum number of cluster for a
  eva = evalclusters(filt.c(:, id450_c) ./ filt.fdom, kmeanfunc, 'silhouette','klist', klist);
  clusters.c = kmeans(filt.c(:, id450_c) ./ filt.fdom, eva.OptimalK);

  % col = distinguishable_colors(3);

  % figure; hold on
  % scatter(lambda.a, filt.a(clusters.a == 1, :) ./ filt.fdom(clusters.a == 1), ...
  %   'MarkerEdgeColor', col(1, :), 'MarkerEdgeAlpha', 0.2, 'Marker', 'x')
  % scatter(lambda.a, filt.a(clusters.a == 2, :) ./ filt.fdom(clusters.a == 2), ...
  %   'MarkerEdgeColor', col(2, :), 'MarkerEdgeAlpha', 0.2, 'Marker', 'o')
  % scatter(lambda.a, filt.a(clusters.a == 3, :) ./ filt.fdom(clusters.a == 3), ...
  %   'MarkerEdgeColor', col(3, :), 'MarkerEdgeAlpha', 0.2, 'Marker', 'd')
  % 
  % figure; hold on
  % scatter(lambda.c, filt.c(clusters.c == 1, :) ./ filt.fdom(clusters.c == 1), ...
  %   'MarkerEdgeColor', col(1, :), 'MarkerEdgeAlpha', 0.2, 'Marker', 'x')
  % scatter(lambda.c, filt.c(clusters.c == 2, :) ./ filt.fdom(clusters.c == 2), ...
  %   'MarkerEdgeColor', col(2, :), 'MarkerEdgeAlpha', 0.2, 'Marker', 'o')
  % scatter(lambda.c, filt.c(clusters.c == 3, :) ./ filt.fdom(clusters.c == 3), ...
  %   'MarkerEdgeColor', col(3, :), 'MarkerEdgeAlpha', 0.2, 'Marker', 'd')

  % regress a&c with fdom
  regress_stats = struct();
  [regress_stats.a, regress_stats.c] = regress_acfilt(filt, lambda, filt, 'fdom', clusters);

  % find wavelength where fdom fit and temperature fit work poorly
  bad_fit_a = all(regress_stats.a.R2 < prctile(regress_stats.a.R2, 5, 1), 2);
  bad_fit_c = all(regress_stats.c.R2 < prctile(regress_stats.c.R2, 5, 1), 2);
  
  %%% apply the regression to get a and c dissolved from fdom
  groups_a = unique(clusters.a(~isnan(clusters.a)));
  groups_c = unique(clusters.c(~isnan(clusters.c)));

  % Weight cluster for each filter event
  cluster_weighted = table(filt_avg.dt, 'VariableNames', {'dt'});
  cluster_weighted.a = NaN(size(cluster_weighted,1), size(groups_a, 1));
  cluster_weighted.c = NaN(size(cluster_weighted,1), size(groups_c, 1));
  for i=1:size(filt_avg, 1)
    sel_filt = filt_avg.start(i) <= filt.dt & filt.dt <= filt_avg.end(i);
    foo = clusters(sel_filt,:);
    for j = 1:size(groups_a,1)
      cluster_weighted.a(i,j) = sum(foo.a == groups_a(j)) ./ sum(~isnan(foo.a));
    end
    for j = 1:size(groups_c,1)
      cluster_weighted.c(i,j) = sum(foo.c == groups_c(j)) ./ sum(~isnan(foo.c));
    end
  end
  % weigth a slope and intercept
  cluster_weighted.slope_a = zeros(size(cluster_weighted, 1), size(regress_stats.a, 1));
  cluster_weighted.intercept_a = zeros(size(cluster_weighted, 1), size(regress_stats.a, 1));
  for j = 1:size(groups_a, 1)
    cluster_weighted.slope_a = cluster_weighted.slope_a + regress_stats.a.slope(:,j)'.*cluster_weighted.a(:,j);
    cluster_weighted.intercept_a = cluster_weighted.intercept_a + regress_stats.a.intercept(:,j)'.*cluster_weighted.a(:,j);
  end
  % weigth c slope and intercept
  cluster_weighted.slope_c = zeros(size(cluster_weighted, 1), size(regress_stats.c, 1));
  cluster_weighted.intercept_c = zeros(size(cluster_weighted, 1), size(regress_stats.c, 1));
  for j = 1:size(groups_c, 1)
    cluster_weighted.slope_c = cluster_weighted.slope_c + regress_stats.c.slope(:,j)'.*cluster_weighted.c(:,j);
    cluster_weighted.intercept_c = cluster_weighted.intercept_c + regress_stats.c.intercept(:,j)'.*cluster_weighted.c(:,j);
  end
  % interpolate slope and intercept cluster specific onto filt_interp.dt
  filt_interp.slope_interp_a = interp1(cluster_weighted.dt, cluster_weighted.slope_a, filt_interp.dt, 'spline');
  filt_interp.slope_interp_a = fillmissing(filt_interp.slope_interp_a, 'nearest');
  filt_interp.intercept_interp_a = interp1(cluster_weighted.dt, cluster_weighted.intercept_a, filt_interp.dt, 'spline');
  filt_interp.intercept_interp_a = fillmissing(filt_interp.intercept_interp_a, 'nearest');
  filt_interp.slope_interp_c = interp1(cluster_weighted.dt, cluster_weighted.slope_c, filt_interp.dt, 'spline');
  filt_interp.slope_interp_c = fillmissing(filt_interp.slope_interp_c, 'nearest');
  filt_interp.intercept_interp_c = interp1(cluster_weighted.dt, cluster_weighted.intercept_c, filt_interp.dt, 'spline');
  filt_interp.intercept_interp_c = fillmissing(filt_interp.intercept_interp_c, 'nearest');
  
  filt_avg.slope_a = cluster_weighted.slope_a;
  filt_avg.intercept_a = cluster_weighted.intercept_a;
  filt_avg.slope_c = cluster_weighted.slope_c;
  filt_avg.intercept_c = cluster_weighted.intercept_c;
  
  % % reconstruct a filtered from fdom
  % % TODO remember here I've removed the bad fit component
  % % filt_interp.a = filt_interp.fdom .* filt_interp.slope_interp_a(:, ~bad_fit_a) + filt_interp.intercept_interp_a(:, ~bad_fit_a);
  % % filt_interp.a(:, ~bad_fit_a) = filt_interp.fdom .* filt_interp.slope_interp_a(:, ~bad_fit_a) + filt_interp.intercept_interp_a(:, ~bad_fit_a);
  % filt_interp.a = filt_interp.fdom .* filt_interp.slope_interp_a + filt_interp.intercept_interp_a;
  % % reconstruct c filtered from fdom
  % % filt_interp.c = filt_interp.fdom .* filt_interp.slope_interp_c(:, ~bad_fit_c) + filt_interp.intercept_interp_c(:, ~bad_fit_c);
  % % filt_interp.c(:, ~bad_fit_c) = filt_interp.fdom .* filt_interp.slope_interp_c(:, ~bad_fit_c) + filt_interp.intercept_interp_c(:, ~bad_fit_c);
  % filt_interp.c = filt_interp.fdom .* filt_interp.slope_interp_c + filt_interp.intercept_interp_c;


  % for i = 1:size(lambda.a, 2)
  %   % reconstruct a filtered from fdom
  %   % TODO remember here I've removed the bad fit component
  %   for j = 1:size(groups_a, 1)
  %     id_grp_a = clusters.a == groups_a(j);
  %     % filt_interp.a(id_grp_a, :) = filt_interp.fdom(id_grp_a) .* regress_stats.a.slope(~bad_fit_a,j)' + repmat(regress_stats.a.intercept(~bad_fit_a,j)', sum(id_grp_a), 1);
  %     % filt_interp.a(id_grp_a, ~bad_fit_a) = filt_interp.fdom(id_grp_a) .* regress_stats.a.slope(~bad_fit_a,j)' + repmat(regress_stats.a.intercept(~bad_fit_a,j)', sum(id_grp_a), 1);
  %     filt_interp.a(id_grp_a, :) = filt_interp.fdom(id_grp_a) .* regress_stats.a.slope(:,j)' + repmat(regress_stats.a.intercept(:,j)', sum(id_grp_a), 1);
  %   end
  %   % reconstruct c filtered from fdom
  %   for j = 1:size(groups_c, 1)
  %     id_grp_c = clusters.c == groups_c(j);
  %     % filt_interp.c(id_grp_c, :) = filt_interp.fdom(id_grp_c) .* regress_stats.c.slope(~bad_fit_c,j)' + repmat(regress_stats.c.intercept(~bad_fit_c,j)', sum(id_grp_c), 1);
  %     % filt_interp.c(id_grp_c, ~bad_fit_c) = filt_interp.fdom(id_grp_c) .* regress_stats.c.slope(~bad_fit_c,j)' + repmat(regress_stats.c.intercept(~bad_fit_c,j)', sum(id_grp_c), 1);
  %     filt_interp.c(id_grp_c, :) = filt_interp.fdom(id_grp_c) .* regress_stats.c.slope(:,j)' + repmat(regress_stats.c.intercept(:,j)', sum(id_grp_c), 1);
  %   end
  % end


  % % % reconstruct a and c filtered from fdom
  % % % TODO remember here I've removed the bad fit component
  % % % filt_interp.a = filt_interp.fdom .* regress_stats.a.slope(~bad_fit)' + repmat(regress_stats.a.intercept(~bad_fit)', size(tot.a, 1), 1);
  % % filt_interp.a = filt_interp.fdom .* regress_stats.a.slope' + repmat(regress_stats.a.intercept', size(tot.a, 1), 1);
  % % % filt_interp.c = filt_interp.fdom .* regress_stats.c.slope(~bad_fit)' + repmat(regress_stats.c.intercept(~bad_fit)', size(tot.c, 1), 1);
  % % filt_interp.c = filt_interp.fdom .* regress_stats.c.slope' + repmat(regress_stats.c.intercept', size(tot.c, 1), 1);
  % % % filt_interp.a(:, ~bad_fit) = filt_interp.fdom .* regress_stats.a.slope(~bad_fit)' + repmat(regress_stats.a.intercept(~bad_fit)', size(tot.a, 1), 1);
  % % % filt_interp.c(:, ~bad_fit) = filt_interp.fdom .* regress_stats.c.slope(~bad_fit)' + repmat(regress_stats.c.intercept(~bad_fit)', size(tot.c, 1), 1);


  % % fill missing modelled data for bad fit wavelength with spline interpolation
  % filt_interp.a(:, bad_fit) = interp1(lambda.a, filt_interp.a', lambda.a(bad_fit), 'spline')';
  % filt_interp.c(:, bad_fit) = interp1(lambda.c, filt_interp.c', lambda.c(bad_fit), 'spline')';


  % % compare filter even measured and modelled
  %     % compute a&c filt avg from best regression (fcdom / temperature / fcdom*temperature)
  %     filt_avg_mod = table(filt_avg.dt, 'VariableNames', {'dt'});
  %     filt_avg_mod.a = NaN(size(filt_avg.a));
  %     filt_avg_mod.c = NaN(size(filt_avg.c));
  %     % fill regression variable missing data with spline interpolation and nearest extrapolation
  %     filt_avg_mod.fdom = interp1(cdom.dt, cdom.fdom, filt_avg.dt, 'spline');
  %     filt_avg_mod.fdom(isnan(filt_avg_mod.fdom)) = interp1(filt_avg_mod.dt, filt_avg_mod.fdom, filt_avg.dt(isnan(filt_avg_mod.fdom)), 'nearest', 'extrap');
  %     % reconstruct a and c filtered from fdom
  %     filt_avg_mod.a = filt_avg_mod.fdom .* regress_stats.a.slope' + repmat(regress_stats.a.intercept', size(filt_avg.a, 1), 1);
  %     filt_avg_mod.c = filt_avg_mod.fdom .* regress_stats.c.slope' + repmat(regress_stats.c.intercept', size(filt_avg.c, 1), 1);
  %     % filt_avg_mod.a(:, ~bad_fit) = filt_avg_mod.fdom .* regress_stats.a.slope(~bad_fit)' + repmat(regress_stats.a.intercept(~bad_fit)', size(filt_avg.a, 1), 1);
  %     % filt_avg_mod.c(:, ~bad_fit) = filt_avg_mod.fdom .* regress_stats.c.slope(~bad_fit)' + repmat(regress_stats.c.intercept(~bad_fit)', size(filt_avg.c, 1), 1);
  %     % 
  %     % % fill missing modelled data for bad fit wavelength with spline interpolation
  %     % filt_avg_mod.a(:, bad_fit) = interp1(lambda.a, filt_avg_mod.a', lambda.a(bad_fit), 'spline')';
  %     % filt_avg_mod.c(:, bad_fit) = interp1(lambda.c, filt_avg_mod.c', lambda.c(bad_fit), 'spline')';
  % 
  %     figure;
  %     hold on
  %     scatter(lambda.a, regress_stats.a.slope, 20, [113	113	198]/255, 'filled')
  %     scatter(lambda.c, regress_stats.c.slope, 20, [113	113	198]/255)
  %     legend('slope of c_{filt} to fCDOM', 'slope of a_{filt} to fCDOM', 'FontSize', 16, 'Location', 'Northwest')
  %     ylabel('slope of a_{filt} and c_{filt} to fCDOM', 'FontSize', 16)
  %     xlabel("wavelength (nm)")
  % 
  %     % % percentage difference between modelled a and c from fdom and a and c
  %     filt_diff = table();
  %     filt_diff.a = abs((filt_avg_mod.a - filt_avg.a) ./ filt_avg.a) * 100;
  %     filt_diff.c = abs((filt_avg_mod.c - filt_avg.c) ./ filt_avg.c) * 100;
  %     avg_diff_a = mean(filt_diff.a, 'omitnan');
  %     avg_diff_c = mean(filt_diff.c, 'omitnan');
  %     std_diff_a = std(filt_diff.a, 'omitnan');
  %     std_diff_c = std(filt_diff.c, 'omitnan');
  % 
  %     figure; hold on
  %     plot(lambda.a, avg_diff_a, '-b')
  %     plot(lambda.c, avg_diff_c, '-r')
  %     ylabel('percentage difference a&c filtered fdom and measured')
  %     legend('percentage difference a', 'percentage difference c')
  % 
  %     visProd3D(lambda.a, filt_avg_mod.dt, filt_avg_mod.a, false, 'Wavelength', false, 23);
  %     visProd3D(lambda.c, filt_avg_mod.dt, filt_avg_mod.c, false, 'Wavelength', false, 24);
  % 
  %     visProd3D(lambda.a, filt_avg.dt, filt_avg.a, false, 'Wavelength', false, 25);
  %     visProd3D(lambda.c, filt_avg.dt, filt_avg.c, false, 'Wavelength', false, 26);
  % 
  %     visProd3D(lambda.a, filt_interp.dt, filt_interp.a, false, 'Wavelength', false, 23);
  %     visProd3D(lambda.c, filt_interp.dt, filt_interp.c, false, 'Wavelength', false, 24);
  % 
  %     visProd3D(lambda.a, filt_interp.dt, filt_interp.a, false, 'Wavelength', false, 25);
  %     visProd3D(lambda.c, filt_interp.dt, filt_interp.c, false, 'Wavelength', false, 26);
end


%% regression between a/c and other variable in filter events
function [a_reg, c_reg] = regress_acfilt(filt, lambda, ancillary, varr, clusters)
  if nargin < 5
    clusters = table();
    clusters.a = ones(size(filt, 1), 1);
    clusters.c = ones(size(filt, 1), 1);
  end
  % Interpolate regression variable over filtered datetime
  filt.([varr '_b']) = interp1(ancillary.dt, ancillary.(varr), ancillary.dt, 'linear');
  % robust linear regression between filt a&c and ancillary variable
  a_reg = array2table(NaN(size(lambda.a, 2), 5), 'VariableNames', {'slope', 'intercept', 'RMSE', 'nRMSE', 'R2'});
  c_reg = array2table(NaN(size(lambda.a, 2), 5), 'VariableNames', {'slope', 'intercept', 'RMSE', 'nRMSE', 'R2'});
  groups_a = unique(clusters.a(~isnan(clusters.a)));
  groups_c = unique(clusters.c(~isnan(clusters.c)));
  for i = 1:size(lambda.a, 2)
    for j = 1:size(groups_a, 1)
      id_grp_a = clusters.a == groups_a(j);
      % regress variable with filtered water absorption
      [stats.b, stats.stats] = robustfit(ancillary.(varr)(id_grp_a), filt.a(id_grp_a, i));
      a_reg.slope(i, j) = stats.b(2);
      a_reg.intercept(i, j) = stats.b(1);
      a_reg.RMSE(i, j) = stats.stats.robust_s;
      % calculated normalize RMSE by average and R2
      id_nonan = all(~isnan([filt.a(:, i) ancillary.(varr)]), 2);
      a_reg.nRMSE(i, j) = a_reg.RMSE(i, j)/abs(mean(filt.a(id_grp_a & id_nonan, i), 'omitnan'))*100;
      a_reg.R2(i, j) = corr(filt.a(id_grp_a & id_nonan, i), stats.b(1) + stats.b(2) * ancillary.(varr)(id_grp_a & id_nonan))^2;
    end
    for j = 1:size(groups_c, 1)
      id_grp_c = clusters.c == groups_c(j);
      % regress variable with filtered water attenuation
      [stats.b, stats.stats] = robustfit(ancillary.(varr)(id_grp_c), filt.c(id_grp_c, i));
      c_reg.slope(i, j) = stats.b(2);
      c_reg.intercept(i, j) = stats.b(1);
      c_reg.RMSE(i, j) = stats.stats.robust_s;
      % calculated normalize RMSE by average and R2
      id_nonan = all(~isnan([filt.c(:, i) ancillary.(varr)]), 2);
      c_reg.nRMSE(i, j) = c_reg.RMSE(i, j)/abs(mean(filt.c(id_grp_c & id_nonan, i), 'omitnan'))*100;
      c_reg.R2(i, j) = corr(filt.c(id_grp_c & id_nonan, i), stats.b(1) + stats.b(2) * ancillary.(varr)(id_grp_c & id_nonan))^2;
    end
  end

  % figure;
  % subplot(2,4,1)
  % plot_linreg(ancillary.(varr), filt.a(:, 76), 'robust', 'linear', false, 0.05);
  % xlabel(varr)
  % ylabel(['a' num2str(lambda.a(76))])
  % subplot(2,4,2)
  % plot_linreg(ancillary.(varr), filt.a(:, 78), 'robust', 'linear', false, 0.05);
  % xlabel(varr)
  % ylabel(['a' num2str(lambda.a(78))])
  % subplot(2,4,3)
  % plot_linreg(ancillary.(varr), filt.a(:, 79), 'robust', 'linear', false, 0.05);
  % xlabel(varr)
  % ylabel(['a' num2str(lambda.a(79))])
  % subplot(2,4,4)
  % plot_linreg(ancillary.(varr), filt.a(:, 80), 'robust', 'linear', false, 0.05);
  % xlabel(varr)
  % ylabel(['a' num2str(lambda.a(80))])
  % subplot(2,4,5)
  % plot_linreg(ancillary.(varr), filt.a(:, 81), 'robust', 'linear', false, 0.05);
  % xlabel(varr)
  % ylabel(['a' num2str(lambda.a(81))])
  % subplot(2,4,6)
  % plot_linreg(ancillary.(varr), filt.a(:, 82), 'robust', 'linear', false, 0.05);
  % xlabel(varr)
  % ylabel(['a' num2str(lambda.a(82))])
  % subplot(2,4,7)
  % plot_linreg(ancillary.(varr), filt.a(:, 84), 'robust', 'linear', false, 0.05);
  % xlabel(varr)
  % ylabel(['a' num2str(lambda.a(84))])
  % subplot(2,4,8)
  % plot_linreg(ancillary.(varr), filt.a(:, 86), 'robust', 'linear', false, 0.05);
  % xlabel(varr)
  % ylabel(['a' num2str(lambda.a(86))])
end

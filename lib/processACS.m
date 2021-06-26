function [p, g, bad] = processACS(lambda, tot, filt, modelG50, modelmphi, di, ...
  cdom, fth, fth_constants, interpolation_method, di_method, compute_ad_aphi)
% NOTE: wavelength of c are interpolated to wavelength of a
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

% interpolate fth.swt onto binned data to fill missing flow data
fth_interp = table([tot.dt; fth.dt; filt.dt], 'VariableNames', {'dt'});
[~,b] = sort(fth_interp.dt); % sort dates
fth_interp.dt = fth_interp.dt(b,:);
fth_interp.swt = interp1(fth.dt, fth.swt, fth_interp.dt, 'previous');%, 'linear', 'extrap');
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
switch interpolation_method
  case 'CDOM'
    % Require both CDOM & Switch position
    % Parameters for CDOM signal (in units of cdom.fdom, by default it's counts);
    param_extrap_threshold = 4;  % counts
    param_min_variability = 1;   % counts
    % Make filtered period median
    filt_avg = table((fth_interp.dt(sel_start) + fth_interp.dt(sel_end)) ./ 2, 'VariableNames', {'dt'});
    filt_avg.cdom = NaN(size(filt_avg,1),1);
    filt_avg.a = NaN(size(filt_avg,1), size(lambda.a, 2));
    filt_avg.c = NaN(size(filt_avg,1), size(lambda.c, 2));
    for i=1:size(sel_start, 1)
      sel_cdom = fth_interp.dt(sel_start(i)) <= cdom.dt & cdom.dt <= fth_interp.dt(sel_end(i));
      filt_avg.cdom(i) = median(cdom.fdom(sel_cdom));
      sel_filt = fth_interp.dt(sel_start(i)) <= filt.dt & filt.dt <= fth_interp.dt(sel_end(i));
      filt_avg.a(i,:) = median(filt.a(sel_filt,:));
      filt_avg.c(i,:) = median(filt.c(sel_filt,:));
    end
    % Remove periods were there is no filtered spectrum of ACS
    filt_avg(any(isnan(filt_avg.a),2),:) = [];
    % use regressions as interpolation method
  %       warning('off', 'stats:regress:RankDefDesignMat');
  %       % Make CDOM function to fill gaps between 2 filtered periods
  %       n_periods = size(filt_avg,1)-1;
  %       cdom_fun = table(cdom.dt, 'VariableNames', {'dt'});
  %       cdom_fun.a = NaN(size(cdom_fun,1),n_wv);
  %       cdom_fun.c = NaN(size(cdom_fun,1),n_wv);
  %       slope_a = NaN(n_periods,n_wv); inter_a = NaN(n_periods,n_wv);
  %       slope_c = NaN(n_periods,n_wv); inter_c = NaN(n_periods,n_wv);
  %       for i=1:n_periods % For each total period
  %         is = i; ie = i + 1;
  %         % Compute regression for each wavelength
  %         for j=1:n_wv
  %           if any(isnan(filt_avg.a(is:ie, j)))
  %             slope_a(i,j) = NaN; inter_a(i,j) = NaN;
  %             slope_c(i,j) = NaN; inter_c(i,j) = NaN;
  %           else
  %             foo = regress(filt_avg.a(is:ie, j), [filt_avg.cdom(is:ie), [1;1]]);
  %             slope_a(i,j) = foo(1); inter_a(i,j) = foo(2);
  %             foo = regress(filt_avg.c(is:ie, j), [filt_avg.cdom(is:ie), [1;1]]);
  %             slope_c(i,j) = foo(1); inter_c(i,j) = foo(2);
  %           end
  %         end
  %         % Build a and c on cdom timestamp
  %         sel = filt_avg.dt(is) <= cdom.dt & cdom.dt <= filt_avg.dt(ie);
  %         cdom_fun.a(sel,:) = slope_a(i,:) .* cdom.fdom(sel) .* ones(1,n_wv) + inter_a(i,:);
  %         cdom_fun.c(sel,:) = slope_c(i,:) .* cdom.fdom(sel) .* ones(1,n_wv) + inter_c(i,:);
  %       end
  %     %   cdom_fun(any(isnan(cdom_fun.a),1),:) = []; % Remove NaN values
  %       % Interpolate cdom_fun on tot
  %       filt_interp = table(tot.dt, 'VariableNames', {'dt'});
  %       filt_interp.a = interp1(cdom_fun.dt, cdom_fun.a, filt_interp.dt);%, 'linear', 'extrap');
  %       filt_interp.c = interp1(cdom_fun.dt, cdom_fun.c, filt_interp.dt);%, 'linear', 'extrap');
  %       warning('on', 'stats:regress:RankDefDesignMat');
    % Use simple mathematical function to interpolate based on CDOM
    n_periods = size(filt_avg,1)-1;
    filt_interp = table(tot.dt, 'VariableNames', {'dt'});
    filt_interp.cdom = interp1(cdom.dt, cdom.fdom, tot.dt); % as independent from tot|filt period
    filt_interp.a = NaN(size(filt_interp,1),n_wv);
    filt_interp.c = NaN(size(filt_interp,1),n_wv);
    % For each period going from t0 to t1, starting and finishing by a filtered time
    for i=1:n_periods
      it0 = i; it1 = i + 1;
      it = filt_avg.dt(it0) <= filt_interp.dt & filt_interp.dt <= filt_avg.dt(it1);
      var_period = abs(filt_avg.cdom(it1) - filt_avg.cdom(it0));
      if max(filt_interp.cdom(it)) - min(filt_interp.cdom(it)) > var_period + param_extrap_threshold
        % If variability during total period is higher than the variability
        % between the two boundaries, extrapolation in needed 
        % Processing is ignored for now.
        filt_interp.a(it,:) = NaN;
        filt_interp.c(it,:) = NaN;
  %       fprintf('processACS: interpolation of filtered period requires an extrapolation of the CDOM signals which is not supported.');
        fprintf('processACS: Unable to interpolate filt %s-%s', datestr(filt_avg.dt(it0)), datestr(filt_avg.dt(it1)));
      elseif var_period > param_min_variability
        % Signficant variability in CDOM between 2 filtered periods
        % Use cdom signal to interpolate filt during tot periods
        Xt = (filt_avg.cdom(it1) - filt_interp.cdom(it)) / (filt_avg.cdom(it1) - filt_avg.cdom(it0));
        filt_interp.a(it,:) = Xt .* filt_avg.a(it0,:) + (1 - Xt) .* filt_avg.a(it1,:);
        filt_interp.c(it,:) = Xt .* filt_avg.c(it0,:) + (1 - Xt) .* filt_avg.c(it1,:);
      else
        % Variability in CDOM is not significant between 2 filtered periods
        % Apply linear interpolation to estimate filt during tot period
        filt_interp.a(it,:) = interp1(filt_avg.dt(it0:it1), filt_avg.a(it0:it1,:), filt_interp.dt(it));
        filt_interp.c(it,:) = interp1(filt_avg.dt(it0:it1), filt_avg.c(it0:it1,:), filt_interp.dt(it));
      end
    end
  %   % Interpolate a_filt and c_filt based on CDOM timestamp to tot timestamp
  %   filt_interp = table(tot.dt, 'VariableNames', {'dt'});
  %   filt_interp.a = interp1(filt_interp.dt, filt_interp.a, filt_interp.dt);%, 'linear', 'extrap');
  %   filt_interp.c = interp1(filt_interp.dt, filt_interp.c, filt_interp.dt);%, 'linear', 'extrap');

    % Note std is interpolated linearly (not using the cdom function)
    filt_interp.a_avg_sd = interp1(filt.dt, filt.a_avg_sd, filt_interp.dt);%, 'linear', 'extrap');
    filt_interp.c_avg_sd = interp1(filt.dt, filt.c_avg_sd, filt_interp.dt);%, 'linear', 'extrap');
    
  case 'linear'
    % Compute filtered period median
    filt_avg = table((fth_interp.dt(sel_start) + fth_interp.dt(sel_end)) ./ 2, 'VariableNames', {'dt'});
    filt_avg.a = NaN(size(filt_avg,1), size(lambda.a, 2));
    filt_avg.c = NaN(size(filt_avg,1), size(lambda.c, 2));
    filt_avg.a_avg_sd = NaN(size(filt_avg,1), size(lambda.a, 2));
    filt_avg.c_avg_sd = NaN(size(filt_avg,1), size(lambda.c, 2));
    filt_avg.a_avg_n = NaN(size(filt_avg,1), 1);
    filt_avg.c_avg_n = NaN(size(filt_avg,1), 1);
    for i=1:size(sel_start, 1)
      sel_filt = fth_interp.dt(sel_start(i)) <= filt.dt & filt.dt <= fth_interp.dt(sel_end(i));
      foo = filt(sel_filt,:);
      foo.a_avg_sd(foo.a > prctile(foo.a, 25, 1)) = NaN;
      foo.c_avg_sd(foo.c > prctile(foo.c, 25, 1)) = NaN;
      foo.a(foo.a > prctile(foo.a, 25, 1)) = NaN;
      foo.c(foo.c > prctile(foo.c, 25, 1)) = NaN;
      % compute average of all values smaller than 25th percentile for each filter event
      filt_avg.a(i,:) = mean(foo.a, 1, 'omitnan');
      filt_avg.c(i,:) = mean(foo.c, 1, 'omitnan');
      filt_avg.a_avg_sd(i,:) = mean(foo.a_avg_sd, 1, 'omitnan');
      filt_avg.c_avg_sd(i,:) = mean(foo.c_avg_sd, 1, 'omitnan');
      filt_avg.a_avg_n(i) = sum(foo.a_avg_n(any(~isnan(foo.a), 2)), 'omitnan');
      filt_avg.c_avg_n(i) = sum(foo.c_avg_n(any(~isnan(foo.c), 2)), 'omitnan');
    end
    filt_avg(all(isnan(filt_avg.a), 2) | all(isnan(filt_avg.c), 2), :) = [];
    % Interpolate filtered on total linearly
    filt_interp = table(tot.dt, 'VariableNames', {'dt'});
    filt_interp.a = interp1(filt_avg.dt, filt_avg.a, filt_interp.dt);%, 'linear', 'extrap');
    filt_interp.c = interp1(filt_avg.dt, filt_avg.c, filt_interp.dt);%, 'linear', 'extrap');
    filt_interp.a_avg_sd = interp1(filt_avg.dt, filt_avg.a_avg_sd, filt_interp.dt);%, 'linear', 'extrap');
    filt_interp.c_avg_sd = interp1(filt_avg.dt, filt_avg.c_avg_sd, filt_interp.dt);%, 'linear', 'extrap');
  otherwise
  error('Method not supported.');
end

if exist('visFlag', 'file') && exist('fth', 'var')
  fh = visFlag([], filt_interp, tot, [], filt_avg, [], 'a', round(size(tot.a, 2)/2), [], fth);
  title('Check filter event interpolation, press q to continue', 'FontSize', 14)
  legend('Filtered interpolated', 'Total', 'Filtered median', 'Flow rate',...
    'AutoUpdate','off', 'FontSize', 12)
  guiSelectOnTimeSeries(fh);
elseif exist('visFlag', 'file')
  fh = visFlag([], filt_interp, tot, [], filt_avg, [], 'a', round(size(tot.a, 2)/2), [], []);
  title('Check filter event interpolation, press q to continue', 'FontSize', 14)
  legend('Filtered interpolated', 'Total', 'Filtered median', 'Flow rate',...
    'AutoUpdate','off', 'FontSize', 12)
  guiSelectOnTimeSeries(fh);
end

% Remove lines of NaNs
sel2rm = any(isnan(tot.a),2) | any(isnan(tot.c),2) |...
         any(isnan(filt_interp.a),2) | any(isnan(filt_interp.c),2);
tot(sel2rm,:) = [];
filt_interp(sel2rm,:) = [];

% Particulate = Total - FSW
p = table(tot.dt, 'VariableNames', {'dt'});
p.ap = tot.a - filt_interp.a;
p.cp = tot.c - filt_interp.c;

if size(lambda.a, 2) > 50 % perform too separate correction for ap and cp only for ACS data, not AC9
  % Interpolate wavelengths for Scattering & Residual temperature correction
  ap_for_cpresiduals_corr = interp1(lambda.a', p.ap', lambda.c', 'linear', 'extrap')';
  cp_for_apresiduals_corr = interp1(lambda.c', p.cp', lambda.a', 'linear', 'extrap')';
  % ap Scattering & Residual temperature correction
  % cp Residual correction (for efficiency use the one computed from ap as it should be the same)
  fprintf('ap residual temperature and scattering correction ... ')
  [p.ap, ~] = ResidualTemperatureAndScatteringCorrection(p.ap, cp_for_apresiduals_corr, lambda.a);
  fprintf('Done\n')
  fprintf('cp residual temperature and scattering correction ... ')
  [~, p.cp] = ResidualTemperatureAndScatteringCorrection(ap_for_cpresiduals_corr, p.cp, lambda.c);
  fprintf('Done\n')
else
  fprintf('ap and cp residual temperature and scattering corrections ... ')
  [p.ap, p.cp] = ResidualTemperatureAndScatteringCorrection(p.ap, p.cp, lambda.ref);
  fprintf('Done\n')
end

% Propagate error (using geometric mean of measurement errors)
%   Note: Error is not propagated through Scattering & Residual temperature
%         correction as required by SeaBASS
p.ap_sd = sqrt(tot.a_avg_sd .* filt_interp.a_avg_sd);
p.cp_sd = sqrt(tot.c_avg_sd .* filt_interp.c_avg_sd);
p.ap_n = tot.a_avg_n;
p.cp_n = tot.c_avg_n;

% delete ap spectrum full of NaNs
p(all(isnan(p.ap),2),:) = [];
p(all(isnan(p.cp),2),:) = [];

% Unsmoothing ACS spectrum
% Ron Zaneveld, WET Labs, Inc., 2005
if size(lambda.a, 2) > 50 % unsmooth only ACS, not AC9
  p = unsmoothACS(p, lambda);
%   p_unsmooth = unsmoothACS(p, lambda);
end

% QC using ap spectrum
wla_430 = lambda.a(find(lambda.a <= 430, 1,'last')); % find lower and closest to 430nm wavelength
wla_700 = lambda.a(find(lambda.a >= 700, 1,'first')); % find higher and closest to 700nm wavelength

% delete absorption values < -0.0015
p.ap_sd(p.ap < -0.0015 & lambda.a < wla_430) = NaN;
p.ap_sd(p.ap < -0.0015 & lambda.a >= wla_700) = NaN;
p.ap(p.ap < -0.0015 & lambda.a < wla_430) = NaN;
p.ap(p.ap < -0.0015 & lambda.a >= wla_700) = NaN;
todelete = any(p.ap < -0.0015 & lambda.a >= wla_430 & lambda.a <= wla_700 | ...
  p.cp < -0.0015 | p.cp > 10, 2);
fprintf('%.2f%% (%i) spectrum failed auto-QC step 1: ap 430-700 | cp < -0.0015 | p.cp > 10\n', ...
  sum(todelete) / size(p, 1) * 100, sum(todelete))
bad = [p(todelete, :) table(repmat({'ap 430-700  | cp < -0.0015 | p.cp > 10'}, ...
  sum(todelete), 1), 'VariableNames', {'QC_failed'})];
p(todelete, :) = [];

% find wavelength below and above which ap and cp are unrealistic and replace by NaNs
if size(lambda.a, 2) > 50 % clean only ACS data, not AC9
  % replace unrealistic ap and cp red wavelenght by NaN
  [~, foo] = min(p.cp(:, lambda.c > 700),[],2);
  minred_wlc = sum(lambda.c <= 700) + foo;
  [~, foo] = min(abs(p.ap(:, lambda.a > 700)),[],2);
  minred_wla = sum(lambda.a <= 700) + foo;
  p.cp(lambda.c > lambda.c(minred_wlc)' & lambda.c > 710) = NaN;
  p.cp_sd(lambda.c > lambda.c(minred_wlc)' & lambda.c > 710) = NaN;
  p.ap(lambda.a > lambda.a(minred_wla)' & lambda.a > wla_700) = NaN;
  p.ap_sd(lambda.a > lambda.a(minred_wla)' & lambda.a > wla_700) = NaN;
  
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
    above_median_d600_ap450 = ratiod600_ap450 > fudge_list(i) * median(ratiod600_ap450, 'omitnan');
    ndel_spec(i) = sum(above_median_d600_ap450);
  end
  fudge_factor = fudge_list(find(abs(diff(ndel_spec)) == min(abs(diff(ndel_spec))) & ...
    ndel_spec(2:end) < 0.05 * max(ndel_spec), 1,  'first')); % threshold on first derivative of number of spectrum deleted
  above_median_d600_ap450 = ratiod600_ap450 > fudge_factor * median(ratiod600_ap450);
  fprintf('%.2f%% (%i) spectrum failed auto-QC step 2: sum(abs(d(ap)/d(lambda(600-650)))) / ap450nm\n', ...
    sum(above_median_d600_ap450) / size(p, 1) * 100, sum(above_median_d600_ap450))

%   delete bad spectrum
  bad = [bad; p(above_median_d600_ap450, :) ...
    table(repmat({'sum(abs(d(ap)/d(lambda(600-650)))) / ap_{450nm}'}, sum(above_median_d600_ap450), 1), ...
    'VariableNames', {'QC_failed'})];
  p(above_median_d600_ap450, :) = [];
end

% Auto QC when a positive first derivatives of ap over
% wavelentght between 460 and 640 nm is larger than 0.2 times ap at 450nm
% or second derivatives of ap over wavelentght between 460 and 640 nm is
% larger than 0.006
ap_450 = p.ap(:, find(lambda.a >= 450, 1,'first'));
d460_640 = diff(p.ap(:, lambda.a > 460 & lambda.a <= 640),[],2);
% delete bad spectrum
todelete = any(d460_640 > 0.4 * ap_450,2) | any(abs(diff(d460_640,[],2)) > 0.05, 2);
fprintf('%.2f%% (%i) spectrum failed auto-QC step 3: d(ap)/d(lambda460-640) > 0.4 * ap_{450nm} | abs(d"(ap)/d(lambda460-640)) > 0.05)\n', ...
  sum(todelete) / size(p, 1) * 100, sum(todelete))
bad = [bad; p(todelete, :) table(repmat({'d(ap)/d(lambda460-640) > 0.4 * ap_{450nm} | abs(d"(ap)/d(lambda460-640)) > 0.05)'}, ...
  sum(todelete), 1), 'VariableNames', {'QC_failed'})];
p(todelete, :) = [];

% Auto QC when spectrum contains 3 consecutive positive first derivatives of ap over
% wavelentght between 485 and 570 nm
d485_570 = diff(p.ap(:, lambda.a > 485 & lambda.a <= 570),[],2);
pos_d485_570 = d485_570 > 0;
todelete = false(size(pos_d485_570, 1), 1);
N = 3; % Required number of consecutive numbers following a first one
for i = 1:size(todelete,1)
  t = [false pos_d485_570(i,:) false];
  if any(find(diff(t)==-1)-find(diff(t)==1)>=N) % First t followed by >=N consecutive numbers
    todelete(i) = true;
  end
end
fprintf('%.2f%% (%i) spectrum failed auto-QC step 4: 3 consecutive d(ap)/d(lambda485-570) > 0\n', ...
  sum(todelete) / size(p, 1) * 100, sum(todelete))
bad = [bad; p(todelete, :) table(repmat({'3 consecutive d(ap)/d(lambda485-570) > 0'}, ...
  sum(todelete), 1), 'VariableNames', {'QC_failed'})];
bad = sortrows(bad, 'dt');
p(todelete, :) = [];

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

% Derive Chl (Line heigh at 676 compared to 650 and 715)
% 	Chlorophyll a (chl) is computed using the particulate absorption line height at 676 nm and the global relationship from Tara Ocean (Boss et al. 2013)
% REFERENCES:
% Emmanuel Boss, Marc Picheral, Thomas Leeuw, Alison Chase, Eric Karsenti, Gabriel Gorsky, Lisa Taylor, Wayne Slade, Josephine Ras, Herve Claustre, 2013.The characteristics of particulate absorption, scattering and attenuation coefficients in the surface ocean; Contribution of the Tara Oceans expedition, Methods in Oceanography.
ap_a = interp1(lambda.a, p.ap',[650 676 715],'linear')';
p(ap_a(:,1) > ap_a(:,2), :) = []; % deleted unrealistic spectrum
ap_a(ap_a(:,1) > ap_a(:,2), :) = []; % deleted unrealistic spectrum
p.ap676_lh = ap_a(:,2)-(39/65*ap_a(:,1)+26/65*ap_a(:,3));
p.chl_ap676lh = 157*p.ap676_lh.^1.22;
p.chl_ap676lh(real(p.chl_ap676lh) ~= p.chl_ap676lh) = NaN;

% 3.3 Derive Gamma (does not support NaN values) (Boss et al. 2001)
% REFERENCES:
% Boss, E., W.S. Pegau, W.D. Gardner, J.R.V. Zaneveld, A.H. Barnard., M.S. Twardowski, G.C. Chang, and T.D. Dickey, 2001. Spectral particulate attenuation and particle size distribution in the bottom boundary layer of a continental shelf. Journal of Geophysical Research, 106, 9509-9516.
% [~,p.gamma] = FitSpectra_HM2(lambda.ref(:,1:end-2),p.cp(:,1:end-2));
% Correct bu on March 5, 2018, FitSpectra_HM2 does not accept NaNs in cp
p(all(isnan(p.cp),2),:) = [];
sel = ~any(isnan(p.cp));
[~,p.gamma] = FitSpectra_HM2(lambda.c(:,sel), p.cp(:,sel));
fprintf('Done\n')

fprintf('Calculating chlorophyll from cp (H_alh) ... ')
% Chlorophyll absorption (alh) and phytoplankton size eigenvectors (P) inferred from cp
% Houskeeper, H.F., Draper, D., Kudela, R.M., Boss, E., 2020. Chlorophyll absorption and phytoplankton size information inferred from hyperspectral particulate beam attenuation. Appl. Opt. 59, 6765. https://doi.org/10.1364/AO.396832
% First compute hskpr P parameters (put link to github of hkpr do not include in github).
[p.Halh, Pr] = houskeeperetal2020(lambda.c, p.cp);
p.chl_Halh = 157*p.Halh.^1.22;
fprintf('Done\n')

fprintf('Estimating G50 and mphi (slope of PSD) ... ')
% Use fit from Haëntjens et al. 2021v22 to get median average
% cross-sectional area (G50) and slope of phytoplankton size distribution
% (in abundance) (mphi):
p.HH_G50 = modelG50.predictFcn(Pr');
p.HH_mphi = modelmphi.predictFcn(Pr');
fprintf('Done\n')

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
  
  % remove NaNs
  filt_avg(all(isnan(filt_avg.a), 2) | all(isnan(filt_avg.c), 2),:) = [];
  
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
  [g.ag, g.cg] = TemperatureAndSalinityDependence(g.ag, g.cg, lambda.ref);
%   [g.ag, g.cg] = ResidualTemperatureAndScatteringCorrection(g.ag, g.cg, lambda.ref);
  fprintf('Done\n')

  % Propagate error
  %   Note: Error is not propagated through Scattering & Residual temperature
  %         correction as required by SeaBASS
  g.ag_sd = sqrt(filt_avg.a_avg_sd + di_interp.a_avg_sd);
  g.cg_sd = sqrt(filt_avg.c_avg_sd + di_interp.c_avg_sd);
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
  
%   % QC with ag and cg spectrums (limited testing on the QC)
%   g(g.ag(:,1) < 0 & g.cg(:,end-3) < -0.005, :) = [];
else
  g = table();
end
end

function [a_corr, c_corr, a_slope, c_slope] = TemperatureAndSalinityDependence(a, c, wl)
% Note that a and c does not have to be on the same wvelength, however the
%   function would need to be edited for that
% Sullivan et al. 2006 values
psi_wl = [400;402;404;406;408;410;412;414;416;418;420;422;424;426;428;430;432;434;436;438;440;442;444;446;448;450;452;454;456;458;460;462;464;466;468;470;472;474;476;478;480;482;484;486;488;490;492;494;496;498;500;502;504;506;508;510;512;514;516;518;520;522;524;526;528;530;532;534;536;538;540;542;544;546;548;550;552;554;556;558;560;562;564;566;568;570;572;574;576;578;580;582;584;586;588;590;592;594;596;598;600;602;604;606;608;610;612;614;616;618;620;622;624;626;628;630;632;634;636;638;640;642;644;646;648;650;652;654;656;658;660;662;664;666;668;670;672;674;676;678;680;682;684;686;688;690;692;694;696;698;700;702;704;706;708;710;712;714;716;718;720;722;724;726;728;730;732;734;736;738;740;742;744;746;748;750]';
psiT = [0.0001;0.0001;0.0001;0.0001;0;0;0;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0002;0.0003;0.0003;0.0004;0.0005;0.0006;0.0006;0.0007;0.0008;0.0009;0.001;0.001;0.001;0.001;0.001;0.0009;0.0009;0.0008;0.0007;0.0006;0.0006;0.0005;0.0004;0.0003;0.0003;0.0002;0.0001;0.0001;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0002;0.0002;0.0002;0.0001;0.0001;0.0001;0;0;-0.0001;-0.0001;-0.0001;-0.0001;-0.0001;-0.0001;0;0;0.0001;0.0002;0.0003;0.0005;0.0007;0.0009;0.0013;0.0017;0.0021;0.0026;0.0032;0.0038;0.0045;0.0054;0.0063;0.0073;0.0083;0.0094;0.0104;0.0113;0.0121;0.0128;0.0133;0.0136;0.0136;0.0133;0.0129;0.0124;0.0116;0.0107]';
c_psiS = [-1.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-4.0e-05;-3.0e-05;-3.0e-05;-2.0e-05;-1.0e-05;0;1.0e-05;1.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;0;0;0;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-3.0e-05;-4.0e-05;-5.0e-05;-6.0e-05;-6.0e-05;-8.0e-05;-9.0e-05;-0.00010;-0.00011;-0.00013;-0.00014;-0.00016;-0.00017;-0.00018;-0.00019;-0.00020;-0.00021;-0.00022;-0.00022;-0.00023;-0.00023;-0.00023;-0.00024;-0.00024;-0.00024;-0.00024;-0.00022;-0.00021;-0.00017;-0.00012;-6.0e-05;2.0e-05;0.00012;0.00022;0.00031;0.00041;0.00049;0.00056;0.00062]';
a_psiS = [3.0e-05;3.0e-05;3.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;0;0;0;0;0;0;0;0;0;0;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;0;1.0e-05;2.0e-05;3.0e-05;3.0e-05;4.0e-05;5.0e-05;5.0e-05;6.0e-05;6.0e-05;6.0e-05;6.0e-05;6.0e-05;6.0e-05;6.0e-05;5.0e-05;5.0e-05;5.0e-05;5.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;3.0e-05;3.0e-05;3.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;1.0e-05;1.0e-05;1.0e-05;0;0;-1.0e-05;-2.0e-05;-3.0e-05;-4.0e-05;-6.0e-05;-7.0e-05;-8.0e-05;-9.0e-05;-0.00011;-0.00012;-0.00014;-0.00015;-0.00016;-0.00017;-0.00018;-0.00019;-0.00020;-0.00020;-0.00020;-0.00021;-0.00021;-0.00021;-0.00021;-0.00021;-0.00020;-0.00017;-0.00013;-8.0e-05;-1.0e-05;7.0e-05;0.00016;0.00026;0.00037;0.00046;0.00054;0.00061;0.00067]';
% Interpolate Sullivan values on the current ACS
psiT=interp1(psi_wl,psiT,wl,'spline'); % PCHIP or SPLINE -> better than linear
a_psiS=interp1(psi_wl,a_psiS,wl,'spline');
c_psiS=interp1(psi_wl,c_psiS,wl,'spline');
% Center psiS on 0 instead of +/- 0.001
a_psiS = a_psiS - median(a_psiS(wl <= 590));
c_psiS = c_psiS - median(c_psiS(wl <= 590));
% tmp = xlsread('Sullivan_etal_2006_instrumentspecific.xls');
% psi_wl = tmp(:,1);
% phiT=interp1(tmp(:,1),tmp(:,2),wl,'linear');
% c_phiS=interp1(tmp(:,1),tmp(:,4),wl,'linear');
% a_phiS=interp1(tmp(:,1),tmp(:,6),wl,'linear');


% Parameters of minization routine
opts = optimset('fminsearch');      
opts = optimset(opts,'MaxIter',20000000); 
opts = optimset(opts,'MaxFunEvals',20000); % 20000
opts = optimset(opts,'TolX',1e-8);
opts = optimset(opts,'TolFun',1e-8);

% Wavelength selection
iwl = wl >= 450;

% Init minimization parameters
deltaT = 10;
deltaS = 20;
amp = a(20);
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
  x = fminsearch(@costfun_TSD, [deltaT, deltaS, amp, slope], opts, a(k,iwl), psiT(iwl), a_psiS(iwl), wl(iwl));
  a_deltaT(k) = x(1); a_deltaS(k) = x(2); a_amp(k) = x(3); a_slope(k) = x(4);
  x = fminsearch(@costfun_TSD, [deltaT, deltaS, amp, slope], opts, c(k,iwl), psiT(iwl), c_psiS(iwl), wl(iwl));
  c_deltaT(k) = x(1); c_deltaS(k) = x(2); c_amp(k) = x(3); c_slope(k) = x(4);
  % Without Salinity forcing
%   x = fminsearch(@costfun_TSD, [deltaT, amp, slope], opts, a(k,iwl), psiT(iwl), a_psiS(iwl), wl(iwl));
%   a_deltaT(k) = x(1); a_amp(k) = x(2); a_slope(k) = x(3);
%   x = fminsearch(@costfun_TSD, [deltaT, amp, slope], opts, c(k,iwl), psiT(iwl), c_psiS(iwl), wl(iwl));
%   c_deltaT(k) = x(1); c_amp(k) = x(2); c_slope(k) = x(3);
end

% Apply correction
a_corr = a - a_deltaT.*psiT - deltaS.*a_psiS;
c_corr = c - c_deltaT.*psiT - deltaS.*c_psiS;

% Display the forcing parameters
% disp([a_deltaT, c_deltaT, a_deltaS, c_deltaS, a_slope, c_slope])
end

function cost = costfun_TSD(x, a, psiT, psiS, wl)
  % Force Temperature, Salinity, and Slope
  cost = sum((a - psiT.*x(1) - psiS.*x(2) - x(3).*exp(-x(4)*(wl-450))).^2);
  % Without Salinity forcing
%   cost = sum((a - psiT.*x(1) - x(2).*exp(-x(3)*(wl-450))).^2);
end

function [ap_corr, cp_corr] = ResidualTemperatureAndScatteringCorrection(ap, cp, wl)
% Function from Emmanuel Boss improved by Nils Haëntjens

% Sullivan et al. 2006 values
psi_wl = [400;402;404;406;408;410;412;414;416;418;420;422;424;426;428;430;432;434;436;438;440;442;444;446;448;450;452;454;456;458;460;462;464;466;468;470;472;474;476;478;480;482;484;486;488;490;492;494;496;498;500;502;504;506;508;510;512;514;516;518;520;522;524;526;528;530;532;534;536;538;540;542;544;546;548;550;552;554;556;558;560;562;564;566;568;570;572;574;576;578;580;582;584;586;588;590;592;594;596;598;600;602;604;606;608;610;612;614;616;618;620;622;624;626;628;630;632;634;636;638;640;642;644;646;648;650;652;654;656;658;660;662;664;666;668;670;672;674;676;678;680;682;684;686;688;690;692;694;696;698;700;702;704;706;708;710;712;714;716;718;720;722;724;726;728;730;732;734;736;738;740;742;744;746;748;750]';
psiT = [0.0001;0.0001;0.0001;0.0001;0;0;0;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0002;0.0003;0.0003;0.0004;0.0005;0.0006;0.0006;0.0007;0.0008;0.0009;0.001;0.001;0.001;0.001;0.001;0.0009;0.0009;0.0008;0.0007;0.0006;0.0006;0.0005;0.0004;0.0003;0.0003;0.0002;0.0001;0.0001;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0002;0.0002;0.0002;0.0001;0.0001;0.0001;0;0;-0.0001;-0.0001;-0.0001;-0.0001;-0.0001;-0.0001;0;0;0.0001;0.0002;0.0003;0.0005;0.0007;0.0009;0.0013;0.0017;0.0021;0.0026;0.0032;0.0038;0.0045;0.0054;0.0063;0.0073;0.0083;0.0094;0.0104;0.0113;0.0121;0.0128;0.0133;0.0136;0.0136;0.0133;0.0129;0.0124;0.0116;0.0107]';
psiT = interp1(psi_wl, psiT, wl);

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
% Find nearest wavelength to greater than 715 nm to use as reference for correction
iref = find(715 <= wl, 1,'first'); % 730
% If ACS spectrum does not go up to 730 nm take the closest wavelength to 730 nm
if isempty(iref); [~, iref] = max(wl); end % works as there is data in iNIR so lowest wavelength is 710

% Initialize output arrays
deltaT = NaN(size(ap,1),1);

% Init routine parameters
bp = cp - ap;

% Run minimization routine on good spectrum only
sel = find(all(isfinite(ap),2));
for k = sel'
  deltaT(k) = fminsearch(@costFun_RTSC, 0, opts, ap(k,:), bp(k,:), psiT, iNIR, iref);         
end
ap_corr = ap - psiT.*deltaT - ((ap(:,iref) - psiT(iref).*deltaT) ./ bp(:,iref)) .* bp; 
cp_corr = cp - psiT.*deltaT;
end

function cost = costFun_RTSC(deltaT, ap, bp, psiT, iNIR, iref)
    cost = sum(abs(ap(iNIR) - psiT(iNIR).*deltaT - ((ap(iref)-psiT(iref).*deltaT)./bp(iref)).*bp(iNIR)));
end

%%%%%%%%%%%%%%%%%
% Ohter Methods %
%%%%%%%%%%%%%%%%%
function [a_ts, c_ts] = TemperatureAndSalinityCorrection(a, c, a_wl, c_wl, delta_t, delta_s)
wl_psi = [400;402;404;406;408;410;412;414;416;418;420;422;424;426;428;430;432;434;436;438;440;442;444;446;448;450;452;454;456;458;460;462;464;466;468;470;472;474;476;478;480;482;484;486;488;490;492;494;496;498;500;502;504;506;508;510;512;514;516;518;520;522;524;526;528;530;532;534;536;538;540;542;544;546;548;550;552;554;556;558;560;562;564;566;568;570;572;574;576;578;580;582;584;586;588;590;592;594;596;598;600;602;604;606;608;610;612;614;616;618;620;622;624;626;628;630;632;634;636;638;640;642;644;646;648;650;652;654;656;658;660;662;664;666;668;670;672;674;676;678;680;682;684;686;688;690;692;694;696;698;700;702;704;706;708;710;712;714;716;718;720;722;724;726;728;730;732;734;736;738;740;742;744;746;748;750]';
psi_t = [0.0001;0.0001;0.0001;0.0001;0;0;0;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0002;0.0003;0.0003;0.0004;0.0005;0.0006;0.0006;0.0007;0.0008;0.0009;0.001;0.001;0.001;0.001;0.001;0.0009;0.0009;0.0008;0.0007;0.0006;0.0006;0.0005;0.0004;0.0003;0.0003;0.0002;0.0001;0.0001;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0002;0.0002;0.0002;0.0001;0.0001;0.0001;0;0;-0.0001;-0.0001;-0.0001;-0.0001;-0.0001;-0.0001;0;0;0.0001;0.0002;0.0003;0.0005;0.0007;0.0009;0.0013;0.0017;0.0021;0.0026;0.0032;0.0038;0.0045;0.0054;0.0063;0.0073;0.0083;0.0094;0.0104;0.0113;0.0121;0.0128;0.0133;0.0136;0.0136;0.0133;0.0129;0.0124;0.0116;0.0107]';
c_psi_s = [-1.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-4.0e-05;-3.0e-05;-3.0e-05;-2.0e-05;-1.0e-05;0;1.0e-05;1.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;0;0;0;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-3.0e-05;-4.0e-05;-5.0e-05;-6.0e-05;-6.0e-05;-8.0e-05;-9.0e-05;-0.00010;-0.00011;-0.00013;-0.00014;-0.00016;-0.00017;-0.00018;-0.00019;-0.00020;-0.00021;-0.00022;-0.00022;-0.00023;-0.00023;-0.00023;-0.00024;-0.00024;-0.00024;-0.00024;-0.00022;-0.00021;-0.00017;-0.00012;-6.0e-05;2.0e-05;0.00012;0.00022;0.00031;0.00041;0.00049;0.00056;0.00062]';
a_psi_s = [3.0e-05;3.0e-05;3.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;0;0;0;0;0;0;0;0;0;0;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;0;1.0e-05;2.0e-05;3.0e-05;3.0e-05;4.0e-05;5.0e-05;5.0e-05;6.0e-05;6.0e-05;6.0e-05;6.0e-05;6.0e-05;6.0e-05;6.0e-05;5.0e-05;5.0e-05;5.0e-05;5.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;3.0e-05;3.0e-05;3.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;1.0e-05;1.0e-05;1.0e-05;0;0;-1.0e-05;-2.0e-05;-3.0e-05;-4.0e-05;-6.0e-05;-7.0e-05;-8.0e-05;-9.0e-05;-0.00011;-0.00012;-0.00014;-0.00015;-0.00016;-0.00017;-0.00018;-0.00019;-0.00020;-0.00020;-0.00020;-0.00021;-0.00021;-0.00021;-0.00021;-0.00021;-0.00020;-0.00017;-0.00013;-8.0e-05;-1.0e-05;7.0e-05;0.00016;0.00026;0.00037;0.00046;0.00054;0.00061;0.00067]';
%sigma_psi_t = [0.0002;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0;0;0;0;0.0001;0.0001;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0;0;0;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0002;0.0002;0.0003;0.0003;0.0004;0.0004;0.0004;0.0004;0.0005;0.0005;0.0006;0.0006;0.0007;0.0007;0.0007;0.0006;0.0005;0.0004;0.0003;0.0003;0.0004;0.0005;0.0006;0.0007;0.0008;0.0009;];
%c_sigma_psi_s = [4e-005;4e-005;4e-005;4e-005;4e-005;4e-005;4e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;0;-1e-005;-2e-005;-3e-005;-4e-005;-6e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;2e-005;2e-005;2e-005;3e-005;3e-005;3e-005;];
%a_sigma_psi_s = [3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;NaN;NaN;NaN;NaN;NaN;NaN;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;2e-005;2e-005;2e-005;3e-005;];
%values from Sullivan et al 2006 (Applied Optics)

%interpolate literature psiT values to acs wavelengths
a_psi_t = interp1(wl_psi,psi_t,a_wl,'linear','extrap');
c_psi_t = interp1(wl_psi,psi_t,c_wl,'linear','extrap');
a_psi_s = interp1(wl_psi,a_psi_s,a_wl,'linear','extrap');
c_psi_s = interp1(wl_psi,c_psi_s,c_wl,'linear','extrap');
%a_sigma_psi_t = interp1(wl_psi,a_sigma_psi_s,wl_a,'linear','extrap');
%c_sigma_psi_t = interp1(wl_psi,c_sigma_psi_s,wl_c,'linear','extrap');

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
  for j = 1:max(size(lambda.(i{:}(1))))
      for jkl = 1:max(size(wavelength))
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

function agaus = GaussDecomp(p, lambda, compute_ad_aphi)
% Gaussian decomposition
% Ron Zaneveld, WET Labs, Inc., 2005
% Ali Chase, University of Maine, 2014
% Adaptation to InLineAnalysis: Guillaume Bourdin, March 2021
%
% Reference:
% Chase, A., et al., Decomposition of in situ particulate absorption
% spectra. Methods in Oceanography (2014), http://dx.doi.org/10.1016/j.mio.2014.02.022

%% delete row and columns full of nans
lambda(all(isnan(p.ap),1)) = [];
%   data(all(isnan(p.ap),2),:) = [];
p.ap_sd(:, all(isnan(p.ap),1)) = [];
p.ap(:, all(isnan(p.ap),1)) = [];

% Peak center values ("peak_loc") determined using a interative
% approach that allows the location to vary (uses the matlab
% function LSQNONLIN), and are rounded to nearest integer.
% Sigma values ("lsqsig") are determined similarly. FWHM = sigma*2.355
peak_loc = [406,434,453,470,492,523,550,584,617,638,660,675];
lsqsig = [16,12,12,13,16,14,14,16,13,11,11,10];
onenm = 400:1:720;

fprintf('Gaussian decomposition ... ')

% extrapolate NaN tails with constant nearest non-NaN value
ap_filled = fillmissing(p.ap, 'nearest', 2);
ap_sd_filled = fillmissing(p.ap_sd, 'nearest', 2);

% interpolate the un-smoothed ap spectra to one nm resolution
acorr2onenm = interp1(lambda, ap_filled', onenm, 'spline');

% visProd3D(onenm, p.dt, ...
%   acorr2onenm', false, 'Wavelength', false, 72);
% zlabel('a_p (m^{-1})'); xlabel('{\lambda} (nm)'); ylabel('Time');


% define the matrix of component Gaussian functions using the peaks
% and widths (sigma) above
coef2 = exp(-0.5 .* (((onenm .* ones(size(peak_loc,2),1))' - peak_loc .* ...
    ones(size(onenm,2),1)) ./ lsqsig) .^ 2);

% define a function for non-algal particles and concatenate this to the Gaussian matrix
coef2nap = exp(-0.01 * (onenm - 400));
coef2 = [coef2nap', coef2];

% normalize both the component functions and the measured ap
% spectrum by the uncertainty (std dev) in the ap spectra
ap_sd_int = interp1(lambda, ap_sd_filled', onenm, 'linear', 'extrap');
acorr2onenm_new = (acorr2onenm ./ ap_sd_int);

amps = NaN(size(p,1), size(coef2,2));
sumspec_temp = NaN(size(coef2,1), size(p,1));
% compspec_temp = NaN(size(coef2,1), size(coef2,2), size(p,1));


% visProd3D(peak_loc, p.dt, ...
%   popo(:, 2:end), false, 'Wavelength', false, 73);
% zlabel('a_p (m^{-1})'); xlabel('{\lambda} (nm)'); ylabel('Time');
% 
% visProd3D(peak_loc, p.dt, ...
%   amps(:, 2:end), false, 'Wavelength', false, 74);
% zlabel('a_p (m^{-1})'); xlabel('{\lambda} (nm)'); ylabel('Time');



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

agaus = array2table([amps uncertainty], 'VariableNames', ...
    [{'ad_model400'} cellfun(@(x) ['agaus' x], cellstr(num2str(peak_loc'))', 'un', 0) {'agaus-mae'}]);

if compute_ad_aphi
  % Compute a non-algal and a phytoplankton from Zheng, G., and D. Stramski (2013) model
  % Reference: 
  % Zheng, G., and D. Stramski (2013), A model based on stacked-constraints approach for partitioning the light absorption coefficient of seawater 
  % into phytoplankton and non-phytoplankton components, J. Geophys. Res. Oceans, 118, 2155?2174, doi:10.1002/jgrc.20115.
  %
  % Zheng, G., and D. Stramski (2013), A model for partitioning the light absorption coefficient of suspended marine particles into phytoplankton 
  % and nonalgal components, J. Geophys. Res. Oceans, 118, 2977?2991, doi:10.1002/jgrc.20206.
  wl_ZS13 = [400 412 420 430 443 450 467 490 500 510 550 555 630 650 670 700];
  qwl = unique([wl_ZS13(:)' 442 676]);
  ap_ZS13 = interp1(lambda, ap_filled', qwl, 'linear', 'extrap'); % need to extrapolate for 400

%   % vectorized partition_ap (DOESN'T WORK)
%   agaus.ad_ZS13 = NaN(size(agaus, 1), size(qwl, 2));
%   agaus.aphi_ZS13 = NaN(size(agaus, 1), size(qwl, 2));
%   [agaus.ad_ZS13, agaus.aphi_ZS13] = partition_ap_vec(ap_ZS13, qwl', 50);  % Must be in colum direction
  
  % Run partition_ap
  ad_ZS13 = NaN(size(agaus, 1), size(qwl, 2));
  aphi_ZS13 = NaN(size(agaus, 1), size(qwl, 2));
%   parfor i = 1:300 %size(ap_ZS13, 2)
  for i = progress(1:300) %size(ap_ZS13, 2)
    [ad_ZS13(i,:), aphi_ZS13(i,:)] = partition_ap(ap_ZS13(:, i), qwl', 50);
    [agaus.ad_ZS13(i, :), agaus.aphi_ZS13(i, :)] = partition_ap(ap_ZS13(:, i), qwl', 50);  % Must be in colum direction
  end
%   figure; hold on
%   plot(qwl, ad)
%   plot(qwl, aphi_ZS13)
  agaus.Properties.VariableUnits = repmat({'1/m'}, 1, size(agaus, 2));
end
fprintf('Done\n')
end



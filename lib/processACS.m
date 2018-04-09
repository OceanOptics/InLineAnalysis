function [p, g] = processACS(lambda, tot, filt, di, cdom, fth, fth_constants)
% NOTE: wavelength of c are interpolated to wavelength of a
%% ap & cp
if nargin > 5
  n_wv = size(lambda.ref,2);
  % Require both CDOM & Switch position
  % Parameters for CDOM signal (in units of cdom.fdom, by default it's counts);
  param_extrap_threshold = 4;  % counts
  param_min_variability = 1;   % counts
  if nargin < 7
    % Assume most recent FlowControl software
    SWITCH_FILTERED = 1;
    SWITCH_TOTAL = 0;
  else
    SWITCH_FILTERED = fth_constants.SWITCH_FILTERED;
    SWITCH_TOTAL = fth_constants.SWITCH_TOTAL;
  end
  % Find switch events from total to filtered
  sel_start = find(fth.swt(1:end-1) == SWITCH_TOTAL & fth.swt(2:end) == SWITCH_FILTERED);
  % Find switch events from filtered to total
  sel_end = find(fth.swt(1:end-1) == SWITCH_FILTERED & fth.swt(2:end) == SWITCH_TOTAL);
  % Verify selections of filtered period
  if sel_start(1) > sel_end(1); sel_end(1) = []; end
  if sel_start(end) > sel_end(end); sel_end(end+1) = size(fth.swt,1); end
  if size(sel_start,1) ~= size(sel_end,1); error('Inconsistent fth data'); end
  % Make filtered period averages
  filt_avg = table((fth.dt(sel_start) + fth.dt(sel_end)) ./ 2, 'VariableNames', {'dt'});
  filt_avg.cdom = NaN(size(filt_avg,1),1);
  filt_avg.a = NaN(size(filt_avg,1),n_wv);
  filt_avg.c = NaN(size(filt_avg,1),n_wv);
  for i=1:size(sel_start, 1)
    sel_cdom = fth.dt(sel_start(i)) <= cdom.dt & cdom.dt <= fth.dt(sel_end(i));
    filt_avg.cdom(i) = median(cdom.fdom(sel_cdom));
    sel_filt = fth.dt(sel_start(i)) <= filt.dt & filt.dt <= fth.dt(sel_end(i));
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
else
  % Interpolate filtered on total linearly
  filt_interp = table(tot.dt, 'VariableNames', {'dt'});
  filt_interp.a = interp1(filt.dt, filt.a, filt_interp.dt);%, 'linear', 'extrap');
  filt_interp.c = interp1(filt.dt, filt.c, filt_interp.dt);%, 'linear', 'extrap');
  filt_interp.a_avg_sd = interp1(filt.dt, filt.a_avg_sd, filt_interp.dt);%, 'linear', 'extrap');
  filt_interp.c_avg_sd = interp1(filt.dt, filt.c_avg_sd, filt_interp.dt);%, 'linear', 'extrap');
end

% Particulate = Total - FSW
p = table(tot.dt, 'VariableNames', {'dt'});
p.ap = tot.a - filt_interp.a;
p.cp = tot.c - filt_interp.c;

% Interpolate wavelengths
% p.cp = cell2mat(arrayfun(@(i) interp1(c_wl, p.cp(i,:), a_wl, 'linear', 'extrap'), 1:size(p,1), 'UniformOutput', false)');
p.ap = interp1(lambda.a', p.ap', lambda.ref', 'linear', 'extrap')';
p.cp = interp1(lambda.c', p.cp', lambda.ref', 'linear', 'extrap')';

% ap Scattering & Residual temperature correction
% cp Residual correction (for efficiency use the one computed from ap as it should be the same)
[p.ap, p.cp] = ResidualTemperatureAndScatteringCorrection(p.ap, p.cp, lambda.ref);

% Propagate error
%   Note: Error is not propagated through Scattering & Residual temperature
%         correction as required by SeaBASS
p.ap_sd = sqrt(tot.a_avg_sd + filt_interp.a_avg_sd);
p.cp_sd = sqrt(tot.c_avg_sd + filt_interp.c_avg_sd);
p.ap_n = tot.a_avg_n;
p.cp_n = tot.c_avg_n;

% QC using ap spectrum
% TODO replace the transpose by std(.., 0, 1 or 2); line 149 for speed
sel_bad = any(p.ap(:,lambda.ref < 430) < 0,2)...
          | any(p.ap(:,:) < -0.0015,2)...
          | std(p.ap(:,lambda.ref < 430)')' > 6 * 10^-3;
p(sel_bad,:) = [];

% Derive standard products from ap and cp
% Derive POC (Specific to region)
cp660 = interp1(lambda.ref,p.cp',660,'linear')';
p.poc = cp660.*380;
% Derive Chl (Line heigh at 676 compared to 650 and 715)
ap_a=interp1(lambda.ref,p.ap',[650 676 715],'linear')';
line_height = (ap_a(:,2)-(39/65*ap_a(:,1)+26/65*ap_a(:,3)));
p.chl=157*line_height.^1.22;
p.chl(real(p.chl) ~= p.chl) = NaN;
% 3.3 Derive Gamma (does not support NaN values)
[~,p.gamma] = FitSpectra_HM2(lambda.ref(:,1:end-2),p.cp(:,1:end-2));

%% ag & cg
if nargin > 3 && nargout > 1
  % Interpolate filtered on Total
  di_interp = table(filt.dt, 'VariableNames', {'dt'});
  di_interp.a = interp1(di.dt, di.a, di_interp.dt, 'linear', 'extrap');
  di_interp.c = interp1(di.dt, di.c, di_interp.dt, 'linear', 'extrap');
  di_interp.a_avg_sd = interp1(di.dt, di.a_avg_sd, di_interp.dt, 'linear', 'extrap');
  di_interp.c_avg_sd = interp1(di.dt, di.c_avg_sd, di_interp.dt, 'linear', 'extrap');

  % Dissolved = Filtered - DI
  g = table(filt.dt, 'VariableNames', {'dt'});
  g.ag = filt.a - di_interp.a;
  g.cg = filt.c - di_interp.c;
  % Interpolate wavelength of c on a 
  % g.cg = cell2mat(arrayfun(@(i) interp1(c_wl, g.cg(i,:), a_wl, 'linear', 'extrap'), 1:size(g,1), 'UniformOutput', false)');
  g.ag = interp1(lambda.a', g.ag', lambda.ref', 'linear', 'extrap')';
  g.cg = interp1(lambda.c', g.cg', lambda.ref', 'linear', 'extrap')';
  % Temperature & Salinity Correction (No Scattering correction needed)
  [g.ag, g.cg] = TemperatureAndSalinityDependence(g.ag, g.cg, lambda.ref);
%   [g.ag, g.cg] = ResidualTemperatureAndScatteringCorrection(g.ag, g.cg, lambda.ref);

  % Propagate error
  %   Note: Error is not propagated through Scattering & Residual temperature
  %         correction as required by SeaBASS
  g.ag_sd = sqrt(filt.a_avg_sd + di_interp.a_avg_sd);
  g.cg_sd = sqrt(filt.c_avg_sd + di_interp.c_avg_sd);
  g.ag_n = filt.a_avg_n;
  g.cg_n = filt.c_avg_n;
  
  % QC with ag and cg spectrums (limited testing on the QC)
  g(g.ag(:,1) < 0 & g.cg(:,end-3) < -0.005, :) = [];
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
opts = optimset(opts,'MaxFunEvals',20000);
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
iref = find(730 <= wl, 1,'first');
% iref = find(720 <= wl, 1,'first'); % TODO: Add cfg variable as instrument dependent

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
c = interp1(c_wl, c, a_wl, 'linear', 'extrap');

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
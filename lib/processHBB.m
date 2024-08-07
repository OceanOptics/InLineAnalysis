function [p, g, FiltStat] = processHBB(param, tot, filt_qc, filt_raw, filt_bad, di, tsg, di_method, filt_method, fth, fth_constants, days2run)
% Note DI is not interpolated as it's assume to be stable in time
% BB3 parameters is a structure
%   param.lambda <1xM double> wavelength (nm)
%   param.theta <1x1 double> scattering angle (degree)
%   param.dark  <1xM double> dark (DEPRECATED)
%   param.slope <1xM double> slope (DEPRECATED)

% make sure lambda are in right dimension
param.lambda = param.lambda(:)';

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
  if islogical(fth_temp.swt); fth_temp.swt = double(fth_temp.swt); end

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
  fth_interp.swt = interp1(fth_temp.dt, fth_temp.swt, fth_interp.dt, 'previous');%, 'linear', 'extrap');
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
  filt_avg.beta = NaN(size(filt_avg,1), size(param.lambda, 2));
  filt_avg.beta_avg_sd = NaN(size(filt_avg,1), size(param.lambda, 2));
  filt_avg.beta_avg_n = NaN(size(filt_avg,1), 1);
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
          else
            perc25 = foo.beta > prctile(foo.beta, 25, 1);
            foo.beta_avg_sd(perc25) = NaN;
            foo.beta(perc25) = NaN;
            % compute average of all values smaller than 25th percentile for each filter event
            filt_avg.dt(i) = mean(foo.dt(any(~perc25, 2)), 'omitnan');
            filt_avg.beta(i,:) = mean(foo.beta, 1, 'omitnan');
            filt_avg.beta_avg_sd(i,:) = mean(foo.beta_avg_sd, 1, 'omitnan');
            filt_avg.beta_avg_n(i) = sum(foo.beta_avg_n(any(~isnan(foo.beta), 2)), 'omitnan');
          end
        end
      end
    case 'exponential_fit'
      % Based on method in: Dall’Olmo, G., Westberry, T.K., Behrenfeld, M.J., Boss, 
      %       E., Slade, W.H., 2009. Direct contribution of phytoplankton-sized particles 
      %       to optical backscattering in the open ocean. Biogeosciences Discuss 6, 291–340. 
      %       https://doi.org/10.5194/bgd-6-291-2009
      fprintf('Fitting exponential to filter events ... ')
      filt_avg.dt = (fth_interp.dt(sel_start) + fth_interp.dt(sel_end)) ./ 2;
      [filt_avg, FiltStat] = FiltExpFit('beta', filt_avg, filt_raw, filt_bad, fth_interp.dt(sel_start), fth_interp.dt(sel_end));
      fprintf('Done\n')      
      % run 25 percentile method on failed exponential fits
      for i=1:size(sel_start, 1)
        sel_filt = fth_interp.dt(sel_start(i)) <= filt_qc.dt & filt_qc.dt <= fth_interp.dt(sel_end(i));
        if sum(sel_filt) > 0
          foo = filt_qc(sel_filt,:);
          if sum(sel_filt) == 1
            filt_avg.dt(i) = foo.dt;
            filt_avg.beta(i,~FiltStat.exitflag(i,:)) = foo.beta(:,~FiltStat.exitflag(i,:));
          else
            perc25 = foo.beta > prctile(foo.beta, 25, 1);
            foo.beta_avg_sd(perc25) = NaN;
            foo.beta(perc25) = NaN;
            % compute average of all values smaller than 25th percentile for each filter event
            filt_avg.dt(i) = mean(foo.dt(any(~perc25, 2)), 'omitnan');
            filt_avg.beta(i,~FiltStat.exitflag(i,:)) = mean(foo.beta(:,~FiltStat.exitflag(i,:)), 1, 'omitnan');
            filt_avg.beta_avg_sd(i,~FiltStat.exitflag(i,:)) = mean(foo.beta_avg_sd(:,~FiltStat.exitflag(i,:)), 1, 'omitnan');
            filt_avg.beta_avg_n(i) = sum(foo.beta_avg_n(any(~isnan(foo.beta(:,~FiltStat.exitflag(i,:))), 2)), 'omitnan');
          end
        else
          filt_avg.beta(i,:) = NaN;
        end
      end
    otherwise
      error('filter event method "filt_method" not supported')
  end
  filt_avg(all(isnan(filt_avg.beta), 2), :) = [];
else
  filt_avg = filt;
end

% Interpolate filtered on total linearly
filt_interp = table(tot.dt, 'VariableNames', {'dt'});
filt_interp.beta = interp1(filt_avg.dt, filt_avg.beta, filt_interp.dt);%, 'linear', 'extrap');
filt_interp.beta_avg_sd = interp1(filt_avg.dt, filt_avg.beta_avg_sd, filt_interp.dt);%, 'linear', 'extrap');

% id only day to run in all tables to plot
filt_interp_id = filt_interp.dt >= min(days2run) & filt_interp.dt < max(days2run)+1;
tot_id = tot.dt >= min(days2run) & tot.dt < max(days2run)+1;
filt_avg_id = filt_avg.dt >= min(days2run) & filt_avg.dt < max(days2run)+1;

% plot
if exist('visFlag', 'file') && exist('fth', 'var')
  fh = visFlag([], filt_interp(filt_interp_id, :), tot(tot_id, :), [], filt_avg(filt_avg_id, :), [], 'beta', round(size(tot.beta, 2)/2), ...
    [], fth_temp, fth.view.spd_variable);
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
p.betap = tot.beta - filt_interp.beta;

% Interpolate X_p with values from Sullivan et al. 2013
theta_ref = 90:10:170;
X_p_ref = [0.709 0.884 1.024 1.121 1.174 1.171 1.138 1.039 0.931]; % Weighted averages
% X_p_ref = [0.684 0.858 1.000 1.097 1.153 1.167 1.156 1.131 1.093];
X_p = interp1(theta_ref, X_p_ref, param.theta, 'spline');

% Compute backscattering
p.bbp = 2 * pi * X_p .* p.betap;

% Propagate error
%   Note: Error is not propagated through Scattering & Residual temperature
%         correction as required by SeaBASS
% p.betap_sd = param.slope .* sqrt(tot.beta_avg_sd.^2 + filt_interp.beta_avg_sd.^2);
p.betap_sd = sqrt(tot.beta_avg_sd.^2 + filt_interp.beta_avg_sd.^2);
p.betap_n = tot.beta_avg_n;

% Derive Gamma_bbp (does not support NaN values)
% Correct bu on March 5, 2018, FitSpectra_HM2 does not accept NaNs in spectra
p(all(isnan(p.bbp),2),:)=[];
sel_lambda = param.lambda < 510 | param.lambda > 650;
sel = ~isnan(p.bbp);
sel = sum(sel, 2) > 5;
fprintf('Computing Gamma_bbp ... ')
[~, p.gamma_bbp(sel)] = FitSpectra_HM2(param.lambda(sel_lambda), p.bbp(sel, sel_lambda));
p.gamma_bbp(p.gamma_bbp < -0.5) = NaN;
fprintf('Done\n')

% Estimate POC and Cphyto from bbp
[p.poc, ~, ~, p.cphyto, ~, ~,] = estimatePOC_Cphyto(p.bbp, param.lambda, 'soccom');

% QC
p(any(p.bbp < 0,2),:) = [];

if nargout > 1 && any(~isempty(di) | strcmp(di_method,'SW_scattering')) && ~isempty(tsg)
  if ~isempty(tsg.prod.a)
    tsg_data = tsg.prod.a;
  elseif ~isempty(tsg.qc.tsw)
    tsg_data = tsg.qc.tsw;
  else
    error('No TSG qc or prod data loaded')
  end
  if ~any(tsg_data.dt >= min([tot.dt; filt_qc.dt]) & tsg_data.dt <= max([tot.dt; filt_qc.dt]))
    fprintf('Warning: TSG dates do not correspond to ACS dates: salinity correction not performed\n')
  else
    % round time stamp and remove time duplicates
    tsg_data = round_timestamp(tsg_data);
  end

  replace_consecutive_nan = 3*60; % 3h
  filt_avg = merge_timeseries(filt_avg, tsg_data, {tsg.temperature_variable, 's'}, '', replace_consecutive_nan);
  switch di_method
    case 'interpolate'
      % Interpolate DI on Filtered
      %     + recommende if sensor drift with time
      di_pp = table(filt_qc.dt, 'VariableNames', {'dt'});
      di_pp.beta = interp1(di.dt, di.beta, di_pp.dt);
      di_pp.beta_avg_sd = interp1(di.dt, di.beta_avg_sd, di_pp.dt);
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
    case 'SW_scattering'
    otherwise
      error('Method not supported.');
  end
  % Get beta salt from Zhang et al. 2009
%   beta_s = NaN(size(filt_qc.beta));
%   for i = 1:max(size(param.lambda))
%     for j = 1:size(t,1)
%       beta_s(j,i) = betasw_ZHH2009(param.lambda(i), t(j), param.theta, s(j)) - betasw_ZHH2009(param.lambda(i), t(j), param.theta, 0);
%     end
%   end

  switch di_method
    case {'interpolate', 'constant'}
      % Get beta salt from Zhang et al. 2009
      beta_s = NaN(size(filt_qc.beta));
      for j = 1:size(filt_avg,1)
        beta_s(j, :) = betasw_ZHH2009(param.lambda, filt_avg.(tsg.temperature_variable)(j), param.theta, filt_avg.s(j)) - ...
          betasw_ZHH2009(param.lambda, filt_avg.(tsg.temperature_variable)(j), param.theta, 0);
      end
      % Compute beta dissolved
      g = table(filt_qc.dt, 'VariableNames', {'dt'});
      g.betag = (filt_qc.beta - di_pp.beta) - beta_s ;
      % Propagate error
      %   Note: Error is not propagated through Scattering & Residual temperature
      %         correction as required by SeaBASS
      % g.betag_sd = param.slope .* sqrt(filt_qc.beta_avg_sd.^2 + di_pp.beta_avg_sd.^2);
      g.betag_sd = sqrt(filt_qc.beta_avg_sd.^2 + di_pp.beta_avg_sd.^2);
    case 'SW_scattering'
      % Get beta salt from Zhang et al. 2009
      beta_sw = NaN(size(filt_qc.beta));
      for j = 1:size(filt_avg,1)
        beta_sw(j, :) = betasw_ZHH2009(param.lambda, filt_avg.(tsg.temperature_variable)(j), param.theta, filt_avg.s(j));
      end
      % Compute beta dissolved
      g = table(filt_qc.dt, 'VariableNames', {'dt'});
      % g.betag = param.slope .* filt_qc.beta - beta_s ;
      g.betag = filt_qc.beta - beta_sw ;
      % Propagate error
      %   Note: Error is not propagated through Scattering & Residual temperature
      %         correction as required by SeaBASS
      % g.betag_sd = param.slope .* sqrt(filt_qc.beta_avg_sd.^2 + di_pp.beta_avg_sd.^2);
      g.betag_sd = filt_qc.beta_avg_sd;
    otherwise
      error('Method not supported.');
  end
  g.betag_n = filt_qc.beta_avg_n;
  
  % Compute bbg and gamma_bbg
  g.bbg = 2 * pi * X_p .* g.betag;
  
  % Derive Gamma_bbp (does not support NaN values)
  % Correct bu on March 5, 2018, FitSpectra_HM2 does not accept NaNs in spectra
  fprintf('Computing beta filtered slope... ')
  g(all(isnan(g.bbg),2),:)=[]; % delete line full of NaNs
  sel = ~any(isnan(g.bbg)); % select wavelenght with no NaN
  sel(param.lambda < 430 | param.lambda > 650) = false;
  [~, g.gamma_bbg] = FitSpectra_HM2(param.lambda(sel), g.bbg(:, sel));  
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

if size(lambda,2) < size(lambda,1)
  lambda = lambda';
end

% Check size of input/content of input
if size(bbp,2) ~= size(lambda,2)
  error('size(bbp,2) ~= numbe of wavelength input');
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


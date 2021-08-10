function [p, g, FiltStat] = processBB3(param, tot, filt_qc, filt_raw, filt_bad, di, tsg, di_method, filt_method, fth, fth_constants)
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
  
  % sort filt_qc data
  filt_qc = sortrows(filt_qc, 'dt');
  fth.swt(fth.swt > 0) = 1;
  % delete duplicats
  [~, L, ~] = unique(fth.dt,'first');
  indexToDump = not(ismember(1:numel(fth.dt),L));
  fth(indexToDump, :) = [];
  
  % merge and sort filt_raw filt_bad data
  if any(~isempty(filt_raw) | ~isempty(filt_bad))
    filt_raw_merged = sortrows([filt_raw; filt_bad], 'dt');
  end
  
  % interpolate fth.swt onto binned data to fill missing flow data
  fth_interp = table([tot.dt; fth.dt; filt_qc.dt], 'VariableNames', {'dt'});
  % delete duplicats
  [~, L, ~] = unique(fth_interp.dt,'first');
  indexToDump = not(ismember(1:numel(fth_interp.dt),L));
  fth_interp(indexToDump, :) = [];
  % sort dates
  [~,b] = sort(fth_interp.dt);
  fth_interp.dt = fth_interp.dt(b,:);
  fth_interp.swt = interp1(fth.dt, fth.swt, fth_interp.dt, 'previous');
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
  filt_avg = table((fth_interp.dt(sel_start) + fth_interp.dt(sel_end)) ./ 2, 'VariableNames', {'dt'});
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
            filt_avg.beta(i,:) = foo.beta;
            filt_avg.beta_avg_sd(i,:) = foo.beta_avg_sd;
            filt_avg.beta_avg_n(i,:) = foo.beta_avg_n;
          else
            foo.beta_avg_sd(foo.beta > prctile(foo.beta, 25, 1)) = NaN;
            foo.beta(foo.beta > prctile(foo.beta, 25, 1)) = NaN;
            % compute average of all values smaller than 25th percentile for each filter event
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
      fprintf('Fitting exponential to filter events ... \n')
      [filt_avg, FiltStat] = FiltExpFit(filt_avg, filt_raw_merged, fth_interp.dt(sel_start), fth_interp.dt(sel_end));
      fprintf('Done\n')
      % run 25 percentile method on failed exponential fits
      for i=1:size(sel_start, 1)
        sel_filt = fth_interp.dt(sel_start(i)) <= filt_qc.dt & filt_qc.dt <= fth_interp.dt(sel_end(i));
        if sum(sel_filt) > 0
          if any(~FiltStat.exitflag(i,:))
            foo = filt_qc(sel_filt,:);
            if sum(sel_filt) == 1
              filt_avg.beta(i,~FiltStat.exitflag(i,:)) = foo.beta(:,~FiltStat.exitflag(i,:));
            else
              foo.beta_avg_sd(foo.beta > prctile(foo.beta, 25, 1)) = NaN;
              foo.beta(foo.beta > prctile(foo.beta, 25, 1)) = NaN;
              % compute average of all values smaller than 25th percentile for each filter event
              filt_avg.beta(i,~FiltStat.exitflag(i,:)) = mean(foo.beta(:,~FiltStat.exitflag(i,:)), 1, 'omitnan');
              filt_avg.beta_avg_sd(i,~FiltStat.exitflag(i,:)) = mean(foo.beta_avg_sd(:,~FiltStat.exitflag(i,:)), 1, 'omitnan');
              filt_avg.beta_avg_n(i) = sum(foo.beta_avg_n(any(~isnan(foo.beta(:,~FiltStat.exitflag(i,:))), 2)), 'omitnan');
            end
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
  filt_avg = filt_qc;
end

% Interpolate filtered on total linearly
filt_interp = table(tot.dt, 'VariableNames', {'dt'});
filt_interp.beta = interp1(filt_avg.dt, filt_avg.beta, filt_interp.dt);%, 'linear', 'extrap');
filt_interp.beta_avg_sd = interp1(filt_avg.dt, filt_avg.beta_avg_sd, filt_interp.dt);%, 'linear', 'extrap');

if exist('visFlag', 'file') && exist('fth', 'var')
  fh = visFlag([], filt_interp, tot, [], filt_avg, [], 'beta', round(size(tot.beta, 2)/2), [], fth);
  title('Check filter event interpolation, press q to continue', 'FontSize', 14)
  legend('Filtered interpolated', 'Total', 'Filtered median', 'Flow rate',...
    'AutoUpdate','off', 'FontSize', 12)
  guiSelectOnTimeSeries(fh);
elseif exist('visFlag', 'file')
  fh = visFlag([], filt_interp, tot, [], filt_avg, [], 'beta', round(size(tot.beta, 2)/2), [], []);
  title('Check filter event interpolation, press q to continue', 'FontSize', 14)
  legend('Filtered interpolated', 'Total', 'Filtered median', 'Flow rate',...
    'AutoUpdate','off', 'FontSize', 12)
  guiSelectOnTimeSeries(fh);
end

% Compute beta particulate
p = table(tot.dt, 'VariableNames', {'dt'});
p.betap = tot.beta - filt_interp.beta;

% Calibrate beta_p (counts to scientific units)
p.betap = param.slope .* p.betap; % Dark independent

% Interpolate X_p with values from Sullivan et al. 2013
% theta_ref = 90:10:170;
% X_p_ref = [0.684 0.858 1.000 1.097 1.153 1.167 1.156 1.131 1.093];
% X_p = interp1(theta_ref, X_p_ref, theta, 'spline');
X_p = 1.076; % Specific to ECO-BB3 with a scattering angle of 124

% Compute backscattering
p.bbp = 2 * pi * X_p .* p.betap;

% Propagate error
%   Note: Error is not propagated through Scattering & Residual temperature
%         correction as required by SeaBASS
p.betap_sd = param.slope .* sqrt(tot.beta_avg_sd + filt_interp.beta_avg_sd);
p.betap_n = tot.beta_avg_n;

% Derive Gamma_bbp (does not support NaN values)
% Correct bu on March 5, 2018, FitSpectra_HM2 does not accept NaNs in spectra
p(all(isnan(p.bbp),2),:)=[];
fprintf('Computing Gamma_bbp ... ')
[~, p.gamma_bbp] = FitSpectra_HM2(param.lambda, p.bbp);
fprintf('Done\n')

% Estimate POC and Cphyto from bbp
[p.poc, ~, ~, p.cphyto, ~, ~,] = estimatePOC_Cphyto(p.bbp, [470 532 650], 'soccom');

% remove negative values
p(any(p.bbp < 0,2),:) = [];

if nargout > 1 && any(~isempty(di) | strcmp(di_method,'SW_scattering'))
  if isempty(tsg)
    error('T/S data required:, no TSG data loaded')
  end
  t = interp1(tsg.dt, tsg.t, filt_avg.dt);
  s = interp1(tsg.dt, tsg.s, filt_avg.dt);
  switch di_method
    case 'interpolate'
      % Interpolate DI on Filtered
      %     + recommende if sensor drift with time
      di_pp = table(filt_avg.dt, 'VariableNames', {'dt'});
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
  
%   % Get beta salt from Zhang et al. 2009
%   beta_s = NaN(size(filt_qc.beta));
%   for i = 1:max(size(param.lambda))
%     for j = 1:size(t,1)
%       beta_s(j,i) = betasw_ZHH2009(param.lambda(i), t(j), param.theta, s(j)) - betasw_ZHH2009(param.lambda(i), t(j), param.theta, 0);
%     end
%   end

  switch di_method
    case {'interpolate', 'constant'}
      % Get beta salt from Zhang et al. 2009
      beta_s = NaN(size(filt_avg.beta));
      for j = 1:size(t,1)
        beta_s(j, :) = betasw_ZHH2009(param.lambda, t(j), param.theta, s(j)) - ...
          betasw_ZHH2009(param.lambda, t(j), param.theta, 0);
      end
      % Compute beta dissolved
      g = table(filt_avg.dt, 'VariableNames', {'dt'});
      g.betag = param.slope .* (filt_avg.beta - di_pp.beta) - beta_s ;
      % Propagate error
      %   Note: Error is not propagated through Scattering & Residual temperature
      %         correction as required by SeaBASS
      g.betag_sd = param.slope .* sqrt(filt_avg.beta_avg_sd + di_pp.beta_avg_sd);
    case 'SW_scattering'
      % Get beta salt from Zhang et al. 2009
      beta_sw = NaN(size(filt_avg.beta));
      for j = 1:size(t,1)
        beta_sw(j, :) = betasw_ZHH2009(param.lambda, t(j), param.theta, s(j));
      end
      % Compute beta dissolved
      g = table(filt_avg.dt, 'VariableNames', {'dt'});
      % g.betag = param.slope .* filt_qc.beta - beta_s ;
      g.betag = param.slope .* (filt_avg.beta - param.dark) - beta_sw ;
      % Propagate error
      %   Note: Error is not propagated through Scattering & Residual temperature
      %         correction as required by SeaBASS
      % g.betag_sd = param.slope .* sqrt(filt_qc.beta_avg_sd + di_pp.beta_avg_sd);
      g.betag_sd = filt_avg.beta_avg_sd;
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
  [~, g.gamma_bbg] = FitSpectra_HM2(param.lambda, g.bbg);  
  fprintf('Done\n')
else
  g = table();
end
end

function [poc, poc_lower, poc_upper, cphyto, cphyto_lower, cphyto_upper] = estimatePOC_Cphyto(bbp, lambda, method, conf_int)
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
if ~exist('conf_int','var') || isempty(conf_int)
  conf_int = true;
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


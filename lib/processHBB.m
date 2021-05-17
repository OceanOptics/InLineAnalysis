function [p, g] = processHBB(param, tot, filt, di, tsg, di_method)
% Note DI is not interpolated as it's assume to be stable in time
% BB3 parameters is a structure
%   param.lambda <1xM double> wavelength (nm)
%   param.theta <1x1 double> scattering angle (degree)
%   param.dark  <1xM double> dark (DEPRECATED)
%   param.slope <1xM double> slope (DEPRECATED)


% Interpolate filtered on Total
filt_interp = table(tot.dt, 'VariableNames', {'dt'});
filt_interp.beta = interp1(filt.dt, filt.beta, filt_interp.dt);
filt_interp.beta_avg_sd = interp1(filt.dt, filt.beta_avg_sd, filt_interp.dt);

% Compute beta particulate
p = table(tot.dt, 'VariableNames', {'dt'});
p.betap = tot.beta - filt_interp.beta;

% Interpolate X_p with values from Sullivan et al. 2013
if param.theta ~= 124
  theta_ref = 90:10:170;
  X_p_ref = [0.684 0.858 1.000 1.097 1.153 1.167 1.156 1.131 1.093];
  X_p = interp1(theta_ref, X_p_ref, param.theta, 'spline');
else
  X_p = 1.076; % Specific to ECO-BB3 with a scattering angle of 124
end

% Compute backscattering
p.bbp = 2 * pi * X_p .* p.betap;

% Propagate error
%   Note: Error is not propagated through Scattering & Residual temperature
%         correction as required by SeaBASS
% p.betap_sd = param.slope .* sqrt(tot.beta_avg_sd + filt_interp.beta_avg_sd);
p.betap_sd = sqrt(tot.beta_avg_sd + filt_interp.beta_avg_sd);
p.betap_n = tot.beta_avg_n;

% Derive Gamma_bbp (does not support NaN values)
% Correct bu on March 5, 2018, FitSpectra_HM2 does not accept NaNs in spectra
p(all(isnan(p.bbp),2),:)=[];
sel = ~any(isnan(p.bbp));
sel(param.lambda < 510 | param.lambda > 650) = false;
fprintf('Computing Gamma_bbp ... ')
[~, p.gamma_bbp] = FitSpectra_HM2(param.lambda(sel), p.bbp(:, sel));
fprintf('Done\n')

% Estimate POC and Cphyto from bbp
[p.poc, ~, ~, p.cphyto, ~, ~,] = estimatePOC_Cphyto(p.bbp, param.lambda, 'soccom');

% QC
p(any(p.bbp < 0,2),:) = [];

if nargout > 1 && nargin > 4
  t = interp1(tsg.dt, tsg.t, filt.dt);
  s = interp1(tsg.dt, tsg.s, filt.dt);
  switch di_method
    case 'interpolate'
      % Interpolate DI on Filtered
      %     + recommende if sensor drift with time
      di_pp = table(filt.dt, 'VariableNames', {'dt'});
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
%   beta_s = NaN(size(filt.beta));
%   for i = 1:max(size(param.lambda))
%     for j = 1:size(t,1)
%       beta_s(j,i) = betasw_ZHH2009(param.lambda(i), t(j), param.theta, s(j)) - betasw_ZHH2009(param.lambda(i), t(j), param.theta, 0);
%     end
%   end


  switch di_method
    case {'interpolate', 'constant'}
      % Get beta salt from Zhang et al. 2009
      beta_s = NaN(size(filt.beta));
      for j = 1:size(t,1)
        beta_s(j, :) = betasw_ZHH2009(param.lambda, t(j), param.theta, s(j)) - ...
          betasw_ZHH2009(param.lambda, t(j), param.theta, 0);
      end
      % Compute beta dissolved
      g = table(filt.dt, 'VariableNames', {'dt'});
      g.betag = (filt.beta - di_pp.beta) - beta_s ;
      % Propagate error
      %   Note: Error is not propagated through Scattering & Residual temperature
      %         correction as required by SeaBASS
      % g.betag_sd = param.slope .* sqrt(filt.beta_avg_sd + di_pp.beta_avg_sd);
      g.betag_sd = sqrt(filt.beta_avg_sd + di_pp.beta_avg_sd);
    case 'SW_scattering'
      % Get beta salt from Zhang et al. 2009
      beta_sw = NaN(size(filt.beta));
      for j = 1:size(t,1)
        beta_sw(j, :) = betasw_ZHH2009(param.lambda, t(j), param.theta, s(j));
      end
      % Compute beta dissolved
      g = table(filt.dt, 'VariableNames', {'dt'});
      % g.betag = param.slope .* filt.beta - beta_s ;
      g.betag = filt.beta - beta_sw ;
      % Propagate error
      %   Note: Error is not propagated through Scattering & Residual temperature
      %         correction as required by SeaBASS
      % g.betag_sd = param.slope .* sqrt(filt.beta_avg_sd + di_pp.beta_avg_sd);
      g.betag_sd = filt.beta_avg_sd;
    otherwise
      error('Method not supported.');
  end
  g.betag_n = filt.beta_avg_n;
  
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


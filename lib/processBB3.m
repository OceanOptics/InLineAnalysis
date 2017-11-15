function [p, g] = processBB3(param, tot, filt, di, tsg)
% Note DI is not interpolated as it's assume to be stable in time
% BB3 parameters is a structure
%   param.lambda <1x3 double> wavelength (nm)
%   param.theta <1x1 double> scattering angle (degree)
%   param.dark  <1x3 double> dark
%   param.slope <1x3 double> slope


% Interpolate filtered on Total
filt_interp = table(tot.dt, 'VariableNames', {'dt'});
filt_interp.beta = interp1(filt.dt, filt.beta, filt_interp.dt);
filt_interp.beta_avg_sd = interp1(filt.dt, filt.beta_avg_sd, filt_interp.dt);

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

if nargout > 1 && nargin > 3
  % Get beta salt from Zhang et al. 2009
  t = interp1(tsg.dt, tsg.t, filt.dt);
  s = interp1(tsg.dt, tsg.s, filt.dt);
  beta_s = NaN(size(filt.beta));
  for i = 1:size(param.lambda,2)
    for j = 1:size(t,1)
      beta_s(j,i) = betasw_ZHH2009(param.lambda(i), t(j), param.theta, s(j)) - betasw_ZHH2009(param.lambda(i), t(j), param.theta, 0);
    end
  end

  % Compute beta dissolved
  g = table(filt.dt, 'VariableNames', {'dt'});
  g.betag = param.slope .* (filt.beta - di.beta) - beta_s ;
  
  % Propagate error
  %   Note: Error is not propagated through Scattering & Residual temperature
  %         correction as required by SeaBASS
  g.betag_sd = param.slope .* sqrt(filt.beta_avg_sd + di.beta_avg_sd);
  g.betag_n = filt.beta_avg_n;
end
end

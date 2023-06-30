function [p] = processFL_old(param, tot, filt)
% FL parameters is a structure
%   param.dark  <1xN double> dark
%   param.slope <1xN double> slope
% Assume only particules >.2 um chla fluoresce

% Interpolate filtered on Total
filt_interp = table(tot.dt, 'VariableNames', {'dt'});
filt_interp.fchl = interp1(filt.dt, filt.fchl, filt_interp.dt);
filt_interp.fchl_avg_sd = interp1(filt.dt, filt.fchl_avg_sd, filt_interp.dt);

% Compute chlorophyll a fluorescence (dark independent)
p = table(tot.dt, 'VariableNames', {'dt'});
p.fchl = param.slope .* (tot.fchl - filt_interp.fchl);

% Propagate error
p.fchl_sd = param.slope .* sqrt(tot.fchl_avg_sd.^2 + filt_interp.fchl_avg_sd.^2);
p.fchl_n = tot.fchl_avg_n;

end
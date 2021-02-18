function [pd] = processCD(param, data)
% FL parameters is a structure
%   param.dark  <1xN double> dark
%   param.slope <1xN double> slope
% Assume no drift in time

% Compute chlorophyll a fluorescence (dark independent)
pd = table(data.dt, 'VariableNames', {'dt'});
pd.fdom = param.slope * ((data.fdom *4.33 - param.dark)/1000);

% Propagate error
pd.fdom_sd = param.slope * ((data.fdom_avg_sd *4.33 - param.dark)/1000);
pd.fdom_n = data.fdom_avg_n;

end
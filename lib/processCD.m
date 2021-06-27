function pd = processCD(param, data, diw)
% FL parameters is a structure
%   param.dark  <1xN double> dark
%   param.slope <1xN double> slope

% set gain
gain = 4.33;

% Compute CDOM fluorescence
pd = table(data.dt, 'VariableNames', {'dt'});

% method without DIW measurement and factory slope and dark
if isempty(diw)
  pd.fdom = param.slope * (data.fdom * gain - param.dark) / 1000;
  % Propagate error
  pd.fdom_sd = param.slope * (data.fdom_avg_sd * gain - param.dark) / 1000;
  pd.fdom_n = data.fdom_avg_n;
  
else % method with DIW measurement (dark not required)  
  % Consider the lowest DIW of the whole period as the dark 
  pd.fdom = (data.fdom * gain - min(diw.fdom)) / 1000;
  % Propagate error
  pd.fdom_sd = (data.fdom_avg_sd * gain - diw.fdom_avg_sd(diw.fdom == min(diw.fdom))) / 1000;
  pd.fdom_n = data.fdom_avg_n;
end
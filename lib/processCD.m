function pd = processCD(param, data, diw)
% FL parameters is a structure
%   param.dark  <1xN double> dark
%   param.slope <1xN double> slope

% set gain OLD INLININO
% gain = 4.33;
gain = 1;

% Compute CDOM fluorescence
pd = table(data.dt, 'VariableNames', {'dt'});

% method without DIW measurement and factory slope and dark
if isempty(diw)
%   pd.fdom = param.slope * (data.(param.varname) * gain - param.dark) / 1000;
  pd.fdom = param.slope * (data.(param.varname) * gain - param.dark);
  % Propagate error
  pd.fdom_sd = param.slope * (data.([param.varname '_avg_sd']) * gain - param.dark);
%   pd.fdom_sd = param.slope * (data.([param.varname '_avg_sd'] * gain - param.dark) / 1000;
  pd.fdom_n = data.([param.varname '_avg_n']);
  
else % method with DIW measurement (dark not required)
  % Consider the lowest DIW of the whole period as the dark
%   pd.fdom = (data.(param.varname) * gain - min(diw.(param.varname))) / 1000;
  pd.fdom = (data.(param.varname) * gain - min(diw.(param.varname)));
  % Propagate error
  pd.fdom_sd = (data.([param.varname '_avg_sd']) * gain - diw.([param.varname '_avg_sd'])(diw.(param.varname) == min(diw.(param.varname))));
%   pd.fdom_sd = (data.([param.varname '_avg_sd']) * gain - diw.([param.varname '_avg_sd'])(diw.(param.varname) == min(diw.(param.varname)))) / 1000;
  pd.fdom_n = data.([param.varname '_avg_n']);
end
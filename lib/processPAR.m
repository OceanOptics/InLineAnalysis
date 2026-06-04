function [a] = processPAR(param, data)
% parameters is a structure
%   param.scale <1xN double> slope
% data: table with fields dt, par, v

% TODO: move autoQC criteria to autoQC section
% specific TaraPacific auto QC based on instrument input voltage
if any(strcmp(data.Properties.VariableNames, 'v'))
  % QC data if input tension is not around 7.5v (± 0.3v) old logger
  data = data(7.2 <= data.v & data.v <= 7.8, :);
end
% Might want to add QC if temperature too high
par_varname = data.Properties.VariableNames{strcmpi(data.Properties.VariableNames, 'par')};

% Calibrate PAR
a = table(data.dt, 'VariableNames', {'dt'});
insitudark = min(data.(par_varname));
if insitudark < param.dark; dark = insitudark; else; dark = param.dark; end
a.par = (data.(par_varname) - dark)./ param.scale;

% Propagate error
a.par_sd = data.([par_varname '_avg_sd']) ./ param.scale;
a.par_n = data.([par_varname '_avg_n']);

% remove obvious bad data
sel_bad = any(a.par > 5000 | a.par < 0);
a(sel_bad,:) = [];
end
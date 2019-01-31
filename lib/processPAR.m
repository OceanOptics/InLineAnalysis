function [a] = processPAR(param, data)
% parameters is a structure
%   param.scale <1xN double> slope
% data: table with fields dt, par, v


% QC data if tension (volts) is not around 7.5v (± 0.3v)
sel = 7.2 <= data.v & data.v <= 7.8;
% Might want to add QC if temperature too high

% Calibrate PAR
a = table(data.dt(sel), 'VariableNames', {'dt'});
a.par = data.par(sel) ./ param.scale;

% Propagate error
a.par_sd = data.par_avg_sd(sel) ./ param.scale;
a.par_n = data.par_avg_n(sel);

end
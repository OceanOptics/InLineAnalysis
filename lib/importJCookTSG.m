function [ data ] = importJCookTSG( filename, verbose )
% import custom made TSG files made by Emmanuel on JCook

if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% Open file
try
  load(filename, 'Time', 'Lat', 'Lon', 'Cond', 'Temp_h', 'Sal', 'Temp_r');
catch
  error('Unable to open file: %s', filename);
end

% Build table
data = table(Time', Lat', Lon', Temp_h', Cond', Sal', Temp_r', 'VariableNames', ...
  {'dt', 'lat', 'lon', 'tcal', 'c', 's', 't'});
data(isnan(data.dt), :) = [];

if verbose; fprintf('Done\n'); end
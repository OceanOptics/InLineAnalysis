function [ data ] = importSBEformatTSG( filename, verbose )
% 20210507-030542-SBE45.sbe45
if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% % Get date from filename
% s = strsplit(filename, '_');
% dt_ref = s{end-1};

% Set parser
parser = '%s%{dd/MM/yy}D%{hh:mm:ss.SSS}T%s%f%f%f%f%f%f';

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end

% Read data
t = textscan(fid, parser, 'delimiter', ',');
% Close file
fclose(fid);

% Build table
data = table(t{2} + t{3}, t{6}, t{7}, t{8}, t{10}, 'VariableNames', ...
  {'dt', 'tcal', 'c', 's', 't'});
data.dt = datenum(data.dt);
data(isnan(data.dt), :) = [];

if verbose; fprintf('Done\n'); end
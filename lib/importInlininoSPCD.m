function [ data ] = importInlininoSPCD( filename, verbose )
% Import SPCD data logged with Inlinino
%
%%
if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% Set parser
parser = '%s%f';

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end

% Get header
hd = strip(strsplit(fgetl(fid), ','));
hd(strcmp(hd, 'time')) = {'dt'};
% get units skipping empty lines (bug in Inlinino)
unit = fgetl(fid);
while isempty(unit)
    unit = fgetl(fid);
end
unit = strip(strsplit(unit, ','));

% Read data
t = textscan(fid, parser, 'delimiter', ',');
% Close file
fclose(fid);

% Build table
data = table(datenum(t{1}), t{2}, 'VariableNames', {'dt', 'fdom'});
data(isnan(data.dt), :) = [];
data.Properties.VariableUnits = unit;

if verbose; fprintf('Done\n'); end
function [ data ] = importInlininoTSG( filename, verbose )
% IMPORTINLININO Import GPSSC701 data from csv files
% Author: Guillaume Bourdin
% Date: Dec 2021
%
% Input:
%   - filename: <char> filename including full path
%   - verbose (optional)
% 
% Example: [ data, lambda] = importInlininoGPSSC701( filename, verbose )
%%
if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% Get date from filename
s = strsplit(filename, '_');
dt_ref = s{end-1};

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end

% Get header
hd = strip(strsplit(fgetl(fid), ','));
% rename datetime
hd{strcmp(hd, 'time')} = 'dt';

% get units skipping empty lines (bug in old Inlinino)
unit = fgetl(fid);
while isempty(unit)
    unit = fgetl(fid);
end
% get units and lambda
unit = strip(strsplit(unit, ','));

% Set parser
parser = ['%s' repmat('%f', 1, size(hd,2)-1)];

% Read data
t = textscan(fid, parser, 'delimiter',',');
% Close file
fclose(fid);

% Build table
if all(contains(t{1}, '/'))
  data = table(datenum(t{1}, 'yyyy/mm/dd HH:MM:SS.FFF'), t{2}, t{3}, t{4}, t{5}, ...
    t{6}, 'VariableNames', hd);
else
  data = table(datenum(cellfun(@(x) [dt_ref x], t{1}, 'UniformOutput', false), 'yyyymmddHH:MM:SS.FFF'), ...
    t{2}, t{3}, t{4}, t{5}, t{6}, 'VariableNames', hd);
end
data.Properties.VariableUnits = unit;

% Remove last line if it's past midnight (bug in old Inlinino)
if ~isempty(data) && size(data,1) > 1
  if data.dt(end-1) > data.dt(end)
    data(end,:) = [];
  end
end

if verbose; fprintf('Done\n'); end




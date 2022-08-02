function [ data ] = importInlininoNMEA( filename, verbose )
% importInlininoNMEA Import NMEA data from csv files
% Author: Guillaume Bourdin
% Date: August 1st, 2022
%
% Input:
%   - filename: <char> filename including full path
%   - verbose (optional)
% 
% Example: [ data, lambda] = importInlininoNMEA( filename, verbose )
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
hd{strcmp(hd, 'time')} = 'dt';
hd{strcmp(hd, 'datetime')} = 'gps_dt';
hd{strcmp(hd, 'latitude')} = 'lat';
hd{strcmp(hd, 'longitude')} = 'lon';
% get units skipping empty lines (bug in old Inlinino)
unit = fgetl(fid);
while isempty(unit)
    unit = fgetl(fid);
end
% get units and lambda
unit = strip(strsplit(unit, ','));

% Set parser
parser = ['%s%s' repmat('%f', 1, size(hd,2)-2)];

% Read data
t = textscan(fid, parser, 'delimiter',',');
% Close file
fclose(fid);

% remove GPS time
unit = unit(~strcmp(hd, 'gps_dt'));
t = t(~strcmp(hd, 'gps_dt'));
hd = hd(~strcmp(hd, 'gps_dt'));

% Build table
dat = [];
for i = 1:size(hd, 2)
  if strcmp(hd{i}, 'dt')
    dat = [dat datenum(t{i}, 'yyyy/mm/dd HH:MM:SS.FFF')];
  elseif contains(hd{i}, 'swt')
    foo2 = t{i};
    foo2(contains(foo2, 'True')) = {'1'};
    foo2(contains(foo2, 'False')) = {'0'};
    dat = [dat logical(cell2mat(foo2))];
  else
    dat = [dat t{i}];
  end
end
data = array2table(dat, 'VariableNames', hd);
data.Properties.VariableUnits = unit;

% Remove last line if it's past midnight (bug in old Inlinino)
if ~isempty(data) && size(data,1) > 1
  if data.dt(end-1) > data.dt(end)
    data(end,:) = [];
  end
end

if verbose; fprintf('Done\n'); end




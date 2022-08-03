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

% filename = '/Volumes/Samsung_T5/Data/TaraMicrobiome/raw/TSG/SBE38+450091_20220324_150229.csv';

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
dat = [];
for i = 1:size(hd, 2)
  if strcmp(hd{i}, 'dt')
    if all(contains(t{1}, '/'))
      dat = [dat datenum(t{i}, 'yyyy/mm/dd HH:MM:SS.FFF')];
    else
      dat = [dat datenum(cellfun(@(x) [dt_ref x], t{1}, 'UniformOutput', false), 'yyyymmddHH:MM:SS.FFF')];
    end
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




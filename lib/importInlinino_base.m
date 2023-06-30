function [ data, lambda ] = importInlinino_base( filename, verbose )
% importInlininoNMEA Import data from basic Inlinino csv files
% Author: Guillaume Bourdin
% Date: June 30st, 2023
%
% Input:
%   - filename: <char> filename including full path
%   - verbose (optional)
% 
% Example: [ data, lambda] = importInlinino_base( filename, verbose )
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
pars = repmat({''}, size(hd));
% rename variables
for i = 1:size(hd, 2)
  if strcmp(hd{i}, 'time')
    hd{i} = 'dt';
    pars{i} = '%s';
  elseif strcmp(hd{i}, 'Switch')
    hd{i} = 'swt';
    pars{i} = '%s';
  elseif strcmp(hd{i}, 'Flow(0)')
    hd{i} = 'spd1';
    pars{i} = '%f';
  elseif strcmp(hd{i}, 'Flow(1)')
    hd{i} = 'spd2';
    pars{i} = '%f';
  elseif strcmp(hd{i}, 'Flow(2)')
    hd{i} = 'spd3';
    pars{i} = '%f';
  elseif strcmp(hd{i}, 'Flow(3)')
    hd{i} = 'spd4';
    pars{i} = '%f';
  else
    pars{i} = '%f';
  end
  hd{i} = strrep(strrep(hd{i}, '(', ''), ')', '');
end
% Set parser
parser = strjoin(pars, '');

% get units skipping empty lines (bug in old Inlinino)
unit = fgetl(fid);
while isempty(unit)
    unit = fgetl(fid);
end
% get units and lambda
unit = strip(strsplit(unit, ','));

% Read data
t = textscan(fid, parser, 'delimiter',',');
% Close file
fclose(fid);

% Build table
dat = [];
for i = 1:size(hd, 2)
  if strcmp(hd{i}, 'dt')
    if contains(t{1}, '/')
      dat = [dat datenum(t{i}, 'yyyy/mm/dd HH:MM:SS.FFF')];
    else
      dat = [dat datenum(cellfun(@(x) [dt_ref x], t{1}, 'UniformOutput', false), 'yyyymmddHH:MM:SS.FFF')];
    end
  elseif contains(hd{i}, 'swt')
    foo2 = t{i};
    foo2(contains(foo2, 'True')) = {'1'};
    foo2(contains(foo2, 'False')) = {'0'};
    dat = [dat logical(str2num(cell2mat(foo2)))];
  else
    dat = [dat t{i}];
  end
end

% detect if lambda to extract
hd_tx = cell(size(hd));
for i = 1:size(hd, 2)
  hd_tx{i} = hd{i}(isstrprop(hd{i},'alpha'));
  hd_digit{i} = hd{i}(isstrprop(hd{i},'digit'));
end
[~, d] = unique(hd_tx, 'first');
% find wavelength and extract
lambda = str2double(hd_digit(contains(hd, unique(hd_tx(not(ismember(1:numel(hd_tx), d)))))));

data = array2table(dat, 'VariableNames', hd);
data.Properties.VariableUnits = unit;

% Remove last line if it's past midnight (bug in old Inlinino)
if ~isempty(data) && size(data,1) > 1
  if data.dt(end-1) > data.dt(end)
    data(end,:) = [];
  end
end

if verbose; fprintf('Done\n'); end

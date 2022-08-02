function [ data ] = importInlininoFlowControl( filename, verbose )
% importInlininoFlowControl Import FlowControl data from csv files
% Author: Guillaume Bourdin
% Date: August 1st, 2022
%
% Input:
%   - filename: <char> filename including full path
%   - verbose (optional)
% 
% Example: [ data, lambda] = importInlininoFlowControl( filename, verbose )
%%
if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% % Get date from filename
% s = strsplit(filename, '_');
% dt_ref = s{end-1};

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end

% Get header
hd = strip(strsplit(fgetl(fid), ','));
for i = 1:size(hd, 2)
  if strcmp(hd{i}, 'time')
    hd{i} = 'dt';
  elseif strcmp(hd{i}, 'Switch')
    hd{i} = 'swt';
  elseif strcmp(hd{i}, 'Flow(0)')
    hd{i} = 'spd1';
  elseif strcmp(hd{i}, 'Flow(1)')
    hd{i} = 'spd2';
  elseif strcmp(hd{i}, 'Flow(2)')
    hd{i} = 'spd3';
  elseif strcmp(hd{i}, 'Flow(3)')
    hd{i} = 'spd4';
  end
end
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
unit = unit(~strcmp(hd, 'datetime'));
hd = hd(~strcmp(hd, 'datetime'));

dat = [];
for i = 1:size(hd, 2)
  if strcmp(hd{i}, 'dt')
    dat = [dat datenum(t{i}, 'yyyy/mm/dd HH:MM:SS.FFF')];
  elseif contains(hd{i}, 'swt')
    foo2 = t{i};
    foo2(contains(foo2, 'True')) = {'1'};
    foo2(contains(foo2, 'False')) = {'0'};
    dat = [dat logical(str2num(cell2mat(foo2)))];
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




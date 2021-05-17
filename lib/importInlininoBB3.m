function [ data ] = importInlininoBB3( filename, verbose )
%IMPORTINLININO Import BB3 data from EXPORTS

if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% Get date from filename
s = strsplit(filename, '_');
dt_ref = s{end-1};

% Set parser
parser = '%s%f%f%f';

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end

% Get header
hd = strip(strsplit(fgetl(fid), ','));
% get units skipping empty lines (bug in Inlinino)
unit = fgetl(fid);
while isempty(unit)
    unit = fgetl(fid);
end
% get units and lambda
unit = strip(strsplit(unit, ','));
unit(3:4) = [];
lambda = str2double(cellfun(@(c) strrep(c, 'beta', ''), hd(2:4), 'un', 0));

% Read data
t = textscan(fid, parser, 'delimiter',',');
% Close file
fclose(fid);

% Build table
if all(contains(t{1}, '/'))
    data = table(datenum(t{1}, 'yyyy/mm/dd HH:MM:SS.FFF'), [t{2:4}], ...
             'VariableNames', {'dt', 'beta'});
else
    data = table(datenum(cellfun(@(x) [dt_ref x], t{1}, 'UniformOutput', false), 'yyyymmddHH:MM:SS.FFF'),...
                 [t{2:4}], 'VariableNames', {'dt', 'beta'});
end
data.Properties.VariableUnits = unit;

% Remove last line if it's past midnight (Bug in Inlinino)
if ~isempty(data)
  if data.dt(end-1) > data.dt(end)
    data(end,:) = [];
  end
end

if verbose; fprintf('Done\n'); end




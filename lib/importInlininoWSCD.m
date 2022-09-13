function [ data ] = importInlininoWSCD( filename, verbose )
%IMPORTINLININO Import BB3 and WSCD data from NAAMES 3

if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% Get date from filename
s = strsplit(filename, '_');
dt_ref = s{end-1};

% Set parser
parser = '%s%f';

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end

% Get header
hd = strip(strsplit(fgetl(fid), ','));
% get units skipping empty lines (bug in old Inlinino)
unit = fgetl(fid);
while isempty(unit)
    unit = fgetl(fid);
end

% Read data
t = textscan(fid, parser, 'delimiter',',');
% Close file
fclose(fid);

% Build table
if all(contains(t{1}, '/'))
    data = table(datenum(t{1}, 'yyyy/mm/dd HH:MM:SS.FFF'), t{2}, ...
             'VariableNames', {'dt', 'fdom'});
else
    data = table(datenum(cellfun(@(x) [dt_ref x], t{1}, 'UniformOutput', false), 'yyyymmddHH:MM:SS.FFF'),...
             t{2}, 'VariableNames', {'dt', 'fdom'});
end


% Remove last line if it's past midnight (Bug in Inlinino)
if ~isempty(data)
  if data.dt(end-1) > data.dt(end)
    data(end,:) = [];
  end
end

if verbose; fprintf('Done\n'); end

end



function [ data ] = importTSGTeraTerm( filename, verbose )
%importBB3TeraTerm Import BB3 data logged with teraterm

if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% % Get date from filename
% s = strsplit(filename, '_');
% dt_ref = s{end-1};

% Set parser
parser = '%s%s%s%s%s';

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end

% Read data
t = textscan(fid, parser, 'delimiter', {',', ']'});
% Close file
fclose(fid);

% reformat data
t{1} = cellfun(@(x) strrep(x, '[', ''), t{1}, 'un', 0);
for i = 2:size(t, 2)
  t{i} = cell2mat(cellfun(@(s) str2double(s{end}), cellfun(@(x) strsplit(x, '= '), t{i}, 'un', 0), 'un', 0));
end

% Build table
n = min(cellfun('size', t, 1));
data = table(datenum(t{1}(1:n), 'yyyy-mm-dd HH:MM:SS.FFF'),...
             NaN(n, 1), NaN(n, 1), t{2}(1:n), t{3}(1:n), t{4}(1:n), t{5}(1:n), ...
             'VariableNames', {'dt', 'lat', 'lon', 'tcal', 'c', 's', 't'});

% Remove last line if it's past midnight (Bug in Inlinino)
if ~isempty(data)
  if data.dt(end-1) > data.dt(end)
    data(end,:) = [];
  end
end

if verbose; fprintf('Done\n'); end
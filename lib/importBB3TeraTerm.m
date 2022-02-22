function [ data ] = importBB3TeraTerm( filename, verbose )
%importBB3TeraTerm Import BB3 data logged with teraterm

if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% Get date from filename
s = strsplit(filename, '_');
dt_ref = s{end-1};

% Set parser
parser = '%s%s%d%f%d%f%d%f%d';

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end
% Skip header
fgetl(fid); fgetl(fid);
% Read data
t = textscan(fid, parser, 'delimiter','\t');
% Close file
fclose(fid);

% Build table
n = min([size(t{4}, 1), size(t{6}, 1), size(t{8}, 1)]);
data = table(datenum(strcat(t{1}(1:n), t{2}(1:n)), 'mm/dd/yyHH:MM:SS'),...
             [t{4}(1:n), t{6}(1:n), t{8}(1:n)], 'VariableNames', {'dt', 'beta'});

% Remove last line if it's past midnight (Bug in Inlinino)
if ~isempty(data)
  if data.dt(end-1) > data.dt(end)
    data(end,:) = [];
  end
end

if verbose; fprintf('Done\n'); end

end
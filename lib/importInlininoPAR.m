function [ data ] = importInlininoPAR( filename, verbose )
%IMPORTPAR Import PAR data from Tara Pacific

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
% Skip header
fgetl(fid); fgetl(fid);
% Read data
t = textscan(fid, parser, 'delimiter',',');
% Close file
fclose(fid);

% Build table
try 
    data = table(datenum(cellfun(@(x) [dt_ref x], t{1}, 'UniformOutput', false), 'yyyymmddHH:MM:SS.FFF'),...
             t{2}, t{3}, t{4}, 'VariableNames', {'dt', 'par', 't', 'v'});
catch
    data = table(datenum(t{1},'yyyy-mm-dd HH:MM:SS'),...
             t{2}, t{3}, t{4}, 'VariableNames', {'dt', 'par', 't', 'v'});
end

% Remove last line if it's past midnight (Bug in Inlinino)
if ~isempty(data)
  if data.dt(end-1) > data.dt(end)
    data(end,:) = [];
  end
end

if verbose; fprintf('Done\n'); end

end



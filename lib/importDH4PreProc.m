function [ data] = importDH4PreProc( filename, verbose )
%IMPORTDH4PREPROC Import data recorded with DH4 and extracted with matlab
%code for both BB3 and WSCD from NAAMES 1

if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% Get date from filename
s = strsplit(filename, '_');
dt_ref = s{end-1};

% Set parser
parser = '%f%f%f%f%f';

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end
% Skip header
fgetl(fid);
% Read data
t = textscan(fid, parser, 'delimiter','\t');
% Close file
fclose(fid);

% Build table
data = table(t{1}, [t{2:4}], t{5}, 'VariableNames', {'dt', 'beta', 'fdom'});

if verbose; fprintf('Done\n'); end
end


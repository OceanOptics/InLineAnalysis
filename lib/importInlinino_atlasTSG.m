function [ data ] = importInlinino_atlasTSG( filename, verbose )
%IMPORTINLININO Import miniTSG data from csv files

if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% Set parser
parser = '%s%f';

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end
% Get header
hd = strip(strsplit(fgetl(fid), ','));
hd(strcmp(hd, 'time')) = {'dt'};
% Skip empty lines (bug in Inlinino)
unit = fgetl(fid);
while isempty(unit)
    unit = fgetl(fid);
end

% Read data
t = textscan(fid, parser, 'delimiter',',');
% Close file
fclose(fid);

% Build table
data = table(datenum(t{1}, 'yyyy/mm/dd HH:MM:SS.FFF'), t{2}, 'VariableNames', hd);
data.Properties.VariableUnits = strip(strsplit(unit, ','));

if verbose; fprintf('Done\n'); end
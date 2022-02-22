function [ data ] = importInlininoLISSTcsv( filename, verbose )
%IMPORTINLININO Import LISST data from csv files logged using Inlinino

if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, filesep);
  fprintf('Importing %s ... ', foo{end});
end

% Set parser
parser = '%s%s%f%f%f%f%f%f%f';

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end
% Get header
hd = strip(strsplit(fgetl(fid), ','));
hd(strcmp(hd, 'time')) = {'dt'};
% get units skipping empty lines (bug in Inlinino)
unit = fgetl(fid);
while isempty(unit)
    unit = fgetl(fid);
end
% get counts angle
cangle = strsplit(unit, {', counts	angle=', ','});
counts_angle = str2double(strsplit(cangle{2}, ' '));

% Read data
t = textscan(fid, parser, 'delimiter',',');
% Format count data into matrix of double
t{2} = cell2mat(cellfun(@(c) str2double(c(2:end-1)), ...
    cellfun(@(x) strsplit(x, {'[', ' ', ']'}), t{2}, 'un', 0), 'un', 0));

% Close file
fclose(fid);

% Build table
data = table(datenum(t{1}, 'yyyy/mm/dd HH:MM:SS.FFF'), [t{2}], ...
             t{3}, t{4}, t{5}, t{6}, t{7}, t{8}, t{9}, 'VariableNames', hd);
data.Properties.VariableUnits = strip(strsplit(unit, ','));

% Remove last line if it's past midnight (Bug in Inlinino)
if ~isempty(data)
  if data.dt(end-1) > data.dt(end)
    data(end,:) = [];
  end
end

if verbose; fprintf('Done\n'); end




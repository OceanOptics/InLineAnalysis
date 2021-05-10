function [ data ] = importInlininoACScsv( filename, verbose )
%IMPORTINLININO Import ACS data from csv files

if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% Set parser
parser = '%s%f%s%s%f%f%s';

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
% get lambda
lambda = strsplit(unit, {', 1/m\tlambda=', ','});
lambda_a = strsplit(lambda{3}, ' ');
lambda_c = strsplit(lambda{4}, ' ');

% Read data
t = textscan(fid, parser, 'delimiter',',');
% Format a and c data into matrix of double
t{3} = cell2mat(cellfun(@(c) str2double(c(2:end-1)), ...
    cellfun(@(x) strsplit(x, {'[', ' ', ']'}), t{3}, 'un', 0), 'un', 0));
t{4} = cell2mat(cellfun(@(c) str2double(c(2:end-1)), ...
    cellfun(@(x) strsplit(x, {'[', ' ', ']'}), t{4}, 'un', 0), 'un', 0));
% Format flag into boolean
t{end} = strcmp(t{end}, 'True');

% Close file
fclose(fid);

% Build table
data = table(datenum(t{1}, 'yyyy/mm/dd HH:MM:SS.FFF'), t{2}, [t{3}], [t{4}], ...
             t{5}, t{6}, t{7}, 'VariableNames', hd);
data.Properties.VariableUnits = strip(strsplit(unit, ','));

% Remove last line if it's past midnight (Bug in Inlinino)
if ~isempty(data)
  if data.dt(end-1) > data.dt(end)
    data(end,:) = [];
  end
end

if verbose; fprintf('Done\n'); end




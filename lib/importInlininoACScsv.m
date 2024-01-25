function [ data, lambda_a, lambda_c ] = importInlininoACScsv( filename, verbose )
% IMPORTINLININO Import ACS data from csv files
% Author: Guillaume Bourdin
% Date: Dec 2020
%
% Input:
%   - filename: <char> filename including full path
%   - verbose (optional)
% 
% Example: [ data, lambda_a, lambda_c ] = importInlininoACScsv( filename, verbose )
%%

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

try
  % Get header
  hd = strip(strsplit(fgetl(fid), ','));
  hd(strcmp(hd, 'time')) = {'dt'};
  % get units skipping empty lines (bug in old Inlinino)
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
catch
  % Close file
  fclose(fid);
  warning('Unable to load data from file: %s\nFile might be corrupted or data line incomplete\nTrying to remove incomplete lines', filename);
  % Open file
  fid=fopen(filename);
  if fid==-1
    error('Unable to open file: %s', filename);
  end

  % Get header
  hd = strip(strsplit(fgetl(fid), ','));
  hd(strcmp(hd, 'time')) = {'dt'};
  % get units skipping empty lines (bug in old Inlinino)
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
  % Format a, search incomplete lines, and delete
  foo = cellfun(@(x) strsplit(x, {'[', ' ', ']'}), t{3}, 'un', 0);
  sz_ln = cell2mat(cellfun(@(s) size(s, 2), foo, 'un', 0));
  corrupted_row_a = sz_ln ~= size(lambda_a, 2)+2;
  if any(corrupted_row_a)
    for i = 1:size(t, 2)
      if size(corrupted_row_a, 1) == size(t{i}, 1)
        t{i}(corrupted_row_a) = [];
      end
    end
    foo(corrupted_row_a) = [];
  end
  t{3} = cell2mat(cellfun(@(c) str2double(c(2:end-1)), foo, 'un', 0));
  % Format c, search incomplete lines, and delete
  foo = cellfun(@(x) strsplit(x, {'[', ' ', ']'}), t{4}, 'un', 0);
  sz_ln = cell2mat(cellfun(@(s) size(s, 2), foo, 'un', 0));
  corrupted_row_c = sz_ln ~= size(lambda_c, 2)+2;
  if any(corrupted_row_c)
    for i = 1:size(t, 2)
      if size(corrupted_row_c, 1) == size(t{i}, 1)
        t{i}(corrupted_row_c) = [];
      end
    end
    foo(corrupted_row_c) = [];
  end
  t{4} = cell2mat(cellfun(@(c) str2double(c(2:end-1)), foo, 'un', 0));
  fprintf('Success!! %i "a" line deleted and %i "c" line deleted', sum(corrupted_row_a), sum(corrupted_row_c))
end
% Format flag into boolean
t{end} = strcmp(t{end}, 'True');

% Close file
fclose(fid);

% Build table
data = table(datenum(t{1}, 'yyyy/mm/dd HH:MM:SS.FFF'), t{2}, [t{3}], [t{4}], ...
           t{5}, t{6}, t{7}, 'VariableNames', hd);
data.Properties.VariableUnits = strip(strsplit(unit, ','));

% Remove last line if it's past midnight (bug in old Inlinino)
if ~isempty(data)
  if data.dt(end-1) > data.dt(end)
    data(end,:) = [];
  end
end

if verbose; fprintf('Done\n'); end




function [ data ] = importTaraTSG( filename, verbose )
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
parser = '%s%f%f%f';

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end
% Parse data
i = 1;
dt = NaN(8640,1); lat = NaN(8640,1); lon = NaN(8640,1); t = NaN(8640,1); s = NaN(8640,1);
c = NaN(8640,1); t2 = NaN(8640,1);
while ~feof(fid)
  l = fgetl(fid);
  l = strsplit(l, ',');
  foo = strsplit(l{1},'='); t(i) = str2double(foo{2});
  foo = strsplit(l{2},'='); c(i) = str2double(foo{2});
  foo = strsplit(l{3},'='); s(i) = str2double(foo{2});
  foo = strsplit(l{4},'='); t2(i) = str2double(foo{2});
  foo = strsplit(l{5},'='); foo = strsplit(foo{2},' '); lat(i) = str2double(foo{1}) + str2double(foo{2}) / 60;
  if strcmp(foo{3}, 'S'); lat(i) = - lat(i); end
  foo = strsplit(l{6},'='); foo = strsplit(foo{2},' '); lon(i) = str2double(foo{1}) + str2double(foo{2}) / 60;
  if strcmp(foo{3}, 'W'); lon(i) = - lon(i); end
  foo = strsplit(l{7},'='); foo2 = strsplit(l{8},'='); dt(i) = datenum([foo{2} foo2{2}], 'HHMMSSddmmyy');
  i = i+1;
end

% Close file
fclose(fid);

%datenum(cellfun(@(x) [dt_ref x], t{1}, 'UniformOutput', false), 'yyyymmddHH:MM:SS.FFF')

% Build table
data = table(dt, lat, lon, t, s, c, t2);
data(isnan(data.dt), :) = [];

if verbose; fprintf('Done\n'); end

end

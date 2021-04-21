function [ data ] = importTaraTSG( filename, verbose )

if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% % Get date from filename
% s = strsplit(filename, '_');
% dt_ref = s{end-1};

% Set parser
parser = 't1=%f, c1=%f, s=%f, t2=%f, lat=%u %f %c, lon=%u %f %c, hms=%6c, dmy=%6c';

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end

% Read data
t = textscan(fid, parser);
% Close file
fclose(fid);

% convert the dates to matlab datenum
dt = datenum(str2num(strcat('20', t{12}(:,5:6))), str2num(t{12}(:,3:4)), str2num(t{12}(:,1:2)), ...
        str2num(t{11}(:,1:2)), str2num(t{11}(:,3:4)), str2num(t{11}(:,5:6)));

% reformat data
lat = double(t{5})+ t{6}/60;
lon = double(t{8})+ t{9}/60;
lat(t{7}=='S') = -lat(t{7}=='S');
lon(t{10}=='W') = -lon(t{10}=='W');

% Build table
data = table(dt, lat, lon, t{1}, t{2}, t{3}, t{4}, 'VariableNames', ...
  {'dt', 'lat', 'lon', 'tcal', 'c', 's', 't'});
data(isnan(data.dt), :) = [];

if verbose; fprintf('Done\n'); end


% % Parse data DEPRECATED
% i = 1;
% dt = NaN(8640,1); lat = NaN(8640,1); lon = NaN(8640,1); t = NaN(8640,1); s = NaN(8640,1);
% c = NaN(8640,1); t2 = NaN(8640,1);
% while ~feof(fid)
%   l = fgetl(fid);
%   l = strsplit(l, ',');
%   foo = strsplit(l{1},'='); t(i) = str2double(foo{2});
%   foo = strsplit(l{2},'='); c(i) = str2double(foo{2});
%   foo = strsplit(l{3},'='); s(i) = str2double(foo{2});
%   foo = strsplit(l{4},'='); t2(i) = str2double(foo{2});
%   foo = strsplit(l{5},'='); foo = strsplit(foo{2},' '); lat(i) = str2double(foo{1}) + str2double(foo{2}) / 60;
%   if strcmp(foo{3}, 'S'); lat(i) = - lat(i); end
%   foo = strsplit(l{6},'='); foo = strsplit(foo{2},' '); lon(i) = str2double(foo{1}) + str2double(foo{2}) / 60;
%   if strcmp(foo{3}, 'W'); lon(i) = - lon(i); end
%   foo = strsplit(l{7},'='); foo2 = strsplit(l{8},'='); dt(i) = datenum([foo{2} foo2{2}], 'HHMMSSddmmyy');
%   i = i+1;
% end
% 
% % Close file
% fclose(fid);
% 
% %datenum(cellfun(@(x) [dt_ref x], t{1}, 'UniformOutput', false), 'yyyymmddHH:MM:SS.FFF')
% 
% % Build table
% data = table(dt, lat, lon, t, s, c, t2);
% data(isnan(data.dt), :) = [];



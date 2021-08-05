function [ data ] = importSBE37_TSG( filename, verbose )
% Author: Guillaume bourdin
% Date: August 4th, 2021
%
% Import TSG data from SBE37 logged with matlab
%%
if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% % Get date from filename
% s = strsplit(filename, '_');
% dt_ref = s{end-1};

% Set parser
parser = 'lat=%u %f %c, lon=%u %f %c, hms=%6c, dmy=%6c, t1=%f, c1=%f, s=%f, sigma=%f';

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
dt = datenum(str2num(strcat('20', t{8}(:,5:6))), str2num(t{8}(:,3:4)), str2num(t{8}(:,1:2)), ...
        str2num(t{7}(:,1:2)), str2num(t{7}(:,3:4)), str2num(t{7}(:,5:6)));

% reformat data
lat = double(t{1})+ t{2}/60;
lon = double(t{4})+ t{5}/60;
lat(t{3}=='S') = -lat(t{3}=='S');
lon(t{6}=='W') = -lon(t{6}=='W');

% Build table
data = table(dt, lat, lon, t{9}, t{10}, t{11}, t{12}, 'VariableNames', ...
  {'dt', 'lat', 'lon', 't', 'c', 's', 'sigma'});
data(isnan(data.dt), :) = [];

if verbose; fprintf('Done\n'); end


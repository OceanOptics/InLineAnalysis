function [data] = importRRevelleUnderway(filename, verbose)
%IMPORTRREVELLEUNDERWAY import underway MET files from R/V Roger Revelle
%   Written during EXPORTS

% filename = '/Users/nils/Data/EXPORTS/InLine/RRUnderway/raw/180810.MET';
% verbose = true;

if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% Init parser
n_column = 91;
% i_str = [1];
parser = repmat('%f',1,n_column);
% parser(i_str*2) = 's';

% Get date from filename
[~, yymmdd] = fileparts(filename);
dt_ref = datenum(['20' yymmdd], 'yyyymmdd');

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end

% Skip header lines
fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);
% Read data
t = textscan(fid, parser, 'delimiter','\t');
% Close file
fclose(fid);

% Get time from first column
HH = floor(t{1} / 10000);
MM = floor((t{1} - HH * 10000) / 100);
SS = (t{1} - HH * 10000 - MM * 100);

% Build table
% 14 PA Surface PAR uE/Second/Meter^2
% 31 TT-2 Temperature (degC)
% 33 SA-2 Salinity (PSU)
% 39 OX Oxygen (ml/l)
% 40 OS Oxygen Saturation (ml/l)
% 41 FL (ug/l)
% 53 LA Latitue (decimal degree)
% 54 LO Longitude (decimal degree)
% 58 ZD GPS DateTime GMT Secs Since 00:00:00 01/01/1970
% 55 GT GPS Time of Day GMT Secs 0-86400
data = table(dt_ref + datenum(0, 0, 0, HH, MM, SS),...
             t{53}, t{54}, t{31}, t{33}, t{39}, t{40}, t{41}, t{14},...
             'VariableNames',...
             {'dt', 'lat', 'lon', 't', 's', 'o2', 'o2_sat', 'fchl', 'par'});
           
% Replace -99 by NaN
for v = data.Properties.VariableNames; v = v{1};
  data.(v)(data.(v) == -99) = NaN;
end

% Remove lines for which temperature & salinity are empty
data(isnan(data.t) & isnan(data.s), :) = [];

if verbose; fprintf('Done\n'); end
end


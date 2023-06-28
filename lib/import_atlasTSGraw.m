function [ data ] = import_atlasTSGraw( filename, verbose )
%IMPORTINLININO Import miniTSG data from csv files

if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end
% Read header line by line
get_header=true;
% file_type = 'Separated';
delimit = ',';
while ~feof(fid) && get_header
  % get header
  hd = strip(strsplit(fgetl(fid), delimit));
  if strcmpi(hd{1}, 'time')
    get_header = false;
    format = strip(strsplit(fgetl(fid), delimit));
  elseif strcmpi(hd{1}, 'datetime') % CTDizzle format
    get_header = false;
    format = {'yyyy-mm-ddTHH:MM:SS', 'uS/cm', 'DegC', 'dBar', 'volts'};
%     format = {'yyyy-MM-dd''T''HH:mm:SS', 'uS/cm', 'DegC', 'dBar', 'volts'}; % input format for datetime
  end
end
% parser = ['%q' repmat(',%f', 1, size(hd{1}, 2)-1) '%[^\n\r]']; % MM-dd-yyyy HH:mm:ss.SSS
parser = ['%q' repmat('%f', 1, size(hd, 2)-1)]; % MM-dd-yyyy HH:mm:ss.SSS
t = textscan(fid, parser, 'Delimiter', delimit);

% Close file
fclose(fid);

% Build table
data = table();
data.dt = datenum(t{1}, format{1});
% data.dt = datetime(t{1}, 'InputFormat', ...
%   strrep(strrep(strrep(strrep(format{1}, ':SS.', ':ss.'), 'fff', 'SSS'), '/mm/', '/MM/'), ':MM:', ':mm:'));
data = [data array2table(horzcat(t{2:end}), 'VariableNames', hd(2:end))];
data.Properties.VariableUnits = format;

if verbose; fprintf('Done\n'); end
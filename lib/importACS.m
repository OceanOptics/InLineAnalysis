function [ data ] = importACS( filename, verbose )
%IMPORTACS Import ACS data recorded with Compass r2.1 in schedule mode

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
flag_header=true;
file_type = 'Manual';
while ~feof(fid) && flag_header
  % Get line
  l = fgetl(fid);
  % Split in parameter and value
  s = regexp(l,'\t','split'); % keep empty cell
  % Get file type
  if ~isempty(strfind(s{1}, 'This file was created from the following files'))
    file_type = 'Scheduled';
  end
  % Get each wavelength
  if strcmp(s{1}, 'Time(ms)')
    flag_header = false;
    i_c = cellfun(@any, strfind(s,'c'));
    if sum(i_c) == 0
      i_c = cellfun(@any, strfind(s,'C'));
    end
%     c_wv = cellfun(@(x) str2double(x(2:end)), {s{i_c}});
    i_a = cellfun(@any, strfind(s,'a'));
    if sum(i_a) == 0
      i_a = cellfun(@any, strfind(s,'A'));
    end
%     a_wv = cellfun(@(x) str2double(x(2:end)), {s{i_a}});
  end
end

% Check number of wavelength
if sum(i_a) ~= sum(i_c)
  error('Unexpected number of wavelength.\n');
else
  n_wv = sum(i_a);
end

% Build data parser
parser = ['%f' repmat('%f',1,n_wv) repmat('%f',1,n_wv) '%[^\n\r]'];

% x=cell2mat(textscan(fid,repmat('%f',1,40), ...
%                         'delimiter',' ', ...
%                         'multipledelimsasone',true));

% Parse data
t = textscan(fid, parser);
dt = t{1};
c = [t{i_c}];
a = [t{i_a}];

fclose(fid);

% build table
data = table(dt, c, a);

% Check if there is a jump in timestmap of compass
[m,i] = max(data.dt);
if i ~= size(data.dt,1)
  % There is a jump in the timestamp, correct for it
  delta = median(data.dt(2:i) - data.dt(1:i-1));
  data.dt(i+1:end) = data.dt(i+1:end) + data.dt(i) - data.dt(i+1) + delta;
end

% Update date & time
s = strsplit(filename, '/');
s = strsplit(s{end}, '_');
data.dt = datenum(s{2}(1:14), 'yyyymmddHHMMSS') + datenum(0,0,0, 0,0,(data.dt - data.dt(1))/1000);
if strcmp(file_type, 'Scheduled')
  % File time stamp is done at the end of the recording
  % Thereafter subtract need to substract the length of the recording to the timestamp
  data.dt = data.dt - (data.dt(end)-data.dt(1));
end

if verbose; fprintf('Done\n'); end

end

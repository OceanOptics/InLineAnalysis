function [ data ] = importACSBin( filename, devicefile, verbose )
%IMPORTACS Import ACS binary data recorded with Compass r2.1 using scheduler

if nargin < 3; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% Check if unpacked binary file already exist (.bin.dat)
if ~isfile([filename '.dat'])
  % Check if there is the file to unpack first  
  if ~isfile(filename)
    error('Unable to open file: %s', filename);
  end
  % Run prepacs.exe
  wk_dir = pwd();
  if ismac || isunix
    % OSX cmd
    cmd = ['!/usr/local/bin/wine ' wk_dir '/packages/prepACS/prepacs.exe ' devicefile ' ' filename ' ' filename '.dat'];
  elseif ispc
    % Windows cmd
    cmd = ['!' wk_dir '\packages\prepACS\prepacs.exe ' devicefile ' ' filename ' ' filename '.dat'];
  else
    error('Platform not supported by prepacs.exe.')
  end
%   fprintf('\n%s\n', cmd);
  eval(cmd);
end

% Open file
fid=fopen([filename '.dat']);
if fid==-1
  warning('Unable to open file: %s', [filename '.dat']);
  data = [];
  return
end

% Read header line by line
flag_header=true;
file_type = 'Scheduled'; % unpacked with prepACS.exe
while ~feof(fid) && flag_header
  % Get line
  l = fgetl(fid);
  % Split in parameter and value
  s = regexp(l,'\t','split'); % keep empty cell
  % Get each wavelength
  if strcmp(strtrim(s{1}), 'Time')
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
% if corrupted file, delete one row to have same number of row to merge
if max(cellfun('size',t,1))~=min(cellfun('size',t,1))
    toobig=cellfun('size',t,1)==max(cellfun('size',t,1));
    tb=t(toobig);
    notoobig=cellfun(@(c) c(1:end-1,:), tb, 'un',0 );
    t=[notoobig,t(:,~toobig)];
    sprintf('%s corrupted: last row deleted', filename);
end

% (cellfun(@(x) iscell(x),t)) = [];
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
s = strsplit(filename, filesep);
s = strsplit(s{end}, '_');
if strcmp(file_type, 'Scheduled'); si = 2; % Scheduled recording
else; si = 3; end% Manual recording (save at prompt when start acquisition with compass)
data.dt = datenum(s{si}(1:14), 'yyyymmddHHMMSS') + datenum(0,0,0, 0,0,(data.dt - data.dt(1))/1000);
if strcmp(file_type, 'Scheduled')
  % Scheduled recording
  % File time stamp is done at the end of the recording
  % Thereafter subtract need to substract the length of the recording to the timestamp
  data.dt = data.dt - (data.dt(end)-data.dt(1));
end

if verbose; fprintf('Done\n'); end

end

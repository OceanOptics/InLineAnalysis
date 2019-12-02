function [ data ] = importAC9wetview( filename, verbose )
%IMPORTACS Import ACS data recorded with Compass r2.1

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
flag_header=true; foo = 1;
% file_type = 'Manual';
while ~feof(fid) && flag_header && foo < 2
  % Get line
  l = fgetl(fid);
  % Split in parameter and value
  s = regexp(l,'\t','split'); % keep empty cell
  % Get start time
    if contains(s{1}, 'WetView')
    dt_st = datenum(s{2},'mm-dd-yyyy') +...
        datenum(0,0,0,str2double(s{3}(1:2)),str2double(s{3}(4:5)),str2double(s{3}(7:8)));
    end
    if foo == 2
        foo = foo +1;
    end
    
  % Get each wavelength
  if contains(l, 'acquisition binsize')
    foo = foo + 1;
  end
  i_c = [17 18 19 5 6 7 11 12 13]; % order wl c
  i_a = [8 9 10 14 15 16 2 3 4]; % order wl a

%   if foo == 3
%     i_c = cellfun(@any, strfind(s,'c'));
%     if sum(i_c) == 0
%       i_c = cellfun(@any, strfind(s,'C'));
%     end
%     i_a = cellfun(@any, strfind(s,'a'));
%     if sum(i_a) == 0
%       i_a = cellfun(@any, strfind(s,'A'));
%     end
%   end
end

% Check number of wavelength
% if sum(i_a) ~= sum(i_c)
%   error('Unexpected number of wavelength.\n');
% else
  n_wv = size(i_c,2);
% end

% Build data parser
parser = ['%f' repmat('%f',1,n_wv) repmat('%f',1,n_wv) '%[^\n\r]'];

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
data.dt = dt_st + datenum(0,0,0,0,0,(data.dt - data.dt(1))/1000);

if verbose; fprintf('Done\n'); end
end

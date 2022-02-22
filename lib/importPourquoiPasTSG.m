function [ data ] = importAtlantisTSG( filename, verbose )
%IMPORTUNDERWAY from R/V Atlantis
% Written for NAAMES 3
% Does not import seconds and milliseconds as data is binned every minute
if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% Init data arrays
% num_line = 1441; % 1 sample / min during 24 hour + 1 min margin
% dt = NaN(num_line, 1);
% lat = NaN(num_line, 1);
% lon = NaN(num_line, 1);
% s = NaN(num_line, 1);
% t = NaN(num_line, 1);
% fchl = NaN(num_line, 1);

% Init parser
n_column = 212;
parser = repmat('%s',1,n_column);

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end

% Skip header lines
fgetl(fid); fgetl(fid);
% Read data
t = textscan(fid, parser, 'delimiter','\t');
% Close file
fclose(fid);
% Build table
data = table(datenum(cell2mat(horzcat(t{1:2})), 'dd/mm/yyyyHH:MM:SS'),...
	pll(t{3}), pll(t{4}), str2double(strrep(t{27}, ',','.')),...
  str2double(strrep(t{25}, ',','.')),...
  str2double(strrep(t{29}, ',','.')),...
  'VariableNames', {'dt', 'lat', 'lon', 's', 't', 'density'});

% Remove last line if mid-night
if day(data.dt(end-1)) ~= day(data.dt(end))
  data(end,:) = [];
end

if verbose; fprintf('Done\n'); end

end

function num=pll(str)
  % Convert Lat or Lon to decimal degree
  num = NaN(size(str));
  for i=1:size(str,1)
    s = strsplit(str{i}, ' ');
    switch s{1}
      case {'N', 'E'}
        f=1;
      case {'S', 'W'}
        f=-1;
      otherwise
        error('Unknow direction of latitude');
    end
    num(i)=dm2degrees([f*str2double(s{2}(1:end-1)) str2double(strrep(s{3}(1:end-1), ',','.'))]);
  end
end


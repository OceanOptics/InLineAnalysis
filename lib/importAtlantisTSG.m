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
n_column = 37;
i_str = [1 2 5 35 36 37];
parser = repmat('%f',1,n_column);
parser(i_str*2) = 's';

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end

% Skip header lines
fgetl(fid); fgetl(fid);
% Read data
t = textscan(fid, parser, 'delimiter',',');
% Close file
fclose(fid);
% Build table
data = table(datenum(cell2mat(horzcat(t{1:2})), 'yyyy/mm/ddHH:MM:SS.FFF'),...
	t{3}, t{4}, t{29}, t{30}, t{33},'VariableNames',...
  {'dt', 'lat', 'lon', 's', 't', 'fchl'});

% % Read file line by line (very slow due to str2double conversion)
% flag_header=0;
% i=1;
% while ~feof(fid)
%   % Get line
%   l = fgetl(fid);
%   
%   % Read first line of header
%   if flag_header < 2
% %     header = s
%     flag_header = flag_header + 1;
%     continue
%   end
%   
%   % Split in parameter and value
%   sl = regexp(l,', ','split'); % keep empty cell
%   
%   % Parse line of data
%   dt(i) = datenum([sl{1} sl{2}(1:5)], 'yyyy/mm/ddHH:MM');
%   lat(i) = str2double(sl{3});
%   lon(i) = str2double(sl{4});
%   s(i) = str2double(sl{29});
%   t(i) = str2double(sl{30});
%   fchl(i) = str2double(sl{33});
%   
%   % Write in next raw of arrays
%   i=i+1;
% end

% fclose(fid);


% remove empty lines
% sel2rm = isnan(data.beta) & isnan(data.fdom);
% sel2rm = isnan(data.s) & isnan(data.t) & isnan(data.fchl);
% data(sel2rm, :) = [];

if verbose; fprintf('Done\n'); end

end



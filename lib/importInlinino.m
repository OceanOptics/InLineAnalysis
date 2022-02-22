function [ data ] = importInlinino( filename, verbose )
%IMPORTINLININO Import BB3 and WSCD data from NAAMES 3

if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% % Init data arrays
% num_line = 3601; % 1 sample / second during 1 hour + 1 second margin
% dt = NaN(num_line, 1);
% beta = NaN(num_line, 3);
% fdom = NaN(num_line, 1);

% Get date from filename
s = strsplit(filename, '_');
dt_ref = s{end-1};

% Set parser
parser = '%s%f%f%f%f';

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end
% Skip header
fgetl(fid); fgetl(fid);
% Read data
t = textscan(fid, parser, 'delimiter',',');
% Close file
fclose(fid);

% Build table
data = table(datenum(cellfun(@(x) [dt_ref x], t{1}, 'UniformOutput', false), 'yyyymmddHH:MM:SS.FFF'),...
             [t{2:4}], t{5}, 'VariableNames', {'dt', 'beta', 'fdom'});

% % Read file line by line (slow)
% flag_header=0;
% i=1;
% while ~feof(fid)
%   % Get line
%   l = fgetl(fid);
%   
%   % Split in parameter and value
%   s = regexp(l,',','split'); % keep empty cell
%   
%   % Read first line of header
%   if flag_header < 2
% %     header = s
%     flag_header = flag_header + 1;
%     continue
%   end
%   
%   % Parse line of data
%   dt(i) = datenum([dt_ref,s{1}], 'yyyymmddHH:MM:SS'); % .FFF
%   beta(i,1) = str2double(s{2});
%   beta(i,2) = str2double(s{3});
%   beta(i,3) = str2double(s{4});
%   fdom(i) = str2double(s{5});
%   
%   % Write in next raw of arrays
%   i=i+1;
% end
% 
% fclose(fid);
% 
% % build table
% data = table(dt, beta, fdom, 'VariableNames', {'dt', 'beta', 'fdom'});
% 
% % remove empty lines
% % sel2rm = isnan(data.beta) & isnan(data.fdom);
% sel2rm = isnan(data.beta(:,1)) & isnan(data.fdom);
% data(sel2rm, :) = [];

% Remove last line if it's past midnight (Bug in Inlinino)
if ~isempty(data)
  if data.dt(end-1) > data.dt(end)
    data(end,:) = [];
  end
end

if verbose; fprintf('Done\n'); end

end



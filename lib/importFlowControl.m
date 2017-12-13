function [ data ] = importFlowControl( filename, verbose )
%IMPORTFLOWCONTROL import data from Sequoia flow through system
% return time stamped switch position and flow speed 1
% 1 is total, 0 is filtered
% file lines
% yyyy-mm-dd HH:MM:SS UTC\t<switch position>\t0000\t<freq flowmeter 1>\t<freq flowmeter 2>\t<speed flowmeter 1(L/min)>\t<speed flowmeter 2(L/min)>\n
% 2017-05-12 17:10:26 UTC	0	0000	12.82	0.00	4.474	0.000

if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% Count number of line to init array
if (~ispc) 
  [status, cmdout]= system(['wc -l ' filename]);
  if(status~=1)
      scan_cell = textscan(cmdout,'%u %s');
      line_count = scan_cell{1} - 1; % for empty last line
  else
      fprintf(1,'\nFailed to find line count of %s\n', filename);
      line_count = -1;
  end
else
  fprintf(1,'\nScript not compatible with Windows\n');
  line_count = -1;
  return
end

% % Init data arrays
% dt = NaN(line_count, 1);
% swt = false(line_count, 1);
% spd = NaN(line_count, 1);

% Set parser
parser = '%s%d%d%f%f%f%f';

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end;

% Read data
t = textscan(fid, parser, 'delimiter','\t');
% Close file
fclose(fid);
% Build table
data = table(datenum(t{1}, 'yyyy-mm-dd HH:MM:SS UTC'), logical(t{2}), t{6},...
             'VariableNames', {'dt', 'swt', 'spd'});
% % Read file line by line (slow)
% i=1;
% while ~feof(fid)
%   % Get line
%   l = fgetl(fid);
%   
%   % Split in parameter and value
%   s = regexp(l,'\t','split'); % keep empty cell
%   
%   % Parse line of data
%   dt(i) = datenum(s{1}, 'yyyy-mm-dd HH:MM:SS UTC');
%   swt(i) = logical(str2double(s{2}));
%   spd(i) = str2double(s{6});
%   
%   % Write in next raw of arrays
%   i=i+1;
% end


% build table
% data = table(dt, swt, spd, 'VariableNames', {'dt', 'swt', 'spd'});

if verbose; fprintf('Done\n'); end

end

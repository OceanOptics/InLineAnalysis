function [ data ] = iRead( fun, dirname_in, dirname_out, prefix, dt, software, force, nowrite, verbose, read_margin, postfix )
%IMPORTALLUNDERWAY import underway data from all files matching regex settings
%   in specified directory. Loaded files are saved as mat files for faster run.
%
% INPUT:
%    - fun <@function> function to be run to import data
%    - dirname <string> directory in which to look for files
%    - prefix <string> prefix of files to look for (e.g. 'Inlinino_')
%    - dt <2x1 datenum> date & time of start and end of dataset imported
%        dt(1) start of data set (1 hour margin is substracted);
%        dt(2) end of data set (1 hour margin is added);
%    - postfix <string> postfix of files to look for (e.g. '.csv');
%    - dtformat <'yyyymmdd'|'yyyydoy'> format of date in filenames to be
%         imported
%    - force <boolean> force to import data (no loading from previous files)
%    - nowrite <boolean> no files will be writen to ./mat/ when import data
%    - verbose <boolean> nothing to write
%
% OUTPUT:
%    - data <format of data output by @fun> all data from files
%    - write files if nowrite is set to false
if nargin < 9; verbose = false; end
if nargin < 10; read_margin = true; end
if nargin < 11; postfix = ''; end

dir_in = dirname_in;
dir_out = dirname_out;

% switch dtformat
%   case 'yyyymmdd'
%     dtformat = @dt_yyyymmdd;
%   case 'yymmdd'
%     dtformat = @dt_yymmdd;
%   case 'yyyydoy'
%     dtformat = @dt_yyyydoy;
%   otherwise
%     error('Add your own dtformat subfunction in iRead');
% end

% Make sure dt is in day (0 hours 0 min and 0 seconds)
dt = floor(dt);

gdata = [];
for i=1:length(dt)
  fn_out = [prefix dt_yyyymmdd(dt(i)) postfix '.mat'];
  if ~force && exist([dir_out fn_out], 'file')
    if verbose; fprintf('Loading %s... ', fn_out); end
    load([dir_out fn_out]);
    if verbose; fprintf('Done\n'); end
    gdata = [gdata; data];
  else
    % List files matching date and prefix
    l = list_files_from_software(software, dir_in, prefix, dt(i), postfix);
    % Check if found files
    if isempty(l)
      fprintf('WARNING: No files found on %s\n', datestr(dt(i)));
    else
      % Import data from selection
      ddata = [];
      for j=1:size(l,1)
        foo = fun([dir_in l{j}], verbose);
        ddata = [ddata; foo];
      end
      % Keep only data of day
      sel = dt(i) <= ddata.dt & ddata.dt < dt(i) + 1;
      ddata = ddata(sel,:);
      % Write data of day
      if ~nowrite
        if ~isdir(dir_out); mkdir(dir_out); end
        data = ddata;
        if verbose; fprintf('Saving %s... ', fn_out); end
        save([dir_out fn_out], 'data');
        if verbose; fprintf('Done\n'); end
      end
      % Add data of day to the global dataset
      gdata = [gdata; ddata];
    end
  end
end

if read_margin
  % Load margin to dataset (calling myself)
  margin = 1/24; % day
  if verbose; fprintf('Reading margin ... \n'); end
  pre_data = iRead( fun, dirname_in, dirname_out, prefix, dt(1)-margin, software, force, nowrite, verbose, false, postfix );
  if ~isempty(pre_data)
    pre_data = pre_data(dt(1)-margin <= pre_data.dt,:);
  end
  post_data = iRead( fun, dirname_in, dirname_out, prefix, dt(end)+1+margin, software, force, nowrite, verbose, false, postfix );
  if ~isempty(post_data)
    post_data = post_data(post_data.dt <= dt(end)+1+margin,:);
  end
  if verbose; fprintf('Reading margin ... [Done]\n'); end
  % Adding margin to dataset
  data = [pre_data; gdata; post_data];
else
  % Export dataset
  data = gdata;
end

end

function [filenames] = list_files_from_software(software, dir_in, prefix, dt, postfix)
% dt <1x1 datenum> day of data to import
  switch software
    case {'Compass_2.1rc_scheduled', 'Compass_2.1rc'}
      % Compass does not reset files at mid-night thereafter some data from the
      % selected day might be in the first file of the following day
      % List all files in directory
      l = dir([dir_in filesep prefix '*' postfix '.dat']);
      if ~isempty(l)
        % Get date of all files
        n = length(prefix);
        l_dt = datenum(cellfun(@(x) x(n+1:n+14), {l.name}, 'UniformOutput', false), 'yyyymmddHHMMSS');
        % Get selection of files to import
        sel = dt <= l_dt & l_dt <= dt + 1 + 1/24; % Add one hour margin
        % Return selected filenames
        filenames = {l(sel).name}';
      else
        warning(['No files found for ' software]);
      end
    case 'Inlinino'
      % List all files in directory
      l = dir([dir_in filesep prefix dt_yyyymmdd(dt) '*' postfix '.csv']);
      filenames = {l.name}';
    case 'FlowControl'
      % List all files in directory
      l = dir([dir_in filesep prefix dt_yyyydoy(dt) '*' postfix '.log']);
      filenames = {l.name}';
    case 'DH4PreProc'
      % Get day of data +/- 1 day
      lp = dir([dir_in filesep prefix dt_doy(dt-1) '*' postfix '.dat']);
      l = dir([dir_in filesep prefix dt_doy(dt) '*' postfix '.dat']);
      la = dir([dir_in filesep prefix dt_doy(dt+1) '*' postfix '.dat']);
      filenames = {lp.name, l.name, la.name}';
    case 'AtlantisTSG'
      % List all files in directory
      l = dir([dir_in filesep prefix dt_yymmdd(dt) '*' postfix '.csv']);
      filenames = {l.name}';
    case 'PourquoiPasTSG'
      % List all files in directory
      l = dir([dir_in filesep dt_yyyymmdd(dt) '*' postfix '.csv']);
      filenames = {l.name}';
    case 'TeraTerm'
      % TeraTerm filenames correspond to the time at which logging started
      % Data from a given day could be in any file preceding that date
      % List all files in directory
      l = dir([dir_in filesep prefix '*' postfix '.log*']);
      % Get date of all files
      n = length(prefix);
      l_dt = floor(datenum(cellfun(@(x) x(n+1:n+15), {l.name}, 'UniformOutput', false), 'yyyymmdd_HHMMSS'));
      % Select the earliest file from the previous date (within 10 days)
      i = 1;
      while ~any(l_dt == dt - i) && i < 10; i = i + 1; end
      % Select previous valid day (within 10 days) and asked dt
      sel = l_dt == dt - i | l_dt == dt;
      % Return selected filenames
      filenames = {l(sel).name}';
    otherwise
      error('Software not supported: %s.', software);
  end

end

function [str] = dt_yyyymmdd(dt)
  str = datestr(dt, 'yyyymmdd');
end

function [str] = dt_yymmdd(dt)
  str = datestr(dt, 'yymmdd');
end

function [str] = dt_yyyydoy(dt)
  str = sprintf('%d%03d',year(dt),datevec2doy(datevec(dt)));
end

function [str] = dt_doy(dt)
  str = sprintf('%03d',datevec2doy(datevec(dt)));
end
function [ data ] = iRead( fun, dirname_in, dirname_out, prefix, dt, software, ...
    force, nowrite, verbose, read_margin, postfix, parallel_flag, otherarg1, otherarg2 )
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
  if nargin < 12; parallel_flag = Inf; end
  if nargin < 13; otherarg1 = {}; end
  if nargin < 14; otherarg2 = {}; end
  
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
  if ~isdatetime(dt)
    dt = datetime(dt, 'ConvertFrom', 'datenum');
  end
  % dt = floor(dt);
  dt = dateshift(dt, 'Start', 'Day');
  
  gdata = [];
  for i=1:length(dt)
    fn_out = [prefix dt_yyyymmdd(dt(i)) postfix '.mat'];
    if ~force && exist(fullfile(dir_out, fn_out), 'file')
      if verbose; fprintf('Loading %s... ', fn_out); end
      load(fullfile(dir_out, fn_out), 'data');
      if verbose; fprintf('Done\n'); end
      if ~isdatetime(data.dt)
        data.dt = datetime(data.dt, 'ConvertFrom', 'datenum');
      end
      % check nb of variable and add data of day to the global dataset
      gdata = check_nb_variables([{gdata}; {data}], dt(i));
    else
      % List files matching date and prefix
      l = list_files_from_software(software, dir_in, prefix, dt(i), postfix);
      % Check if found files
      if isempty(l)
        fprintf('WARNING: No files found on %s\n', dt(i));
      else
        % % Import data from selection
        % ddata = [];
        % % for j=1:size(l,1)
        % parfor (j=1:size(l,1), parallel_flag)
        %   if dir(fullfile(dir_in, l{j})).bytes > 0
        %     if isempty(otherarg1) && isempty(otherarg2)
        %       foo = fun(fullfile(dir_in, l{j}), verbose);
        %     elseif isempty(otherarg2)
        %       foo = fun(fullfile(dir_in, l{j}), otherarg1, verbose);
        %     else
        %       foo = fun(fullfile(dir_in, l{j}), otherarg1, otherarg2, verbose);
        %     end
        %     [ddata, foo] = check_nb_variables(ddata, foo, dt(i));
        %     ddata = [ddata; foo];
        %   else
        %     warning('Empty file: %s - ignored\n', fullfile(dir_in, l{j}))
        %   end
        % end
  
        % Get file info once before the loop
        file_paths = fullfile(dir_in, l);
        file_info = cellfun(@dir, file_paths, 'UniformOutput', false);
        file_sizes = cellfun(@(f) ~isempty(f) && f(1).bytes > 0, file_info);
        % Import data from selection
        nFiles = size(l, 1);
        ddata_cell = cell(nFiles, 1); % Use a cell array to store output from each iteration
        % for j=1:nFiles
        parfor (j = 1:nFiles, parallel_flag)
          if file_sizes(j)
            if isempty(otherarg1) && isempty(otherarg2)
              foo = fun(fullfile(dir_in, l{j}), verbose);
            elseif isempty(otherarg2)
              foo = fun(fullfile(dir_in, l{j}), otherarg1, verbose);
            else
              foo = fun(fullfile(dir_in, l{j}), otherarg1, otherarg2, verbose);
            end
            ddata_cell{j} = foo;
          else
            warning('Empty or missing file: %s - ignored\n', file_paths{j})
            ddata_cell{j} = []; % Empty placeholder
          end
        end
        % check nb of variable
        if any(~cellfun('isempty', ddata_cell))
          ddata = check_nb_variables(ddata_cell, dt(i));
        else
          ddata = [];
        end
        % Keep only data of day
        if ~any(strcmp(software, {'internal_logger', 'TeraTerm', 'WetView', 'Compass_2.1rc'}))
          if ~isempty(ddata)
            % sel = ddata.dt >= dt(i) & ddata.dt < dt(i) + 1;
            sel = ddata.dt >= dt(i) & ddata.dt < dt(i) + days(1);
            ddata = ddata(sel,:);
          end
        end
        if ~isempty(ddata)
          % Write data of day
          if ~nowrite
            if ~isfolder(dir_out); mkdir(dir_out); end
            data = ddata;
            if verbose; fprintf('Saving %s... ', fn_out); end
            save(fullfile(dir_out, fn_out), 'data');
            if verbose; fprintf('Done\n'); end
          end
        else
          fprintf('WARNING: No data found on %s\n', dt(i));
        end
        % check nb of variable and add data of day to the global dataset
        gdata = check_nb_variables([{gdata}; {ddata}], dt(i));
      end
    end
  end
  
  if read_margin
    % Load margin to dataset (calling myself)
    % margin = 1/24; % 1 hour
    margin = hours(1)+minutes(1); % 1 hour + 1 minute
    if verbose; fprintf('Reading margin ... \n'); end
    if isempty(dt); error("'days2run' is empty."); end
    % pre_data = iRead( fun, dirname_in, dirname_out, prefix, dt(1)-1, ...
    %     software, force, nowrite, verbose, false, postfix, parallel_flag, otherarg1, otherarg2 );
    pre_data = iRead( fun, dirname_in, dirname_out, prefix, dt(1)-days(1), ...
        software, force, nowrite, verbose, false, postfix, parallel_flag, otherarg1, otherarg2 );
    if ~isempty(pre_data)
      pre_data = pre_data(dt(1)-margin <= pre_data.dt,:);
    end
    % post_data = iRead( fun, dirname_in, dirname_out, prefix, dt(end)+1+margin, ...
    %     software, force, nowrite, verbose, false, postfix, parallel_flag, otherarg1, otherarg2 );
    post_data = iRead( fun, dirname_in, dirname_out, prefix, dt(end)+days(1)+margin, ...
        software, force, nowrite, verbose, false, postfix, parallel_flag, otherarg1, otherarg2 );
    if ~isempty(post_data)
      % post_data = post_data(post_data.dt <= dt(end)+1+margin,:);
      post_data = post_data(post_data.dt <= dt(end)+days(1)+margin,:);
    end
    if verbose; fprintf('Reading margin ... [Done]\n'); end
    % check nb of variable and add margin data to the global dataset
    gdata = check_nb_variables([{pre_data}; {gdata}], dt(i));
    data = check_nb_variables([{gdata}; {post_data}], dt(i));
    if ~isempty(data)
      data = sortrows(data, 'dt');
    end
  else
    % Export dataset
    data = gdata;
  end
end

function [filenames] = list_files_from_software(software, dir_in, prefix, dt, postfix)
  % dt <1x1 datetime> day of data to import
  if ~isdatetime(dt)
    dt = datetime(dt, 'ConvertFrom', 'datenum');
  end
  switch software
    case {'WetView', 'Compass_2.1rc_scheduled', 'Compass_2.1rc', ...
            'Compass_2.1rc_scheduled_bin'}
      % Compass does not reset files at mid-night thereafter some data from the
      % selected day might be in the first file of the following day
      % List all files in directory
      switch software
        case 'Compass_2.1rc_scheduled_bin'
          fprintf('Looking for %s\n', fullfile(dir_in, [prefix '*' postfix '.bin']))
          l = dir(fullfile(dir_in, [prefix '*' postfix '.bin']));
        otherwise
          fprintf('Looking for %s\n', fullfile(dir_in, [prefix '*' postfix '.dat']))
          l = dir(fullfile(dir_in, [prefix '*' postfix '.dat']));
      end
      if ~isempty(l)
        % Get date of all files
        n = length(prefix);
        % l_dt = datenum(cellfun(@(x) x(n+1:n+14), {l.name}, 'UniformOutput', false), 'yyyymmddHHMMSS');
        l_dt = datetime(cellfun(@(x) x(n+1:n+14), {l.name}, 'UniformOutput', false), 'InputFormat', 'yyyyMMddHHmmss');
        % Get selection of files to import
        switch software
            case 'WetView'
                % sel = dt - 7 <= l_dt & l_dt <= dt + 1 + 7; % Add 7 day margin
                sel = dt - days(7) <= l_dt & l_dt <= dt + days(8); % Add 7 day margin
            case 'Compass_2.1rc'
                % sel = dt - 6 <= l_dt & l_dt <= dt + 1 + 6; % Add 6 day margin
                sel = dt - days(6) <= l_dt & l_dt <= dt + days(7); % Add 6 day margin
            otherwise
                % sel = dt <= l_dt & l_dt <= dt + 1 + 1/24; % Add one hour margin
                sel = dt <= l_dt & l_dt <= dt + days(1) + hours(1); % Add one hour margin
        end
        % Return selected filenames
        filenames = {l(sel).name}';
      else
        warning(['No files found for ' software]);
        filenames = [];
      end
    case {'Inlinino', 'Inlinino_atlasTSG', 'InlininoADU100'}
      % List all files in directory
      filenames = struct2table(dir(fullfile(dir_in, [prefix dt_yyyymmdd(dt) '*' postfix '.csv']))).name;
      fprintf('Looking for %s\n', fullfile(dir_in, [prefix dt_yyyymmdd(dt) '*' postfix '.csv']))
      if isempty(filenames)
        filenames = struct2table(dir(fullfile(dir_in, [prefix dt_yyyymmdd(dt) '*' postfix '.raw']))).name;
        fprintf('No file found: looking for %s instead\n', fullfile(dir_in, [prefix dt_yyyymmdd(dt) '*' postfix '.raw']))
      end
      if isempty(filenames)
        filenames = struct2table(dir(fullfile(dir_in, [prefix dt_yyyymmdd(dt) '*' postfix '.log']))).name;
        fprintf('No file found: looking for %s instead\n', fullfile(dir_in, [prefix dt_yyyymmdd(dt) '*' postfix '.log']))
      end
    case 'SBE45software'
      % List all files in directory
      fprintf('Looking for %s\n', fullfile(dir_in, [prefix dt_yyyymmdd(dt) '*' postfix '.sbe45']))
      l = dir(fullfile(dir_in, [prefix dt_yyyymmdd(dt) '*' postfix '.sbe45']));
      filenames = {l.name}';
    case 'matlab_Emmanuel'
      % List all files in directory
      fprintf('Looking for %s\n', fullfile(dir_in, [prefix dt_yyyymmdd(dt) '*' postfix '.mat']))
      l = dir(fullfile(dir_in, [prefix dt_yyyymmdd(dt) '*' postfix '.mat']));
      filenames = {l.name}';
    case 'FlowControl'
      % List all files in directory
      fprintf('Looking for %s\n', fullfile(dir_in, [prefix dt_yyyydoy(dt) '*' postfix '.log']))
      l = dir(fullfile(dir_in, [prefix dt_yyyydoy(dt) '*' postfix '.log']));
      filenames = {l.name}';
    case 'DH4PreProc'
      fprintf('Looking for %s\n', fullfile(dir_in, [prefix dt_doy(dt) '*' postfix '.dat']))
      l = dir(fullfile(dir_in, [prefix dt_doy(dt) '*' postfix '.dat']));
      % Get day of data +/- 1 day
      % fprintf('Looking for %s\n', fullfile(dir_in, [prefix dt_doy(dt-1) '*' postfix '.dat']))
      % lp = dir(fullfile(dir_in, [prefix dt_doy(dt-1) '*' postfix '.dat']));
      % fprintf('Looking for %s\n', fullfile(dir_in, [prefix dt_doy(dt+1) '*' postfix '.dat']))
      % la = dir(fullfile(dir_in, [prefix dt_doy(dt+1) '*' postfix '.dat']));
      fprintf('Looking for %s\n', fullfile(dir_in, [prefix dt_doy(dt-days(1)) '*' postfix '.dat']))
      lp = dir(fullfile(dir_in, [prefix dt_doy(dt-days(1)) '*' postfix '.dat']));
      fprintf('Looking for %s\n', fullfile(dir_in, [prefix dt_doy(dt+days(1)) '*' postfix '.dat']))
      la = dir(fullfile(dir_in, [prefix dt_doy(dt+days(1)) '*' postfix '.dat']));
      filenames = {lp.name, l.name, la.name}';
    case 'AtlantisTSG'
      % List all files in directory
      fprintf('Looking for %s\n', fullfile(dir_in, [prefix dt_yymmdd(dt) '*' postfix '.csv']))
      l = dir(fullfile(dir_in, [prefix dt_yymmdd(dt) '*' postfix '.csv']));
      filenames = {l.name}';
    case 'PourquoiPasTSG'
      % List all files in directory
      fprintf('Looking for %s\n', fullfile(dir_in, [dt_yyyymmdd(dt) '*' postfix '.csv']))
      l = dir(fullfile(dir_in, [dt_yyyymmdd(dt) '*' postfix '.csv']));
      filenames = {l.name}';
    case 'RRevelleUnderway'
      % List all files in directory
      fprintf('Looking for %s\n', fullfile(dir_in, [prefix dt_yymmdd(dt) '*' postfix '.MET']))
      l = dir(fullfile(dir_in, [prefix dt_yymmdd(dt) '*' postfix '.MET']));
      filenames = {l.name}';
    case 'MatlabTSG'
      % List all files in directory
      fprintf('Looking for %s\n', fullfile(dir_in, [prefix dt_yyyymmdd(dt) '*' postfix '.txt']))
      l = dir(fullfile(dir_in, [prefix dt_yyyymmdd(dt) '*' postfix '.txt']));
      filenames = {l.name}';
    case 'internal_logger'
      % List all files in directory
      fprintf('Looking for %s\n', fullfile(dir_in, [prefix dt_yyyymmdd(dt) '*' postfix '.RBN']))
      l = dir(fullfile(dir_in, [prefix dt_yyyymmdd(dt) '*' postfix '.RBN']));
      filenames = {l.name}';
    case 'TeraTerm'
      % TeraTerm filenames correspond to the time at which logging started
      % Data from a given day could be in any file preceding that date
      % List all files in directory
      fprintf('Looking for %s\n', fullfile(dir_in, [prefix '*' postfix '.log*']))
      l = dir(fullfile(dir_in, [prefix '*' postfix '.log*']));
      % Get date of all files
      n = length(prefix);
      % l_dt = floor(datenum(cellfun(@(x) x(n+1:n+15), {l.name}, 'UniformOutput', false), 'yyyymmdd_HHMMSS'));
      l_dt = dateshift(datetime(cellfun(@(x) x(n+1:n+15), {l.name}, 'UniformOutput', false), ...
        'InputFormat', 'yyyyMMdd_HHmmss'), 'Start', 'Day');
      % Find adjacent date before the asked date (within 10 days)
      i = 1;
      % Select previous valid day (within 10 days) and asked dt
      % while ~any(l_dt == dt - i) && i < 10; i = i + 1; end
      % sel = l_dt == dt - i | l_dt == dt;
      while ~any(l_dt == dt - days(i)) && i < 10; i = i + 1; end
      sel = l_dt == dt - days(i) | l_dt == dt;
      % Return selected filenames
      filenames = {l(sel).name}';
    case 'ALFA_LabView_m'
      % ALFA LabView filenames correspond to the time at which logging started
      % Data from a given day could be in any file preceding that date
      % List all files in directory
      fprintf('Looking for %s\n', fullfile(dir_in, [prefix '*' postfix '_m.txt']))
      l = dir(fullfile(dir_in, [prefix '*' postfix '_m.txt']));
      % Get date of all files
      n = length(prefix);
      % l_dt = floor(datenum(cellfun(@(x) x(n+1:n+15), {l.name}, 'UniformOutput', false), 'yyyymmdd_HHMMSS'));
      l_dt = dateshift(datetime(cellfun(@(x) x(n+1:n+15), {l.name}, 'UniformOutput', false), ...
        'InputFormat', 'yyyyMMdd_HHmmss'), 'Start', 'Day');
      % Find adjacent date before the asked date (within 10 days)
      i = 1;
      % Select previous valid day (within 10 days) and asked dt
      % while ~any(l_dt == dt - i) && i < 10; i = i + 1; end
      % sel = l_dt == dt - i | l_dt == dt;
      while ~any(l_dt == dt - days(i)) && i < 10; i = i + 1; end
      sel = l_dt == dt - days(i) | l_dt == dt;
      % Return selected filenames
      filenames = {l(sel).name}';
    otherwise
      error('Software not supported: %s.', software);
  end
  if ~isempty(filenames)
    if ~iscell(filenames); filenames = cellstr(filenames); end
  end
end

function str = dt_yyyymmdd(dt)
  % str = datestr(dt, 'yyyymmdd');
  str = char(datetime(dt, 'Format', 'yyyyMMdd'));
end

function str = dt_yymmdd(dt)
  % str = datestr(dt, 'yymmdd');
  str = char(datetime(dt, 'Format', 'yyMMdd'));
end

function str = dt_yyyydoy(dt)
%   str = sprintf('%d%03d',year(dt),datevec2doy(datevec(dt)));
  dtvec = datevec(dt);
  str = sprintf('%d%03d', dtvec(1), datevec2doy(dtvec));
end

function [str] = dt_doy(dt)
  str = sprintf('%03d', datevec2doy(datevec(dt)));
end

function merged_data = check_nb_variables(gdata, dt)
  all_table = table();
  all_table.sz = cell2mat(cellfun(@(x) size(x,2), gdata, 'un', 0));
  gdata(all_table.sz == 0) = [];
  all_table(all_table.sz == 0, :) = [];
  all_table.varnames = cellfun(@(c) c.Properties.VariableNames, gdata, 'un', 0);
  all_table.var_combi = cellfun(@(x) cell2mat(join(x, ',')), all_table.varnames, 'un', 0);
  % detect all unique combinaisons of variables
  foo_var = table();
  foo_var.uvc = unique(all_table.var_combi);
  if size(foo_var.uvc,1) > 1
    % count occurence of each combinason
    foo_var.uvc_ct = NaN(size(foo_var.uvc));
    for i = 1:size(foo_var.uvc,1)
      foo_var.uvc_ct(i) = sum(strcmp(all_table.var_combi, foo_var.uvc(i)));
    end
    % find reference table
    varmax = all_table.sz == max(all_table.sz);
    if size(all_table, 1) == 2
      foosz = cell2mat(cellfun(@(x) size(x,1), gdata, 'un', 0));
      tabref = find(varmax & foosz == max(foosz(varmax)), 1, 'first');
    else
      usel_varmax = logical(sum(categorical(all_table.var_combi(varmax))' == categorical(foo_var.uvc), 2));
      id_max = find(usel_varmax & foo_var.uvc_ct == max(foo_var.uvc_ct(usel_varmax)), 1, 'first');
      tabref = find(varmax & strcmp(all_table.var_combi, foo_var.uvc(id_max)), 1, 'first');
    end
    % re-order all table variable names to match the max number of variable and the order as the most frequent combinaison
    for i = 1:size(gdata,1)
      if ~strcmp(all_table.var_combi{i}, all_table.var_combi(tabref))
        missing_var = all_table.varnames{tabref}(~ismember(all_table.varnames{tabref}, all_table.varnames{i}));
        if size(missing_var, 2) == all_table.sz(i) - 1
          error('Impossible to consolidating files with different variable names: %s', dt)
        elseif size(missing_var, 2) > 0
          warning('Consolidating files with different number of variables. %s missing in %s files: empty variable added', ...
            ['"' cell2mat(join(missing_var, '", "')) '"'], dt)
        else
          warning('Consolidating %s files with variables in different order: variables re-ordered', ...
            ['"' cell2mat(join(missing_var, '", "')) '"'], dt)
        end
        foo = gdata{i};
        for j = missing_var
          if isnumeric(gdata{tabref}.(j{:}))
            foo.(j{:}) = NaN(size(gdata{i}.dt));
          elseif islogical(gdata{tabref}.(j{:}))
            foo.(j{:}) = false(size(gdata{i}.dt));
          elseif iscell(gdata{tabref}.(j{:}))
            foo.(j{:}) = cell(size(gdata{i}.dt));
          elseif isdatetime(gdata{tabref}.(j{:}))
            foo.(j{:}) = NaT(size(gdata{i}.dt));
          end
        end
        % force order
        [~,Y] = ismember(foo.Properties.VariableNames, all_table.varnames{tabref});
        [~,ord] = sort(Y);
        gdata{i} = foo(:,ord);
      end
    end
  end
  merged_data = vertcat(gdata{:});
% function [gdata, data] = check_nb_variables(gdata, data, dt)
%   if all(size(gdata, 2) ~= size(data, 2) & ~isempty(gdata) & ~isempty(data))
%     if all(ismember(gdata.Properties.VariableNames, data.Properties.VariableNames)) && ...
%         ~all(ismember(data.Properties.VariableNames, gdata.Properties.VariableNames))
%       data = data(:, ismember(data.Properties.VariableNames, gdata.Properties.VariableNames));
%       missing_var = data.Properties.VariableNames(~ismember(data.Properties.VariableNames, gdata.Properties.VariableNames));
%       before_or_after = 'before';
%     elseif ~all(ismember(gdata.Properties.VariableNames, data.Properties.VariableNames)) && ...
%         all(ismember(data.Properties.VariableNames, gdata.Properties.VariableNames))
%       gdata = gdata(:, ismember(gdata.Properties.VariableNames, data.Properties.VariableNames));
%       missing_var = data.Properties.VariableNames(~ismember(gdata.Properties.VariableNames, data.Properties.VariableNames));
%       before_or_after = 'after';
%     else
%       error('Impossible to consolidate table with different variable names, check files before and after %s', dt)
%     end
%     warning('Consolidating files with different number of variables. %s missing in files %s %s: variable ignored', ...
%       ['"' cell2mat(join(missing_var, '", "')) '"'], before_or_after, dt)
%   end
end


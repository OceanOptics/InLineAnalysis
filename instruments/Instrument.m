classdef (Abstract) Instrument < handle
  %INSTRUMENTS Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    % Instrument caracteristics
    model = '';
    sn = '';
    logger = '';
    prefix = '';
    di_cfg = struct('logger', '', 'prefix', '', 'postfix', '');
    
    % Path to data
    path = struct('raw', '', 'di', '', 'wk', '', 'prod', '', 'ui', '');
    
    % Instrument-specific processing parameters
    split = struct('mode', 'split');
    bin_method = 'SB_ALL';
    
    % Instrument-specific view parametes
    view = struct('varname', '', 'varcol', 1);
    
    % Instrument data
    %%% NOTE: For DIW QC is done before the Binning %%%
    % Raw & sync
    data = table();
    % Split
    raw = struct('tsw', table(), 'fsw', table(), 'bad', table(), 'diw', table());
    % Bin
    bin = struct('tsw', table(), 'fsw', table(), 'diw', table());
    % QC
    qc = struct('tsw', table(), 'fsw', table(), 'diw', table());
    suspect = struct('tsw', table(), 'fsw', table(), 'diw', table());
    bad = struct('tsw', table(), 'fsw', table(), 'diw', table());
    % Calibrated
    prod = struct();
  end
  properties (SetAccess = private, GetAccess = public)
    % Processing parameters modiwfied by instrument methods only
    sync_delay = 0; % days
    stretch_delta = 0; % days (can also be array of size of dt)
  end
  
  methods
    % Constructor
    function obj = Instrument(cfg, ~)
      
      % Object Initilization
      % obj = obj@handle();
      
      % Post initialization
      if nargin ~= 0
        % Load required fields
        if isfield(cfg, 'model'); obj.model = cfg.model;
        else; error('Missing field model.'); end
        if isfield(cfg, 'path') && isfield(cfg.path, 'raw') &&...
          isfield(cfg.path, 'wk') && isfield(cfg.path, 'prod') && isfield(cfg.path, 'ui')
          for f=fieldnames(cfg.path)'; obj.path.(f{1}) = cfg.path.(f{1}); end
        else
          error('Missing field path, path.raw, path.wk, path.prod, and/or path.ui.');
        end
        
        % Load optional fields
        if isfield(cfg, 'sn'); obj.sn = cfg.sn; end
        if isfield(cfg, 'logger'); obj.logger = cfg.logger; end
        if isfield(cfg, 'ila_prefix'); obj.prefix = cfg.ila_prefix;
        else; obj.prefix = [obj.model obj.sn '_']; end
        if isempty(obj.prefix); obj.prefix = [obj.model obj.sn '_']; end

        if isfield(cfg, 'view')
          if isfield(cfg.view, 'varname'); obj.view.varname = cfg.view.varname; end
          if isfield(cfg.view, 'varcol'); obj.view.varcol = cfg.view.varcol; end
        end
        if isfield(cfg, 'di')
          if isfield(cfg.di, 'prefix'); obj.di_cfg.prefix = cfg.di.prefix; end
          if isfield(cfg.di, 'prefix'); obj.di_cfg.prefix = cfg.di.prefix; end
          if isfield(cfg.di, 'postfix'); obj.di_cfg.postfix = cfg.di.postfix; end
%           if isfield(cfg.di, 'logger'); obj.di_cfg.logger = cfg.di.logger; end
          obj.di_cfg.logger = cfg.logger;
        end
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PRE-PROCESSING METHODS %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Sync(obj, delay_in_seconds)
      % SYNC synchronize time of instrument with other insrtruments
      %   Note: if function is run several time the delay applied is the
      %         always in reference to the original data (not the one
      %         already synchronized).
      
      % Convert in delay in day
      delay_in_days = delay_in_seconds / 3600 / 24;
      % Apply synchronisation
      obj.data.dt = obj.data.dt + obj.sync_delay - delay_in_days;
      % Update delay stored in instrument properties
      obj.sync_delay = delay_in_days;
    end
    
    function Stretch(obj, delta_in_seconds)
      % STRETCH stretch the time of instrument
      %   this function is mainly intendended for instrument recorded with
      %   the DH4 which clock shift if in time about 3 min/ 30 days
      
      % Convert delay in day
      delta_in_days = delta_in_seconds / 3600 / 24;
      % Linearly interpolate delay
      %   no delay at the beginning, full delay at the end.
      delta_interp = interp1([obj.data.dt(1) obj.data.dt(end)],...
                             [0 delta_in_days],obj.data.dt, 'linear');
      % Apply Stretch
      obj.data.dt = obj.data.dt + obj.stretch_delta - delta_interp;
      % Update delta stored in instrument properties
      obj.stretch_delta = delta_interp;
    end
    
    function Split(obj, ref, buffer)
      reference_data = ref.data;
      reference_constants = struct('SWITCH_FILTERED', ref.SWITCH_FILTERED, 'SWITCH_TOTAL', ref.SWITCH_TOTAL);
      [obj.raw.tsw, obj.raw.fsw, obj.raw.bad] = splitTable(reference_data,...
        obj.data, buffer, obj.split.mode, reference_constants, false);
    end
    
%     function Bin(obj, bin_size_minutes, method, prctile_detection, prctile_average, parallel)
%       bin_size_days = bin_size_minutes / 60 / 24;
%       if isempty(obj.raw.tsw)
%         fprintf('WARNING: No raw.tsw data to bin\n');
%       else
%         obj.bin.tsw = binTable(obj.raw.tsw, bin_size_days, method, prctile_detection, prctile_average, false, parallel, false);
%       end
%       if isempty(obj.raw.fsw)
%         if ~(strcmp(obj.split.mode, 'rmBuffer') || strcmp(obj.split.mode, 'None'))
%           fprintf('WARNING: No raw.fsw data to bin\n');
%         end
%       else
%         obj.bin.fsw = binTable(obj.raw.fsw, bin_size_days, method, prctile_detection, prctile_average, false, parallel, false);
%       end
%     end
    
    function Bin(obj, bin_size_minutes, prctile_detection, prctile_average, parallel, mode)
      bin_size_days = bin_size_minutes / 60 / 24;
      if isempty(obj.raw.tsw)
        fprintf('WARNING: No raw.tsw data to bin\n');
      else
        fprintf('\tTSW\n');
        switch mode
          case 'OneShot'
            obj.bin.tsw = binTable(obj.raw.tsw, bin_size_days, obj.bin_method, ...
              prctile_detection, prctile_average, false, parallel, false);
          case 'ByDay'
            for d = floor(min(obj.raw.tsw.dt)):floor(max(obj.raw.tsw.dt))
              fprintf('\t\t%s', datestr(d)); tic;
              sel = d <= obj.raw.tsw.dt & obj.raw.tsw.dt < d + 1;
              if sum(sel) == 0
                fprintf('  No total data to bin\n');
                continue
              end
              obj.bin.tsw = [obj.bin.tsw; binTable(obj.raw.tsw(sel,:), bin_size_days, ...
                obj.bin_method, prctile_detection, prctile_average, false, parallel, false)];
              t = toc; fprintf('  %1.3f s\n', t);
            end
          otherwise
            error('Binning mode not supported.');
        end
      end
      if isempty(obj.raw.fsw)
        if ~(strcmp(obj.split.mode, 'rmBuffer') || strcmp(obj.split.mode, 'None'))
          fprintf('WARNING: No raw.fsw data to bin\n');
        end
      else
        fprintf('\tFSW\n');
        switch mode
          case 'OneShot'
            obj.bin.fsw = binTable(obj.raw.fsw, bin_size_days/3, obj.bin_method, ...
              prctile_detection, prctile_average, false, parallel, false);
          case 'ByDay'
            for d = floor(min(obj.raw.tsw.dt)):floor(max(obj.raw.tsw.dt))
              fprintf('\t\t%s', datestr(d)); tic;
              sel = d <= obj.raw.fsw.dt & obj.raw.fsw.dt < d + 1;
              if sum(sel) == 0
                fprintf('  No filtered data to bin\n');
                continue
              end
              obj.bin.fsw = [obj.bin.fsw; binTable(obj.raw.fsw(sel,:), bin_size_days/3, ...
                obj.bin_method, prctile_detection, prctile_average, false, parallel, false)];
              t = toc; fprintf('  %1.3f s\n', t);
            end
          otherwise
            error('Binning mode not supported.');
        end
      end
    end
    
    function BinDI(obj, bin_size_minutes, prctile_detection, prctile_average, parallel)
      %%% NOTE: For DIW QC is done before the Binning %%%
%       bin_size_minutes = 60;
      % BinDI is only in mode OneShot, no mode day by day (as in a typical setup they won't many samples
      bin_size_days = bin_size_minutes / 60 / 24;
      if isempty(obj.qc.diw)
        fprintf('WARNING: No qc.diw data to bin\n');
      else
        obj.bin.diw = binTable(obj.qc.diw, bin_size_days, obj.bin_method, prctile_detection, ....
          prctile_average, true, parallel, false);
      end
    end
    
    function Flag(obj, params)
      % Make flags
      flags_tsw = flagTable(obj.bin.tsw, params.tot);
      % Make flag selection
      [sel_good_tsw, sel_suspect_tsw] = selFlag(flags_tsw, params.tot.primary_varname, ...
        params.tot.min_flag_n);
      % Keep data
      obj.qc.tsw = obj.bin.tsw(sel_good_tsw,:);
      obj.suspect.tsw = obj.bin.tsw(sel_suspect_tsw,:);
      % Same for filtered sea water (fsw)
      if ~isempty(obj.bin.fsw)
        flags_fsw = flagTable(obj.bin.fsw, params.filt);
        [sel_good_fsw, sel_suspect_fsw] = selFlag(flags_fsw, params.filt.primary_varname, ...
          params.filt.min_flag_n);
        obj.suspect.fsw = obj.bin.fsw(sel_suspect_fsw,:);
        obj.qc.fsw = obj.bin.fsw(sel_good_fsw,:);
      end
    end
    
    function DeleteUserSelection(obj, user_selection, level, chan)
      if nargin < 3
        level = 'qc';
        chan = {'all'};
      elseif nargin < 4
        chan{2} = {'all'};
      end
      if size(chan, 2) < 2 || any(cellfun('isempty', chan))
        chan = {'all'};
      end
      if size(user_selection, 2) == 2
        if any(strcmp(chan{1}, 'all'))
          fieldn = fieldnames(obj.(level))';
          for j = fieldn; j = j{1};
            if ~isempty(obj.(level).(j))
              for i=1:size(user_selection, 1)
                obj.(level).(j)(user_selection(i,1) <= obj.(level).(j).dt & ...
                  obj.(level).(j).dt <= user_selection(i,2), :) = []; 
              end
            end
          end
        elseif any(strcmp(chan{2}, 'all'))
          if ~isempty(obj.(level).(chan{1}))
            for i=1:size(user_selection, 1)
              obj.(level).(chan{1})(user_selection(i,1) <= obj.(level).(chan{1}).dt & ...
                obj.(level).(chan{1}).dt <=  user_selection(i,2), :) = [];
            end
          end
        else
          if ~isempty(obj.(level).(chan{1}))
            for i=1:size(user_selection, 1)
              obj.(level).(chan{1}).(chan{2})(user_selection(i,1) <= obj.(level).(chan{1}).dt & ...
                obj.(level).(chan{1}).dt <=  user_selection(i,2), :) = NaN;
            end
          end
        end
      elseif size(user_selection, 2) == 1
        if any(strcmp(chan{1}, 'all'))
          fieldn = fieldnames(obj.(level))';
          for j = fieldn; j = j{1};
            if ~isempty(obj.(level).(j))
              obj.(level).(j)(ismember(round(obj.(level).(j).dt, 9), ...
                round(user_selection, 9)), :) = [];
            end
          end
        elseif any(strcmp(chan{2}, 'all'))
          if isfield(obj.(level), chan{1})
            if ~isempty(obj.(level).(chan{1}))
              obj.(level).(chan{1})(ismember(round(obj.(level).(chan{1}).dt, 9), ...
                round(user_selection, 9)), :) = [];
            end
          end
        else
          if isfield(obj.(level), chan{1})
            if ~isempty(obj.(level).(chan{1}))
              obj.(level).(chan{1}).(chan{2})(ismember(round(obj.(level).(chan{1}).dt, 9), ...
                round(user_selection, 9)), :) = NaN;
            end
          end
        end
      end
    end
    
    function Write(obj, filename_prefix, days2write, level, part_or_diw)
      if nargin < 4
        level = 'prod';
        part_or_diw = 'all';
      elseif nargin < 5
        part_or_diw = 'all';
      end
      part_or_diw = lower(part_or_diw);
      fieldna = fieldnames(obj.(level))';
      switch part_or_diw
        case {'part', 'particulate'}
          fieldna(contains(fieldna, {'diw','g'})) = [];
        case {'diw', 'dissolved'}
          fieldna(contains(fieldna, {'tsw','fsw','p','FiltStat','QCfailed'})) = [];
        otherwise
          error("Type of data to write not recognized: 'part' or 'diw' or 'all'")
      end
%       if isstruct(obj.(level))
        % For each product type (particulate, dissoved...)
        for f = fieldna; f = f{1};
          filename = [filename_prefix '_' level '_' f '.mat'];
          if isempty(obj.(level).(f)); continue; end
          days2write = floor(days2write); % force days2write to entire day
          sel = min(days2write) <= obj.(level).(f).dt & obj.(level).(f).dt < max(days2write) + 1;
          if ~any(sel); fprintf('WRITE: %s_%s_%s No data.\n', filename_prefix, level, f); continue; end
          data = obj.(level).(f)(sel,:);
          if ~isfolder(obj.path.wk); mkdir(obj.path.(level)); end
          tic
          save(fullfile(obj.path.wk, filename), 'data');
          fprintf('WRITE: %s saved', filename);
          t = toc; fprintf('  %1.3f s\n', t);
        end
%       else
%         % One Table at the level
%         filename = [filename_prefix '_' level '.mat'];
%         if isempty(obj.(level)); fprintf('WRITE: %s_%s No data.\n', filename_prefix, level); return; end
%         days2write = floor(days2write); % force days2write to entire day
%         sel = min(days2write) <= obj.(level).dt & obj.(level).dt < max(days2write) + 1;
%         if ~any(sel); fprintf('WRITE: %s_%s No data.\n', filename_prefix, level); return; end
%         data = obj.(level)(sel,:);
%         if ~isfolder(obj.path.wk); mkdir(obj.path.(level)); end
%         save(fullfile(obj.path.wk, filename), 'data');
%       end
    end
    
    function Read(obj, filename_prefix, days2read, level)
      % For each product type (particulate, dissoved...)
      % This will simply add data at the end of the current table
      %   (if data was already in memory it could duplicate timestamps)
      if contains(obj.model, 'AC')
        obj.ReadDeviceFile()
      end
      if nargin < 4; level = 'prod'; end
      if strcmp(level, 'data')
        l = dir(fullfile(obj.path.wk, [filename_prefix '.mat']));
      else
        l = dir(fullfile(obj.path.wk, [filename_prefix '_' level '*.mat']));
      end
      if isempty(l)
        fprintf('%s: %s_%s No data.\n', datestr(days2read), filename_prefix, level);
      else
        for f = {l.name}; f = f{1};
          fprintf('\t\t%s', f); tic;
          load(fullfile(obj.path.wk, f), 'data'); % data variable is created
          data_temp = sortrows(data, 1);
          % add column of spd to make sure all files have equal number of colums
          if strcmp(obj.model, 'FTH')
            foospd = {'spd1', 'spd2'};
            if strcmp(level, 'raw')
              if size(data_temp, 2) < 4
                data_temp = renamevars(data_temp, 'spd', obj.view.spd_variable);
                data_temp.(foospd{~strcmp(foospd, obj.view.spd_variable)}) = ...
                  NaN(size(data_temp.(obj.view.spd_variable)));
                data_temp = movevars(data_temp,'spd2','After','spd1');
              end
            else
              if size(data_temp, 2) < 7
                data_temp = renamevars(data_temp, 'spd', obj.view.spd_variable);
                data_temp = renamevars(data_temp, 'spd_avg_sd', [obj.view.spd_variable '_avg_sd']);
                data_temp = renamevars(data_temp, 'spd_avg_n', [obj.view.spd_variable '_avg_n']);
                data_temp.(foospd{~strcmp(foospd, obj.view.spd_variable)}) = ...
                  NaN(size(data_temp.(obj.view.spd_variable)));
                data_temp.([foospd{~strcmp(foospd, obj.view.spd_variable)} '_avg_sd']) = ...
                  NaN(size(data_temp.(obj.view.spd_variable)));
                data_temp.([foospd{~strcmp(foospd, obj.view.spd_variable)} '_avg_n']) = ...
                  NaN(size(data_temp.(obj.view.spd_variable)));
                data_temp = movevars(data_temp,'spd2','After','spd1_avg_n');
                data_temp = movevars(data_temp,'spd2_avg_sd','After','spd2');
                data_temp = movevars(data_temp,'spd2_avg_n','After','spd2_avg_sd');
              end
            end
          end
          if isempty(data_temp)
            warning('%s is empty, the file was deleted', f)
            delete(fullfile(obj.path.wk, f))
          else
            sel = min(days2read) <= data_temp.dt & data_temp.dt < max(days2read) + 1;
            fn = strsplit(f, {'_','.'}); fn = fn{end-1};%(1:end-4);
            if strcmp(fn, level)
              if ~isempty(obj.(level))
                [obj.(level), data_temp] = check_nb_variables(obj.(level), data_temp(sel,:), f);
              end
              obj.(level)(end+1:end+sum(sel),:) = data_temp(sel,:);
              obj.(level).Properties.CustomProperties = data_temp.Properties.CustomProperties;
            else
              if isfield(obj.(level), fn)
                if ~isempty(obj.(level).(fn))
                  [obj.(level).(fn), data_temp] = check_nb_variables(obj.(level).(fn), data_temp(sel,:), f);
                end
                obj.(level).(fn)(end+1:end+sum(sel),:) = data_temp(sel,:);
                obj.(level).(fn).Properties.CustomProperties = data_temp.Properties.CustomProperties;
              else
                if strcmp(level, 'prod')
                  obj.(level).(fn) = data_temp(sel,:);
                  obj.(level).(fn).Properties.CustomProperties = data_temp.Properties.CustomProperties;
                else
                  if ~isempty(obj.(level))
                    [obj.(level), data_temp] = check_nb_variables(obj.(level), data_temp(sel,:), f);
                  end
                  obj.(level)(end+1:end+sum(sel),:) = data_temp(sel,:);
                  obj.(level).Properties.CustomProperties = data_temp.Properties.CustomProperties;
                end
              end
            end
            t = toc; fprintf('  %1.3f s\n', t);
          end
        end
      end
    end
    
%     function Write(obj, filename_prefix, days2write)
%       % For each product type (particulate, dissoved...)
%       for f = fieldnames(obj.prod)'; f = f{1};
%         filename = [filename_prefix '_' f '_prod.mat'];
%         sel = min(days2write) <= obj.prod.(f).dt & obj.prod.(f).dt < max(days2write) + 1;
%         data = obj.prod.(f)(sel,:);
%         if ~isdir(obj.path.prod); mkdir(obj.path.prod); end
%         save(fullfile(obj.path.prod, filename), 'data');
%       end
%     end
    
%     function LoadProducts(obj, filename_prefix, days2read)
%       % For each product type (particulate, dissoved...)
%       % This will simply add data at the end of the current table
%       %   (if data was already in memory it could duplicate timestamps)
%       l = dir(fullfile(obj.path.prod, [filename_prefix '_*_prod.mat']);
%       if isempty(l)
%         fprintf('%s: %s No data.\n', datestr(days2read), filename_prefix);
%       else
%         for f = {l.name}; f = f{1};
%           load(fullfile(obj.path.prod, f), 'data'); % data variable is created
%           sel = min(days2read) <= data.dt & data.dt < max(days2read) + 1;
%           fn = strsplit(f, '_'); fn = fn{end-1};
%           if isfield(obj.prod, fn)
%             obj.prod.(fn)(end+1:end+sum(sel),:) = data(sel,:);
%           else
%             obj.prod.(fn) = data(sel,:);
%           end
%         end
%       end
%     end
    
  end
  
%   methods (Abstract=true)
%     % Methods without implementation. Subclasses will define these methods
%     function ReadRaw(obj, days2run, force_import, write)      
%       % READRAW import | load raw data
% %       error('Function not implemented.\n');
%     end
%     
%     function ApplyUserInput(obj, user_selection, mode)
%       % Only used by FTH for now
%     end
%   end
  
end

function [gdata, data] = check_nb_variables(gdata, data, dt)
  if all(size(gdata, 2) ~= size(data, 2) & ~isempty(gdata) & ~isempty(data))
    if all(ismember(gdata.Properties.VariableNames, data.Properties.VariableNames)) && ...
        ~all(ismember(data.Properties.VariableNames, gdata.Properties.VariableNames))
      missing_var = data.Properties.VariableNames(~ismember(data.Properties.VariableNames, gdata.Properties.VariableNames));
      data = data(:, ismember(data.Properties.VariableNames, gdata.Properties.VariableNames));
      before_or_after = 'before';
    elseif ~all(ismember(gdata.Properties.VariableNames, data.Properties.VariableNames)) && ...
        all(ismember(data.Properties.VariableNames, gdata.Properties.VariableNames))
      missing_var = data.Properties.VariableNames(~ismember(gdata.Properties.VariableNames, data.Properties.VariableNames));
      gdata = gdata(:, ismember(gdata.Properties.VariableNames, data.Properties.VariableNames));
      before_or_after = 'after';
    end
    warning('Consolidating files with different number of variables. %s missing in files %s %s: variable ignored', ...
      ['"' cell2mat(join(missing_var, '", "')) '"'], before_or_after, dt)
  end
end


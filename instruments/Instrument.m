classdef (Abstract) Instrument < handle
  %INSTRUMENTS Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    % Instrument caracteristics
    model = '';
    sn = '';
    logger = '';
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
    function obj = Instrument(cfg, instrument_id)
      
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
        if isfield(cfg, 'view')
          if isfield(cfg.view, 'varname'); obj.view.varname = cfg.view.varname; end
          if isfield(cfg.view, 'varcol'); obj.view.varcol = cfg.view.varcol; end
        end
        if isfield(cfg, 'di')
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
            obj.bin.tsw = binTable(obj.raw.tsw, bin_size_days, obj.bin_method, prctile_detection, prctile_average, false, parallel, false);
          case 'ByDay'
            for d = floor(min(obj.raw.tsw.dt)):floor(max(obj.raw.tsw.dt))
              fprintf('\t\t%s', datestr(d)); tic;
              sel = d <= obj.raw.tsw.dt & obj.raw.tsw.dt < d + 1;
              if sum(sel) == 0
                fprintf('  No total data to bin\n');
                continue
              end
              obj.bin.tsw = [obj.bin.tsw; binTable(obj.raw.tsw(sel,:), bin_size_days, obj.bin_method, prctile_detection, prctile_average, false, parallel, false)];
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
            obj.bin.fsw = binTable(obj.raw.fsw, bin_size_days, obj.bin_method, prctile_detection, prctile_average, false, parallel, false);
          case 'ByDay'
            for d = floor(min(obj.raw.tsw.dt)):floor(max(obj.raw.tsw.dt))
              fprintf('\t\t%s', datestr(d)); tic;
              sel = d <= obj.raw.fsw.dt & obj.raw.fsw.dt < d + 1;
              if sum(sel) == 0
                fprintf('  No filtered data to bin\n');
                continue
              end
              obj.bin.fsw = [obj.bin.fsw; binTable(obj.raw.fsw(sel,:), bin_size_days, obj.bin_method, prctile_detection, prctile_average, false, parallel, false)];
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
        obj.bin.diw = binTable(obj.qc.diw, bin_size_days, obj.bin_method, prctile_detection, prctile_average, true, parallel, false);
      end
    end
    
    function Flag(obj, params)
      % Make flags
      flags_tsw = flagTable(obj.bin.tsw, params.tot);
      % Make flag selection
      [sel_good_tsw, sel_suspect_tsw] = selFlag(flags_tsw, params.tot.primary_varname, params.tot.min_flag_n);
      % Keep data
      obj.qc.tsw = obj.bin.tsw(sel_good_tsw,:);
      obj.suspect.tsw = obj.bin.tsw(sel_suspect_tsw,:);
      % Same for filtered sea water (fsw)
      if ~isempty(obj.bin.fsw)
        flags_fsw = flagTable(obj.bin.fsw, params.filt);
        [sel_good_fsw, sel_suspect_fsw] = selFlag(flags_fsw, params.filt.primary_varname, params.filt.min_flag_n);
        obj.suspect.fsw = obj.bin.fsw(sel_suspect_fsw,:);
        obj.qc.fsw = obj.bin.fsw(sel_good_fsw,:);
      end
    end
    
    function DeleteUserSelection(obj, user_selection)
      for i=1:size(user_selection, 1)
%         obj.bad.tsw = obj.qc.tsw(user_selection(i,1) <= obj.qc.tsw.dt & obj.qc.tsw.dt <= user_selection(i,2),:);
%         obj.bad.fsw(user_selection(i,1) <= obj.qc.fsw.dt & obj.qc.fsw.dt <= user_selection(i,2),:);
        if ~isempty(obj.qc.tsw)
          obj.qc.tsw(user_selection(i,1) <= obj.qc.tsw.dt & obj.qc.tsw.dt <= user_selection(i,2),:) = []; 
        end
        if ~isempty(obj.qc.fsw)
          obj.qc.fsw(user_selection(i,1) <= obj.qc.fsw.dt & obj.qc.fsw.dt <= user_selection(i,2),:) = []; 
        end
        if ~isempty(obj.qc.diw)
          obj.qc.diw(user_selection(i,1) <= obj.qc.diw.dt & obj.qc.diw.dt <= user_selection(i,2),:) = []; 
        end
      end
    end
    
    function Write(obj, filename_prefix, days2write, level)
      if nargin < 4; level = 'prod'; end
      if isstruct(obj.(level))
        % For each product type (particulate, dissoved...)
        for f = fieldnames(obj.(level))'; f = f{1};
          filename = [filename_prefix '_' level '_' f '.mat'];
          if isempty(obj.(level).(f)); continue; end
          days2write = floor(days2write); % force days2write to entire day
          sel = min(days2write) <= obj.(level).(f).dt & obj.(level).(f).dt < max(days2write) + 1;
          if ~any(sel); fprintf('WRITE: %s_%s_%s No data.\n', filename_prefix, level, f); continue; end
          data = obj.(level).(f)(sel,:);
          if ~isdir(obj.path.wk); mkdir(obj.path.(level)); end
          save([obj.path.wk filename], 'data');
        end
      else
        % One Table at the level
        filename = [filename_prefix '_' level '.mat'];
        if isempty(obj.(level)); fprintf('WRITE: %s_%s No data.\n', filename_prefix, level); return; end
        days2write = floor(days2write); % force days2write to entire day
        sel = min(days2write) <= obj.(level).dt & obj.(level).dt < max(days2write) + 1;
        if ~any(sel); fprintf('WRITE: %s_%s No data.\n', filename_prefix, level); return; end
        data = obj.(level)(sel,:);
        if ~isdir(obj.path.wk); mkdir(obj.path.(level)); end
        save([obj.path.wk filename], 'data');
      end
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
        l = dir([obj.path.wk filename_prefix '.mat']);
      else
        l = dir([obj.path.wk filename_prefix '_' level '*.mat']);
      end
      if isempty(l)
        fprintf('%s: %s_%s No data.\n', datestr(days2read), filename_prefix, level);
      else
        for f = {l.name}; f = f{1};
          fprintf('\t\t%s', f); tic;
          load([obj.path.wk f], 'data'); % data variable is created
          if isempty(data)
            warning('%s is empty, the file was deleted', f)
            delete([obj.path.wk f])
          else
            sel = min(days2read) <= data.dt & data.dt < max(days2read) + 1;
            fn = strsplit(f, {'_','.'}); fn = fn{end-1};%(1:end-4);
            if strcmp(fn, level)
              obj.(level)(end+1:end+sum(sel),:) = data(sel,:);
            else
              if isfield(obj.(level), fn)
                obj.(level).(fn)(end+1:end+sum(sel),:) = data(sel,:);
              else
                if strcmp(level, 'prod')
                  obj.(level).(fn) = data(sel,:);
                else
                  obj.(level)(end+1:end+sum(sel),:) = data(sel,:);
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
%         save([obj.path.prod filename], 'data');
%       end
%     end
    
%     function LoadProducts(obj, filename_prefix, days2read)
%       % For each product type (particulate, dissoved...)
%       % This will simply add data at the end of the current table
%       %   (if data was already in memory it could duplicate timestamps)
%       l = dir([obj.path.prod filename_prefix '_*_prod.mat']);
%       if isempty(l)
%         fprintf('%s: %s No data.\n', datestr(days2read), filename_prefix);
%       else
%         for f = {l.name}; f = f{1};
%           load([obj.path.prod f], 'data'); % data variable is created
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


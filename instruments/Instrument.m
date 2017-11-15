classdef (Abstract) Instrument < handle
  %INSTRUMENTS Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    % Instrument caracteristics
    model = '';
    sn = '';
    logger = '';
    
    % Path to data
    path = struct('raw', '', 'wk', '', 'prod', '');
    
    % Instrument-specific processing parameters
    split = struct('mode', 'split');
    
    % Instrument-specific view parametes
    view = struct('varname', '', 'varcol', 1);
    
    % Instrument data
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
  end
  
  methods
    % Constructor
    function obj = Instrument(cfg)
      
      % Object Initilization
      % obj = obj@handle();
      
      % Post initialization
      if nargin ~= 0
        % Load required fields
        if isfield(cfg, 'model'); obj.model = cfg.model;
        else; error('Missing field model.'); end
        if isfield(cfg, 'path') && isfield(cfg.path, 'raw') &&...
          isfield(cfg.path, 'wk') && isfield(cfg.path, 'prod')
          obj.path = cfg.path;
        else
          error('Missing field path, path.raw, path.wk, or path.prod.');
        end
        
        % Load optional fields
        if isfield(cfg, 'sn'); obj.sn = cfg.sn; end
        if isfield(cfg, 'logger'); obj.logger = cfg.logger; end
        if isfield(cfg, 'view'); 
          if isfield(cfg.view, 'varname'); obj.view.varname = cfg.view.varname; end
          if isfield(cfg.view, 'varcol'); obj.view.varcol = cfg.view.varcol; end
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
    
    function Split(obj, ref, buffer)
      reference_data = ref.data;
      reference_constants = struct('SWITCH_FILTERED', ref.SWITCH_FILTERED, 'SWITCH_TOTAL', ref.SWITCH_TOTAL);
      [obj.raw.tsw, obj.raw.fsw, obj.raw.bad] = splitTable(reference_data,...
        obj.data, buffer, obj.split.mode, reference_constants, false);
    end
    
    function Bin(obj, bin_size_minutes, prctile_detection, prctile_average, parallel)
      bin_size_days = bin_size_minutes / 60 / 24;
      if isempty(obj.raw.tsw)
        fprintf('WARNING: No raw.tsw data to bin\n');
      else
        obj.bin.tsw = binTable(obj.raw.tsw, bin_size_days, '4flag', prctile_detection, prctile_average, parallel, false);
      end
      if isempty(obj.raw.fsw)
        if ~strcmp(obj.split.mode, 'rmBuffer')
          fprintf('WARNING: No raw.fsw data to bin\n');
        end
      else
        obj.bin.fsw = binTable(obj.raw.fsw, bin_size_days, '4flag', prctile_detection, prctile_average, parallel, false);
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
        obj.qc.tsw(user_selection(i,1) <= obj.qc.tsw.dt & obj.qc.tsw.dt <= user_selection(i,2),:) = []; 
        if ~isempty(obj.qc.fsw)
          obj.qc.fsw(user_selection(i,1) <= obj.qc.fsw.dt & obj.qc.fsw.dt <= user_selection(i,2),:) = []; 
        end
      end
    end
    
    function Write(obj, filename_prefix, days2write)
      % For each product type (particulate, dissoved...)
      for f = fieldnames(obj.prod); f = f{1};
        filename = [filename_prefix '_' f '_prod.mat'];
        sel = min(days2write) <= obj.prod.(f).dt & obj.prod.(f).dt < max(days2write) + 1;
        data = obj.prod.(f)(sel,:);
        save([obj.path.prod filename], 'data');
      end
    end
    
    function LoadProducts(obj, filename_prefix, days2read)
      % For each product type (particulate, dissoved...)
      l = dir([obj.path.prod filename_prefix '_*_prod.mat']);
      for f = {l.name}'; f = f{1};
        load([obj.path.prod f], 'data'); % data variable is created
        sel = min(days2read) <= data.dt & data.dt < max(days2read) + 1;
        fn = strsplit(f, '_'); fn = fn{end-1};
        if isfield(obj.prod, fn)
          obj.prod.(fn)(end+1:end+sum(sel),:) = data(sel,:);
        else
          obj.prod.(fn) = data(sel,:);
        end
      end
    end
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


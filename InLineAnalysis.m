classdef InLineAnalysis < handle
  %INLINEANALYSIS Summary of this class goes here
  %   Detailed explanation goes here
  
 % The handle class as representing objects whose identity is independent of the values of their properties
  
  properties
    meta=struct();
    cfg=struct();
    instrument = struct();
  end
  
  methods
    % Constructor
    function obj = InLineAnalysis(cfg_file_name)
      % Pre-initialization
      % Add path to library of functions
      addpath('lib', 'instruments', 'packages', 'packages/datetick2_doy', 'packages/spectral_color');
      
      % Object Initilization
      obj = obj@handle();
      
      % Post initialization
      if nargin ~= 0
        % Load configuration
        cfg = loadjson(cfg_file_name, 'SimplifyCell', 0);
        obj.meta = cfg.meta;
        obj.cfg = cfg.process;
        
        % Reformat string cell arrays in cfg
        if isfield(obj.cfg,'instruments2run')
          obj.cfg.instruments2run = cellfun(@(x) char(x), obj.cfg.instruments2run, 'UniformOutput', false);
        end
        tocheck1 = {'sync', 'split', 'bin', 'flag', 'calibrate', 'write'};
        tocheck2 = {'skip', 'skip', 'skip', 'skip', 'skip', 'skip'};
        for i=1:size(tocheck1,2)
          if isfield(obj.cfg,tocheck1{i}) && isfield(obj.cfg.(tocheck1{i}),tocheck2{i})
            obj.cfg.(tocheck1{i}).(tocheck2{i}) = cellfun(@(x) char(x), obj.cfg.(tocheck1{i}).(tocheck2{i}), 'UniformOutput', false);
          end
        end
        tocheck1 = {'qc', 'qc'};
        tocheck2 = {'global', 'specific'};
        tocheck3 = {'apply', 'run'};
        for i=1:size(tocheck1,2)
          if isfield(obj.cfg,tocheck1{i}) && isfield(obj.cfg.(tocheck1{i}),tocheck2{i}) && isfield(obj.cfg.(tocheck1{i}).(tocheck2{i}),tocheck3{i})
            obj.cfg.(tocheck1{i}).(tocheck2{i}).(tocheck3{i}) = cellfun(@(x) char(x), obj.cfg.(tocheck1{i}).(tocheck2{i}).(tocheck3{i}), 'UniformOutput', false);
          end
        end
        
        % Reformat parallel flag
        if ischar(obj.cfg.parallel)
          obj.cfg.parallel = str2double(obj.cfg.parallel);
          if isnan(obj.cfg.parallel)
            obj.cfg.parallel = 0;
            fprintf('WARNING: process.parallel forced to 0.\n');
          end
        end
        % Format cfg in flag
        cfg_flag_buf = obj.cfg.flag;
        obj.cfg.flag = struct();
        obj.cfg.flag.skip = cfg_flag_buf.skip;
        % Set all instruments with default parameters
        for i = fieldnames(cfg.instruments)'; i = i{1};
          obj.cfg.flag.(i) = struct();
          for p = fieldnames(cfg_flag_buf.default)'; p = p{1};
            switch p
              case {'tot', 'filt'}
                % Overwrite tot|filt specific parameters
                for sp = fieldnames(cfg_flag_buf.default.(p))'; sp = sp{1};
                  obj.cfg.flag.(i).(p).(sp) = cfg_flag_buf.default.(p).(sp);
                end
              otherwise
                % Overwrite with instrument specific parameters (both filt and tot)
                obj.cfg.flag.(i).tot.(p) = cfg_flag_buf.default.(p);
                obj.cfg.flag.(i).filt.(p) = cfg_flag_buf.default.(p);
            end
          end
        end
        % Set instruments with specific parameters
        % Note: orders of the parameters in the cfg files matters
        for i = fieldnames(cfg.instruments)'; i = i{1};
          if isfield(cfg_flag_buf, i)
            for p = fieldnames(cfg_flag_buf.(i))'; p = p{1};
              switch p
                case {'tot', 'filt'}
                  % Overwrite tot|filt specific parameters
                  for sp = fieldnames(cfg_flag_buf.(i).(p))'; sp = sp{1};
                    obj.cfg.flag.(i).(p).(sp) = cfg_flag_buf.(i).(p).(sp);
                  end
                otherwise
                  % Overwrite with instrument specific parameters (both filt and tot)
                  obj.cfg.flag.(i).tot.(p) = cfg_flag_buf.(i).(p);
                  obj.cfg.flag.(i).filt.(p) = cfg_flag_buf.(i).(p);
              end
            end
          end
        end
        
        % Initialize each instrument
        for i = fieldnames(cfg.instruments)'; i = i{1};
          switch cfg.instruments.(i).model
            case 'TSG'
              obj.instrument.(i) = TSG(cfg.instruments.(i));
            case 'FTH'
              obj.instrument.(i) = FTH(cfg.instruments.(i));
            case 'ACS'
              obj.instrument.(i) = ACS(cfg.instruments.(i));
            case 'BB3'
              obj.instrument.(i) = BB3(cfg.instruments.(i));
            case 'WSCD'
              obj.instrument.(i) = WSCD(cfg.instruments.(i));
            case 'LISST'
              obj.instrument.(i) = LISST(cfg.instruments.(i));
            case 'ECO'
              obj.instrument.(i) = ECO(cfg.instruments.(i));
            case 'FL'
              obj.instrument.(i) = FL(cfg.instruments.(i));
            case {'BB', 'BB3'}
              obj.instrument.(i) = BB(cfg.instruments.(i));
            case {'CD', 'WSCD'}
              obj.instrument.(i) = CD(cfg.instruments.(i));
            otherwise
              error('Instrument not supported: %s.', cfg.instruments.(i).model);
          end
        end
        
        % Set optional default parameters
        if ~isfield(obj.cfg,'instruments2run')
          obj.cfg.instruments2run = fieldnames(obj.instrument)';
        end
      end
    end
    
    % Pre-Process
    function Read(obj)
      for i=obj.cfg.instruments2run; i = i{1};
        obj.instrument.(i).ReadRaw(obj.cfg.days2run, obj.cfg.force_import, true);
      end
    end
    
    function Sync(obj)
      % Note: Run all days loaded (independent of days2run)
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(strcmp(i,obj.cfg.sync.skip))
          fprintf('SYNC: Skip %s\n', i);
        else
          obj.instrument.(i).Sync(obj.cfg.sync.delay.(i));
        end
      end
    end
    
    function QCRef(obj)
      switch obj.cfg.qcref.mode
        case 'ui'
          % Fresh selection does not take into account previous QC
          % TOTAL Sections
          fh = fig(31);
          title('Select Total Section'); fprintf('Select total section\n');
          yyaxis('left');
          plot(obj.instrument.(obj.cfg.qcref.reference).data.dt,...
               obj.instrument.(obj.cfg.qcref.reference).data.(obj.instrument.(obj.cfg.qcref.reference).view.varname), 'k', 'LineWidth', obj.instrument.(obj.cfg.qcref.reference).view.varcol);
          yyaxis('right'); 
          plot(obj.instrument.(obj.cfg.qcref.view).data.dt,...
               obj.instrument.(obj.cfg.qcref.view).data.(obj.instrument.(obj.cfg.qcref.view).view.varname)(:,obj.instrument.(obj.cfg.qcref.view).view.varcol),'.');
          ylabel(obj.instrument.(obj.cfg.qcref.view).view.varname);
          user_selection_total = guiSelectOnTimeSeries(fh);
          obj.instrument.(obj.cfg.qcref.reference).ApplyUserInput(user_selection_total, 'total');
          
          % FILTERED Sections
          fh = fig(31);
          title('Select Filtered Section'); fprintf('Select filtered section\n');
          yyaxis('left');
          plot(obj.instrument.(obj.cfg.qcref.reference).data.dt,...
               obj.instrument.(obj.cfg.qcref.reference).data.(obj.instrument.(obj.cfg.qcref.reference).view.varname), 'k', 'LineWidth', 1);
          yyaxis('right'); 
          plot(obj.instrument.(obj.cfg.qcref.view).data.dt,...
               obj.instrument.(obj.cfg.qcref.view).data.(obj.instrument.(obj.cfg.qcref.view).view.varname)(:,obj.instrument.(obj.cfg.qcref.view).view.varcol),'.');
          ylabel(obj.instrument.(obj.cfg.qcref.view).view.varname);
          user_selection_filtered = guiSelectOnTimeSeries(fh);
          obj.instrument.(obj.cfg.qcref.reference).ApplyUserInput(user_selection_filtered, 'filtered');
          
          filename = [obj.instrument.(obj.cfg.qcref.reference).path.wk, 'QCRef_UserSelection.json'];
          if exist(filename, 'file')
            % Load file
            file_selection = loadjson(filename);
            % Remove old (days2run) selections
            sel = min(obj.cfg.days2run) <= file_selection.total(:,1) & file_selection.total(:,1) < max(obj.cfg.days2run) + 1;
            file_selection.total(sel,:) = [];
            sel = min(obj.cfg.days2run) <= file_selection.filtered(:,1) & file_selection.filtered(:,1) < max(obj.cfg.days2run) + 1;
            file_selection.filtered(sel,:) = [];
            % Add new user selection
            file_selection.total = [file_selection.total; user_selection_total];
            file_selection.filtered = [file_selection.filtered; user_selection_filtered];
          else
            file_selection = struct('total', user_selection_total, 'filtered', user_selection_filtered);
          end
          % Save user selection
          savejson('',file_selection,filename);  
        case 'load'
          % Load previous QC and apply it
          file_selection = loadjson([obj.instrument.(obj.cfg.qcref.reference).path.wk, 'QCRef_UserSelection.json']);
          % Remove selection from days before & after days2run
          sel = file_selection.total(:,2) < min(obj.cfg.days2run) | max(obj.cfg.days2run) + 1 < file_selection.total(:,1);
          file_selection.total(sel,:) = [];
          sel = file_selection.filtered(:,2) < min(obj.cfg.days2run) | max(obj.cfg.days2run) + 1 < file_selection.filtered(:,1);
          file_selection.filtered(sel,:) = [];
          % Apply selection
          obj.instrument.(obj.cfg.qcref.reference).ApplyUserInput(file_selection.total, 'total');
          obj.instrument.(obj.cfg.qcref.reference).ApplyUserInput(file_selection.filtered, 'filtered');
        case 'skip'
          fprintf('WARNING: Reference is not QC.\n');
        otherwise
          error('Unknown mode.');
      end
    end
    
    function Split(obj)
       % Note: Run all days loaded (independent of days2run)
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(strcmp(i,obj.cfg.split.skip))
          fprintf('SPLIT: Skip %s (copy data to next level)\n', i);
          obj.instrument.(i).raw.tsw = obj.instrument.(i).data;
        elseif strcmp(obj.instrument.(i).split.mode, 'None')
          fprintf('SPLIT: Not available for %s\n', i);
        else
          obj.instrument.(i).Split(obj.instrument.(obj.cfg.split.reference),...
                                   obj.cfg.split.buffer.(i));
        end
      end
    end
    
    function Bin(obj)
      % Note: Run all days loaded (independent of days2run)
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(strcmp(i,obj.cfg.bin.skip))
          fprintf('BIN: Skip %s (copy data to next level)\n', i);
          obj.instrument.(i).bin.tsw = obj.instrument.(i).raw.tsw;
          obj.instrument.(i).bin.fsw = obj.instrument.(i).raw.fsw;
%           obj.instrument.(i).bin.diw = obj.instrument.(i).raw.diw;
        else
          obj.instrument.(i).Bin(obj.cfg.bin.bin_size.(i),...
                                 obj.cfg.bin.prctile_detection,...
                                 obj.cfg.bin.prctile_average,...
                                 obj.cfg.parallel);
        end
      end
    end
    
    function Flag(obj)
      % Automatic selection of good and bad data
      % Note: Run all days loaded (independent of days2run)
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(strcmp(i,obj.cfg.flag.skip))
          fprintf('FLAG: Skip %s (copy data to next level)\n', i);
          obj.instrument.(i).qc.tsw = obj.instrument.(i).bin.tsw;
          obj.instrument.(i).qc.fsw = obj.instrument.(i).bin.fsw;
        else
          obj.instrument.(i).Flag(obj.cfg.flag.(i))
        end
      end
    end
    
    function QC(obj)
      % Manual quality check of the data resulting in good and bad data
      switch obj.cfg.qc.mode
        case 'ui'
          if obj.cfg.qc.global.active
            % Display interactive figure
            foo = obj.instrument.(obj.cfg.qc.global.view);
            fh=visFlag(foo.raw.tsw, foo.raw.fsw,...
                       foo.qc.tsw, foo.suspect.tsw,...
                       foo.qc.fsw, foo.suspect.fsw,...
                       foo.view.varname, foo.view.varcol);
            title('Global QC');
            user_selection = guiSelectOnTimeSeries(fh);
            % For each instrument 
            for i=obj.cfg.qc.global.apply; i = i{1};
              if ~any(strcmp(obj.cfg.instruments2run, i)); continue; end
              % Apply user selection
              obj.instrument.(i).DeleteUserSelection(user_selection);
              % Save user selection
              filename = [obj.instrument.(i).path.wk i '_QCGlobal_UserSelection.json'];
              obj.updatejson_userselection_bad(filename, user_selection);
            end
          end
          if obj.cfg.qc.specific.active
            % For each instrument
            for i=obj.cfg.qc.specific.run; i = i{1};
              if ~any(strcmp(obj.cfg.instruments2run, i)); continue; end
              % Display interactive figure
              foo = obj.instrument.(i);
              fh=visFlag(foo.raw.tsw, foo.raw.fsw,...
                         foo.qc.tsw, foo.suspect.tsw,...
                         foo.qc.fsw, foo.suspect.fsw,...
                         foo.view.varname, foo.view.varcol);
              title([i ' QC']);
              user_selection = guiSelectOnTimeSeries(fh);
              % Apply user selection
              obj.instrument.(i).DeleteUserSelection(user_selection);
              % Save user selection
              filename = [obj.instrument.(i).path.wk i '_QCSpecific_UserSelection.json'];
              obj.updatejson_userselection_bad(filename, user_selection);
            end
          end
          if ~obj.cfg.qc.global.active && ~obj.cfg.qc.specific.active
            fprintf('WARNING: Quality check is NOT performed.\n');
          end
        case 'load'
          % Load previous QC files and apply them
          if obj.cfg.qc.global.active
            for i=obj.cfg.qc.global.apply; i = i{1};
              if ~any(strcmp(obj.cfg.instruments2run, i)); continue; end
              file_selection = loadjson([obj.instrument.(i).path.wk i '_QCGlobal_UserSelection.json']);
              obj.instrument.(i).DeleteUserSelection(file_selection.bad);
            end
          end
          if obj.cfg.qc.specific.active
            for i=obj.cfg.qc.specific.run; i = i{1};
              if ~any(strcmp(obj.cfg.instruments2run, i)); continue; end
              file_selection = loadjson([obj.instrument.(i).path.wk i '_QCSpecific_UserSelection.json']);
              obj.instrument.(i).DeleteUserSelection(file_selection.bad);
            end
          end
          if ~obj.cfg.qc.global.active && ~obj.cfg.qc.specific.active
            fprintf('WARNING: Quality check is NOT performed.\n');
          end
        case 'skip'
          fprintf('WARNING: Quality Check is NOT performed.\n');
        otherwise
          error('Unknown mode.');
      end
    end
    
    % Process
    function Calibrate(obj)
      % Note: Run all days loaded (independent of days2run)
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(strcmp(i,obj.cfg.calibrate.skip))
          fprintf('CALIBRATE: Skip %s (copy data to next level)\n', i);
          obj.instrument.(i).prod.a = obj.instrument.(i).qc.tsw;
        else
          switch obj.instrument.(i).model
            case 'ACS'
              obj.instrument.(i).Calibrate(obj.cfg.calibrate.(i).compute_dissolved,...
                                           obj.cfg.calibrate.(i).interpolation_method,...
                                           obj.instrument.(obj.cfg.calibrate.(i).CDOM_source),...
                                           obj.instrument.(obj.cfg.calibrate.(i).FTH_source));
            case 'BB'
              obj.instrument.(i).Calibrate(obj.cfg.calibrate.(i).compute_dissolved,...
                                           obj.instrument.(obj.cfg.calibrate.(i).TSG_source))
            otherwise
              obj.instrument.(i).Calibrate()
          end
        end
      end
    end
    
    % Write
    function Write(obj)
      % for each instrument
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(strcmp(i,obj.cfg.write.skip))
          fprintf('WRITE: Skip %s\n', i);
        else
          switch obj.cfg.write.mode
            case 'One file'
              % Save all days2run in one file
              obj.instrument.(i).Write([i 'ALL'], obj.cfg.days2run);
            case 'One day one file'
              % Save each day from days2run in independent files
              for d=obj.cfg.days2run
                obj.instrument.(i).Write([i datestr(d,'yyyymmdd')], d);
              end
            otherwise
              error('Unknow writing mode.');
          end
        end
      end
    end
    
    % Load
    function LoadProducts(obj)
      % for each instrument
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(strcmp(i,obj.cfg.write.skip))
          fprintf('WRITE: Skip %s\n', i);
        else
          switch obj.cfg.write.mode
            case 'One file'
              % Save all days2run in one file
              obj.instrument.(i).LoadProducts([i 'ALL'], obj.cfg.days2run);
            case 'One day one file'
              % Save each day from days2run in independent files
              for d=obj.cfg.days2run
                obj.instrument.(i).LoadProducts([i datestr(d,'yyyymmdd')], d);
              end
            otherwise
              error('Unknow writing mode.');
          end
        end
      end
    end
  end
  
  methods (Access=private)
    function updatejson_userselection_bad(obj, filename, user_selection)
      if exist(filename, 'file')
        % Load file
        file_selection = loadjson(filename);
        if isfield(file_selection, 'bad') && ~isempty(file_selection.bad)
          % Remove old (days2run) selections
          sel = min(obj.cfg.days2run) <= file_selection.bad(:,1) & file_selection.bad(:,1) < max(obj.cfg.days2run) + 1;
          file_selection.bad(sel,:) = [];
          % Add new user selection
          file_selection.bad = [file_selection.bad; user_selection];
        else
          file_selection = struct('bad', user_selection);
        end
      else
        file_selection = struct('bad', user_selection);
      end
      % Save user selection
      savejson('',file_selection,filename); 
    end
  end
  
end


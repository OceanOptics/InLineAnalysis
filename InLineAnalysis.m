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
        [~, ~, ext] = fileparts(cfg_file_name);
        switch ext
          case '.json'
            cfg = obj.ReadCfgJSON(cfg_file_name);
          case '.m'
            cfg = obj.ReadCfgM(cfg_file_name);
          otherwise
            error('Unknown configuration file type.');
        end
        
        obj.meta = cfg.meta;
        obj.cfg = cfg.process;
        
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
            case 'AC9'
              obj.instrument.(i) = AC9(cfg.instruments.(i));
%             case 'BB3'
%               obj.instrument.(i) = BB3(cfg.instruments.(i));
%             case 'WSCD'
%               obj.instrument.(i) = WSCD(cfg.instruments.(i));
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
            case 'ALFA'
              obj.instrument.(i) = ALFA(cfg.instruments.(i));
           case 'PAR'
              obj.instrument.(i) = PAR(cfg.instruments.(i));
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
    function ReadRaw(obj)
      % Read was renamed to ReadRaw on Oct 19, 2018
      for i=obj.cfg.instruments2run; i = i{1};
        fprintf('READ RAW: %s\n', i);
            obj.instrument.(i).ReadRaw(obj.cfg.days2run, obj.cfg.force_import, true);
      end
    end
    
    function ReadRawDI(obj)
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(strcmp(i,obj.cfg.di.skip))
          fprintf('READ DI: Skip %s\n', i);
        else
          fprintf('READ DI: %s\n', i);
          obj.instrument.(i).ReadRawDI(obj.cfg.days2run, obj.cfg.force_import, true);
        end
      end
    end
    
    function Sync(obj)
      % Note: Run all days loaded (independent of days2run)
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(strcmp(i,obj.cfg.sync.skip))
          fprintf('SYNC: Skip %s\n', i);
        else
          fprintf('SYNC: %s\n', i);
          obj.instrument.(i).Sync(obj.cfg.sync.delay.(i));
        end
      end
    end
    
    function SplitDetect (obj, MinFiltPeriod)
        fprintf('Detecting %s filter events...\n', obj.cfg.qcref.view);
        obj.instrument.FTH.data = SplitDetect(obj.cfg.qcref.view,...
            obj.instrument.(obj.cfg.qcref.view).data, obj.instrument.FTH.data, MinFiltPeriod);
        fprintf('Done\n');
    end
    
    function StepQC (obj, fudge_factor, bb_threshold)
      if any(nargin < 2 | nargin < 3)
          fudge_factor.filtered.a = 3;
          fudge_factor.filtered.c = 3;
          fudge_factor.total.a = 3;
          fudge_factor.total.c = 3;
      end
      instru = fieldnames(obj.instrument);
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(contains(i,{'ACS','BB','PAR'}))
            if ~isempty(obj.instrument.(i).raw.fsw)
                fprintf('Deleting bad values from %s filtered data...', i);
                if  any(contains(i,'ACS'))
                    lambda.a = obj.instrument.(i).lambda_a;
                    lambda.c = obj.instrument.(i).lambda_c;
                elseif  any(contains(i,'BB'))
                    lambda.bb = obj.instrument.(i).lambda;
                end
                if any(contains(i,{'ACS','BB'}))
                [obj.instrument.(i).raw.fsw, Nbad]= StepQC(obj.instrument.(i).raw.fsw,...
                    lambda, fudge_factor.filtered, obj.instrument.(instru{contains(instru, 'BB3')}).dark,...
                    bb_threshold);
                end
                if  any(contains(i,'ACS'))
                    fprintf('Done\n%4.2f%% of absorption and %4.2f%% of attenuation spectrum deleted from %s filtered data\n',...
                        Nbad.a, Nbad.c, i);
                elseif  any(contains(i,'BB'))
                    for ii = 1:size(Nbad.bb,2)
                        fprintf('Done\n%4.2f%% of beta%i deleted from %s filtered data\n',...
                            Nbad.bb(ii), lambda.bb(ii), i);
                    end
                end
                if isempty(obj.instrument.(i).raw.tsw); fprintf('StepQC [Done]\n');end
            end
            if ~isempty(obj.instrument.(i).raw.tsw)
                fprintf('Deleting bad values from %s total data... ', i);
                if  any(contains(i,'ACS'))
                    lambda.a = obj.instrument.(i).lambda_a;
                    lambda.c = obj.instrument.(i).lambda_c;
                elseif  any(contains(i,'BB'))
                    lambda.bb = obj.instrument.(i).lambda;
                elseif any(contains(i,'PAR'))
                    foo = obj.instrument.(i).raw.tsw.par./obj.instrument.PAR.scale > 4500 ...
                        | obj.instrument.(i).raw.tsw.par./obj.instrument.PAR.scale < 0;
                    obj.instrument.(i).raw.tsw(foo,:) = [];
                    foo = obj.instrument.(i).data.par./obj.instrument.PAR.scale > 4500 ...
                        | obj.instrument.(i).data.par./obj.instrument.PAR.scale < 0;
                    obj.instrument.(i).data(foo,:) = [];
                end
                if any(contains(i,{'ACS','BB'}))
                [obj.instrument.(i).raw.tsw, Nbad]= StepQC(obj.instrument.(i).raw.tsw,...
                    lambda, fudge_factor.total, obj.instrument.(instru{contains(instru, 'BB3')}).dark,...
                    bb_threshold);
                end
                if  any(contains(i,'ACS'))
                    fprintf('Done\n%4.2f%% of absorption and %4.2f%% of attenuation spectrum deleted from %s total data\n',...
                        Nbad.a, Nbad.c, i);
                elseif any(contains(i,'BB'))
                    for ii = 1:size(Nbad.bb,2)
                        fprintf('Done\n%4.2f%% of beta%i deleted from %s total data\n',...
                            Nbad.bb(ii), lambda.bb(ii), i);
                    end
                elseif any(contains(i,'PAR'))
                    fprintf('Done\n%i raw %s values deleted\n', sum(foo), i);
                end
            fprintf('StepQC [Done]\n');
            end
        end
      end
    end
    
    function DiagnosticPlot (obj, instrument, level)
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(contains(i,instrument))
          fprintf('%s Diagnostic plots\n', i);
          DiagnosticPlot(obj.instrument.(i), i, level);
        end
      end
    end
      
    function Stretch(obj)
      % Note: Run all days loaded (independent of days2run)
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(strcmp(i,obj.cfg.stretch.skip))
          fprintf('STRETCH: Skip %s\n', i);
        else
          fprintf('STRETCH: %s\n', i);
          obj.instrument.(i).Stretch(obj.cfg.stretch.delta.(i));
        end
      end
    end
    
    function QCRef(obj)
      switch obj.cfg.qcref.mode
        case 'ui'
          % Fresh selection does not take into account previous QC
          % TOTAL and FILTERED Sections
          fh = fig(31);
          title('Select total (t; red) and filtered (f; green) sections (q to save)'); fprintf('Select total (t; red) and filtered (f; green) sections (q to save)\n');
          yyaxis('left');
          plot(obj.instrument.(obj.cfg.qcref.reference).data.dt,...
               obj.instrument.(obj.cfg.qcref.reference).data.(obj.instrument.(obj.cfg.qcref.reference).view.varname), 'k', 'LineWidth', obj.instrument.(obj.cfg.qcref.reference).view.varcol);
          ylim([-0.1 1.1]);
          yyaxis('right');
          plot(obj.instrument.(obj.cfg.qcref.view).data.dt,...
               obj.instrument.(obj.cfg.qcref.view).data.(obj.instrument.(obj.cfg.qcref.view).view.varname)(:,obj.instrument.(obj.cfg.qcref.view).view.varcol),'.');
          if contains(obj.cfg.qcref.view, 'AC')
            ylabel([obj.instrument.(obj.cfg.qcref.view).view.varname ' ' ...
                num2str(round(obj.instrument.(obj.cfg.qcref.view).('lambda_ref')(obj.instrument.(obj.cfg.qcref.view).view.varcol),0)) 'nm'])
          elseif contains(obj.cfg.qcref.view, 'BB')
            ylabel([obj.instrument.(obj.cfg.qcref.view).view.varname ' ' ...
                num2str(round(obj.instrument.(obj.cfg.qcref.view).('lambda')(obj.instrument.(obj.cfg.qcref.view).view.varcol),0)) 'nm'])
          else
              ylabel([obj.instrument.(obj.cfg.qcref.view).view.varname obj.instrument.(obj.cfg.qcref.view).view.varcol]);
          end
          [user_selection_total, user_selection_filtered] = guiSelectOnTimeSeries(fh);
          obj.instrument.(obj.cfg.qcref.reference).ApplyUserInput(user_selection_total, 'total');
          obj.instrument.(obj.cfg.qcref.reference).ApplyUserInput(user_selection_filtered, 'filtered');
          
          filename = [obj.instrument.(obj.cfg.qcref.reference).path.ui, 'QCRef_UserSelection.json'];
          if exist(filename, 'file')
            % Load file
            file_selection = loadjson(filename);
            % Convert datestr to datenum for newer format
            try
              if ~isempty(file_selection.total)
                file_selection.total = [datenum(file_selection.total(1)), datenum(file_selection.total(2))];
              end
            catch
              if ~isempty(file_selection.total)
                file_selection.total = [datenum(cellfun(@(x) char(x), file_selection.total{1}', 'UniformOutput', false)),...
                    datenum(cellfun(@(x) char(x), file_selection.total{2}', 'UniformOutput', false))];
              end
            end
            try
              if ~isempty(file_selection.filtered)
                file_selection.filtered = [datenum(file_selection.filtered(1)), datenum(file_selection.filtered(2))];
              end
            catch
              if ~isempty(file_selection.filtered)
                file_selection.filtered = [datenum(cellfun(@(x) char(x), file_selection.filtered{1}', 'UniformOutput', false)),...
                    datenum(cellfun(@(x) char(x), file_selection.filtered{2}', 'UniformOutput', false))];
              end
            end
            % Remove old (days2run) selections
            if ~isempty(file_selection.total)
              sel = min(obj.cfg.days2run) <= file_selection.total(:,1) & file_selection.total(:,1) < max(obj.cfg.days2run) + 1;
              file_selection.total(sel,:) = [];
            end
            if ~isempty(file_selection.filtered)
              sel = min(obj.cfg.days2run) <= file_selection.filtered(:,1) & file_selection.filtered(:,1) < max(obj.cfg.days2run) + 1;
              file_selection.filtered(sel,:) = [];
            end
            % Add new user selection
            file_selection.total = [file_selection.total; user_selection_total];
            file_selection.filtered = [file_selection.filtered; user_selection_filtered];
          else
            file_selection = struct('total', user_selection_total, 'filtered', user_selection_filtered);
          end
          % Convert datenum to datestr for newer format
          if ~isempty(file_selection.total)
            file_selection.total = {datestr(file_selection.total(:,1)), datestr(file_selection.total(:,2))};
          end
          if ~isempty(file_selection.filtered)
            file_selection.filtered = {datestr(file_selection.filtered(:,1)), datestr(file_selection.filtered(:,2))};
          end
          % Save user selection
          if ~isdir(obj.instrument.(obj.cfg.qcref.reference).path.ui)
            mkdir(obj.instrument.(obj.cfg.qcref.reference).path.ui);
          end
          savejson('',file_selection,filename);  
        case 'load'
          fprintf('QCRef LOAD: %s\n', obj.cfg.qcref.reference);
          % Load previous QC and apply it
          file_selection = loadjson([obj.instrument.(obj.cfg.qcref.reference).path.ui, 'QCRef_UserSelection.json']);
          % Convert datestr to datenum for newer format
          try
            if ~isempty(file_selection.total)
              file_selection.total = [datenum(file_selection.total(1)), datenum(file_selection.total(2))];
            end
          catch
            if ~isempty(file_selection.total)
              file_selection.total = [datenum(cellfun(@(x) char(x), file_selection.total{1}', 'UniformOutput', false)),...
                  datenum(cellfun(@(x) char(x), file_selection.total{2}', 'UniformOutput', false))];
            end
          end
          try
            if ~isempty(file_selection.filtered)
              file_selection.filtered = [datenum(file_selection.filtered(1)), datenum(file_selection.filtered(2))];
            end
          catch
            if ~isempty(file_selection.filtered)
              file_selection.filtered = [datenum(cellfun(@(x) char(x), file_selection.filtered{1}', 'UniformOutput', false)),...
                  datenum(cellfun(@(x) char(x), file_selection.filtered{2}', 'UniformOutput', false))];
            end
          end
          % Remove selection from days before & after days2run
          if ~isempty(file_selection.total)
            sel = file_selection.total(:,2) < min(obj.cfg.days2run) | max(obj.cfg.days2run) + 1 < file_selection.total(:,1);
            file_selection.total(sel,:) = [];
          end
          if ~isempty(file_selection.filtered)
            sel = file_selection.filtered(:,2) < min(obj.cfg.days2run) | max(obj.cfg.days2run) + 1 < file_selection.filtered(:,1);
            file_selection.filtered(sel,:) = [];
          end
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
          fprintf('SPLIT: %s\n', i);
          obj.instrument.(i).Split(obj.instrument.(obj.cfg.split.reference),...
                                   obj.cfg.split.buffer.(i));
        end
      end
    end
    
    function Bin(obj)
      % Note: Run all days loaded (independent of days2run)
      for i=obj.cfg.instruments2run; i = i{1};
        obj.instrument.(i).bin.tsw = table();
        obj.instrument.(i).bin.fsw = table();
        if  any(strcmp(i,obj.cfg.bin.skip))
          fprintf('BIN: Skip %s (copy data to next level)\n', i);
          obj.instrument.(i).bin.tsw = obj.instrument.(i).raw.tsw;
          obj.instrument.(i).bin.fsw = obj.instrument.(i).raw.fsw;
        else
          fprintf('BIN: %s\n', i);
          obj.instrument.(i).Bin(obj.cfg.bin.bin_size.(i),...
                                 obj.cfg.bin.prctile_detection,...
                                 obj.cfg.bin.prctile_average,...
                                 obj.cfg.parallel,...
                                 obj.cfg.bin.mode);
%                                  obj.cfg.bin.method,...
        end
      end
    end
    
    function BinDI(obj)
      %%% NOTE: For DIW QC is done before the Binning %%%
      % Note: Run all days loaded (independent of days2run)
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(strcmp(i,obj.cfg.di.skip))
          fprintf('BIN DI: Skip %s\n', i);
        else
          fprintf('BIN DI: %s\n', i);
          obj.instrument.(i).BinDI(obj.cfg.di.bin.bin_size,...
                                   obj.cfg.bin.prctile_detection,...
                                   obj.cfg.bin.prctile_average,...
                                   obj.cfg.parallel);
%                                    obj.cfg.bin.method,...
                                   
        end
      end
    end
    
    % Flag is DEPRECATED: Use skip to move data to next level
    function Flag(obj)
      % Automatic selection of good and bad data
      % Note: Run all days loaded (independent of days2run)
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(strcmp(i,obj.cfg.flag.skip))
          fprintf('FLAG: Skip %s (copy data to next level)\n', i);
          obj.instrument.(i).qc.tsw = obj.instrument.(i).bin.tsw;
          obj.instrument.(i).qc.fsw = obj.instrument.(i).bin.fsw;
        else
          fprintf('FLAG: %s\n', i);
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
            fooflow = obj.instrument.FTH.bin.tsw;

            fh=visFlag(foo.raw.tsw, foo.raw.fsw,...
                       foo.qc.tsw, foo.suspect.tsw,...
                       foo.qc.fsw, foo.suspect.fsw,...
                       foo.view.varname, foo.view.varco,...
                       foo.raw.bad,fooflow);
            title('Global QC');
            user_selection = guiSelectOnTimeSeries(fh);
            % For each instrument 
            for i=obj.cfg.qc.global.apply; i = i{1};
              if ~any(strcmp(obj.cfg.instruments2run, i)); continue; end
              % Apply user selection
              obj.instrument.(i).DeleteUserSelection(user_selection);
              % Save user selection
              filename = [obj.instrument.(i).path.ui i '_QCGlobal_UserSelection.json'];
              if ~isdir(obj.instrument.(i).path.ui); mkdir(obj.instrument.(i).path.ui); end
              obj.updatejson_userselection_bad(filename, user_selection);
            end
          end
          if obj.cfg.qc.specific.active
            % For each instrument
            for i=obj.cfg.qc.specific.run; i = i{1};
              if ~any(strcmp(obj.cfg.instruments2run, i)); continue; end
              % Display interactive figure
              foo = obj.instrument.(i);
              fooflow = obj.instrument.FTH.bin.tsw;

              if ~isempty(foo.raw.tsw)
                fh=visFlag(foo.raw.tsw, foo.raw.fsw,...
                           foo.qc.tsw, foo.suspect.tsw,...
                           foo.qc.fsw, foo.suspect.fsw,...
                           foo.view.varname, foo.view.varcol,...
                           foo.raw.bad,fooflow);
              else
                fh=visFlag([], [],...
                           foo.qc.tsw, foo.suspect.tsw,...
                           foo.qc.fsw, foo.suspect.fsw,...
                           foo.view.varname, foo.view.varcol,...
                           foo.raw.bad,fooflow);
              end
              title([i ' QC' newline 'Trash section pressing t (q to save)']);
              fprintf('Trash section pressing t (q to save)\n');
              user_selection = guiSelectOnTimeSeries(fh);
              % Apply user selection
              obj.instrument.(i).DeleteUserSelection(user_selection);
              % Save user selection
              filename = [obj.instrument.(i).path.ui i '_QCSpecific_UserSelection.json'];
              if ~isdir(obj.instrument.(i).path.ui); mkdir(obj.instrument.(i).path.ui); end
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
              fprintf('QC LOAD Global: %s\n', i);
              file_selection = loadjson([obj.instrument.(i).path.ui i '_QCGlobal_UserSelection.json']);
              % Convert datestr to datenum for newer format
              if ~isempty(file_selection.bad)
                  file_selection.bad = [datenum(cellfun(@(x) char(x), file_selection.bad{1}', 'UniformOutput', false)),...
                      datenum(cellfun(@(x) char(x), file_selection.bad{2}', 'UniformOutput', false))];
              end
              obj.instrument.(i).DeleteUserSelection(file_selection.bad);
            end
          end
          if obj.cfg.qc.specific.active
            for i=obj.cfg.qc.specific.run; i = i{1};
              if ~any(strcmp(obj.cfg.instruments2run, i)); continue; end
              fprintf('QC LOAD Specific: %s\n', i);
              file_selection = loadjson([obj.instrument.(i).path.ui i '_QCSpecific_UserSelection.json']);
              % Convert datestr to datenum for newer format
              try
                  if ~isempty(file_selection.bad)
                      file_selection.bad = [datenum(file_selection.bad(1)), datenum(file_selection.bad(2))]; end
              catch
                  if ~isempty(file_selection.bad)
                      file_selection.bad = [datenum(cellfun(@(x) char(x), file_selection.bad{1}', 'UniformOutput', false)),...
                          datenum(cellfun(@(x) char(x), file_selection.bad{2}', 'UniformOutput', false))];
                  end
              end
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
    
    function QCDI(obj)
      %%% NOTE: For DIW QC is done before the Binning %%%
      % Manual quality check of the data resulting in good and bad data
      % Check if qc level is empty for each instrument
      for i=obj.cfg.instruments2run; i = i{1};
        obj.instrument.(i).qc.diw = obj.instrument.(i).raw.diw;
      end
      
      switch obj.cfg.di.qc.mode
        case 'ui'
          % For each instrument
          for i=obj.cfg.instruments2run; i = i{1};
            if ~any(strcmp(obj.cfg.instruments2run, i)) || any(strcmp(obj.cfg.di.skip, i)); continue; end
            % Display interactive figure
            foo = obj.instrument.(i);
            if isempty(foo.raw.diw); error('Empty raw diw\n'); end
            ColorSet = lines(2);
            fh = fig(52); hold('on');
            if iscell(foo.view.varname)
              plot(foo.raw.diw.dt, foo.raw.diw.(foo.view.varname{1})(:,foo.view.varcol), '.', 'Color', ColorSet(1,:));
              plot(foo.raw.diw.dt, foo.raw.diw.(foo.view.varname{2})(:,foo.view.varcol), '.', 'Color', ColorSet(2,:));
            else
              plot(foo.raw.diw.dt, foo.raw.diw.(foo.view.varname)(:,foo.view.varcol), '.', 'Color', ColorSet(1,:));
            end
            ylabel(foo.view.varname); title([i ' QC DI']);
            datetick2_doy();
            set(datacursormode(fh),'UpdateFcn',@data_cursor_display_date);
            % Get user selection
            user_selection = guiSelectOnTimeSeries(fh);
            % Apply user selection
            obj.instrument.(i).DeleteUserSelection(user_selection);
            % Save user selection
            filename = [obj.instrument.(i).path.ui i '_QCDI_UserSelection.json'];
            if ~isdir(obj.instrument.(i).path.ui); mkdir(obj.instrument.(i).path.ui); end
            obj.updatejson_userselection_bad(filename, user_selection);
          end
        case 'load'
          % Load previous QC files and apply them
          for i=obj.cfg.instruments2run; i = i{1};
            if ~any(strcmp(obj.cfg.instruments2run, i)) || any(strcmp(obj.cfg.di.skip, i)); continue; end
            fprintf('QC DI LOAD: %s\n', i);
            file_selection = loadjson([obj.instrument.(i).path.ui i '_QCDI_UserSelection.json']);
            % Convert datestr to datenum for newer format
            if ~isempty(file_selection.bad)
                file_selection.bad = [datenum(cellfun(@(x) char(x), file_selection.bad{1}', 'UniformOutput', false)),...
                    datenum(cellfun(@(x) char(x), file_selection.bad{2}', 'UniformOutput', false))];
            end
            obj.instrument.(i).DeleteUserSelection(file_selection.bad);
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
          fprintf('CALIBRATE: %s\n', i);
          switch obj.instrument.(i).model
            case 'AC9'
              obj.instrument.(i).Calibrate(obj.cfg.calibrate.(i).compute_dissolved,...
                                           obj.cfg.calibrate.(i).interpolation_method,...
                                           obj.instrument.(obj.cfg.calibrate.(i).CDOM_source),...
                                           obj.instrument.(obj.cfg.calibrate.(i).FTH_source));
            case 'ACS'
              obj.instrument.(i).Calibrate(obj.cfg.calibrate.(i).compute_dissolved,...
                                           obj.cfg.calibrate.(i).interpolation_method,...
                                           obj.instrument.(obj.cfg.calibrate.(i).CDOM_source),...
                                           obj.instrument.(obj.cfg.calibrate.(i).FTH_source));
            case {'BB', 'BB3'}
              obj.instrument.(i).Calibrate(obj.cfg.calibrate.(i).compute_dissolved,...
                                           obj.instrument.(obj.cfg.calibrate.(i).TSG_source),...
                                           obj.cfg.calibrate.(i).di_method)
            otherwise
              obj.instrument.(i).Calibrate()
          end
        end
      end
    end
    
    % Write
    function Write(obj, level)
      if nargin < 2; level = 'prod'; end
      % for each instrument
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(strcmp(i,obj.cfg.write.skip))
          fprintf('WRITE %s: Skip %s\n', level, i);
        else
          fprintf('WRITE %s: %s\n', level, i);
          switch obj.cfg.write.mode
            case 'One file'
              % Save all days2run in one file
              obj.instrument.(i).Write([i '_ALL'], obj.cfg.days2run, level);
            case 'One day one file'
              % Save each day from days2run in independent files
              for d=obj.cfg.days2run
                obj.instrument.(i).Write([i '_' datestr(d,'yyyymmdd')], d, level);
              end
            otherwise
              error('Unknow writing mode.');
          end
        end
      end
    end
    
    % Load
    function Read(obj, level)
      % LoadProducts is renamed to Read on Oct 19, 2018
      if nargin < 2; level = 'prod'; end
      % for each instrument
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(strcmp(i,obj.cfg.write.skip))
          fprintf('LOAD: Skip %s\n', i);
        else
          fprintf('LOAD: %s\n', i);
          switch obj.cfg.write.mode
            case 'One file'
              % Read all days2run in one file
              obj.instrument.(i).Read([i '_ALL'], obj.cfg.days2run, level);
            case 'One day one file'
              % Read each day from days2run in independent files
              for d=obj.cfg.days2run
                obj.instrument.(i).Read([i '_' datestr(d,'yyyymmdd')], d, level);
              end
            otherwise
              error('Unknow loading mode.');
          end
        end
      end
    end
    
    % Check size of data in each instrument
    function CheckDataStatus(obj)
      fprintf('---------------------------------------------------------------------------------------------+\n');
      fprintf('Instrument |   Data   |   Raw    |   Bin    |   QC     | Suspect  |   Bad    |    Prod       |\n');
      fprintf('-----------+----------+----------+----------+----------+----------+----------+---------------+\n');
      for i=obj.cfg.instruments2run; i = i{1};
        sdata = size(obj.instrument.(i).data);
        sraw = size(obj.instrument.(i).raw.tsw);
        sbin = size(obj.instrument.(i).bin.tsw);
        sqc = size(obj.instrument.(i).qc.tsw);
        ssuspect = size(obj.instrument.(i).suspect.tsw);
        sbad = size(obj.instrument.(i).bad.tsw);
        sprod = []; nprod = {};
        if ~isempty(fieldnames(obj.instrument.(i).prod))
          for t=fieldnames(obj.instrument.(i).prod); t = t{1};
            sprod(end+1,:) = size(obj.instrument.(i).prod.(t));
            nprod{end+1} = t;
          end
        end
        fprintf('%10s | %5dx%2d | %5dx%2d | %5dx%2d | %5dx%2d | %5dx%2d | %5dx%2d | ',...
                i, sdata(1), sdata(2), sraw(1), sraw(2), sbin(1), sbin(2), sqc(1), sqc(2),...
                ssuspect(1), ssuspect(2), sbad(1), sbad(2));
        if isempty(sprod)
          fprintf('None          | ');
        else
          for k = 1:size(sprod,1)
            fprintf('%3s(%5dx%2d) | ', nprod{k}, sprod(k,1), sprod(k,2));
          end
        end
        fprintf('\n');
      end
      fprintf('---------------------------------------------------------------------------------------------+\n');
    end
    
    % Merge products of two instruments (same model)
    function MergeProducts(obj, primary_instrument, secondary_instrument)
      % Data (at prod level) from the secondary instrument is copied to the
      % primary_instrument. The product table of the primary instrument is
      % then sorted by date & time (dt).

      fprintf('MERGING: %s << %s\n', primary_instrument, secondary_instrument);
      % For each product type (particulate, dissoved...)
      for f = fieldnames(obj.instrument.(primary_instrument).prod)'; f = f{1};
        ns = size(obj.instrument.(secondary_instrument).prod.(f),1);
        obj.instrument.(primary_instrument).prod.(f)(end+1:end+ns,:) = ...
          obj.instrument.(secondary_instrument).prod.(f);
        obj.instrument.(primary_instrument).prod.(f) = ...
          sortrows(obj.instrument.(primary_instrument).prod.(f));
      end
    end
    
  end
  
  
  
  methods (Static)
    % Load configuration file
    function cfg = ReadCfgJSON(cfg_file_name)
      % Read JSON configuration file (original format)
      fprintf('ReadCfgJSON method is DEPRECATED\n');
      
      % Load JSON file
      try
        cfg = loadjson(cfg_file_name, 'SimplifyCell', 0);
      catch
        error('Invalid configuration file: "%s"\n', cfg_file_name);
      end
      
      % Reformat string cell arrays in cfg
      if isfield(cfg.process,'instruments2run')
        cfg.process.instruments2run = cellfun(@(x) char(x), cfg.process.instruments2run, 'UniformOutput', false);
      end
      tocheck1 = {'di', 'sync', 'split', 'stretch', 'bin', 'flag', 'calibrate', 'write'};
      tocheck2 = {'skip', 'skip', 'skip', 'skip', 'skip', 'skip', 'skip', 'skip'};
      for i=1:size(tocheck1,2)
        if isfield(cfg.process,tocheck1{i}) && isfield(cfg.process.(tocheck1{i}),tocheck2{i})
          cfg.process.(tocheck1{i}).(tocheck2{i}) = cellfun(@(x) char(x), cfg.process.(tocheck1{i}).(tocheck2{i}), 'UniformOutput', false);
        end
      end
      tocheck1 = {'qc', 'qc'};
      tocheck2 = {'global', 'specific'};
      tocheck3 = {'apply', 'run'};
      for i=1:size(tocheck1,2)
        if isfield(cfg.process,tocheck1{i}) && isfield(cfg.process.(tocheck1{i}),tocheck2{i}) && isfield(cfg.process.(tocheck1{i}).(tocheck2{i}),tocheck3{i})
          cfg.process.(tocheck1{i}).(tocheck2{i}).(tocheck3{i}) = cellfun(@(x) char(x), cfg.process.(tocheck1{i}).(tocheck2{i}).(tocheck3{i}), 'UniformOutput', false);
        end
      end

      % Reformat parallel flag
      if ischar(cfg.process.parallel)
        cfg.process.parallel = str2double(cfg.process.parallel);
        if isnan(cfg.process.parallel)
          cfg.process.parallel = 0;
          fprintf('WARNING: process.parallel forced to 0.\n');
        end
      end
    end
    
    function cfg = ReadCfgM(filename)
      % Read matlab configuration file (current format)
      run(filename)
    end
  end
  
  methods (Access=private)
    function updatejson_userselection_bad(obj, filename, user_selection)
      if exist(filename, 'file')
        % Load file
        file_selection = loadjson(filename);
        % Convert datestr to datenum for newer format
        try
            if ~isempty(file_selection.bad)
                file_selection.bad = [datenum(file_selection.bad(1)), datenum(file_selection.bad(2))]; end
        catch
            if ~isempty(file_selection.bad)
                file_selection.bad = [datenum(cellfun(@(x) char(x), file_selection.bad{1}', 'UniformOutput', false)),...
                    datenum(cellfun(@(x) char(x), file_selection.bad{2}', 'UniformOutput', false))];
            end
        end
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
      % Convert datenum to datestr for newer format
      if ~isempty(file_selection.bad); file_selection.bad = {datestr(file_selection.bad(:,1)), datestr(file_selection.bad(:,2))}; end
      % Save user selection
      savejson('',file_selection,filename); 
    end
  end
  
end


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
            case 'HBB'
              obj.instrument.(i) = HBB(cfg.instruments.(i));
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
    
    function SplitDetect (obj, MinFiltPeriod, szFilt)
      if nargin < 2
        MinFiltPeriod = 65;
        szFilt = 10;
      elseif nargin < 3
        szFilt = 10;
      end
      fprintf('Detecting %s filter events...\n', obj.cfg.qcref.view);
      obj.instrument.FLOW.data = SplitDetect(obj.cfg.qcref.view,...
        obj.instrument.(obj.cfg.qcref.view).data, obj.instrument.FLOW.data, MinFiltPeriod, szFilt);
      fprintf('Done\n');
    end
    
    function RawAutoQC (obj, fudge_factor, bb_threshold, level)
      if nargin < 2
        fudge_factor.filtered.a = 3;
        fudge_factor.filtered.c = 3;
        fudge_factor.total.a = 3;
        fudge_factor.total.c = 3;
      elseif nargin < 3
        bb_threshold = 4000;
      elseif nargin < 4
        level = 'raw';
      end
      instru = fieldnames(obj.instrument);
      for i=obj.cfg.instruments2run; i = i{1};
        if any(contains(i,{'AC','BB','PAR'}))
          if ~isempty(obj.instrument.(i).(level).fsw)
%             fprintf('Deleting bad values from %s filtered data ...\n', i);
            if any(contains(i,'AC'))
              lambda.a = obj.instrument.(i).lambda_a;
              lambda.c = obj.instrument.(i).lambda_c;
            elseif  any(contains(i,'BB'))
              lambda.bb = obj.instrument.(i).lambda;
            end
            if any(contains(i,{'AC','BB'}))
              [obj.instrument.(i).(level).fsw, Nbad]= RawAutoQC(i, obj.instrument.(i).(level).fsw,...
                lambda, fudge_factor.filtered, obj.instrument.(instru{contains(instru, 'BB3')}).dark,...
                bb_threshold);
            end
            if any(contains(i,'AC'))
              fprintf('%4.2f%% of absorption and %4.2f%% of attenuation spectrum deleted from %s filtered %s data\n',...
                Nbad.a, Nbad.c, i, level);
            elseif any(contains(i,'BB'))
              fprintf('\n')
              fprintf('\n')
              fprintf('Data percentage deleted from %s filtered %s data:\n', i, level)
              lineseg = '-------+';
              headseg = '  %i  |';
              dataseg = ' %2.2f%% |';
              lin = '-----------+';
              head = '|  lambda  |';
              dat = ' %% deleted |';
              for ii = 1:size(Nbad.bb,2)
                lin = [lin lineseg];
                head = [head headseg];
                dat = [dat dataseg];
              end
              lin = [lin '+\n'];
              head = [head '\n'];
              dat = [dat '\n'];
              fprintf(head, lambda.bb)
              fprintf(lin);
              fprintf(dat, Nbad.bb)
%               for ii = 1:size(Nbad.bb,2)
%                 fprintf('%4.2f%% of beta%i deleted from %s filtered data\n',...
%                   Nbad.bb(ii), lambda.bb(ii), i);
%               end
            end
          else
            fprintf('No filtered data loaded: Skip\n');
          end
          if ~isempty(obj.instrument.(i).(level).tsw)
% %             fprintf('Deleting bad values from %s total data ...\n', i);
            if  any(contains(i,'AC'))
              lambda.a = obj.instrument.(i).lambda_a;
              lambda.c = obj.instrument.(i).lambda_c;
            elseif any(contains(i,'BB'))
              lambda.bb = obj.instrument.(i).lambda;
            elseif any(contains(i,'PAR'))
              foo = obj.instrument.(i).(level).tsw.par./obj.instrument.PAR.scale > 4500 ...
                | obj.instrument.(i).(level).tsw.par./obj.instrument.PAR.scale < 0;
              obj.instrument.(i).(level).tsw(foo,:) = [];
              foo = obj.instrument.(i).data.par./obj.instrument.PAR.scale > 4500 ...
                | obj.instrument.(i).data.par./obj.instrument.PAR.scale < 0;
              obj.instrument.(i).data(foo,:) = [];
            end
            if any(contains(i,{'AC','BB'}))
              [obj.instrument.(i).(level).tsw, Nbad]= RawAutoQC(i, obj.instrument.(i).(level).tsw,...
                lambda, fudge_factor.total, obj.instrument.(instru{contains(instru, 'BB3')}).dark,...
                bb_threshold);
            end
            if any(contains(i,'AC'))
              fprintf('%4.2f%% of absorption and %4.2f%% of attenuation spectrum deleted from %s total %s data\n',...
                Nbad.a, Nbad.c, i, level);
            elseif any(contains(i,'BB'))
              fprintf('\n')
              fprintf('\n')
              fprintf('Data percentage deleted from %s total %s data:\n', i, level)
              lineseg = '-------+';
              headseg = '  %i  |';
              dataseg = ' %2.2f%% |';
              lin = '-----------+';
              head = '|  lambda  |';
              dat = ' %% deleted |';
              for ii = 1:size(Nbad.bb,2)
                lin = [lin lineseg];
                head = [head headseg];
                dat = [dat dataseg];
              end
              lin = [lin '+\n'];
              head = [head '\n'];
              dat = [dat '\n'];
              fprintf(head, lambda.bb)
              fprintf(lin);
              fprintf(dat, Nbad.bb)
%               for ii = 1:size(Nbad.bb,2)
%                 fprintf('%4.2f%% of beta%i deleted from %s total data\n',...
%                   Nbad.bb(ii), lambda.bb(ii), i);
%               end
            elseif any(contains(i,'PAR'))
              fprintf('%i %s %s values deleted\n', sum(foo), level, i);
            end
          else
            fprintf('No total data loaded: Skip\n');
          end
          if ~isempty(obj.instrument.(i).(level).diw)
%             fprintf('Deleting bad values from %s dissolved data...\n', i);
            if any(contains(i,'AC'))
              lambda.a = obj.instrument.(i).lambda_a;
              lambda.c = obj.instrument.(i).lambda_c;
            elseif  any(contains(i,'BB'))
              lambda.bb = obj.instrument.(i).lambda;
            end
            if any(contains(i,{'AC','BB'}))
              [obj.instrument.(i).(level).diw, Nbad]= RawAutoQC(i, obj.instrument.(i).(level).diw,...
                lambda, fudge_factor.dissolved, obj.instrument.(instru{contains(instru, 'BB3')}).dark,...
                bb_threshold, true);
            end
            if any(contains(i,'AC'))
              fprintf('%4.2f%% of absorption and %4.2f%% of attenuation spectrum deleted from %s dissolved %s data\n',...
                Nbad.a, Nbad.c, i, level);
            elseif any(contains(i,'BB'))
              fprintf('\n')
              fprintf('\n')
              fprintf('Data percentage deletred from %s dissolved %s data:\n', i, level)
              lineseg = '-------+';
              headseg = '  %i  |';
              dataseg = ' %2.2f%% |';
              lin = '-----------+';
              head = '|  lambda  |';
              dat = ' %% deleted |';
              for ii = 1:size(Nbad.bb,2)
                lin = [lin lineseg];
                head = [head headseg];
                dat = [dat dataseg];
              end
              lin = [lin '+\n'];
              head = [head '\n'];
              dat = [dat '\n'];
              fprintf(head, lambda.bb)
              fprintf(lin);
              fprintf(dat, Nbad.bb)
%               for ii = 1:size(Nbad.bb,2)
%                 fprintf('%4.2f%% of beta%i deleted from %s dissolved data\n',...
%                   Nbad.bb(ii), lambda.bb(ii), i);
%               end
            end
          end
          fprintf('RawAutoQC [Done]\n')
        else
          fprintf('No RawAutoQC for %s [Done]\n', i)
        end
      end
    end
    
    function DiagnosticPlot (obj, instru, level, save_figure, toClean)
      if nargin < 3
        error('Not enough input argument')
      elseif nargin == 3
        save_figure = false;
        toClean = {'',''};
      elseif nargin == 4
        toClean = {'',''};
      elseif nargin > 5
        error('Too many input argument')
      end
      if size(toClean,2) < 2 || ~iscell(toClean)
        error("Indicate the table and variable names in cell array, e.g. {'p', 'ap'} to clean ACS product visualising ap spectrum")
      end
      for i=obj.cfg.instruments2run; i = i{1};
        if any(contains(i,instru))
          fprintf('%s Diagnostic plots\n', i);
          [user_selection] = DiagnosticPlot(obj.instrument.(i), i, level, ...
            save_figure, obj.meta.cruise, toClean);
          % Apply user selection
          if ~isempty(user_selection)
            if ~isfolder(obj.instrument.(i).path.ui); mkdir(obj.instrument.(i).path.ui); end
            if contains(i, 'AC')
%               obj.instrument.(i).DeleteUserSelection(user_selection, toClean{2}(1));
              obj.instrument.(i).DeleteUserSelection(user_selection, level{:}, toClean);
              % Save user selection
              if strcmp(toClean{1}, 'diw')
                filename = [obj.instrument.(i).path.ui i '_QCDIpickSpecific_UserSelection.json'];
              else
                filename = [obj.instrument.(i).path.ui i '_QCpickSpecific_UserSelection.json'];
              end
              obj.updatejson_userselection_bad(filename, user_selection, level{:}, toClean);
            else
              obj.instrument.(i).DeleteUserSelection(user_selection);
              % Save user selection
              if strcmp(toClean{1}, 'diw')
                filename = [obj.instrument.(i).path.ui i '_QCDIpickSpecific_UserSelection.json'];
              else
                filename = [obj.instrument.(i).path.ui i '_QCpickSpecific_UserSelection.json'];
              end
              obj.updatejson_userselection_bad(filename, user_selection);
            end
          end
        end
      end
    end
    
    function visProd_timeseries (obj)
      for i=obj.cfg.instruments2run; i = i{1};
        if  any(~contains(i, {'FLOW', 'FTH', 'Flow'}))
          fprintf('%s products time series plots\n', i);
          ifieldn = fieldnames(obj.instrument.(i).prod);
          for j=1:size(ifieldn,1)
            if ~strcmp(ifieldn{j}, 'QCfailed') && ~isempty(obj.instrument.(i).prod.(ifieldn{j}))
              if contains(i, {'BB'})
                visProd_timeseries(obj.instrument.(i).prod.(ifieldn{j}), i, ...
                  obj.instrument.(i).lambda);
              else
                visProd_timeseries(obj.instrument.(i).prod.(ifieldn{j}), i);
              end
            end
          end
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
                file_selection.total = [datenum(cellfun(@(x) char(x), file_selection.total{1}', 'un', 0)),...
                    datenum(cellfun(@(x) char(x), file_selection.total{2}', 'un', 0))];
              end
            end
            try
              if ~isempty(file_selection.filtered)
                file_selection.filtered = [datenum(file_selection.filtered(1)), datenum(file_selection.filtered(2))];
              end
            catch
              if ~isempty(file_selection.filtered)
                file_selection.filtered = [datenum(cellfun(@(x) char(x), file_selection.filtered{1}', 'un', 0)),...
                    datenum(cellfun(@(x) char(x), file_selection.filtered{2}', 'un', 0))];
              end
            end
%             % Remove old (days2run) selections REMOVED TO KEEP ALL HISTORY OF USER SELECTION
%             if ~isempty(file_selection.total)
%               sel = min(obj.cfg.days2run) <= file_selection.total(:,1) & file_selection.total(:,1) < max(obj.cfg.days2run) + 1;
%               file_selection.total(sel,:) = [];
%             end
%             if ~isempty(file_selection.filtered)
%               sel = min(obj.cfg.days2run) <= file_selection.filtered(:,1) & file_selection.filtered(:,1) < max(obj.cfg.days2run) + 1;
%               file_selection.filtered(sel,:) = [];
%             end
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
          if ~isfolder(obj.instrument.(obj.cfg.qcref.reference).path.ui)
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
              file_selection.total = [datenum(cellfun(@(x) char(x), file_selection.total{1}', 'un', 0)),...
                  datenum(cellfun(@(x) char(x), file_selection.total{2}', 'un', 0))];
            end
          end
          try
            if ~isempty(file_selection.filtered)
              file_selection.filtered = [datenum(file_selection.filtered(1)), datenum(file_selection.filtered(2))];
            end
          catch
            if ~isempty(file_selection.filtered)
              file_selection.filtered = [datenum(cellfun(@(x) char(x), file_selection.filtered{1}', 'un', 0)),...
                  datenum(cellfun(@(x) char(x), file_selection.filtered{2}', 'un', 0))];
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
            if isempty(obj.instrument.FLOW.qc.tsw)
              fooflow = obj.instrument.FLOW.bin.tsw;
            else
              fooflow = obj.instrument.FLOW.qc.tsw;
            end

            fh=visFlag(foo.raw.tsw, foo.raw.fsw, foo.qc.tsw, foo.suspect.tsw,...
                       foo.qc.fsw, foo.suspect.fsw, foo.view.varname, foo.view.varcol,...
                       foo.raw.bad, fooflow);
            title('Global QC');
            user_selection = guiSelectOnTimeSeries(fh);
            % For each instrument 
            for i=obj.cfg.qc.global.apply; i = i{1};
              if ~any(strcmp(obj.cfg.instruments2run, i)); continue; end
              % Apply user selection
              obj.instrument.(i).DeleteUserSelection(user_selection);
              % Save user selection
              filename = [obj.instrument.(i).path.ui i '_QCGlobal_UserSelection.json'];
              if ~isfolder(obj.instrument.(i).path.ui); mkdir(obj.instrument.(i).path.ui); end
              obj.updatejson_userselection_bad(filename, user_selection);
            end
          end
          if obj.cfg.qc.specific.active
            % For each instrument
            for i=obj.cfg.qc.specific.run; i = i{1};
              if ~any(strcmp(obj.cfg.instruments2run, i)); continue; end
              % Display interactive figure
              foo = obj.instrument.(i);
              if ~isempty(obj.instrument.FLOW.qc.tsw)
                fooflow = obj.instrument.FLOW.qc.tsw;
              elseif ~isempty(obj.instrument.FLOW.bin.tsw)
                fooflow = obj.instrument.FLOW.bin.tsw;
              elseif ~isempty(obj.instrument.FLOW.raw.tsw)
                fooflow = obj.instrument.FLOW.raw.tsw;
              else
                warning('Flow data not loaded')
                fooflow = obj.instrument.FLOW.qc.tsw;
              end
              if ~isfolder(obj.instrument.(i).path.ui); mkdir(obj.instrument.(i).path.ui); end
              if contains(i, 'AC') && ~obj.cfg.qc.qc_once_for_AandC  
                channel = {'a', 'c'};
                for j = channel
                  fh=visFlag(foo.raw.tsw, foo.raw.fsw, foo.qc.tsw, foo.suspect.tsw,...
                             foo.qc.fsw, foo.suspect.fsw, j{:}, foo.view.varcol,...
                             foo.raw.bad, fooflow);
                  title([i ' QC of "' j{:} '" only' newline 'Trash section pressing t (q to save)']);
                  % Get user selection
                  user_selection = guiSelectOnTimeSeries(fh);
                  % Apply user selection
                  obj.instrument.(i).DeleteUserSelection(user_selection, 'qc', ['tsw' j]);
                  obj.instrument.(i).DeleteUserSelection(user_selection, 'qc', ['fsw' j]);
                  % Save user selection
                  filename = [obj.instrument.(i).path.ui i '_QCSpecific_UserSelection.json'];
                  obj.updatejson_userselection_bad(filename, user_selection, 'qc', ['tsw' j]);
                  obj.updatejson_userselection_bad(filename, user_selection, 'qc', ['fsw' j]);
                  clf(52)
                end
              else
                if ~isempty(foo.raw.tsw)
                  fh=visFlag(foo.raw.tsw, foo.raw.fsw, foo.qc.tsw, foo.suspect.tsw,...
                             foo.qc.fsw, foo.suspect.fsw, foo.view.varname, foo.view.varcol,...
                             foo.raw.bad, fooflow);
                else
                  fh=visFlag([], [],...
                             foo.qc.tsw, foo.suspect.tsw, foo.qc.fsw, foo.suspect.fsw,...
                             foo.view.varname, foo.view.varcol, foo.raw.bad, fooflow);
                end
                title([i ' QC all' newline 'Trash both a & c section pressing t (q to save)']);
                fprintf('Trash section pressing t (q to save)\n');
                user_selection = guiSelectOnTimeSeries(fh);
                % Apply user selection
                obj.instrument.(i).DeleteUserSelection(user_selection);
                % Save user selection
                filename = [obj.instrument.(i).path.ui i '_QCSpecific_UserSelection.json'];
                obj.updatejson_userselection_bad(filename, user_selection);
              end
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
              if isfile([obj.instrument.(i).path.ui i '_QCSpecific_UserSelection.json'])
                file_selection = loadjson([obj.instrument.(i).path.ui i '_QCSpecific_UserSelection.json']);
                sel_toload = fieldnames(file_selection);
                % Convert datestr to datenum for newer format
                for j = 1:size(sel_toload, 1)
                  if ~isempty(file_selection.(sel_toload{j}))
                    if iscell(file_selection.(sel_toload{j}){1}) && iscell(file_selection.(sel_toload{j}){1})
                      file_selection.(sel_toload{j}) = [datenum(cellfun(@(x) char(x), ...
                        file_selection.(sel_toload{j}){1}, 'UniformOutput', false)),...
                        datenum(cellfun(@(x) char(x), file_selection.(sel_toload{j}){2}, 'UniformOutput', false))];
                    else
                      file_selection.(sel_toload{j}) = [datenum(file_selection.(sel_toload{j})(:, 1)),...
                        datenum(file_selection.(sel_toload{j})(:, 2))];
                    end
                  end
                  if contains(i, 'AC')
                    channel = strsplit(sel_toload{j}, 'bad_');
                    foo = strsplit(channel{end}, '_');
                    if size(foo, 2) == 1
                      level = 'qc';
                      channel = channel(end);
                    elseif size(foo, 2) == 2
                      level = 'qc';
                      channel = foo;
                    elseif size(foo, 2) == 3
                      level = foo{1};
                      channel = foo(2:3);
                    end
                    obj.instrument.(i).DeleteUserSelection(file_selection.(sel_toload{j}), ...
                      level, channel);
                  else
                    obj.instrument.(i).DeleteUserSelection(file_selection.(sel_toload{j}));
                  end
                end
              else
                fprintf(['Warning: ' obj.instrument.(i).path.ui i '_QCSpecific_UserSelection.json not found\n'])
              end
              if isfile([obj.instrument.(i).path.ui i '_QCpickSpecific_UserSelection.json'])
                file_pick_selection = loadjson([obj.instrument.(i).path.ui i '_QCpickSpecific_UserSelection.json']);
                sel_picktoload = fieldnames(file_pick_selection);
                for j = 1:size(sel_picktoload, 1)
                  % load hand picked bad values
                  if ~isempty(file_pick_selection.(sel_picktoload{j}))
                    if iscell(file_pick_selection.(sel_picktoload{j}){1})
                      file_pick_selection.(sel_picktoload{j}) = datenum(cellfun(@(x) char(x), ...
                        file_pick_selection.(sel_picktoload{j}){1}, 'UniformOutput', false));
                    else
                      file_pick_selection.(sel_picktoload{j}) = datenum(file_pick_selection.(sel_picktoload{j})(:, 1));
                    end
                  end
                  if contains(i, 'AC')
                    channel = strsplit(sel_picktoload{j}, 'bad_');
                    foo = strsplit(channel{end}, '_');
                    if size(foo, 2) == 1
                      level = 'qc';
                      channel = channel(end);
                    elseif size(foo, 2) == 2
                      level = 'qc';
                      channel = foo;
                    elseif size(foo, 2) == 3
                      level = foo{1};
                      channel = foo(2:3);
                    end
                    obj.instrument.(i).DeleteUserSelection(file_pick_selection.(sel_picktoload{j}), ...
                      level, channel);
                  else
                    obj.instrument.(i).DeleteUserSelection(file_pick_selection.(sel_picktoload{j}));
                  end
                end
              else
                fprintf(['Warning: ' obj.instrument.(i).path.ui i '_QCpickSpecific_UserSelection.json not found\n'])
              end
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
        if isempty(obj.instrument.(i).qc.diw)
          obj.instrument.(i).qc.diw = obj.instrument.(i).raw.diw;
        end
      end
      switch obj.cfg.di.qc.mode
        case 'ui'
          % For each instrument
          for i=obj.cfg.instruments2run; i = i{1};
            if ~any(strcmp(obj.cfg.instruments2run, i)) || any(strcmp(obj.cfg.di.skip, i)); continue; end
            % Display interactive figure
            foo = obj.instrument.(i);
            if isempty(foo.qc.diw)
              error('Empty qc diw \n');
            end
            % create folder for user input
            if ~isfolder(obj.instrument.(i).path.ui); mkdir(obj.instrument.(i).path.ui); end
            ColorSet = lines(2);
            fh = fig(52); hold('on');
            if contains(i, 'AC') && ~obj.cfg.di.qc.qc_once_for_AandC  
              channel = {'a', 'c'};
              for j = 1:size(channel,2)
                plot(foo.qc.diw.dt, foo.qc.diw.(channel{j})(:,foo.view.varcol), '.', 'Color', ColorSet(j,:));
                title([i ' QC of "' channel{j} '" only' newline 'Trash section pressing t (q to save)']);
                ylabel(channel{j});
                datetick2_doy();
                set(datacursormode(fh), 'UpdateFcn', @data_cursor_display_date);
                % Get user selection
                user_selection = guiSelectOnTimeSeries(fh);
                % Apply user selection
                obj.instrument.(i).DeleteUserSelection(user_selection, 'qc', ['diw' channel(j)]);
                % Save user selection
                filename = [obj.instrument.(i).path.ui i '_QCDI_UserSelection.json'];
                obj.updatejson_userselection_bad(filename, user_selection, 'qc', ['diw' channel(j)]);
                clf(52)
              end
            else
              plot(foo.qc.diw.dt, foo.qc.diw.(foo.view.varname)(:,foo.view.varcol), '.', 'Color', ColorSet(1,:));
              ylabel(foo.view.varname);
              title([i ' QC all' newline 'Trash both a & c section pressing t (q to save)']);
              datetick2_doy();
              set(datacursormode(fh), 'UpdateFcn', @data_cursor_display_date);
              % Get user selection
              user_selection = guiSelectOnTimeSeries(fh);
              % Apply user selection
              obj.instrument.(i).DeleteUserSelection(user_selection);
              % Save user selection
              filename = [obj.instrument.(i).path.ui i '_QCDI_UserSelection.json'];
              obj.updatejson_userselection_bad(filename, user_selection);
            end
          end
        case 'load'
          % Load previous QC DI files and apply them
          for i=obj.cfg.instruments2run; i = i{1};
            if ~any(strcmp(obj.cfg.instruments2run, i)) || any(strcmp(obj.cfg.di.skip, i)); continue; end
            fprintf('QC DI LOAD: %s\n', i);
            if isfile([obj.instrument.(i).path.ui i '_QCDI_UserSelection.json'])
              file_selection = loadjson([obj.instrument.(i).path.ui i '_QCDI_UserSelection.json']);
              sel_toload = fieldnames(file_selection);
              % Convert datestr to datenum for newer format
              for j = 1:size(sel_toload, 1)
                if ~isempty(file_selection.(sel_toload{j}))
                  if iscell(file_selection.(sel_toload{j}){1}) && iscell(file_selection.(sel_toload{j}){1})
                    file_selection.(sel_toload{j}) = [datenum(cellfun(@(x) char(x), ...
                      file_selection.(sel_toload{j}){1}, 'UniformOutput', false)),...
                      datenum(cellfun(@(x) char(x), file_selection.(sel_toload{j}){2}, 'UniformOutput', false))];
                  else
                    file_selection.(sel_toload{j}) = [datenum(file_selection.(sel_toload{j})(:, 1)),...
                      datenum(file_selection.(sel_toload{j})(:, 2))];
                  end
                end
                if contains(i, 'AC')
                  channel = strsplit(sel_toload{j}, 'bad_');
                  foo = strsplit(channel{end}, '_');
                  if size(foo, 2) == 1
                    level = 'qc';
                    channel = channel(end);
                  elseif size(foo, 2) == 2
                    level = 'qc';
                    channel = foo;
                  elseif size(foo, 2) == 3
                    level = foo{1};
                    channel = foo(2:3);
                  end
                  obj.instrument.(i).DeleteUserSelection(file_selection.(sel_toload{j}), ...
                    level, channel);
                else
                  obj.instrument.(i).DeleteUserSelection(file_selection.(sel_toload{j}));
                end
              end
            else
              fprintf(['Warning: ' obj.instrument.(i).path.ui i '_QCDI_UserSelection.json not found\n'])
            end
            if isfile([obj.instrument.(i).path.ui i '_QCDI_pickSpecific_UserSelection.json'])
              file_pick_selection = loadjson([obj.instrument.(i).path.ui i '_QCDI_pickSpecific_UserSelection.json']);
              sel_picktoload = fieldnames(file_pick_selection);
              for j = 1:size(sel_picktoload, 1)
                % load hand picked bad values
                if ~isempty(file_pick_selection.(sel_picktoload{j}))
                  if iscell(file_pick_selection.(sel_picktoload{j}){1}) && iscell(file_pick_selection.(sel_picktoload{j}){1})
                    file_pick_selection.(sel_picktoload{j}) = datenum(cellfun(@(x) char(x), file_pick_selection.(sel_picktoload{j}){1}, 'UniformOutput', false));
                  else
                    file_pick_selection.(sel_picktoload{j}) = datenum(file_pick_selection.(sel_picktoload{j})(:, 1));
                  end
                end
                if contains(i, 'AC')
                  channel = strsplit(file_pick_selection{j}, 'bad_');
                  foo = strsplit(channel{end}, '_');
                  if size(foo, 2) == 1
                    level = 'qc';
                    channel = channel(end);
                  elseif size(foo, 2) == 2
                    level = 'qc';
                    channel = foo;
                  elseif size(foo, 2) == 3
                    level = foo{1};
                    channel = foo(2:3);
                  end
                  obj.instrument.(i).DeleteUserSelection(file_pick_selection.(file_pick_selection{j}), ...
                    level, channel);
                else
                  obj.instrument.(i).DeleteUserSelection(file_pick_selection.(file_pick_selection{j}));
                end
              end
            else
              fprintf(['Warning: ' obj.instrument.(i).path.ui i '_QCDI_pickSpecific_UserSelection.json not found\n'])
            end
          end
        case 'skip'
          fprintf('WARNING: Quality Check is NOT performed.\n');
        otherwise
          error('Unknown mode.');
      end
    end
    
    function QCSwitchPosition(obj, shift_flow)
      if isempty(obj.instrument.(obj.cfg.qcref.view).qc.tsw)
        error('No %s qc tsw data loaded', obj.cfg.qcref.view)
      end
      if isempty(obj.instrument.(obj.cfg.qcref.view).qc.fsw)
        error('No %s qc fsw data loaded', obj.cfg.qcref.view)
      end
      if isempty(obj.instrument.FLOW.qc.tsw)
        error('No FLOW qc data loaded')
      end
      filtdt = datetime(obj.instrument.(obj.cfg.qcref.view).qc.fsw.dt, 'ConvertFrom', 'datenum');
      filtdt = dateshift(filtdt, 'start', 'minute');
      flowdt = datetime(obj.instrument.FLOW.qc.tsw.dt, 'ConvertFrom', 'datenum');
      flowdt = dateshift(flowdt, 'start', 'minute');
      fprintf('%.2f%% of filtered data is included into filtered switch position before QC\n', ...
        sum(ismember(filtdt, flowdt)) / size(filtdt,1) * 100)
      if nargin == 1
        % Automatically detect how much to shift flow data time so that switch filtered
        % position include the most filtered data in number of minutes
        shift_flow_list = (-20:20)';
        matchs = NaN(size(shift_flow_list));
        for sh = 1:size(shift_flow_list, 1)
          popoflow = obj.instrument.FLOW.qc.tsw(obj.instrument.FLOW.qc.tsw.swt == 1, :);
          if shift_flow_list(sh) > 0
            linetoadd_dt = (datetime(min(popoflow.dt), 'ConvertFrom', 'datenum') - ...
              minutes(shift_flow_list(sh)):minutes(1):datetime(min(popoflow.dt), 'ConvertFrom', 'datenum') - ...
              minutes(1))';
            linetoadd = table(datenum(linetoadd_dt), zeros(size(linetoadd_dt)), zeros(size(linetoadd_dt)), ...
              NaN(size(linetoadd_dt)), NaN(size(linetoadd_dt)), NaN(size(linetoadd_dt)), NaN(size(linetoadd_dt)), ...
              'VariableNames', popoflow.Properties.VariableNames);
            popoflow = [repmat(linetoadd, shift_flow_list(sh), 1); popoflow(1:end-shift_flow_list(sh),:)];
          elseif shift_flow_list(sh) < 0
            linetoadd_dt = (datetime(max(popoflow.dt), 'ConvertFrom', 'datenum') + ...
              minutes(shift_flow_list(sh)):minutes(1):datetime(min(popoflow.dt), 'ConvertFrom', 'datenum') + ...
              minutes(1))';
            linetoadd = table(datenum(linetoadd_dt), zeros(size(linetoadd_dt)), zeros(size(linetoadd_dt)), ...
              NaN(size(linetoadd_dt)), NaN(size(linetoadd_dt)), NaN(size(linetoadd_dt)), NaN(size(linetoadd_dt)), ...
              'VariableNames', popoflow.Properties.VariableNames);
            popoflow = [popoflow(abs(shift_flow_list(sh))+1:end, :); repmat(linetoadd, abs(shift_flow_list(sh)), 1)];
          end
          popoflow.dt = datenum(datetime(popoflow.dt, 'ConvertFrom', 'datenum') + minutes(shift_flow_list(sh)));
          filtdt = datetime(obj.instrument.(obj.cfg.qcref.view).qc.fsw.dt, 'ConvertFrom', 'datenum');
          filtdt = dateshift(filtdt, 'start', 'minute');
          flowdt = datetime(popoflow.dt, 'ConvertFrom', 'datenum');
          flowdt = dateshift(flowdt, 'start', 'minute');
          matchs(sh) = sum(ismember(filtdt, flowdt));
        end
        shift_flow = shift_flow_list(matchs == max(matchs));
        shift_flow = min(shift_flow);
        if shift_flow ~= 0
          fprintf('Switch position shifted automatically to best match filtered data: %i minute(s)\n', shift_flow)
        end
      elseif nargin > 2
        error('Too many input arguments')
      end
      % Apply shift to flow data
      shift_flow = round(shift_flow);
      if shift_flow > 0
        linetoadd_dt = (datetime(min(obj.instrument.FLOW.qc.tsw.dt), 'ConvertFrom', 'datenum') - ...
          minutes(shift_flow):minutes(1):datetime(min(obj.instrument.FLOW.qc.tsw.dt), 'ConvertFrom', 'datenum') - ...
          minutes(1))';
        linetoadd = table(datenum(linetoadd_dt), zeros(size(linetoadd_dt)), zeros(size(linetoadd_dt)), ...
          NaN(size(linetoadd_dt)), NaN(size(linetoadd_dt)), NaN(size(linetoadd_dt)), NaN(size(linetoadd_dt)), ...
          'VariableNames', obj.instrument.FLOW.qc.tsw.Properties.VariableNames);
        obj.instrument.FLOW.qc.tsw = [repmat(linetoadd, shift_flow, 1); obj.instrument.FLOW.qc.tsw(1:end-shift_flow,:)];
      elseif shift_flow < 0
        linetoadd_dt = (datetime(max(obj.instrument.FLOW.qc.tsw.dt), 'ConvertFrom', 'datenum') + ...
          minutes(shift_flow):minutes(1):datetime(min(obj.instrument.FLOW.qc.tsw.dt), 'ConvertFrom', 'datenum') + ...
          minutes(1))';
        linetoadd = table(datenum(linetoadd_dt), zeros(size(linetoadd_dt)), zeros(size(linetoadd_dt)), ...
          NaN(size(linetoadd_dt)), NaN(size(linetoadd_dt)), NaN(size(linetoadd_dt)), NaN(size(linetoadd_dt)), ...
          'VariableNames', obj.instrument.FLOW.qc.tsw.Properties.VariableNames);
        obj.instrument.FLOW.qc.tsw = [obj.instrument.FLOW.qc.tsw(abs(shift_flow)+1:end, :); repmat(linetoadd, abs(shift_flow), 1)];
      end
      obj.instrument.FLOW.qc.tsw.dt = datenum(datetime(obj.instrument.FLOW.qc.tsw.dt, 'ConvertFrom', 'datenum') + minutes(shift_flow));
      % create filter event duplicate when long period without filter event 
      fh = visFlag([], [], obj.instrument.(obj.cfg.qcref.view).qc.tsw, [], ...
        obj.instrument.(obj.cfg.qcref.view).qc.fsw, [], obj.instrument.(obj.cfg.qcref.view).view.varname, ...
        obj.instrument.(obj.cfg.qcref.view).view.varcol, [], obj.instrument.FLOW.qc.tsw);
      plot(obj.instrument.FLOW.qc.tsw.dt, obj.instrument.FLOW.qc.tsw.swt, '-k')
      title(['Select filter event to duplicate (press f)' newline 'Select new time slot for filter event duplicated (press s)' newline 'Change switch position (press t)'], ...
        'FontSize', 14)
      legend('Total', 'Filtered', 'Flow rate','switch position (1=filtered | 0=total)',...
        'AutoUpdate','off', 'FontSize', 12)
      [toswitch_t, toduplicate_f, newtime_s] = guiSelectOnTimeSeries(fh);
      % check number of entries
      if size(toduplicate_f,1) ~= size(newtime_s,1)
        fprintf('Warning: Inconsistent number of entries (f and t), new filter event ignored\n')
        toduplicate_f = [];
        newtime_s = [];
      end
      % duplicate filter events f and s commands
      for j = 1:size(toduplicate_f, 1)
        idx_acs = obj.instrument.(obj.cfg.qcref.view).qc.fsw.dt >= toduplicate_f(j, 1) & ...
          obj.instrument.(obj.cfg.qcref.view).qc.fsw.dt <= toduplicate_f(j, 2);
        if sum(idx_acs) > 0
          % copy closest filter event to specific new time slot
          new_filt = obj.instrument.(obj.cfg.qcref.view).qc.fsw(idx_acs, :);
          nts = dateshift(datetime(newtime_s(j), 'ConvertFrom', 'datenum'), 'start', 'minute');
          new_filt.dt = new_filt.dt - (new_filt.dt(1) - datenum(nts));
          obj.instrument.(obj.cfg.qcref.view).qc.fsw = [obj.instrument.(obj.cfg.qcref.view).qc.fsw; new_filt];
          % sort by date
          [~,b] = sort(obj.instrument.(obj.cfg.qcref.view).qc.fsw.dt);
          obj.instrument.(obj.cfg.qcref.view).qc.fsw = obj.instrument.(obj.cfg.qcref.view).qc.fsw(b,:);
          idx_flow_newfilt = obj.instrument.FLOW.qc.tsw.dt >= new_filt.dt(1) & obj.instrument.FLOW.qc.tsw.dt <= new_filt.dt(end);
          if sum(idx_flow_newfilt) > 0
            obj.instrument.FLOW.qc.tsw.swt(idx_flow_newfilt) = 1;
            obj.instrument.FLOW.qc.tsw.swt_avg_sd(idx_flow_newfilt) = 1;
          end
          remaining_idx = size(new_filt,1) - sum(idx_flow_newfilt);
          if remaining_idx > 0
            % create flow table
            new_flow = table();
            new_flow.dt = new_filt.dt(end-remaining_idx+1:end);
            new_flow.swt = ones(remaining_idx,1);
            new_flow.swt_avg_sd = zeros(remaining_idx,1);
            new_flow.swt_avg_n = NaN(remaining_idx,1);
            new_flow.spd = NaN(remaining_idx,1);
            new_flow.spd_avg_sd = NaN(remaining_idx,1);
            new_flow.spd_avg_n = NaN(remaining_idx,1);
            obj.instrument.FLOW.qc.tsw = [obj.instrument.FLOW.qc.tsw; new_flow];
          end

          flow_swt_off = table(datenum(datetime(new_filt.dt(end), 'ConvertFrom', 'datenum') + ...
            minutes(1)'), 0, 0, NaN, NaN, NaN, NaN, 'VariableNames', ...
            obj.instrument.FLOW.qc.tsw.Properties.VariableNames);
          obj.instrument.FLOW.qc.tsw = [obj.instrument.FLOW.qc.tsw; flow_swt_off];
           % sort by date
          [~,b] = sort(obj.instrument.FLOW.qc.tsw.dt);
          obj.instrument.FLOW.qc.tsw = obj.instrument.FLOW.qc.tsw(b,:);
        end
      end
      % change switch position t command
      for j = 1:size(toswitch_t, 1)
        idx_flow = obj.instrument.FLOW.qc.tsw.dt > toswitch_t(j, 1) & ...
          obj.instrument.FLOW.qc.tsw.dt < toswitch_t(j, 2);
        nbtoswitch = floor(minutes(abs(datetime(toswitch_t(j, 2), 'ConvertFrom', 'datenum') - ...
          datetime(toswitch_t(j, 1), 'ConvertFrom', 'datenum'))));
        if sum(idx_flow) < nbtoswitch % replace all FLOW data with new table
            obj.instrument.FLOW.qc.tsw(idx_flow, :) = [];
            closest_position = obj.instrument.FLOW.qc.tsw.swt(abs(obj.instrument.FLOW.qc.tsw.dt - mean(toswitch_t(j, :))) == ...
              min(abs(obj.instrument.FLOW.qc.tsw.dt - mean(toswitch_t(j, :))))); % get the closest switch position in time
            new_flow = table();
            new_flow.dt = datenum(datetime(toswitch_t(j, 1), 'ConvertFrom', 'datenum'):...
              minutes(1):datetime(toswitch_t(j, 2), 'ConvertFrom', 'datenum'))';
            if closest_position > 0 % set the switch position to the opposite to the closest switch position
              new_flow.swt = zeros(size(new_flow,1),1);
            else
              new_flow.swt = ones(size(new_flow,1),1);
            end
            new_flow.swt_avg_sd = zeros(size(new_flow,1),1);
            new_flow.swt_avg_n = NaN(size(new_flow,1),1);
            new_flow.spd = NaN(size(new_flow,1),1);
            new_flow.spd_avg_sd = NaN(size(new_flow,1),1);
            new_flow.spd_avg_n = NaN(size(new_flow,1),1);
            obj.instrument.FLOW.qc.tsw = [obj.instrument.FLOW.qc.tsw; new_flow];
             % sort by date
            [~,b] = sort(obj.instrument.FLOW.qc.tsw.dt);
            obj.instrument.FLOW.qc.tsw = obj.instrument.FLOW.qc.tsw(b,:);
        else
    %       obj.instrument.FLOW.qc.tsw.swt(idx_flow) = 1;
          if sum(obj.instrument.FLOW.qc.tsw.swt(idx_flow) == 0) > sum(obj.instrument.FLOW.qc.tsw.swt(idx_flow) > 0)
            obj.instrument.FLOW.qc.tsw.swt(idx_flow) = 1;
          else
            obj.instrument.FLOW.qc.tsw.swt(idx_flow) = 0;
          end
        end
      end
      % check for duplicats in flow data and delete
      [~, L, ~] = unique(obj.instrument.(obj.cfg.qcref.view).qc.fsw.dt,'first');
      indexToDump = not(ismember(1:numel(obj.instrument.(obj.cfg.qcref.view).qc.fsw.dt),L));
      if sum(indexToDump) > 0
        fprintf('Warning: %i identical dates in AC data => deleted\n', sum(indexToDump))
        obj.instrument.(obj.cfg.qcref.view).qc.fsw(indexToDump, :) = [];
      end
      [~, L, ~] = unique(obj.instrument.FLOW.qc.tsw.dt,'first');
      indexToDump = not(ismember(1:numel(obj.instrument.FLOW.qc.tsw.dt),L));
      if sum(indexToDump) > 0
        fprintf('Warning: %i identical dates in FLOW data => deleted\n', sum(indexToDump))
        obj.instrument.FLOW.qc.tsw(indexToDump, :) = [];
      end
      filtdt = datetime(obj.instrument.(obj.cfg.qcref.view).qc.fsw.dt, 'ConvertFrom', 'datenum');
      filtdt = dateshift(filtdt, 'start', 'minute');
      flowdt = datetime(obj.instrument.FLOW.qc.tsw.dt, 'ConvertFrom', 'datenum');
      flowdt = dateshift(flowdt, 'start', 'minute');
      fprintf('%.2f%% of filtered data is included into filtered switch position after QC\n', ...
        sum(ismember(filtdt, flowdt)) / size(filtdt,1) * 100)
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
                                           obj.instrument.(obj.cfg.calibrate.(i).FLOW_source),...
                                           obj.cfg.calibrate.(i).di_method, ...
                                           obj.cfg.calibrate.(i).compute_ad_aphi);
            case 'ACS'
              obj.instrument.(i).Calibrate(obj.cfg.calibrate.(i).compute_dissolved,...
                                           obj.cfg.calibrate.(i).interpolation_method,...
                                           obj.instrument.(obj.cfg.calibrate.(i).CDOM_source),...
                                           obj.instrument.(obj.cfg.calibrate.(i).FLOW_source),...
                                           obj.cfg.calibrate.(i).di_method, ...
                                           obj.cfg.calibrate.(i).compute_ad_aphi);
            case {'BB', 'BB3', 'HBB'}
              obj.instrument.(i).Calibrate(obj.cfg.calibrate.(i).compute_dissolved,...
                                           obj.instrument.(obj.cfg.calibrate.(i).TSG_source),...
                                           obj.instrument.(obj.cfg.calibrate.(i).FLOW_source),...
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
        fprintf('%s %s flushed\n', i, level);
        switch level
          case 'data'
            obj.instrument.(i).(level) = table();
          case 'bin'
            obj.instrument.(i).(level).tsw = table();
            obj.instrument.(i).(level).fsw = table();
          case 'qc'
            obj.instrument.(i).(level).tsw = table();
            obj.instrument.(i).(level).fsw = table();
            obj.instrument.(i).(level).diw = table();
          case 'prod'
            fna = fieldnames(obj.instrument.(i).(level));
            for j = 1:size(fna, 1)
              if ~isempty(fna{j})
                obj.instrument.(i).(level).(fna{j}) = table();
              end
            end
          otherwise
            error('Level unknown')
        end
        if  any(strcmp(i,obj.cfg.write.skip))
          fprintf('LOAD: Skip %s\n', i);
        else
          day2read = [min(obj.cfg.days2run)-1 obj.cfg.days2run max(obj.cfg.days2run)+1];
          fprintf('LOAD: %s\n', i);
          switch obj.cfg.write.mode
            case 'One file'
              % Read all days2run in one file
              obj.instrument.(i).Read([i '_ALL'], day2read, level);
            case 'One day one file'
              % Read each day from days2run in independent files
              for d=day2read
                obj.instrument.(i).Read([i '*_' datestr(d,'yyyymmdd')], d, level);
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
        sraw = size(obj.instrument.(i).raw.tsw) + size(obj.instrument.(i).raw.diw);
        sbin = size(obj.instrument.(i).bin.tsw) + size(obj.instrument.(i).bin.diw);
        sqc = size(obj.instrument.(i).qc.tsw) + size(obj.instrument.(i).qc.diw);
        ssuspect = size(obj.instrument.(i).suspect.tsw) + size(obj.instrument.(i).suspect.diw);
        sbad = size(obj.instrument.(i).bad.tsw) + size(obj.instrument.(i).bad.diw);
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
    function updatejson_userselection_bad(~, filename, user_selection, level, channel) % obj, 
      if nargin < 4
        level = 'qc';
        channel = '';
      elseif nargin < 5
        channel = ['_' level];
      else
        channel = join(['_' level '_' strjoin(channel, '_')],'');
      end
      if exist(filename, 'file')
        % Load file
        file_selection = loadjson(filename);
        fiedna = fieldnames(file_selection);
        for i = 1:size(fiedna,1)
          if ~isempty(file_selection.(fiedna{i}))
            % Convert datestr to datenum for newer format
            if iscell(file_selection.(fiedna{i}){1})
              if size(file_selection.(fiedna{i}), 2) == 2
                file_selection.(fiedna{i}) = [datenum(cellfun(@(x) char(x), ...
                  file_selection.(fiedna{i}){1}', 'UniformOutput', false)),...
                  datenum(cellfun(@(x) char(x), file_selection.(fiedna{i}){2}', ...
                  'UniformOutput', false))];
              elseif size(file_selection.(fiedna{i}), 2) == 1
                file_selection.(fiedna{i}) = datenum(cellfun(@(x) char(x), ...
                  file_selection.(fiedna{i}){1}', 'UniformOutput', false));
              else
                error('Date/time size in .json file not supported, check %\n', filename)
              end
            else
              if size(file_selection.(fiedna{i}), 2) == 2
                file_selection.(fiedna{i}) = [datenum(file_selection.(fiedna{i})(1)), ...
                  datenum(file_selection.(fiedna{i})(2))];
              elseif size(file_selection.(fiedna{i}), 2) == 1
                file_selection.(fiedna{i}) = datenum(file_selection.(fiedna{i})(1)');
              else
                error('Date/time size in .json file not supported, check %\n', filename)
              end
            end
          end
        end
        if isfield(file_selection, ['bad' channel]) && ~isempty(file_selection.(['bad' channel]))
%           % Remove old (days2run) selections REMOVED TO KEEP ALL HISTORY OF USER SELECTION
%           sel = min(obj.cfg.days2run) <= file_selection.(['bad' channel])(:,1) & ...
%             file_selection.(['bad' channel])(:,1) < max(obj.cfg.days2run) + 1;
%           file_selection.(['bad' channel])(sel,:) = [];
          % Add new user selection
          file_selection.(['bad' channel]) = [file_selection.(['bad' channel]); user_selection];
        else
          file_selection.(['bad' channel]) = user_selection;
          fiedna = fieldnames(file_selection);
        end
      else
        file_selection = struct(['bad' channel], user_selection);
        fiedna = fieldnames(file_selection);
      end
      % Convert datenum to datestr for newer format
      for i = 1:size(fiedna,1)
        if ~isempty(file_selection.(fiedna{i}))
          if size(file_selection.(fiedna{i}), 2) == 2
            file_selection.(fiedna{i}) = {datestr(file_selection.(fiedna{i})(:,1), 'dd-mmm-yyyy HH:MM:SS.FFF'), ...
              datestr(file_selection.(fiedna{i})(:,2), 'dd-mmm-yyyy HH:MM:SS.FFF')};
          elseif size(file_selection.(fiedna{i}), 2) == 1
            file_selection.(fiedna{i}) = {datestr(file_selection.(fiedna{i}), 'dd-mmm-yyyy HH:MM:SS.FFF')};
          end
        end
      end
      % Save user selection
      savejson('',file_selection,filename); 
    end
  end
end


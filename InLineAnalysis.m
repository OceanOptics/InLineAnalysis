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
      addpath('lib', 'instruments', 'packages', 'packages/datetick2_doy', ...
        'packages/spectral_color', 'packages/TEOS-10_subset');
      
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
            case {'NMEA','SC701','GP32','GPS32Tara','GP32Tara','GPSSC701Tara','19X','GPS19X','GPSCOMPASSAT'}
              obj.instrument.(i) = NMEA(cfg.instruments.(i));
            case {'TSG', 'SBE45', 'SBE3845'}
              obj.instrument.(i) = TSG(cfg.instruments.(i));
            case 'atlasTSG'
              obj.instrument.(i) = atlasTSG(cfg.instruments.(i));
            case {'FTH', 'ADU100'}
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
            case {'LISSTTau', 'LISSTTAU', 'TAU'}
              obj.instrument.(i) = LISSTTau(cfg.instruments.(i));
            case 'LISST200X'
              obj.instrument.(i) = LISST200X(cfg.instruments.(i));
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
      fprintf('%s config file loaded in ILA structure\n', cfg_file_name)
    end
    
    % Pre-Process
    function ReadRaw(obj)
      % Read was renamed to ReadRaw on Oct 19, 2018
      for i=obj.cfg.instruments2run; i = i{1};
        if ~isfield(obj.instrument, i)
          error('Instrument2run "%s" does not match any instrument name in the cfg file:  %s', ...
            i, strjoin(fieldnames(obj.instrument), ' / '))
        else
          fprintf('READ RAW: %s\n', i);
          obj.instrument.(i).ReadRaw(obj.cfg.days2run, obj.cfg.force_import, true);
        end
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
    
    function AutoQC (obj, level)
      tolerance = obj.cfg.qc.AutoQC_tolerance;
      saturation_threshold = obj.cfg.qc.AutoQC_Saturation_Threshold;
      if isempty(tolerance)
        tolerance.dissolved.a = 3;
        tolerance.dissolved.c = 3;
        tolerance.filtered.a = 3;
        tolerance.filtered.c = 3;
        tolerance.total.a = 3;
        tolerance.total.c = 3;
      elseif isempty(saturation_threshold)
        saturation_threshold.a = 50;
        saturation_threshold.c = 50;
        saturation_threshold.bb = 4000;
      elseif nargin < 2
        level = 'raw';
      end
      instru = fieldnames(obj.instrument);
      for i=obj.cfg.instruments2run; i = i{1};
        if any(contains(i,{'AC','BB','PAR','ALFA'}))
          if ~isempty(obj.instrument.(i).(level).fsw)
%             fprintf('Deleting bad values from %s filtered data ...\n', i);
            if any(contains(i,'AC'))
              lambda.a = obj.instrument.(i).lambda_a;
              lambda.c = obj.instrument.(i).lambda_c;
              bb_dark = [];
            elseif  any(contains(i,'BB3'))
              lambda.bb = obj.instrument.(i).lambda;
              bb_dark = obj.instrument.(instru{contains(instru, 'BB3')}).dark;
            elseif  any(contains(i,{'HBB', 'HyperBB'}))
              lambda.bb = obj.instrument.(i).lambda;
              bb_dark = [];
            end
            if any(contains(i,{'AC','BB'}))
              [obj.instrument.(i).(level).fsw, bad_spec, Nbad] = AutoQC(i, obj.instrument.(i).(level).fsw,...
                lambda, tolerance.filtered, bb_dark, saturation_threshold);
              obj.instrument.(i).(level).bad = [obj.instrument.(i).(level).bad; bad_spec];
            end
            if any(contains(i,'AC'))
              fprintf('%4.2f%% of absorption and %4.2f%% of attenuation spectra deleted from %s filtered %s data\n',...
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
            warning('No filtered data loaded: Skip');
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
              % get percentage of data deleted
              foo = sum(foo) / size(obj.instrument.(i).(level).tsw, 1) * 100;
            elseif any(contains(i,'ALFA'))
              alf_ar = table2array(obj.instrument.(i).(level).tsw);
              toqc = repmat(~contains(obj.instrument.(i).(level).tsw.Properties.VariableNames, ...
                {'dt', 'WL'}), size(alf_ar,1),1);
              % clean with derivative
              clean_cycles = 50;
              foo = 0;
              for ii = progress(1:clean_cycles)
                alfa_deriv = [zeros(1, size(alf_ar, 2)); diff(alf_ar, [], 1)];
                prc_deriv = repmat(prctile(alfa_deriv, 95), size(alf_ar,1),1);
                toclean = toqc & alfa_deriv > 1.5*prc_deriv & alfa_deriv(:,1) < 0.6;
                alf_ar(toclean) = NaN;
                foo = foo + sum(toclean(:));
              end
              % clean with absolute values
              alfa_prc = repmat(prctile(alf_ar, 95), size(alf_ar,1),1);
              toclean = toqc & alfa_deriv > 1.5*alfa_prc;
              alf_ar(toclean) = NaN;
              foo = foo + sum(toclean(:));
              % rebuild table
              obj.instrument.(i).(level).tsw = array2table(alf_ar, 'VariableNames', ...
                obj.instrument.(i).(level).tsw.Properties.VariableNames);
              % get percentage of data deleted
              foo = foo / sum(toqc(:));
            end
            if any(contains(i,{'AC','BB'}))
              [obj.instrument.(i).(level).tsw, bad_spec, Nbad] = AutoQC(i, obj.instrument.(i).(level).tsw,...
                lambda, tolerance.total, bb_dark, saturation_threshold);
              obj.instrument.(i).(level).bad = [obj.instrument.(i).(level).bad; bad_spec];
            end
            if any(contains(i,'AC'))
              fprintf('%4.2f%% of absorption and %4.2f%% of attenuation spectra deleted from %s total %s data\n',...
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
            elseif any(contains(i,{'PAR', 'ALFA'}))
              fprintf('%.3f%% of %s %s values deleted\n', foo, level, i);
            end
          else
            warning('No total data loaded: Skip');
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
              [obj.instrument.(i).(level).diw, bad_spec, Nbad] = AutoQC(i, obj.instrument.(i).(level).diw,...
                lambda, tolerance.dissolved, obj.instrument.(instru{contains(instru, 'BB3')}).dark,...
                saturation_threshold, true);
              obj.instrument.(i).(level).bad = [obj.instrument.(i).(level).bad; bad_spec];
            end
            if any(contains(i,'AC'))
              fprintf('%4.2f%% of absorption and %4.2f%% of attenuation spectra deleted from %s dissolved %s data\n',...
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
            end
          end
          fprintf('AutoQC [Done]\n')
        else
          fprintf('No AutoQC for %s [Done]\n', i)
        end
      end
    end
    
    function SpectralQC (obj, instru, level, save_figure, toClean)
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
        error("Indicate the table and variable names in cell array, e.g. {'p', 'ap'} to clean ACS product visualising ap spectra")
      end
      for i=obj.cfg.instruments2run; i = i{1};
        if any(contains(i,instru))
          fprintf('%s Spectral QCs\n', i);
          [user_selection] = SpectralQC(obj.instrument.(i), obj.cfg.days2run, i, level, ...
            save_figure, obj.meta.cruise, toClean);
          % Apply user selection
          if ~isempty(user_selection)
            if ~isfolder(obj.instrument.(i).path.ui); mkdir(obj.instrument.(i).path.ui); end
            obj.instrument.(i).DeleteUserSelection(user_selection, level{:}, toClean);
            if strcmp(level{:}, 'prod')
              obj.instrument.(i).DeleteUserSelection(user_selection, 'qc', {'tsw', strrep(toClean{2}, 'p', '')});
            end
            % Save user selection
            if strcmp(toClean{1}, 'diw')
              filename = fullfile(obj.instrument.(i).path.ui, [i '_QCDI_Specific_UserSelection.mat']);
            else
              filename = fullfile(obj.instrument.(i).path.ui, [i '_QCpickSpecific_UserSelection.mat']);
            end
            obj.update_userselection_bad(filename, user_selection, false, level{:}, toClean);
          end
        end
      end
    end
    
    function visProd_timeseries(obj)
      for i=obj.cfg.instruments2run; i = i{1};
        if any(~contains(i, {'FLOW', 'FTH', 'Flow'}))
          fprintf('%s products time series plots\n', i);
          ifieldn = fieldnames(obj.instrument.(i).prod);
          for j=1:size(ifieldn,1)
            if ~strcmp(ifieldn{j}, 'FiltStat') && ~strcmp(ifieldn{j}, 'QCfailed') && ...
                ~isempty(obj.instrument.(i).prod.(ifieldn{j}))
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
      if isempty(obj.instrument.(obj.cfg.qcref.view).data.dt)
        error('%s data table is empty, make sure to view the instrument you are trying to process or load your data before running QCref', ...
          obj.cfg.qcref.view)
      end
      switch obj.cfg.qcref.mode
        case 'ui'
          %% create new FTH data for missing data
          % round date time to second
          obj.instrument.FLOW.data.dt = datenum(floor(datevec(obj.instrument.FLOW.data.dt)));
          % delete duplicates (bug in flowcontrol software)
          [~, L, ~] = unique(obj.instrument.FLOW.data.dt,'first');
          indexToDump = not(ismember(1:numel(obj.instrument.FLOW.data.dt),L));
          obj.instrument.FLOW.data(indexToDump, :) = [];
          obj.instrument.FLOW.data(~isfinite(obj.instrument.FLOW.data.dt),:) = [];
          % Fresh selection does not take into account previous QC
          % TOTAL and FILTERED Sections
          fh = fig(31);
          title(['\fontsize{22}Switch position QC:' newline '\fontsize{18}Select total (t; \color{red}red\color{black}) and filtered (f; \color{green}green\color{black}) sections' newline '\fontsize{14}Press q to save and quit (close graph to cancel and quit)'], 'interpreter', 'tex');
          fprintf('Select total (t; red) and filtered (f; green) sections (q to save)\n');
          yyaxis('left');
          plot(obj.instrument.(obj.cfg.qcref.reference).data.dt,...
               obj.instrument.(obj.cfg.qcref.reference).data.(obj.instrument.(obj.cfg.qcref.reference).view.varname), ...
               'k', 'LineWidth', obj.instrument.(obj.cfg.qcref.reference).view.varcol);
          ylim([-0.1 1.1]);
          ax = gca; ax.YColor = 'k';
          ylabel('Switch position')
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
          datetick2_doy();
          legend('switch position (1=filtered | 0=total)', obj.instrument.(obj.cfg.qcref.view).view.varname, 'FontSize', 14, 'AutoUpdate','off')
          [user_selection.total, user_selection.filtered] = guiSelectOnTimeSeries(fh);
          obj.instrument.(obj.cfg.qcref.reference).ApplyUserInput(user_selection.total, 'total');
          obj.instrument.(obj.cfg.qcref.reference).ApplyUserInput(user_selection.filtered, 'filtered');
          filename = fullfile(obj.instrument.(obj.cfg.qcref.reference).path.ui, 'QCRef_UserSelection.mat');
          if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
            file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
            save(filename, 'file_selection')
          end
          if isfile(filename)
            % Load file
            load(filename, 'file_selection');
            if obj.cfg.qcref.remove_old
              % Remove old (days2run) selections
              if ~isempty(file_selection.total)
                sel = min(obj.cfg.days2run) <= file_selection.total(:,1) & file_selection.total(:,1) < max(obj.cfg.days2run) + 1;
                file_selection.total(sel,:) = [];
              end
              if ~isempty(file_selection.filtered)
                sel = min(obj.cfg.days2run) <= file_selection.filtered(:,1) & file_selection.filtered(:,1) < max(obj.cfg.days2run) + 1;
                file_selection.filtered(sel,:) = [];
              end
            end
            % Add new user selection
            file_selection.total = [file_selection.total; user_selection.total];
            file_selection.filtered = [file_selection.filtered; user_selection.filtered];
          else
            file_selection = user_selection;
          end
          % Save user selection
          if ~isfolder(obj.instrument.(obj.cfg.qcref.reference).path.ui)
            mkdir(obj.instrument.(obj.cfg.qcref.reference).path.ui);
          end
          save(filename, 'file_selection');  
        case 'load'
          fprintf('QCRef LOAD: %s\n', obj.cfg.qcref.reference);
          % Load previous QC and apply it
          filename = fullfile(obj.instrument.(obj.cfg.qcref.reference).path.ui, 'QCRef_UserSelection.mat');
          if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
            file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
            save(filename, 'file_selection')
          end
          if isfile(filename)
            load(filename, 'file_selection');
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
            if ~isempty(file_selection.total)
              obj.instrument.(obj.cfg.qcref.reference).ApplyUserInput(file_selection.total, 'total');
            end
            if ~isempty(file_selection.filtered)
              obj.instrument.(obj.cfg.qcref.reference).ApplyUserInput(file_selection.filtered, 'filtered');
            end
            fh = fig(31);
            title(['\fontsize{22}Switch position QC:' newline '\fontsize{18}Previous selection applied' newline 'Make sure filtered and total are properly seleceted'], 'interpreter', 'tex');
            fprintf('Select total (t; red) and filtered (f; green) sections (q to save)\n');
            yyaxis('left');
            plot(obj.instrument.(obj.cfg.qcref.reference).data.dt,...
                 obj.instrument.(obj.cfg.qcref.reference).data.(obj.instrument.(obj.cfg.qcref.reference).view.varname), ...
                 'k', 'LineWidth', obj.instrument.(obj.cfg.qcref.reference).view.varcol);
            ylim([-0.1 1.1]);
            ax = gca; ax.YColor = 'k';
            ylabel('Switch position')
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
            datetick2_doy();
            legend('switch position (1=filtered | 0=total)', obj.instrument.(obj.cfg.qcref.view).view.varname, 'FontSize', 14, 'AutoUpdate','off')
          else
            fprintf(['Warning: ' filename ' not found\n'])
          end
        case 'skip'
          fprintf('WARNING: Reference is not QC.\n');
        otherwise
          error('Unknown mode.');
      end
    end
    
    function Split(obj)
       % Note: Run all days loaded (independent of days2run)
      for i=obj.cfg.instruments2run; i = i{1};
        if isempty(obj.instrument.(i).data)
          error('%s data table is empty', i)
        end
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
      % add all instrument from instrument2run into the skip so that flag
      % never bugs
      obj.cfg.flag.skip = [obj.cfg.flag.skip; obj.cfg.instruments2run(:)];
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
            foo = obj.instrument.(obj.cfg.qc.global.view{:});
            if isempty(obj.instrument.FLOW.qc.tsw)
              fooflow = obj.instrument.FLOW.bin.tsw;
            else
              fooflow = obj.instrument.FLOW.qc.tsw;
            end
            fh=visFlag(foo.raw.tsw, foo.raw.fsw, foo.qc.tsw, foo.suspect.tsw,...
                       foo.qc.fsw, foo.suspect.fsw, foo.view.varname, foo.view.varcol,...
                       foo.raw.bad, fooflow, obj.instrument.FLOW.view.spd_variable);
            title('\fontsize{22}\color{red}Global QC: \fontsize{18}\color{black}Press t to trash section (press q to save)', 'interpreter', 'tex');
            fprintf('Global QC: Press t to trash section (press q to save)\n');
            user_selection = guiSelectOnTimeSeries(fh);
            % For each instrument 
            for i=obj.cfg.qc.global.apply(:)'; i = i{1};
              if ~any(strcmp(obj.cfg.instruments2run, i)); continue; end
              % Apply user selection
              obj.instrument.(i).DeleteUserSelection(user_selection);
              % Save user selection
              filename = fullfile(fileparts(obj.instrument.(i).path.ui), 'QCGlobal_UserSelection.mat');
              if ~isfolder(fileparts(fileparts(obj.instrument.(i).path.ui)))
                mkdir(fileparts(fileparts(obj.instrument.(i).path.ui)));
              end
              obj.update_userselection_bad(filename, user_selection, obj.cfg.qc.remove_old);
            end
          else
            % Load previous globalQC if file exists and apply
            for i=obj.cfg.qc.global.apply(:)'; i = i{1};
              if ~any(strcmp(obj.cfg.instruments2run, i)); continue; end
              filename = fullfile(fileparts(obj.instrument.(i).path.ui), 'QCGlobal_UserSelection.mat');
              if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
                file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
                save(filename, 'file_selection')
              end
              if isfile(filename)
                fprintf('QC LOAD Global selection: %s ... ', i);
                load(filename, 'file_selection');
                fprintf('done\n')
                obj.instrument.(i).DeleteUserSelection(file_selection.bad);
              end
            end
          end
          if obj.cfg.qc.specific.active
            % For each instrument
            for i=obj.cfg.qc.specific.run(:)'; i = i{1};
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
              if contains(i, 'AC')
                channel = {'a', 'c'};
              elseif contains(i, {'TSG', 'SBE45', 'SBE3845'})
                if ~strcmp(foo.qc.tsw.Properties.VariableNames, 's')
                  channel = {foo.temperature_variable, 'c'};
                else
                  channel = {foo.temperature_variable, 's'};
                end
              end
              if contains(i, {'AC', 'TSG', 'SBE45', 'SBE3845'}) && ~obj.cfg.qc.qc_once_for_all
                for j = channel
                  fh=visFlag(foo.raw.tsw, foo.raw.fsw, foo.qc.tsw, foo.suspect.tsw,...
                             foo.qc.fsw, foo.suspect.fsw, j{:}, foo.view.varcol,...
                             foo.raw.bad, fooflow, obj.instrument.FLOW.view.spd_variable);
                  title(['\fontsize{22}\color{red}' i ' QC of "' j{:} '" only:' newline '\fontsize{18}\color{black}Press t to trash section (press q to save)'], 'interpreter', 'tex');
                  % title([i ' QC of "' j{:} '" only' newline 'Press t to trash section (press q to save)']);
                  fprintf([i ' QC of "' j{:} '" only: Press t to trash section (press q to save)\n']);
                  % Get user selection
                  user_selection = guiSelectOnTimeSeries(fh);
                  % Apply user selection
                  obj.instrument.(i).DeleteUserSelection(user_selection, 'qc', ['tsw' j]);
                  obj.instrument.(i).DeleteUserSelection(user_selection, 'qc', ['fsw' j]);
                  % Save user selection
                  filename = fullfile(obj.instrument.(i).path.ui, [i '_QCSpecific_UserSelection.mat']);
                  obj.update_userselection_bad(filename, user_selection, obj.cfg.qc.remove_old, ...
                          'qc', ['tsw' j]);
                  obj.update_userselection_bad(filename, user_selection, obj.cfg.qc.remove_old, ...
                          'qc', ['fsw' j]);
                  clf(52)
                end
              elseif contains(i, 'ALFA') && ~obj.cfg.qc.qc_once_for_all
                channel = foo.qc.tsw.Properties.VariableNames;
                channel(contains(channel, {'dt', '_avg_sd', '_avg_n'})) = [];
                for j = channel
                  % delete crazy values
                  if contains(j{:}, 'WL')
                    foo.qc.tsw.(j{:})(foo.qc.tsw.(j{:}) > 1000 | foo.qc.tsw.(j{:}) < 540) = NaN;
                  else
                    foo.qc.tsw.(j{:})(foo.qc.tsw.(j{:}) > 100) = NaN;
                  end
                  fh=visFlag(foo.raw.tsw, foo.raw.fsw, foo.qc.tsw, foo.suspect.tsw,...
                             foo.qc.fsw, foo.suspect.fsw, j{:}, foo.view.varcol,...
                             foo.raw.bad, fooflow, obj.instrument.FLOW.view.spd_variable);
                  title(['\fontsize{22}\color{red}' i ' QC of "' j{:} '" only:' newline '\fontsize{18}\color{black}Press t to trash section (press q to save)'], 'interpreter', 'tex');
                  % title([i ' QC of "' j{:} '" only' newline 'Press t to trash section (press q to save)']);
                  fprintf([i ' QC of "' j{:} '" only: Press t to trash section (press q to save)\n']);
                  % Get user selection
                  user_selection = guiSelectOnTimeSeries(fh);
                  % Apply user selection
                  obj.instrument.(i).DeleteUserSelection(user_selection, 'qc', ['tsw' j]);
                  % Save user selection
                  filename = fullfile(obj.instrument.(i).path.ui, [i '_QCSpecific_UserSelection.mat']);
                  obj.update_userselection_bad(filename, user_selection, obj.cfg.qc.remove_old, ...
                          'qc', ['tsw' j]);
                  clf(52)
                end
              else
                if ~isempty(foo.raw.tsw)
                  fh=visFlag(foo.raw.tsw, foo.raw.fsw, foo.qc.tsw, foo.suspect.tsw,...
                             foo.qc.fsw, foo.suspect.fsw, foo.view.varname, foo.view.varcol,...
                             foo.raw.bad, fooflow, obj.instrument.FLOW.view.spd_variable);
                else
                  fh=visFlag([], [],...
                             foo.qc.tsw, foo.suspect.tsw, foo.qc.fsw, foo.suspect.fsw,...
                             foo.view.varname, foo.view.varcol, foo.raw.bad, fooflow,...
                             obj.instrument.FLOW.view.spd_variable);
                end
                title(['\fontsize{22}\color{red}' i ' QC all variables:' newline '\fontsize{18}\color{black}Press t to trash section (press q to save)'], 'interpreter', 'tex');
                % title([i ' specific QC all' newline 'Trash full section pressing t (q to save)']);
                fprintf([i ' QC all: Press t to trash section (press q to save)\n']);
                user_selection = guiSelectOnTimeSeries(fh);
                % Apply user selection
                obj.instrument.(i).DeleteUserSelection(user_selection);
                % Save user selection
                filename = fullfile(obj.instrument.(i).path.ui, [i '_QCSpecific_UserSelection.mat']);
                obj.update_userselection_bad(filename, user_selection, obj.cfg.qc.remove_old);
              end
            end
          end
          if ~obj.cfg.qc.global.active && ~obj.cfg.qc.specific.active
            fprintf('WARNING: Quality check is NOT performed.\n');
          end
        case 'load'
          % Load previous QC files and apply them
          % if obj.cfg.qc.global.active
          for i=obj.cfg.qc.global.apply(:)'; i = i{1};
            if ~any(strcmp(obj.cfg.instruments2run, i)); continue; end
            fprintf('QC LOAD Global selection: %s ... ', i);
            filename = fullfile(fileparts(obj.instrument.(i).path.ui), 'QCGlobal_UserSelection.mat');
            if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
              file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
              save(filename, 'file_selection')
            end
            if isfile(filename)
              load(filename, 'file_selection');
              fprintf('done\n')
              obj.instrument.(i).DeleteUserSelection(file_selection.bad);
            else
              fprintf(['Warning: ' filename ' not found\n'])
            end
          end
          % end
          if obj.cfg.qc.specific.active
            for i=obj.cfg.qc.specific.run(:)'; i = i{1};
              if ~any(strcmp(obj.cfg.instruments2run, i)); continue; end
              fprintf('QC LOAD Specific slection: %s ... ', i);
              filename = fullfile(obj.instrument.(i).path.ui, [i '_QCSpecific_UserSelection.mat']);
              if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
                file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
                save(filename, 'file_selection')
              end
              if isfile(filename)
                fprintf('done\n')
                load(filename, 'file_selection');
                sel_toload = fieldnames(file_selection);
                for j = 1:size(sel_toload, 1)
                  % Keep only selection of the day2run
                  file_selection.(sel_toload{j})(any(file_selection.(sel_toload{j}) < min(obj.cfg.days2run) | ...
                    file_selection.(sel_toload{j}) > max(obj.cfg.days2run)+1, 2)) = [];
                  if ~isempty(file_selection.(sel_toload{j}))
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
                end
              else
                fprintf(['Warning: ' filename ' not found\n'])
              end
              fprintf('QC LOAD Specific pick selection: %s ... ', i);
              filename = fullfile(obj.instrument.(i).path.ui, [i '_QCpickSpecific_UserSelection.mat']);
              if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
                file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
                save(filename, 'file_selection')
              end
              if isfile(filename)
                % load hand picked bad values
                fprintf('done\n')
                load(filename, 'file_selection');
                sel_picktoload = fieldnames(file_selection);
                for j = progress(1:size(sel_picktoload, 1))
                  % Keep only selection of the day2run
                  file_selection.(sel_picktoload{j})(any(file_selection.(sel_picktoload{j}) < min(obj.cfg.days2run) | ...
                    file_selection.(sel_picktoload{j}) > max(obj.cfg.days2run)+1, 2)) = [];
                  if ~isempty(file_selection.(sel_picktoload{j}))
                    channel = strsplit(sel_picktoload{j}, 'bad_');
                    foo = strsplit(channel{end}, '_');
                    if size(foo, 2) == 1
                      level = 'qc';
                      channel = {'all'};
                    elseif size(foo, 2) == 2
                      level = 'qc';
                      channel = foo;
                    elseif size(foo, 2) == 3
                      level = foo{1};
                      channel = foo(2:3);
                    end
                    if contains(i, 'AC')
                      obj.instrument.(i).DeleteUserSelection(file_selection.(sel_picktoload{j}), ...
                        level, channel);
                      if strcmp(level, 'prod')
                        obj.instrument.(i).DeleteUserSelection(file_selection.(sel_picktoload{j}), ...
                          'qc', {'tsw', strrep(channel{2}, 'p', '')});
                      end
                    else
                      obj.instrument.(i).DeleteUserSelection(file_selection.(sel_picktoload{j}), level);
                      if strcmp(level, 'prod')
                        obj.instrument.(i).DeleteUserSelection(file_selection.(sel_picktoload{j}), 'qc');
                      end
                    end
                  end
                end
              else
                fprintf(['Warning: ' filename ' not found\n'])
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
            foo = obj.instrument.(i); %(obj.instrument.(i).dt, ;
            if isempty(foo.qc.diw)
              error('Empty qc diw \n');
            end
            % create folder for user input
            filename = fullfile(obj.instrument.(i).path.ui, [i '_QCDI_UserSelection.mat']);
            if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
              file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
              save(filename, 'file_selection')
            end
            if ~isfolder(obj.instrument.(i).path.ui); mkdir(obj.instrument.(i).path.ui); end
            ColorSet = lines(2);
            fh = fig(52); hold('on');
            if contains(i, 'AC') && ~obj.cfg.di.qc.qc_once_for_all  
              channel = {'a', 'c'};
              for j = 1:size(channel,2)
                plot(foo.qc.diw.dt, foo.qc.diw.(channel{j})(:,foo.view.varcol), 'o', 'Color', ColorSet(j,:));
                title([i ' QC of "' channel{j} '" only' newline 'Press t to trash section (press q to save)']);
                ylabel(channel{j});
                datetick2_doy();
                set(datacursormode(fh), 'UpdateFcn', @data_cursor_display_date);
                % Get user selection
                user_selection = guiSelectOnTimeSeries(fh);
                % Apply user selection
                obj.instrument.(i).DeleteUserSelection(user_selection, 'qc', ['diw' channel(j)]);
                % Save user selection
                obj.update_userselection_bad(filename, user_selection, obj.cfg.di.qc.remove_old, ...
                        'qc', ['diw' channel(j)]);
                clf(52)
              end
            else
              plot(foo.qc.diw.dt, foo.qc.diw.(foo.view.varname)(:,foo.view.varcol), 'o', 'Color', ColorSet(1,:));
              ylabel(foo.view.varname);
              title([i ' QC all' newline 'Trash full section pressing t (q to save)']);
              datetick2_doy();
              set(datacursormode(fh), 'UpdateFcn', @data_cursor_display_date);
              % Get user selection
              user_selection = guiSelectOnTimeSeries(fh);
              % Apply user selection
              obj.instrument.(i).DeleteUserSelection(user_selection);
              % Save user selection
              obj.update_userselection_bad(filename, user_selection, obj.cfg.di.qc.remove_old);
            end
          end
        case 'load'
          % Load previous QC DI files and apply them
          for i=obj.cfg.instruments2run; i = i{1};
            if ~any(strcmp(obj.cfg.instruments2run, i)) || any(strcmp(obj.cfg.di.skip, i)); continue; end
            fprintf('QC DI LOAD: %s\n', i);
            filename = fullfile(obj.instrument.(i).path.ui, [i '_QCDI_UserSelection.mat']);
            if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
              file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
              save(filename, 'file_selection')
            end
            if isfile(filename)
              % load bad DI values
              load(filename, 'file_selection');
              sel_toload = fieldnames(file_selection);
              for j = 1:size(sel_toload, 1)
                % Keep only selection of the day2run
                file_selection.(sel_toload{j})(~any(min(obj.cfg.days2run) < file_selection.(sel_toload{j}) & ...
                  file_selection.(sel_toload{j}) < max(obj.cfg.days2run), 2)) = [];
                if ~isempty(file_selection.(sel_toload{j}))
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
              end
            else
              fprintf(['Warning: ' filename ' not found\n'])
            end
            filename = fullfile(obj.instrument.(i).path.ui, [i '_QCDI_pickSpecific_UserSelection.mat']);
            if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
              file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
              save(filename, 'file_selection')
            end
            if isfile(filename)
              % load hand picked bad DI values
              load(filename, 'file_selection');
              sel_picktoload = fieldnames(file_selection);
              for j = 1:size(sel_picktoload, 1)
                % Keep only selection of the day2run
                file_selection.(sel_picktoload{j})(~any(min(obj.cfg.days2run) < file_selection.(sel_picktoload{j}) & ...
                  file_selection.(sel_picktoload{j}) < max(obj.cfg.days2run), 2)) = [];
                if ~isempty(file_selection.(sel_picktoload{j}))
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
                  if contains(i, 'AC')
                    obj.instrument.(i).DeleteUserSelection(file_selection.(sel_picktoload{j}), ...
                      level, channel);
                    if strcmp(level, 'prod')
                      obj.instrument.(i).DeleteUserSelection(file_selection.(sel_picktoload{j}), ...
                        'qc', {'fsw', strrep(channel{2}, 'g', '')});
                    end
                  else
                    obj.instrument.(i).DeleteUserSelection(file_selection.(sel_picktoload{j}));
                  end
                end
              end
            else
              fprintf(['Warning: ' filename ' not found\n'])
            end
          end
        case 'skip'
          fprintf('WARNING: Quality Check is NOT performed.\n');
        otherwise
          error('Unknown mode.');
      end
    end
    
    function QCSwitchPosition(obj)
      if isfield(obj.cfg.calibrate.(obj.cfg.qcref.view), 'filt_method')
        if strcmp(obj.cfg.calibrate.(obj.cfg.qcref.view).filt_method, 'exponential_fit')
          QCSwitchPosition(obj.instrument.(obj.cfg.qcref.view), obj.instrument.FLOW, 'raw')
        else
          QCSwitchPosition(obj.instrument.(obj.cfg.qcref.view), obj.instrument.FLOW, 'qc')
        end
      else
        QCSwitchPosition(obj.instrument.(obj.cfg.qcref.view), obj.instrument.FLOW, 'qc')
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
          if contains(i, {'AC','LISSTTau','LISSTTAU','LISST-Tau','TAU','CSTAR'})
            if isempty(obj.cfg.calibrate.(i).CDOM_source)
              cdom_source = [];
            else
              cdom_source = obj.instrument.(obj.cfg.calibrate.(i).CDOM_source);
            end
          end
          switch obj.instrument.(i).model
            case 'AC9'
              obj.instrument.(i).Calibrate(obj.cfg.days2run, ...
                                           obj.cfg.calibrate.(i).compute_dissolved,...
                                           obj.cfg.calibrate.(i).interpolation_method,...
                                           cdom_source,...
                                           obj.instrument.(obj.cfg.calibrate.(i).FLOW_source),...
                                           obj.cfg.calibrate.(i).di_method, ...
                                           obj.cfg.calibrate.(i).scattering_correction, ...
                                           obj.cfg.calibrate.(i).compute_ad_aphi, ...
                                           obj.instrument.(obj.cfg.calibrate.(i).TSG_source), ...
                                           obj.cfg.min_nb_pts_per_cluster, ...
                                           obj.cfg.time_weight_for_cluster);
            case 'ACS'
              obj.instrument.(i).Calibrate(obj.cfg.days2run, ...
                                           obj.cfg.calibrate.(i).compute_dissolved,...
                                           obj.cfg.calibrate.(i).interpolation_method,...
                                           cdom_source,...
                                           obj.instrument.(obj.cfg.calibrate.(i).FLOW_source),...
                                           obj.cfg.calibrate.(i).di_method, ...
                                           obj.cfg.calibrate.(i).scattering_correction, ...
                                           obj.cfg.calibrate.(i).compute_ad_aphi, ...
                                           obj.instrument.(obj.cfg.calibrate.(i).TSG_source), ...
                                           obj.cfg.min_nb_pts_per_cluster, ...
                                           obj.cfg.time_weight_for_cluster);
            case {'BB', 'BB3', 'HBB'}
              obj.instrument.(i).Calibrate(obj.cfg.days2run, ...
                                           obj.cfg.calibrate.(i).compute_dissolved,...
                                           obj.instrument.(obj.cfg.calibrate.(i).TSG_source),...
                                           obj.instrument.(obj.cfg.calibrate.(i).FLOW_source),...
                                           obj.cfg.calibrate.(i).di_method,...
                                           obj.cfg.calibrate.(i).filt_method)
            case {'FL'}
              obj.instrument.(i).Calibrate(obj.cfg.days2run, ...
                                           obj.cfg.calibrate.(i).compute_dissolved,...
                                           obj.instrument.(obj.cfg.calibrate.(i).FLOW_source),...
                                           obj.cfg.calibrate.(i).di_method,...
                                           obj.cfg.calibrate.(i).filt_method)
            case 'CD'
              obj.instrument.(i).Calibrate(obj.cfg.calibrate.(i).compute_dissolved)
            case {'LISST', 'LISST100X', 'LISST100x', 'LISST200X', 'LISST200x'}
              obj.instrument.(i).Calibrate(obj.cfg.days2run, ...
                                           obj.cfg.calibrate.(i).compute_dissolved,...
                                           obj.instrument.(obj.cfg.calibrate.(i).FLOW_source),...
                                           obj.cfg.calibrate.(i).di_method)
            case {'LISSTTau','LISSTTAU','LISST-Tau','TAU'}
              obj.instrument.(i).Calibrate(oobj.cfg.days2run, ...
                                           bj.cfg.calibrate.(i).compute_dissolved,...
                                           obj.cfg.calibrate.(i).interpolation_method,...
                                           cdom_source,...
                                           obj.instrument.(obj.cfg.calibrate.(i).FLOW_source),...
                                           obj.cfg.calibrate.(i).di_method);
            otherwise
              obj.instrument.(i).Calibrate()
          end
        end
      end
    end
    
    % Write
    function Write(obj, level, part_or_diw)
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
              obj.instrument.(i).Write([i '_ALL'], obj.cfg.days2run, level, part_or_diw);
            case 'One day one file'
              % Save each day from days2run in independent files
              for d=obj.cfg.days2run
                obj.instrument.(i).Write([i '_' datestr(d,'yyyymmdd')], d, level, part_or_diw);
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
          case 'raw'
            obj.instrument.(i).(level).tsw = table();
            obj.instrument.(i).(level).fsw = table();
            obj.instrument.(i).(level).bad = table();
            obj.instrument.(i).(level).diw = table();
          case {'bin', 'qc'}
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
          % extract custom properties
          non_empty = find(structfun(@(x) ~isempty(x), obj.instrument.(i).(level)));
          fnam = fieldnames(obj.instrument.(i).(level));
          if ~isempty(non_empty)
            cprop = fieldnames(obj.instrument.(i).(level).(fnam{non_empty(1)}).Properties.CustomProperties);
            if ~isempty(cprop)
              for j = 1:size(cprop, 1)
                obj.instrument.(i).(cprop{j}) = obj.instrument.(i).(level).(fnam{non_empty(1)}).Properties.CustomProperties.(cprop{j});
              end
            end
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
    function update_userselection_bad(obj, filename, user_selection, remove_old, level, channel)
      if nargin < 5
        level = 'qc';
        channel = '';
      elseif nargin < 6
        channel = ['_' level];
      else
        channel = join(['_' level '_' strjoin(channel, '_')],'');
      end
      if isfile(filename)
        % Load file
        if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
          file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
          save(filename, 'file_selection')
        end
        load(filename, 'file_selection');
        if isfield(file_selection, ['bad' channel]) && ~isempty(file_selection.(['bad' channel]))
          if remove_old
            % Remove old (days2run) selections
            sel = min(obj.cfg.days2run) <= file_selection.(['bad' channel])(:,1) & ...
              file_selection.(['bad' channel])(:,1) < max(obj.cfg.days2run) + 1;
            file_selection.(['bad' channel])(sel,:) = [];
          end
          % Add new user selection
          file_selection.(['bad' channel]) = [file_selection.(['bad' channel]); user_selection];
        else
          file_selection.(['bad' channel]) = user_selection;
        end
      else
        file_selection = struct(['bad' channel], user_selection);
      end
      % Save user selection
      save(filename, 'file_selection'); 
    end
  end
end


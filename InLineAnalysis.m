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
      addpath('lib', 'instruments', 'packages', 'packages/addaxis6', ... % 'packages/datetick2_doy', 
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
        
        % Make sure Pass2QC is loaded and not the former Flag function now deprecated
        if isfield(obj.cfg, 'flag') && ~isfield(obj.cfg, 'Pass2QC')
          obj.cfg.Pass2QC = obj.cfg.flag;
        end
        % cfg_flag_buf = obj.cfg.flag;
        % obj.cfg.flag = struct();
        % obj.cfg.flag.skip = cfg_flag_buf.skip;
        % % Set all instruments with default parameters
        % for i = fieldnames(cfg.instruments)'; i = i{1};
        %   obj.cfg.flag.(i) = struct();
        %   for p = fieldnames(cfg_flag_buf.default)'; p = p{1};
        %     switch p
        %       case {'tot', 'filt'}
        %         % Overwrite tot|filt specific parameters
        %         for sp = fieldnames(cfg_flag_buf.default.(p))'; sp = sp{1};
        %           obj.cfg.flag.(i).(p).(sp) = cfg_flag_buf.default.(p).(sp);
        %         end
        %       otherwise
        %         % Overwrite with instrument specific parameters (both filt and tot)
        %         obj.cfg.flag.(i).tot.(p) = cfg_flag_buf.default.(p);
        %         obj.cfg.flag.(i).filt.(p) = cfg_flag_buf.default.(p);
        %     end
        %   end
        % end
        % % Set instruments with specific parameters
        % % Note: orders of the parameters in the cfg files matters
        % for i = fieldnames(cfg.instruments)'; i = i{1};
        %   if isfield(cfg_flag_buf, i)
        %     for p = fieldnames(cfg_flag_buf.(i))'; p = p{1};
        %       switch p
        %         case {'tot', 'filt'}
        %           % Overwrite tot|filt specific parameters
        %           for sp = fieldnames(cfg_flag_buf.(i).(p))'; sp = sp{1};
        %             obj.cfg.flag.(i).(p).(sp) = cfg_flag_buf.(i).(p).(sp);
        %           end
        %         otherwise
        %           % Overwrite with instrument specific parameters (both filt and tot)
        %           obj.cfg.flag.(i).tot.(p) = cfg_flag_buf.(i).(p);
        %           obj.cfg.flag.(i).filt.(p) = cfg_flag_buf.(i).(p);
        %       end
        %     end
        %   end
        % end
        
        % Initialize each instrument
        for i = fieldnames(cfg.instruments)'%; i = i{1};
          switch cfg.instruments.(i{:}).model
            case {'NMEA','SC701','GP32','GPS32Tara','GP32Tara','GPSSC701Tara','GPSSC701','19X','GPS19X','GPSCOMPASSAT','GPS','GPSraw'}
              obj.instrument.(i{:}) = NMEA(cfg.instruments.(i{:}));
            case {'TSG', 'SBE45', 'SBE3845'}
              obj.instrument.(i{:}) = TSG(cfg.instruments.(i{:}));
            case {'atlasTSG' ,'ATLASECRTD'}
              obj.instrument.(i{:}) = atlasTSG(cfg.instruments.(i{:}));
            case {'FTH', 'ADU100', 'ADU200'}
              obj.instrument.(i{:}) = FTH(cfg.instruments.(i{:}));
            case 'ACS'
              obj.instrument.(i{:}) = ACS(cfg.instruments.(i{:}));
            case 'AC9'
              obj.instrument.(i{:}) = AC9(cfg.instruments.(i{:}));
            case 'HBB'
              obj.instrument.(i{:}) = HBB(cfg.instruments.(i{:}));
%             case 'BB3'
%               obj.instrument.(i{:}) = BB3(cfg.instruments.(i{:}));
%             case 'WSCD'
%               obj.instrument.(i{:}) = WSCD(cfg.instruments.(i{:}));
            case {'LISST', 'LISST100X'}
              obj.instrument.(i{:}) = LISST100X(cfg.instruments.(i{:}));
            case 'LISST200X'
              obj.instrument.(i{:}) = LISST200X(cfg.instruments.(i{:}));
            case {'LISSTTau', 'LISSTTAU', 'TAU'}
              obj.instrument.(i{:}) = LISSTTau(cfg.instruments.(i{:}));
            case 'ECO'
              obj.instrument.(i{:}) = ECO(cfg.instruments.(i{:}));
            case {'FL', 'WS3S'}
              obj.instrument.(i{:}) = FL(cfg.instruments.(i{:}));
            case {'BB','BB2','BB3'} % 'BB9' TODO
              obj.instrument.(i{:}) = BB(cfg.instruments.(i{:}));
            case {'CD', 'WSCD', 'SUVF'}
              obj.instrument.(i{:}) = CD(cfg.instruments.(i{:}));
            case 'ALFA'
              obj.instrument.(i{:}) = ALFA(cfg.instruments.(i{:}));
            case {'PAR','QSP2150A','QCR2150A'}
              obj.instrument.(i{:}) = PAR(cfg.instruments.(i{:}));
            otherwise
              error('Instrument not supported: %s.', cfg.instruments.(i{:}).model);
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
      if ~isdatetime(obj.cfg.days2run)
        obj.cfg.days2run = datetime(obj.cfg.days2run, 'datetime', 'datenum');
      end
      for i=obj.cfg.instruments2run%; i = i{1};
        if ~isfield(obj.instrument, i{:})
          error('Instrument2run "%s" does not match any instrument name in the cfg file:  %s', ...
            i{:}, strjoin(fieldnames(obj.instrument), ' / '))
        else
          fprintf('READ RAW: %s\n', i{:});
          obj.instrument.(i{:}).ReadRaw(obj.cfg.days2run, obj.cfg.force_import, true);
        end
      end
    end
    
    function ReadRawDI(obj)
      if ~isdatetime(obj.cfg.days2run)
        obj.cfg.days2run = datetime(obj.cfg.days2run, 'datetime', 'datenum');
      end
      for i=obj.cfg.instruments2run%; i = i{1};
        if  any(strcmp(i{:},obj.cfg.di.skip))
          fprintf('READ DI: Skip %s\n', i{:});
        else
          fprintf('READ DI: %s\n', i{:});
          obj.instrument.(i{:}).ReadRawDI(obj.cfg.days2run, obj.cfg.force_import, true);
        end
      end
    end
    
    function Sync(obj)
      % Note: Run all days loaded (independent of days2run)
      for i=obj.cfg.instruments2run%; i = i{1};
        if  any(strcmp(i{:},obj.cfg.sync.skip))
          fprintf('SYNC: Skip %s\n', i{:});
        else
          fprintf('SYNC: %s\n', i{:});
          obj.instrument.(i{:}).Sync(obj.cfg.sync.delay.(i{:}));
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
      for i=obj.cfg.instruments2run%; i = i{1};
        if any(contains(lower(i{:}),{'ac','bb','lisst','par','alfa','qsp','qcr'}))
          if any(contains(lower(i{:}),'ac'))
            lambda.a = obj.instrument.(i{:}).lambda_a;
            lambda.c = obj.instrument.(i{:}).lambda_c;
            bb_dark = [];
          elseif any(contains(lower(i{:}),'bb3'))
            lambda.bb = obj.instrument.(i{:}).lambda;
            bb_dark = obj.instrument.(instru{contains(lower(instru), 'bb3')}).dark;
          elseif any(contains(lower(i{:}),{'hbb', 'hyperbb'}))
            lambda.bb = obj.instrument.(i{:}).lambda;
            bb_dark = [];
          elseif any(contains(lower(i{:}),'lisst'))
            lambda.theta = obj.instrument.(i{:}).theta;
            bb_dark = [];
          end
          if ~isempty(obj.instrument.(i{:}).(level).fsw)
%             fprintf('Deleting bad values from %s filtered data ...\n', i);
            if any(contains(lower(i{:}),{'ac','bb','lisst'}))
              [obj.instrument.(i{:}).(level).fsw, bad_spec, Nbad] = AutoQC(i{:}, obj.instrument.(i{:}).(level).fsw,...
                lambda, tolerance.filtered, bb_dark, saturation_threshold, 'FSW');
              if strcmp(level, 'raw')
                obj.instrument.(i{:}).raw.bad = [obj.instrument.(i{:}).raw.bad; bad_spec];
              elseif strcmp(level, 'qc')
                obj.instrument.(i{:}).bad.fsw = [obj.instrument.(i{:}).bad.fsw; bad_spec];
              else
                error("Level %s not supported, AutoQC available only on 'raw' and 'qc' levels")
              end
            end
            if any(contains(lower(i{:}),'ac'))
              fprintf('%4.2f%% of absorption and %4.2f%% of attenuation spectra deleted from %s filtered %s data\n',...
                Nbad.a, Nbad.c, i{:}, level);
            elseif any(contains(lower(i{:}),'bb')) || any(contains(lower(i{:}),'lisst'))
              fprintf('\n')
              fprintf('\n')
              fprintf('Data percentage deleted from %s filtered %s data:\n', i{:}, level)
              lineseg = '-------+';
              headseg = '  %i  |';
              dataseg = ' %2.2f%% |';
              lin = '-----------+';
              dat = ' %% deleted |';
              if any(contains(lower(i{:}),'bb'))
                head = '|  lambda  |';
                nn = size(Nbad.bb,2);
              elseif any(contains(lower(i{:}),'lisst'))
                head = '|  theta  |';
                nn = size(Nbad.lisst,2);
              end
              for ii = 1:nn
                lin = [lin lineseg];
                head = [head headseg];
                dat = [dat dataseg];
              end
              lin = [lin '+\n'];
              head = [head '\n'];
              dat = [dat '\n'];
              if any(contains(lower(i{:}),'bb'))
                fprintf(head, lambda.bb)
                fprintf(lin);
                fprintf(dat, Nbad.bb)
              elseif any(contains(lower(i{:}),'lisst'))
                fprintf(head, lambda.theta)
                fprintf(lin);
                fprintf(dat, Nbad.lisst)
              end
            end
          else
            warning('No filtered data loaded: Skip');
          end
          if ~isempty(obj.instrument.(i{:}).(level).tsw)
% %             fprintf('Deleting bad values from %s total data ...\n', i);
            if any(contains(lower(i{:}),{'par','qcr','qsp'}))
              foo = obj.instrument.(i{:}).(level).tsw.(obj.instrument.(i{:}).view.varname)./obj.instrument.(i{:}).scale > 4500 ...
                | obj.instrument.(i{:}).(level).tsw.(obj.instrument.(i{:}).view.varname)./obj.instrument.(i{:}).scale < 0;
              obj.instrument.(i{:}).(level).tsw(foo,:) = [];
              % get percentage of data deleted
              foo = sum(foo) / size(obj.instrument.(i{:}).(level).tsw, 1) * 100;
            elseif any(contains(lower(i{:}),'alfa'))
              alf_ar = table2array(obj.instrument.(i{:}).(level).tsw);
              toqc = repmat(~contains(obj.instrument.(i{:}).(level).tsw.Properties.VariableNames, ...
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
              obj.instrument.(i{:}).(level).tsw = array2table(alf_ar, 'VariableNames', ...
                obj.instrument.(i{:}).(level).tsw.Properties.VariableNames);
              % get percentage of data deleted
              foo = foo / sum(toqc(:));
            end
            if any(contains(lower(i{:}),{'ac','bb','lisst'}))
              [obj.instrument.(i{:}).(level).tsw, bad_spec, Nbad] = AutoQC(i{:}, obj.instrument.(i{:}).(level).tsw,...
                lambda, tolerance.total, bb_dark, saturation_threshold, 'TSW');
              if strcmp(level, 'raw')
                obj.instrument.(i{:}).raw.bad = [obj.instrument.(i{:}).raw.bad; bad_spec];
              elseif strcmp(level, 'qc')
                obj.instrument.(i{:}).bad.tsw = [obj.instrument.(i{:}).bad.tsw; bad_spec];
              else
                error("Level %s not supported, AutoQC available only on 'raw' and 'qc' levels")
              end
            end
            if any(contains(lower(i{:}),'ac'))
              fprintf('%4.2f%% of absorption and %4.2f%% of attenuation spectra deleted from %s total %s data\n',...
                Nbad.a, Nbad.c, i{:}, level);
            elseif any(contains(lower(i{:}),'BB')) || any(contains(lower(i{:}),'LISST'))
              fprintf('\n')
              fprintf('\n')
              fprintf('Data percentage deleted from %s total %s data:\n', i{:}, level)
              lineseg = '-------+';
              headseg = '  %i  |';
              dataseg = ' %2.2f%% |';
              lin = '-----------+';
              dat = ' %% deleted |';
              if any(contains(lower(i{:}),'bb'))
                head = '|  lambda  |';
                nn = size(Nbad.bb,2);
              elseif any(contains(lower(i{:}),'lisst'))
                head = '|  theta  |';
                nn = size(Nbad.lisst,2);
              end
              for ii = 1:nn
                lin = [lin lineseg];
                head = [head headseg];
                dat = [dat dataseg];
              end
              lin = [lin '+\n'];
              head = [head '\n'];
              dat = [dat '\n'];
              if any(contains(lower(i{:}),'bb'))
                fprintf(head, lambda.bb)
                fprintf(lin);
                fprintf(dat, Nbad.bb)
              elseif any(contains(lower(i{:}),'lisst'))
                fprintf(head, lambda.theta)
                fprintf(lin);
                fprintf(dat, Nbad.lisst)
              end
            elseif any(contains(lower(i{:}),{'par','qcr','qsp','alfa'}))
              fprintf('%.3f%% of %s %s values deleted\n', foo, level, i{:});
            end
          else
            warning('No total data loaded: Skip');
          end
          if ~isempty(obj.instrument.(i{:}).(level).diw)
%             fprintf('Deleting bad values from %s dissolved data...\n', i);
            if any(contains(lower(i{:}),'ac'))
              lambda.a = obj.instrument.(i{:}).lambda_a;
              lambda.c = obj.instrument.(i{:}).lambda_c;
            elseif  any(contains(lower(i{:}),'bb'))
              lambda.bb = obj.instrument.(i{:}).lambda;
            end
            if any(contains(lower(i{:}),{'ac','bb'}))
              [obj.instrument.(i{:}).(level).diw, bad_spec, Nbad] = AutoQC(i{:}, obj.instrument.(i{:}).(level).diw,...
                lambda, tolerance.dissolved, bb_dark, saturation_threshold, 'DI');
              if strcmp(level, 'raw')
                obj.instrument.(i{:}).raw.bad = [obj.instrument.(i{:}).raw.bad; bad_spec];
              elseif strcmp(level, 'qc')
                obj.instrument.(i{:}).bad.diw = [obj.instrument.(i{:}).bad.diw; bad_spec];
              else
                error("Level %s not supported, AutoQC available only on 'raw' and 'qc' levels")
              end
            end
            if any(contains(lower(i{:}),'ac'))
              fprintf('%4.2f%% of absorption and %4.2f%% of attenuation spectra deleted from %s dissolved %s data\n',...
                Nbad.a, Nbad.c, i{:}, level);
            elseif any(contains(lower(i{:}),'bb'))
              fprintf('\n')
              fprintf('\n')
              fprintf('Data percentage deletred from %s dissolved %s data:\n', i{:}, level)
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
          fprintf('No AutoQC for %s [Done]\n', i{:})
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
      for i=obj.cfg.instruments2run%; i = i{1};
        if any(contains(i{:},instru))
          fprintf('%s Spectral QCs\n', i{:});
          [user_selection] = SpectralQC(obj.instrument.(i{:}), obj.cfg.days2run, i{:}, level, ...
            save_figure, obj.meta.cruise, toClean);
          % Apply user selection
          if ~isempty(user_selection)
            if ~isfolder(obj.instrument.(i{:}).path.ui); mkdir(obj.instrument.(i{:}).path.ui); end
            obj.instrument.(i{:}).DeleteUserSelection(user_selection, level{:}, toClean);
            if strcmp(level{:}, 'prod')
              obj.instrument.(i{:}).DeleteUserSelection(user_selection, 'qc', {'tsw', strrep(toClean{2}, 'p', '')});
            end
            % Save user selection
            if strcmp(toClean{1}, 'diw')
              filename = fullfile(obj.instrument.(i{:}).path.ui, [i{:} '_QCDI_pickSpecific_UserSelection.mat']);
            else
              filename = fullfile(obj.instrument.(i{:}).path.ui, [i{:} '_QCpickSpecific_UserSelection.mat']);
            end
            obj.update_userselection_bad(filename, user_selection, false, level{:}, toClean);
          end
        end
      end
    end
    
    function visProd_timeseries(obj)
      if ~isdatetime(obj.cfg.days2run)
        obj.cfg.days2run = datetime(obj.cfg.days2run, 'ConvertFrom', 'datenum');
      end
      for i=obj.cfg.instruments2run%; i = i{1};
        if any(~contains(lower(i{:}), {'flow', 'fth'}))
          fprintf('%s products time series plots\n', i{:});
          ifieldn = fieldnames(obj.instrument.(i{:}).prod);
          for j=1:size(ifieldn,1)
            if ~isempty(obj.instrument.(i{:}).prod.(ifieldn{j}))
              if ~isdatetime(obj.instrument.(i{:}).prod.(ifieldn{j}).dt)
                obj.instrument.(i{:}).prod.(ifieldn{j}).dt = datetime(obj.instrument.(i{:}).prod.(ifieldn{j}).dt, 'ConvertFrom', 'datenum');
              end
              id_day = obj.instrument.(i{:}).prod.(ifieldn{j}).dt >= min(obj.cfg.days2run) & obj.instrument.(i{:}).prod.(ifieldn{j}).dt < max(obj.cfg.days2run)+days(1);
              if ~strcmp(ifieldn{j}, 'FiltStat') && ~strcmp(ifieldn{j}, 'QCfailed') && ...
                  ~isempty(obj.instrument.(i{:}).prod.(ifieldn{j}))
                if contains(lower(i{:}), {'bb'})
                  if ~strcmp(ifieldn{j}, 'g') & contains(lower(i{:}), {'hyperbb'}) & ~isempty(obj.cfg.AC_source)
                    wrong_ac = CheckAncillary(obj, obj.cfg.AC_source);
                    % if AC source valid, interpolate Chl, cp, and bp onto hyperBB table
                    if ~wrong_ac
                      % interpolate Chl
                      obj.instrument.(i{:}).prod.(ifieldn{j}).chl_ap676lh = interp1(obj.instrument.(obj.cfg.AC_source).prod.p.dt, ...
                        obj.instrument.(obj.cfg.AC_source).prod.p.chl_ap676lh, obj.instrument.(i{:}).prod.(ifieldn{j}).dt, "linear");
                      % interpolate cp660
                      idc660 = min(abs(obj.instrument.(obj.cfg.AC_source).lambda_c - 660)) == abs(obj.instrument.(obj.cfg.AC_source).lambda_c - 660);
                      obj.instrument.(i{:}).prod.(ifieldn{j}).cp660 = interp1(obj.instrument.(obj.cfg.AC_source).prod.p.dt, ...
                        obj.instrument.(obj.cfg.AC_source).prod.p.cp(:,idc660), obj.instrument.(i{:}).prod.(ifieldn{j}).dt, "linear");
                      % interpolate bp550
                      idc550 = min(abs(obj.instrument.(obj.cfg.AC_source).lambda_c - 550)) == abs(obj.instrument.(obj.cfg.AC_source).lambda_c - 550);
                      ida550 = min(abs(obj.instrument.(obj.cfg.AC_source).lambda_a - 550)) == abs(obj.instrument.(obj.cfg.AC_source).lambda_a - 550);
                      obj.instrument.(i{:}).prod.(ifieldn{j}).bp550 = interp1(obj.instrument.(obj.cfg.AC_source).prod.p.dt, ...
                        obj.instrument.(obj.cfg.AC_source).prod.p.cp(:,idc550) - obj.instrument.(obj.cfg.AC_source).prod.p.ap(:,ida550), ...
                        obj.instrument.(i{:}).prod.(ifieldn{j}).dt, "linear");
                    end
                  end
                  visProd_timeseries(obj.instrument.(i{:}).prod.(ifieldn{j})(id_day, :), i{:}, ...
                    obj.instrument.(i{:}).lambda);
                  % remove inteporlated variables
                  % if any(strcmp(obj.instrument.(i{:}).prod.(ifieldn{j}).Properties.VariableNames, 'chl_ap676lh'))
                  %   obj.instrument.(i{:}).prod.(ifieldn{j}) = removeVars(obj.instrument.(i{:}).prod.(ifieldn{j}), 'chl_ap676lh');
                  % end
                  % if any(strcmp(obj.instrument.(i{:}).prod.(ifieldn{j}).Properties.VariableNames, 'cp550'))
                  %   obj.instrument.(i{:}).prod.(ifieldn{j}) = removeVars(obj.instrument.(i{:}).prod.(ifieldn{j}), 'cp550');
                  % end
                  % if any(strcmp(obj.instrument.(i{:}).prod.(ifieldn{j}).Properties.VariableNames, 'bp550'))
                  %   obj.instrument.(i{:}).prod.(ifieldn{j}) = removeVars(obj.instrument.(i{:}).prod.(ifieldn{j}), 'bp550');
                  % end
                elseif contains(lower(i{:}), {'ac'})
                  visProd_timeseries(obj.instrument.(i{:}).prod.(ifieldn{j})(id_day, :), i{:}, ...
                    obj.instrument.(i{:}).lambda_c);
                else
                  visProd_timeseries(obj.instrument.(i{:}).prod.(ifieldn{j})(id_day, :), i{:});
                end
              end
            end
          end
        end
      end
    end
    
    function Stretch(obj)
      % Note: Run all days loaded (independent of days2run)
      for i=obj.cfg.instruments2run%; i = i{1};
        if  any(strcmp(i{:},obj.cfg.stretch.skip))
          fprintf('STRETCH: Skip %s\n', i{:});
        else
          fprintf('STRETCH: %s\n', i{:});
          obj.instrument.(i{:}).Stretch(obj.cfg.stretch.delta.(i{:}));
        end
      end
    end
    
    function QCRef(obj)
      if contains(lower(obj.cfg.qcref.view), {'par','qcr','qsp'})
        error("QCRef not required for %s data: skip this step and run 'Split'", obj.cfg.qcref.view)
      end
      if isempty(obj.instrument.(obj.cfg.qcref.view).data)
        error('%s data table is empty, make sure to view the instrument you are trying to process or load your data before running QCref', ...
          obj.cfg.qcref.view)
      end
      if ~isdatetime(obj.cfg.days2run)
        obj.cfg.days2run = datetime(obj.cfg.days2run, 'ConvertFrom', 'datenum');
      end
      switch obj.cfg.qcref.mode
        case 'ui'
          if isempty(obj.instrument.FLOW.data)
            error('No flow/switch data loaded.')
          end
          % create new FTH data for missing data
          % round date time to second
          if ~isdatetime(obj.instrument.FLOW.data.dt)
            obj.instrument.FLOW.data.dt = datetime(obj.instrument.FLOW.data.dt, 'ConvertFrom', 'datenum');
          end
          % obj.instrument.FLOW.data.dt = datenum(floor(datevec(obj.instrument.FLOW.data.dt)));
          obj.instrument.FLOW.data.dt = dateshift(obj.instrument.FLOW.data.dt, 'Start', 'second');
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
            obj.instrument.(obj.cfg.qcref.reference).data.(obj.instrument.(obj.cfg.qcref.reference).view.swt_variable), ...
            'k', 'LineWidth', 1);
          ylim([-0.1 1.1]);
          ax = gca; ax.YColor = 'k';
          ylabel('Switch position')
          yyaxis('right');
          plot(obj.instrument.(obj.cfg.qcref.view).data.dt, ...
            obj.instrument.(obj.cfg.qcref.view).data.(obj.instrument.(obj.cfg.qcref.view).view.varname)(:,obj.instrument.(obj.cfg.qcref.view).view.varcol),'.');
          ylim(prctile(obj.instrument.(obj.cfg.qcref.view).data.(obj.instrument.(obj.cfg.qcref.view).view.varname)(:,obj.instrument.(obj.cfg.qcref.view).view.varcol), [0.1 99.9]));
          if contains(obj.cfg.qcref.view, 'AC')
            ylabel([obj.instrument.(obj.cfg.qcref.view).view.varname ' ' ...
                num2str(round(obj.instrument.(obj.cfg.qcref.view).('lambda_ref')(obj.instrument.(obj.cfg.qcref.view).view.varcol),0)) 'nm'])
          elseif contains(obj.cfg.qcref.view, 'BB')
            ylabel([obj.instrument.(obj.cfg.qcref.view).view.varname ' ' ...
                num2str(round(obj.instrument.(obj.cfg.qcref.view).('lambda')(obj.instrument.(obj.cfg.qcref.view).view.varcol),0)) 'nm'])
          else
            ylabel([obj.instrument.(obj.cfg.qcref.view).view.varname obj.instrument.(obj.cfg.qcref.view).view.varcol]);
          end
          % datetick2_doy();
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
            % convert to datetime if datenum
            fnames = fieldnames(file_selection);
            for i = fnames'
              if ~isdatetime(file_selection.(i{:})) && ~isempty(file_selection.(i{:}))
                file_selection.(i{:}) = datetime(file_selection.(i{:}), 'ConvertFrom', 'datenum', 'Format','yyyy-MM-dd HH:mm:ss.SSS');
              end
            end
            if obj.cfg.qcref.remove_old
              % Remove old (days2run) selections
              if ~isempty(file_selection.total)
                sel = min(obj.cfg.days2run) <= file_selection.total(:,1) & file_selection.total(:,1) < max(obj.cfg.days2run) + days(1);
                file_selection.total(sel,:) = [];
              end
              if ~isempty(file_selection.filtered)
                sel = min(obj.cfg.days2run) <= file_selection.filtered(:,1) & file_selection.filtered(:,1) < max(obj.cfg.days2run) + days(1);
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
            % convert to datetime if datenum
            file_upd = false;
            fnames = fieldnames(file_selection);
            for i = fnames'
              if ~isdatetime(file_selection.(i{:})) && ~isempty(file_selection.(i{:}))
                file_selection.(i{:}) = datetime(file_selection.(i{:}), 'ConvertFrom', 'datenum', 'Format','yyyy-MM-dd HH:mm:ss.SSS');
                file_upd = true;
              end
            end
            if file_upd
              save(filename, 'file_selection')
            end
            % Remove selection from days before & after days2run
            if ~isempty(file_selection.total)
              sel = file_selection.total(:,2) < min(obj.cfg.days2run) | max(obj.cfg.days2run) + days(1) < file_selection.total(:,1);
              file_selection.total(sel,:) = [];
            end
            if ~isempty(file_selection.filtered)
              sel = file_selection.filtered(:,2) < min(obj.cfg.days2run) | max(obj.cfg.days2run) + days(1) < file_selection.filtered(:,1);
              file_selection.filtered(sel,:) = [];
            end
            % Apply selection
            if ~isempty(file_selection.total)
              obj.instrument.(obj.cfg.qcref.reference).ApplyUserInput(file_selection.total, 'total');
            end
            if ~isempty(file_selection.filtered)
              obj.instrument.(obj.cfg.qcref.reference).ApplyUserInput(file_selection.filtered, 'filtered');
            end
            fig(31);
            title(['\fontsize{22}Switch position QC:' newline '\fontsize{18}Previous selection applied' newline 'Make sure filtered and total are properly seleceted'], 'interpreter', 'tex');
            fprintf('Select total (t; red) and filtered (f; green) sections (q to save)\n');
            yyaxis('left');
            plot(obj.instrument.(obj.cfg.qcref.reference).data.dt,...
                 obj.instrument.(obj.cfg.qcref.reference).data.(obj.instrument.(obj.cfg.qcref.reference).view.swt_variable), ...
                 'k', 'LineWidth', 1);
            ylim([-0.1 1.1]);
            ax = gca; ax.YColor = 'k';
            ylabel('Switch position')
            yyaxis('right');
            plot(obj.instrument.(obj.cfg.qcref.view).data.dt, ...
              obj.instrument.(obj.cfg.qcref.view).data.(obj.instrument.(obj.cfg.qcref.view).view.varname)(:,obj.instrument.(obj.cfg.qcref.view).view.varcol),'.');
            ylim(prctile(obj.instrument.(obj.cfg.qcref.view).data.(obj.instrument.(obj.cfg.qcref.view).view.varname)(:,obj.instrument.(obj.cfg.qcref.view).view.varcol), [0.1 99.9]));
            if contains(obj.cfg.qcref.view, 'AC')
              ylabel([obj.instrument.(obj.cfg.qcref.view).view.varname ' ' ...
                  num2str(round(obj.instrument.(obj.cfg.qcref.view).('lambda_ref')(obj.instrument.(obj.cfg.qcref.view).view.varcol),0)) 'nm'])
            elseif contains(obj.cfg.qcref.view, 'BB')
              ylabel([obj.instrument.(obj.cfg.qcref.view).view.varname ' ' ...
                  num2str(round(obj.instrument.(obj.cfg.qcref.view).('lambda')(obj.instrument.(obj.cfg.qcref.view).view.varcol),0)) 'nm'])
            else
              ylabel([obj.instrument.(obj.cfg.qcref.view).view.varname obj.instrument.(obj.cfg.qcref.view).view.varcol]);
            end
            % datetick2_doy();
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
      for i=obj.cfg.instruments2run%; i = i{1};
        if isempty(obj.instrument.(i{:}).data)
          error('%s data table is empty', i{:})
        end
        if any(strcmp(i{:}, obj.cfg.split.skip))
          fprintf('SPLIT: Skip %s (copy data to next level)\n', i{:});
          obj.instrument.(i{:}).raw.tsw = obj.instrument.(i{:}).data;
        elseif strcmp(obj.instrument.(i{:}).split.mode, 'None')
          fprintf('SPLIT: Not available for %s\n', i{:});
        else
          fprintf('SPLIT: %s\n', i{:});
          obj.instrument.(i{:}).Split(obj.instrument.(obj.cfg.split.reference), obj.cfg.split.buffer.(i{:}));
        end
      end
    end
    
    function Bin(obj)
      % Note: Run all days loaded (independent of days2run)
      for i=obj.cfg.instruments2run%; i = i{1};
        obj.instrument.(i{:}).bin.tsw = table();
        obj.instrument.(i{:}).bin.fsw = table();
        if  any(strcmp(i{:},obj.cfg.bin.skip))
          fprintf('BIN: Skip %s (copy data to next level)\n', i{:});
          obj.instrument.(i{:}).bin.tsw = obj.instrument.(i{:}).raw.tsw;
          obj.instrument.(i{:}).bin.fsw = obj.instrument.(i{:}).raw.fsw;
        else
          fprintf('BIN: %s\n', i{:});
          obj.instrument.(i{:}).Bin(obj.cfg.bin.bin_size.(i{:}),...
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
      for i=obj.cfg.instruments2run%; i = i{1};
        if  any(strcmp(i{:},obj.cfg.di.skip))
          fprintf('BIN DI: Skip %s\n', i{:});
        else
          fprintf('BIN DI: %s\n', i{:});
          obj.instrument.(i{:}).BinDI(obj.cfg.di.bin.bin_size,...
                                   obj.cfg.bin.prctile_detection,...
                                   obj.cfg.bin.prctile_average,...
                                   obj.cfg.parallel);
%                                    obj.cfg.bin.method,...
                                   
        end
      end
    end

    % Flag is DEPRECATED: Use skip to move data to next level: renamed Pass2QC
    function Flag(obj,level)
      warning('on')
      error('Flag function DEPRECATED: Use command "ila.Pass2QC()" instead to copy data to QC level')
      % warning('off')
    end

    function Pass2QC(obj,level)
      % Pass data to QC level
      % Note: Run all days loaded (independent of days2run)
      % add all instrument from instrument2run into the skip so that Pass2QC never bugs
      if nargin == 1
        level = 'all';
      end
      if isfield(obj.cfg, 'flag') && ~isfield(obj.cfg, 'Pass2QC')
        obj.cfg.Pass2QC.skip = obj.cfg.flag.skip;
      end
      obj.cfg.Pass2QC.skip = [obj.cfg.Pass2QC.skip; obj.cfg.instruments2run(:)];
      for i=obj.cfg.instruments2run%; i = i{1};
        if  any(strcmp(i{:},obj.cfg.Pass2QC.skip))
          fprintf('Pass2QC: Copy %s data to QC level)\n', i{:});
          if strcmp(level, 'particulate')
            obj.instrument.(i{:}).qc.tsw = obj.instrument.(i{:}).bin.tsw;
            obj.instrument.(i{:}).qc.fsw = obj.instrument.(i{:}).bin.fsw;
          elseif strcmp(level, 'dissolved')
            obj.instrument.(i{:}).qc.diw = obj.instrument.(i{:}).raw.diw;
          else
            obj.instrument.(i{:}).qc.tsw = obj.instrument.(i{:}).bin.tsw;
            obj.instrument.(i{:}).qc.fsw = obj.instrument.(i{:}).bin.fsw;
            obj.instrument.(i{:}).qc.diw = obj.instrument.(i{:}).raw.diw;
          end
        end
      end
    end
    
    function QC(obj)
      if ischar(obj.cfg.qc.specific.run)
        obj.cfg.qc.specific.run = {obj.cfg.qc.specific.run};
      end
      if ischar(obj.cfg.qc.global.view)
        obj.cfg.qc.global.view = {obj.cfg.qc.global.view};
      end
      if any(~contains(obj.cfg.instruments2run, {'PAR','QSP','QCR'}))
        % Check if remove data when flow below threshold
        if isfield(obj.cfg.qc, 'remove_when_flow_below')
          if islogical(obj.cfg.qc.remove_when_flow_below)
            if obj.cfg.qc.remove_when_flow_below
              flow_threshold = 0.5;
              rm_with_flow = true;
            else
              rm_with_flow = false;
            end
          elseif isnumeric(obj.cfg.qc.remove_when_flow_below)
            flow_threshold = obj.cfg.qc.remove_when_flow_below;
            rm_with_flow = true;
          else
            rm_with_flow = false;
          end
        else
          rm_with_flow = false;
        end
        % select flow while rounding timestamps
        if ~isempty(obj.instrument.FLOW.raw.tsw)
          flow_level = 'raw';
          obj.instrument.FLOW.(flow_level).tsw = round_timestamp(obj.instrument.FLOW.(flow_level).tsw, seconds(1));
        elseif ~isempty(obj.instrument.FLOW.bin.tsw)
          flow_level = 'bin';
          obj.instrument.FLOW.(flow_level).tsw = round_timestamp(obj.instrument.FLOW.(flow_level).tsw, minutes(1));
        elseif ~isempty(obj.instrument.FLOW.qc.tsw)
          flow_level = 'qc';
          obj.instrument.FLOW.(flow_level).tsw = round_timestamp(obj.instrument.FLOW.(flow_level).tsw, minutes(1));
        else
          warning('Flow data not loaded')
          flow_level = 'raw';
        end
        % remove sections with low flow
        if ~isempty(obj.instrument.FLOW.(flow_level).tsw)
          for i=obj.cfg.qc.specific.run(:)'%; i = i{1};
            if ~any(strcmp(obj.cfg.instruments2run, i{:})); continue; end
            % round timestamps
            if ~isempty(obj.instrument.(i{:}).qc.fsw)
              obj.instrument.(i{:}).qc.fsw = round_timestamp(obj.instrument.(i{:}).qc.fsw, seconds(1));
            end
            if ~isempty(obj.instrument.(i{:}).qc.tsw)
              obj.instrument.(i{:}).qc.tsw = round_timestamp(obj.instrument.(i{:}).qc.tsw, minutes(1));
            end
            if rm_with_flow
              remove_low_flow(obj, i{:}, flow_level, flow_threshold)
            end
          end
          fooflow = obj.instrument.FLOW.(flow_level).tsw;
        else
          fooflow = [];
        end
      else
        fooflow = [];
      end
      % Manual quality check of the data resulting in good and bad data
      switch obj.cfg.qc.mode
        case 'ui'
          if obj.cfg.qc.global.active
            % Display interactive figure
            foo = obj.instrument.(obj.cfg.qc.global.view{:});
            fh=visFlag(foo.raw.tsw, foo.raw.fsw, foo.qc.tsw, foo.suspect.tsw,...
                       foo.qc.fsw, foo.suspect.fsw, foo.view.varname, foo.view.varcol,...
                       foo.raw.bad, fooflow, obj.instrument.FLOW.view.spd_variable);
            title('\fontsize{22}\color{red}Global QC: \fontsize{18}\color{black}Press t to select section to trash (press q to save)', 'interpreter', 'tex');
            fprintf('Global QC: Press t to select section to trash (press q to save)\n');
            user_selection = guiSelectOnTimeSeries(fh);
            % For each instrument 
            for i=obj.cfg.qc.global.apply(:)'%; i = i{1};
              if ~any(strcmp(obj.cfg.instruments2run, i{:})); continue; end
              % Apply user selection
              obj.instrument.(i{:}).DeleteUserSelection(user_selection);
              % Save user selection
              filename = fullfile(fileparts(obj.instrument.(i{:}).path.ui), 'QCGlobal_UserSelection.mat');
              if ~isfolder(fileparts(fileparts(obj.instrument.(i{:}).path.ui)))
                mkdir(fileparts(fileparts(obj.instrument.(i{:}).path.ui)));
              end
              obj.update_userselection_bad(filename, user_selection, obj.cfg.qc.remove_old);
            end
          else
            % Load previous globalQC if file exists and apply
            for i=obj.cfg.qc.global.apply(:)'%; i = i{1};
              if ~any(strcmp(obj.cfg.instruments2run, i{:})); continue; end
              filename = fullfile(fileparts(obj.instrument.(i{:}).path.ui), 'QCGlobal_UserSelection.mat');
              if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
                file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
                save(filename, 'file_selection')
              end
              if isfile(filename)
                fprintf('QC LOAD Global selection: %s ... ', i{:});
                load(filename, 'file_selection');
                % convert to datetime if datenum
                file_upd = false;
                fnames = fieldnames(file_selection);
                for f = fnames'
                  if ~isdatetime(file_selection.(f{:})) && ~isempty(file_selection.(f{:}))
                    file_selection.(f{:}) = datetime(file_selection.(f{:}), 'ConvertFrom', 'datenum', 'Format','yyyy-MM-dd HH:mm:ss.SSS');
                    file_upd = true;
                  end
                end
                if file_upd
                  save(filename, 'file_selection')
                end
                fprintf('done\n')
                obj.instrument.(i{:}).DeleteUserSelection(file_selection.bad);
              end
            end
          end
          if obj.cfg.qc.specific.active
            % For each instrument
            for i=obj.cfg.qc.specific.run(:)'%; i = i{1};
              if ~any(strcmp(obj.cfg.instruments2run, i{:})); continue; end
              if ~isfolder(obj.instrument.(i{:}).path.ui); mkdir(obj.instrument.(i{:}).path.ui); end
              if isempty(obj.instrument.(i{:}).qc.tsw) && ~isempty(obj.instrument.(i{:}).bin.tsw)
                error("Data found in 'bin' level but not in 'qc' level: run 'ila.Pass2QC()' command before 'ila.QC()' command to copy bin data into qc")
              end
              if contains(lower(i{:}), 'ac')
                channel = {'a', 'c'};
              elseif contains(lower(i{:}), {'tsg', 'sbe45', 'sbe3845',})
                channel = {obj.instrument.(i{:}).temperature_variable, obj.instrument.(i{:}).salinity_variable};
              elseif contains(lower(i{:}), {'atlas', 'ecrtd'})
                channel = {obj.instrument.(i{:}).temperature_variable, obj.instrument.(i{:}).conductivity_variable};
              else
                channel = obj.instrument.(i{:}).qc.tsw.Properties.VariableNames(~contains(obj.instrument.(i{:}).qc.tsw.Properties.VariableNames,{'dt','lat','lon','_avg_sd','_avg_n'}));
              end
              if contains(lower(i{:}), {'ac', 'tsg', 'sbe45', 'sbe3845', 'atlas', 'ecrtd'}) && ~obj.cfg.qc.qc_once_for_all
                for j = channel
                  % Display interactive figure
                  fh=visFlag(obj.instrument.(i{:}).raw.tsw, obj.instrument.(i{:}).raw.fsw, ...
                        obj.instrument.(i{:}).qc.tsw, obj.instrument.(i{:}).suspect.tsw, obj.instrument.(i{:}).qc.fsw, ...
                        obj.instrument.(i{:}).suspect.fsw, j{:}, obj.instrument.(i{:}).view.varcol,...
                        obj.instrument.(i{:}).raw.bad, fooflow, obj.instrument.FLOW.view.spd_variable);
                  title(['\fontsize{22}\color{red}' i{:} ' QC of "' j{:} '" only:' newline '\fontsize{18}\color{black}Press t to select section to trash (press q to save)'], 'interpreter', 'tex');
                  fprintf([i{:} ' QC of "' j{:} '" only: Press t to select section to trash (press q to save)\n']);
                  % Get user selection
                  user_selection = guiSelectOnTimeSeries(fh);
                  % Apply user selection
                  obj.instrument.(i{:}).DeleteUserSelection(user_selection, 'qc', ['tsw' j]);
                  obj.instrument.(i{:}).DeleteUserSelection(user_selection, 'qc', ['fsw' j]);
                  % Save user selection
                  filename = fullfile(obj.instrument.(i{:}).path.ui, [i{:} '_QCSpecific_UserSelection.mat']);
                  obj.update_userselection_bad(filename, user_selection, obj.cfg.qc.remove_old, ...
                          'qc', ['tsw' j]);
                  obj.update_userselection_bad(filename, user_selection, obj.cfg.qc.remove_old, ...
                          'qc', ['fsw' j]);
                  clf(52)
                end
              elseif contains(lower(i{:}), 'alfa') && ~obj.cfg.qc.qc_once_for_all
                channel = obj.instrument.(i{:}).qc.tsw.Properties.VariableNames;
                channel(contains(channel, {'dt', '_avg_sd', '_avg_n'})) = [];
                for j = channel
                  % delete crazy values
                  if contains(j{:}, 'WL')
                    obj.instrument.(i{:}).qc.tsw.(j{:})(obj.instrument.(i{:}).qc.tsw.(j{:}) > 1000 | obj.instrument.(i{:}).qc.tsw.(j{:}) < 540) = NaN;
                  else
                    obj.instrument.(i{:}).qc.tsw.(j{:})(obj.instrument.(i{:}).qc.tsw.(j{:}) > 100) = NaN;
                  end
                  fh=visFlag(obj.instrument.(i{:}).raw.tsw, obj.instrument.(i{:}).raw.fsw, ...
                        obj.instrument.(i{:}).qc.tsw, obj.instrument.(i{:}).suspect.tsw, obj.instrument.(i{:}).qc.fsw, ...
                        obj.instrument.(i{:}).suspect.fsw, j{:}, obj.instrument.(i{:}).view.varcol,...
                        obj.instrument.(i{:}).raw.bad, fooflow, obj.instrument.FLOW.view.spd_variable);
                  title(['\fontsize{22}\color{red}' i{:} ' QC of "' j{:} '" only:' newline '\fontsize{18}\color{black}Press t to select section to trash (press q to save)'], 'interpreter', 'tex');
                  fprintf([i{:} ' QC of "' j{:} '" only: Press t to select section to trash (press q to save)\n']);
                  % Get user selection
                  user_selection = guiSelectOnTimeSeries(fh);
                  % Apply user selection
                  obj.instrument.(i{:}).DeleteUserSelection(user_selection, 'qc', ['tsw' j]);
                  % Save user selection
                  filename = fullfile(obj.instrument.(i{:}).path.ui, [i{:} '_QCSpecific_UserSelection.mat']);
                  obj.update_userselection_bad(filename, user_selection, obj.cfg.qc.remove_old, ...
                          'qc', ['tsw' j]);
                  clf(52)
                end
              else
                if ~isempty(obj.instrument.(i{:}).raw.tsw)
                  fh=visFlag(obj.instrument.(i{:}).raw.tsw, obj.instrument.(i{:}).raw.fsw, obj.instrument.(i{:}).qc.tsw, ...
                        obj.instrument.(i{:}).suspect.tsw, obj.instrument.(i{:}).qc.fsw, obj.instrument.(i{:}).suspect.fsw, ...
                        obj.instrument.(i{:}).view.varname, obj.instrument.(i{:}).view.varcol,...
                        obj.instrument.(i{:}).raw.bad, fooflow, obj.instrument.FLOW.view.spd_variable);
                else
                  fh=visFlag([], [], obj.instrument.(i{:}).qc.tsw, obj.instrument.(i{:}).suspect.tsw, ...
                        obj.instrument.(i{:}).qc.fsw, obj.instrument.(i{:}).suspect.fsw,...
                        obj.instrument.(i{:}).view.varname, obj.instrument.(i{:}).view.varcol, ...
                        obj.instrument.(i{:}).raw.bad, fooflow, obj.instrument.FLOW.view.spd_variable);
                end
                title(['\fontsize{22}\color{red}' i{:} ' QC all variables:' newline '\fontsize{18}\color{black}Press t to select section to trash (press q to save)'], 'interpreter', 'tex');
                fprintf([i{:} ' QC all: Press t to select section to trash (press q to save)\n']);
                user_selection = guiSelectOnTimeSeries(fh);
                % Apply user selection
                obj.instrument.(i{:}).DeleteUserSelection(user_selection);
                % Save user selection
                filename = fullfile(obj.instrument.(i{:}).path.ui, [i{:} '_QCSpecific_UserSelection.mat']);
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
          for i=obj.cfg.qc.global.apply(:)'%; i = i{1};
            if ~any(strcmp(obj.cfg.instruments2run, i{:})); continue; end
            fprintf('QC LOAD Global selection: %s ... ', i{:});
            filename = fullfile(fileparts(obj.instrument.(i{:}).path.ui), 'QCGlobal_UserSelection.mat');
            if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
              file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
              save(filename, 'file_selection')
            end
            if isfile(filename)
              load(filename, 'file_selection');
              % convert to datetime if datenum
              file_upd = false;
              fnames = fieldnames(file_selection);
              for f = fnames'
                if ~isdatetime(file_selection.(f{:})) && ~isempty(file_selection.(f{:}))
                  file_selection.(f{:}) = datetime(file_selection.(f{:}), 'ConvertFrom', 'datenum', 'Format','yyyy-MM-dd HH:mm:ss.SSS');
                  file_upd = true;
                end
              end
              if file_upd
                save(filename, 'file_selection')
              end
              fprintf('done\n')
              obj.instrument.(i{:}).DeleteUserSelection(file_selection.bad);
            else
              fprintf(['Warning: ' filename ' not found\n'])
            end
          end
          % end
          if obj.cfg.qc.specific.active
            for i=obj.cfg.qc.specific.run(:)'%; i = i{1};
              if ~any(strcmp(obj.cfg.instruments2run, i{:})); continue; end
              fprintf('QC LOAD Specific slection: %s ... ', i{:});
              filename = fullfile(obj.instrument.(i{:}).path.ui, [i{:} '_QCSpecific_UserSelection.mat']);
              if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
                file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
                save(filename, 'file_selection')
              end
              if isfile(filename)
                fprintf('done\n')
                load(filename, 'file_selection');
                % convert to datetime if datenum
                file_upd = false;
                fnames = fieldnames(file_selection);
                for f = fnames'
                  if ~isdatetime(file_selection.(f{:})) && ~isempty(file_selection.(f{:}))
                    file_selection.(f{:}) = datetime(file_selection.(f{:}), 'ConvertFrom', 'datenum', 'Format','yyyy-MM-dd HH:mm:ss.SSS');
                    file_upd = true;
                  end
                end
                if file_upd
                  save(filename, 'file_selection')
                end
                sel_toload = fieldnames(file_selection);
                for j = 1:size(sel_toload, 1)
                  % Keep only selection of the day2run
                  if ~isempty(file_selection.(sel_toload{j}))
                    file_selection.(sel_toload{j})(all(file_selection.(sel_toload{j}) < min(obj.cfg.days2run)-hours(1), 2) | ...
                      all(file_selection.(sel_toload{j}) > max(obj.cfg.days2run)+days(1)+hours(1), 2), :) = [];
                    if contains(lower(i{:}), 'ac')
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
                      obj.instrument.(i{:}).DeleteUserSelection(file_selection.(sel_toload{j}), level, channel);
                    else
                      obj.instrument.(i{:}).DeleteUserSelection(file_selection.(sel_toload{j}));
                    end
                  end
                end
              else
                fprintf(['Warning: ' filename ' not found\n'])
              end
              fprintf('QC LOAD Specific pick selection: %s ... ', i{:});
              filename = fullfile(obj.instrument.(i{:}).path.ui, [i{:} '_QCpickSpecific_UserSelection.mat']);
              if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
                file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
                save(filename, 'file_selection')
              end
              if isfile(filename)
                % load hand picked bad values
                fprintf('done\n')
                load(filename, 'file_selection');
                % convert to datetime if datenum
                file_upd = false;
                fnames = fieldnames(file_selection);
                for f = fnames'
                  if ~isdatetime(file_selection.(f{:})) && ~isempty(file_selection.(f{:}))
                    file_selection.(f{:}) = datetime(file_selection.(f{:}), 'ConvertFrom', 'datenum', 'Format','yyyy-MM-dd HH:mm:ss.SSS');
                    file_upd = true;
                  end
                  if size(file_selection, 2) > 1
                    file_selection = file_selection';
                    file_upd = true;
                  end
                end
                if file_upd
                  save(filename, 'file_selection')
                end
                sel_picktoload = fieldnames(file_selection);
                for j = progress(1:size(sel_picktoload, 1))
                  % Keep only selection of the day2run
                  if ~isdatetime(file_selection.(sel_picktoload{j}))
                    file_selection.(sel_picktoload{j}) = datetime(file_selection.(sel_picktoload{j}), 'ConvertFrom', 'datenum', 'Format','yyyy-MM-dd HH:mm:ss.SSS');
                  end
                  if ~isempty(file_selection.(sel_picktoload{j}))
                    file_selection.(sel_picktoload{j})(file_selection.(sel_picktoload{j}) < min(obj.cfg.days2run)-hours(1) | ...
                      file_selection.(sel_picktoload{j}) > max(obj.cfg.days2run)+days(1)+hours(1)) = [];
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
                    if contains(lower(i{:}), 'ac')
                      obj.instrument.(i{:}).DeleteUserSelection(file_selection.(sel_picktoload{j}), ...
                        level, channel);
                      if strcmp(level, 'prod')
                        obj.instrument.(i{:}).DeleteUserSelection(file_selection.(sel_picktoload{j}), ...
                          'qc', {'tsw', strrep(channel{2}, 'p', '')});
                      end
                    else
                      obj.instrument.(i{:}).DeleteUserSelection(file_selection.(sel_picktoload{j}), level);
                      if strcmp(level, 'prod')
                        obj.instrument.(i{:}).DeleteUserSelection(file_selection.(sel_picktoload{j}), 'qc');
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
      switch obj.cfg.di.qc.mode
        case 'ui'
          % For each instrument
          for i=obj.cfg.instruments2run%; i = i{1};
            if ~any(strcmp(obj.cfg.instruments2run, i{:})) || any(strcmp(obj.cfg.di.skip, i{:})); continue; end
            % Display interactive figure
            foo = obj.instrument.(i{:}); %(obj.instrument.(i).dt, ;
            if isempty(foo.qc.diw)
              error('Empty qc diw \n');
            end
            % create folder for user input
            filename = fullfile(obj.instrument.(i{:}).path.ui, [i{:} '_QCDI_UserSelection.mat']);
            if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
              file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
              save(filename, 'file_selection')
            end
            if ~isfolder(obj.instrument.(i{:}).path.ui); mkdir(obj.instrument.(i{:}).path.ui); end
            ColorSet = lines(2);
            fh = fig(52); hold('on');
            if contains(lower(i{:}), 'ac') && ~obj.cfg.di.qc.qc_once_for_all  
              channel = {'a', 'c'};
              for j = 1:size(channel,2)
                plot(foo.qc.diw.dt, foo.qc.diw.(channel{j})(:,foo.view.varcol), 'o', 'Color', ColorSet(j,:));
                title([i{:} ' QC of "' channel{j} '" only' newline 'Press t to select section to trash (press q to save)']);
                ylabel(channel{j});
                % datetick2_doy();
                set(datacursormode(fh), 'UpdateFcn', @data_cursor_display_date);
                % Get user selection
                user_selection = guiSelectOnTimeSeries(fh);
                % Apply user selection
                obj.instrument.(i{:}).DeleteUserSelection(user_selection, 'qc', ['diw' channel(j)]);
                % Save user selection
                obj.update_userselection_bad(filename, user_selection, obj.cfg.di.qc.remove_old, ...
                        'qc', ['diw' channel(j)]);
                clf(52)
              end
            else
              plot(foo.qc.diw.dt, foo.qc.diw.(foo.view.varname)(:,foo.view.varcol), 'o', 'Color', ColorSet(1,:));
              ylabel(foo.view.varname);
              title([i{:} ' QC all' newline 'Trash full section pressing t (q to save)']);
              % datetick2_doy();
              set(datacursormode(fh), 'UpdateFcn', @data_cursor_display_date);
              % Get user selection
              user_selection = guiSelectOnTimeSeries(fh);
              % Apply user selection
              obj.instrument.(i{:}).DeleteUserSelection(user_selection);
              % Save user selection
              obj.update_userselection_bad(filename, user_selection, obj.cfg.di.qc.remove_old);
            end
          end
        case 'load'
          % Load previous QC DI files and apply them
          for i=obj.cfg.instruments2run%; i = i{1};
            if ~any(strcmp(obj.cfg.instruments2run, i{:})) || any(strcmp(obj.cfg.di.skip, i{:})); continue; end
            fprintf('QC DI LOAD: %s\n', i{:});
            filename = fullfile(obj.instrument.(i{:}).path.ui, [i{:} '_QCDI_UserSelection.mat']);
            if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
              file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
              save(filename, 'file_selection')
            end
            if isfile(filename)
              % load bad DI values
              load(filename, 'file_selection');
              % convert to datetime if datenum
              file_upd = false;
              fnames = fieldnames(file_selection);
              for f = fnames'
                if ~isdatetime(file_selection.(f{:})) && ~isempty(file_selection.(f{:}))
                  file_selection.(f{:}) = datetime(file_selection.(f{:}), 'ConvertFrom', 'datenum', 'Format','yyyy-MM-dd HH:mm:ss.SSS');
                  file_upd = true;
                end
              end
              if file_upd
                save(filename, 'file_selection')
              end
              sel_toload = fieldnames(file_selection);
              for j = 1:size(sel_toload, 1)
                % Keep only selection of the day2run
                file_selection.(sel_toload{j})(all(file_selection.(sel_toload{j}) < min(obj.cfg.days2run)-hours(1), 2) | ...
                  all(file_selection.(sel_toload{j}) > max(obj.cfg.days2run)+days(1)+hours(1), 2), :) = [];
                if ~isempty(file_selection.(sel_toload{j}))
                  if contains(lower(i{:}), 'ac')
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
                    obj.instrument.(i{:}).DeleteUserSelection(file_selection.(sel_toload{j}), ...
                      level, channel);
                  else
                    obj.instrument.(i{:}).DeleteUserSelection(file_selection.(sel_toload{j}));
                  end
                end
              end
            else
              fprintf(['Warning: ' filename ' not found\n'])
            end
            filename = fullfile(obj.instrument.(i{:}).path.ui, [i{:} '_QCDI_pickSpecific_UserSelection.mat']);
            if all(isfile(strrep(filename, '.mat', '.json')) & ~isfile(filename))
              file_selection = json_to_mat(strrep(filename, '.mat', '.json'));
              save(filename, 'file_selection')
            end
            if isfile(filename)
              % load hand picked bad DI values
              load(filename, 'file_selection');
              % convert to datetime if datenum
              file_upd = false;
              fnames = fieldnames(file_selection);
              for f = fnames'
                if ~isdatetime(file_selection.(f{:})) && ~isempty(file_selection.(f{:}))
                  file_selection.(f{:}) = datetime(file_selection.(f{:}), 'ConvertFrom', 'datenum', 'Format','yyyy-MM-dd HH:mm:ss.SSS');
                  file_upd = true;
                end
                if size(file_selection, 2) > 1
                  file_selection = file_selection';
                  file_upd = true;
                end
              end
              if file_upd
                save(filename, 'file_selection')
              end
              sel_picktoload = fieldnames(file_selection);
              for j = 1:size(sel_picktoload, 1)
                % Keep only selection of the day2run
                file_selection.(sel_picktoload{j})(file_selection.(sel_picktoload{j}) < min(obj.cfg.days2run)-hours(1) | ...
                  file_selection.(sel_picktoload{j}) > max(obj.cfg.days2run)+days(1)+hours(1)) = [];
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
                  if contains(lower(i{:}), 'ac')
                    obj.instrument.(i{:}).DeleteUserSelection(file_selection.(sel_picktoload{j}), ...
                      level, channel);
                    if strcmp(level, 'prod')
                      obj.instrument.(i{:}).DeleteUserSelection(file_selection.(sel_picktoload{j}), ...
                        'qc', {'fsw', strrep(channel{2}, 'g', '')});
                    end
                  else
                    obj.instrument.(i{:}).DeleteUserSelection(file_selection.(sel_picktoload{j}));
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
    
    function FindAncillarySource(obj)
      warning('on','all')
      for i=obj.cfg.instruments2run
        % extract instrument classes
        if contains(lower(i{:}), {'ac','lissttau','lisst-tau','tau','cstar','bb','bb3','hbb'})
          instrument_class = table(fieldnames(obj.instrument), 'VariableNames', {'instruments'});
          instrument_class.class = cell(size(instrument_class.instruments));
          for co = 1:size(instrument_class, 1)
            instrument_class.class{co} = class(obj.instrument.(instrument_class.instruments{co}));
          end
        end
        % Check FDOM and TSG source for 'ac','lissttau','lisst-tau','tau','cstar' instruments
        if contains(lower(i{:}), {'ac','lissttau','lisst-tau','tau','cstar'})
          CheckFDOMSource(obj, instrument_class, i{:})
          CheckTSGSource(obj, instrument_class, i{:})
        end
        % Check TSG and AC source for 'bb','bb3','hbb' instruments, if no
        % AC source check FDOM source and FDOM ag parameters for ag reconstruct from FDOM
        if contains(lower(i{:}), {'bb','bb3','hbb'})
          CheckTSGSource(obj, instrument_class, i{:})
          AC_ok = CheckACSource(obj, instrument_class, i{:});
          % Check FDOM source and fdom_ag parameters if AC source not loaded or no fdom variable in AC source for attenuation correction
          if ~AC_ok
            CheckFDOMSource(obj, instrument_class, i{:})
          end
        end
        % % make sure the right AC source was input, and if not, select the good one
        % if contains(lower(i{:}), {'bb', 'bb3', 'hbb'})
        %   CheckACSource(obj, instrument_class, i{:})
        % end
      end
    end

    function wrong_source = CheckAncillary(obj, instrument_source)
      % find instrument table name
      tblname = fieldnames(obj.instrument.(instrument_source).prod);
      wrong_source = false;
      if isempty(tblname)
        wrong_source = true;
      elseif ~isempty(obj.instrument.(instrument_source).prod.(tblname{1}))
        time_match = obj.instrument.(instrument_source).prod.(tblname{1}).dt >= min(obj.cfg.days2run) & ...
          obj.instrument.(instrument_source).prod.(tblname{1}).dt < max(obj.cfg.days2run) + days(1);
        if ~any(time_match)
          wrong_source = true;
        end
      else
        wrong_source = true;
      end
    end
    
    function CheckFDOMSource(obj, instrument_class, i)
      if ~isempty(obj.cfg.calibrate.(i).CDOM_source)
        wrong_cdom = CheckAncillary(obj, obj.cfg.calibrate.(i).CDOM_source);
        if wrong_cdom
          list_cdom = instrument_class.instruments(strcmp(instrument_class.class, 'CD') | strcmp(instrument_class.class, 'WSCD') | strcmp(instrument_class.class, 'SUVF'));
          wrong_cdom = true;
          c = 1;
          while wrong_cdom && c <= size(list_cdom, 1)
            wrong_cdom = CheckAncillary(obj, list_cdom{c});
            if wrong_cdom && c == size(list_cdom, 1)
              if contains(lower(i), {'ac','lissttau','lisst-tau','tau','cstar'}) && strcmp(obj.cfg.calibrate.(i).interpolation_method, 'CDOM')
                warning('%s processing: no FDOM data available for CDOM interpolation.', i)
              elseif contains(lower(i), {'bb','bb3','hbb'})
                warning('%s processing: no FDOM data available for attenuation correction.', i)
              end
              c = c + 1;
            elseif wrong_cdom
              c = c + 1;
            else
              warning('%s processing: wrong FDOM source selected: %s\n%s selected automatically instead', i, obj.cfg.calibrate.(i).CDOM_source, list_cdom{c})
              obj.cfg.calibrate.(i).CDOM_source = list_cdom{c};
              obj.cfg.CDOM_source = list_cdom{c};
            end
          end
        end
      end
    end

    function CheckTSGSource(obj, instrument_class, i)
      if ~isempty(obj.cfg.calibrate.(i).TSG_source)
        wrong_tsg = CheckAncillary(obj, obj.cfg.calibrate.(i).TSG_source);
        if wrong_tsg
          list_tsg = instrument_class.instruments(strcmp(instrument_class.class, 'TSG'));
          wrong_tsg = true;
          c = 1;
          while wrong_tsg && c <= size(list_tsg, 1)
            wrong_tsg = CheckAncillary(obj, list_tsg{c});
            if wrong_tsg && c == size(list_tsg, 1)
              if contains(lower(i), {'ac','lissttau','lisst-tau','tau','cstar'})
                warning(['%s processing: no TSG data available, T/S correction not applied before removing dissolved part from total.' newline ...
                  'T/S assumed to be constant between filter and total events.' newline ...
                  'Remaining T variation will be corrected using the residual temperature correction\n'], i)
              elseif contains(lower(i), {'bb','bb3','hbb'})
                warning('%s processing: "compute_dissolved" on but no TSG data available for required seawater scattering correction.', i)
              end
              c = c + 1;
            elseif wrong_tsg
              c = c + 1;
            else
              warning('%s processing: wrong TSG source selected: %s\n%s selected automatically instead', i, obj.cfg.calibrate.(i).TSG_source, list_tsg{c})
              obj.cfg.calibrate.(i).TSG_source = list_tsg{c};
              obj.cfg.TSG_source = list_tsg{c};
            end
          end
        end
      end
    end

    function AC_ok = CheckACSource(obj, instrument_class, i)
      if ~isempty(obj.cfg.calibrate.(i).AC_source)
        wrong_ac = CheckAncillary(obj, obj.cfg.calibrate.(i).AC_source);
        if wrong_ac
          list_ac = instrument_class.instruments(strcmp(instrument_class.class, 'ACS') | strcmp(instrument_class.class, 'AC9'));
          wrong_ac = true;
          c = 1;
          while wrong_ac && c <= size(list_ac, 1)
            wrong_ac = CheckAncillary(obj, list_ac{c});
            if wrong_ac
              c = c + 1;
            else
              warning('%s processing: wrong AC source selected: %s\n%s selected automatically instead', i, obj.cfg.calibrate.(i).AC_source, list_ac{c})
              obj.cfg.calibrate.(i).AC_source = list_ac{c};
              obj.cfg.AC_source = list_ac{c};
            end
          end
        end
        AC_ok = false;
        % find instrument table name
        ac_tblname = fieldnames(obj.instrument.(obj.cfg.AC_source).prod);
        if isempty(ac_tblname)
          if ~isprop(obj.instrument.(i), 'fdom_ag_parameters')
            warning('%s processing: no AC data available for attenuation correction.', i)
          elseif isempty(obj.instrument.(i).fdom_ag_parameters)
            warning('%s processing: no AC data available for attenuation correction.', i)
          end
        elseif ~isempty(ac_tblname)
          p_ok = false;
          p_fdom_ok = false;
          g_ok = false;
          g_fdom_ok = false;
          fdom_ag_parameters_ok = false;
          fdom_ok = false;
          if any(strcmp(ac_tblname, 'p'))
            if ~isempty(obj.instrument.(obj.cfg.calibrate.(i).AC_source).prod.p)
              p_ok = true;
              if any(strcmp(obj.instrument.(obj.cfg.calibrate.(i).AC_source).prod.p.Properties.VariableNames, 'fdom'))
                p_fdom_ok = true;
              end
            end
          end
          if any(strcmp(ac_tblname, 'g'))
            if ~isempty(obj.instrument.(obj.cfg.calibrate.(i).AC_source).prod.g)
              g_ok = true;
              if any(strcmp(obj.instrument.(obj.cfg.calibrate.(i).AC_source).prod.g.Properties.VariableNames, 'fdom'))
                g_fdom_ok = true;
              end
            end
          elseif any(strcmp(obj.instrument.(obj.cfg.calibrate.(i).AC_source).prod.p.Properties.VariableNames, 'ag_modelled'))
            g_ok = true;
            if any(strcmp(obj.instrument.(obj.cfg.calibrate.(i).AC_source).prod.p.Properties.VariableNames, 'fdom'))
              g_fdom_ok = true;
            end
          end
          if isprop(obj.instrument.(i), 'fdom_ag_parameters')
            if ~isempty(obj.instrument.(i).fdom_ag_parameters)
              fdom_ag_parameters_ok = true;
            end
          end
          if ~isempty(obj.cfg.calibrate.(i).CDOM_source)
            fdom_tblname = fieldnames(obj.instrument.(obj.cfg.calibrate.(i).CDOM_source).prod);
            if ~isempty(obj.instrument.(obj.cfg.calibrate.(i).CDOM_source).prod.(fdom_tblname{:}))
              fdom_ok = true;
            end
          end
          if p_ok && g_ok && p_fdom_ok && g_fdom_ok
            AC_ok = true;
          elseif p_ok && g_ok && fdom_ok && (~p_fdom_ok || ~g_fdom_ok)
            if ~p_fdom_ok
              obj.instrument.(obj.cfg.AC_source).prod.p = merge_timeseries(...
                obj.instrument.(obj.cfg.AC_source).prod.p, ...
                obj.instrument.(obj.cfg.calibrate.(i).CDOM_source).prod.(fdom_tblname{:}), 'fdom', '', 30);
            end
            if ~g_fdom_ok
              obj.instrument.(obj.cfg.AC_source).prod.g = merge_timeseries(...
                obj.instrument.(obj.cfg.AC_source).prod.g, ...
                obj.instrument.(obj.cfg.calibrate.(i).CDOM_source).prod.(fdom_tblname{:}), 'fdom', '', 30);
            end
          elseif p_ok && fdom_ag_parameters_ok && fdom_ok
            AC_ok = true;
          elseif p_ok && g_ok && ~fdom_ag_parameters_ok && ~fdom_ok
            warning('%s processing: no FDOM data available for attenuation correction. ag will be linearly interpolated between filter events.', i)
          elseif p_ok && ~fdom_ag_parameters_ok && ~fdom_ok
            warning('%s processing: no ag data available for attenuation correction.', i)
          else
            warning('%s processing: no AC nor FDOM data available for attenuation correction.', i)
          end
        end
      end
    end

    function QCSwitchPosition(obj)
      % Check ancillary data instrument source
      FindAncillarySource(obj)
      if all(~strcmp(obj.cfg.qcref.view, obj.cfg.split.skip)) && ~contains(obj.cfg.qcref.view, {'SUVF', 'WSCD'})
        if contains(lower(obj.cfg.qcref.view), 'bb') && isfield(obj.cfg.calibrate.(obj.cfg.qcref.view), 'filt_method')
          if strcmp(obj.cfg.calibrate.(obj.cfg.qcref.view).filt_method, 'exponential_fit')
            [obj.instrument.(obj.cfg.qcref.view), obj.instrument.FLOW] = QCSwitchPosition(obj.instrument.(obj.cfg.qcref.view), obj.instrument.FLOW, obj.cfg.days2run, 'raw');
          else
            [obj.instrument.(obj.cfg.qcref.view), obj.instrument.FLOW] = QCSwitchPosition(obj.instrument.(obj.cfg.qcref.view), obj.instrument.FLOW, obj.cfg.days2run, 'qc');
          end
        elseif contains(lower(obj.cfg.qcref.view), 'ac') && isfield(obj.cfg.calibrate.(obj.cfg.qcref.view), 'interpolation_method') && ~isempty(obj.cfg.CDOM_source)
          if strcmp(obj.cfg.calibrate.(obj.cfg.qcref.view).interpolation_method, 'CDOM') && ~isempty(obj.instrument.(obj.cfg.CDOM_source).qc)
            [obj.instrument.(obj.cfg.qcref.view), obj.instrument.FLOW, obj.instrument.(obj.cfg.CDOM_source)] = ...
              QCSwitchPosition(obj.instrument.(obj.cfg.qcref.view), obj.instrument.FLOW, obj.cfg.days2run, 'qc', [], obj.instrument.(obj.cfg.CDOM_source));
          else
            [obj.instrument.(obj.cfg.qcref.view), obj.instrument.FLOW] = QCSwitchPosition(obj.instrument.(obj.cfg.qcref.view), obj.instrument.FLOW, obj.cfg.days2run, 'qc');
          end
        else
          [obj.instrument.(obj.cfg.qcref.view), obj.instrument.FLOW] = QCSwitchPosition(obj.instrument.(obj.cfg.qcref.view), obj.instrument.FLOW, obj.cfg.days2run, 'qc');
        end
      else
        fprintf('QCSwitchPosition %s: skip\n', obj.cfg.qcref.view)
      end
    end
    
    % Process
    function Calibrate(obj)
      % Check ancillary data instrument source
      FindAncillarySource(obj)
      % Calibrate (Runs all days loaded: independent of days2run)
      for i=obj.cfg.instruments2run%; i = i{1};
        if any(contains(lower(obj.instrument.(i{:}).model), {'ac9','acs','bb','bb3','hbb','tau'}))
          if ~isempty(obj.cfg.calibrate.(i{:}).TSG_source)
            TSG_source = obj.instrument.(obj.cfg.calibrate.(i{:}).TSG_source);
          else
            TSG_source = [];
          end
          if ~isempty(obj.cfg.calibrate.(i{:}).CDOM_source)
            CDOM_source = obj.instrument.(obj.cfg.calibrate.(i{:}).CDOM_source);
          else
            CDOM_source = [];
          end
        end
        if any(contains(lower(obj.instrument.(i{:}).model), {'bb','bb3','hbb'}))
          if ~isempty(obj.cfg.calibrate.(i{:}).AC_source)
            AC_source = obj.instrument.(obj.cfg.calibrate.(i{:}).AC_source);
          else
            AC_source = [];
          end
        end
        if ~isfolder(obj.instrument.(i{:}).path.prod)
          mkdir(obj.instrument.(i{:}).path.prod)
        end
        if any(strcmp(i{:},obj.cfg.calibrate.skip))
          fprintf('CALIBRATE: Skip %s (copy data to next level)\n', i{:});
          obj.instrument.(i{:}).prod.a = obj.instrument.(i{:}).qc.tsw;
        else
          fprintf('CALIBRATE: %s\n', i{:});
          % Calibrate
          switch obj.instrument.(i{:}).model
            case 'AC9'
              obj.instrument.(i{:}).Calibrate(obj.cfg.days2run, ...
                                           obj.cfg.calibrate.(i{:}).compute_dissolved,...
                                           obj.cfg.calibrate.(i{:}).interpolation_method,...
                                           CDOM_source,...
                                           obj.instrument.(obj.cfg.calibrate.(i{:}).FLOW_source),...
                                           obj.cfg.calibrate.(i{:}).di_method, ...
                                           obj.cfg.calibrate.(i{:}).scattering_correction, ...
                                           obj.cfg.calibrate.(i{:}).compute_ad_aphi, ...
                                           TSG_source); % obj.cfg.min_nb_pts_per_cluster, obj.cfg.time_weight_for_cluster
              if ishandle(36)
                savefig(36, fullfile(obj.instrument.(i{:}).path.prod, sprintf('%s_%s_%s-%s_da_dfdom_slopes', obj.meta.cruise, i{:}, ...
                  min(datetime(obj.cfg.days2run, 'Format', 'yyyyMMdd')), max(datetime(obj.cfg.days2run, 'Format', 'yyyyMMdd')))))
              end
              if ishandle(52)
                savefig(52, fullfile(obj.instrument.(i{:}).path.prod, sprintf('%s_%s_%s-%s_interpolation', obj.meta.cruise, i{:}, ...
                  min(datetime(obj.cfg.days2run, 'Format', 'yyyyMMdd')), max(datetime(obj.cfg.days2run, 'Format', 'yyyyMMdd')))))
              end
            case 'ACS'
              obj.instrument.(i{:}).Calibrate(obj.cfg.days2run, ...
                                           obj.cfg.calibrate.(i{:}).compute_dissolved,...
                                           obj.cfg.calibrate.(i{:}).interpolation_method,...
                                           CDOM_source,...
                                           obj.instrument.(obj.cfg.calibrate.(i{:}).FLOW_source),...
                                           obj.cfg.calibrate.(i{:}).di_method, ...
                                           obj.cfg.calibrate.(i{:}).scattering_correction, ...
                                           obj.cfg.calibrate.(i{:}).compute_ad_aphi, ...
                                           TSG_source); % obj.cfg.min_nb_pts_per_cluster, obj.cfg.time_weight_for_cluster
              if ishandle(36)
                savefig(36, fullfile(obj.instrument.(i{:}).path.prod, sprintf('%s_%s_%s-%s_da_dfdom_slopes', obj.meta.cruise, i{:}, ...
                  min(datetime(obj.cfg.days2run, 'Format', 'yyyyMMdd')), max(datetime(obj.cfg.days2run, 'Format', 'yyyyMMdd')))))
              end
              if ishandle(52)
                savefig(52, fullfile(obj.instrument.(i{:}).path.prod, sprintf('%s_%s_%s-%s_interpolation', obj.meta.cruise, i{:}, ...
                  min(datetime(obj.cfg.days2run, 'Format', 'yyyyMMdd')), max(datetime(obj.cfg.days2run, 'Format', 'yyyyMMdd')))))
              end
            case {'BB', 'BB3', 'HBB'}
              obj.instrument.(i{:}).Calibrate(obj.cfg.days2run,...
                                           obj.cfg.calibrate.(i{:}).compute_dissolved,...
                                           TSG_source,...
                                           obj.instrument.(obj.cfg.calibrate.(i{:}).FLOW_source),...
                                           AC_source,...
                                           CDOM_source,...
                                           obj.cfg.calibrate.(i{:}).di_method,...
                                           obj.cfg.calibrate.(i{:}).filt_method)
              if ishandle(52)
                savefig(52, fullfile(obj.instrument.(i{:}).path.prod, sprintf('%s_%s_%s-%s_interpolation', obj.meta.cruise, i{:}, ...
                  min(datetime(obj.cfg.days2run, 'Format', 'yyyyMMdd')), max(datetime(obj.cfg.days2run, 'Format', 'yyyyMMdd')))))
              end
            case {'FL', 'WS3S'}
              obj.instrument.(i{:}).Calibrate(obj.cfg.days2run, ...
                                           obj.cfg.calibrate.(i{:}).compute_dissolved,...
                                           obj.instrument.(obj.cfg.calibrate.(i{:}).FLOW_source),...
                                           obj.cfg.calibrate.(i{:}).di_method,...
                                           obj.cfg.calibrate.(i{:}).filt_method)
              if ishandle(52)
                savefig(52, fullfile(obj.instrument.(i{:}).path.prod, sprintf('%s_%s_%s-%s_interpolation', obj.meta.cruise, i{:}, ...
                  min(datetime(obj.cfg.days2run, 'Format', 'yyyyMMdd')), max(datetime(obj.cfg.days2run, 'Format', 'yyyyMMdd')))))
              end
            case {'CD', 'WSCD', 'SUVF'}
              obj.instrument.(i{:}).Calibrate(obj.cfg.calibrate.(i{:}).compute_dissolved)
            case {'LISST', 'LISST100X', 'LISST100x', 'LISST200X', 'LISST200x'}
              obj.instrument.(i{:}).Calibrate(obj.cfg.days2run, ...
                                           obj.cfg.calibrate.(i{:}).compute_dissolved,...
                                           obj.instrument.(obj.cfg.calibrate.(i{:}).FLOW_source),...
                                           obj.cfg.calibrate.(i{:}).di_method)
              if ishandle(52)
                savefig(52, fullfile(obj.instrument.(i{:}).path.prod, sprintf('%s_%s_%s-%s_interpolation', obj.meta.cruise, i{:}, ...
                  min(datetime(obj.cfg.days2run, 'Format', 'yyyyMMdd')), max(datetime(obj.cfg.days2run, 'Format', 'yyyyMMdd')))))
              end
            case {'LISSTTau','LISSTTAU','LISST-Tau','TAU'}
              obj.instrument.(i{:}).Calibrate(oobj.cfg.days2run, ...
                                           bj.cfg.calibrate.(i{:}).compute_dissolved,...
                                           obj.cfg.calibrate.(i{:}).interpolation_method,...
                                           CDOM_source,...
                                           obj.instrument.(obj.cfg.calibrate.(i{:}).FLOW_source),...
                                           obj.cfg.calibrate.(i{:}).di_method);
              if ishandle(52)
                savefig(52, fullfile(obj.instrument.(i{:}).path.prod, sprintf('%s_%s_%s-%s_interpolation', obj.meta.cruise, i{:}, ...
                  min(datetime(obj.cfg.days2run, 'Format', 'yyyyMMdd')), max(datetime(obj.cfg.days2run, 'Format', 'yyyyMMdd')))))
              end
            otherwise
              obj.instrument.(i{:}).Calibrate()
          end
        end
      end
    end
    
    % Write
    function Write(obj, level, part_or_diw)
      if nargin < 2; level = 'prod'; end
      % for each instrument
      for i=obj.cfg.instruments2run%; i = i{1};
        if  any(strcmp(i{:},obj.cfg.write.skip))
          fprintf('WRITE %s: Skip %s\n', level, i{:});
        else
          fprintf('WRITE %s: %s\n', level, i{:});
          switch obj.cfg.write.mode
            case 'One file'
              % Save all days2run in one file
              obj.instrument.(i{:}).Write([i{:} '_ALL'], obj.cfg.days2run, level, part_or_diw);
            case 'One day one file'
              % Save each day from days2run in independent files
              for d=obj.cfg.days2run
                obj.instrument.(i{:}).Write([i{:} '_' char(datetime(d, 'Format','yyyyMMdd'))], d, level, part_or_diw);
              end
            otherwise
              error('Unknow writing mode.');
          end
        end
      end
    end
    
    % Load
    function Read(obj, level)
      if ~isdatetime(obj.cfg.days2run)
        obj.cfg.days2run = datetime(obj.cfg.days2run, 'ConvertFrom', 'datenum');
      end
      % LoadProducts is renamed to Read on Oct 19, 2018
      if nargin < 2; level = 'prod'; end
      % for each instrument
      for i=obj.cfg.instruments2run(:)'%; i = i{1};
        fprintf('%s %s flushed\n', i{:}, level);
        switch level
          case 'raw'
            obj.instrument.(i{:}).(level).tsw = table();
            obj.instrument.(i{:}).(level).fsw = table();
            obj.instrument.(i{:}).(level).bad = table();
            obj.instrument.(i{:}).(level).diw = table();
          case {'bin', 'qc'}
            obj.instrument.(i{:}).(level).tsw = table();
            obj.instrument.(i{:}).(level).fsw = table();
            obj.instrument.(i{:}).(level).diw = table();
          case 'prod'
            fna = fieldnames(obj.instrument.(i{:}).(level));
            for j = 1:size(fna, 1)
              if ~isempty(fna{j})
                obj.instrument.(i{:}).(level).(fna{j}) = table();
              end
            end
          otherwise
            error('Level unknown')
        end
        if  any(strcmp(i{:},obj.cfg.write.skip))
          fprintf('LOAD: Skip %s\n', i{:});
        else
          day2read = [min(obj.cfg.days2run)-days(1) obj.cfg.days2run max(obj.cfg.days2run)+days(1)];
          fprintf('LOAD: %s\n', i{:});
          switch obj.cfg.write.mode
            case 'One file'
              % Read all days2run in one file
              obj.instrument.(i{:}).Read([i{:} '_ALL'], day2read, level);
            case 'One day one file'
              % Read each day from days2run in independent files
              for d=day2read
                if isstruct(obj.instrument.(i{:}))
                  fnam = fieldnames(obj.instrument);
                  instru_loaded = false(size(fnam, 1),1);
                  for s = 1:size(fnam, 1)
                    if ~isstruct(obj.instrument.(fnam{s}))
                      instru_loaded(s) = true;
                    end
                  end
                  warning('Instrument loaded from cfg file: %s', strjoin(fnam(instru_loaded), ', '))
                  error("%s entered in instruments2run doesn't correspond to any intrument loaded from cfg file", i)
                else
                  obj.instrument.(i{:}).Read([i{:} '*_' char(datetime(d,'Format','yyyyMMdd'))], d, level);
                end
              end
            otherwise
              error('Unknow loading mode.');
          end
          % extract custom properties
          non_empty = find(structfun(@(x) ~isempty(x), obj.instrument.(i{:}).(level)));
          fnam = fieldnames(obj.instrument.(i{:}).(level));
          if ~isempty(non_empty)
            cprop = fieldnames(obj.instrument.(i{:}).(level).(fnam{non_empty(1)}).Properties.CustomProperties);
            if ~isempty(cprop)
              for j = 1:size(cprop, 1)
                obj.instrument.(i{:}).(cprop{j}) = obj.instrument.(i{:}).(level).(fnam{non_empty(1)}).Properties.CustomProperties.(cprop{j});
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
      for i=obj.cfg.instruments2run%; i = i{1};
        sdata = size(obj.instrument.(i{:}).data);
        sraw = size(obj.instrument.(i{:}).raw.tsw) + size(obj.instrument.(i{:}).raw.diw);
        sbin = size(obj.instrument.(i{:}).bin.tsw) + size(obj.instrument.(i{:}).bin.diw);
        sqc = size(obj.instrument.(i{:}).qc.tsw) + size(obj.instrument.(i{:}).qc.diw);
        ssuspect = size(obj.instrument.(i{:}).suspect.tsw) + size(obj.instrument.(i{:}).suspect.diw);
        sbad = size(obj.instrument.(i{:}).bad.tsw) + size(obj.instrument.(i{:}).bad.diw);
        sprod = []; nprod = {};
        if ~isempty(fieldnames(obj.instrument.(i{:}).prod))
          for t=fieldnames(obj.instrument.(i{:}).prod)'%; t = t{1};
            sprod(end+1,:) = size(obj.instrument.(i{:}).prod.(t{1})); %#ok<AGROW>
            nprod{end+1} = t{:}; %#ok<AGROW>
          end
        end
        fprintf('%10s | %5dx%2d | %5dx%2d | %5dx%2d | %5dx%2d | %5dx%2d | %5dx%2d | ',...
                i{:}, sdata(1), sdata(2), sraw(1), sraw(2), sbin(1), sbin(2), sqc(1), sqc(2),...
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

    % Remove sections with low flow
    function remove_low_flow(obj, i, flow_level, flow_threshold)
      % get flow data
      fooflow = obj.instrument.FLOW.(flow_level).tsw;
      if isempty(fooflow)
        error('Trying to delete data when low flow but FLOW data not loaded, either load FLOW data or choose "remove_when_flow_below = false"')
      else
        if ~isdatetime(fooflow.dt)
          flow_dt = datetime(fooflow.dt, 'ConvertFrom', 'datenum');
        else
          flow_dt = fooflow.dt;
        end
        % speed up the process
        if median(diff(flow_dt)) == seconds(1)
          flow_dt = flow_dt(1:60:end);
        end
      end
      % find flow speed variable
      if isfield(obj.instrument.FLOW.view, 'spd_variable')
        spd_var = obj.instrument.FLOW.view.spd_variable;
      else
        if ~any(strcmp(fooflow.Properties.VariableNames, 'spd'))
          spd_vid = contains(fooflow.Properties.VariableNames, 'spd') & ...
            ~contains(fooflow.Properties.VariableNames, 'avg');
          spd_var = fooflow.Properties.VariableNames{spd_vid & any(~isnan(table2array(fooflow)) & ...
            table2array(fooflow) > 0)};
        end
      end
      % ID low flow period if fsw
      if ~isempty(obj.instrument.(i).qc.fsw) && ~isempty(fooflow)
        % FSW
        fprintf('Deleting %s filtered data when flow <= %.1f LPM ... \n', i, flow_threshold)
        interp_flow = interp1(fooflow.dt, fooflow.(spd_var), obj.instrument.(i).qc.fsw.dt, 'previous');
        if ~isdatetime(obj.instrument.(i).qc.fsw.dt)
          obj.instrument.(i).qc.fsw.dt = datetime(obj.instrument.(i).qc.fsw.dt, 'ConvertFrom', 'datenum');
        end
        % var_dt = obj.instrument.(i).qc.fsw.dt;
        % has_flow = false(size(var_dt));
        % for f = progress(1:size(var_dt, 1))
        %   if any(min(abs(var_dt(f) - flow_dt)) < minutes(2))
        %     has_flow(f) = true;
        %   end
        % end
        low_flow = interp_flow <= flow_threshold;
        % low_flow_id_fsw = obj.instrument.(i).qc.fsw.dt(has_flow & low_flow);
        low_flow_id_fsw = obj.instrument.(i).qc.fsw.dt(low_flow);
        obj.instrument.(i).DeleteUserSelection(low_flow_id_fsw, 'qc', ['fsw' {'all'}]);
        fprintf(' done\n')
      end
      if ~isempty(obj.instrument.(i).qc.tsw) && ~isempty(fooflow)
        % TSW
        fprintf('Deleting %s total data when flow <= %.1f LPM ...', i, flow_threshold)
        interp_flow = interp1(fooflow.dt, fooflow.(spd_var), obj.instrument.(i).qc.tsw.dt, 'previous');
        if ~isdatetime(obj.instrument.(i).qc.tsw.dt)
          obj.instrument.(i).qc.tsw.dt = datetime(obj.instrument.(i).qc.tsw.dt, 'ConvertFrom', 'datenum');
        end
        % var_dt = obj.instrument.(i).qc.tsw.dt;
        % has_flow = false(size(var_dt));
        % for f = progress(1:size(var_dt, 1))
        %   if any(min(abs(var_dt(f) - flow_dt)) < minutes(2))
        %     has_flow(f) = true;
        %   end
        % end
        low_flow = interp_flow <= flow_threshold;
        % low_flow_id_tsw = obj.instrument.(i).qc.tsw.dt(has_flow & low_flow);
        low_flow_id_tsw = obj.instrument.(i).qc.tsw.dt(low_flow);
        obj.instrument.(i).DeleteUserSelection(low_flow_id_tsw, 'qc', ['tsw' {'all'}]);
        fprintf(' done\n')
      end
    end
                          

    % Merge products of two instruments (same model)
    function MergeProducts(obj, primary_instrument, secondary_instrument)
      % Data (at prod level) from the secondary instrument is copied to the
      % primary_instrument. The product table of the primary instrument is
      % then sorted by date & time (dt).

      fprintf('MERGING: %s << %s\n', primary_instrument, secondary_instrument);
      % For each product type (particulate, dissoved...)
      for f = fieldnames(obj.instrument.(primary_instrument).prod)'%; f = f{1};
        ns = size(obj.instrument.(secondary_instrument).prod.(f{:}),1);
        obj.instrument.(primary_instrument).prod.(f{:})(end+1:end+ns,:) = ...
          obj.instrument.(secondary_instrument).prod.(f{:});
        obj.instrument.(primary_instrument).prod.(f{:}) = ...
          sortrows(obj.instrument.(primary_instrument).prod.(f{:}));
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
    
    function cfg = ReadCfgM(filename) %#ok<STOUT>
      % Read matlab configuration file (current format)
      run(filename)
    end
  end
  
  methods (Access=private)
    function update_userselection_bad(obj, filename, user_selection, remove_old, level, channel)
      if nargin < 5
        level = 'qc'; %#ok<NASGU>
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
        % convert to datetime if datenum
        fnames = fieldnames(file_selection);
        for f = fnames'
          if ~isdatetime(file_selection.(f{:})) && ~isempty(file_selection.(f{:}))
            file_selection.(f{:}) = datetime(file_selection.(f{:}), 'ConvertFrom', 'datenum', 'Format','yyyy-MM-dd HH:mm:ss.SSS');
          end
        end
        if isfield(file_selection, ['bad' channel]) && ~isempty(file_selection.(['bad' channel]))
          if remove_old
            % Remove old (days2run) selections
            sel = min(obj.cfg.days2run) <= file_selection.(['bad' channel])(:,1) & ...
              file_selection.(['bad' channel])(:,1) < max(obj.cfg.days2run) + days(1);
            file_selection.(['bad' channel])(sel,:) = [];
          end
          % Add new user selection
          if ~isdatetime(file_selection.(['bad' channel]))
            file_selection.(['bad' channel]) = datetime(file_selection.(['bad' channel]), 'ConvertFrom', 'datenum', 'Format','yyyy-MM-dd HH:mm:ss.SSS');
          end
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


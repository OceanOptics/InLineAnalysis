classdef HBB < Instrument
  % HBB Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    hbb_cal = '';
    hbb_tcal = '';
    fdom_ag_parameters = '';
    lambda = [];
    k_exp = [];
    theta = [];
%     muFactors = [];
  end
  
  methods
    function obj = HBB(cfg)
      % HBB Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Change default processing method
      obj.bin_method = 'SB_IN_PRCTL';
      
      % Post initialization
      if isempty(obj.view.varname); obj.view.varname = 'beta'; end

      % Load HBB plaque calibration file
      load_legacy_cal = false;
      if isfield(cfg,'PlaqueCal') && ~isfield(cfg,'hbb_cal')
        cfg.hbb_cal = cfg.PlaqueCal;
        cfg = rmfield(cfg,'PlaqueCal');     
      end
      if isfield(cfg, 'hbb_cal')
        if ~isempty(cfg.hbb_cal) & isfile(cfg.hbb_cal)
          [calfolder, calfname, ext] = fileparts(cfg.hbb_cal);
          if strcmp(ext, '.hbb_cal')
            obj.hbb_cal = Hbb_ReadBinaryCalFile(cfg.hbb_cal);
          elseif isfile(fullfile(calfolder, ['Hbb_Cal_Plaque_SN' cfg.sn '.hbb_cal'])) || isfile(fullfile(calfolder, [calfname '.hbb_cal']))
            error('The Cal Plaque file extension in cfg is ".mat" but a newer format ".hbb_cal" was found. Make sure to use the correct Cal Plaque file.')
          else
            load_legacy_cal = true;
          end
        else
          error('Cal Plaque file not found: %s', cfg.hbb_cal)
        end
      else
        error('Filedname "hbb_cal" not found in HyperBB cfg.')
      end
      % Load HBB temperature calibration file
      load_legacy_tcal = false;
      if isfield(cfg,'TemperatureCal') && ~isfield(cfg,'hbb_tcal')
        cfg.hbb_cal = cfg.TemperatureCal;
        cfg = rmfield(cfg,'TemperatureCal');     
      end
      if isfield(cfg, 'hbb_tcal')
        if ~isempty(cfg.hbb_tcal) & isfile(cfg.hbb_tcal)
          [calfolder, calfname, ext] = fileparts(cfg.hbb_tcal);
          if strcmp(ext, '.hbb_tcal')
            obj.hbb_tcal = Hbb_ReadBinaryTempCalFile(cfg.hbb_tcal);
          elseif isfile(fullfile(calfolder, ['Hbb_Cal_Temp_SN' cfg.sn '.hbb_tcal'])) || isfile(fullfile(calfolder, [calfname '.hbb_tcal']))
            error('The Cal Temp file extension in cfg is ".mat" but a newer format ".hbb_tcal" was found. Make sure to use the correct Cal Temp file.')
          else
            load_legacy_tcal = true;
          end
        else
          error('Cal Temp file not found: %s', cfg.hbb_cal)
        end
      else
        error('Filedname "hbb_tcal" not found in HyperBB cfg.')
      end
      % Convert cal files into new format
      if load_legacy_cal || load_legacy_tcal
        Hbb_ConvertCalibrations(cfg.hbb_cal, cfg.hbb_tcal)
        if load_legacy_cal
          movefile(strrep(cfg.hbb_cal,'.mat','.hbb_cal'), fullfile(calfolder, ['Hbb_Cal_Plaque_SN' cfg.sn '.hbb_cal']))
          cfg.hbb_cal = fullfile(calfolder, ['Hbb_Cal_Plaque_SN' cfg.sn '.hbb_cal']);
          obj.hbb_cal = Hbb_ReadBinaryCalFile(cfg.hbb_cal);
        else
          delete(strrep(cfg.hbb_cal,'.mat','.hbb_cal'))
        end
        if load_legacy_tcal
          movefile(strrep(cfg.hbb_tcal,'.mat','.hbb_tcal'), fullfile(calfolder, ['Hbb_Cal_Temp_SN' cfg.sn '.hbb_tcal']))
          cfg.hbb_tcal = fullfile(calfolder, ['Hbb_Cal_Temp_SN' cfg.sn '.hbb_tcal']);
          obj.hbb_tcal = Hbb_ReadBinaryTempCalFile(cfg.hbb_tcal);
        else
          delete(strrep(cfg.hbb_tcal,'.mat','.hbb_tcal'))
        end
      end

      if isfield(cfg, 'fdom_ag_correlation')
        [~, ~, ext] = fileparts(cfg.fdom_ag_correlation);
        if any(strcmp(ext, {'.csv', 'xlsx', 'xls'}))
          obj.fdom_ag_parameters = readtable(cfg.fdom_ag_correlation);
        elseif strcmp(ext, '.mat')
          foo = load(cfg.fdom_ag_correlation);
          fname = fieldnames(foo);
          obj.fdom_ag_parameters = foo.(fname{1});
        else
          error('File extension not supported: %s', cfg.fdom_ag_correlation)
        end
        % check variable names
        if ~all(any(strcmpi(obj.fdom_ag_parameters.Properties.VariableNames, 'wl')) & ...
            any(strcmpi(obj.fdom_ag_parameters.Properties.VariableNames, 'slope')) & ...
            any(strcmpi(obj.fdom_ag_parameters.Properties.VariableNames, 'intercept')))
          error("Cannot find all necessary variables, fdom_ag parameters table must contain variable: 'wl', 'slope', and 'intercept'")
        end
      else
        obj.fdom_ag_parameters = [];
        warning('Missing field fdom_ag_correlation for attenuation correction.');
      end
    
      if isfield(cfg, 'theta'); obj.theta = cfg.theta;
      else; error('Missing field theta.'); end
      if isfield(cfg, 'k_exp'); obj.k_exp = cfg.k_exp; % HBB pathlength 
      else; obj.k_exp = 0.1058; fprintf("Missing field pathlength 'k_exp' in *.cfg, using HyperBB's default k_exp = 0.1058.\n"); end
%       if isfield(cfg, 'muFactors'); obj.muFactors = cfg.muFactors;
%       else; error('Missing field muFactors.'); end
      if isempty(obj.logger)
        fprintf('WARNING: Logger set to InlininoHBB.\n');
        obj.logger = 'InlininoHBB';
      end
    end
    
    function ReadRaw(obj, days2run, force_import, write)
      % Get wavelengths from calibration file
      % create wk directory if doesn't exist
      if ~isfolder(obj.path.wk); mkdir(obj.path.wk); end
      % display calibrations
      DisplayCalibration(obj, days2run)
      % Read raw data
      switch obj.logger
        case 'InlininoHBB'
          obj.data = iRead(@importInlininoHBB, obj.path.raw, obj.path.wk, ['HyperBB' obj.sn '_'], days2run, ...
            'Inlinino', force_import, ~write, true, true, '', Inf, obj.hbb_cal, obj.hbb_tcal);
        otherwise
          error('HBB: Unknown logger.');
      end
      % get lambda from data imported directly in case calibration wasn't recorded with same lambda
      obj.lambda = obj.data.Properties.CustomProperties.lambda;
    end
    
    function ReadRawDI(obj, days2run, force_import, write)
      % Get wavelengths from calibration file
      % Set default parameters
      if isempty(obj.path.di)
        fprintf('WARNING: DI Path is same as raw.\n');
        obj.path.di = obj.path.raw;
      end
      if isempty(obj.di_cfg.logger)
        fprintf('WARNING: DI Logger set to InlininoHBB.\n');
        obj.di_cfg.logger = 'InlininoHBB';
      end
      if isempty(obj.di_cfg.postfix)
        fprintf('WARNING: DI Postfix isempty \n');
%         fprintf('WARNING: DI Postfix set to "_DI" \n'); DEPRECATED with Inlinino
%         obj.di_cfg.postfix = '_DI'; DEPRECATED with Inlinino
      end
      if isempty(obj.di_cfg.prefix); obj.di_cfg.prefix = ['HyperBB' obj.sn '_']; end
      % display calibrations
      DisplayCalibration(obj, days2run)
      switch obj.di_cfg.logger
        case 'InlininoHBB'
          obj.raw.diw = iRead(@importInlininoHBB, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'Inlinino', force_import, ~write, true, true, ...
                         obj.di_cfg.postfix, Inf, obj.hbb_cal, obj.hbb_tcal);
        otherwise
          error('HBB: Unknown logger.');
      end
    end

    function DisplayCalibration(obj, days2run)
      % Display selected calibration based on days2run datetime
      idhbb_cal = find(min(days2run) > cell2mat({obj.hbb_cal.date})',1,'last'):find(max(days2run) > cell2mat({obj.hbb_cal.date})',1,'last');
      % plot calibration mufactors
      figure(600); clf; hold on
      title(sprintf('HyperBB%s calibration history (selected calibration(s) in bold)', obj.sn))
      % convert dates to numeric
      numeric_dates = posixtime(cell2mat({obj.hbb_cal.date})');
      % normalize dates between 0 and 1
      norm_dates = (numeric_dates - min(numeric_dates)) / (max(numeric_dates) - min(numeric_dates));
      % select a base colormap and extract colors at exact normalized date positions
      base_colormap = colormap(jet(256));
      % linearly interpolate colors matching date spacing
      col = interp1(linspace(0, 1, 256), base_colormap, norm_dates);
      for c = 1:size(obj.hbb_cal,2)
        if any(c == idhbb_cal)
          plot(obj.hbb_cal(c).muFactorWl,obj.hbb_cal(c).muFactors,'Color',col(c,:),'LineWidth',3)
        else
          plot(obj.hbb_cal(c).muFactorWl,obj.hbb_cal(c).muFactors,'Color',col(c,:))
        end
      end
      cb = colorbar;
      clim([min(numeric_dates), max(numeric_dates)]); % lock color limits to timestamps
      % place ticks based on calibration dates
      cb.Ticks = numeric_dates;
      cb.TickLabels = string(datetime(cell2mat({obj.hbb_cal.date})','Format','yyyy-MM-dd'));
      % write selected calibration in bold and larger font
      ax = cb.Ruler;
      for i = 1:length(ax.TickLabels)
        if any(i == idhbb_cal)
          ax.TickLabels{i} = ['\bf\fontsize{12}' ax.TickLabels{i}]; 
        end
      end
      cb.Label.String = 'Calibration dates (selected calibration(s) in bold)';
      cb.Label.FontSize = 12;
      xlabel('\lambda')
      ylabel('\mu factor')
      drawnow
    end
    
    function Calibrate(obj, days2run, compute_dissolved, TSG, SWT, AC, CDOM, di_method, filt_method)
      SWT_constants = struct('SWITCH_FILTERED', SWT.SWITCH_FILTERED, 'SWITCH_TOTAL', SWT.SWITCH_TOTAL);
      param = struct('lambda', obj.lambda, 'theta', obj.theta, 'k_exp', ...
        obj.k_exp, 'fdom_ag_parameters', obj.fdom_ag_parameters);
%       param = struct('lambda', obj.lambda, 'theta', obj.theta, 'muFactors', obj.muFactors);
      % ProcessHBB
      if compute_dissolved
        switch filt_method
          case '25percentil'
            [obj.prod.p, obj.prod.g] = processHBB(param, obj.qc.tsw, obj.qc.fsw, [], [], ...
              obj.bin.diw, TSG, di_method, filt_method, SWT, SWT_constants, AC, CDOM, days2run);
          case 'exponential_fit'
            [obj.prod.p, obj.prod.g, obj.prod.FiltStat] = processHBB(param, obj.qc.tsw, ...
              obj.qc.fsw, obj.raw.fsw, obj.raw.bad, obj.bin.diw, TSG, ...
              di_method, filt_method, SWT, SWT_constants, AC, CDOM, days2run);
        end
      else
        switch filt_method
          case '25percentil'
            obj.prod.p = processHBB(param, obj.qc.tsw, obj.qc.fsw, [], [], [], TSG, [], ...
              filt_method, SWT, SWT_constants, AC, CDOM, days2run);
          case 'exponential_fit'
            [obj.prod.p, obj.prod.g, obj.prod.FiltStat] = processHBB(param, obj.qc.tsw, obj.qc.fsw, ...
              obj.raw.fsw, obj.raw.bad, [], TSG, [], filt_method, SWT, SWT_constants, AC, CDOM, days2run);
        end
      end
    end
  end
end
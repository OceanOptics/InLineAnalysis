classdef HBB < Instrument
  % HBB Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    PlaqueCal = '';
    TemperatureCal = '';
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

      % Load HBB calibration files and lambda
      if isfield(cfg, 'PlaqueCal')
        % obj.PlaqueCal = cfg.PlaqueCal;
        load(cfg.PlaqueCal, 'cal');
        obj.PlaqueCal = cal;
      else
        error('Missing field PlaqueCal.')
      end

      if isfield(cfg, 'TemperatureCal')
        load(cfg.TemperatureCal, 'cal_temp');
        obj.TemperatureCal = cal_temp;
        obj.lambda = cal_temp.wl;
      else
        error('Missing field TemperatureCal.');
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
      % Read raw data
      switch obj.logger
        case 'InlininoHBB'
          obj.data = iRead(@importInlininoHBB, obj.path.raw, obj.path.wk, ['HyperBB' obj.sn '_'], days2run, ...
            'Inlinino', force_import, ~write, true, true, '', Inf, obj.PlaqueCal, obj.TemperatureCal);
        otherwise
          error('HBB: Unknown logger.');
      end
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
      switch obj.di_cfg.logger
        case 'InlininoHBB'
          obj.raw.diw = iRead(@importInlininoHBB, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'Inlinino', force_import, ~write, true, true, ...
                         obj.di_cfg.postfix, Inf, obj.PlaqueCal, obj.TemperatureCal);
        otherwise
          error('HBB: Unknown logger.');
      end
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
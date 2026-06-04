classdef BB < ECO
  %BB Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    fdom_ag_parameters = '';
    lambda = NaN;
    k_exp = NaN;
    theta = NaN;
  end
  
  methods
    function obj = BB(cfg)
      %BB Construct an instance of this class
      
      % Object Initilization
      obj = obj@ECO(cfg);
      
      % Post initialization
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
          error("%s: Cannot find all necessary variables, fdom_ag parameters table must contain variable: 'wl', 'slope', and 'intercept'", [obj.model obj.SN])
        end
      else
        obj.fdom_ag_parameters = [];
        warning('%s: Missing field fdom_ag_correlation for attenuation correction.', [obj.model obj.sn]);
      end

      if isempty(obj.varname); obj.varname = 'beta'; end % Required for ECO class
      if isfield(cfg, 'lambda'); obj.lambda = cfg.lambda;
      else; error('%s: Missing field lambda.', [obj.model obj.sn]); end
      if isfield(cfg, 'k_exp'); obj.k_exp = cfg.k_exp; % HBB pathlength 
      else; obj.k_exp = 0.01635; fprintf("%s: Missing field pathlength 'k_exp' in *.cfg, using BB3's default k_exp = 0.01635.\n", [obj.model obj.sn]); end
      if isfield(cfg, 'theta'); obj.theta = cfg.theta;
      else; error('%s: Missing field theta.', [obj.model obj.sn]); end
      if isnan(obj.slope); error('%s: Missing field slope.', [obj.model obj.sn]); end
      if isnan(obj.dark); error('%s: Missing field dark.', [obj.model obj.sn]); end
      
    end

    function Calibrate(obj, days2run, compute_dissolved, TSG, SWT, AC, CDOM, di_method, filt_method)
      SWT_constants = struct('SWITCH_FILTERED', SWT.SWITCH_FILTERED, 'SWITCH_TOTAL', SWT.SWITCH_TOTAL);
      param = struct('lambda', obj.lambda, 'theta', obj.theta, 'slope', obj.slope, 'dark', ...
        obj.dark, 'k_exp', obj.k_exp, 'fdom_ag_parameters', obj.fdom_ag_parameters);
      % ProcessBB3
      if compute_dissolved
        switch filt_method
          case '25percentil'
            [obj.prod.p, obj.prod.g] = processBB3(param, obj.qc.tsw, obj.qc.fsw, [], [], ...
              obj.bin.diw, TSG, di_method, filt_method, SWT, SWT_constants, AC, CDOM, days2run);
          case 'exponential_fit'
            [obj.prod.p, obj.prod.g, obj.prod.FiltStat] = processBB3(param, obj.qc.tsw, ...
              obj.qc.fsw, obj.raw.fsw, obj.raw.bad, obj.bin.diw, TSG, di_method, ...
              filt_method, SWT, SWT_constants, AC, CDOM, days2run);
        end
      else
        switch filt_method
          case '25percentil'
            obj.prod.p = processBB3(param, obj.qc.tsw, obj.qc.fsw, [], [], [], TSG, [], ...
              filt_method, SWT, SWT_constants, AC, CDOM, days2run);
          case 'exponential_fit'
            [obj.prod.p, obj.prod.g, obj.prod.FiltStat] = processBB3(param, obj.qc.tsw, obj.qc.fsw, ...
              obj.raw.fsw, obj.raw.bad, [], TSG, [], filt_method, SWT, SWT_constants, AC, CDOM, days2run);
        end
      end
    end
  end
end
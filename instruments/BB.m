classdef BB < ECO
  %BB Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    lambda = NaN;
    theta = NaN;
  end
  
  methods
    function obj = BB(cfg)
      %BB Construct an instance of this class
      
      % Object Initilization
      obj = obj@ECO(cfg);
      
      % Post initialization
      if isempty(obj.varname); obj.varname = 'beta'; end % Required for ECO class
      if isfield(cfg, 'lambda'); obj.lambda = cfg.lambda;
      else; error('Missing field lambda.'); end
      if isfield(cfg, 'theta'); obj.theta = cfg.theta;
      else; error('Missing field theta.'); end
      if isnan(obj.slope); error('Missing field slope.'); end
      if isnan(obj.dark); error('Missing field dark.'); end
      
    end

    function Calibrate(obj, compute_dissolved, TSG, SWT, di_method)
      SWT_constants = struct('SWITCH_FILTERED', SWT.SWITCH_FILTERED, 'SWITCH_TOTAL', SWT.SWITCH_TOTAL);
      param = struct('lambda', obj.lambda, 'theta', obj.theta, 'slope', obj.slope, 'dark', obj.dark);
      % linear interpolation only, CDOM interpolation is not yet available
      if compute_dissolved
        [obj.prod.p, obj.prod.g] = processBB3(param, obj.qc.tsw, obj.qc.fsw, ...
          obj.bin.diw, TSG.qc.tsw, di_method, SWT.qc.tsw, SWT_constants);
      else
        [obj.prod.p] = processBB3(param, obj.qc.tsw, obj.qc.fsw, [], [], [], ...
          SWT.qc.tsw, SWT_constants);
      end
    end
  end
end
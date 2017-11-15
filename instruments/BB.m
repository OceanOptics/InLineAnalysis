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
      if isempty(obj.varname); obj.varname = 'beta'; end
      if isfield(cfg, 'lambda'); obj.lambda = cfg.lambda;
      else; error('Missing field lambda.'); end
      if isfield(cfg, 'theta'); obj.theta = cfg.theta;
      else; error('Missing field theta.'); end
      if isnan(obj.slope); error('Missing field slope.'); end
      if isnan(obj.dark); error('Missing field dark.'); end
      
    end

    function Calibrate(obj, compute_dissolved, TSG)
      param = struct('lambda', obj.lambda, 'theta', obj.theta, 'slope', obj.slope);
      % linear interpolation only, CDOM interpolation is not yet available
      if compute_dissolved
        [obj.prod.p, obj.prod.g] = processBB3(param, obj.qc.tsw, obj.qc.fsw, obj.qc.diw, TSG);
      else
        [obj.prod.p] = processBB3(param, obj.qc.tsw, obj.qc.fsw, [], TSG);
      end
    end
  end
end
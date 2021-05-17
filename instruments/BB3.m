classdef BB3 < Instrument
  %BB3 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    lambda = NaN(1,3);
    theta = NaN;
    slope = NaN(1,3);
    dark = NaN(1,3);
  end
  
  methods
    function obj = BB3(cfg)
      %BB3 Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Change default processing method
      obj.bin_method = 'SB_IN_PRCTL';
      
      % Post initialization
      if isfield(cfg, 'lambda'); obj.lambda = cfg.lambda;
      else; error('Missing field lambda.'); end
      if isfield(cfg, 'theta'); obj.theta = cfg.theta;
      else; error('Missing field theta.'); end
      if isfield(cfg, 'slope'); obj.slope = cfg.slope;
      else; error('Missing field slope.'); end
      if isfield(cfg, 'dark'); obj.dark = cfg.dark;
      else; warning('Missing field dark'); end
      
      if isempty(obj.logger); obj.logger = 'Inlinino'; end
      warning('BB3 Instrument is deprecated.');
    end
    
    function ReadRaw(obj, days2run, force_import, write)
      switch obj.logger
        case 'Inlinino'
          obj.data = iRead(@importInlinino, obj.path.raw, obj.path.wk, 'Inlinino_',...
                         days2run, 'Inlinino', force_import, ~write, true);
          % Remove fdom field if present in Inlinino
          if any(strcmp(obj.data.Properties.VariableNames, 'fdom'))
            obj.data.fdom = [];
          end
          % Remove empty lines
          obj.data(all(isnan(obj.data.beta),2),:) = [];
        otherwise
          error('BB3: Unknown logger.');
      end
    end

    function Calibrate(obj)
      param = struct('lambda', obj.lambda, 'theta', obj.theta, ...
        'slope', obj.slope, 'dark', obj.dark);
      [obj.prod.p, obj.prod.g] = processBB3(param, obj.qc.tsw, obj.qc.fsw, obj.qc.diw, obj.qc.tsg);
    end
  end
end
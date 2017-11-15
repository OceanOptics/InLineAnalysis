classdef ACS < Instrument
  %ACS Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    lambda_ref = [];
    lambda_a = [];
    lambda_c = [];
  end
  
  methods
    function obj = ACS(cfg)
      %ACS Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Post initialization
      if isfield(cfg, 'lambda_a'); obj.lambda_a = cfg.lambda_a;
      else; error('Missing field lambda_a.'); end
      if isfield(cfg, 'lambda_c'); obj.lambda_c = cfg.lambda_c;
      else; error('Missing field lambda_c.'); end
      if isfield(cfg, 'lambda_reference'); obj.lambda_ref = cfg.lambda_reference;
      else; error('Missing field lambda_reference.'); end
      if isempty(obj.logger) 
        fprintf('WARNING: Logger set to Compass_2.1rc_scheduled.\n');
        obj.logger = 'Compass_2.1rc_scheduled';
      end
    end    
    
    function ReadRaw(obj, days2run, force_import, write)
      switch obj.logger
        case 'Compass_2.1rc_scheduled'
          obj.data = iRead(@importACS, obj.path.raw, obj.path.wk, ['acs' obj.sn '_'],...
                         days2run, 'Compass_2.1rc_scheduled', force_import, ~write, true);
        otherwise
          error('ACS: Unknown logger.');
      end
    end
    
    function Calibrate(obj, compute_dissolved, interpolation_method, CDOM, SWT)
      lambda = struct('ref', obj.lambda_ref, 'a', obj.lambda_a, 'c', obj.lambda_c);
      SWT_constants = struct('SWITCH_FILTERED', SWT.SWITCH_FILTERED, 'SWITCH_TOTAL', SWT.SWITCH_TOTAL);
      switch interpolation_method
        case 'linear'
          if compute_dissolved
            [obj.prod.p, obj.prod.g] = processACS(lambda, obj.qc.tsw, obj.qc.fsw, obj.qc.diw);
          else
            [obj.prod.p] = processACS(lambda, obj.qc.tsw, obj.qc.fsw);
          end
        case 'CDOM'
          if compute_dissolved
            [obj.prod.p, obj.prod.g] = processACS(lambda, obj.qc.tsw, obj.qc.fsw, obj.qc.diw, CDOM.qc.tsw, SWT.qc.tsw, SWT_constants);
          else
            [obj.prod.p] = processACS(lambda, obj.qc.tsw, obj.qc.fsw, [], CDOM.qc.tsw, SWT.qc.tsw, SWT_constants);
          end
        otherwise
          error('Method not supported.');
      end
    end
  end
end
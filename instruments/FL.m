classdef FL < ECO
  %FL Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
  end
  
  methods
    function obj = FL(cfg)
      %FL Construct an instance of this class
      
      % Object Initilization
      obj = obj@ECO(cfg);
      
      % Post initialization
      if isempty(obj.varname); obj.varname = 'fchl'; end % Required for ECO class
      if isnan(obj.slope); error('Missing field slope.'); end
      if isnan(obj.dark); error('Missing field dark.'); end
      
    end

    function Calibrate(obj, SWT)
      SWT_constants = struct('SWITCH_FILTERED', SWT.SWITCH_FILTERED, 'SWITCH_TOTAL', SWT.SWITCH_TOTAL);
      param = struct('slope', obj.slope, 'dark', obj.dark);
      obj.prod.p = processFL(param, obj.qc.tsw, obj.qc.fsw, SWT, SWT_constants);
    end
  end
end
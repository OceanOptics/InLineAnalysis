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
      if isempty(obj.varname); obj.varname = 'fchl'; end
      if isnan(obj.slope); error('Missing field slope.'); end
      if isnan(obj.dark); error('Missing field dark.'); end
      
    end

    function Calibrate(obj)
      param = struct('slope', obj.slope);
      obj.prod.p = processFL(param, obj.qc.tsw, obj.qc.fsw);
    end
  end
end
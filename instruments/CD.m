classdef CD < ECO
  %CD Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
  end
  
  methods
    function obj = CD(cfg)
      %CD Construct an instance of this class
      
      % Object Initilization
      obj = obj@ECO(cfg);
      
      % Post initialization
      if isempty(obj.varname); obj.varname = 'fdom'; end % Required for ECO class
      if isempty(obj.view.varname); obj.view.varname = 'fdom'; end
      
      % Change default Split method
      obj.split.mode = 'rmBuffer';
    end

    function Calibrate(obj)
      param = struct('slope', obj.slope, 'dark', obj.dark);
      obj.prod.pd = processCD(param, obj.qc.tsw);
    end
  end
end
classdef PAR < Instrument
  %PAR Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    scale = '';
    dark = '';
  end
  
  methods
    function obj = PAR(cfg)
      %PAR Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Post initialization
      if isfield(cfg, 'logger'); obj.logger = cfg.logger;
      else; error('Missing field logger.'); end
      if isfield(cfg, 'scale'); obj.scale = cfg.scale;
      else; error('Missing field scale.'); end
      if isfield(cfg, 'dark'); obj.dark = cfg.dark;
      else; error('Missing field dark.'); end
      
      % Change default Split method
      obj.split.mode = 'None';
    end
    
    function ReadRaw(obj, days2run, force_import, write)
      switch obj.logger
        case 'Inlinino'
          obj.data = iRead(@importInlininoPAR, obj.path.raw, obj.path.wk, 'Inlinino_',...
                         days2run, 'Inlinino', force_import, ~write, true);
        otherwise
          error('PAR: Unknown logger.');
      end
    end
    
    function Calibrate(obj)
      param = struct('scale', obj.scale, 'dark', obj.dark);
      obj.prod.a = processPAR(param, obj.qc.tsw);
    end
    
  end
end
classdef NMEA < Instrument
  %NMEA Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
  end
  
  methods
    function obj = NMEA(cfg)
      %TSG Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Post initialization
      if isfield(cfg, 'prefix')
        obj.prefix = cfg.prefix;
      else
        obj.prefix = [obj.model obj.sn];
        warning('Missing field "prefix": Default used "[model serial_number]".');
      end
      
      % Change default Split method
      obj.split.mode = 'None';
    end
    
    function ReadRaw(obj, days2run, force_import, write)
      switch obj.logger
        case {'Inlinino'}
            obj.data = iRead(@importInlininoNMEA, obj.path.raw, obj.path.wk, [obj.prefix '_'],...
                           days2run, 'Inlinino', force_import, ~write, true);
        otherwise
          error('NMEA: Unknown logger.');
      end
    end
  end
end


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
      
      % Change default Split method
      obj.split.mode = 'None';
    end
    
    function ReadRaw(obj, days2run, force_import, write)
      switch obj.logger
        case {'Inlinino'}
            obj.data = iRead(@importInlininoGPSSC701, obj.path.raw, obj.path.wk, 'SC701Tara_',...
                           days2run, 'Inlinino', force_import, ~write, true);
        otherwise
          error('NMEA: Unknown logger.');
      end
    end
  end
end


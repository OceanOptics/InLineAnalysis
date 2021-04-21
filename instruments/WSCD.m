classdef WSCD < Instrument
  %WSCD Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
  end
  
  methods
    function obj = WSCD(cfg)
      %WSCD Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Post initialization
      if isempty(obj.logger); obj.logger = 'Inlinino'; end
      
      % Change default processing methods
      obj.split.mode = 'rmBuffer';
      obj.bin_method = 'SB_IN_PRCTL';
      
      warning('WSCD Instrument is deprecated.');
    end
    
    function ReadRaw(obj, days2run, force_import, write)
      switch obj.logger
        case 'Inlinino'
          obj.data = iRead(@importInlinino, obj.path.raw, obj.path.wk, 'Inlinino_',...
                         days2run, 'Inlinino', force_import, ~write, true);
          % Remove beta field if present in Inlinino
          if any(strcmp(obj.data.Properties.VariableNames, 'beta'))
            obj.data.beta = [];
          end
          % Remove empty lines
          obj.data(isnan(obj.data.fdom),:) = [];
        otherwise
          error('WSCD: Unknown logger.');
      end
    end
  end
end
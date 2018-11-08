classdef ALFA < Instrument
  
  properties
    
  end
  
  methods
    function obj = ALFA(cfg)
      %ALFA Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Post initialization
      if isempty(obj.view.varname); obj.view.varname = 'FvFm'; end
      if isempty(obj.logger); obj.logger = 'TeraTerm'; end
    end    
    
    function ReadRaw(obj, days2run, force_import, write)
      switch obj.logger
        case 'ALFA_LabView_m'
          obj.data = iRead(@importALFAm, obj.path.raw, obj.path.wk, ['ALFA' obj.sn '_'],...
                         days2run, 'ALFA_LabView_m', force_import, ~write, true);
        otherwise
          error('ALFA: Unknown logger.');
      end
    end
    
    function ReadRawDI(obj, days2run, force_import, write)
      if isempty(obj.path.di)
        fprintf('WARNING: DI Path is same as raw.\n');
        obj.path.di = obj.path.raw;
      end
      if isempty(obj.di_cfg.postfix) 
        fprintf('WARNING: DI Postfix set to "_DI" \n');
        obj.di_cfg.postfix = '_DI';
      end
      switch obj.logger
        case 'ALFA_LabView_m'
          obj.raw.diw = iRead(@importALFATeraTerm, obj.path.di, obj.path.wk, ['ALFA' obj.sn '_'],...
                         days2run, 'ALFA_LabView_m', force_import, ~write, true, false, obj.di_cfg.postfix);
        otherwise
          error('ALFA: Unknown logger.');
      end
    end

    function Calibrate(obj)
      fprintf('Not Implemented');
    end

  end
end
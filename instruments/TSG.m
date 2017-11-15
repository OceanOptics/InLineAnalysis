classdef TSG < Instrument
  %TSG Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    boat = '';
  end
  
  methods
    function obj = TSG(cfg)
      %TSG Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Post initialization
      if isfield(cfg, 'boat'); obj.boat = cfg.boat;
      else; error('Missing field boat.'); end
      
      % Change default Split method
      obj.split.mode = 'None';
    end
    
    function ReadRaw(obj, days2run, force_import, write)
      switch obj.boat
        case 'Atlantis'
          obj.data = iRead(@importAtlantisTSG, obj.path.raw, obj.path.wk, 'AT',...
                         days2run, 'AtlantisTSG', force_import, ~write, true);
        case 'Pourquoi Pas ?'
          obj.data = iRead(@importPourquoiPasTSG, obj.path.raw, obj.path.wk, 'PP',...
                         days2run, 'PourquoiPasTSG', force_import, ~write, true);
        otherwise
          error('TSG: Unknown boat.');
      end
    end
    
  end
end


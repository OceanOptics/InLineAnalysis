classdef ECO < Instrument
  %ECO Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    varname = '';
    slope = NaN;
    dark = NaN;
  end
  
  methods
    function obj = ECO(cfg)
      %ECO Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Post initialization
      if isfield(cfg, 'varname'); obj.varname = cfg.varname; end
      if isfield(cfg, 'slope'); obj.slope = cfg.slope; end
      if isfield(cfg, 'dark'); obj.dark = cfg.dark; end
      
      if isempty(obj.logger); obj.logger = 'Inlinino'; end
    end
    
    function ReadRaw(obj, days2run, force_import, write)
      switch obj.logger
        case 'Inlinino'
          obj.data = iRead(@importInlinino, obj.path.raw, obj.path.wk, 'Inlinino_',...
                         days2run, 'Inlinino', force_import, ~write, true);
        case 'InlininoPourquoiPas'
          obj.data = iRead(@importInlininoPourquoiPas, obj.path.raw, obj.path.wk, 'Inlinino_',...
                         days2run, 'Inlinino', force_import, ~write, true);
        otherwise
          error('ECO: Unknown logger.');
      end
      % Keep only varname field
      foo = find(strcmp(obj.data.Properties.VariableNames, obj.varname));
      if ~isempty(foo); obj.data = obj.data(:,[1,foo]); end
      % Remove empty lines
      obj.data(all(isnan(obj.data.(obj.varname)),2),:) = [];
    end

  end
end
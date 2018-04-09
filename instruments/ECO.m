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
        case 'DH4PreProc'
          obj.data = iRead(@importDH4PreProc, obj.path.raw, obj.path.wk, 'DH4_bb_CDOM_day',...
                         days2run, 'DH4PreProc', force_import, ~write, true);
        otherwise
          error(['ECO: Unknown logger: ' obj.logger]);
      end
      if ~isempty(obj.data)
        % Keep only varname field
        foo = find(strcmp(obj.data.Properties.VariableNames, obj.varname));
        if ~isempty(foo); obj.data = obj.data(:,[1,foo]); end
        % Remove empty lines
        obj.data(all(isnan(obj.data.(obj.varname)),2),:) = [];
      end
    end
    
    function ReadDI(obj, days2run, force_import, write)
      % DI parameters
      if isempty(obj.path.di)
        fprintf('WARNING: DI Path is same as raw.\n');
        obj.path.di = obj.path.raw;
      end
      if isempty(obj.di_cfg.prefix) 
        fprintf('WARNING: DI Prefix set to "DI_" \n');
        obj.di_cfg.prefix = 'DI_';
      end
      % Other parameters
      switch obj.logger
        case 'Inlinino'
          obj.raw.diw = iRead(@importInlinino, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'Inlinino', force_import, ~write, true, false);
        case 'InlininoPourquoiPas'
          obj.raw.diw = iRead(@importInlininoPourquoiPas, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'Inlinino', force_import, ~write, true, false);
        case 'DH4PreProc'
          obj.raw.diw = iRead(@importDH4PreProc, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'DH4PreProc', force_import, ~write, true, false);
        otherwise
          error(['ECO: Unknown logger: ' obj.logger]);
      end
      if ~isempty(obj.data)
        % Keep only varname field
        foo = find(strcmp(obj.data.Properties.VariableNames, obj.varname));
        if ~isempty(foo); obj.data = obj.data(:,[1,foo]); end
        % Remove empty lines
        obj.data(all(isnan(obj.data.(obj.varname)),2),:) = [];
      end
    end

  end
end
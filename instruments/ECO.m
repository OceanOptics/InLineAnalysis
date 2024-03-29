classdef ECO < Instrument
  %ECO Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    varname = '';
    slope = NaN;
    dark = NaN;
    analog_channel = '';
  end
  
  methods
    function obj = ECO(cfg)
      %ECO Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Change default processing method
      obj.bin_method = 'SB_IN_PRCTL';
      
      % Post initialization
      if isfield(cfg, 'varname'); obj.varname = cfg.varname; end
      if isfield(cfg, 'slope'); obj.slope = cfg.slope; end
      if isfield(cfg, 'dark'); obj.dark = cfg.dark; end
      if isfield(cfg, 'analog_channel'); obj.analog_channel = strrep(strrep(cfg.analog_channel, '(', ''), ')', ''); end
      
      if isempty(obj.logger); obj.logger = 'Inlinino'; end
    end
    
    function ReadRaw(obj, days2run, force_import, write)
      switch obj.logger
        case 'Inlinino'
          obj.data = iRead(@importInlinino, obj.path.raw, obj.path.wk, 'Inlinino_',...
                         days2run, 'Inlinino', force_import, ~write, true);
        case 'Inlinino_base'
          obj.data = iRead(@importInlinino_base, obj.path.raw, obj.path.wk, [obj.prefix '_'],...
                         days2run, 'Inlinino', force_import, ~write, true);
        case 'InlininoBB3'
          obj.data = iRead(@importInlininoBB3, obj.path.raw, obj.path.wk, 'Inlinino_',...
                         days2run, 'Inlinino', force_import, ~write, true);
        case 'InlininoBB3SN'
          obj.data = iRead(@importInlininoBB3, obj.path.raw, obj.path.wk, ['BB3' obj.sn '_'],...
                         days2run, 'Inlinino', force_import, ~write, true);
        case 'InlininoWSCD'
          obj.data = iRead(@importInlininoWSCD, obj.path.raw, obj.path.wk, 'Inlinino_',...
                         days2run, 'Inlinino', force_import, ~write, true);
        case 'InlininoWSCDSN'
          obj.data = iRead(@importInlininoWSCD, obj.path.raw, obj.path.wk, ['WSCD' obj.sn '_'],...
                         days2run, 'Inlinino', force_import, ~write, true);
        case 'InlininoSUVFSN'
          obj.data = iRead(@importInlininoSUVF, obj.path.raw, obj.path.wk, ['SUVF' obj.sn '_'],...
                         days2run, 'Inlinino', force_import, ~write, true);
        case 'InlininoADU100'
          obj.data = iRead(@importInlinino_base, obj.path.raw, obj.path.wk, obj.prefix,...
                         days2run, 'Inlinino', force_import, ~write, true); % @importInlininoFlowControl InlininoADU100
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
        if ~isempty(obj.analog_channel)
          obj.data = renamevars(obj.data, obj.analog_channel, obj.varname);
        end
        % Keep only varname field
        foo = find(strcmp(obj.data.Properties.VariableNames, obj.varname));
        if ~isempty(foo); obj.data = obj.data(:,[1,foo]); end
        % Remove empty lines
        obj.data(all(isnan(table2array(obj.data(:, ~strcmp(obj.data.Properties.VariableNames, 'dt')))), 2), :) = [];
      end
    end
    
    function ReadRawDI(obj, days2run, force_import, write)
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
                         days2run, 'Inlinino', force_import, ~write, true, true, obj.di_cfg.postfix);
        case {'InlininoBB3', 'InlininoBB3SN'}
          obj.raw.diw = iRead(@importInlininoBB3, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'Inlinino', force_import, ~write, true, true, obj.di_cfg.postfix);
        case {'InlininoWSCD', 'InlininoWSCDSN'}
          obj.raw.diw = iRead(@importInlininoWSCD, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'Inlinino', force_import, ~write, true, true, obj.di_cfg.postfix);
        case 'InlininoSUVFSN'
          obj.raw.diw = iRead(@importInlininoSUVF, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'Inlinino', force_import, ~write, true, true, obj.di_cfg.postfix);
        case 'InlininoPourquoiPas'
          obj.raw.diw = iRead(@importInlininoPourquoiPas, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'Inlinino', force_import, ~write, true, true, obj.di_cfg.postfix);
        case 'DH4PreProc'
          obj.raw.diw = iRead(@importDH4PreProc, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'DH4PreProc', force_import, ~write, true, true, obj.di_cfg.postfix);
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
classdef TSG < Instrument
  %TSG Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    boat = '';
    temperature_variable = '';
  end
  
  methods
    function obj = TSG(cfg)
      %TSG Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Post initialization
      if isfield(cfg, 'boat'); obj.boat = cfg.boat;
      else; error('Missing field boat.'); end
      if isfield(cfg, 'temperature_variable'); obj.temperature_variable = cfg.temperature_variable;
      else; error('Missing field temperature variable.'); end
      
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
        case 'Tara'
          switch obj.logger
            case {'Matlab', 'matlab', 'MATLAB'}
              obj.data = iRead(@importTaraTSG, obj.path.raw, obj.path.wk, 'tsg_',...
                             days2run, 'MatlabTSG', force_import, ~write, true);
            case {'TeraTerm', 'teraterm', 'TERATERM'}
              obj.data = iRead(@importTSGTeraTerm, obj.path.raw, obj.path.wk, 'TeraTerm_tsg_',...
                             days2run, 'TeraTerm', force_import, ~write, true);
            case 'Inlinino_base'
              obj.data = iRead(@importInlinino_base, obj.path.raw, obj.path.wk, obj.prefix,...
                             days2run, 'Inlinino', force_import, ~write, true);
            % case {'Inlinino','inlinino'}
            %   obj.data = iRead(@importInlininoTSG, obj.path.raw, obj.path.wk, [obj.model obj.sn '_'],...
            %                  days2run, 'Inlinino', force_import, ~write, true);
%               obj.data = renamevars(obj.data, obj.temperature_variable, 'sst');
%               if strcmp(obj.view.varname, obj.temperature_variable)
%                 obj.view.varname = 'sst';
%               end
%               obj.temperature_variable = 'sst';
%               other_temp = obj.data.Properties.VariableNames{strcmp(obj.data.Properties.VariableNames, 't1') | ...
%                 strcmp(obj.data.Properties.VariableNames, 't2') & ~strcmp(obj.data.Properties.VariableNames, obj.temperature_variable)};
%               obj.data = renamevars(obj.data, other_temp, 't_cal');
%               obj.data = renamevars(obj.data, 's', 'sss');
%               if strcmp(obj.data.Properties.VariableNames, 'c1')
%                 obj.data = renamevars(obj.data, 'c1', 'c');
%               end
%             case {'Inlinino'}
%               obj.data = iRead(@importInlinino_miniTSG, obj.path.raw, obj.path.wk, 'T01_',...
%                              days2run, 'Inlinino', force_import, ~write, true);
            otherwise
              error('TSG: Unknown logger for %s.', obj.boat);
          end
        case 'RRevelle'
          obj.data = iRead(@importRRevelleUnderway, obj.path.raw, obj.path.wk, '',...
                         days2run, 'RRevelleUnderway', force_import, ~write, true);
        case 'IRA-C'
          obj.data = iRead(@importSBE37_TSG, obj.path.raw, obj.path.wk, 'tsg_',...
                         days2run, 'MatlabTSG', force_import, ~write, true);
        case 'JCook'
          switch obj.logger
            case {'SBE45software'}
              obj.data = iRead(@importSBEformatTSG, obj.path.raw, obj.path.wk, '',...
                             days2run, 'SBE45software', force_import, ~write, true);
            case {'matlab_Emmanuel'}
              obj.data = iRead(@importJCookTSG, obj.path.raw, obj.path.wk, '',...
                             days2run, 'matlab_Emmanuel', force_import, ~write, true,...
                             [], 'GPS_TSG');
            otherwise
              error('TSG: Unknown logger for %s.', obj.boat);
          end                           
        otherwise
          error('TSG: Unknown boat.');
      end
    end
    
  end
end


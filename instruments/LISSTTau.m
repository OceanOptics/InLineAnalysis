classdef LISSTTau < Instrument
  %NMEA Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
  end
  
  methods
    function obj = LISSTTau(cfg)
      %TSG Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Post initialization
      
      % Change default Split method
%       obj.split.mode = 'None';
    end
    
    function ReadRaw(obj, days2run, force_import, write)
      switch obj.logger
        case {'Inlinino'}
            obj.data = iRead(@importInlininoTAU, obj.path.raw, obj.path.wk, [obj.prefix '_'],...
                           days2run, 'Inlinino', force_import, ~write, true);
        otherwise
          error('LISST-Tau: Unknown logger.');
      end
    end
    
    function ReadRawDI(obj, days2run, force_import, write)
      % Set default parameters
      if isempty(obj.path.di)
        fprintf('WARNING: DI Path is same as raw.\n');
        obj.path.di = obj.path.raw;
      end
      if isempty(obj.di_cfg.logger)
        fprintf('WARNING: DI Logger set to Inlinino.\n');
        obj.di_cfg.logger = 'Inlinino';
      end
      if isempty(obj.di_cfg.postfix)
        fprintf('WARNING: DI Postfix isempty\n');
%         fprintf('WARNING: DI Postfix set to "_DI" \n'); DEPRECATED with Inlinino
%         obj.di_cfg.postfix = '_DI'; DEPRECATED with Inlinino
      end
      if isempty(obj.di_cfg.prefix); obj.di_cfg.prefix = ['LISSTTau' obj.sn '_']; end
      switch obj.di_cfg.logger
        case 'Inlinino'
          obj.raw.diw = iRead(@importInlininoTAU, obj.path.di, obj.path.wk, obj.di_cfg.prefix, ...
                         days2run, 'Inlinino', force_import, ~write, true, false, obj.di_cfg.postfix);
        otherwise
          error('LISST-Tau: Unknown logger.');
      end
        end
    
    function Calibrate(obj, compute_dissolved, interpolation_method, CDOM, SWT, di_method)
      SWT_constants = struct('SWITCH_FILTERED', SWT.SWITCH_FILTERED, 'SWITCH_TOTAL', SWT.SWITCH_TOTAL);
      switch interpolation_method
        case 'linear'
          if compute_dissolved
            [obj.prod.p, obj.prod.g] = processTAU(obj.qc.tsw, obj.qc.fsw, ...
              obj.bin.diw, [], SWT, SWT_constants, interpolation_method, di_method);
          else
            [obj.prod.p, ~] = processTAU(obj.qc.tsw, obj.qc.fsw, [], [], SWT, ...
              SWT_constants, interpolation_method, []);
          end
        case 'CDOM'
          if ~isfield(CDOM.prod, 'pd') && isempty(CDOM.qc.tsw)
            error('No CDOM data loaded');
          end
          if isempty(CDOM.qc.tsw)
            cdom = CDOM.prod.pd;
          else
            cdom = CDOM.qc.tsw;            
          end
          if compute_dissolved
            [obj.prod.p, obj.prod.g] = processTAU(obj.qc.tsw, obj.qc.fsw, obj.bin.diw, ...
              cdom, SWT, SWT_constants, interpolation_method, di_method);
          else
            [obj.prod.p, ~] = processTAU(obj.qc.tsw, obj.qc.fsw, [], cdom, SWT, ...
              SWT_constants, interpolation_method, []);
          end
        otherwise
          error('Method not supported.');
      end
    end
  end
end


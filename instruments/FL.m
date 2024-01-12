classdef FL < ECO
  %FL Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
  end
  
  methods
    function obj = FL(cfg)
      %FL Construct an instance of this class
      
      % Object Initilization
      obj = obj@ECO(cfg);
      
      % Post initialization
      if isempty(obj.varname); obj.varname = 'fchl'; end % Required for ECO class
      if isempty(obj.view.varname); obj.view.varname = 'fdom'; end
      if isfield(cfg, 'analog_channel'); obj.analog_channel = strrep(strrep(cfg.analog_channel, '(', ''), ')', ''); end
      if isnan(obj.slope); warning('Missing field slope.'); end
      if isnan(obj.dark); warning('Missing field dark.'); end
      
    end

    % function Calibrate(obj)
    %   param = struct('slope', obj.slope, 'dark', obj.dark);
    %   [obj.prod.p, obj.prod.g] = processFL(param, obj.qc.tsw, obj.qc.fsw, obj.qc.diw);
    % end


    function Calibrate(obj, days2run, compute_dissolved, SWT, di_method, filt_method)
      SWT_constants = struct('SWITCH_FILTERED', SWT.SWITCH_FILTERED, 'SWITCH_TOTAL', SWT.SWITCH_TOTAL);
      param = struct('slope', obj.slope, 'dark', obj.dark);
      % linear interpolation only, CDOM interpolation is not yet available (TODO)
      if compute_dissolved
        switch filt_method
          case '25percentil'
            [obj.prod.p, obj.prod.g] = processFL(param, obj.qc.tsw, obj.qc.fsw, [], [], ...
              obj.bin.diw, di_method, filt_method, SWT, SWT_constants, days2run);
          case 'exponential_fit'
            [obj.prod.p, obj.prod.g, obj.prod.FiltStat] = processFL(param, obj.qc.tsw, ...
              obj.qc.fsw, obj.raw.fsw, obj.raw.bad, obj.bin.diw, di_method, ...
              filt_method, SWT, SWT_constants, days2run);
        end
      else
        switch filt_method
          case '25percentil'
            obj.prod.p = processFL(param, obj.qc.tsw, obj.qc.fsw, [], [], [], [], ...
              filt_method, SWT, SWT_constants, days2run);
          case 'exponential_fit'
            [obj.prod.p, ~, obj.prod.FiltStat] = processFL(param, obj.qc.tsw, ...
              obj.qc.fsw, obj.raw.fsw, obj.raw.bad, [], [], filt_method, SWT, ...
              SWT_constants, days2run);
        end
      end
    end

  end
end
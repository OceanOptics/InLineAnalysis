classdef LISST200X < Instrument
  % HBB Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    zsc = [];
    fzs = [];
    dcal = [];
    config = [];
    housek = [];
    inversion = '';
  end
  
  methods
    function obj = LISST200X(cfg)
      % HBB Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Change default processing method
      obj.bin_method = 'SB_IN_PRCTL';
      
      % Post initialization
      if isfield(cfg, 'inversion'); obj.inversion = cfg.inversion;
      else; error('Missing field "inversion".'); end
%       if isempty(obj.view.varname); obj.view.varname = 'beta'; end
%       if isfield(cfg, 'calfile_plaque'); obj.calfile_plaque = cfg.calfile_plaque;
%       else; error('Missing field calfile_plaque.'); end
%       if isfield(cfg, 'calfile_temp'); obj.calfile_temp = cfg.calfile_temp;
%       else; error('Missing field calfile_temp.'); end
%       if isfield(cfg, 'theta'); obj.theta = cfg.theta;
%       else; error('Missing field theta.'); end
% %       if isfield(cfg, 'muFactors'); obj.muFactors = cfg.muFactors;
% %       else; error('Missing field muFactors.'); end

      if isempty(obj.logger)
        fprintf('WARNING: Logger set to "internal_logger".\n');
        obj.logger = 'internal_logger';
      end
    end
    
    function ReadRaw(obj, days2run, force_import, write)
      % create wk directory if doesn't exist
      if ~isfolder(obj.path.wk); mkdir(obj.path.wk); end
      % Read raw data
      switch obj.logger
        case 'internal_logger'
          obj.data = iRead(@importLISST200x_raw, obj.path.raw, obj.path.wk, obj.prefix,...
                         days2run, 'internal_logger', force_import, ~write, true, true, '', Inf);
        otherwise
          error('LISST200X: Unknown logger.');
      end
      obj.zsc = obj.data.Properties.CustomProperties.zsc;
      obj.fzs = obj.data.Properties.CustomProperties.fzs;
      obj.dcal = obj.data.Properties.CustomProperties.dcal;
      obj.config = obj.data.Properties.CustomProperties.config;
      obj.housek = obj.data.Properties.CustomProperties.housek;
    end
    
    function ReadRawDI(obj, days2run, force_import, write)
      % Get wavelengths from calibration file
      obj.ReadHBBCalFiles()
      % Set default parameters
      if isempty(obj.path.di)
        fprintf('WARNING: DI Path is same as raw.\n');
        obj.path.di = obj.path.raw;
      end
      if isempty(obj.di_cfg.logger)
        fprintf('WARNING: DI Logger set to "internal_logger".\n');
        obj.di_cfg.logger = 'internal_logger';
      end
      if isempty(obj.di_cfg.postfix)
        fprintf('WARNING: DI Postfix isempty \n');
%         fprintf('WARNING: DI Postfix set to "_DI" \n'); DEPRECATED with Inlinino
%         obj.di_cfg.postfix = '_DI'; DEPRECATED with Inlinino
      end
      if isempty(obj.di_cfg.prefix); obj.di_cfg.prefix = ['LISST200X' obj.sn '_']; end
      switch obj.di_cfg.logger
        case 'internal_logger'
          obj.raw.diw = iRead(@importLISST200x_raw, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'internal_logger', force_import, ~write, true, true, ...
                         obj.di_cfg.postfix, Inf, obj.calfile_plaque, obj.calfile_temp);
        otherwise
          error('LISST200X: Unknown logger.');
      end
    end
    
    function Calibrate(obj, days2run, compute_dissolved, SWT, di_method)
      SWT_constants = struct('SWITCH_FILTERED', SWT.SWITCH_FILTERED, 'SWITCH_TOTAL', SWT.SWITCH_TOTAL);
      param = struct('zsc', obj.zsc, 'fzs', obj.fzs, 'dcal', obj.dcal, 'config', obj.config,...
                     'housek', obj.housek, 'inversion', obj.inversion);
      if compute_dissolved
        obj.prod.p = ProcessLISST200X(param, obj.qc.tsw, obj.qc.fsw, [], ...
          SWT, SWT_constants, di_method, days2run);
      else
        obj.prod.p = ProcessLISST200X(param, obj.qc.tsw, obj.qc.fsw, obj.bin.diw, ...
          SWT, SWT_constants, di_method, days2run);
      end
    end
  end
end
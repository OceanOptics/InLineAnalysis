classdef HBB < Instrument
  % HBB Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    calfile_plaque = '';
    calfile_temp = '';
    lambda = [];
    theta = [];
%     muFactors = [];
  end
  
  methods
    function obj = HBB(cfg)
      % HBB Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Change default processing method
      obj.bin_method = 'SB_IN_PRCTL';
      
      % Post initialization
      if isempty(obj.view.varname); obj.view.varname = 'beta'; end
      if isfield(cfg, 'calfile_plaque'); obj.calfile_plaque = cfg.calfile_plaque;
      else; error('Missing field calfile_plaque.'); end
      if isfield(cfg, 'calfile_temp'); obj.calfile_temp = cfg.calfile_temp;
      else; error('Missing field calfile_temp.'); end
      if isfield(cfg, 'theta'); obj.theta = cfg.theta;
      else; error('Missing field theta.'); end
%       if isfield(cfg, 'muFactors'); obj.muFactors = cfg.muFactors;
%       else; error('Missing field muFactors.'); end

      if isempty(obj.logger)
        fprintf('WARNING: Logger set to InlininoHBB.\n');
        obj.logger = 'InlininoHBB';
      end
    end
    
    function ReadHBBCalFiles(obj)
      % Load HBB lambda from Calibration files
      load(obj.calfile_temp, 'cal_temp');
      obj.lambda = cal_temp.wl;
    end
    
    function ReadRaw(obj, days2run, force_import, write)
      % Get wavelengths from calibration file
      obj.ReadHBBCalFiles()
      % create wk directory if doesn't exist
      if ~isfolder(obj.path.wk); mkdir(obj.path.wk); end
      % Read raw data
      switch obj.logger
        case 'InlininoHBB'
          obj.data = iRead(@importInlininoHBB, obj.path.raw, obj.path.wk, ['HyperBB' obj.sn '_'],...
                         days2run, 'Inlinino', force_import, ~write, true, true, '', Inf, ...
                         obj.calfile_plaque, obj.calfile_temp);
        otherwise
          error('HBB: Unknown logger.');
      end
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
        fprintf('WARNING: DI Logger set to InlininoHBB.\n');
        obj.di_cfg.logger = 'InlininoHBB';
      end
      if isempty(obj.di_cfg.postfix)
        fprintf('WARNING: DI Postfix isempty \n');
%         fprintf('WARNING: DI Postfix set to "_DI" \n'); DEPRECATED with Inlinino
%         obj.di_cfg.postfix = '_DI'; DEPRECATED with Inlinino
      end
      if isempty(obj.di_cfg.prefix); obj.di_cfg.prefix = ['HyperBB' obj.sn '_']; end
      switch obj.di_cfg.logger
        case 'InlininoHBB'
          obj.raw.diw = iRead(@importInlininoHBB, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'Inlinino', force_import, ~write, true, true, '', Inf, ...
                         obj.calfile_plaque, obj.calfile_temp);
        otherwise
          error('HBB: Unknown logger.');
      end
    end
    
    function Calibrate(obj, compute_dissolved, TSG, di_method)
      param = struct('lambda', obj.lambda, 'theta', obj.theta);
%       param = struct('lambda', obj.lambda, 'theta', obj.theta, 'muFactors', obj.muFactors);
      % linear interpolation only, CDOM interpolation is not yet available
      if compute_dissolved
        [obj.prod.p, obj.prod.g] = processHBB(param, obj.qc.tsw, obj.qc.fsw, obj.bin.diw, TSG.qc.tsw, di_method);
      else
        [obj.prod.p] = processHBB(param, obj.qc.tsw, obj.qc.fsw);
      end
    end
  end
end
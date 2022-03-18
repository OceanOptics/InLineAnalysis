classdef AC9 < Instrument
  %AC9 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    device_file = '';
    lambda_ref = [];
    lambda_c = [];
    lambda_a = [];
    modelG50 = [];
    modelmphi = [];
  end
  
  methods
    function obj = AC9(cfg)
      %AC9 Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Change default processing method
      obj.bin_method = 'SB_IN_PRCTL';
      
      % Post initialization
%       if isfield(cfg, 'device_file'); obj.device_file = cfg.device_file;
%       else; error('Missing field device_file.'); end
%       % Set lambda reference to c by default
%       if isfield(cfg, 'lambda_reference'); obj.lambda_ref = cfg.lambda_reference; end
      % DEPRECATED AS will be read in device file now
%       if isfield(cfg, 'lambda_a'); obj.lambda_a = cfg.lambda_a;
%       else; error('Missing field lambda_a.'); end
%       if isfield(cfg, 'lambda_c'); obj.lambda_c = cfg.lambda_c;
%       else; error('Missing field lambda_c.'); end

      if isempty(obj.logger) 
        fprintf('WARNING: Logger set to WetView.\n');
        obj.logger = 'WetView';
      end
    end    
    
%     function ReadDeviceFile(obj)
%         % Read Device file to set wavelength
%       [obj.lambda_c, obj.lambda_a] = importAC9DeviceFile(obj.device_file);
%       if isempty(obj.lambda_ref); obj.lambda_ref = obj.lambda_c; end
%     end
    function ReadDeviceFile(obj)
        % Read Device file to set wavelength
      obj.lambda_ref = [412 440 488 510 532 555 650 676 715];
      obj.lambda_a = [412 440 488 510 532 555 650 676 715];
      obj.lambda_c = [412 440 488 510 532 555 650 676 715];
%       if isempty(obj.lambda_ref); obj.lambda_ref = obj.lambda_c; end
    end
    
    function obj = load_HTJSetal2021_model(obj)
      % Load model from Haëntjens et al. 2021v22 to estimate cross-sectional
      % area (G50) and slope of phytoplankton size distribution (in abundance) (mphi):
      load('HTJS20_LinearRegression_5features.mat', 'model_G50');
      load('HTJS20_LinearRegression_5P-mphi.mat', 'model_mphi');
      obj.modelG50 = model_G50;
      obj.modelmphi = model_mphi;
    end
    
    function ReadRaw(obj, days2run, force_import, write)
      % Get wavelengths from device file
      obj.ReadDeviceFile()
      % Read raw data
      switch obj.logger
        case 'WetView'
          obj.data = iRead(@importAC9wetview, obj.path.raw, obj.path.wk, ['ac9_' obj.sn '_'],...
                         days2run, 'WetView', force_import, ~write, true);
        otherwise
          error('AC9: Unknown logger.');
      end
    end
    
    function ReadRawDI(obj, days2run, force_import, write)
      % Set default parameters
      if isempty(obj.path.di)
        fprintf('WARNING: DI Path is same as raw.\n');
        obj.path.di = obj.path.raw;
      end
      if isempty(obj.di_cfg.logger) 
        fprintf('WARNING: DI Logger set to WetView.\n');
        obj.di_cfg.logger = 'WetView';
      end
      if isempty(obj.di_cfg.postfix) 
        fprintf('WARNING: DI Postfix set to "_DI" \n');
        obj.di_cfg.postfix = '_DI';
      end
      switch obj.di_cfg.logger
        case 'WetView'
          if isempty(obj.di_cfg.prefix); obj.di_cfg.prefix = ['ac9' obj.sn '_']; end
          obj.raw.diw = iRead(@importAC9wetview, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'WetView', force_import, ~write, true, false, obj.di_cfg.postfix);
        otherwise
          error('AC9: Unknown logger.');
      end
    end
    
    function Calibrate(obj, compute_dissolved, interpolation_method, CDOM, SWT, di_method, compute_ad_aphi)
      lambda = struct('ref', obj.lambda_ref, 'a', obj.lambda_ref, 'c', obj.lambda_ref);
      SWT_constants = struct('SWITCH_FILTERED', SWT.SWITCH_FILTERED, 'SWITCH_TOTAL', SWT.SWITCH_TOTAL);
      % Load model from Haëntjens et al. 2021v22 to estimate cross-sectional
      % area (G50) and slope of phytoplankton size distribution (in abundance) (mphi):
      obj.load_HTJSetal2021_model();
      switch interpolation_method
        case 'linear'
          if compute_dissolved
            [obj.prod.p, obj.prod.g, obj.prod.QCfailed] = processACS(lambda, ...
              obj.qc.tsw, obj.qc.fsw, obj.modelG50, obj.modelmphi, obj.bin.diw, ...
              [], SWT, SWT_constants, interpolation_method, di_method, ...
              compute_ad_aphi);
          else
            [obj.prod.p, ~, obj.prod.QCfailed] = processACS(lambda, ...
              obj.qc.tsw, obj.qc.fsw, obj.modelG50, obj.modelmphi, [], ...
              [], SWT, SWT_constants, interpolation_method, [],...
              compute_ad_aphi);
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
            [obj.prod.p, obj.prod.g, obj.prod.QCfailed] = processACS(lambda, ...
              obj.qc.tsw, obj.qc.fsw, obj.modelG50, obj.modelmphi, obj.bin.diw, ...
              cdom, SWT, SWT_constants, interpolation_method, di_method, ...
              compute_ad_aphi);
          else
            [obj.prod.p, ~, obj.prod.QCfailed] = processACS(lambda, ...
              obj.qc.tsw, obj.qc.fsw, obj.modelG50, obj.modelmphi, [], ...
              cdom, SWT, SWT_constants, interpolation_method, [],...
              compute_ad_aphi);
          end
        otherwise
          error('Method not supported.');
      end
    end
  end
end
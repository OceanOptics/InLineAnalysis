classdef ACS < Instrument
  %ACS Summary of this class goes here
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
    function obj = ACS(cfg)
      %ACS Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Change default processing method
      obj.bin_method = 'SB_IN_PRCTL';
      
      % Post initialization
      if isfield(cfg, 'device_file'); obj.device_file = cfg.device_file;
      else; error('Missing field device_file.'); end
      % Set lambda reference to c by default
      if isfield(cfg, 'lambda_reference'); obj.lambda_ref = cfg.lambda_reference; end
      % DEPRECATED AS will be red in device file now
%       if isfield(cfg, 'lambda_a'); obj.lambda_a = cfg.lambda_a;
%       else; error('Missing field lambda_a.'); end
%       if isfield(cfg, 'lambda_c'); obj.lambda_c = cfg.lambda_c;
%       else; error('Missing field lambda_c.'); end

      if isempty(obj.logger) 
        fprintf('WARNING: Logger set to Compass_2.1rc_scheduled.\n');
        obj.logger = 'Compass_2.1rc_scheduled';
      end
    end    
    
    function ReadDeviceFile(obj)
      % Read Device file to set wavelength
      [obj.lambda_c, obj.lambda_a] = importACSDeviceFile(obj.device_file);
      if isempty(obj.lambda_ref); obj.lambda_ref = obj.lambda_a; end
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
      % create wk directory if doesn't exist
      if ~isfolder(obj.path.wk); mkdir(obj.path.wk); end
      % Read raw data
      switch obj.logger
        case 'InlininoACScsv'
          obj.data = iRead(@importInlininoACScsv, obj.path.raw, obj.path.wk, ['acs' obj.sn '_'],...
                         days2run, 'InlininoACScsv', force_import, ~write, true);
        case 'Compass_2.1rc_scheduled'
          warning('DEPRECATED: as wavelength are shifted when saved to .dat by the scheduler of compass.');
          obj.data = iRead(@importACS, obj.path.raw, obj.path.wk, ['acs' obj.sn '_'],...
                         days2run, 'Compass_2.1rc_scheduled', force_import, ~write, true);
        case 'Compass_2.1rc_scheduled_bin'
          obj.data = iRead(@importACSBin, obj.path.raw, obj.path.wk, ['acs' obj.sn '_'],...
                         days2run, 'Compass_2.1rc_scheduled_bin', force_import, ~write, true, true, '', Inf, obj.device_file);
        case 'Compass_2.1rc'
          obj.data = iRead(@importACS, obj.path.raw, obj.path.wk, ['acs' obj.sn '_'],...
                         days2run, 'Compass_2.1rc', force_import, ~write, true);
        case 'WetView'
          obj.data = iRead(@importACSwetview, obj.path.raw, obj.path.wk, ['acs' obj.sn '_'],...
                         days2run, 'WetView', force_import, ~write, true);
        otherwise
          error('ACS: Unknown logger.');
      end
    end
    
    function ReadRawDI(obj, days2run, force_import, write)
      % Get wavelengths from device file
      obj.ReadDeviceFile()
      % Set default parameters
      if isempty(obj.path.di)
        fprintf('WARNING: DI Path is same as raw.\n');
        obj.path.di = obj.path.raw;
      end
      if isempty(obj.di_cfg.logger)
        fprintf('WARNING: DI Logger set to Compass_2.1rc.\n');
        obj.di_cfg.logger = 'Compass_2.1rc';
      end
      if isempty(obj.di_cfg.postfix)
        fprintf('WARNING: DI Postfix isempty \n');
%         fprintf('WARNING: DI Postfix set to "_DI" \n'); DEPRECATED with Inlinino
%         obj.di_cfg.postfix = '_DI'; DEPRECATED with Inlinino
      end
      if isempty(obj.di_cfg.prefix); obj.di_cfg.prefix = ['acs' obj.sn '_']; end
      switch obj.di_cfg.logger
        case 'InlininoACScsv'
          obj.raw.diw = iRead(@importInlininoACScsv, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'InlininoACScsv', force_import, ~write, true, false, obj.di_cfg.postfix);
        case 'Compass_2.1rc_scheduled'
          obj.raw.diw = iRead(@importACS, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'Compass_2.1rc_scheduled', force_import, ~write, true, false, obj.di_cfg.postfix);
        case 'Compass_2.1rc_scheduled_bin'
          obj.raw.diw = iRead(@importACSBin, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'Compass_2.1rc_scheduled_bin', force_import, ~write, true, false, obj.di_cfg.postfix);
        case 'Compass_2.1rc'
          obj.raw.diw = iRead(@importACS, obj.path.di, obj.path.wk, obj.di_cfg.prefix,...
                         days2run, 'Compass_2.1rc', force_import, ~write, true, false, obj.di_cfg.postfix);
        otherwise
          error('ACS: Unknown logger.');
      end
    end
    
    function Calibrate(obj, compute_dissolved, interpolation_method, CDOM, SWT)
      lambda = struct('ref', obj.lambda_ref, 'a', obj.lambda_a, 'c', obj.lambda_c);
      SWT_constants = struct('SWITCH_FILTERED', SWT.SWITCH_FILTERED, 'SWITCH_TOTAL', SWT.SWITCH_TOTAL);
      % Load model from Haëntjens et al. 2021v22 to estimate cross-sectional
      % area (G50) and slope of phytoplankton size distribution (in abundance) (mphi):
      obj.load_HTJSetal2021_model();
      switch interpolation_method
        case 'linear'
          if compute_dissolved
            [obj.prod.p, obj.prod.g, obj.prod.QCfailed] = processACS(lambda, obj.qc.tsw, obj.qc.fsw, ...
              obj.modelG50, obj.modelmphi, obj.bin.diw, [], SWT.qc.tsw, SWT_constants, interpolation_method);
          else
            [obj.prod.p, obj.prod.g, obj.prod.QCfailed] = processACS(lambda, obj.qc.tsw, obj.qc.fsw, obj.modelG50, ...
              obj.modelmphi, [], [], SWT.qc.tsw, SWT_constants, interpolation_method);
          end
        case 'CDOM'
          if compute_dissolved
            [obj.prod.p, obj.prod.g, obj.prod.QCfailed] = processACS(lambda, obj.qc.tsw, obj.qc.fsw, ...
              obj.modelG50, obj.modelmphi, obj.bin.diw, CDOM.qc.tsw, SWT.qc.tsw, SWT_constants, interpolation_method);
          else
            [obj.prod.p, obj.prod.g, obj.prod.QCfailed] = processACS(lambda, obj.qc.tsw, obj.qc.fsw, obj.modelG50, obj.modelmphi, ...
              [], CDOM.qc.tsw, SWT.qc.tsw, SWT_constants, interpolation_method);
          end
        otherwise
          error('Method not supported.');
      end
    end
  end
end
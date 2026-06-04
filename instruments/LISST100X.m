classdef LISST100X < Instrument
  % LISST100X Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    type = '';
    % Calibration parameters
    zsc = NaN(1,32);
    dcal = NaN(1,32);
    vcc = NaN;
    % Type of inversion
    inversion = 'spherical';
    ds = NaN(1,33); % 1.25*1.18.^(0:1:32); Type B & Spherical
                    % 1*1.18.^(0:1:32);    Type B & Non-Spherical
    diameters = NaN(1,32);
    % VSF Angles
    theta = NaN(1,32);
  end
  
  methods
    function obj = LISST100X(cfg)
      % LISST100X Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Change default processing method
      obj.bin_method = 'SB_IN_PRCTL';
      
      % Load LISST calibration data
      if isfield(cfg, 'zsc')
        if isnumeric(cfg.zsc); obj.zsc = cfg.zsc; else; obj.zsc = importLISST100XDeviceFile(cfg.zsc); end
      else
        error('Missing field zsc.')
      end
      if isfield(cfg, 'dcal')
        if isnumeric(cfg.dcal); obj.dcal = cfg.dcal; else; obj.dcal = importLISST100XDeviceFile(cfg.dcal); end
      else
        error('Missing field dcal.')
      end
      if isfield(cfg, 'vcc')
        if isnumeric(cfg.vcc)
          obj.vcc = cfg.vcc;
        else
          instrument_data = importLISST100XDeviceFile(cfg.vcc);
          obj.vcc = instrument_data{4};
          if ~strcmp(num2str(instrument_data{1}), cfg.sn)
            error("Serial number in calibration file '%s' is different than the one entered in .cfg file", cfg.vcc)
          end
          if isfield(cfg, 'type')
            if ~isempty(cfg.type)
              if ~strcmpi(cfg.type, cell2mat(instrument_data{2}))
                error("LISST100X type in calibration file: '%s', is different than the one entered in .cfg file: '%s'", ...
                  lower(cfg.type), lower(cell2mat(instrument_data{2})))
              end
              obj.type = cfg.type;
            else
              obj.type = instrument_data{2};
            end
          end
        end
        if ~any(strcmpi(obj.type, {'b','c'}))
          error("LISST100X type '%s' not supported: must be either 'b' or 'c'", lower(obj.type))
        end
      else
        warning("LISST100X import: Missing field 'vcc', set vcc to 13000 (default).")
        obj.vcc = 13000;
        if isfield(cfg, 'type')
          warning("LISST100X import: Missing field 'type', set 'type' to 'b' (default).")
          obj.type = 'b';
        end
      end
      if isfield(cfg, 'inversion')
        if strcmpi(obj.inversion, 'spherical') || strcmpi(obj.inversion, 'non-spherical')
          obj.inversion = cfg.inversion;
        else
          error('Unknown inversion type.');
        end
      else
        fprintf('WARNING: Set inversion to spherical (default).'); 
      end
      if isfield(cfg, 'ds')
        obj.ds = cfg.ds;
      elseif strcmpi(obj.inversion, 'spherical')
        if strcmpi(obj.type, 'b')
          X = 1.25;
          fprintf("'ds' set to default (spherical inversion, LISST100X Type B).\n");
        elseif strcmpi(obj.type, 'c')
          X = 2.5;
          fprintf("'ds' set to default (spherical inversion, LISST100X Type C).\n");
        end
      elseif strcmpi(obj.inversion, 'non-spherical')
        if strcmpi(obj.type, 'b')
          X = 1;
          fprintf("'ds' set to default (non-spherical inversion, LISST100X Type B).\n");
        elseif strcmpi(obj.type, 'c')
          X = 1.9;
          fprintf("'ds' set to default (non-spherical inversion, LISST100X Type C).\n");
        end
      end
      if ~isfield(cfg, 'ds')
        obj.ds = X*200^(1/32).^(0:1:32);
      end

      if isfield(cfg, 'theta')
        obj.theta = cfg.theta;
      elseif strcmpi(obj.type, 'b')
        start_angle = 0.1;
        fprintf("'theta' (VSF Angle) set to default (LISST100X Type B).\n");
      elseif strcmpi(obj.type, 'c')
        start_angle = 0.05;
        fprintf("'theta' (VSF Angle) set to default (LISST100X Type C).\n");
      end
      if ~isfield(cfg, 'theta')
        angles_foo = (logspace(0,log10(200),33)*start_angle)';
        angles_O2(:,1) = angles_foo(1:32); % The lower limits are the first 32 (of 33)
        angles_O2(:,2) = angles_foo(2:33); % The upper limits are the last 32 (of 33)
        obj.theta = sqrt(angles_O2(:,1).*angles_O2(:,2)) ./ 1.33; % Midpoints / water refractive index to get angles in water
      end
      % Compute diameters from initialized parameters
      %   avoid special method for writting and reading data
      obj.ComputeDiameters();
      
      if isempty(obj.logger); obj.logger = 'TeraTerm'; end
      if isempty(obj.view.varname); obj.view.varname = 'beta'; end
    end    
    
    function ReadRaw(obj, days2run, force_import, write)
      switch obj.logger
        case 'InlininoLISSTcsv'
          obj.data = iRead(@importInlininoLISSTcsv, obj.path.raw, obj.path.wk, obj.prefix,...
                         days2run, 'Inlinino', force_import, ~write, true);
        case 'TeraTerm'
          obj.data = iRead(@importLISSTTeraTerm, obj.path.raw, obj.path.wk, obj.prefix,...
                         days2run, 'TeraTerm', force_import, ~write, true);
        otherwise
          error('LISST: Unknown logger.');
      end
    end

    function ReadRawDI(obj, days2run, force_import, write)
      if isempty(obj.path.di)
        fprintf('WARNING: DI Path is same as raw.\n');
        obj.path.di = obj.path.raw;
      end
      if isempty(obj.di_cfg.postfix) 
        fprintf('WARNING: DI Postfix set to "_DI" \n');
        obj.di_cfg.postfix = '_DI';
      end
      switch obj.logger
        case 'InlininoLISSTcsv'
          obj.raw.diw = iRead(@importInlininoLISSTcsv, obj.path.di, obj.path.wk, obj.prefix,...
                         days2run, 'Inlinino', force_import, ~write, true, false, obj.di_cfg.postfix);
        case 'TeraTerm'
          obj.raw.diw = iRead(@importLISSTTeraTerm, obj.path.di, obj.path.wk, obj.prefix,...
                         days2run, 'TeraTerm', force_import, ~write, true, false, obj.di_cfg.postfix);
        otherwise
          error('LISST100X: Unknown logger.');
      end
    end

    function ComputeDiameters(obj)
      obj.diameters = sqrt(obj.ds(1:end-1) .* obj.ds(2:end));
    end
    
    function Calibrate(obj, days2run, compute_dissolved, SWT, di_method)
      SWT_constants = struct('SWITCH_FILTERED', SWT.SWITCH_FILTERED, 'SWITCH_TOTAL', SWT.SWITCH_TOTAL);
      param = struct('zsc', obj.zsc, 'dcal', obj.dcal, 'vcc', obj.vcc,...
                     'inversion', obj.inversion, 'ds', obj.ds,...
                     'diameters', obj.diameters, 'theta', obj.theta);
      if compute_dissolved
        obj.prod.p = processLISST100X(param, obj.qc.tsw, obj.qc.fsw, [], ...
          SWT, SWT_constants, di_method, days2run);
      else
        obj.prod.p = processLISST100X(param, obj.qc.tsw, obj.qc.fsw, obj.bin.diw, ...
          SWT, SWT_constants, di_method, days2run);
      end
    end
    
    % Old deprecated method overload (can now use standard methods)
%     function Write(obj, filename_prefix, days2write)
%       % Overload instrument write class
%       % Call superclass Write class
% %       Write@Instrument(obj, filename_prefix, days2write) % Not needed as
% %       modify the full fprocess
%       % For each product type (particulate, dissoved...)
%       for f = fieldnames(obj.prod); f = f{1};
%         filename = [filename_prefix '_' f '_prod.mat'];
%         sel = min(days2write) <= obj.prod.(f).dt & obj.prod.(f).dt < max(days2write) + 1;
%         data = obj.prod.(f)(sel,:);
%         diameters = obj.diameters;
%         if ~isdir(obj.path.prod); mkdir(obj.path.prod); end
%         save([obj.path.prod filename], 'data', 'diameters');
%       end
%     end
% 
%     function LoadProducts(obj, filename_prefix, days2read)
%       % Overload instrument LoadProducts class
%       % For each product type (particulate, dissoved...)
%       l = dir([obj.path.prod filename_prefix '_*_prod.mat']);
%       for f = {l.name}'; f = f{1};
%         load([obj.path.prod f]); % data variable is create
%         if exist('diameters', 'var'); obj.diameters = diameters; end
%         sel = min(days2read) <= data.dt & data.dt < max(days2read) + 1;
%         fn = strsplit(f, '_'); fn = fn{end-1};
%         if isfield(obj.prod, fn)
%           obj.prod.(fn)(end+1:end+sum(sel),:) = data(sel,:);
%         else
%           obj.prod.(fn) = data(sel,:);
%         end
%       end
%     end

  end
end
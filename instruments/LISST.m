classdef LISST < Instrument
  %LISST Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    type = 'B';
    % Calibration parameters
    zsc = NaN(1,32);
    dcal = NaN(1,32);
    vcc = NaN;
    % Type of inversion
    inversion = 'spherical';
    ds = NaN(1,33); % 1.25*1.18.^(0:1:32); Type B & Spherical
                    % 1*1.18.^(0:1:32);    Type B & Non-Spherical
    diameters = NaN(1,32);
  end
  
  methods
    function obj = LISST(cfg)
      %LISST Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Post initialization
      if isfield(cfg, 'zsc'); obj.zsc = cfg.zsc;
      else; error('Missing field zsc.'); end
      if isfield(cfg, 'dcal'); obj.dcal = cfg.dcal;
      else; error('Missing field dcal.'); end
      if isfield(cfg, 'type'); obj.type = cfg.type;
      else; fprintf('WARNING: Set type to B (default).'); end
      if isfield(cfg, 'vcc'); obj.vcc = cfg.vcc;
      else; fprintf('WARNING: Set vcc to 13000 (default).'); end
      if isfield(cfg, 'inversion'); obj.inversion = cfg.inversion;
      else; fprintf('WARNING: Set inversion to spherical (default).'); end
      if isfield(cfg, 'ds'); obj.ds = cfg.ds;
      else; fprintf('WARNING: Set ds to default (spherical inversion, LISST Type B).'); end
      
%       obj.diameters = NaN(1,33);
      
      if isempty(obj.logger); obj.logger = 'TeraTerm'; end
    end    
    
    function ReadRaw(obj, days2run, force_import, write)
      switch obj.logger
        case 'TeraTerm'
          obj.data = iRead(@importLISSTTeraTerm, obj.path.raw, obj.path.wk, ['LISST' obj.sn '_'],...
                         days2run, 'TeraTerm', force_import, ~write, true);
        otherwise
          error('LISST: Unknown logger.');
      end
    end
    
    function ReadDI(obj, days2run, force_import, write)
      if isempty(obj.path.di)
        fprintf('WARNING: DI Path is same as raw.\n');
        obj.path.di = obj.path.raw;
      end
      if isempty(obj.di_cfg.postfix) 
        fprintf('WARNING: DI Postfix set to "_DI" \n');
        obj.di_cfg.postfix = '_DI';
      end
      switch obj.logger
        case 'TeraTerm'
          obj.raw.diw = iRead(@importLISSTTeraTerm, obj.path.di, obj.path.wk, ['LISST' obj.sn '_'],...
                         days2run, 'TeraTerm', force_import, ~write, true, false, obj.di_cfg.postfix);
        otherwise
          error('LISST: Unknown logger.');
      end
    end

    function Calibrate(obj)
      param = struct('zsc', obj.zsc, 'dcal', obj.dcal, 'vcc', obj.vcc, 'non_spherical', NaN, 'ds', obj.ds);
      switch obj.inversion
        case 'spherical'
          param.non_spherical = 0;
        case 'non-spherical'
          param.non_spherical = 1;
      end
      [obj.prod.p, obj.diameters] = processLISST(param, obj.qc.tsw, obj.qc.fsw);
    end
    
    function Write(obj, filename_prefix, days2write)
      % Overload instrument write class
      % Call superclass Write class
%       Write@Instrument(obj, filename_prefix, days2write) % Not needed as
%       modify the full fprocess
      % For each product type (particulate, dissoved...)
      for f = fieldnames(obj.prod); f = f{1};
        filename = [filename_prefix '_' f '_prod.mat'];
        sel = min(days2write) <= obj.prod.(f).dt & obj.prod.(f).dt < max(days2write) + 1;
        data = obj.prod.(f)(sel,:);
        diameters = obj.diameters;
        if ~isdir(obj.path.prod); mkdir(obj.path.prod); end
        save([obj.path.prod filename], 'data', 'diameters');
      end
    end

    function LoadProducts(obj, filename_prefix, days2read)
      % Overload instrument LoadProducts class
      % For each product type (particulate, dissoved...)
      l = dir([obj.path.prod filename_prefix '_*_prod.mat']);
      for f = {l.name}'; f = f{1};
        load([obj.path.prod f]); % data variable is create
        if exist('diameters', 'var'); obj.diameters = diameters; end
        sel = min(days2read) <= data.dt & data.dt < max(days2read) + 1;
        fn = strsplit(f, '_'); fn = fn{end-1};
        if isfield(obj.prod, fn)
          obj.prod.(fn)(end+1:end+sum(sel),:) = data(sel,:);
        else
          obj.prod.(fn) = data(sel,:);
        end
      end
    end

  end
end
classdef LISST < Instrument
  %LISST Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    type = 'B';
    zsc = NaN(1,32);
    dcal = NaN(1,32);
    vcc = 13000;
    inversion = 'spherical';
    ds = NaN(1,33); % 1.25*1.18.^(0:1:32); Type B & Spherical
                    % 1*1.18.^(0:1:32);    Type B & Non-Spherical
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
      
      if isempty(obj.logger); obj.logger = 'TeraTerm'; end
    end    
    
    function ReadRaw(obj, days2run, force_import, write)
      switch obj.logger
        case 'TeraTerm'
          obj.data = iRead(@importLISSTTeraTerm, obj.path.raw, obj.path.wk, 'LISST1183_',...
                         days2run, 'TeraTerm', force_import, ~write, true);
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
  end
end
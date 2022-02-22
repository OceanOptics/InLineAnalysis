classdef atlasTSG < Instrument
  %TSG Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    TSseparated = '';
    coef = struct('t', struct('slope','','interc',''), ...
      'c', struct('slope','','interc',''));
%       'c', struct('p1','','p2','','p3',''));
  end
  
  methods
    function obj = atlasTSG(cfg)
      %TSG Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Post initialization
      if isfield(cfg, 'TSseparated'); obj.TSseparated = cfg.TSseparated;
      else; error('Missing field "TSseparated".'); end
      if isfield(cfg, 'coef'); obj.coef = cfg.coef;
      else; error('Missing field "coef": Strucutre containing calibration coefficients.'); end
      
      % Change default Split method
      obj.split.mode = 'None';
    end
    
    function ReadRaw(obj, days2run, force_import, write)
      switch obj.logger
        case {'Inlinino', 'Inlinino_atlasTSG'}
          if obj.TSseparated
            obj.data = iRead(@importInlinino_atlasTSG, obj.path.raw, obj.path.wk, 'AST01_',...
                           days2run, 'Inlinino_atlasTSG', force_import, ~write, true);
            foo = iRead(@importInlinino_atlasTSG, obj.path.raw, obj.path.wk, 'ASC01_',...
                           days2run, 'Inlinino_atlasTSG', force_import, ~write, true);
            % interpolate C onto T table and merge
            obj.data.C = interp1(foo.dt, foo.C, obj.data.dt, 'linear', 'extrap'); % extrap needed for first minute of data
          else
            obj.data = iRead(@importInlinino_atlasTSG, obj.path.raw, obj.path.wk, 'AS_',...
                           days2run, 'Inlinino_atlasTSG', force_import, ~write, true);
          end
        otherwise
          error('atlasTSG: Unknown logger.');
      end
      obj.data = renamevars(obj.data, 'T' , 't');
      obj.data = renamevars(obj.data, 'C' , 'c');
    end
    
    function Calibrate(obj)
      % apply calibration coefficients
      obj.prod.a = obj.qc.tsw;
%       obj.prod.a.c = obj.coef.c.p1 .* obj.prod.a.c.^2 + obj.coef.c.p2 .* obj.prod.a.c + obj.coef.c.p3;
      obj.prod.a.c = obj.prod.a.c / obj.coef.c.slope + obj.coef.c.interc;
      obj.prod.a.t = obj.prod.a.t / obj.coef.t.slope + obj.coef.t.interc;
      % compute salinity from conductivity and temperature using average pressure at 5
       % C * 10 to get mS/cm to input GSW | pressure
      obj.prod.a.s = gsw_SP_from_C(obj.prod.a.c * 10, obj.prod.a.t, 36*11.1094 - 10.1325);
      obj.prod.a.s = gsw_SA_from_SP(obj.prod.a.s, 36*11.1094 - 10.1325, -14.8, 49.1); % absolute salinity from practical salinity
    end
  end
end


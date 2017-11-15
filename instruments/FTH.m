classdef FTH < Instrument
  %FTH Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Constant=true, Hidden=true)
    SAMPLING_FREQUENCY = 1; % Hz
  end
  
  properties (Hidden=true)
    SWITCH_FILTERED = 1;
    SWITCH_TOTAL = 0;
  end
  
  methods
    function obj = FTH(cfg)
      %FTH Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Post initialization
      if isempty(obj.logger); obj.logger = 'FlowControl'; end
      
      switch obj.logger
        case 'FlowControl_old'
          obj.SWITCH_FILTERED = 0;
          obj.SWITCH_TOTAL = 1;
      end
      % Change default Split method
      obj.split.mode = 'None';
    end
    
    function ReadRaw(obj, days2run, force_import, write)
      switch obj.logger
        case {'FlowControl', 'FlowControl_old'}
          obj.data = iRead(@importFlowControl, obj.path.raw, obj.path.wk, 'Flow_',...
                         days2run, 'FlowControl', force_import, ~write, true);
        otherwise
          error('FTH: Unknown logger.');
      end
    end
    
    function ApplyUserInput(obj, user_selection, mode)
      % Correct part of switch data
      % Note: spd data is lost when corrected
      for i=1:size(user_selection,1)
        % Remove existing data from fth
        obj.data(user_selection(i,1) <= obj.data.dt &...
          obj.data.dt <= user_selection(i,2),:) = [];
        % Add user selection in fth
        dt = [user_selection(i,1):1/obj.SAMPLING_FREQUENCY/3600/24:user_selection(i,2)]';
        switch mode
          case 'total'
            obj.data = [obj.data; table(dt, ones(size(dt)) * obj.SWITCH_TOTAL, NaN(size(dt)), 'VariableNames', {'dt', 'swt', 'spd'})];
          case 'filtered'
            obj.data = [obj.data; table(dt, ones(size(dt)) * obj.SWITCH_FILTERED, NaN(size(dt)), 'VariableNames', {'dt', 'swt', 'spd'})];
          otherwise
            error('Unknown mode.');
        end
      end
      obj.data = sortrows(obj.data);
    end
  end
end
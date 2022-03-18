classdef FTH < Instrument
  %FTH Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Constant=true, Hidden=true)
    SAMPLING_FREQUENCY = 1; % Hz
  end
  
  properties (Hidden=true)
    SWITCH_FILTERED = 1;
    SWITCH_TOTAL = 0;
    spd_variable = [];
  end
  
  methods
    function obj = FTH(cfg)
      %FTH Construct an instance of this class
      
      % Object Initilization
      obj = obj@Instrument(cfg);
      
      % Post initialization
      if isempty(cfg.logger)
        obj.logger = 'FlowControl';
      else
        obj.logger = cfg.logger;
      end
      if isempty(cfg.view.spd_variable)
        obj.view.spd_variable = 'spd';
      else
        obj.view.spd_variable = cfg.view.spd_variable;
      end
      
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
      fprintf('User input ');
      for i=progress(1:size(user_selection,1))
        % Add user selection in fth
        % round user input to the second and create continuous time vector
        dt_st = datenum(floor(datevec(user_selection(i,1))));
        dt_end = datenum(floor(datevec(user_selection(i,2))));%, 'ConvertFrom','datenum', 'Format', 'yyyy/MM/dd hh:mm:ss') ...
        dt = (dt_st:1/obj.SAMPLING_FREQUENCY/3600/24:dt_end)';
        
        % delete data with duplicats timestamp
        if ~ isempty(obj.data.dt)
          % round data time vector to the second
          obj.data.dt = datenum(floor(datevec(obj.data.dt)));
          % delete duplicates
          [~, L, ~] = unique(obj.data.dt,'first');
          indexToDump = not(ismember(1:numel(obj.data.dt),L));
          obj.data(indexToDump, :) = [];
        else
          error('Raw data not loaded')
        end
        % interpolate flow rate data over continuous time vector
        spd_var = obj.data.Properties.VariableNames(contains(obj.data.Properties.VariableNames, 'spd'));
%         interp_spd = interp1(obj.data.dt, obj.data.(obj.view.spd_variable), dt, 'linear');
        interp_spd = NaN(size(dt,1),2);
        for j = 1:size(spd_var, 2)
          interp_spd(:,j) = interp1(obj.data.dt, obj.data.(spd_var{j}), dt, 'linear');
        end
        % delete interpolated data to keep only 'true' data
        ism = ~ismember(dt, obj.data.dt);
        interp_spd(ism) = NaN;

        % Remove existing data from fth
        obj.data(user_selection(i,1) <= obj.data.dt & obj.data.dt <= user_selection(i,2),:) = [];
        
        % get variable names
        varnam = obj.data.Properties.VariableNames;
        if size(varnam,2) == 3
          varnam = {'dt','swt','spd1','spd2'};
          obj.data.spd2 = NaN(size(obj.data.dt));
          obj.data = renamevars(obj.data, 'spd', 'spd1');
        end
        switch mode
          case 'total'
            obj.data = [obj.data; array2table([dt ones(size(dt))*obj.SWITCH_TOTAL ...
              interp_spd], 'VariableNames', varnam)];
          case 'filtered'
            obj.data = [obj.data; array2table([dt ones(size(dt))*obj.SWITCH_FILTERED ...
              interp_spd], 'VariableNames', varnam)];
          otherwise
            error('Unknown mode.');
        end
      end
      obj.data = sortrows(obj.data, 'dt');
      fprintf('Done\n');
    end
  end
end
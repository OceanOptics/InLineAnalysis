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
    swt_variable = [];
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
      if isempty(cfg.view.spd_variable)
        obj.view.swt_variable = 'swt';
      else
        obj.view.swt_variable = cfg.view.swt_variable;
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
        case 'Inlinino'
          obj.data = iRead(@importInlininoFlowControl, obj.path.raw, obj.path.wk, [obj.model obj.sn '_'],...
                       days2run, 'Inlinino', force_import, ~write, true);
        case 'Inlinino_base'
          obj.data = iRead(@importInlinino_base, obj.path.raw, obj.path.wk, obj.prefix,...
                       days2run, 'Inlinino', force_import, ~write, true);
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
        if ~isdatetime(obj.data.dt)
          obj.data.dt = datetime(obj.data.dt, 'ConvertFrom', 'datenum');
        end
        if ~isdatetime(obj.data.dt)
          obj.data.dt = datetime(obj.data.dt, 'ConvertFrom', 'datenum');
        end
        % dt_st = datenum(floor(datevec(user_selection(i,1))));
        % dt_end = datenum(floor(datevec(user_selection(i,2))));%, 'ConvertFrom','datenum', 'Format', 'yyyy/MM/dd hh:mm:ss') ...
        % dt = (dt_st:1/obj.SAMPLING_FREQUENCY/3600/24:dt_end)';

        % % remove duplicates
        % [~, L, ~] = unique(obj.data.dt,'first');
        % indexToDump = not(ismember(1:numel(obj.data.dt), L));
        % if any(indexToDump)
        %   obj.data(indexToDump, :) = [];
        % end

        obj.data = round_timestamp(obj.data, seconds(0.001));
        user_selection = dateshift(user_selection, 'Start', 'seconds');
        dt = (user_selection(i,1):seconds(1):user_selection(i,2))';
        
        % switch filtered/total
        foo = table();
        foo.dt = dt;
        % id interpolated data to keep only 'true' data
        isnm = ~ismember(foo.dt, obj.data.dt);
        for j = 2:size(obj.data, 2)
          varname = obj.data.Properties.VariableNames{j};
          if contains(varname, 'swt')
            foo.(varname) = interp1(obj.data.dt, double(obj.data.(varname)), foo.dt, 'previous');
          else
            foo.(varname) = interp1(obj.data.dt, obj.data.(varname), foo.dt, 'linear');
          end
          foo.(varname)(isnm) = NaN;
        end
        % Remove existing data from fth
        obj.data(user_selection(i,1) <= obj.data.dt & obj.data.dt <= user_selection(i,2),:) = [];
        switch mode
          case 'total'
            foo.(obj.view.swt_variable) = ones(size(foo.dt))*obj.SWITCH_TOTAL;
          case 'filtered'
            foo.(obj.view.swt_variable) = ones(size(foo.dt))*obj.SWITCH_FILTERED;
          otherwise
            error('Unknown switch change query.');
        end
        obj.data = [obj.data; foo];
        obj.data = sortrows(obj.data, 'dt');
      end
      fprintf('Done\n');
    end
  end
end
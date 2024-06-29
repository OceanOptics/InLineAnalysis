function QCSwitchPosition(instru, flow, level, shift_flow, fdom)
  % Automatically sync flow data to filt data, change switch position and duplicate filter events
  % Author: Guillaume Bourdin
  % Date: June 2021, updated June 2024
  %%
  if strcmp(level, 'raw') && isempty(instru.raw.tsw)
    error('raw data must be loaded to process using expenential fit method for filter calibration')
  end
  if strcmp(level, 'qc') && isempty(instru.qc.tsw)
    error('qc data must be loaded for QCSwitchPosition')
  end
  if isempty(instru.(level).tsw)
    error('No %s %s tsw data loaded', instru.model, level)
  end
  if isempty(instru.(level).fsw)
    error('No %s %s fsw data loaded', instru.model, level)
  end
  if isempty(flow.(level).tsw)
    error('No FLOW %s data loaded', level)
  end
  if strcmp(level, 'raw')
    shif = secondes(1);
  else
    shif = minutes(1);
  end
  flow.(level).tsw = round_timestamp(flow.(level).tsw, shif);
  
  % get sampling period
  flow_prd = median(diff(flow.(level).tsw.dt), 'omitnan');
  flow_prd_dt = median(diff(datetime(flow.(level).tsw.dt, 'ConvertFrom', 'datenum')), 'omitnan');
  % add one line of total at the beginning and the end of timseries in case a filter event is cut
  if flow.(level).tsw.swt(1) > 0
    flow.(level).tsw = [flow.(level).tsw(1,:); flow.(level).tsw];
    flow.(level).tsw.dt(1) = flow.(level).tsw.dt(1) - flow_prd;
    flow.(level).tsw.swt(1) = 0;
  end
  if flow.(level).tsw.swt(end) > 0
    flow.(level).tsw = [flow.(level).tsw; flow.(level).tsw(end,:)];
    flow.(level).tsw.dt(end) = flow.(level).tsw.dt(end) + flow_prd;
    flow.(level).tsw.swt(end) = 0;
  end
  % prepare switch position
  swt = fillmissing(flow.(level).tsw.swt, 'previous');%, 'linear', 'extrap');
  swt = swt > 0;
  % Find switch events from total to filtered
  sel_start = find(swt(1:end-1) == flow.SWITCH_TOTAL & swt(2:end) == flow.SWITCH_FILTERED);
  % Find switch events from filtered to total
  sel_end = find(swt(1:end-1) == flow.SWITCH_FILTERED & swt(2:end) == flow.SWITCH_TOTAL);
  % Verify selections of filtered period
  if sel_start(1) > sel_end(1); sel_end(1) = []; end
  if sel_start(end) > sel_end(end); sel_end(end+1) = size(swt,1); end
  if size(sel_start,1) ~= size(sel_end,1); error('Inconsistent fth data'); end
  % id filt data in filter event
  id_infilt = false(size(instru.(level).fsw.dt));
  for i = 1:size(sel_start,1)
    id_infilt(instru.(level).fsw.dt >= flow.(level).tsw.dt(sel_start(i)) & ...
      instru.(level).fsw.dt <= flow.(level).tsw.dt(sel_end(i))) = true;
  end
  fprintf('%.2f%% of filtered data is included into filtered switch position before QC\n', ...
    sum(id_infilt) / size(instru.(level).fsw.dt, 1) * 100)
  
  if nargin < 2
    error('Not enough input argument')
  elseif nargin == 2
    level = 'qc';
    autoshift = true;
    fdom = [];
  elseif nargin == 3
    autoshift = true;
    fdom = [];
  elseif nargin == 4
    fdom = [];
    if isempty(shift_flow)
      autoshift = true;
    else
      autoshift = false;
    end
  elseif nargin == 5
    if isempty(shift_flow)
      autoshift = true;
    else
      autoshift = false;
    end
  elseif nargin > 5
    error('Too many input arguments')
  end
  % automatic detect of time shift
  if autoshift
    if sum(id_infilt) / size(instru.(level).fsw.dt, 1) < 1
      % Automatically detect how much to shift flow data time so that switch filtered
      % position include the most filtered data in number of samples
      shift_flow_list = (-20*shif:shif:20*shif)';
      matchs = NaN(size(shift_flow_list));
      for sh = 1:size(shift_flow_list, 1)
        % shift flow timeseries in time
        popoflow = shift_timeseries(flow.(level).tsw, shift_flow_list(sh), flow_prd);
        % add one line of total at the beginning and the end of timseries in case a filter event is cut
        if popoflow.swt(1) > 0
          popoflow = [popoflow(1,:); popoflow];
          popoflow.dt(1) = popoflow.dt(1) - flow_prd;
          popoflow.swt(1) = 0;
        end
        if popoflow.swt(end) > 0
          popoflow = [popoflow; popoflow(end,:)];
          popoflow.dt(end) = popoflow.dt(end) + flow_prd;
          popoflow.swt(end) = 0;
        end
        % prepare switch position
        swt = popoflow.swt;
        swt = swt > 0;
        % Find switch events from total to filtered
        sel_start = find(swt(1:end-1) == flow.SWITCH_TOTAL & swt(2:end) == flow.SWITCH_FILTERED);
        % Find switch events from filtered to total
        sel_end = find(swt(1:end-1) == flow.SWITCH_FILTERED & swt(2:end) == flow.SWITCH_TOTAL);
        % id filt data in filter event
        id_infilt = false(size(instru.(level).fsw.dt));
        for i = 1:size(sel_start,1)
          id_infilt(instru.(level).fsw.dt >= popoflow.dt(sel_start(i)) & instru.(level).fsw.dt <= popoflow.dt(sel_end(i))) = true;
        end
        matchs(sh) = sum(id_infilt);
      end
      % get mininum shift
      shift_flow = shift_flow_list(matchs == max(matchs));
      shift_flow = shift_flow(abs(shift_flow) == min(abs(shift_flow)));
      if shift_flow ~= 0
        fprintf('Switch position shifted automatically to best match filtered data: %s(s)\n', shift_flow)
      end
    else
      shift_flow = 0;
    end
  end
  % Apply shift to flow data
  flow.(level).tsw = shift_timeseries(flow.(level).tsw, shift_flow, flow_prd);

  % create filter event duplicate when long period without filter event
  fh = visFlag(instru.raw.tsw, instru.raw.fsw, instru.qc.tsw, [], instru.qc.fsw, [], ...
    instru.view.varname, instru.view.varcol, instru.raw.bad, flow.qc.tsw, flow.view.spd_variable);
  plot(flow.(level).tsw.dt, flow.(level).tsw.swt, '-k');

  if ~isempty(fdom)
    ax1 = gca; % current axes
    ax1_pos = ax1.Position;
    ax2 = axes('Position',ax1_pos, 'YAxisLocation','right',...
        'xColor','none', 'yColor','none','Color', 'none');
    scatter(ax2, fdom.(level).tsw.dt, fdom.(level).tsw.fdom, 30, 'MarkerFaceColor', ...
      [0	205	205]/255, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.8)
    ax2.Color = 'none';
    ax2.XColor = 'none';
    ax2.YColor = 'none';
    linkaxes([ax1 ax2],'x')
    leg = [findobj('Type','Line'); findobj('Type','Scatter')];
    legend(leg, 'switch position (1=filtered | 0=total)', 'Flow speed', 'Binned filtered', ...
      'Binned total', 'Binned fDOM', 'AutoUpdate','off', 'FontSize', 12)
  else
    leg = findobj('Type','Line');
    legend(leg, 'switch position (1=filtered | 0=total)', 'Flow speed', 'Binned filtered', 'Binned total',...
      'AutoUpdate','off', 'FontSize', 12)
  end
  title(['Select filter event to duplicate (press x)' newline 'Select new time slot for filter event duplicated (press s)' newline 'Change switch position to filtered (press f)'  newline 'Change switch position to total (press t)'], ...
    'FontSize', 14)

  [totalswitch, filterswitch, newtime_s, ~, toduplicate_x] = guiSelectOnTimeSeries(fh);
  % check number of entries
  if size(toduplicate_x,1) ~= size(newtime_s,1)
    fprintf('Warning: Inconsistent number of entries (x and s), filter events duplication ignored\n')
    toduplicate_x = [];
    newtime_s = [];
  end
  % duplicate filter events x and s commands
  if size(toduplicate_x, 1) > 0
    toduplicate_x = round_timestamp(toduplicate_x, shif);
  end
  nts = round_timestamp(newtime_s, shif);
  for j = 1:size(toduplicate_x, 1)
    if strcmp(level, 'raw')
      idx_inst_qc = instru.qc.fsw.dt >= toduplicate_x(j, 1) & instru.qc.fsw.dt <= toduplicate_x(j, 2);
      idx_inst_bad = instru.raw.bad.dt >= toduplicate_x(j, 1) & instru.raw.bad.dt <= toduplicate_x(j, 2);
    end
    idx_inst = instru.(level).fsw.dt >= toduplicate_x(j, 1) & instru.(level).fsw.dt <= toduplicate_x(j, 2);
    if any(idx_inst)
      % copy filter event to specific new time slot
      new_filt = instru.(level).fsw(idx_inst, :);
      lag = min(new_filt.dt) - toduplicate_x(j, 1);
      if strcmp(level, 'raw')
        foo_newfilt_dt = [instru.qc.fsw.dt(idx_inst_qc); instru.(level).fsw.dt(idx_inst); ...
          instru.raw.bad.dt(idx_inst_bad)];
      else
        foo_newfilt_dt = instru.(level).fsw.dt(idx_inst);
      end
      delta_dt = min(foo_newfilt_dt) - nts(j);
      new_filt.dt = new_filt.dt - delta_dt + lag;
      % delete filt data at time of new filter event
      id_newfilt = instru.(level).fsw.dt >= min(new_filt.dt) & instru.(level).fsw.dt <= max(new_filt.dt);
      instru.(level).fsw(id_newfilt, :) = [];
      % add new filt
      instru.(level).fsw = [instru.(level).fsw; new_filt];
      % sort by date
      instru.(level).fsw = sortrows(instru.(level).fsw, 'dt');
      if strcmp(level, 'raw')
        % copy qc if level raw to appear on the plot even if not used
        new_filt_qc = instru.qc.fsw(idx_inst_qc, :);
        lag = min(new_filt_qc.dt) - toduplicate_x(j, 1);
        new_filt_qc.dt = new_filt_qc.dt - delta_dt + lag;
        % delete filt data at time of new filter event
        id_newfilt = instru.qc.fsw.dt >= min(new_filt_qc.dt) & instru.qc.fsw.dt <= max(new_filt_qc.dt);
        instru.qc.fsw(id_newfilt, :) = [];
        % add new filt
        instru.qc.fsw = [instru.qc.fsw; new_filt_qc];
        instru.qc.fsw = sortrows(instru.qc.fsw, 'dt');
  
        % copy bad if level raw
        new_filt_bad = instru.raw.bad(idx_inst_bad, :);
        lag = min(new_filt_bad.dt) - toduplicate_x(j, 1);
        new_filt_bad.dt = new_filt_bad.dt - delta_dt + lag;
        % delete filt data at time of new filter event
        id_newfilt = instru.raw.fsw.dt >= min(new_filt_bad.dt) & instru.raw.fsw.dt <= max(new_filt_bad.dt);
        instru.raw.fsw(id_newfilt, :) = [];
        % add new filt
        instru.raw.bad = [instru.raw.bad; new_filt_bad];
        instru.raw.bad = sortrows(instru.raw.bad, 'dt');
      end
  
      % duplicate fdom if fdom not empty
      if ~isempty(fdom)
        idx_inst = fdom.(level).tsw.dt >= toduplicate_x(j, 1) & fdom.(level).tsw.dt <= toduplicate_x(j, 2);
        if any(idx_inst)
          % copy filter event to specific new time slot
          new_fdom = fdom.(level).tsw(idx_inst, :);
          lag = min(new_fdom.dt) - toduplicate_x(j, 1);
          if strcmp(level, 'raw')
            foo_newfdom_dt = [fdom.qc.tsw.dt(idx_inst_qc); fdom.(level).tsw.dt(idx_inst)];
          else
            foo_newfdom_dt = fdom.(level).tsw.dt(idx_inst);
          end
          delta_dt = min(foo_newfdom_dt) - nts(j);
          new_fdom.dt = new_fdom.dt - delta_dt + lag;
          % delete filt data at time of new filter event
          id_newfdom = fdom.(level).tsw.dt >= min(new_fdom.dt) & fdom.(level).tsw.dt <= max(new_fdom.dt);
          fdom.(level).fsw(id_newfdom, :) = [];
          % add new fdom
          fdom.(level).tsw = [fdom.(level).tsw; new_fdom];
          % sort by date
          fdom.(level).tsw = sortrows(fdom.(level).tsw, 'dt');
        end
      end
      % create flow data at new filter time
      if strcmp(level, 'raw')
        newflow_st = round_timestamp(datetime(min([new_filt_qc.dt; new_filt.dt; new_filt_bad.dt]), ...
          'ConvertFrom', 'datenum'), shif);
        newflow_end = round_timestamp(datetime(max([new_filt_qc.dt; new_filt.dt; new_filt_bad.dt]), ...
          'ConvertFrom', 'datenum'), shif);
      else
        newflow_st = round_timestamp(datetime(min(new_filt.dt), 'ConvertFrom', 'datenum'), shif);
        newflow_end = round_timestamp(datetime(max(new_filt.dt), 'ConvertFrom', 'datenum'), shif);
      end
      newflow_dt = (newflow_st-flow_prd_dt:flow_prd_dt:newflow_end+flow_prd_dt)';
      idx_flow_newfilt = flow.(level).tsw.dt >= datenum(min(newflow_dt)-flow_prd_dt) & flow.(level).tsw.dt <= datenum(max(newflow_dt)+flow_prd_dt);
      old_flow = flow.(level).tsw(idx_flow_newfilt,:);
      old_flow.dt = datetime(old_flow.dt, 'ConvertFrom', 'datenum');
      old_flow = round_timestamp(old_flow, shif);
      % create flow table with swt == 0 as first and end row
      new_flow = table('Size', [size(newflow_dt,1)+2 size(flow.(level).tsw, 2)], ...
        'VariableNames', flow.(level).tsw.Properties.VariableNames, ...
        'VariableTypes', ['datetime', repmat({'doubleNaN'}, 1, size(flow.(level).tsw(:,2:end), 2))]);
      new_flow.dt = [min(newflow_dt)-flow_prd_dt; newflow_dt; max(newflow_dt)+flow_prd_dt];
      new_flow.swt = [0; ones(size(newflow_dt,1), 1); 0];
      % put back already existing flow data within the newflow table
      new_flow(ismember(new_flow.dt, old_flow.dt), :) = old_flow;
      new_flow.swt(ismember(new_flow.dt, old_flow.dt)) = 1;
      new_flow.swt([1 end]) = 0;
      new_flow.dt = datenum(new_flow.dt);
      % merge newflow with flow data
      flow.(level).tsw(idx_flow_newfilt,:) = [];
      flow.(level).tsw = [flow.(level).tsw; new_flow];
      % sort by date
      flow.(level).tsw = sortrows(flow.(level).tsw, 'dt');
    end
  end
  % change switch position to total (t)
  for j = 1:size(totalswitch, 1)
    idx_flow = flow.(level).tsw.dt > totalswitch(j, 1) & flow.(level).tsw.dt < totalswitch(j, 2);
    linetoadd_dt = (totalswitch(j, 1):flow_prd:totalswitch(j, 2))';
    new_flow = array2table([linetoadd_dt zeros(size(linetoadd_dt, 1), 1) ...
        NaN(size(linetoadd_dt, 1), size(popoflow, 2) - 2)], ...
        'VariableNames', flow.(level).tsw.Properties.VariableNames);
    if sum(idx_flow) >= 2
      new_flow.(flow.view.spd_variable) = interp1(flow.(level).tsw.dt(idx_flow), flow.(level).tsw.(flow.view.spd_variable)(idx_flow), new_flow.dt, 'linear');
    end
    flow.(level).tsw(idx_flow, :) = [];
    flow.(level).tsw = [flow.(level).tsw; new_flow];
    % sort by date
    flow.(level).tsw = sortrows(flow.(level).tsw, 'dt');
  end
  % change switch position to filtered (f)
  for j = 1:size(filterswitch, 1)
    idx_flow = flow.(level).tsw.dt > filterswitch(j, 1) & flow.(level).tsw.dt < filterswitch(j, 2);
    linetoadd_dt = (filterswitch(j, 1):flow_prd:filterswitch(j, 2))';
    new_flow = array2table([linetoadd_dt ones(size(linetoadd_dt, 1), 1) ...
        NaN(size(linetoadd_dt, 1), size(popoflow, 2) - 2)], ...
        'VariableNames', flow.(level).tsw.Properties.VariableNames);
    if sum(idx_flow) >= 2
      % check for duplicats in flow data and delete
      [~, L, ~] = unique(flow.(level).tsw.dt,'first');
      indexToDump = not(ismember(1:numel(flow.(level).tsw.dt),L));
      if sum(indexToDump) > 0
        fprintf('Warning: %i identical dates in flow data => deleted\n', sum(indexToDump))
        flow.(level).tsw(indexToDump, :) = [];
      end
      new_flow.(flow.view.spd_variable) = interp1(flow.(level).tsw.dt(idx_flow), flow.(level).tsw.(flow.view.spd_variable)(idx_flow), new_flow.dt, 'linear');
    end
    flow.(level).tsw(idx_flow, :) = [];
    flow.(level).tsw = [flow.(level).tsw; new_flow];
    % sort by date
    flow.(level).tsw = sortrows(flow.(level).tsw, 'dt');
  end
  
  % check for duplicats in flow data and delete
  [~, L, ~] = unique(instru.(level).fsw.dt,'first');
  indexToDump = not(ismember(1:numel(instru.(level).fsw.dt),L));
  if sum(indexToDump) > 0
    fprintf('Warning: %i identical dates in instrument data => deleted\n', sum(indexToDump))
    instru.(level).fsw(indexToDump, :) = [];
  end
  [~, L, ~] = unique(flow.(level).tsw.dt,'first');
  indexToDump = not(ismember(1:numel(flow.(level).tsw.dt),L));
  if sum(indexToDump) > 0
    fprintf('Warning: %i identical dates in FLOW data => deleted\n', sum(indexToDump))
    flow.(level).tsw(indexToDump, :) = [];
  end

  % check amount if filt data matching with flow filt events
  swt = fillmissing(flow.(level).tsw.swt, 'previous');%, 'linear', 'extrap');
  swt = swt > 0;
  % Find switch events from total to filtered
  sel_start = find(swt(1:end-1) == flow.SWITCH_TOTAL & swt(2:end) == flow.SWITCH_FILTERED);
  % Find switch events from filtered to total
  sel_end = find(swt(1:end-1) == flow.SWITCH_FILTERED & swt(2:end) == flow.SWITCH_TOTAL);
  % Verify selections of filtered period
  if sel_start(1) > sel_end(1); sel_end(1) = []; end
  if sel_start(end) > sel_end(end); sel_end(end+1) = size(swt,1); end
  if size(sel_start,1) ~= size(sel_end,1); error('Inconsistent fth data'); end
  % id filt data in filter event
  id_infilt = false(size(instru.(level).fsw.dt));
  for i = 1:size(sel_start,1)
    id_infilt(instru.(level).fsw.dt >= flow.(level).tsw.dt(sel_start(i)) & ...
      instru.(level).fsw.dt <= flow.(level).tsw.dt(sel_end(i))) = true;
  end
  fprintf('%.2f%% of filtered data is included into filtered switch position before QC\n', ...
    sum(id_infilt) / size(instru.(level).fsw.dt, 1) * 100)
end


function data_out = shift_timeseries(data_in, shift, flow_prd)
  if shift == 0 
    data_out = data_in;
  else
    if shift > 0
      linetoadd_dt = (min(data_in.dt)-datenum(shift):flow_prd:min(data_in.dt)-flow_prd)';
    elseif shift < 0
      linetoadd_dt = (max(data_in.dt)+flow_prd:flow_prd:max(data_in.dt)-datenum(shift))';
    end
    linetoadd = table('Size', [size(linetoadd_dt,1) size(data_in, 2)], ...
      'VariableNames', data_in.Properties.VariableNames, ...
      'VariableTypes', repmat({'doubleNaN'}, 1, size(data_in, 2)));
    linetoadd.dt = linetoadd_dt;
    linetoadd.swt = zeros(size(linetoadd_dt,1), 1);
    % merge new line to data_in timeseries
    [h, m, s] = hms(shift);
    shift_sz = sum([h, m, s]);
    if shift_sz > 0
      data_out = [data_in(1:end-abs(shift_sz), :); linetoadd];
    else
      data_out = [data_in(abs(shift_sz)+1:end, :); linetoadd];
    end
    data_out = sortrows(data_out, 'dt');
  end
end

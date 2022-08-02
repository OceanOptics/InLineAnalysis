function QCSwitchPosition(instru, flow, level, shift_flow)
% Automatically sync flow data to filt data, change switch position and
% duplicate filter events
%
% Author: Guillaume Bourdin
% Date: Jun 2021
%
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
  shif = 'second';
else
  shif = 'minute';
end
filtdt = datetime(instru.(level).fsw.dt, 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
filtdt = datenum(dateshift(filtdt, 'start', shif));
flowdt = datetime(flow.(level).tsw.dt(flow.(level).tsw.swt == 1), 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
flowdt = datenum(dateshift(flowdt, 'start', shif));
fprintf('%.2f%% of filtered data is included into filtered switch position before QC\n', ...
  sum(ismember(filtdt, flowdt)) / size(filtdt,1) * 100)
d_filtdt = median(diff(filtdt), 'omitnan');
d_flowdt = median(diff(flowdt), 'omitnan');
if nargin < 2
  error('Not enough input argument')
elseif nargin == 2
  level = 'qc';
elseif nargin < 4
  % Automatically detect how much to shift flow data time so that switch filtered
  % position include the most filtered data in number of minutes
  shift_flow_list = (-20:20)';
  matchs = NaN(size(shift_flow_list));
  for sh = 1:size(shift_flow_list, 1)
    popoflow = flow.(level).tsw(flow.(level).tsw.swt == 1, :);
    if shift_flow_list(sh) > 0
      linetoadd_dt = (min(popoflow.dt) - datenum(minutes(shift_flow_list(sh))):d_flowdt:...
        min(popoflow.dt) - d_flowdt)';
      linetoadd = array2table([linetoadd_dt zeros(size(linetoadd_dt, 1), size(popoflow, 2) - 1)], ...
        'VariableNames', popoflow.Properties.VariableNames);
      popoflow = [repmat(linetoadd, shift_flow_list(sh), 1); popoflow(1:end-shift_flow_list(sh),:)];
    elseif shift_flow_list(sh) < 0
      linetoadd_dt = (max(popoflow.dt) + d_flowdt:d_flowdt:...
        max(popoflow.dt) - datenum(minutes(shift_flow_list(sh))))';
      linetoadd = array2table([linetoadd_dt zeros(size(linetoadd_dt, 1), size(popoflow, 2) - 1)], ...
        'VariableNames', popoflow.Properties.VariableNames);
      popoflow = [popoflow(abs(shift_flow_list(sh))+1:end, :); repmat(linetoadd, abs(shift_flow_list(sh)), 1)];
    end
    popoflow.dt = popoflow.dt + datenum(minutes(shift_flow_list(sh)));
    filtdt = datetime(instru.(level).fsw.dt, 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
    filtdt = datenum(dateshift(filtdt, 'start', shif));
    flowdt = datetime(popoflow.dt(popoflow.swt == 1), 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
    flowdt = datenum(dateshift(flowdt, 'start', shif));
    matchs(sh) = sum(ismember(filtdt, flowdt));
  end
  shift_flow = shift_flow_list(matchs == max(matchs));
  shift_flow = min(shift_flow);
  if shift_flow ~= 0
    fprintf('Switch position shifted automatically to best match filtered data: %i minute(s)\n', shift_flow)
  end
elseif nargin > 4
  error('Too many input arguments')
end
% Apply shift to flow data
shift_flow = round(shift_flow);
if shift_flow > 0
  linetoadd_dt = (min(flow.(level).tsw.dt) - datenum(minutes(shift_flow)):d_flowdt:...
    min(flow.(level).tsw.dt) - d_flowdt)';
  linetoadd = array2table([linetoadd_dt zeros(size(linetoadd_dt, 1), size(flow.(level).tsw, 2) - 1)], ...
    'VariableNames', flow.(level).tsw.Properties.VariableNames);
  flow.(level).tsw = [repmat(linetoadd, shift_flow, 1); flow.(level).tsw(1:end-shift_flow,:)];
elseif shift_flow < 0
  linetoadd_dt = (max(flow.(level).tsw.dt) + d_flowdt:d_flowdt:...
    max(flow.(level).tsw.dt) - datenum(minutes(shift_flow)))';
  linetoadd = array2table([linetoadd_dt zeros(size(linetoadd_dt, 1), size(flow.(level).tsw, 2) - 1)], ...
    'VariableNames', flow.(level).tsw.Properties.VariableNames);
  flow.(level).tsw = [flow.(level).tsw(abs(shift_flow)+1:end, :); repmat(linetoadd, abs(shift_flow), 1)];
end
flow.(level).tsw.dt = flow.(level).tsw.dt + datenum(minutes(shift_flow));
% create filter event duplicate when long period without filter event
fh = visFlag(instru.raw.tsw, instru.raw.fsw, instru.qc.tsw, [], instru.qc.fsw, [], ...
  instru.view.varname, instru.view.varcol, instru.raw.bad, flow.qc.tsw, flow.view.spd_variable );
leg = plot(flow.(level).tsw.dt, flow.(level).tsw.swt, '-k');
title(['Select filter event to duplicate (press x)' newline 'Select new time slot for filter event duplicated (press s)' newline 'Change switch position to filtered (press f)'  newline 'Change switch position to total (press t)'], ...
  'FontSize', 14)
legend(leg, 'switch position (1=filtered | 0=total)','AutoUpdate','off', 'FontSize', 12)
[totalswitch, filterswitch, newtime_s, ~, toduplicate_x] = guiSelectOnTimeSeries(fh);
% check number of entries
if size(toduplicate_x,1) ~= size(newtime_s,1)
  fprintf('Warning: Inconsistent number of entries (f and s), filter events duplication ignored\n')
  toduplicate_x = [];
  newtime_s = [];
end
% duplicate filter events f and s commands
for j = 1:size(toduplicate_x, 1)
  if strcmp(level, 'raw')
    idx_inst_qc = instru.qc.fsw.dt >= toduplicate_x(j, 1) & instru.qc.fsw.dt <= toduplicate_x(j, 2);
    idx_inst_bad = instru.raw.bad.dt >= toduplicate_x(j, 1) & instru.raw.bad.dt <= toduplicate_x(j, 2);
  end
  idx_inst = instru.(level).fsw.dt >= toduplicate_x(j, 1) & instru.(level).fsw.dt <= toduplicate_x(j, 2);
  if any(idx_inst)
    % copy filter event to specific new time slot
    new_filt = instru.(level).fsw(idx_inst, :);
    nts = dateshift(datetime(newtime_s(j), 'ConvertFrom', 'datenum'), 'start', shif);
    if strcmp(level, 'raw')
      foo_newfilt_dt = [instru.qc.fsw.dt(idx_inst_qc); instru.(level).fsw.dt(idx_inst); ...
        instru.raw.bad.dt(idx_inst_bad)];
    else
      foo_newfilt_dt = instru.(level).fsw.dt(idx_inst);
    end
    delta_dt = min(foo_newfilt_dt) - datenum(nts);
    new_filt.dt = new_filt.dt - delta_dt;
    instru.(level).fsw = [instru.(level).fsw; new_filt];
    % sort by date
    instru.(level).fsw = sortrows(instru.(level).fsw, 'dt');
    if strcmp(level, 'raw')
      % copy qc if level raw
      new_filt_qc = instru.qc.fsw(idx_inst_qc, :);
      new_filt_qc.dt = new_filt_qc.dt - delta_dt;
      instru.qc.fsw = [instru.qc.fsw; new_filt_qc];
      instru.qc.fsw = sortrows(instru.qc.fsw, 'dt');
      % copy bad if level raw
      new_filt_bad = instru.raw.bad(idx_inst_bad, :);
      new_filt_bad.dt = new_filt_bad.dt - delta_dt;
      instru.raw.bad = [instru.raw.bad; new_filt_bad];
      instru.raw.bad = sortrows(instru.raw.bad, 'dt');
      % get end switch
      end_swt = max(instru.(level).fsw.dt(idx_inst)) - delta_dt;
    end
    idx_flow_newfilt = flow.(level).tsw.dt >= min(foo_newfilt_dt - delta_dt)  - d_filtdt & ...
      flow.(level).tsw.dt <= max(foo_newfilt_dt - delta_dt) + d_filtdt;
    flow.(level).tsw(idx_flow_newfilt,:) = [];
    % create flow table
    linetoadd_dt = (min(foo_newfilt_dt - delta_dt) - d_filtdt:d_flowdt:...
      max(foo_newfilt_dt - delta_dt) + d_filtdt)';
    flow_swt_offst = array2table([linetoadd_dt(1) - d_flowdt', 0, NaN(1, size(popoflow, 2) - 2)], ...
      'VariableNames', flow.(level).tsw.Properties.VariableNames);
    flow_swt_offend = array2table([linetoadd_dt(end) + d_flowdt', 0, NaN(1, size(popoflow, 2) - 2)], ...
      'VariableNames', flow.(level).tsw.Properties.VariableNames);
    new_flow = array2table([linetoadd_dt ones(size(linetoadd_dt, 1), 1) ...
      NaN(size(linetoadd_dt, 1), size(popoflow, 2) - 2)], ...
      'VariableNames', flow.(level).tsw.Properties.VariableNames);
    if strcmp(level, 'raw')
      new_flow.swt(new_flow.dt > end_swt) = 0;
    end
    flow.(level).tsw = [flow.(level).tsw; flow_swt_offst; new_flow; flow_swt_offend];
    % sort by date
    flow.(level).tsw = sortrows(flow.(level).tsw, 'dt');
  end
end
% change switch position to total (t)
for j = 1:size(totalswitch, 1)
  idx_flow = flow.(level).tsw.dt > totalswitch(j, 1) & flow.(level).tsw.dt < totalswitch(j, 2);
  linetoadd_dt = (totalswitch(j, 1):d_flowdt:totalswitch(j, 2))';
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
  linetoadd_dt = (filterswitch(j, 1):d_flowdt:filterswitch(j, 2))';
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
%   nbtoswitch = size(totalswitch(j, 1):d_flowdt:totalswitch(j, 2),2);
%   if sum(idx_flow) < nbtoswitch % replace all FLOW data with new table
%     closest_position = flow.(level).tsw.swt(abs(flow.(level).tsw.dt(~idx_flow) - mean(totalswitch(j, :))) == ...
%       min(abs(flow.(level).tsw.dt(~idx_flow) - mean(totalswitch(j, :))))); % get the closest switch position in time
%     linetoadd_dt = (totalswitch(j, 1):d_flowdt:totalswitch(j, 2))';
%     new_flow = array2table([linetoadd_dt repmat(closest_position, size(linetoadd_dt, 1), 1) ...
%       NaN(size(linetoadd_dt, 1), size(popoflow, 2) - 2)], ...
%       'VariableNames', flow.(level).tsw.Properties.VariableNames);
%     new_flow.(flow.view.spd_variable) = interp1(flow.(level).tsw.dt(idx_flow), flow.(level).tsw.(flow.view.spd_variable )(idx_flow), new_flow.dt, 'linear');
%     flow.(level).tsw(idx_flow, :) = [];
%     flow.(level).tsw = [flow.(level).tsw; new_flow];
%     % sort by date
%     flow.(level).tsw = sortrows(flow.(level).tsw, 'dt');
%   else
% %     flow.(level).tsw.swt(idx_flow) = 1;
%     if sum(flow.(level).tsw.swt(idx_flow) == 0) > sum(flow.(level).tsw.swt(idx_flow) > 0)
%       flow.(level).tsw.swt(idx_flow) = 1;
%     else
%       flow.(level).tsw.swt(idx_flow) = 0;
%     end
%   end
% end

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
filtdt = datetime(instru.(level).fsw.dt, 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
filtdt = datenum(dateshift(filtdt, 'start', shif));
flowdt = datetime(flow.(level).tsw.dt(flow.(level).tsw.swt == 1), 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
flowdt = datenum(dateshift(flowdt, 'start', shif));
fprintf('%.2f%% of filtered data is included into filtered switch position after QC\n', ...
  sum(ismember(filtdt, flowdt)) / size(filtdt,1) * 100)
end
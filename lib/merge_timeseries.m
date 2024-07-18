function merged_data = merge_timeseries(data, data_tomerge, vars, suffix, replace_consecutive_nan)
  if ~isdatetime(data.dt)
    data_wasdatenum = true;
    data.dt = datetime(data.dt, 'ConvertFrom', 'datenum');
  end
  if ~isdatetime(data_tomerge.dt)
    data_tomerge_wasdatenum = true;
    data_tomerge.dt = datetime(data_tomerge.dt, 'ConvertFrom', 'datenum');
  end
  
  % order data
  data = sortrows(data, 'dt');
  data_tomerge = sortrows(data_tomerge, 'dt');
  % merge variable on timeseries
  data = clean_dt(data);
  data_tomerge = clean_dt(data_tomerge);
  % create table without gap in dt to remove correct number of consecutive NaNs
  merged_data = table();
  merged_data.dt = (min(data.dt):minutes(1):max(data.dt))';
  idex = ismember(merged_data.dt, data.dt);
  var_todo = data.Properties.VariableNames(~strcmp(data.Properties.VariableNames, 'dt'));
  for i = var_todo
    if iscell(data.(i{:}))
      merged_data.(i{:}) = repmat({''}, size(merged_data.dt, 1), size(data.(i{:}), 2));
    elseif isdatetime(data.(i{:}))
      merged_data.(i{:}) = NaT(size(merged_data.dt, 1), size(data.(i{:}), 2));
    else
      merged_data.(i{:}) = NaN(size(merged_data.dt, 1), size(data.(i{:}), 2));
    end
    merged_data.(i{:})(idex, :) = data.(i{:});
  end
  merged_data.Properties.VariableUnits = data.Properties.VariableUnits;
  merged_data.Properties.VariableDescriptions = data.Properties.VariableDescriptions;

  if nargin == 3
    suffix = '';
    replace_consecutive_nan = '';
  elseif nargin == 4
    if ~isempty(suffix)
      suffix = ['_' suffix];
    end
    replace_consecutive_nan = '';
  end

  % flip variable to add to match the order input
  vars = flip(vars);
  % find id of data matching
  id = ismember(merged_data.dt, data_tomerge.dt);
  % merge variables
  for j = 1:size(vars, 2)
    % create variable if not already in target table
    if ~any(strcmp(merged_data.Properties.VariableNames, [vars{j} suffix]))
      if iscell(data_tomerge.(vars{j}))
        merged_data = addvars(merged_data, repmat({''}, size(merged_data.dt)), 'NewVariableNames', [vars{j} suffix], 'After', 'dt');
      elseif isdatetime(data_tomerge.(vars{j}))
        merged_data = addvars(merged_data, NaT(size(merged_data.dt)), 'NewVariableNames', [vars{j} suffix], 'After', 'dt');
      else
        merged_data = addvars(merged_data, NaN(size(merged_data.dt)), 'NewVariableNames', [vars{j} suffix], 'After', 'dt');
      end
    end
    % find only id of variable with nan to replace
    if iscell(data_tomerge.(vars{j}))
      idnan = cellfun('isempty', merged_data.([vars{j} suffix]));
    else
      idnan = isnan(merged_data.([vars{j} suffix]));
    end
    id_tomerge = ismember(data_tomerge.dt, merged_data.dt(idnan));
    % merge data
    merged_data.([vars{j} suffix])(id & idnan, :) = data_tomerge.(vars{j})(id_tomerge);
    % id merged variable location
    idvar = strcmp(merged_data.Properties.VariableNames, [vars{j} suffix]);
    % merge variable unit
    if size(data_tomerge.Properties.VariableUnits, 2) == size(data_tomerge.Properties.VariableNames, 2)
      if isempty(merged_data.Properties.VariableUnits{idvar})
        merged_data.Properties.VariableUnits(idvar) = data_tomerge.Properties.VariableUnits(strcmp(data_tomerge.Properties.VariableNames, vars{j}));
      end
    end
    % merge variable description
    if size(data_tomerge.Properties.VariableDescriptions, 2) == size(data_tomerge.Properties.VariableNames, 2)
      if isempty(merged_data.Properties.VariableDescriptions{idvar})
        merged_data.Properties.VariableDescriptions(idvar) = data_tomerge.Properties.VariableDescriptions(strcmp(data_tomerge.Properties.VariableNames, vars{j}));
      end
    end
    % linearly interpolate missing lat/lon when consecutive missing data < replace_consecutive_nan
    if ~isempty(replace_consecutive_nan)
      missing_data = isnan(merged_data.([vars{j} suffix]));
      t = [true; diff(missing_data) ~= 0];
      k = diff(find([t; true])) .* missing_data(t);
      long_nan = k(cumsum(t)) > replace_consecutive_nan;
      if all(long_nan)
        long_nan = false(size(long_nan));
      end
      merged_data.([vars{j} suffix]) = fillmissing(merged_data.([vars{j} suffix]), 'linear', 'SamplePoints', merged_data.dt);
      merged_data.([vars{j} suffix])(long_nan) = NaN;
    end
  end
  merged_data = merged_data(idex, :);
  if data_wasdatenum && data_tomerge_wasdatenum
    merged_data.dt = datenum(merged_data.dt);
  end
end

function cleaned_tbl = clean_dt(data)
  if isempty(data.Properties.VariableUnits)
    data.Properties.VariableUnits = repmat({''}, size(data.Properties.VariableNames));
  end
  if isempty(data.Properties.VariableDescriptions)
    data.Properties.VariableDescriptions = repmat({''}, size(data.Properties.VariableNames));
  end
  data = sortrows(data, 'dt');
  cleaned_tbl = data;
  % round timestamp to minute start
  cleaned_tbl.dt = dateshift(cleaned_tbl.dt, 'start', 'minute');
  % remove duplicates
  [~, L, ~] = unique(cleaned_tbl.dt,'first');
  indexToDump = not(ismember(1:numel(cleaned_tbl.dt), L));
  if sum(indexToDump) > 0
    % fprintf('Warning: %i identical dates => ignored\n', sum(indexToDump))
    cleaned_tbl(indexToDump, :) = [];
  end
end
function data_out = round_timestamp(data_in)
  % force timestamp to be synched
  data_out = table();
  % make sure data_in.dt in rounded to the time binning frequency
  if isdatetime(data_in.dt)
    datetime_data_in_dt = data_in.dt;
    datenum_in = false;
  else
    datetime_data_in_dt = datetime(data_in.dt, 'ConvertFrom', 'datenum');
    datenum_in = true;
  end
  % get time binning frequency
  Tbin_data_in = median(diff(datetime_data_in_dt));
  if Tbin_data_in >= hours(1)
    data_out.dt = dateshift(datetime_data_in_dt, 'start', 'hours');
  elseif Tbin_data_in >= minutes(1)
    % round start/end time to minute
    data_out.dt = dateshift(datetime_data_in_dt, 'start', 'minutes');
  elseif Tbin_data_in >= seconds(1)
    % round start/end time to seconds
    data_out.dt = dateshift(datetime_data_in_dt, 'start', 'seconds');
  else
    error('automatic detection of sampling rate detected a frequency not supported: check round_timestamp function in processACS.m')
  end
  if datenum_in
    data_out.dt = datenum(data_out.dt);
  end 
  % remove duplicates
  [~, L, ~] = unique(data_in.dt,'first');
  indexToDump = not(ismember(1:numel(data_in.dt), L));
  if sum(indexToDump) > 0
    data_in(indexToDump, :) = [];
  end
  % remove duplicates
  [~, L, ~] = unique(data_out.dt,'first');
  indexToDump = not(ismember(1:numel(data_out.dt), L));
  if sum(indexToDump) > 0
    data_out(indexToDump, :) = [];
  end
  % interpolate data on rounded datetime
  vars = data_in.Properties.VariableNames;
  vars(strcmp(vars, 'dt')) = [];
  for v = vars
    data_out.(v{:}) = interp1(data_in.dt, double(data_in.(v{:})), data_out.dt, 'linear');
    if islogical(data_in.(v{:}))
      data_out.(v{:}) = data_out.(v{:}) > 0;
    end
  end
  if datenum_in
    data_in.dt = datenum(data_in.dt);
  end
end
function data_out = round_timestamp(data_in, Tbin_data_in)
  % force timestamp to be synched
  data_out = table();
  % make sure data_in.dt in rounded to the time binning frequency
  if istable(data_in)
    istab = true;
  else
    istab = false;
    data_in = table(data_in, 'VariableNames', {'dt'});
  end
  if isdatetime(data_in.dt)
    datetime_data_in_dt = data_in.dt;
    datenum_in = false;
  else
    datetime_data_in_dt = datetime(data_in.dt, 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
    datenum_in = true;
  end
  if nargin < 2
    % get time binning frequency
    Tbin_data_in = median(diff(datetime_data_in_dt));
  end
  if Tbin_data_in >= hours(1)
    data_out.dt = dateshift(datetime_data_in_dt, 'start', 'hours', 'nearest');
  elseif Tbin_data_in >= minutes(1)
    % round start/end time to minute
    data_out.dt = dateshift(datetime_data_in_dt, 'start', 'minutes', 'nearest');
  elseif Tbin_data_in >= seconds(1)
    % round start/end time to seconds
    data_out.dt = dateshift(datetime_data_in_dt, 'start', 'seconds', 'nearest');
  else
    error('automatic detection of sampling rate detected a frequency not supported: check round_timestamp function in "lib" directory')
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
    if sum(~isnan(data_in.(v{:}))) >= 2
      data_out.(v{:}) = interp1(data_in.dt, double(data_in.(v{:})), data_out.dt, 'linear');
      if islogical(data_in.(v{:}))
        data_out.(v{:}) = data_out.(v{:}) > 0;
      end
    else
      data_out.(v{:}) = data_in.(v{:});
    end
  end
  if datenum_in
    data_in.dt = datenum(data_in.dt);
  end
  if ~istab
    data_out = data_out.dt;
  end
end
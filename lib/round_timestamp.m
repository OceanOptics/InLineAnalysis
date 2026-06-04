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
    % get time binning frequency if not input
    foo = diff(datetime_data_in_dt);
    foo(foo == duration(0,0,0)) = [];
    Tbin_data_in = median(foo);
  end
  if Tbin_data_in >= years(1)
    % round timestamp to year or multiple years
    if Tbin_data_in == years(1)
      data_out.dt = dateshift(datetime_data_in_dt, 'start', 'year', 'nearest');
    else
      t_rounded = datetime_data_in_dt;
      t_rounded.Year = round(year(datetime_data_in_dt) / years(Tbin_data_in)) * years(Tbin_data_in);
      data_out.dt = dateshift(t_rounded, 'start', 'year');
    end
  elseif iscalendarduration(Tbin_data_in)
    % round timestamp to month or multiple months
    if Tbin_data_in == calmonths(1)
      data_out.dt = dateshift(datetime_data_in_dt, 'start', 'month', 'nearest');
    else
      baseDate = datetime(2000, 1, 1);
      data_out.dt = baseDate + calmonths(calmonths(Tbin_data_in) * ...
        round(calmonths(between(baseDate, datetime_data_in_dt, 'months')) / ...
        calmonths(Tbin_data_in)));
    end
  elseif Tbin_data_in >= days(1)
    % round timestamp to days or multiple days
    if Tbin_data_in == days(1)
      data_out.dt = dateshift(datetime_data_in_dt, 'start', 'days', 'nearest');
    else
      data_out.dt = datetime(round(posixtime(datetime_data_in_dt) / ...
        (days(Tbin_data_in) * 24*60*60)) * (days(Tbin_data_in) * 24*60*60), ...
        'ConvertFrom', 'posixtime','Format', datetime_data_in_dt.Format);
    end
  elseif Tbin_data_in >= hours(1)
    % round timestamp to hours or multiple hours
    if Tbin_data_in == hours(1)
      data_out.dt = dateshift(datetime_data_in_dt, 'start', 'hours', 'nearest');
    else
      data_out.dt = datetime(round(posixtime(datetime_data_in_dt) / ...
        (hours(Tbin_data_in) * 60*60)) * (hours(Tbin_data_in) * 60*60), ...
        'ConvertFrom', 'posixtime','Format', datetime_data_in_dt.Format);
    end
  elseif Tbin_data_in >= minutes(1)
    % round timestamp to minute or multiple minutes
    if Tbin_data_in == minutes(1)
      data_out.dt = dateshift(datetime_data_in_dt, 'start', 'minutes', 'nearest');
    else
      data_out.dt = datetime(round(posixtime(datetime_data_in_dt) / ...
        (minutes(Tbin_data_in) * 60)) * (minutes(Tbin_data_in) * 60), ...
        'ConvertFrom', 'posixtime','Format', datetime_data_in_dt.Format);
    end
  elseif Tbin_data_in >= seconds(1)
    % round timestamp to seconds or multiple minutes
    if Tbin_data_in == seconds(1)
      data_out.dt = dateshift(datetime_data_in_dt,'start','minute') + seconds(round(second(datetime_data_in_dt),1));
    else
      data_out.dt = datetime(round(posixtime(datetime_data_in_dt) / ...
        seconds(Tbin_data_in)) * seconds(Tbin_data_in), ...
        'ConvertFrom', 'posixtime','Format', datetime_data_in_dt.Format);
    end
  elseif Tbin_data_in >= seconds(0.1)
    % round timestamp to 0.1 seconds
    data_out.dt = dateshift(datetime_data_in_dt, 'start', 'minute') + seconds(round(second(datetime_data_in_dt), 2));
  elseif Tbin_data_in >= seconds(0.01)
    % round timestamp to 0.01 seconds
    data_out.dt = dateshift(datetime_data_in_dt, 'start', 'minute') + seconds(round(second(datetime_data_in_dt), 3));
  elseif Tbin_data_in >= seconds(0.001)
    % round timestamp to 0.001 seconds
    data_out.dt = dateshift(datetime_data_in_dt, 'start', 'minute') + seconds(round(second(datetime_data_in_dt), 4));
  elseif Tbin_data_in >= seconds(0.0001)
    % round timestamp to 0.0001 seconds
    data_out.dt = dateshift(datetime_data_in_dt, 'start', 'minute') + seconds(round(second(datetime_data_in_dt), 5));
  elseif Tbin_data_in >= seconds(0.00001)
    % round timestamp to 0.00001 seconds
    data_out.dt = dateshift(datetime_data_in_dt, 'start', 'minute') + seconds(round(second(datetime_data_in_dt), 6));
  else
    error('automatic detection of sampling rate detected a frequency not supported: check round_timestamp function in "lib" directory')
  end
  if datenum_in
    data_out.dt = datenum(data_out.dt);
  end
  % remove duplicates
  [~, L, ~] = unique(data_in.dt,'first');
  indexToDump = not(ismember(1:numel(data_in.dt), L));
  if any(indexToDump)
    data_in(indexToDump, :) = [];
  end
  data_in = sortrows(data_in, 'dt');
  % remove duplicates
  [~, L, ~] = unique(data_out.dt,'first');
  indexToDump = not(ismember(1:numel(data_out.dt), L));
  if any(indexToDump)
    data_out(indexToDump, :) = [];
  end
  data_out = sortrows(data_out, 'dt');
  % remove near duplicates
  data_in([diff(data_in.dt,[],1); repmat(Tbin_data_in, 1, size(data_in.dt,2))] < seconds(0.00001)/2, :) = [];
  data_out([diff(data_out.dt,[],1); repmat(Tbin_data_in, 1, size(data_out.dt,2))] < seconds(0.00001)/2, :) = [];
  % interpolate data on rounded datetime
  vars = data_in.Properties.VariableNames;
  vars(strcmp(vars, 'dt')) = [];
  for v = vars
    if isdatetime(data_in.(v{:}))
      data_out.(v{:}) = NaT(size(data_out, 1), size(data_in.(v{:}), 2));
      idint = sum(~isnat(data_in.(v{:}))) >= 2;
      % interpolate all columns with at least 2 non-nan values
      if any(idint)
        data_out.(v{:})(:,idint) = interp1(data_in.dt, data_in.(v{:})(:,idint), data_out.dt, 'linear');
      end
      % merge columns with 1 non-nan data
      if any(~idint)
        idjoin = sum(~isnat(data_in.(v{:}))) == 1;
        % if only one row just copy data
        if size(data_in.(v{:}),1) == 1
          data_out.(v{:})(:,idjoin) = data_in.(v{:})(:,idjoin);
        elseif any(idjoin)
          % if more than one row, find the closest data_out.dt to the data_in.dt of the non-nan data
          for y = find(idjoin)
            id = abs(data_in.dt(~isnat(data_in.(v{:})(:,y))) - data_out.dt) == min(abs(data_in.dt(~isnat(data_in.(v{:})(:,y))) - data_out.dt));
            data_out.(v{:})(id,y) = data_in.(v{:})(~isnat(data_in.(v{:})(:,y)),y);
          end
        end
      end
    else
      data_out.(v{:}) = NaN(size(data_out, 1), size(data_in.(v{:}), 2));
      idint = sum(~isnan(data_in.(v{:}))) >= 2;
      % interpolate all columns with at least 2 non-nan data
      if any(idint)
        if any(contains(v{:}, {'wind_dir','heading','cog','course_over_ground','wind_direction','longitude'})) | strcmp(v{:}, 'lon')
          data_out.(v{:})(:,idint) = interp1_circular(data_in.dt, double(data_in.(v{:})(:,idint)), data_out.dt, 'linear');
        else
          data_out.(v{:})(:,idint) = interp1(data_in.dt, double(data_in.(v{:})(:,idint)), data_out.dt, 'linear');
        end
      end
      % merge columns with 1 non-nan value
      if any(~idint)
        idjoin = sum(~isnan(data_in.(v{:}))) == 1;
        % if only one row just copy data
        if size(data_in.(v{:}),1) == 1
          data_out.(v{:})(:,idjoin) = data_in.(v{:})(:,idjoin);
        elseif any(idjoin)
          % if more than one row, find the closest data_out.dt to the data_in.dt of the non-nan data
          for y = find(idjoin)
            id = abs(data_in.dt(~isnan(data_in.(v{:})(:,y))) - data_out.dt) == min(abs(data_in.dt(~isnan(data_in.(v{:})(:,y))) - data_out.dt));
            data_out.(v{:})(id,y) = data_in.(v{:})(~isnan(data_in.(v{:})(:,y)),y);
          end
        end
      end
      % convert data_out to logical if data_in is logical
      if islogical(data_in.(v{:}))
        data_out.(v{:}) = data_out.(v{:}) > 0;
      end
    end
  end
  if datenum_in
    data_in.dt = datenum(data_in.dt);
  end
  if ~istab
    data_out = data_out.dt;
  end
end
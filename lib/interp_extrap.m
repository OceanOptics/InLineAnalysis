function data_out = interp_extrap(data_in, dt_vector, var_tointerp, max_missing_length, extrap_bool, interp_method, extrap_method)
  % Function to interpolate only gaps smaller than "max_missing_length" and
  % to interpolate and extrapolate data using different methods.
  % Author: Guillaume Bourdin
  % Date: 2024-08-29
  %%
  if nargin < 4
    max_missing_length = 30;
    extrap_bool = true;
    interp_method = 'linear';
    extrap_method = 'nearest';
  elseif nargin < 5
    extrap_bool = true;
    interp_method = 'linear';
    extrap_method = 'nearest';
  elseif nargin < 6
    interp_method = 'linear';
    extrap_method = 'nearest';
  elseif nargin < 7
    extrap_method = 'nearest';
  end
  % remove row full of NaNs
  data_in(all(isnan(data_in.(var_tointerp)), 2), :) = [];
  % convert dt_vector to datetime
  if ~isdatetime(dt_vector)
    dt_vector = datetime(dt_vector, 'ConvertFrom', 'datenum');
  end
  dt = (min(dt_vector):median(diff(dt_vector)):max(dt_vector))';
  % id extrapolation
  extrapolated_id = isnan(interp1(data_in.dt, data_in.(var_tointerp), dt, 'nearest'));
  % find consecutive NaN longer than max_missing_length
  if ~isempty(max_missing_length) && any(any(isnan(data_in.(var_tointerp)), 2))
    missing_data = ~ismember(dt, data_in.dt);
    t = [true; diff(missing_data) ~= 0];
    k = diff(find([t; true])) .* missing_data(t);
    long_nan = k(cumsum(t)) > max_missing_length;
  end
  % interpolate
  data_out = interp1(data_in.dt, data_in.(var_tointerp), dt, interp_method);
  % replace extrapolated data by NaN
  data_out(extrapolated_id) = NaN;
  % fill missing data with extrapolation method
  if extrap_bool
    data_out = fillmissing(data_out, extrap_method, 'SamplePoints', dt);
  end
  % replace interpolated values over gaps > max_missing_length by NaN
  if ~isempty(max_missing_length) && any(any(isnan(data_in.(var_tointerp)), 2))
    if any(long_nan)
      data_out(long_nan, :) = NaN;
    end
  end
  data_out = interp1(dt, data_out, dt_vector, 'nearest');
end
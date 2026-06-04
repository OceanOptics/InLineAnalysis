function data_filled = fillmissing_circular(data_deg, method, sample_points)
  if nargin < 4
    method = 'linear';
  end
  % Unroll the longitude into X and Y vector coordinates
  x = cosd(data_deg);
  y = sind(data_deg);
  
  % Perform standard linear interpolation on both components
  if nargin < 5
    x_filled = fillmissing(x, method);
    y_filled = fillmissing(y, method);
  else
    x_filled = fillmissing(x, method, 'SamplePoints', sample_points);
    y_filled = fillmissing(y, method, 'SamplePoints', sample_points);
  end
  
  % Reconstruct the angles from the interpolated components
  % and force the output back into the -180 to 180 degree range
  data_filled = atan2d(y_filled, x_filled);
end
function data_degq = interp1_circular(x, data_deg, xq, method)
  if nargin < 4
    method = 'linear';
  end
  % Unroll the longitude into X and Y vector coordinates
  x_coords = cosd(data_deg);
  y_coords = sind(data_deg);
  
  % Perform standard linear interpolation on both components
  xq_coords = interp1(x, x_coords, xq, method);
  yq_coords = interp1(x, y_coords, xq, method);
  
  % Reconstruct the angles from the interpolated components
  % and force the output back into the -180 to 180 degree range
  data_degq = atan2d(yq_coords, xq_coords);
end
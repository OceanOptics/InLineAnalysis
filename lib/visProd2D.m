function fh = visProd2D(x, dt, data, smooth, figid)
% INPUT:
%   x <1xM double> wavelength (lambda) or scattering angle (theta)
%   dt <Nx1 datenum> date & time
%   data <NxM double> attenuation or absorption spectrum
%   smooth <boolean> smooth data
%
% Require: colorval

% Unselec all NaN
sel = any(~isnan(data),2);

% Smooth
if nargin < 4;  smooth = false; end
if smooth
  Y = filtfilt(ones(10, 1), 10, data(sel,:)); 
else
  Y = data(sel,:);
end

% Figure id
if nargin < 5; figid = 61; end

fh = fig(figid);
if size(dt, 1) > 2
  % Save current defaultAxesColorOrder
  defaultAxesColorOrder = get(groot, 'defaultAxesColorOrder');
  set(groot,'defaultAxesColorOrder',colorval(dt(sel), parula));
  plot(x, Y)
else
  h = plot(x, Y, 'LineWidth',  3);
  col = reshape(spectrumRGB(x), max(size(x)),  3);
  col = uint8(col'*255); % need a 4xN uint8 array
  col(4,:) = 255; % last column is transparency
  pause(0.001)
  set(h.Edge,'ColorBinding','interpolated','ColorData',col)
end
pause(0.001)
xlim([min(x) max(x)])
if size(dt, 1) > 1
  colormap(parula);
  cb = colorbar();
  caxis([min(dt(sel)) max(dt(sel))]);
  datetick(cb, 'y', 'HH:MM mmm dd', 'keepticks');
  xlabel('Wavelength (nm)');
end
% ylabel('a_p (m^{-1})');

if size(dt, 1) > 1
  % Reset defaultAxesColorOrder
  set(groot,'defaultAxesColorOrder', defaultAxesColorOrder)
end
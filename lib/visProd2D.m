function fh = visProd2D(x, dt, data, smooth)
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
if smooth; Y = filtfilt(ones(10, 1), 10, data(sel,:)); 
else Y = data(sel,:); end

fh = fig(61);
if size(dt, 1) > 1
  % Save current defaultAxesColorOrder
  defaultAxesColorOrder = get(groot, 'defaultAxesColorOrder');
  set(groot,'defaultAxesColorOrder',colorval(dt(sel), parula));
end
plot(x, Y);
if size(dt, 1) > 1
  colormap(parula);
  cb = colorbar();
  caxis([min(dt(sel)) max(dt(sel))]);
  datetick(cb, 'y', 'HH:MM mmm dd', 'keepticks');
  xlabel('Wavelength (nm)');
end
% ylabel('a_p (m^{-1})');

% Reset defaultAxesColorOrder
set(groot,'defaultAxesColorOrder', defaultAxesColorOrder)

end
function fh = visProd3D(x, dt, data, smooth, color, autorotate)
%
% INPUT:
%   x <1xM double> wavelength (lambda) or scattering angle (theta)
%   dt <Nx1 datenum> date & time
%   data <NxM double> attenuation or absorption spectrum
%   smooth <boolean> smooth data
%   color <'Wavelength'|'Intensity'> Color data
%   autorotate <boolean> animate figure (will stop execution of program)
%
% EXAMPLE:
%   visProd3D(lambda.ref, ACS.p.dt, ACS.p.ap, false); zlabel('a_p (m^{-1})');

% Unselect NaN
sel = any(~isnan(data),2);
% Smooth
if nargin < 4;  smooth = False; end
if smooth; Z = filtfilt(ones(10, 1), 10, data(sel,:)); 
else Z = data(sel,:); end
% Color
if nargin < 5; color='Wavelength'; end
switch color
  case 'Wavelength'
    C = spectrumRGB(x) .* ones(size(dt(sel),1),1);
  case 'Intensity'
    C = Z;
  case 'Log'
    C = log10(Z);
  otherwise
    error('Color mode not supported');
end
% Rotate figure
if nargin < 6; autorotate=false; end;

 %% 3D Mesh Plot
fh = fig(72);
%   waterfall(x, dt, Z, C);
mesh(x, dt(sel), Z, C);

datetick('y');
view(20,20);

set(datacursormode(fh),'UpdateFcn',@data_cursor_display_date_y);

% xlabel('Wavelength (nm)');
% zlabel('a_p (m^{-1})');


%% Add light and shaddow on graph (looks nicer)
light

%% Rotate the figure
if autorotate
  while true
    for i=1:360
      view(i+20,20); 
      drawnow();
    end
  end
end

end
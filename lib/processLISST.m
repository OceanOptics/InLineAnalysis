function [p, diameter, bin_size] = processLISST(param, tot, filt)
% PROCESSLISST process LISST from flow through system with both
%   filtered and total periods
%
% OUTPUT
%     p.betap   % <Nx32 double> counts
%     p.beamC   % <Nx1 double> 1/m
%     p.VD      % <Nx32 double> volume distribution (\muL/L)
%     p.PSD     % <Nx32 double> numbers / mL / \mum
%     p.VSD     % <Nx32 double> numbers x 10^-6  / \mum
%


% Interpolate filtered on Total
filt_interp = table(tot.dt, 'VariableNames', {'dt'});
filt_interp.beta = interp1(filt.dt, filt.beta, filt_interp.dt, 'linear', 'extrap');
filt_interp.laser_reference = interp1(filt.dt, filt.laser_reference, filt_interp.dt, 'linear', 'extrap');
filt_interp.laser_power = interp1(filt.dt, filt.laser_power, filt_interp.dt, 'linear', 'extrap');

% Estimate Particulate Scattering
p = table(tot.dt, 'VariableNames', {'dt'});
p.betap = tot.beta/10 - filt_interp.beta/10 .* tot.laser_reference ./ filt_interp.laser_reference;
% p.laser_power = tot.laser_power - filt_interp.laser_power .* tot.laser_reference ./ filt_interp.laser_reference;

% Calibrate Particulate Scattering based on detector area
p.betap = (ones(size(p,1),1) * param.dcal) .* p.betap;
% interp on scat_angle 30

% Compute Beam C
% Rings from file (reference for instrument)
r = param.zsc(33)/param.zsc(36);
tau = tot.laser_power ./ r ./ tot.laser_reference;
p.beamc = -log(tau)/0.05;
flref=param.zsc(36);

% Compute Volume Distribution using spherical|non-spherical inversion model
% Spherical parameter for inversion
%    0 -> spherical
%   ~0 -> non-spherical
% non_spherical = 0;
p.VD = NaN(size(p.betap));
for i=1:size(p,1)
  p.VD(i,:) = invert_2014b(p.betap(i,:),2,0,param.non_spherical,0,0,0) / param.vcc .* flref./(tot.laser_reference(i).*ones(1,32)); % unit of microl/l, or 10^-6 micron^3/micron^3
end

% Get bin size specific to inversion model and instrument type
% Diameter in microns
%   Constant specific to type B LISST and spherical inversion
% ds=1.25*1.18.^(0:1:32);
%   Constant specific to type B LISST and non-spherical inversion
% ds=1*1.18.^(0:1:32);
% Compute bin width
bin_size = param.ds(2:end) - param.ds(1:end-1);
bs = bin_size .* ones(size(p,1),1);
% Compute diameters
diameter = sqrt(param.ds(1:end-1).*param.ds(2:end));
d = diameter .* ones(size(p,1),1);

% Convert to number distribution
p.PSD = (p.VD ./ (pi*d.^3/6) ./ bs)  * 10^6; % units: # / \mum^3 / \mum
p.VSD = (p.VD ./ bs) * 10^6;                 % units: # / \mum 

end
function [p, bin_size] = processLISST(param, tot, filt, fth, fth_constants)
% PROCESSLISST process LISST from flow through system with both
%   filtered and total periods
%
% OUTPUT
%     p.betap   % <Nx32 double> particulate VSF (counts)
%     p.cp      % <Nx1 double>  particulate beam attenuation (1/m)
%     p.VD      % <Nx32 double> volume distribution (\muL/L)
%     p.PSD     % <Nx32 double> numbers / mL / m
%     p.VSD     % <Nx32 double> numbers x 10^-6  / m
%

% Instrument constants
% Path Length (m)
path_length = 0.05;
% Fraction of circle covered by detectors
phi = 1/6;

if ~exist('fth', 'var')
  % check FTH data
  if ~exist('fth_constants', 'var')
    % Assume most recent FlowControl software
    SWITCH_FILTERED = 1;
    SWITCH_TOTAL = 0;
  else
    SWITCH_FILTERED = fth_constants.SWITCH_FILTERED;
    SWITCH_TOTAL = fth_constants.SWITCH_TOTAL;
  end

  % interpolate fth.swt onto binned data to fill missing flow data
  fth_interp = table([tot.dt; fth.dt; filt.dt], 'VariableNames', {'dt'});
  [~,b] = sort(fth_interp.dt); % sort dates
  fth_interp.dt = fth_interp.dt(b,:);
  fth_interp.swt = interp1(fth.dt, fth.swt, fth_interp.dt, 'previous');%, 'linear', 'extrap');
  fth_interp.swt = fth_interp.swt > 0;
  % Find switch events from total to filtered
  sel_start = find(fth_interp.swt(1:end-1) == SWITCH_TOTAL & fth_interp.swt(2:end) == SWITCH_FILTERED);
  % Find switch events from filtered to total
  sel_end = find(fth_interp.swt(1:end-1) == SWITCH_FILTERED & fth_interp.swt(2:end) == SWITCH_TOTAL);
  % Verify selections of filtered period
  if sel_start(1) > sel_end(1); sel_end(1) = []; end
  if sel_start(end) > sel_end(end); sel_end(end+1) = size(fth_interp.swt,1); end
  if size(sel_start,1) ~= size(sel_end,1); error('Inconsistent fth data'); end

  % Compute filtered period median
  filt_avg = table((fth_interp.dt(sel_start) + fth_interp.dt(sel_end)) ./ 2, 'VariableNames', {'dt'});
  for i=1:size(sel_start, 1)
    sel_filt = fth_interp.dt(sel_start(i)) <= filt.dt & filt.dt <= fth_interp.dt(sel_end(i));
    foo = filt(sel_filt,:);
    % compute 5 percentile for each filter event
    filt_avg.beta(i,:) = prctile(foo.beta, 5, 1);
    filt_avg.beta_avg_sd(i,:) = prctile(foo.beta_avg_sd, 5, 1);
    filt_avg.laser_reference(i,:) = prctile(foo.laser_reference, 5, 1);
    filt_avg.laser_transmission(i,:) = prctile(foo.laser_transmission, 5, 1);
  end
  filt_avg(all(isnan(filt_avg.beta), 2), :) = [];
else
  filt_avg = filt;
end

% Interpolate filtered on Total
filt_interp = table(tot.dt, 'VariableNames', {'dt'});
filt_interp.beta = interp1(filt_avg.dt, filt_avg.beta, filt_interp.dt, 'linear', 'extrap');
filt_interp.beta_avg_sd = interp1(filt_avg.dt, filt_avg.beta_avg_sd, filt_interp.dt, 'linear', 'extrap');
filt_interp.laser_reference = interp1(filt_avg.dt, filt_avg.laser_reference, filt_interp.dt, 'linear', 'extrap');
filt_interp.laser_transmission = interp1(filt_avg.dt, filt_avg.laser_transmission, filt_interp.dt, 'linear', 'extrap');

if exist('visFlag', 'file') && exist('fth', 'var')
  visFlag(tot, filt_interp, [], [], filt_avg, [], 'beta', round(size(tot.beta, 2)/2), [], fth);
  title('Check filter event interpolation')
  legend('Total', 'Filtered interpolated', 'Filtered median', 'Flow rate')
elseif exist('visFlag', 'file')
  visFlag(tot, filt_interp, [], [], filt_avg, [], 'beta', round(size(tot.beta, 2)/2), [], []);
  title('Check filter event interpolation')
  legend('Total', 'Filtered interpolated', 'Filtered median')
end
  

% Set output table
p = table(tot.dt, 'VariableNames', {'dt'});

% Code use in Boss et al. 2018
% % % % Capture change in laser reference between filtered and total periods
% % % r = param.zsc(33)/param.zsc(36); % Rings from file (reference for instrument)
% % % tau_tot = tot.laser_transmission ./ r ./ tot.laser_reference;
% % % tau_filt = filt_interp.laser_transmission ./ r ./ filt_interp.laser_reference;
% % % 
% % % % Estimate particulate scattering in counts
% % % % p.betap = tot.beta/10 - filt_interp.beta/10 .* tot.laser_reference ./ filt_interp.laser_reference; % Does not take into account variation in laser power
% % % p.betap = tot.beta / 10 ./ tau_tot - filt_interp.beta / 10 ./ tau_filt;
% % % % p.laser_transmission = tot.laser_transmission - filt_interp.laser_transmission .* tot.laser_reference ./ filt_interp.laser_reference;
% % % p.betap_sd = sqrt(tot.beta_avg_sd.^2 + filt_interp.beta_avg_sd.^2) / 10;
% % % 
% % % % Correct particulate scattering for detector area and responsivness
% % % p.betap = (ones(size(p,1),1) * param.dcal) .* p.betap;
% % % p.betap_sd = (ones(size(p,1),1) * param.dcal) .* p.betap_sd;
% % % 
% % % % Correct for laser power in TSW relative to FSW
% % % % p.betap = p.betap ./ tot.laser_reference ./ filt_interp.laser_reference;
% % % 
% % % % Compute Beam C
% % % beamc_tot = -log(tau_tot) / path_length;
% % % beamc_filt = -log(tau_filt) / path_length;
% % % p.cp = beamc_tot - beamc_filt;
% % % 
% % % % Compute Volume Distribution using spherical|non-spherical inversion model
% % % % Spherical parameter for inversion
% % % %    0 -> spherical
% % % %   ~0 -> non-spherical
% % % zsc_laser_reference = param.zsc(36);
% % % p.invert = NaN(size(p.betap));
% % % p.VD = NaN(size(p.betap));
% % % for i=1:size(p,1)
% % %   % Use factory zsc taken at the same time as the VCC for the inversion.
% % %   % Need betap in counts (NOT Scientific Units)
% % %   p.VD(i,:) = invert_2014b(p.betap(i,:),2,0,param.non_spherical,0,0,0) ./ param.vcc ...
% % %               .* zsc_laser_reference ./ (tot.laser_reference(i).*ones(1,32)); % unit of microl/l, or 10^-6 micron^3/micron^3
% % %   p.VD_sd(i,:) = invert_2014b(p.betap_sd(i,:),2,0,param.non_spherical,0,0,0) ./ param.vcc ...
% % %               .* zsc_laser_reference ./ (tot.laser_reference(i).*ones(1,32)); % standard deviation (same units)
% % % end

% Code updated on Jan 16, 2019
%   dcal is independent of zsc so removed zsc from code and normalize directly to FSW
% Capture change in laser reference
r = filt_interp.laser_transmission ./ filt_interp.laser_reference;
tau = tot.laser_transmission ./ tot.laser_reference ./ r;
% Compute Beam C
p.cp = -log(tau) / path_length;
% Estimate particulate scattering in counts
p.betap = tot.beta / 10 ./ tau - filt_interp.beta / 10;
p.betap_sd = sqrt((tot.beta_avg_sd/10./tau).^2 + (filt_interp.beta_avg_sd/10).^2);
% Correct for laser power in TSW relative to FSW
p.betap = p.betap .* tot.laser_reference ./ filt_interp.laser_reference;
% Correct particulate scattering for detector area and responsivness
p.betap = (ones(size(p,1),1) * param.dcal) .* p.betap;
p.betap_sd = (ones(size(p,1),1) * param.dcal) .* p.betap_sd;

% DCal correction for ring 30 in low gain
% Interpolate betap at ring 30
betap_interp_30 = interp1(param.theta([1:29 31:32])', p.betap(:,[1:29 31:32])', param.theta(30))';
% Select all betap from ring 30 that are more than 15% away from interpolated value
sel_LOW_GAIN = (betap_interp_30 - p.betap(:,30)) ./ ((betap_interp_30 + p.betap(:,30)) / 2) > 0.15;
% Replace bad value by interpolated values as can't do inversion without values
p.betap(sel_LOW_GAIN,30) = betap_interp_30(sel_LOW_GAIN);

% Compute Volume Distribution using spherical|non-spherical inversion model
% Spherical parameter for inversion
%    0 -> spherical
%   ~0 -> non-spherical
p.invert = NaN(size(p.betap));
p.VD = NaN(size(p.betap));
for i=1:size(p,1)
  % Use factory zsc taken at the same time as the VCC for the inversion.
  % Need betap in counts (NOT Scientific Units)
  p.VD(i,:) = invert_2014b(p.betap(i,:),2,0,param.non_spherical,0,0,0) ./ param.vcc ...
              .* filt_interp.laser_reference(i) ./ (tot.laser_reference(i).*ones(1,32)); % unit of microl/l, or 10^-6 micron^3/micron^3
  p.VD_sd(i,:) = invert_2014b(p.betap_sd(i,:),2,0,param.non_spherical,0,0,0) ./ param.vcc ...
              .* filt_interp.laser_reference(i) ./ (tot.laser_reference(i).*ones(1,32)); % standard deviation (same units)
end

% Convert betap from counts to scientific units (1/m)
p.betap =  p.betap ./ (ones(height(p),1) .* (pi * phi * path_length * (param.ds(2:end).^2 - param.ds(1:end-1).^2)));

% Get bin size specific to inversion model and instrument type
% Diameter in microns
%   Constant specific to type B LISST and spherical inversion
% ds=1.25*1.18.^(0:1:32);
%   Constant specific to type B LISST and non-spherical inversion
% ds=1*1.18.^(0:1:32);
% Compute bin width
bin_size = param.ds(2:end) - param.ds(1:end-1);
bs = bin_size .* ones(size(p,1),1);
% Compute diameters (deprecated now taken as argument)
% diameter = sqrt(param.ds(1:end-1).*param.ds(2:end));
% d = diameter .* ones(size(p,1),1);
% Load diameters
d = param.diameters .* ones(size(p,1),1);

% Convert to number distribution
p.PSD = (p.VD ./ (pi*d.^3/6) ./ bs) * 10^6; % units: # / \mum^3 / m
p.VSD = (p.VD ./ bs) * 10^6;                % units: ppm / m 

% Perform simple QC
p(any(p.PSD < 0,2),:) = [];

end
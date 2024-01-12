function [p, pbin_size, g, gbin_size] = processLISST(param, tot, filt, di, fth, fth_constants, di_method, days2run)
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%% DIW NOT IMPLEMENTED BUT READY TO BE %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Instrument constants
% Path Length (m)
path_length = 0.05;
% Fraction of circle covered by detectors
phi = 1/6;

% detect if laser_power or laser_transmission
totvar = tot.Properties.VariableNames;
LASER_POWER_VAR = totvar{contains(totvar, {'laser_power', 'laser_transmission'}) & ...
  ~contains(totvar, {'_avg_sd', '_avg_n'})};

if ~isempty(filt)
  if exist('fth', 'var')
    % check FTH data
    if ~exist('fth_constants', 'var')
      % Assume most recent FlowControl software
      SWITCH_FILTERED = 1;
      SWITCH_TOTAL = 0;
    else
      SWITCH_FILTERED = fth_constants.SWITCH_FILTERED;
      SWITCH_TOTAL = fth_constants.SWITCH_TOTAL;
    end
    % remove duplicates
    [~, L, ~] = unique(fth.qc.tsw.dt,'first');
    indexToDump = not(ismember(1:numel(fth.qc.tsw.dt), L));
    if sum(indexToDump) > 0
      fprintf('Warning: %i identical dates in FLOW data => deleted\n', sum(indexToDump))
      fth.qc.tsw(indexToDump, :) = [];
    end
    % remove duplicates
    [~, L, ~] = unique(tot.dt,'first');
    indexToDump = not(ismember(1:numel(tot.dt), L));
    if sum(indexToDump) > 0
      fprintf('Warning: %i identical dates in total data => deleted\n', sum(indexToDump))
      tot(indexToDump, :) = [];
    end
    % remove duplicates
    [~, L, ~] = unique(filt.dt,'first');
    indexToDump = not(ismember(1:numel(filt.dt), L));
    if sum(indexToDump) > 0
      fprintf('Warning: %i identical dates in filtered data => deleted\n', sum(indexToDump))
      filt(indexToDump, :) = [];
    end

    % interpolate fth.qc.tsw.swt onto binned data to fill missing flow data
    fth_interp = table([tot.dt; fth.qc.tsw.dt; filt.dt], 'VariableNames', {'dt'});
    [~,b] = sort(fth_interp.dt); % sort dates
    fth_interp.dt = fth_interp.dt(b,:);
    fth_interp.swt = interp1(fth.qc.tsw.dt, fth.qc.tsw.swt, fth_interp.dt, 'previous');%, 'linear', 'extrap');
    fth_interp.swt = fth_interp.swt > 0;
    % Find switch events from total to filtered
    sel_start = find(fth_interp.swt(1:end-1) == SWITCH_TOTAL & fth_interp.swt(2:end) == SWITCH_FILTERED);
    % Find switch events from filtered to total
    sel_end = find(fth_interp.swt(1:end-1) == SWITCH_FILTERED & fth_interp.swt(2:end) == SWITCH_TOTAL);
    % Verify selections of filtered period
    if sel_start(1) > sel_end(1); sel_end(1) = []; end
    if sel_start(end) > sel_end(end); sel_end(end+1) = size(fth_interp.swt,1); end
    if size(sel_start,1) ~= size(sel_end,1); error('Inconsistent fth data'); end
    % interpolate filtered event
  
    % Compute filtered period median
    filt_avg = table(NaN(size(sel_start)), 'VariableNames', {'dt'});
    filt_avg.beta = NaN(size(filt_avg,1), size(filt.beta, 2));
    filt_avg.beta_avg_sd = NaN(size(filt_avg,1), size(filt.beta_avg_sd, 2));
    filt_avg.laser_reference = NaN(size(filt_avg,1), 1);
    filt_avg.(LASER_POWER_VAR) = NaN(size(filt_avg,1), 1);
    for i=1:size(sel_start, 1)
      sel_filt = fth_interp.dt(sel_start(i)) <= filt.dt & filt.dt <= fth_interp.dt(sel_end(i));
      foo = filt(sel_filt,:);
      if sum(sel_filt) == 1
        filt_avg.dt(i) = foo.dt;
        filt_avg.beta(i,:) = foo.beta;
        filt_avg.beta_avg_sd(i,:) = foo.beta_avg_sd;
        filt_avg.laser_reference(i,:) = foo.laser_reference;
        filt_avg.(LASER_POWER_VAR)(i,:) = foo.(LASER_POWER_VAR);
      else
        perc25 = foo.beta > prctile(foo.beta, 25, 1);
        foo.beta_avg_sd(perc25) = NaN;
        foo.beta(perc25) = NaN;
        foo.laser_reference(perc25) = NaN;
        foo.(LASER_POWER_VAR)(perc25) = NaN;
        % compute average of all values smaller than 25th percentile for each filter event
        filt_avg.dt(i) = mean(foo.dt(any(~perc25, 2)), 'omitnan');
        filt_avg.beta(i,:) = mean(foo.beta, 1, 'omitnan');
        filt_avg.beta_avg_sd(i,:) = mean(foo.beta_avg_sd, 1, 'omitnan');
        filt_avg.laser_reference(i) = mean(foo.laser_reference, 1, 'omitnan');
        filt_avg.(LASER_POWER_VAR)(i) = mean(foo.(LASER_POWER_VAR), 1, 'omitnan');
      end
    end
    filt_avg(all(isnan(filt_avg.beta), 2), :) = [];
  else
    filt_avg = filt;
  end

  % Interpolate filtered on Total
  filt_interp = table(tot.dt, 'VariableNames', {'dt'});
  filt_interp.beta = interp1(filt_avg.dt, filt_avg.beta, filt_interp.dt);
  filt_interp.beta_avg_sd = interp1(filt_avg.dt, filt_avg.beta_avg_sd, filt_interp.dt);
  filt_interp.laser_reference = interp1(filt_avg.dt, filt_avg.laser_reference, filt_interp.dt);
  filt_interp.(LASER_POWER_VAR) = interp1(filt_avg.dt, filt_avg.(LASER_POWER_VAR), filt_interp.dt);
  
  % id only day to run in all tables to plot
  filt_interp_id = filt_interp.dt >= min(days2run) & filt_interp.dt < max(days2run)+1;
  tot_id = tot.dt >= min(days2run) & tot.dt < max(days2run)+1;
  filt_avg_id = filt_avg.dt >= min(days2run) & filt_avg.dt < max(days2run)+1;
  
  % plot
  if exist('visFlag', 'file') && exist('fth', 'var')
    fh = visFlag([], filt_interp(filt_interp_id, :), tot(tot_id, :), [], filt_avg(filt_avg_id, :), [], 'beta', round(size(tot.beta, 2)/2), ...
      [], fth.qc.tsw, fth.view.spd_variable);
    title('Check filter event interpolation, press q to continue', 'FontSize', 14)
    legend('Filtered interpolated', 'Total', 'Filtered median', 'Flow rate',...
      'AutoUpdate','off', 'FontSize', 12)
    guiSelectOnTimeSeries(fh);
  elseif exist('visFlag', 'file')
    fh = visFlag([], filt_interp(filt_interp_id, :), tot(tot_id, :), [], filt_avg(filt_avg_id, :), [], 'beta', round(size(tot.beta, 2)/2), [], []);
    title('Check filter event interpolation, press q to continue', 'FontSize', 14)
    legend('Filtered interpolated', 'Total', 'Filtered median', 'Flow rate',...
      'AutoUpdate','off', 'FontSize', 12)
    guiSelectOnTimeSeries(fh);
  end
else
  warning('Filtered measurements not provided files, instrument zscat was used.')
  filt_interp = [];
end

% % Set output table
% p = table(tot.dt, 'VariableNames', {'dt'});

% Code use in Boss et al. 2018
% % % % Capture change in laser reference between filtered and total periods
% % % r = param.zsc(33)/param.zsc(36); % Rings from file (reference for instrument)
% % % tau_tot = tot.(LASER_POWER_VAR) ./ r ./ tot.laser_reference;
% % % tau_filt = filt_interp.(LASER_POWER_VAR) ./ r ./ filt_interp.laser_reference;
% % % 
% % % % Estimate particulate scattering in counts
% % % % p.betap = tot.beta/10 - filt_interp.beta/10 .* tot.laser_reference ./ filt_interp.laser_reference; % Does not take into account variation in laser power
% % % p.betap = tot.beta / 10 ./ tau_tot - filt_interp.beta / 10 ./ tau_filt;
% % % % p.(LASER_POWER_VAR) = tot.(LASER_POWER_VAR) - filt_interp.(LASER_POWER_VAR) .* tot.laser_reference ./ filt_interp.laser_reference;
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

[p, pbin_size] = compute_product_LISST100X(tot, filt_interp, param, phi, path_length, LASER_POWER_VAR);

% % Code updated on Jan 16, 2019
% %   dcal is independent of zsc so removed zsc from code and normalize directly to FSW
% % Capture change in laser reference
% r = filt_interp.(LASER_POWER_VAR) ./ filt_interp.laser_reference;
% tau = tot.(LASER_POWER_VAR) ./ tot.laser_reference ./ r;
% % Compute Beam C
% p.cp = -log(tau) / path_length;
% % Estimate particulate scattering in counts
% p.betap = tot.beta / 10 ./ tau - filt_interp.beta / 10;
% p.betap_sd = sqrt((tot.beta_avg_sd/10./tau).^2 + (filt_interp.beta_avg_sd/10).^2);
% % Correct for laser power in TSW relative to FSW
% p.betap = p.betap .* tot.laser_reference ./ filt_interp.laser_reference;
% % Correct particulate scattering for detector area and responsivness
% p.betap = (ones(size(p,1),1) * param.dcal) .* p.betap;
% p.betap_sd = (ones(size(p,1),1) * param.dcal) .* p.betap_sd;
% 
% % DCal correction for ring 30 in low gain
% % Interpolate betap at ring 30
% betap_interp_30 = interp1(param.theta([1:29 31:32])', p.betap(:,[1:29 31:32])', param.theta(30))';
% % Select all betap from ring 30 that are more than 15% away from interpolated value
% sel_LOW_GAIN = (betap_interp_30 - p.betap(:,30)) ./ ((betap_interp_30 + p.betap(:,30)) / 2) > 0.15;
% % Replace bad value by interpolated values as can't do inversion without values
% p.betap(sel_LOW_GAIN,30) = betap_interp_30(sel_LOW_GAIN);
% 
% % Compute Volume Distribution using spherical|non-spherical inversion model
% % Spherical parameter for inversion
% %    0 -> spherical
% %   ~0 -> non-spherical
% p.invert = NaN(size(p.betap));
% p.VD = NaN(size(p.betap));
% for i=1:size(p,1)
%   % Use factory zsc taken at the same time as the VCC for the inversion.
%   % Need betap in counts (NOT Scientific Units)
%   p.VD(i,:) = invert_2014b(p.betap(i,:),2,0,param.non_spherical,0,0,0) ./ param.vcc ...
%               .* filt_interp.laser_reference(i) ./ (tot.laser_reference(i).*ones(1,32)); % unit of microl/l, or 10^-6 micron^3/micron^3
%   p.VD_sd(i,:) = invert_2014b(p.betap_sd(i,:),2,0,param.non_spherical,0,0,0) ./ param.vcc ...
%               .* filt_interp.laser_reference(i) ./ (tot.laser_reference(i).*ones(1,32)); % standard deviation (same units)
% end
% 
% % Convert betap from counts to scientific units (1/m)
% p.betap =  p.betap ./ (ones(height(p),1) .* (pi * phi * path_length * (param.ds(2:end).^2 - param.ds(1:end-1).^2)));
% 
% % Get bin size specific to inversion model and instrument type
% % Diameter in microns
% %   Constant specific to type B LISST and spherical inversion
% % ds=1.25*1.18.^(0:1:32);
% %   Constant specific to type B LISST and non-spherical inversion
% % ds=1*1.18.^(0:1:32);
% % Compute bin width
% bin_size = param.ds(2:end) - param.ds(1:end-1);
% bs = bin_size .* ones(size(p,1),1);
% % Compute diameters (deprecated now taken as argument)
% % diameter = sqrt(param.ds(1:end-1).*param.ds(2:end));
% % d = diameter .* ones(size(p,1),1);
% % Load diameters
% d = param.diameters .* ones(size(p,1),1);
% 
% % Convert to number distribution
% p.PSD = (p.VD ./ (pi*d.^3/6) ./ bs) * 10^6; % units: # / \mum^3 / m
% p.VSD = (p.VD ./ bs) * 10^6;                % units: ppm / m 
% 
% % Perform simple QC
% p(any(p.PSD < 0,2),:) = [];

if nargout > 2 && ~isempty(di)
  switch di_method
    case 'interpolate'
      % Interpolate DI on Filtered
      %     + recommende if sensor drift with time
      di_pp = table(filt_avg.dt, 'VariableNames', {'dt'});
      di_pp.beta = interp1(di.dt, di.beta, di_pp.dt);
      di_pp.beta_avg_sd = interp1(di.dt, di.beta_avg_sd, di_pp.dt);
    case 'constant'
      % Average all given DI samples
      %     + recommended if no drift are observed with sensor
      di_pp = table(NaN, 'VariableNames', {'dt'});
      % Get not NaN DI values
      di_beta_sel = di.beta(all(~isnan(di.beta),2));
      di_beta_avg_sd_sel = di.beta_avg_sd(all(~isnan(di.beta),2));
      % Select only DI values within 5th and 75th percentile
      foo = prctile(di_beta_sel,[5, 75],1);
      avg_pl(1,:) = foo(1,:); % low percentile
      avg_ph(1,:) = foo(2,:); % high percentile
      avg_sel = any(avg_pl(1,:) <= di_beta_sel & di_beta_sel <= avg_ph(1,:),2);
      % Average values
      di_pp.beta = mean(di_beta_sel(avg_sel,:));
      di_pp.beta_avg_sd = mean(di_beta_avg_sd_sel(avg_sel,:));
    otherwise
      error('Method not supported.');
  end
  [g, gbin_size] = compute_product_LISST100X(filt_avg, di_pp, param, phi, path_length, LASER_POWER_VAR);
else
  g = table();
end
end


function [prod, bin_size] = compute_product_LISST100X(tot, filt_interp, param, phi, path_length, LASER_POWER_VAR)
% Set output table
prod = table(tot.dt, 'VariableNames', {'dt'});

% Code updated on Jan 16, 2019
% dcal is independent of zsc so removed zsc from code and normalize directly to FSW
% Capture change in laser reference
r = filt_interp.(LASER_POWER_VAR) ./ filt_interp.laser_reference;
tau = tot.(LASER_POWER_VAR) ./ tot.laser_reference ./ r;
% Compute Beam C
prod.cp = -log(tau) / path_length;
% Estimate particulate scattering in counts
prod.betap = tot.beta / 10 ./ tau - filt_interp.beta / 10;
prod.betap_sd = sqrt((tot.beta_avg_sd/10./tau).^2 + (filt_interp.beta_avg_sd/10).^2);
% Correct for laser power in TSW relative to FSW or FSW relative to DIW 
prod.betap = prod.betap .* tot.laser_reference ./ filt_interp.laser_reference;
% Correct particulate scattering for detector area and responsivness
prod.betap = (ones(size(prod,1),1) * param.dcal) .* prod.betap;
prod.betap_sd = (ones(size(prod,1),1) * param.dcal) .* prod.betap_sd;

% DCal correction for ring 30 in low gain
% Interpolate betap at ring 30
betap_interp_30 = interp1(param.theta([1:29 31:32])', prod.betap(:,[1:29 31:32])', param.theta(30))';
% Select all betap from ring 30 that are more than 15% away from interpolated value
sel_LOW_GAIN = (betap_interp_30 - prod.betap(:,30)) ./ ((betap_interp_30 + prod.betap(:,30)) / 2) > 0.15;
% Replace bad value by interpolated values as can't do inversion without values
prod.betap(sel_LOW_GAIN,30) = betap_interp_30(sel_LOW_GAIN);

% Compute Volume Distribution using spherical|non-spherical inversion model
% Spherical parameter for inversion
%    0 -> spherical
%   ~0 -> non-spherical
prod.invert = NaN(size(prod.betap));
prod.VD = NaN(size(prod.betap));
for i=1:size(prod,1)
  % Use factory zsc taken at the same time as the VCC for the inversion.
  % Need betap in counts (NOT Scientific Units)
  prod.VD(i,:) = invert_2014b(prod.betap(i,:),2,0,param.non_spherical,0,0,0) ./ param.vcc ...
              .* filt_interp.laser_reference(i) ./ (tot.laser_reference(i).*ones(1,32)); % unit of microl/l, or 10^-6 micron^3/micron^3
  prod.VD_sd(i,:) = invert_2014b(prod.betap_sd(i,:),2,0,param.non_spherical,0,0,0) ./ param.vcc ...
              .* filt_interp.laser_reference(i) ./ (tot.laser_reference(i).*ones(1,32)); % standard deviation (same units)
end

% Convert betap from counts to scientific units (1/m)
prod.betap =  prod.betap ./ (ones(height(prod),1) .* (pi * phi * path_length * (param.ds(2:end).^2 - param.ds(1:end-1).^2)));

% Get bin size specific to inversion model and instrument type
% Diameter in microns
%   Constant specific to type B LISST and spherical inversion
% ds=1.25*1.18.^(0:1:32);
%   Constant specific to type B LISST and non-spherical inversion
% ds=1*1.18.^(0:1:32);
% Compute bin width
bin_size = param.ds(2:end) - param.ds(1:end-1);
bs = bin_size .* ones(size(prod,1),1);
% Compute diameters (deprecated now taken as argument)
% diameter = sqrt(param.ds(1:end-1).*param.ds(2:end));
% d = diameter .* ones(size(p,1),1);
% Load diameters
d = param.diameters .* ones(size(prod,1),1);

% Convert to number distribution
prod.PSD = (prod.VD ./ (pi*d.^3/6) ./ bs) * 10^6; % units: # / \mum^3 / m
prod.VSD = (prod.VD ./ bs) * 10^6;                % units: ppm / m 

prod.Properties.VariableUnits = {'','1/m','1/m/Sr','1/m/Sr', '','uL/L','uL/L','nb/mL/micron', 'ppm/m'};
prod.Properties.UserData = param;

% Perform simple QC
prod(any(prod.PSD < 0,2),:) = [];
end
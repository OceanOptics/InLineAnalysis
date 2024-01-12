function [p, pbin_size, g, gbin_size] = ProcessLISST200X(param, tot, filt, di, fth, fth_constants, di_method, days2run)
  % PROCESSLISST200X process LISST200X from flow through system
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
      filt_avg.RingValues = NaN(size(filt_avg,1), size(filt.RingValues, 2));
      filt_avg.RingValues_avg_sd = NaN(size(filt_avg,1), size(filt.RingValues_avg_sd, 2));
      filt_avg.RingValues_avg_n = NaN(size(filt_avg,1), size(filt.RingValues_avg_sd, 2));
      filt_avg.LaserReference = NaN(size(filt_avg,1), 1);
      filt_avg.LaserReference_avg_sd = NaN(size(filt_avg,1), 1);
      filt_avg.LaserReference_avg_n = NaN(size(filt_avg,1), 1);
      filt_avg.LaserTransmission = NaN(size(filt_avg,1), 1);
      filt_avg.LaserTransmission_avg_sd = NaN(size(filt_avg,1), 1);
      filt_avg.LaserTransmission_avg_n = NaN(size(filt_avg,1), 1);
      filt_avg.TotalVolumeConcentration = NaN(size(filt_avg,1), 1);
      filt_avg.TotalVolumeConcentration_avg_sd = NaN(size(filt_avg,1), 1);
      filt_avg.TotalVolumeConcentration_avg_n = NaN(size(filt_avg,1), 1);
      for i=1:size(sel_start, 1)
        sel_filt = fth_interp.dt(sel_start(i)) <= filt.dt & filt.dt <= fth_interp.dt(sel_end(i));
        foo = filt(sel_filt,:);
        if sum(sel_filt) == 1
          filt_avg.dt(i) = foo.dt;
          filt_avg.RingValues(i,:) = foo.RingValues;
          filt_avg.RingValues_avg_sd(i,:) = foo.RingValues_avg_sd;
          filt_avg.RingValues_avg_n(i,:) = foo.RingValues_avg_n;
          filt_avg.LaserReference(i,:) = foo.LaserReference;
          filt_avg.LaserReference_avg_sd(i,:) = foo.LaserReference_avg_sd;
          filt_avg.LaserReference_avg_n(i,:) = foo.LaserReference_avg_n;
          filt_avg.LaserTransmission(i,:) = foo.LaserTransmission;
          filt_avg.LaserTransmission_avg_sd(i,:) = foo.LaserTransmission_avg_sd;
          filt_avg.LaserTransmission_avg_n(i,:) = foo.LaserTransmission_avg_n;
          filt_avg.TotalVolumeConcentration(i,:) = foo.TotalVolumeConcentration;
          filt_avg.TotalVolumeConcentration_avg_sd(i,:) = foo.TotalVolumeConcentration_avg_sd;
          filt_avg.TotalVolumeConcentration_avg_n(i,:) = foo.TotalVolumeConcentration_avg_n;
        else
          perc25 = foo.LaserTransmission > prctile(foo.LaserTransmission, 25, 1);
          foo.RingValues_avg_sd(perc25) = NaN;
          foo.RingValues(perc25) = NaN;
          foo.LaserReference(perc25) = NaN;
          foo.LaserReference_avg_sd(perc25) = NaN;
          foo.LaserTransmission_avg_sd(perc25) = NaN;
          foo.LaserTransmission(perc25) = NaN;
          foo.TotalVolumeConcentration_avg_sd(perc25) = NaN;
          foo.TotalVolumeConcentration(perc25) = NaN;
          % compute average of all values smaller than 25th percentile for each filter event
          filt_avg.dt(i) = mean(foo.dt(any(~perc25, 2)), 'omitnan');
          filt_avg.RingValues(i,:) = mean(foo.RingValues, 1, 'omitnan');
          filt_avg.RingValues_avg_sd(i,:) = mean(foo.RingValues_avg_sd, 1, 'omitnan');
          filt_avg.LaserReference(i) = mean(foo.LaserReference, 1, 'omitnan');
          filt_avg.LaserReference_avg_sd(i) = mean(foo.LaserReference_avg_sd, 1, 'omitnan');
          filt_avg.LaserTransmission(i) = mean(foo.LaserTransmission, 1, 'omitnan');
          filt_avg.LaserTransmission_avg_sd(i) = mean(foo.LaserTransmission_avg_sd, 1, 'omitnan');
          filt_avg.TotalVolumeConcentration(i) = mean(foo.TotalVolumeConcentration, 1, 'omitnan');
          filt_avg.TotalVolumeConcentration_avg_sd(i) = mean(foo.TotalVolumeConcentration_avg_sd, 1, 'omitnan');
          filt_avg.RingValues_avg_n(i) = sum(foo.RingValues_avg_n(any(~isnan(foo.RingValues), 2)), 'omitnan');
          filt_avg.LaserReference_avg_n(i) = sum(foo.LaserReference_avg_n(any(~isnan(foo.LaserReference), 2)), 'omitnan');
          filt_avg.LaserTransmission_avg_n(i) = sum(foo.LaserTransmission_avg_n(any(~isnan(foo.LaserTransmission), 2)), 'omitnan');
          filt_avg.TotalVolumeConcentration_avg_n(i) = sum(foo.TotalVolumeConcentration_avg_n(any(~isnan(foo.TotalVolumeConcentration), 2)), 'omitnan');
        end
      end
      filt_avg(all(isnan(filt_avg.RingValues), 2), :) = [];
    else
      filt_avg = filt;
    end
    % Interpolate filtered on Total
    filt_interp = table(tot.dt, 'VariableNames', {'dt'});
    filt_interp.RingValues = interp1(filt_avg.dt, filt_avg.RingValues, filt_interp.dt);
    filt_interp.RingValues_avg_sd = interp1(filt_avg.dt, filt_avg.RingValues_avg_sd, filt_interp.dt);
    filt_interp.LaserReference = interp1(filt_avg.dt, filt_avg.LaserReference, filt_interp.dt);
    filt_interp.LaserReference_avg_sd = interp1(filt_avg.dt, filt_avg.LaserReference_avg_sd, filt_interp.dt);
    filt_interp.LaserTransmission = interp1(filt_avg.dt, filt_avg.LaserTransmission, filt_interp.dt);
    filt_interp.LaserTransmission_avg_sd = interp1(filt_avg.dt, filt_avg.LaserTransmission_avg_sd, filt_interp.dt);
    filt_interp.TotalVolumeConcentration = interp1(filt_avg.dt, filt_avg.TotalVolumeConcentration, filt_interp.dt);
    filt_interp.TotalVolumeConcentration_avg_sd = interp1(filt_avg.dt, filt_avg.TotalVolumeConcentration_avg_sd, filt_interp.dt);
  
    % id only day to run in all tables to plot
    filt_interp_id = filt_interp.dt >= min(days2run) & filt_interp.dt < max(days2run)+1;
    tot_id = tot.dt >= min(days2run) & tot.dt < max(days2run)+1;
    filt_avg_id = filt_avg.dt >= min(days2run) & filt_avg.dt < max(days2run)+1;
  
    % plot
    if exist('visFlag', 'file') && exist('fth', 'var')
      fh = visFlag([], filt_interp(filt_interp_id, :), tot(tot_id, :), [], filt_avg(filt_avg_id, :), [], 'LaserTransmission', round(size(tot.LaserTransmission, 2)/2), ...
        [], fth.qc.tsw, fth.view.spd_variable);
      title('Check filter event interpolation, press q to continue', 'FontSize', 14)
      legend('Filtered interpolated', 'Total', 'Filtered median', 'Flow rate',...
        'AutoUpdate','off', 'FontSize', 12)
      guiSelectOnTimeSeries(fh);
    elseif exist('visFlag', 'file')
      fh = visFlag([], filt_interp(filt_interp_id, :), tot(tot_id, :), [], filt_avg(filt_avg_id, :), [], 'LaserTransmission', round(size(tot.LaserTransmission, 2)/2), [], []);
      title('Check filter event interpolation, press q to continue', 'FontSize', 14)
      legend('Filtered interpolated', 'Total', 'Filtered median', 'Flow rate',...
        'AutoUpdate','off', 'FontSize', 12)
      guiSelectOnTimeSeries(fh);
    end
  else
    warning('Filtered measurements not provided files, instrument zscat was used.')
    filt_interp = [];
  end
  
  [p, pbin_size] = compute_product_LISST200X(tot, filt_interp, param);
  
  if nargout > 2 && ~isempty(di)
    switch di_method
      case 'interpolate'
        % Interpolate DI on Filtered
        %     + recommende if sensor drift with time
        di_pp = table(filt_avg.dt, 'VariableNames', {'dt'});
        di_pp.RingValues = interp1(di.dt, di.RingValues, di_pp.dt);
        di_pp.RingValues_avg_sd = interp1(di.dt, di.RingValues_avg_sd, di_pp.dt);
      case 'constant'
        % Average all given DI samples
        %     + recommended if no drift are observed with sensor
        di_pp = table(NaN, 'VariableNames', {'dt'});
        % Get not NaN DI values
        di_RingValues_sel = di.RingValues(all(~isnan(di.RingValues),2));
        di_RingValues_avg_sd_sel = di.RingValues_avg_sd(all(~isnan(di.RingValues),2));
        % Select only DI values within 5th and 75th percentile
        foo = prctile(di_RingValues_sel,[5, 75],1);
        avg_pl(1,:) = foo(1,:); % low percentile
        avg_ph(1,:) = foo(2,:); % high percentile
        avg_sel = any(avg_pl(1,:) <= di_RingValues_sel & di_RingValues_sel <= avg_ph(1,:),2);
        % Average values
        di_pp.RingValues = mean(di_RingValues_sel(avg_sel,:));
        di_pp.RingValues_avg_sd = mean(di_RingValues_avg_sd_sel(avg_sel,:));
      otherwise
        error('Method not supported.');
    end
    [g, gbin_size] = compute_product_LISST200X(filt_avg, di_pp, param);
  else
    g = table();
  end
end


function [prod, bin_size] = compute_product_LISST200X(tot, filt_interp, param)
  % if no filtered data, use instrument zscat
  if isempty(filt_interp)
    zsc = table();
    zsc.RingValues = repmat(param.zsc(1:36), size(tot, 1), 1);
    zsc.LaserTransmission = repmat(param.zsc(37), size(tot, 1), 1);
    zsc.LaserReference = repmat(param.zsc(40), size(tot, 1), 1);
  else
    zsc = filt_interp;
  end
  % Set output table
  prod = table(tot.dt, 'VariableNames', {'dt'});
  
  % rawData, param.zsc, param.fzs, param.dcal, param.config, param.housek, param.InversionType
  
  % Divide rings 1-36 by 10 to match L100 convention
  ringData = tot.RingValues ./ 10;
  param.fzs(1:36) = param.fzs(1:36) ./ 10;
  
  if isempty(filt_interp) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO
    zsc.RingValues = zsc.RingValues ./ 10;
  end
  
  %**************************************************************************
  %2: Compute optical transmission, raw scattering, corrected scattering and VSF
  %**************************************************************************
  
  LaserRatio = zsc.LaserTransmission ./ zsc.LaserReference;                       % ratio of transmitted power / laser ref
  prod.Transmission = zsc.LaserTransmission ./ LaserRatio ./ zsc.LaserReference;  % compute optical transmission, taking the eventual drift in laser power into account
  prod.Beamc = -log(prod.Transmission) ./ (param.config.effPath / 1000);
  
  scat = ringData ./ prod.Transmission;                                           % correct for attenuation
  scat = scat - zsc.RingValues .* zsc.LaserReference ./ zsc.LaserReference;       % subtract the background
  cscat = scat .* repmat(param.dcal, size(zsc.RingValues, 1), 1);                 % apply ring area file
  [angles, prod.VSF] = getVSF_L200X(cscat, param.config, zsc.LaserReference);     % calculate uncalibrated VSF
  cscat = cscat .* repmat(param.fzs(40) ./ zsc.LaserReference, 1, 36);            % normilize to factory LREF
  cscat = cscat ./ param.config.VCC;                                              % apply concentration calibration
  
  cscat(cscat < 0) = 0; % negative cscats are not possible, so set them to 0.
  
  if any(strcmp(param.inversion, {'Spherical', 'spherical'}))
    InvType = 0;
  elseif any(strcmp(param.inversion, {'Irregular', 'irregular'}))
    InvType = 1;
  else
    error('Unknown LISST-200X Inversion Type')
  end
  
  [prod.PSD, bin_size] = LISST_ComputePSD(cscat, InvType, 1, 0);
  
  prod = addprop(prod, {'SizeBins', 'VSFAngles'}, {'table', 'table'});
    
  prod.Properties.CustomProperties.SizeBins = bin_size;
  prod.Properties.CustomProperties.VSFAngles = angles;
  
  %**************************************************************************
  % TODO 3: Apply corrections to auxiliary data
  %**************************************************************************


end


function [angles,vsf] = getVSF_L200X(cscat, config, Lref)

  [rows,~] = size(cscat);
  
  % calculate angles in water in radians (120mm focal length)
  rho = 1.18;
  theta0air = (0.1/120);
  edge_angles = theta0air*rho.^(0:36);
  edge_angles = asin(sin(edge_angles)/1.33); % convert air angles to in water angles
  
  % find solid angle
  dOmega=cos(edge_angles(1:36))-cos(edge_angles(2:37));  
  dOmega=dOmega*2*pi/6; % factor 6 takes care of rings covering only 1/6th circle
  
  % calculate detector center angles in degrees
  angles = (180./pi) .* (sqrt(edge_angles(1:36).*edge_angles(2:37)));
  
  % compute light on rings
  Watt_per_count_on_rings = 1.9e-10; % assumed the same for all detectors
  light_on_rings = cscat.*Watt_per_count_on_rings;
  
  % calculate incident laser power from LREF
  Watt_per_count_laser_ref = 1;                 % THIS IS UNKNOWN! 
  laser_incident_power = Lref*Watt_per_count_laser_ref; 
  
  % compute VSF (no calibration for Lref, so magnitude is not correct)
  vsf = light_on_rings ./ repmat(dOmega,rows,1) ./ repmat(laser_incident_power,1,36) ./ (config.effPath/1000);
end
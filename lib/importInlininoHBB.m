function [ data, lambda ] = importInlininoHBB( filename, hbb_cal, hbb_tcal, verbose )
  % Import HyperBB data logged with Inlinino and apply SEQUOIA's processing code
  % using hbb_cal and hbb_tcal
  % Authors: SEQUOIA SCIENTIFIC
  % Modified: Guillaume Bourdin
  % Date : May 2021
  %
  %%
  % verbose = true;
  
  if nargin < 4; verbose = false; end
  if verbose
    foo = strsplit(filename, '/');
    fprintf('Importing %s ... ', foo{end});
  end
  
  p.RemoveMultiplePmtGains = true;
  try
    % Open file
    fid=fopen(filename);
    if fid==-1
      error('Unable to open file: %s', filename);
    end
  
    % Get header
    hd = strip(strsplit(strrep(fgetl(fid), ',,', ', ,'), ','));
    hd(strcmp(hd, 'time')) = {'dt'};
    % get units skipping empty lines (bug in Inlinino)
    unit = fgetl(fid);
    while isempty(unit)
      unit = fgetl(fid);
    end
    unit = strip(strsplit(strrep(unit, ',,', ', ,'), ','));
  
    % Set parser & read data
    parser = ['%s%f%f%s%s' repmat('%f', 1, size(hd(6:end),2))];
    t = textscan(fid, parser, 'delimiter', ',');
    % Close file
    fclose(fid);
  
    % Build table
    dat = table();
    for i = 1:size(hd, 2)
      if strcmp(hd{i}, 'dt')
        if contains(t{1}(1), '/')
          dat.(hd{i}) = datetime(t{i}, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
          % dat.(hd{i}) = datenum(t{i}, 'yyyy/mm/dd HH:MM:SS.FFF');
        else
          dat.(hd{i}) = datetime(cellfun(@(x) [dt_ref x], t{1}, 'UniformOutput', false), 'InputFormat', 'yyyyMMddHH:mm:ss.SSS');
          % dat.(hd{i}) = datenum(cellfun(@(x) [dt_ref x], t{1}, 'UniformOutput', false), 'yyyymmddHH:MM:SS.FFF');
        end
      else
        dat.(hd{i}) = t{i};
      end
    end
    if any(strcmp(dat.Properties.VariableNames, 'Date'))
      dat.Date = [];
    end
    if any(strcmp(dat.Properties.VariableNames, 'Time'))
      dat.Time = [];
    end
  
  catch
    warning('Error loading %s file\n', filename)
    % Open file
    fid=fopen(filename);
    if fid==-1
      error('Unable to open file: %s', filename);
    end
  
    % Get header
    hd = strip(strsplit(fgetl(fid), ','));
    hd(strcmp(hd, 'time')) = {'dt'};
    % get units skipping empty lines (bug in Inlinino)
    unit = fgetl(fid);
    while isempty(unit)
        unit = fgetl(fid);
    end
    unit = strip(strsplit(unit, ','));
    
    % get file line by line skipping end of scan count lines (instrument end of scanning cycle)
    tline = fgetl(fid);
    
    fprintf('Checking %s line by line ...\n', filename)
    t = [];
    i = 0;
    while all(~contains(tline, 'Date, TIME, POS, WL, PWR, GAIN, SIG, AVG, STD,') & ~strcmp(tline, '-1'))
      t = [t; strrep(strsplit(tline, ','), ' ', '')];
      tline = fgetl(fid);
      if isnumeric(tline)
        tline = num2str(tline);
      end
      i = i + 1;
    end
    if i >= round(i, -3)
      fprintf('Warning: %s lines checked ...\n', num2str(i))
    end
    % Read data
    tend = textscan(fid, parser, 'delimiter', ',');
    % Close file
    fclose(fid);
    
    % Build table
    if ~isempty(t)
      for i = 1:size(tend, 2)
        if isnumeric(tend{1, i})
          tend(1, i) = {[str2double(t(:,i)); tend{1, i}]};
        else
          tend(1, i) = {[t(:,i); tend{1, i}]};
        end
      end
    end
  
    % Build table
    dat = table();
    for i = 1:size(hd, 2)
      if strcmp(hd{i}, 'dt')
        if contains(tend{1}(1), '/')
          dat.(hd{i}) = datetime(tend{i}, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
          % dat.(hd{i}) = datenum(tend{i}, 'yyyy/mm/dd HH:MM:SS.FFF');
        else
          dat.(hd{i}) = datetime(cellfun(@(x) [dt_ref x], tend{1}, 'UniformOutput', false), 'InputFormat', 'yyyyMMddHH:mm:ss.SSS');
          % dat.(hd{i}) = datenum(cellfun(@(x) [dt_ref x], tend{1}, 'UniformOutput', false), 'yyyymmddHH:MM:SS.FFF');
        end
      else
        dat.(hd{i}) = tend{i};
      end
    end
    if any(strcmp(dat.Properties.VariableNames, 'Date'))
      dat.Date = [];
    end
    if any(strcmp(dat.Properties.VariableNames, 'Time'))
      dat.Time = [];
    end
  end
    
  dat(isnat(dat.dt), :) = [];
  % dat(isnan(dat.dt), :) = [];
  % dat.Properties.VariableUnits = unit(~contains(hd, {'Date', 'Time'}));
  
  mg = 0;
  for ks = unique(dat.ScanIdx')
    scanSel = dat.ScanIdx==ks;
    if length(unique(dat.PmtGain(scanSel)))~=1
      if p.RemoveMultiplePmtGains
        dat(scanSel,:) = [];
      end
      mg = mg + 1;
    end
  end
  if p.RemoveMultiplePmtGains && verbose && mg > 0
    warning(sprintf('  %i scans (%.1f%%) with multiple PmtGain: deleted', mg, mg/size(unique(dat.ScanIdx),1)*100))
  elseif verbose && mg > 0
    warning(sprintf('  %i scans (%.1f%%) with multiple PmtGain', mg, mg/size(unique(dat.ScanIdx),1)*100))
  end
  
  satLevel = 4000;
  dat.SigOn3(dat.SigOn3>satLevel) = NaN;
  dat.SigOn2(dat.SigOn2>satLevel) = NaN;
  dat.SigOn1(dat.SigOn1>satLevel) = NaN;
  dat.SigOff3(dat.SigOff3>satLevel) = NaN;
  dat.SigOff2(dat.SigOff2>satLevel) = NaN;
  dat.SigOff1(dat.SigOff1>satLevel) = NaN;
  % dat.SigOn3(dat.SigOn3>satLevel | dat.SigOn3<0) = NaN;
  % dat.SigOn2(dat.SigOn2>satLevel | dat.SigOn2<0) = NaN;
  % dat.SigOn1(dat.SigOn1>satLevel | dat.SigOn1<0) = NaN;
  
  % figure, hold on, plot(dat.SigOn3,'.-')    
  % plot(dat.SigOn1,'.-')
  % plot(dat.SigOn1,'.-')
  
  % Calculate the net high and low gain, net ref
  dat = addvars(dat, dat.SigOn3 - dat.SigOff3, 'After', 'NetSig1', 'NewVariableNames', 'NetSig3');
  dat = addvars(dat, dat.SigOn2 - dat.SigOff2, 'After', 'NetSig1', 'NewVariableNames', 'NetSig2');
  dat = addvars(dat, dat.RefOn - dat.RefOff, 'Before', 'NetSig1', 'NewVariableNames', 'NetRef');
  
  
  dat = addvars(dat, dat.NetSig1 ./ dat.NetRef, 'Before', 'NetSig1', 'NewVariableNames', 'Scat1');
  dat = addvars(dat, dat.NetSig2 ./ dat.NetRef, 'Before', 'NetSig1', 'NewVariableNames', 'Scat2');
  dat = addvars(dat, dat.NetSig3 ./ dat.NetRef, 'Before', 'NetSig1', 'NewVariableNames', 'Scat3');
  
  % ScatX is the highest raw net signal from the front end that is not
  % saturated. GainX indicates the front end channel (1,2,3) assigned to
  % ScatX.
  % Start by assigning ScatX and GainX to highest front end gain level (3)
  ScatX = dat.Scat3;
  GainX = repmat(3, size(ScatX));
  if any(isnan(ScatX)) % If high gain is saturated, replace with low gain (2)
    GainX(isnan(ScatX)) = 2;
    ScatX(isnan(ScatX)) = dat.Scat2(isnan(ScatX)); 
  end
  if any(isnan(ScatX)) % it low gain is still saturated, replace with raw pmt no gain (1)
    GainX(isnan(ScatX)) = 1;
    ScatX(isnan(ScatX)) = dat.Scat1(isnan(ScatX)); 
  end
  
  dat = addvars(dat, ScatX, GainX, 'After', 'Scat3', 'NewVariableNames', {'ScatX', 'GainX'});
  
  % read cal files
  if ~all(isfield(hbb_cal, {'gain1_2','gain2_3','PMTGamma','PMTReferenceGain', ...
      'muFactorWl','muFactors','darkOffsetPMTGain','darkOffsetWl',...
      'darkOffsetScat1','darkOffsetScat2','darkOffsetScat3','muFactorLEDTemp'})) %  'H', 'rho', 
    error('Input cal struct does not contain required fields.')
  end
  % select calibration based on data dt
  idhbb_cal = find(min(dat.dt) > cell2mat({hbb_cal.date})',1,'last');
  idhbb_tcal = cell2mat({hbb_tcal.date})' == hbb_cal(idhbb_cal).tempCalDate;

  % % Display calibration mufactors
  % figure(600); clf; hold on
  % % convert dates to numeric
  % numeric_dates = posixtime(cell2mat({hbb_cal.date})');
  % % normalize dates between 0 and 1
  % norm_dates = (numeric_dates - min(numeric_dates)) / (max(numeric_dates) - min(numeric_dates));
  % % select a base colormap and extract colors at exact normalized date positions
  % base_colormap = colormap(jet(256));
  % % linearly interpolate colors matching date spacing
  % col = interp1(linspace(0, 1, 256), base_colormap, norm_dates);
  % for c = 1:size(hbb_cal,2)
  %   if idhbb_cal == c
  %     plot(hbb_cal(c).muFactorWl,hbb_cal(c).muFactors,'Color',col(c,:),'LineWidth',3)
  %   else
  %     plot(hbb_cal(c).muFactorWl,hbb_cal(c).muFactors,'Color',col(c,:))
  %   end
  % end
  % cb = colorbar;
  % clim([min(numeric_dates), max(numeric_dates)]); % lock color limits to timestamps
  % % place ticks based on calibration dates
  % cb.Ticks = numeric_dates;
  % cb.TickLabels = string(datetime(cell2mat({hbb_cal.date})','Format','yyyy-MM-dd'));
  % % write selected calibration in bold and larger font
  % ax = cb.Ruler;
  % for i = 1:length(ax.TickLabels)
  %   if idhbb_cal == i
  %     ax.TickLabels{i} = ['\bf\fontsize{12}' ax.TickLabels{i}]; 
  %   end
  % end
  % cb.Label.String = 'Calibration dates (selected calibration in bold)';
  % cb.Label.FontSize = 12;
  % xlabel('\lambda')
  % ylabel('\mu factor')
  % drawnow

  % Apply calibration
  if istable(dat)
    if sum(ismember(dat.Properties.VariableNames, {'Scat1', 'Scat2', 'Scat3', 'PmtGain'})) == 4
      dat = processDataTable(dat, hbb_cal(idhbb_cal), hbb_tcal(idhbb_tcal));
    else
      if verbose
        warning(['dat{' num2str(kd) '} does not contain proper data.'])
      end
    end
  else
    if verbose
      warning(['dat{' num2str(kd) '} is not a data table.'])
    end
  end
  
  % get wavelength from file and count occurences 
  [wl_counts, wlDat] = groupcounts(dat.wl);
  % interpolate hbb_cal to file wavelengths if not recorded at the same wavelength
  sel_hbb_cal = hbb_cal(idhbb_cal);
  sel_hbb_tcal = hbb_tcal(idhbb_tcal);
  % if at least two occurences of each wavelength => scan full: file lambda are all there
  % => interpolate mufactor on file lambda before storing calibration data in .mat files
  if all(wl_counts > 1) % different lambda in temperature calibration not tested, to be added here eventually
    sel_hbb_cal.muFactors = interp1(sel_hbb_cal.muFactorWl, sel_hbb_cal.muFactors, wlDat, 'pchip');
    sel_hbb_cal.muFactorTempCorr = interp1(sel_hbb_cal.muFactorWl, sel_hbb_cal.muFactorTempCorr, wlDat, 'pchip');
    sel_hbb_cal.muFactorLEDTemp = interp1(sel_hbb_cal.muFactorWl, sel_hbb_cal.muFactorLEDTemp, wlDat, 'pchip');
    sel_hbb_cal.muFactorWl = wlDat;
    lambda = wlDat;
  elseif mean(diff(sel_hbb_tcal.wl)) == mean(diff(wlDat))
    lambda = sel_hbb_tcal.wl;
  elseif mean(diff(sel_hbb_cal.muFactorWl)) == mean(diff(wlDat))
    lambda = sel_hbb_cal.muFactorWl;
  else
    lambda = wlDat;
    warning('Scan incomplete, complete lambda vector could not be retrieved')
  end

  % reshape data into clean table
  [data, lambda] = reformatHBB(dat, lambda);
  data(all(isnan(data.beta), 2), :) = [];
  
  % save lambda and calibration information in table properties
  data = addprop(data, {'lambda','hbb_cal','hbb_tcal'}, {'table','table','table'});
  data.Properties.CustomProperties.lambda = lambda;
  data.Properties.UserData = lambda;
  data.Properties.CustomProperties.hbb_cal = sel_hbb_cal;
  data.Properties.CustomProperties.hbb_tcal = sel_hbb_tcal;
  
  if verbose; fprintf('Done\n'); end
end

function dat = processDataTable(dat, hbb_cal, hbb_tcal)

  % Interpolate wl and pmt to find dark offset
  darkOffset_scat1 = interp2(hbb_cal.darkOffsetPMTGain, hbb_cal.darkOffsetWl, ...
    hbb_cal.darkOffsetScat1, dat.PmtGain, dat.wl, 'linear');
  darkOffset_scat2 = interp2(hbb_cal.darkOffsetPMTGain, hbb_cal.darkOffsetWl, ...
    hbb_cal.darkOffsetScat2, dat.PmtGain, dat.wl, 'linear');
  darkOffset_scat3 = interp2(hbb_cal.darkOffsetPMTGain, hbb_cal.darkOffsetWl, ...
    hbb_cal.darkOffsetScat3, dat.PmtGain,dat.wl, 'linear');
  
  % Subtract Dark Offset
  scat1_darkRemoved = dat.Scat1 - darkOffset_scat1;
  scat2_darkRemoved = dat.Scat2 - darkOffset_scat2;
  scat3_darkRemoved = dat.Scat3 - darkOffset_scat3;
  
  % Apply PMT and front end gain factors
  Gpmt = (dat.PmtGain ./ hbb_cal.PMTReferenceGain) .^ hbb_cal.PMTGamma;
  scat1_gainCorrected = scat1_darkRemoved .* hbb_cal.gain1_2 .* hbb_cal.gain2_3 .* Gpmt;
  scat2_gainCorrected = scat2_darkRemoved .* hbb_cal.gain2_3 .* Gpmt;
  scat3_gainCorrected = scat3_darkRemoved .* Gpmt;
  
  % Add dat to output table
  dat = addvars(dat, scat1_gainCorrected, 'After', 'Scat3', 'NewVariableNames', 'ScatCor1');
  dat = addvars(dat, scat2_gainCorrected, 'After', 'ScatCor1', 'NewVariableNames', 'ScatCor2');
  dat = addvars(dat, scat3_gainCorrected, 'After', 'ScatCor2', 'NewVariableNames', 'ScatCor3');
  
  % Apply temperature correction
  tempCoeff = GetTemperatureCoefficients(hbb_tcal, dat.wl, dat.LedTemp);
  scat1_tempCorrected = dat.ScatCor1 .* tempCoeff;
  scat2_tempCorrected = dat.ScatCor2 .* tempCoeff;
  scat3_tempCorrected = dat.ScatCor3 .* tempCoeff;
  
  % Add data to output table
  dat = addvars(dat, scat1_tempCorrected, 'After','ScatCor3', 'NewVariableNames', 'ScatTempCor1');
  dat = addvars(dat, scat2_tempCorrected, 'After','ScatTempCor1', 'NewVariableNames', 'ScatTempCor2');
  dat = addvars(dat, scat3_tempCorrected, 'After','ScatTempCor2', 'NewVariableNames', 'ScatTempCor3');
  dat = addvars(dat, tempCoeff, 'After','ScatTempCor3', 'NewVariableNames', 'TempCorrCoeff'); % save the calculated temp coeff.
  
  % Select highest non-saturated gain channel
  ScatTempCorX = dat.ScatTempCor3;
  % If high gain is saturated, replace with low gain
  ScatTempCorX(isnan(ScatTempCorX)) = dat.ScatTempCor2(isnan(ScatTempCorX));
  % it low gain is saturated, replace with raw pmt
  ScatTempCorX(isnan(ScatTempCorX)) = dat.ScatTempCor1(isnan(ScatTempCorX));
  dat = addvars(dat, ScatTempCorX, 'After', 'ScatTempCor3', 'NewVariableNames', 'ScatTempCorX');
  
  % Temperature correct mu calibration
  tempCoeff_mu = GetTemperatureCoefficients(hbb_tcal, hbb_cal.muFactorWl, hbb_cal.muFactorLEDTemp);
  muFactors_tempCorrected = hbb_cal.muFactors .* tempCoeff_mu;
  
  % Calculate Beta
  wlDat = sort(unique(dat.wl))';
  muFactors = interp1(hbb_cal.muFactorWl, muFactors_tempCorrected, wlDat, 'pchip');
  dat = addvars(dat, NaN(height(dat),1), 'After','ScatTempCorX', 'NewVariableNames', 'beta'); % preallocate array
  for kw = 1:length(wlDat)
    dat.beta(dat.wl == wlDat(kw)) = dat.ScatTempCorX(dat.wl == wlDat(kw)) .* muFactors(kw);
  end        
end

function tempCoeff = GetTemperatureCoefficients(hbb_tcal, wavelength, temperature) 
  % Generate temperature correction grid
  LEDTempRange = min(temperature):0.1:max(temperature) + 0.1; % need to make sure the max value is included
  TempCorrGrid = NaN(length(hbb_tcal.wl), length(LEDTempRange));
  for n = 1:length(hbb_tcal.wl)
    TempCorrGrid(n,:) = polyval(hbb_tcal.coeff(n,:), LEDTempRange);
  end
  
  % 2D interpolate (wavelength x temperature) to find temperature
  % correction factor
  tempCoeff = interp2(LEDTempRange, hbb_tcal.wl, TempCorrGrid, temperature, wavelength, 'linear');
end

function [data, lambda] = reformatHBB(dat, lambda)
  % Reformat HBB data to get lambda as column
  uscan = unique(dat.ScanIdx);
  
  foo = dat.wl;
  for i = 1:max(size(lambda))
    foo(foo == lambda(i)) = i;
  end
  
  data = table();
  data.dt = NaT(size(uscan, 1), 1);
  data.ScanIdx = uscan;
  data.WaterTemp = NaN(size(uscan, 1), 1);
  data.Depth = NaN(size(uscan, 1), 1);
  data.beta = NaN(size(uscan, 1), max(size(lambda)));
  data.beta_u = NaN(size(uscan, 1), max(size(lambda)));
  data.bb_inlinino = NaN(size(uscan, 1), max(size(lambda)));
  
  for i = 1:size(uscan, 1)
    data.dt(i) = median(dat.dt(dat.ScanIdx == uscan(i)),'omitnan');
    data.WaterTemp(i) = median(dat.WaterTemp(dat.ScanIdx == uscan(i)),'omitnan');
    data.Depth(i) = median(dat.Depth(dat.ScanIdx == uscan(i)),'omitnan');
    data.beta(i, foo(dat.ScanIdx == uscan(i))) = dat.beta(dat.ScanIdx == uscan(i))';
    data.beta_u(i, foo(dat.ScanIdx == uscan(i))) = dat.beta_u(dat.ScanIdx == uscan(i))';
    data.bb_inlinino(i, foo(dat.ScanIdx == uscan(i))) = dat.bb(dat.ScanIdx == uscan(i))';
  end
  data = sortrows(data, 'dt');
end

function [ data, lambda ] = importInlininoHBB( filename, calfile_plaque, calfile_temp, verbose )
% Import HyperBB data logged with Inlinino and apply SEQUOIA's processing code
% using calfile_plaque and calfile_temp
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

p.RemoveMultiplePmtGains = false;
try
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

  % Set parser & read data
%   readtable(filename)
  try
    parser = '%s%f%f%{yyyy/MM/dd}D%{hh:mm:ss}T%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
    t = textscan(fid, parser, 'delimiter', ',');
  catch
    fid=fopen(filename);
    % Get header
    hd = strip(strsplit(fgetl(fid), ','));
    hd(strcmp(hd, 'time')) = {'dt'};
    % get units skipping empty lines (bug in Inlinino)
    unit = fgetl(fid);
    while isempty(unit)
        unit = fgetl(fid);
    end
    unit = strip(strsplit(unit, ','));
    parser = '%s%f%f%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
    t = textscan(fid, parser, 'delimiter', ',');
    t{4} = NaT(size(t{1}));
    t{5} = NaT(size(t{1}));
  end
  % Close file
  fclose(fid);

  % Build table
  dat = table(datenum(t{1}), t{2}, t{3}, t{6}, t{7}, t{8}, t{9}, t{10}, ...
    t{11}, t{12}, t{13}, t{14}, t{15}, t{16}, t{17}, t{18}, t{19}, t{20}, ...
    t{21}, t{22}, t{23}, t{24}, t{25}, t{26}, t{27}, t{28}, t{29}, t{30}, t{31}, ...
    'VariableNames', hd(~contains(hd, {'Date', 'Time'})));
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
  
  t = [];
  i = 0;
  while all(~contains(tline, 'Date, TIME, POS, WL, PWR, GAIN, SIG, AVG, STD,') & ~strcmp(tline, '-1'))
    t = [t; strsplit(tline, ', ')];
    tline = fgetl(fid);
    if isnumeric(tline)
      tline = num2str(tline);
    end
    if i >= round(i, -3)
      warning('%s lines checked ...\n', num2str(i))
    end
    i = i + 1;
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
  dat = table(datenum(tend{1}), tend{2}, tend{3}, tend{6}, tend{7}, tend{8}, tend{9}, tend{10}, ...
    tend{11}, tend{12}, tend{13}, tend{14}, tend{15}, tend{16}, tend{17}, tend{18}, tend{19}, tend{20}, ...
    tend{21}, tend{22}, tend{23}, tend{24}, tend{25}, tend{26}, tend{27}, tend{28}, tend{29}, tend{30}, tend{31}, ...
    'VariableNames', hd(~contains(hd, {'Date', 'Time'})));
end
  
dat(isnan(dat.dt), :) = [];
% dat.Properties.VariableUnits = unit(~contains(hd, {'Date', 'Time'}));

for ks = unique(dat.ScanIdx')
  scanSel = dat.ScanIdx==ks;
  if length(unique(dat.PmtGain(scanSel)))~=1
    if p.RemoveMultiplePmtGains
      if verbose
        warning(['  Scan ' num2str(ks) ' has multiple PmtGain, REMOVING'])
      end
      dat(scanSel,:) = [];
    else
      if verbose
        warning(['  Scan ' num2str(ks) ' has multiple PmtGain'])
      end
    end
  end
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
if ~all(isfield(calfile_plaque, {'gain12', 'gain23', 'pmtGamma', 'pmtRefGain', 'H', 'rho', 'muWavelengths', 'muFactors'}))
  error('Input cal struct does not contain required fields.')
end

if istable(dat)
  if sum(ismember(dat.Properties.VariableNames, {'Scat1', 'Scat2', 'Scat3', 'PmtGain'})) == 4
    dat = processDataTable(dat, calfile_plaque, calfile_temp);
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

[data, lambda] = reformatHBB(dat, calfile_temp.wl);
data(all(isnan(data.beta), 2), :) = [];

% save calibration information in table properties
data = addprop(data, {'PlaqueCal', 'TemperatureCal'}, ...
              {'table', 'table'});

data.Properties.CustomProperties.PlaqueCal = calfile_plaque;
data.Properties.CustomProperties.TemperatureCal = calfile_temp;

if verbose; fprintf('Done\n'); end
end 

function dat = processDataTable(dat, cal_plaque, cal_temp)

% Interpolate wl and pmt to find dark offset
darkOffset_scat1 = interp2(cal_plaque.darkCalPmtGain, cal_plaque.darkCalWavelength, ...
  cal_plaque.darkCalScat1, dat.PmtGain, dat.wl, 'linear');
darkOffset_scat2 = interp2(cal_plaque.darkCalPmtGain, cal_plaque.darkCalWavelength, ...
  cal_plaque.darkCalScat2, dat.PmtGain, dat.wl, 'linear');
darkOffset_scat3 = interp2(cal_plaque.darkCalPmtGain, cal_plaque.darkCalWavelength, ...
  cal_plaque.darkCalScat3, dat.PmtGain,dat.wl, 'linear');

% Subtract Dark Offset
scat1_darkRemoved = dat.Scat1 - darkOffset_scat1;
scat2_darkRemoved = dat.Scat2 - darkOffset_scat2;
scat3_darkRemoved = dat.Scat3 - darkOffset_scat3;

% Apply PMT and front end gain factors
Gpmt = (dat.PmtGain ./ cal_plaque.pmtRefGain) .^ cal_plaque.pmtGamma;
scat1_gainCorrected = scat1_darkRemoved .* cal_plaque.gain12 .* cal_plaque.gain23 .* Gpmt;
scat2_gainCorrected = scat2_darkRemoved .* cal_plaque.gain23 .* Gpmt;
scat3_gainCorrected = scat3_darkRemoved .* Gpmt;

% Add dat to output table
dat = addvars(dat, scat1_gainCorrected, 'After', 'Scat3', 'NewVariableNames', 'ScatCor1');
dat = addvars(dat, scat2_gainCorrected, 'After', 'ScatCor1', 'NewVariableNames', 'ScatCor2');
dat = addvars(dat, scat3_gainCorrected, 'After', 'ScatCor2', 'NewVariableNames', 'ScatCor3');

% Apply temperature correction
tempCoeff = GetTemperatureCoefficients(cal_temp, dat.wl, dat.LedTemp);
scat1_tempCorrected = dat.ScatCor1 .* tempCoeff;
scat2_tempCorrected = dat.ScatCor2 .* tempCoeff;
scat3_tempCorrected = dat.ScatCor3 .* tempCoeff;

% Add data to output table
dat = addvars(dat, scat1_tempCorrected, 'After','ScatCor3', ...
  'NewVariableNames', 'ScatTempCor1');
dat = addvars(dat, scat2_tempCorrected, 'After','ScatTempCor1', ...
  'NewVariableNames', 'ScatTempCor2');
dat = addvars(dat, scat3_tempCorrected, 'After','ScatTempCor2', ...
  'NewVariableNames', 'ScatTempCor3');
dat = addvars(dat, tempCoeff, 'After','ScatTempCor3', 'NewVariableNames', ...
  'TempCorrCoeff'); % save the calculated temp coeff.

% Select highest non-saturated gain channel
ScatTempCorX = dat.ScatTempCor3;
% If high gain is saturated, replace with low gain
ScatTempCorX(isnan(ScatTempCorX)) = dat.ScatTempCor2(isnan(ScatTempCorX));
% it low gain is saturated, replace with raw pmt
ScatTempCorX(isnan(ScatTempCorX)) = dat.ScatTempCor1(isnan(ScatTempCorX));
dat = addvars(dat, ScatTempCorX, 'After', 'ScatTempCor3', 'NewVariableNames', 'ScatTempCorX');

% Temperature correct mu calibration
tempCoeff_mu = GetTemperatureCoefficients(cal_temp, cal_plaque.muWavelengths, ...
  cal_plaque.muLedTemp);
muFactors_tempCorrected = cal_plaque.muFactors .* tempCoeff_mu;

% Calculate Beta
wlDat = sort(unique(dat.wl))';
muFactors = interp1(cal_plaque.muWavelengths, muFactors_tempCorrected, wlDat, 'pchip');
dat = addvars(dat, NaN(height(dat),1), 'After','ScatTempCorX', 'NewVariableNames', 'beta'); % preallocate array
for kw = 1:length(wlDat)
  dat.beta(dat.wl == wlDat(kw)) = dat.ScatTempCorX(dat.wl == wlDat(kw)) .* muFactors(kw);
end        
end

function tempCoeff = GetTemperatureCoefficients(cal_temp, wavelength, temperature) 
% Generate temperature correction grid
LEDTempRange = min(temperature):0.1:max(temperature) + 0.1; % need to make sure the max value is included
TempCorrGrid = NaN(length(cal_temp.wl), length(LEDTempRange));
for n = 1:length(cal_temp.wl)
  TempCorrGrid(n,:) = polyval(cal_temp.coeff(n,:), LEDTempRange);
end

% 2D interpolate (wavelength x temperature) to find temperature
% correction factor
tempCoeff = interp2(LEDTempRange, cal_temp.wl, TempCorrGrid, temperature, wavelength, 'linear');
end

function [data, lambda] = reformatHBB(dat, lambda)
% Reformat HBB data to get lambda as column
uscan = unique(dat.ScanIdx);

foo = dat.wl;
for i = 1:max(size(lambda))
  foo(foo == lambda(i)) = i;
end

data = table();
data.dt = NaN(size(uscan, 1), 1);
data.ScanIdx = uscan;
data.WaterTemp = NaN(size(uscan, 1), 1);
data.Depth = NaN(size(uscan, 1), 1);
data.beta = NaN(size(uscan, 1), max(size(lambda)));

for i = 1:size(uscan, 1)
  data.dt(i) = median(dat.dt(dat.ScanIdx == uscan(i)),'omitnan');
  data.WaterTemp(i) = median(dat.WaterTemp(dat.ScanIdx == uscan(i)),'omitnan');
  data.Depth(i) = median(dat.Depth(dat.ScanIdx == uscan(i)),'omitnan');
  data.beta(i, foo(dat.ScanIdx == uscan(i))) = dat.beta(dat.ScanIdx == uscan(i))';
end
data = sortrows(data, 'dt');
end

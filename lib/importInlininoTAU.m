function [ data ] = importInlininoTAU( filename, verbose )
% Import LISST-Tau data logged with Inlinino
%
%%
if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

data = readtable(filename,'FileType','text','Delimiter',{'\t', ','});
data.Properties.VariableNames = {'dt', 'Name','InstrTime','Seconds','RefNet','RefOn','RefOff','SigNet','SigOn','SigOff','SigSat','RefNetStd','SigNetStd', ...
    'Tr','CorrFunTemp','TrCorr','Tau','Beamc','BeamcStd','NumSamp','DacLedOn','DacLedOff', ...
    'TempLedRaw','TempMainRaw','TempRcvrRaw', 'VinRaw','TempLed','TempMain','TempRcvr','Vin', ...
    'Analog1Raw','Analog2Raw','Analog1','Analog2', ...
    'InstrFirmware', 'CalTime','TrCal','TempCal','CorrFunCal'};

data.Name = [];
data.dt = datenum(data.dt);
data(isnan(data.dt), :) = [];
data.InstrTime = datenum(data.InstrTime);
data.CalTime = datenum(data.CalTime);

if verbose; fprintf('Done\n'); end
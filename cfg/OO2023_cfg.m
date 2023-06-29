% OO2023 Configuration file
% author: Guilaume Bourdin
% created: June 27, 2023

cfg = struct('meta', struct(), 'instruments', struct(), 'process', struct());

%%%%%%%%%%%%%%
%% METADATA %%
%%%%%%%%%%%%%%
cfg.meta.investigators = 'OO2023_participants';
cfg.meta.affiliations = 'University_of_Maine';
cfg.meta.emails = 'emmanuel.boss@maine.edu';
cfg.meta.experiment = 'OO2023';
cfg.meta.cruise = 'OO2023';
cfg.meta.station = 'NA';
cfg.meta.documents = 'NA';
cfg.meta.calibration_files = 'NA';
cfg.meta.data_type = 'flow_thru';
cfg.meta.data_status = 'preliminary';
cfg.meta.measurement_depth = 1.5;


%%%%%%%%%%%%%%%%%
%% INSTRUMENTS %%
%%%%%%%%%%%%%%%%%

PATH_ROOT = '/Volumes/Samsung_T5/Data/OO2023';

%%% TSG %%%
model = 'SBE45';
SN = '36073';
cfg.instruments.([model SN]) = struct();
cfg.instruments.([model SN]).model = model;
cfg.instruments.([model SN]).boat = 'Tara';
cfg.instruments.([model SN]).logger = 'Inlinino'; % TeraTerm Matlab Inlinino
cfg.instruments.([model SN]).sn = SN;
cfg.instruments.([model SN]).path = struct('raw',  fullfile(PATH_ROOT, 'raw', [model SN]),...
                                  'wk',   fullfile(PATH_ROOT, 'wk', [model SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', [model SN]));
cfg.instruments.([model SN]).view = struct('varname', 't');
cfg.instruments.([model SN]).temperature_variable = 't';

%%% atlasTSG002 %%%
model = 'atlasTSG';
SN = '002';
cfg.instruments.([model SN]) = struct();
cfg.instruments.([model SN]).model = 'atlasTSG';
cfg.instruments.([model SN]).TSseparated = true;
cfg.instruments.([model SN]).coef = struct('t', struct('slope',1,'interc',0), ...
                                  'c', struct('slope',1,'interc',0));
cfg.instruments.([model SN]).logger = 'Inlinino_atlasTSG';
cfg.instruments.([model SN]).sn = SN;
cfg.instruments.([model SN]).prefix = 'ATLAS';
cfg.instruments.([model SN]).path = struct('raw',  [PATH_ROOT 'raw' filesep 'atlasTSG' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'atlasTSG' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'atlasTSG' filesep]);
cfg.instruments.([model SN]).view = struct('varname', 't');

%%% atlasTSG003 %%%
model = 'atlasTSG';
SN = '003';
cfg.instruments.([model SN]) = struct();
cfg.instruments.([model SN]).model = 'atlasTSG';
cfg.instruments.([model SN]).TSseparated = false;
cfg.instruments.([model SN]).coef = struct('t', struct('slope',1,'interc',0), ...
                                  'c', struct('slope',1,'interc',0));
cfg.instruments.([model SN]).logger = 'Inlinino_atlasTSG';
cfg.instruments.([model SN]).sn = SN;
cfg.instruments.([model SN]).prefix = 'ATLASECRTD';
cfg.instruments.([model SN]).path = struct('raw',  [PATH_ROOT 'raw' filesep 'atlasTSG' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'atlasTSG' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'atlasTSG' filesep]);
cfg.instruments.([model SN]).view = struct('varname', 't');

%%% atlasTSG004 %%%
model = 'atlasTSG';
SN = '004';
cfg.instruments.([model SN]) = struct();
cfg.instruments.([model SN]).model = 'atlasTSG';
cfg.instruments.([model SN]).TSseparated = false;
cfg.instruments.([model SN]).coef = struct('t', struct('slope',1,'interc',0), ...
                                  'c', struct('slope',1,'interc',0));
cfg.instruments.([model SN]).logger = 'Inlinino_atlasTSG';
cfg.instruments.([model SN]).sn = SN;
cfg.instruments.([model SN]).prefix = 'ATLASECRTD';
cfg.instruments.([model SN]).path = struct('raw',  [PATH_ROOT 'raw' filesep 'atlasTSG' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'atlasTSG' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'atlasTSG' filesep]);
cfg.instruments.([model SN]).view = struct('varname', 't');

%%% NMEA %%%
model = '19X'; % GPSSC701 GPSGP32
SN = '2J5055713';
cfg.instruments.NMEA = struct();
cfg.instruments.NMEA.model = model;
cfg.instruments.NMEA.boat = 'Ira-C';
cfg.instruments.NMEA.logger = 'Inlinino';
cfg.instruments.NMEA.sn = SN;
cfg.instruments.NMEA.prefix = 'GPS19x';
cfg.instruments.NMEA.path = struct('raw',  fullfile(PATH_ROOT, 'raw', ['GPS' cfg.instruments.NMEA.model]),...
                                  'wk',   fullfile(PATH_ROOT, 'wk', ['GPS' cfg.instruments.NMEA.model]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', ['GPS' cfg.instruments.NMEA.model]));
cfg.instruments.NMEA.view = struct('varname', 'lat');

%%% FLOW (FlowControl) %%%
SN = 'B02782';
cfg.instruments.FLOW = struct();
cfg.instruments.FLOW.model = 'ADU100';
cfg.instruments.FLOW.logger = 'Inlinino';
cfg.instruments.FLOW.sn = SN;
cfg.instruments.FLOW.LoadPrevious = true;
cfg.instruments.FLOW.analog1 = '';
cfg.instruments.FLOW.analog2 = 'SUVF6254';
cfg.instruments.FLOW.path = struct('raw',  fullfile(PATH_ROOT, 'raw', 'FlowControl'),...
                                  'wk',   fullfile(PATH_ROOT, 'wk', 'FlowControl'),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', 'FlowControl'));
cfg.instruments.FLOW.view = struct('varname', 'swt','spd_variable','spd1'); % spd1 spd2

%%% ACS412 %%%
SN = '412';
cfg.instruments.(['ACS' SN]) = struct();
cfg.instruments.(['ACS' SN]).di = struct();
cfg.instruments.(['ACS' SN]).di.prefix = ['DIW_ACS' SN '_'];
cfg.instruments.(['ACS' SN]).di.postfix = '';
cfg.instruments.(['ACS' SN]).model = 'ACS';
cfg.instruments.(['ACS' SN]).sn = SN;
cfg.instruments.(['ACS' SN]).logger = 'InlininoACScsv';
cfg.instruments.(['ACS' SN]).device_file = fullfile(PATH_ROOT, 'DeviceFiles', 'ACS0412_20230510.DEV');
cfg.instruments.(['ACS' SN]).path = struct('raw',  fullfile(PATH_ROOT, 'raw', ['ACS' SN]),...
                                    'di',  fullfile(PATH_ROOT, 'raw', ['ACS' SN], 'DI'),...
                                    'wk',   fullfile(PATH_ROOT, 'wk', ['ACS' SN]),...
                                    'prod', fullfile(PATH_ROOT, 'prod'),...
                                    'ui', fullfile(PATH_ROOT, 'ui', ['ACS' SN]));
cfg.instruments.(['ACS' SN]).view = struct('varname', 'a', 'varcol', 40);

%%% HBB %%%
SN = '8002';
cfg.instruments.(['HyperBB' SN]) = struct();
cfg.instruments.(['HyperBB' SN]).di = struct();
cfg.instruments.(['HyperBB' SN]).di.prefix = ['DIW_HyperBB' SN '_'];
cfg.instruments.(['HyperBB' SN]).di.postfix = '';
cfg.instruments.(['HyperBB' SN]).model = 'HBB';
cfg.instruments.(['HyperBB' SN]).sn = SN;
cfg.instruments.(['HyperBB' SN]).ila_prefix = ['HyperBB' SN];
cfg.instruments.(['HyperBB' SN]).logger = 'InlininoHBB';
cfg.instruments.(['HyperBB' SN]).theta = 135;
cfg.instruments.(['HyperBB' SN]).PlaqueCal = fullfile(PATH_ROOT, 'DeviceFiles', 'Hbb_Cal_Plaque_20230606_120731.mat');
cfg.instruments.(['HyperBB' SN]).TemperatureCal = fullfile(PATH_ROOT, 'DeviceFiles', 'Hbb_Cal_Temp_20230601_145659.mat');
cfg.instruments.(['HyperBB' SN]).path = struct('raw',  fullfile(PATH_ROOT, 'raw', ['HyperBB' SN]),...
                                  'di',  fullfile(PATH_ROOT, 'raw', ['HyperBB' SN], 'DI'),...
                                  'wk',   fullfile(PATH_ROOT, 'wk', ['HyperBB' SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', ['HyperBB' SN]));
cfg.instruments.(['HyperBB' SN]).view = struct('varname', 'beta', 'varcol', 23);

%%% BB3 %%%
SN = '349';
cfg.instruments.(['BB3' SN]) = struct();
cfg.instruments.(['BB3' SN]).di = struct();
cfg.instruments.(['BB3' SN]).di.prefix = ['DIW_BB3' SN '_'];
cfg.instruments.(['BB3' SN]).di.postfix = '';
cfg.instruments.(['BB3' SN]).model = 'BB';
cfg.instruments.(['BB3' SN]).sn = SN;
cfg.instruments.(['BB3' SN]).ila_prefix = 'BB3';
cfg.instruments.(['BB3' SN]).logger = 'InlininoBB3SN';
cfg.instruments.(['BB3' SN]).lambda = [468,528,655];
cfg.instruments.(['BB3' SN]).theta = 120;
cfg.instruments.(['BB3' SN]).slope = [8.407E-06,4.624E-06,4.090E-06]; % Emmanuel cal 2021
cfg.instruments.(['BB3' SN]).dark = [50, 44, 46];
cfg.instruments.(['BB3' SN]).path = struct('raw',  fullfile(PATH_ROOT, 'raw', ['BB3' SN]),...
                                  'di',  fullfile(PATH_ROOT, 'raw', ['BB3' SN]),...
                                  'wk',   fullfile(PATH_ROOT, 'wk', ['BB3' SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', ['BB3' SN]));
cfg.instruments.(['BB3' SN]).view = struct('varname', 'beta', 'varcol', 2);

%%% SUVF %%%
SN = '6254';
cfg.instruments.(['SUVF' SN]) = struct();
cfg.instruments.(['SUVF' SN]).model = 'CD';
cfg.instruments.(['SUVF' SN]).sn = SN;
cfg.instruments.(['SUVF' SN]).slope = 1;
cfg.instruments.(['SUVF' SN]).dark = 0;
cfg.instruments.(['SUVF' SN]).logger = 'InlininoADU100';
cfg.instruments.(['SUVF' SN]).analog_channel = 'Analog(2)';
if strcmp(cfg.instruments.(['SUVF' SN]).logger, 'InlininoADU100')
  cfg.instruments.(['SUVF' SN]).ila_prefix = [cfg.instruments.FLOW.model cfg.instruments.FLOW.sn '_'];
  [~, raw_folder_name] = fileparts(cfg.instruments.FLOW.path.raw);
else
  cfg.instruments.(['SUVF' SN]).ila_prefix = ['SUVF' SN];
  raw_folder_name = ['SUVF' SN];
end
cfg.instruments.(['SUVF' SN]).di = struct();
cfg.instruments.(['SUVF' SN]).di.prefix = ['DIW_' cfg.instruments.(['SUVF' SN]).ila_prefix SN '_'];
cfg.instruments.(['SUVF' SN]).di.postfix = '';
cfg.instruments.(['SUVF' SN]).path = struct('raw',  fullfile(PATH_ROOT, 'raw', raw_folder_name),...
                                  'di',  fullfile(PATH_ROOT, 'raw', ['SUVF' SN], 'DI'),...
                                  'wk',   fullfile(PATH_ROOT, 'wk', ['SUVF' SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', ['SUVF' SN]));
cfg.instruments.(['SUVF' SN]).view = struct('varname', 'fdom', 'varcol', 1);

%%% WS3S %%%
SN = '1081';
cfg.instruments.(['WS3S' SN]) = struct();
cfg.instruments.(['WS3S' SN]).sn = SN;
cfg.instruments.(['WS3S' SN]).model = 'FL';
cfg.instruments.(['WS3S' SN]).ila_prefix = ['WS3S' SN];
cfg.instruments.(['WS3S' SN]).logger = 'InlininoWSCDSN';
cfg.instruments.(['WS3S' SN]).slope = 5.6;
cfg.instruments.(['WS3S' SN]).dark = 0.068;
cfg.instruments.(['WS3S' SN]).path = struct('raw',  fullfile(PATH_ROOT, 'raw', ['WS3S' SN]),...
                                  'wk',   fullfile(PATH_ROOT, 'wk', ['WS3S' SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', ['WS3S' SN]));
cfg.instruments.(['WS3S' SN]).view = struct('varname', 'fchl');

%%% LISST %%%
SN = '1183';
cfg.instruments.(['LISST' SN]) = struct();
cfg.instruments.(['LISST' SN]).di = struct();
cfg.instruments.(['LISST' SN]).di.prefix = ['DIW_LISST' SN '_'];
cfg.instruments.(['LISST' SN]).di.postfix = '';
cfg.instruments.(['LISST' SN]).model = 'LISST';
cfg.instruments.(['LISST' SN]).type = 'B'; 
cfg.instruments.(['LISST' SN]).sn = SN;
cfg.instruments.(['LISST' SN]).logger = 'InlininoLISSTcsv';
cfg.instruments.(['LISST' SN]).zsc = [2.203500e+001, 2.568500e+001, 2.503000e+001, 2.986000e+001, 2.842500e+001, 3.283000e+001, 3.077000e+001, 3.659500e+001, 2.978000e+001, 3.552000e+001, 3.198000e+001, 4.216000e+001, 3.916500e+001, 4.662500e+001, 3.974000e+001, 4.454000e+001, 4.403500e+001, 4.604500e+001, 4.430000e+001, 4.510500e+001, 4.719500e+001, 3.850000e+001, 5.373000e+001, 2.664000e+001, 3.180500e+001, 1.655500e+001, 2.205500e+001, 1.554000e+001, 1.422000e+001, 1.123000e+001, 8.780000e+000, 8.555000e+000, 1.515000e+003, 1.167900e+003, 6.410000e+001, 1.055150e+003, 7.700000e+001, 2.116600e+003, 1.807000e+003, 5.476500e+003];
% Original dcal file (ring area)
% cfg.instruments.(['LISST' SN]).dcal = [1.0000000e+000, 1.0038000e+000, 9.9360000e-001, 1.0027000e+000, 9.9720000e-001, 9.9570000e-001, 9.9030000e-001, 9.9430000e-001, 9.9290000e-001, 9.9000000e-001, 9.9290000e-001, 9.9300000e-001, 9.9150000e-001, 9.9300000e-001, 9.9230000e-001, 9.9090000e-001, 1.1032000e+000, 1.1123000e+000, 1.2430000e+000, 1.1562000e+000, 1.3273000e+000, 1.1999000e+000, 1.0740000e+000, 1.7489000e+000, 1.5382000e+000, 2.5109000e+000, 2.5468000e+000, 3.5504000e+000, 3.9338000e+000, 5.1747342e+000, 7.5143548e+000, 1.2528083e+001];
% dcal adhoc from email of Wayne Slade Dec 8, 2017
cfg.instruments.(['LISST' SN]).dcal = [ 1.0179083e+00	   9.9213489e-01	   1.0108161e+00	   9.9492883e-01	   1.0043707e+00	   9.9891840e-01	   9.9859055e-01	   1.0042049e+00	   1.0000763e+00	   9.9889997e-01	   1.0009497e+00	   1.0004019e+00	   1.0011130e+00	   1.0004677e+00	   1.0213554e+00	   9.9990262e-01	   1.1115630e+00	   1.1206668e+00	   1.2493699e+00	   1.1643199e+00	   1.3355657e+00	   1.2090892e+00	   1.0781540e+00	   1.7620752e+00	   1.5508563e+00	   2.5304119e+00	   2.5638592e+00	   3.5757212e+00	   3.9631987e+00	   5.0166411e+00	   5.6381118e+00	   8.6881539e+00];
% customize dcal of 1183 on Jan 16, 2018 due to bump in VSF at ring 30
cfg.instruments.(['LISST' SN]).dcal(30) = 2.6; % Good for low values but bad for high values ??
cfg.instruments.(['LISST' SN]).vcc = 48493;
cfg.instruments.(['LISST' SN]).inversion = 'spherical';
cfg.instruments.(['LISST' SN]).ds = [1.2500,1.4750,1.7405,2.0538,2.4235,2.8597,3.3744,3.9818,4.6986,5.5443,6.5423,7.7199,9.1095,10.7492,12.6841,14.9672,17.6613,20.8403,24.5916,29.0180,34.2413,40.4047,47.6776,56.2595,66.3863,78.3358,92.4362,109.0747,128.7082,151.8757,179.2133,211.4717,249.5366];
cfg.instruments.(['LISST' SN]).theta = [0.082, 0.096, 0.114, 0.134, 0.158, 0.187, 0.221, 0.260, 0.307, 0.362, 0.428, 0.505, 0.596, 0.703, 0.829, 0.979, 1.155, 1.363, 1.609, 1.898, 2.240, 2.643, 3.119, 3.681, 4.344, 5.126, 6.049, 7.138, 8.424, 9.941, 11.73, 13.84];
cfg.instruments.(['LISST' SN]).path = struct('raw',  fullfile(PATH_ROOT, 'raw', ['LISST' SN]),...
                                  'wk',   fullfile(PATH_ROOT, 'wk', ['LISST' SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', ['LISST' SN]));
cfg.instruments.(['LISST' SN]).view = struct('varname', 'beta', 'varcol', 15);


%% %%%%%%% %%
%  PROCESS  %
%%%%%%%%%%%%%

%%% General parameters %%%
cfg.process.days2run = datenum(2020,12,12):datenum(2022,10,15);
% cfg.process.instruments2run = {'FLOW', 'NMEA', 'ACS57', 'BB31502', 'WSCD859', ...
%   'SBE4536073', 'SUVF6244', 'LISST1183', 'HyperBB8005'};
cfg.process.instruments2run = fieldnames(cfg.instruments);
cfg.process.write = true;
cfg.process.force_import = false;
cfg.process.parallel = Inf; % 0: disable parallel or Inf: as many thread available
cfg.process.di = struct();
cfg.process.di.skip = cfg.process.instruments2run(contains(cfg.process.instruments2run, ...
  {'FLOW','TSG','SBE45','SBE3845','NMEA','PAR'}));
cfg.process.di.qc = struct('mode', 'ui',... % ui or load
                           'qc_once_for_all', false,... % true = QC all variables | false = QC variables separately);
                           'remove_old', false); % remove old selection of the same period
cfg.process.di.bin = struct('bin_size', 30);

%%% Synchronization %%%
cfg.process.sync = struct();
cfg.process.sync.delay = struct();
cfg.process.sync.skip = cfg.process.instruments2run(contains(cfg.process.instruments2run, ...
  {'FLOW','TSG','SBE45','SBE3845','ALFA','NMEA', 'PAR'}));
% Set default sync delay.
% To customize sync delay, uncomment section below
for i = 1:size(cfg.process.instruments2run)
  cfg.process.sync.delay.(cfg.process.instruments2run{i}) = 0;
end
% % Manually customize sync delay
% cfg.process.sync.delay.FLOW = 0;
% cfg.process.sync.delay.ACS412 = 0;
% cfg.process.sync.delay.BB3349 = 0;
% cfg.process.sync.delay.WS3S1081 = 0;
% cfg.process.sync.delay.NMEA = 0;
% cfg.process.sync.delay.HyperBB8002 = 0;
% cfg.process.sync.delay.LISST1183 = 0;
% cfg.process.sync.delay.SUVF6244 = 0;

%%% QC Reference (Flow Control/FLOW) %%%
cfg.process.qcref = struct();
cfg.process.qcref.reference = 'FLOW';
cfg.process.qcref.view = cfg.process.instruments2run(find(contains(cfg.process.instruments2run, ...
  {'ACS', 'AC9'}),1, 'first'));
cfg.process.qcref.mode = 'ui'; % load or ui
cfg.process.qcref.remove_old = false; % remove old selection of the same period
cfg.process.qcref.MinFiltPeriod = 50; % filter even period in minute
cfg.process.qcref.szFilt = 10; % filter even length in minute

%%% Split total/filtered %%%
cfg.process.split = struct();
cfg.process.split.reference = 'FLOW';
cfg.process.split.buffer = struct();
cfg.process.split.skip = cfg.process.instruments2run(contains(cfg.process.instruments2run, ...
  {'FLOW','TSG','SBE45','SBE3845','ALFA','NMEA', 'PAR'}));
% Set buffer length depending on instrument type (default).
% To customize buffer length, uncomment section below
for i = 1:size(cfg.process.instruments2run)
  if any(contains(cfg.process.instruments2run{i}, {'ACS', 'AC9'}))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [180, 60]; % [180, 60] for AC meters
  elseif any(contains(cfg.process.instruments2run{i}, {'BB', 'BB'}) & ~contains(cfg.process.instruments2run{i}, {'HyperBB','HBB','hbb'}))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [200, 100]; % [420, 220] for ECO-BB
  elseif any(contains(cfg.process.instruments2run{i}, 'WSCD'))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [200, 100]; % [540, 100] for ECO-fluo
  elseif any(contains(cfg.process.instruments2run{i}, 'WS3S'))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [200, 150]; % [420, 220] for ECO-fluo
  elseif any(contains(cfg.process.instruments2run{i}, 'SUVF'))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [200, 100]; % [240, 100] for Seapoint fluo
  elseif any(contains(cfg.process.instruments2run{i}, {'HyperBB', 'HBB', 'hbb'}))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [200, 100]; % [240, 140] for HyperBB
  elseif any(contains(cfg.process.instruments2run{i}, 'LISST') & ~contains(cfg.process.instruments2run{i}, {'TAU','tau','Tau','200X','200x'}))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [180, 60]; % [540, 360] for LISST
  elseif any(contains(cfg.process.instruments2run{i}, {'LISST200X', 'LISST200x'}))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [180, 60]; % [540, 360] for LISST
  elseif any(contains(cfg.process.instruments2run{i}, {'LISSTTAU','LISSTTau','LISST-TAU'}))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [180, 60]; % [180, 60] for LISST-Tau
  end
end
% % Manually customize buffer length
% cfg.process.split.buffer.ACS57 = [180, 60];
% cfg.process.split.buffer.ACS348 = [180, 60];
% cfg.process.split.buffer.LISSTTau1002G = [180, 60];
% cfg.process.split.buffer.BB31502 = [420, 220];
% cfg.process.split.buffer.WSCD859 = [540, 100];
% cfg.process.split.buffer.SUVF6244 = [240, 100]; % [660, 100]
% cfg.process.split.buffer.HyperBB8005 = [240, 140]; % [540, 340]
% cfg.process.split.buffer.LISST1183 = [540, 360];

%%% Binning %%%
cfg.process.bin = struct('bin_size', struct());
cfg.process.bin.prctile_detection = [2.5, 97.5];
% Ff StepQC with ACS: prctile_average = [2.5, 97.5]; otherwise prctile_average = [5, 75];
cfg.process.bin.prctile_average = [2.5, 97.5];
% Bin mode does not affect the outcome of the data but just the way the data is presented to the computer
% cfg.process.bin.mode = 'OneShot'; % Faster for small dataset fiting in the memory of the computer
cfg.process.bin.mode = 'ByDay'; % Slightly slower but can handle a lot more data at once as it will be binned one day at a time
% Set binning length depending on instrument type (default).
% To customize bin sizes, uncomment section below
for i = 1:size(cfg.process.instruments2run)
  if contains(cfg.process.instruments2run{i}, 'FLOW')
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for FLOW
  elseif any(contains(cfg.process.instruments2run{i}, {'GPS', 'NMEA'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for NMEA
  elseif any(contains(cfg.process.instruments2run{i}, {'TSG', 'SBE38', 'SBE45'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for TSG
  elseif any(contains(cfg.process.instruments2run{i}, {'ACS', 'AC9'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for AC meters
  elseif any(contains(cfg.process.instruments2run{i}, {'BB', 'BB'}) & ~contains(cfg.process.instruments2run{i}, {'HyperBB','HBB','hbb'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for ECO-BB
  elseif any(contains(cfg.process.instruments2run{i}, 'WSCD'))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for ECO-fluo
  elseif any(contains(cfg.process.instruments2run{i}, 'WS3S'))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for ECO-fluo
  elseif any(contains(cfg.process.instruments2run{i}, 'SUVF'))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for Seapoint fluo
  elseif any(contains(cfg.process.instruments2run{i}, {'HyperBB', 'HBB', 'hbb'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 5; % 5 min for HyperBB
  elseif any(contains(cfg.process.instruments2run{i}, 'LISST') & ~contains(cfg.process.instruments2run{i}, {'TAU','tau','Tau','200X','200x'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 10; % 10 min for LISST
  elseif any(contains(cfg.process.instruments2run{i}, 'LISST200X'))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for LISST
  elseif any(contains(cfg.process.instruments2run{i}, {'LISSTTAU','LISSTTau','LISST-TAU'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for LISST-Tau
  elseif any(contains(cfg.process.instruments2run{i}(1:3), 'ALFA'))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 10; % 10 min for ALFA
  end
end
% % Manually customize bin sizes
% cfg.process.bin.bin_size.FLOW = 1;
% cfg.process.bin.bin_size.ACS57 = 1;
% cfg.process.bin.bin_size.ACS348 = 1;
% cfg.process.bin.bin_size.LISSTTau1002G = 1;
% cfg.process.bin.bin_size.BB31502 = 1;
% cfg.process.bin.bin_size.WSCD859 = 1;
% cfg.process.bin.bin_size.SUVF6244 = 1;
% cfg.process.bin.bin_size.SBE38450091 = 1;
% cfg.process.bin.bin_size.NMEA = 1;
% cfg.process.bin.bin_size.HyperBB8005 = 5;
% cfg.process.bin.bin_size.LISST1183 = 10;
% cfg.process.bin.bin_size.SUVF6244 = 1;
% cfg.process.bin.bin_size.ALFA = 10;
% cfg.process.bin.skip = {};

%%% Automatically flagging %%%
cfg.process.flag = struct();
cfg.process.flag.skip = cfg.process.instruments2run;
% Default: parameters set to all instruments if not specific parameters set
cfg.process.flag.default = struct();
% cfg.process.flag.default.maximum_fudge_factor = 4;
% cfg.process.flag.default.variance_fudge_factor = 3;
% cfg.process.flag.default.avg_sensitivity = 1;
% cfg.process.flag.default.unc1_sensitivity = 1;
% cfg.process.flag.default.unc2_sensitivity = 2;
% cfg.process.flag.default.smooth_threshold = 60;
% cfg.process.flag.default.min_flag_n = 1;
% cfg.process.flag.default.filt = struct('smooth_threshold', 2);
  
%%% Auto QC %%%
cfg.process.qc = struct();
cfg.process.qc.RawAutoQCLim.filtered.a = 3;
cfg.process.qc.RawAutoQCLim.filtered.c = 3;
cfg.process.qc.RawAutoQCLim.total.a = 3;
cfg.process.qc.RawAutoQCLim.total.c = 3;
cfg.process.qc.RawAutoQCLim.dissolved.a = 3;
cfg.process.qc.RawAutoQCLim.dissolved.c = 3;
cfg.process.qc.RawAutoQCLim.filtered.bb = 3;
cfg.process.qc.RawAutoQCLim.total.bb = 3;
cfg.process.qc.RawAutoQCLim.dissolved.bb = 3;
cfg.process.qc.Saturation_Threshold_bb = 4100; % (counts) max being 4130

%%% Manually QC %%%
cfg.process.qc.mode = 'ui';
cfg.process.qc.qc_once_for_all = false;  % true = QC all variables | false = QC variables separately);
cfg.process.qc.remove_old = false; % remove old selection of the same period
cfg.process.qc.global = struct();
cfg.process.qc.global.active = false;
cfg.process.qc.global.view = cfg.process.qcref.view;
cfg.process.qc.global.apply = cfg.process.instruments2run(~contains(cfg.process.instruments2run, ...
  {'FLOW','NMEA', 'PAR'}));
cfg.process.qc.specific = struct();
cfg.process.qc.specific.active = true;
cfg.process.qc.specific.run = cfg.process.qcref.view;

%%% Calibrate %%%
cfg.process.calibrate = struct();
cfg.process.calibrate.skip = cfg.process.instruments2run(contains(cfg.process.instruments2run, ...
  {'FLOW','TSG','SBE45','SBE3845','ALFA','NMEA'}));
% Set calibrate options depending on instrument type (default).
for i = 1:size(cfg.process.instruments2run)
  % AC meter options
  if any(contains(cfg.process.instruments2run{i}, {'ACS', 'AC9'}))
    cfg.process.calibrate.(cfg.process.instruments2run{i}) = struct('compute_dissolved', false, ...
                                      'interpolation_method', 'linear', ... % linear CDOM
                                      'CDOM_source', 'SUVF6254', ... %
                                      'FLOW_source', 'FLOW', ...
                                      'di_method', 'normal', ... % best_di normal
                                      'compute_ad_aphi', false); % VERY SLOW: compute ad and aphi from Zheng and Stramski 2013
  % ECO-BB options
  elseif any(contains(cfg.process.instruments2run{i}, {'BB', 'BB'}) & ~contains(cfg.process.instruments2run{i}, {'HyperBB','HBB','hbb'}))
    cfg.process.calibrate.(cfg.process.instruments2run{i}) = struct('compute_dissolved', true, ...
                                      'TSG_source', 'SBE4536073', ...
                                      'FLOW_source', 'FLOW', ...
                                      'di_method', 'SW_scattering', ... % interpolate constant SW_scattering
                                      'filt_method', 'exponential_fit'); % 25percentil exponential_fit
  % HyperBB options
  elseif any(contains(cfg.process.instruments2run{i}, {'HyperBB', 'HBB', 'hbb'}))
    cfg.process.calibrate.(cfg.process.instruments2run{i}) = struct('compute_dissolved', false, ...
                                      'TSG_source', 'SBE4536073', ...
                                      'FLOW_source', 'FLOW', ...
                                      'di_method', 'SW_scattering', ... % interpolate constant SW_scattering
                                      'filt_method', 'exponential_fit'); % 25percentil exponential_fit
  % LISST100X options
  elseif any(contains(cfg.process.instruments2run{i}, 'LISST') & ~contains(cfg.process.instruments2run{i}, {'TAU','tau','Tau'}))
    cfg.process.calibrate.(cfg.process.instruments2run{i}) = struct('compute_dissolved', false, ...
                                      'FLOW_source', 'FLOW', ...
                                      'di_method', 'interpolate'); % interpolate constant
  % LISST200X options
  elseif any(contains(cfg.process.instruments2run{i}, {'LISST200X', 'LISST200x'}))
    cfg.process.calibrate.(cfg.process.instruments2run{i}) = struct('compute_dissolved', false, ...
                                      'FLOW_source', 'FLOW', ...
                                      'di_method', 'interpolate'); % interpolate constant
  % LISST-Tau options
  elseif any(contains(cfg.process.instruments2run{i}, {'LISSTTau','LISSTTAU','LISST-Tau','TAU'}))
    cfg.process.calibrate.(cfg.process.instruments2run{i}) = struct('compute_dissolved', false, ...
                                      'interpolation_method', 'linear', ... % linear CDOM
                                      'CDOM_source', 'SUVF6254', ... % 
                                      'FLOW_source', 'FLOW', ...
                                      'di_method', 'normal'); % VERY SLOW: compute ad and aphi from Zheng and Stramski 2013
  end
end
% cfg.process.calibrate.ACS57 = struct('compute_dissolved', false, ...
%                                   'interpolation_method', 'CDOM', ... % 
%                                   'CDOM_source', 'SUVF6244', ... % WSCD859
%                                   'FLOW_source', 'FLOW', ...
%                                   'di_method', 'normal', ... % best_di normal
%                                   'compute_ad_aphi', false); % VERY SLOW: compute ad and aphi from Zheng and Stramski 2013
% cfg.process.calibrate.ACS348 = struct('compute_dissolved', false, ...
%                                   'interpolation_method', 'CDOM', ... % 
%                                   'CDOM_source', 'SUVF6244', ... % WSCD859
%                                   'FLOW_source', 'FLOW', ...
%                                   'di_method', 'normal', ... % best_di normal
%                                   'compute_ad_aphi', false); % VERY SLOW: compute ad and aphi from Zheng and Stramski 2013
% cfg.process.calibrate.BB31502 = struct('compute_dissolved', true, ...
%                                   'TSG_source', 'SBE38450091', ...
%                                   'FLOW_source', 'FLOW', ...
%                                   'di_method', 'SW_scattering', ... % interpolate constant SW_scattering
%                                   'filt_method', 'exponential_fit'); % 25percentil exponential_fit
% cfg.process.calibrate.HyperBB8005 = struct('compute_dissolved', false, ...
%                                   'TSG_source', 'SBE38450091', ...
%                                   'FLOW_source', 'FLOW', ...
%                                   'di_method', 'SW_scattering', ... % interpolate constant SW_scattering
%                                   'filt_method', 'exponential_fit'); % 25percentil exponential_fit
% cfg.process.calibrate.LISST1183 = struct('compute_dissolved', false, ...
%                                   'FLOW_source', 'FLOW', ...
%                                   'di_method', 'interpolate'); % interpolate constant

%%% Write %%%
cfg.process.write = struct();
cfg.process.write.mode = 'One day one file'; % 'One file' or 'One day one file'
cfg.process.write.skip = {};







% %%% General parameters %%%
% cfg.process.days2run = datenum(2020,12,12):datenum(2021,01,05);
% cfg.process.instruments2run = {'FLOW', 'NMEA', 'ACS412', 'BB31052', 'WSCD859', ...
%   'SBE4536073', 'SUVF6254', 'LISST1183', 'HBB'}; % 'TSG', 'TSG36073', 'BB31052', 'WSCD859'
% cfg.process.write = true;
% cfg.process.force_import = false;
% cfg.process.parallel = Inf; % 0: disable parallel or Inf: as many thread available
% cfg.process.di = struct();
% cfg.process.di.skip = cfg.process.instruments2run(contains(cfg.process.instruments2run, ...
%   {'FLOW','TSG','SBE45','SBE38+45','NMEA','PAR'}));
% cfg.process.di.qc = struct('mode', 'ui',... % ui or load
%                            'qc_once_for_all', false,... % true = QC all variables | false = QC variables separately);
%                            'remove_old', false); % remove old selection of the same period
% cfg.process.di.bin = struct('bin_size', 30);
% 
% %%% Synchronization %%%
% cfg.process.sync = struct();
% cfg.process.sync.delay = struct();
% cfg.process.sync.delay.FLOW = 0;
% cfg.process.sync.delay.ACS412 = 0;
% cfg.process.sync.delay.BB31052 = 0;
% cfg.process.sync.delay.WSCD859 = 0;
% cfg.process.sync.delay.NMEA = 0;
% cfg.process.sync.delay.HyperBB8005 = 0;
% cfg.process.sync.delay.LISST1183 = 0;
% cfg.process.sync.delay.SUVF6254 = 0;
% cfg.process.sync.skip = cfg.process.instruments2run(contains(cfg.process.instruments2run, ...
%   {'FLOW','TSG','SBE45','SBE38+45','ALFA','NMEA', 'PAR'}));
% 
% %%% QC Reference (Flow Control/FLOW) %%%
% cfg.process.qcref = struct();
% cfg.process.qcref.reference = 'FLOW';
% cfg.process.qcref.view = 'ACS57';
% cfg.process.qcref.mode = 'ui'; % load or ui
% cfg.process.qcref.remove_old = false; % remove old selection of the same period
% cfg.process.qcref.MinFiltPeriod = 50; % filter even period in minute
% cfg.process.qcref.szFilt = 10; % filter even length in minute
% 
% %%% Split total/filtered %%%
% cfg.process.split = struct();
% cfg.process.split.reference = 'FLOW';
% cfg.process.split.buffer = struct();
% cfg.process.split.buffer.ACS412 = [180, 60];
% cfg.process.split.buffer.BB31052 = [420, 220];
% cfg.process.split.buffer.WSCD859 = [540, 100];
% cfg.process.split.buffer.SUVF6254 = [100, 50]; % [660, 100]
% cfg.process.split.buffer.HyperBB8005 = [240, 140]; % [540, 340]
% cfg.process.split.buffer.LISST1183 = [540, 360];
% cfg.process.split.skip = cfg.process.instruments2run(contains(cfg.process.instruments2run, ...
%   {'FLOW','TSG','SBE45','SBE38+45','ALFA','NMEA', 'PAR'}));
% 
% %%% Binning %%%
% cfg.process.bin = struct('bin_size', struct());
% cfg.process.bin.prctile_detection = [2.5, 97.5];
% % if StepQC with ACS: prctile_average = [2.5, 97.5]; otherwise prctile_average = [5, 75];
% cfg.process.bin.prctile_average = [2.5, 97.5];
% % Bin mode does not affect the outcome of the data but just the way the data is presented to the computer
% % cfg.process.bin.mode = 'OneShot'; % Faster for small dataset fiting in the memory of the computer
% cfg.process.bin.mode = 'ByDay'; % Slightly slower but can handle a lot more data at once as it will be binned one day at a time
% cfg.process.bin.bin_size.FLOW = 1;
% cfg.process.bin.bin_size.ACS412 = 1;
% cfg.process.bin.bin_size.BB31052 = 1;
% cfg.process.bin.bin_size.WSCD859 = 1;
% cfg.process.bin.bin_size.SUVF6254 = 1;
% cfg.process.bin.bin_size.SBE4536073 = 1;
% cfg.process.bin.bin_size.NMEA = 1;
% cfg.process.bin.bin_size.HyperBB8005 = 5;
% cfg.process.bin.bin_size.LISST1183 = 10;
% cfg.process.bin.bin_size.ALFA = 10;
% cfg.process.bin.skip = {};
% 
% %%% Automatically flagging %%%
% cfg.process.flag = struct();
% cfg.process.flag.skip = {'FLOW', 'TSG', 'NMEA', 'HBB', 'ACS57', 'BB31052', ...
%   'WSCD859', 'AC9', 'ACS007', 'ACS348', 'ACS091', 'ACS111', 'ACS279', 'LISST1183', ...
%   'WSCD1082P', 'ALFA','PAR','ACS298'};
% cfg.process.flag.skip = [cfg.process.flag.skip cfg.process.instruments2run];
% % Default: parameters set to all instruments if not specific parameters set
% cfg.process.flag.default = struct();
% cfg.process.flag.default.maximum_fudge_factor = 4;
% cfg.process.flag.default.variance_fudge_factor = 3;
% cfg.process.flag.default.avg_sensitivity = 1;
% cfg.process.flag.default.unc1_sensitivity = 1;
% cfg.process.flag.default.unc2_sensitivity = 2;
% cfg.process.flag.default.smooth_threshold = 60;
% cfg.process.flag.default.min_flag_n = 1;
% cfg.process.flag.default.filt = struct('smooth_threshold', 2);
%   
% %%% Auto QC %%%
% cfg.process.qc = struct();
% cfg.process.qc.RawAutoQCLim.filtered.a = 3;
% cfg.process.qc.RawAutoQCLim.filtered.c = 3;
% cfg.process.qc.RawAutoQCLim.total.a = 3;
% cfg.process.qc.RawAutoQCLim.total.c = 3;
% cfg.process.qc.RawAutoQCLim.dissolved.a = 3;
% cfg.process.qc.RawAutoQCLim.dissolved.c = 3;
% cfg.process.qc.RawAutoQCLim.filtered.bb = 3;
% cfg.process.qc.RawAutoQCLim.total.bb = 3;
% cfg.process.qc.RawAutoQCLim.dissolved.bb = 3;
% cfg.process.qc.Saturation_Threshold_bb = 4100; % (counts) max being 4130
% 
% %%% Manually QC %%%
% cfg.process.qc.mode = 'ui';
% cfg.process.qc.qc_once_for_all = false;  % true = QC all variables | false = QC variables separately);
% cfg.process.qc.remove_old = false; % remove old selection of the same period
% cfg.process.qc.global = struct();
% cfg.process.qc.global.active = false;
% cfg.process.qc.global.view = cfg.process.qcref.view;
% cfg.process.qc.global.apply = cfg.process.instruments2run(contains(cfg.process.instruments2run, ...
%   {'FLOW','NMEA', 'PAR'}));
% cfg.process.qc.specific = struct();
% cfg.process.qc.specific.active = true;
% cfg.process.qc.specific.run = cfg.process.qcref.view;
% 
% %%% Calibrate %%%
% cfg.process.calibrate = struct();
% cfg.process.calibrate.ACS298 = struct('compute_dissolved', false, ...
%                                   'interpolation_method', 'CDOM', ... % 
%                                   'CDOM_source', 'SUVF6254', ... % WSCD859
%                                   'FLOW_source', 'FLOW', ...
%                                   'di_method', 'normal', ... % best_di normal
%                                   'compute_ad_aphi', false); % VERY SLOW: compute ad and aphi from Zheng and Stramski 2013
% cfg.process.calibrate.BB31052 = struct('compute_dissolved', false, ...
%                                   'TSG_source', 'SBE4536073', ...
%                                   'FLOW_source', 'FLOW', ...
%                                   'di_method', 'constant', ... % interpolate constant SW_scattering
%                                   'filt_method', 'exponential_fit'); % 25percentil exponential_fit
% cfg.process.calibrate.HyperBB8005 = struct('compute_dissolved', false, ...
%                                   'TSG_source', 'SBE4536073', ...
%                                   'FLOW_source', 'FLOW', ...
%                                   'di_method', 'constant', ... % interpolate constant SW_scattering
%                                   'filt_method', 'exponential_fit'); % 25percentil exponential_fit
% cfg.process.calibrate.LISST1183 = struct('compute_dissolved', false, ...
%                                   'FLOW_source', 'FLOW', ...
%                                   'di_method', 'interpolate'); % interpolate constant
% cfg.process.calibrate.skip = cfg.process.instruments2run(contains(cfg.process.instruments2run, ...
%   {'FLOW','TSG','SBE45','SBE38+45','ALFA','NMEA'}));
% 
% %%% Write %%%
% cfg.process.write = struct();
% cfg.process.write.mode = 'One day one file'; % 'One file' or 'One day one file'
% cfg.process.write.skip = {};

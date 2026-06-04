% SOPACE-TN444 Configuration file
% author: Guilaume Bourdin
% created: May 2025

cfg = struct('meta', struct(), 'instruments', struct(), 'process', struct());

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                   METADATA (ENTER YOUR METADATA HERE)                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg.meta.investigators = 'Alison_Chase,Guillaume_Bourdin,Charlotte_Begouen_Demeaux';
cfg.meta.affiliations = 'University_of_Washington';
cfg.meta.emails = 'alichase@uw.edu';
cfg.meta.experiment = 'TN444SOPACE';
cfg.meta.cruise = 'TN444SOPACE';
cfg.meta.station = 'NA';
cfg.meta.documents = 'NA';
cfg.meta.calibration_files = 'NA';
cfg.meta.data_type = 'flow_thru';
cfg.meta.data_status = 'preliminary';
cfg.meta.measurement_depth = 2;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                INSTRUMENTS (SETUP YOUR INSTRUMENTS HERE)                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PATH_ROOT = '/Volumes/KM24-SOPACE2';
PATH_ROOT = '/Volumes/Chase-SSD1/TN444-SOPACE';

%%% TSG %%%
model = 'SBE3845';
SN = 'TN444';
cfg.instruments.([model SN]) = struct();
cfg.instruments.([model SN]).model = model;
cfg.instruments.([model SN]).TSG_source = true;
cfg.instruments.([model SN]).boat = 'Tara';
cfg.instruments.([model SN]).logger = 'Inlinino_base'; % TeraTerm Matlab Inlinino
cfg.instruments.([model SN]).sn = SN;
cfg.instruments.([model SN]).path = struct('raw', fullfile(PATH_ROOT, 'raw', [model SN]),...
                                  'wk', fullfile(PATH_ROOT, 'wk', [model SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', [model SN]));
cfg.instruments.([model SN]).view = struct('varname', 'sss');
cfg.instruments.([model SN]).temperature_variable = 'sst';
cfg.instruments.([model SN]).salinity_variable = 'sss';

% %%% TSG %%%
% model = 'SBE3845';
% SN = '04970269';
% cfg.instruments.([model SN]) = struct();
% cfg.instruments.([model SN]).model = model;
% cfg.instruments.([model SN]).TSG_source = true;
% cfg.instruments.([model SN]).boat = 'Tara';
% cfg.instruments.([model SN]).logger = 'Inlinino_base'; % TeraTerm Matlab Inlinino
% cfg.instruments.([model SN]).sn = SN;
% cfg.instruments.([model SN]).path = struct('raw', fullfile(PATH_ROOT, 'raw', [model SN]),...
%                                   'wk', fullfile(PATH_ROOT, 'wk', [model SN]),...
%                                   'prod', fullfile(PATH_ROOT, 'prod'),...
%                                   'ui', fullfile(PATH_ROOT, 'ui', [model SN]));
% cfg.instruments.([model SN]).view = struct('varname', 't2');
% cfg.instruments.([model SN]).temperature_variable = 't2';

%%% NMEA %%%
model = 'GPS'; % GPSSC701Tara GPS32Tara GPSCOMPASSAT
SN = 'TN444';
cfg.instruments.NMEA = struct();
cfg.instruments.NMEA.model = model;
% cfg.instruments.NMEA.boat = 'Tara';
cfg.instruments.NMEA.logger = 'Inlinino';
cfg.instruments.NMEA.sn = SN;
cfg.instruments.NMEA.prefix = [model SN];
cfg.instruments.NMEA.path = struct('raw', fullfile(PATH_ROOT, 'raw', [model SN]),...
                                  'wk', fullfile(PATH_ROOT, 'wk', [model SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', [model SN]));
cfg.instruments.NMEA.view = struct('varname', 'lat');

%%% FLOW (FlowControl) %%%
SN = 'A54704';
model='ADU200';
cfg.instruments.FLOW = struct();
cfg.instruments.FLOW.model = model;
cfg.instruments.FLOW.logger = 'Inlinino_base';% 'FlowControl';
cfg.instruments.FLOW.sn = SN;
cfg.instruments.FLOW.LoadPrevious = true;
cfg.instruments.FLOW.analog1 = '';
cfg.instruments.FLOW.analog2 = '';
cfg.instruments.FLOW.path = struct('raw', fullfile(PATH_ROOT, 'raw', [model SN]),...
                                  'wk', fullfile(PATH_ROOT, 'wk', [model SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', [model SN]));
cfg.instruments.FLOW.view = struct('varname','spd1','swt_variable','swt1','spd_variable','spd1'); % swt1 swt2 spd1

%%% ACS111 %%%
SN = '111';
cfg.instruments.(['ACS' SN]) = struct();
cfg.instruments.(['ACS' SN]).di = struct();
cfg.instruments.(['ACS' SN]).di.prefix = ['DIW_ACS' SN '_'];
cfg.instruments.(['ACS' SN]).di.postfix = '';
cfg.instruments.(['ACS' SN]).model = 'ACS';
cfg.instruments.(['ACS' SN]).AC_source = true;
cfg.instruments.(['ACS' SN]).sn = SN;
cfg.instruments.(['ACS' SN]).logger = 'InlininoACScsv';
cfg.instruments.(['ACS' SN]).device_file = fullfile(PATH_ROOT, 'DeviceFiles', 'ACS111_20231006.dev');
cfg.instruments.(['ACS' SN]).path = struct('raw', fullfile(PATH_ROOT, 'raw', ['ACS' SN]),...
                                    'di', fullfile(PATH_ROOT, 'raw', ['ACS' SN]),...
                                    'wk', fullfile(PATH_ROOT, 'wk', ['ACS' SN]),...
                                    'prod', fullfile(PATH_ROOT, 'prod'),...
                                    'ui', fullfile(PATH_ROOT, 'ui', ['ACS' SN]));
cfg.instruments.(['ACS' SN]).view = struct('varname', 'a', 'varcol', 40);

% %%% HBB %%%
SN = '8005';
cfg.instruments.(['HyperBB' SN]) = struct();
cfg.instruments.(['HyperBB' SN]).di = struct();
cfg.instruments.(['HyperBB' SN]).di.prefix = ['DIW_HyperBB' SN '_'];
cfg.instruments.(['HyperBB' SN]).di.postfix = '';
cfg.instruments.(['HyperBB' SN]).model = 'HBB';
cfg.instruments.(['HyperBB' SN]).sn = SN;
cfg.instruments.(['HyperBB' SN]).k_exp = 0.1058;
cfg.instruments.(['HyperBB' SN]).ila_prefix = ['HyperBB' SN];
cfg.instruments.(['HyperBB' SN]).logger = 'InlininoHBB';
cfg.instruments.(['HyperBB' SN]).theta = 135;
cfg.instruments.(['HyperBB' SN]).fdom_ag_correlation = fullfile(PATH_ROOT, 'DeviceFiles', 'parameter_to_convert_FDOM_to_ag_TaraEuropa_2023_2024.csv');
cfg.instruments.(['HyperBB' SN]).PlaqueCal = fullfile(PATH_ROOT, 'DeviceFiles', 'Hbb_Cal_Plaque_20241010_142447.mat');
cfg.instruments.(['HyperBB' SN]).TemperatureCal = fullfile(PATH_ROOT, 'DeviceFiles', 'Hbb_Cal_Temp_20241008_130132.mat');
cfg.instruments.(['HyperBB' SN]).path = struct('raw', fullfile(PATH_ROOT, 'raw', ['HyperBB' SN]),...
                                  'di', fullfile(PATH_ROOT, 'raw', ['HyperBB' SN]),...
                                  'wk', fullfile(PATH_ROOT, 'wk', ['HyperBB' SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', ['HyperBB' SN]));
cfg.instruments.(['HyperBB' SN]).view = struct('varname', 'beta', 'varcol', 23);


%%% SUVF %%%
SN = '6266';
interface = 'DataQDI2108999';
cfg.instruments.(['SUVF' SN]) = struct();
cfg.instruments.(['SUVF' SN]).model = 'CD';
cfg.instruments.(['SUVF' SN]).CDOM_source = true;
cfg.instruments.(['SUVF' SN]).sn = SN;
cfg.instruments.(['SUVF' SN]).slope = 1;
cfg.instruments.(['SUVF' SN]).dark = 0.076;
cfg.instruments.(['SUVF' SN]).logger = 'Inlinino_base';
cfg.instruments.(['SUVF' SN]).analog_channel = ['C5_SUVF' SN];
if contains(interface, 'ADU100') || contains(interface, 'DataQ')
  cfg.instruments.(['SUVF' SN]).ila_prefix = interface;
  raw_path = fullfile(PATH_ROOT, 'raw', interface);
else
  cfg.instruments.(['SUVF' SN]).ila_prefix = ['SUVF' SN];
  raw_path = ['SUVF' SN];
end
cfg.instruments.(['SUVF' SN]).di = struct();
cfg.instruments.(['SUVF' SN]).di.prefix = ['DIW_' cfg.instruments.(['SUVF' SN]).ila_prefix SN '_'];
cfg.instruments.(['SUVF' SN]).di.postfix = '';
cfg.instruments.(['SUVF' SN]).path = struct('raw', raw_path,...
                                  'di', raw_path,...
                                  'wk', fullfile(PATH_ROOT, 'wk', ['SUVF' SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', ['SUVF' SN]));
cfg.instruments.(['SUVF' SN]).view = struct('varname', 'fdom', 'varcol', 1);

%%% WS3S %%%
SN = '1081P';
interface = 'DataQDI2108999';
cfg.instruments.(['WS3S' SN]) = struct();
cfg.instruments.(['WS3S' SN]).model = 'FL';
cfg.instruments.(['WS3S' SN]).sn = SN;
cfg.instruments.(['WS3S' SN]).slope = 1;
cfg.instruments.(['WS3S' SN]).dark = 0.076;
cfg.instruments.(['WS3S' SN]).logger = 'Inlinino_base';
cfg.instruments.(['WS3S' SN]).analog_channel = ['C1_WS3S' SN];
if contains(interface, 'ADU100') || contains(interface, 'DataQ')
  cfg.instruments.(['WS3S' SN]).ila_prefix = [interface '_'];
  raw_path = fullfile(PATH_ROOT, 'raw', interface);
else
  cfg.instruments.(['WS3S' SN]).ila_prefix = ['WS3S' SN];
  raw_path = ['WS3S' SN];
end
cfg.instruments.(['WS3S' SN]).di = struct();
cfg.instruments.(['WS3S' SN]).di.prefix = ['DIW_' cfg.instruments.(['WS3S' SN]).ila_prefix SN '_'];
cfg.instruments.(['WS3S' SN]).di.postfix = '';
cfg.instruments.(['WS3S' SN]).path = struct('raw',  raw_path,...
                                  'di',  raw_path,...
                                  'wk', fullfile(PATH_ROOT, 'wk', ['WS3S' SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', ['WS3S' SN]));
cfg.instruments.(['WS3S' SN]).view = struct('varname', 'fdom', 'varcol', 1);

%%% LISST100X %%%
model = 'LISST100X'; % 100X 200X
SN = '1183';
cfg.instruments.([model SN]) = struct();
cfg.instruments.([model SN]).di = struct();
cfg.instruments.([model SN]).di.prefix = ['DIW_LISST' SN '_'];
cfg.instruments.([model SN]).di.postfix = '';
cfg.instruments.([model SN]).model = model;
cfg.instruments.([model SN]).type = 'B'; 
cfg.instruments.([model SN]).sn = SN;
cfg.instruments.([model SN]).logger = 'InlininoLISSTcsv';
cfg.instruments.([model SN]).zsc = fullfile(PATH_ROOT, 'DeviceFiles', 'LISST100X1183_20230215_factory_zsc.asc');
cfg.instruments.([model SN]).dcal = fullfile(PATH_ROOT, 'DeviceFiles', 'LISST100X1183_20230215_Ringarea.asc');
cfg.instruments.([model SN]).vcc = fullfile(PATH_ROOT, 'DeviceFiles', 'LISST100X1183_20230215_InstrumentData.txt');
% cfg.instruments.([model SN]).zsc = [2.203500e+001, 2.568500e+001, 2.503000e+001, 2.986000e+001, 2.842500e+001, 3.283000e+001, 3.077000e+001, 3.659500e+001, 2.978000e+001, 3.552000e+001, 3.198000e+001, 4.216000e+001, 3.916500e+001, 4.662500e+001, 3.974000e+001, 4.454000e+001, 4.403500e+001, 4.604500e+001, 4.430000e+001, 4.510500e+001, 4.719500e+001, 3.850000e+001, 5.373000e+001, 2.664000e+001, 3.180500e+001, 1.655500e+001, 2.205500e+001, 1.554000e+001, 1.422000e+001, 1.123000e+001, 8.780000e+000, 8.555000e+000, 1.515000e+003, 1.167900e+003, 6.410000e+001, 1.055150e+003, 7.700000e+001, 2.116600e+003, 1.807000e+003, 5.476500e+003];
% % Original dcal file (ring area)
% % cfg.instruments.([model SN]).dcal = [1.0000000e+000, 1.0038000e+000, 9.9360000e-001, 1.0027000e+000, 9.9720000e-001, 9.9570000e-001, 9.9030000e-001, 9.9430000e-001, 9.9290000e-001, 9.9000000e-001, 9.9290000e-001, 9.9300000e-001, 9.9150000e-001, 9.9300000e-001, 9.9230000e-001, 9.9090000e-001, 1.1032000e+000, 1.1123000e+000, 1.2430000e+000, 1.1562000e+000, 1.3273000e+000, 1.1999000e+000, 1.0740000e+000, 1.7489000e+000, 1.5382000e+000, 2.5109000e+000, 2.5468000e+000, 3.5504000e+000, 3.9338000e+000, 5.1747342e+000, 7.5143548e+000, 1.2528083e+001];
% % dcal adhoc from email of Wayne Slade Dec 8, 2017
% cfg.instruments.([model SN]).dcal = [ 1.0179083e+00	   9.9213489e-01	   1.0108161e+00	   9.9492883e-01	   1.0043707e+00	   9.9891840e-01	   9.9859055e-01	   1.0042049e+00	   1.0000763e+00	   9.9889997e-01	   1.0009497e+00	   1.0004019e+00	   1.0011130e+00	   1.0004677e+00	   1.0213554e+00	   9.9990262e-01	   1.1115630e+00	   1.1206668e+00	   1.2493699e+00	   1.1643199e+00	   1.3355657e+00	   1.2090892e+00	   1.0781540e+00	   1.7620752e+00	   1.5508563e+00	   2.5304119e+00	   2.5638592e+00	   3.5757212e+00	   3.9631987e+00	   5.0166411e+00	   5.6381118e+00	   8.6881539e+00];
% % customize dcal of 1183 on Jan 16, 2018 due to bump in VSF at ring 30
% cfg.instruments.([model SN]).dcal(30) = 2.6; % Good for low values but bad for high values ??
% cfg.instruments.([model SN]).vcc = 47447; % 2018 cal: 48493 | 2023 cal 47447
cfg.instruments.([model SN]).inversion = 'spherical';
% cfg.instruments.([model SN]).ds = [1.2500,1.4750,1.7405,2.0538,2.4235,2.8597,3.3744,3.9818,4.6986,5.5443,6.5423,7.7199,9.1095,10.7492,12.6841,14.9672,17.6613,20.8403,24.5916,29.0180,34.2413,40.4047,47.6776,56.2595,66.3863,78.3358,92.4362,109.0747,128.7082,151.8757,179.2133,211.4717,249.5366];
% cfg.instruments.([model SN]).theta = [0.082, 0.096, 0.114, 0.134, 0.158, 0.187, 0.221, 0.260, 0.307, 0.362, 0.428, 0.505, 0.596, 0.703, 0.829, 0.979, 1.155, 1.363, 1.609, 1.898, 2.240, 2.643, 3.119, 3.681, 4.344, 5.126, 6.049, 7.138, 8.424, 9.941, 11.73, 13.84];
cfg.instruments.([model SN]).path = struct('raw', fullfile(PATH_ROOT, 'raw', [model SN]),...
                                  'di', fullfile(PATH_ROOT, 'raw', [model SN]),...
                                  'wk', fullfile(PATH_ROOT, 'wk', [model SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', [model SN]));
cfg.instruments.([model SN]).view = struct('varname', 'beta', 'varcol', 15);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%              PROCESS (NO CHANGE REQUIRED BEYOND THIS POINT)             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% General parameters %%%
cfg.process.days2run = datetime(2020,12,12):datetime(2022,10,15);
% cfg.process.instruments2run = {'FLOW', 'NMEA', 'ACS57', 'BB31502', 'WSCD859', ...
%   'SBE4536073', 'SUVF6244', 'LISST1183', 'HyperBB8005'};
cfg.process.instruments2run = fieldnames(cfg.instruments);
cfg.process.write = true;
cfg.process.force_import = false;
cfg.process.parallel = Inf; % 0: disable parallel or Inf: as many thread available
cfg.process.di = struct();
cfg.process.di.skip = cfg.process.instruments2run(contains(lower(cfg.process.instruments2run), ...
  {'flow','tsg','sbe45','sbe3845','nmea','par'}));
cfg.process.di.qc = struct('mode', 'ui',... % ui or load
                           'qc_once_for_all', false,... % true = QC all variables | false = QC variables separately);
                           'remove_old', false,... % remove old selection of the same period
                           'remove_when_flow_below', false); % true = remove data when flow <= 0.5 | false = no data removal data depending on flow | number = remove data when flow <= number)
cfg.process.di.bin = struct('bin_size', 30);

%%% Synchronization %%%
cfg.process.sync = struct();
cfg.process.sync.delay = struct();
cfg.process.sync.skip = cfg.process.instruments2run(contains(lower(cfg.process.instruments2run), ...
  {'flow','tsg','sbe45','sbe3845','nmea','par', 'alfa'}));
% Set default sync delay.
% To customize sync delay, uncomment section below
for i = 1:length(cfg.process.instruments2run)
  cfg.process.sync.delay.(cfg.process.instruments2run{i}) = seconds(0);
end
% % Manually customize sync delay
% cfg.process.sync.delay.FLOW = seconds(0);
% cfg.process.sync.delay.ACS412 = seconds(0);
% cfg.process.sync.delay.BB3349 = seconds(0);
% cfg.process.sync.delay.WS3S1081 = seconds(0);
% cfg.process.sync.delay.NMEA = seconds(0);
% cfg.process.sync.delay.HyperBB8002 = seconds(0);
% cfg.process.sync.delay.LISST1183 = seconds(0);
% cfg.process.sync.delay.SUVF6244 = seconds(0);

%%% QC Reference (Flow Control/FLOW) %%%
cfg.process.qcref = struct();
cfg.process.qcref.reference = 'FLOW';
cfg.process.qcref.view = cfg.process.instruments2run{find(contains(lower(cfg.process.instruments2run), ...
  {'acs', 'ac9'}),1, 'first')};
cfg.process.qcref.mode = 'ui'; % load or ui
cfg.process.qcref.remove_old = false; % remove old selection of the same period
cfg.process.qcref.MinFiltPeriod = minutes(50); % filter even period in minute
cfg.process.qcref.szFilt = minutes(10); % filter even length in minute

%%% Split total/filtered %%%
cfg.process.split = struct();
cfg.process.split.reference = 'FLOW';
cfg.process.split.buffer = struct();
cfg.process.split.skip = cfg.process.instruments2run(contains(lower(cfg.process.instruments2run), ...
  {'flow','tsg','sbe45','sbe3845','nmea','par', 'alfa'}));
% Set buffer length depending on instrument type (default).
% To customize buffer length, uncomment section below
for i = 1:length(cfg.process.instruments2run)
  if any(contains(lower(cfg.process.instruments2run{i}), {'acs', 'ac9'}))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = seconds([180, 60]); % [180, 60] for AC meters
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'bb') & ~contains(lower(cfg.process.instruments2run{i}), {'hyperbb','hbb'}))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = seconds([200, 100]); % [420, 220] for ECO-BB
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'wscd'))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = seconds([200, 100]); % [540, 100] for ECO-fluo
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'ws3s'))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = seconds([200, 150]); % [420, 220] for ECO-fluo
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'suvf'))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = seconds([180, 60]); % [240, 100] for Seapoint fluo
  elseif any(contains(lower(cfg.process.instruments2run{i}), {'hyperbb','hbb'}))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = seconds([200, 100]); % [240, 140] for HyperBB
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'lisst100x'))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = seconds([180, 60]); % [540, 360] for LISST
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'lisst200x'))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = seconds([180, 60]); % [540, 360] for LISST
  elseif any(contains(lower(cfg.process.instruments2run{i}), {'lissttau','lisst-tau'}))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = seconds([180, 60]); % [180, 60] for LISST-Tau
  else
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = seconds([180, 60]);
  end
end
% % Manually customize buffer length
% cfg.process.split.buffer.ACS57 = seconds([180, 60]);
% cfg.process.split.buffer.ACS348 = seconds([180, 60]);
% cfg.process.split.buffer.LISSTTau1002G = seconds([180, 60]);
% cfg.process.split.buffer.BB31502 = seconds([420, 220]);
% cfg.process.split.buffer.WSCD859 = seconds([540, 100]);
% cfg.process.split.buffer.SUVF6244 = seconds([180, 60]); % [660, 100]
% cfg.process.split.buffer.HyperBB8005 = seconds([240, 140]); % [540, 340]
% cfg.process.split.buffer.LISST1183 = seconds([540, 360]);

%%% Binning %%%
cfg.process.bin = struct('bin_size', struct());
cfg.process.bin.skip = {};
cfg.process.bin.prctile_detection = [2.5, 97.5];
% StepQC with ACS: prctile_average = [2.5, 97.5]; otherwise prctile_average = [5, 75];
cfg.process.bin.prctile_average = [2.5, 97.5];
% Bin mode does not affect the outcome of the data but just the way the data is presented to the computer
% cfg.process.bin.mode = 'OneShot'; % Faster for small dataset fiting in the memory of the computer
cfg.process.bin.mode = 'ByDay'; % Slightly slower but can handle a lot more data at once as it will be binned one day at a time
% Set binning length depending on instrument type (default).
% To customize bin sizes, uncomment section below
for i = 1:length(cfg.process.instruments2run)
  if contains(lower(cfg.process.instruments2run{i}), 'flow')
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = minutes(1); % 1 min for FLOW
  elseif any(contains(lower(cfg.process.instruments2run{i}), {'gps', 'nmea'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = minutes(1); % 1 min for NMEA
  elseif any(contains(lower(cfg.process.instruments2run{i}), {'tsg', 'sbe38', 'sbe45'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = minutes(1); % 1 min for TSG
  elseif any(contains(lower(cfg.process.instruments2run{i}), {'acs', 'ac9'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = minutes(1); % 1 min for AC meters
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'bb') & ~contains(lower(cfg.process.instruments2run{i}), {'hyperbb','hbb'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = minutes(1); % 1 min for ECO-BB
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'wscd'))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = minutes(1); % 1 min for ECO-fluo
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'ws3s'))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = minutes(1); % 1 min for ECO-fluo
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'suvf'))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = minutes(1); % 1 min for Seapoint fluo
  elseif any(contains(lower(cfg.process.instruments2run{i}), {'hyperbb','hbb'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = minutes(5); % 5 min for HyperBB
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'lisst100x'))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = minutes(2); % 2 min for LISST100X
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'lisst200x'))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = minutes(2); % 1 min for LISST200X
  elseif any(contains(lower(cfg.process.instruments2run{i}), {'lissttau','lisst-tau'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = minutes(1); % 1 min for LISST-Tau
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'alfa'))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = minutes(10); % 10 min for ALFA
  else
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = minutes(1);
  end
end
% % Manually customize bin sizes
% cfg.process.bin.bin_size.FLOW = minutes(1);
% cfg.process.bin.bin_size.ACS57 = minutes(1);
% cfg.process.bin.bin_size.ACS348 = minutes(1);
% cfg.process.bin.bin_size.LISSTTau1002G = minutes(1);
% cfg.process.bin.bin_size.BB31502 = minutes(1);
% cfg.process.bin.bin_size.WSCD859 = minutes(1);
% cfg.process.bin.bin_size.SUVF6244 = minutes(1);
% cfg.process.bin.bin_size.SBE38450091 = minutes(1);
% cfg.process.bin.bin_size.NMEA = minutes(1);
% cfg.process.bin.bin_size.HyperBB8005 = minutes(5);
% cfg.process.bin.bin_size.LISST1183 = minutes(2);
% cfg.process.bin.bin_size.SUVF6244 = minutes(1);
% cfg.process.bin.bin_size.ALFA = minutes(10);
% cfg.process.bin.skip = {};

%%% Pass data to QC level %%% (FLag function deprecated
cfg.process.Pass2QC = struct();
cfg.process.Pass2QC.skip = cfg.process.instruments2run;
% Default: parameters set to all instruments if not specific parameters set
% cfg.process.flag.default = struct();
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
cfg.process.qc.AutoQC_tolerance.filtered.a = 3;
cfg.process.qc.AutoQC_tolerance.filtered.c = 3;
cfg.process.qc.AutoQC_tolerance.total.a = 3;
cfg.process.qc.AutoQC_tolerance.total.c = 3;
cfg.process.qc.AutoQC_tolerance.dissolved.a = 3;
cfg.process.qc.AutoQC_tolerance.dissolved.c = 3;
cfg.process.qc.AutoQC_tolerance.filtered.bb = 3;
cfg.process.qc.AutoQC_Saturation_Threshold.a = 50; % in uncalibrated m^-1
cfg.process.qc.AutoQC_Saturation_Threshold.c = 50; % in uncalibrated m^-1
cfg.process.qc.AutoQC_tolerance.total.bb = 3;
cfg.process.qc.AutoQC_tolerance.dissolved.bb = 3;
cfg.process.qc.AutoQC_Saturation_Threshold.bb = 4100; % (counts) max being 4130

%%% Manually QC %%%
cfg.process.qc.mode = 'ui';
cfg.process.qc.qc_once_for_all = false;  % true = QC all variables | false = QC variables separately);
cfg.process.qc.remove_old = false; % remove old selection of the same period
cfg.process.qc.remove_when_flow_below = false; % true = remove data when flow <= 0.5 | false = no data removal data depending on flow | number = remove data when flow <= number)
cfg.process.qc.global = struct();
cfg.process.qc.global.active = false;
cfg.process.qc.global.view = cfg.process.qcref.view;
cfg.process.qc.global.apply = cfg.process.instruments2run(~contains(lower(cfg.process.instruments2run), ...
  {'flow','nmea', 'par'}));
cfg.process.qc.specific = struct();
cfg.process.qc.specific.active = true;
cfg.process.qc.specific.run = {cfg.process.qcref.view};

%%% Calibrate %%%
cfg.process.calibrate = struct();
cfg.process.calibrate.skip = cfg.process.instruments2run(contains(lower(cfg.process.instruments2run), ...
  {'flow','tsg','sbe45','sbe3845','nmea','alfa'}));
% cfg.process.min_nb_pts_per_cluster = 200;
% look for TSG, SUVF, and C, and automatically turn off source not available
cfg.process.TSG_source = '';
cfg.process.CDOM_source = '';
cfg.process.AC_source = '';
for i = fieldnames(cfg.instruments)'
  if isfield(cfg.instruments.(i{:}), 'TSG_source')
    if cfg.instruments.(i{:}).TSG_source
      cfg.process.TSG_source = i{:};
    end
  end
  if isfield(cfg.instruments.(i{:}), 'CDOM_source')
    if cfg.instruments.(i{:}).CDOM_source
      cfg.process.CDOM_source = i{:};
    end
  end
  if isfield(cfg.instruments.(i{:}), 'AC_source')
    if cfg.instruments.(i{:}).AC_source
      cfg.process.AC_source = i{:};
    end
  end
end
% if no TSG, CDOM, and C source indicated in cfg: find TSG, CDOM, and C instruments automatically
if isempty(cfg.process.TSG_source)
  if any(contains(lower(cfg.process.instruments2run), {'sbe3845', 'sbe45', 'atlasecrtd'}))
    cfg.process.TSG_source = cfg.process.instruments2run{find(contains(lower(cfg.process.instruments2run), {'sbe3845', 'sbe45', 'atlasecrtd'}), 1, 'first')};
  end
end
if isempty(cfg.process.CDOM_source)
  if any(contains(lower(cfg.process.instruments2run), {'wscd', 'suvf'}))
    cfg.process.CDOM_source = cfg.process.instruments2run{find(contains(lower(cfg.process.instruments2run), {'wscd', 'suvf'}), 1, 'first')};
  end
end
if isempty(cfg.process.AC_source)
  if any(contains(lower(cfg.process.instruments2run), {'acs', 'ac9'}))
    cfg.process.AC_source = cfg.process.instruments2run{find(contains(lower(cfg.process.instruments2run), {'acs', 'ac9'}), 1, 'first')};
  end
end

% Set calibrate options depending on instrument type (default).
for i = 1:length(cfg.process.instruments2run)
  % AC meter options
  if any(contains(lower(cfg.process.instruments2run{i}), {'acs', 'ac9'}))
    cfg.process.calibrate.(cfg.process.instruments2run{i}) = struct('compute_dissolved', false, ...
                                      'TSG_source', cfg.process.TSG_source, ...
                                      'interpolation_method', 'CDOM', ... % choose one: linear CDOM
                                      'CDOM_source', cfg.process.CDOM_source, ... 
                                      'FLOW_source', 'FLOW', ...
                                      'di_method', 'best_di', ... % best_di normal
                                      'scattering_correction', 'Semiempirical_blended2', ... % Zaneveld1994_proportional Rottgers2013_semiempirical Semiempirical_blended1 Semiempirical_blended2 Semiempirical_blended3
                                      'compute_ad_aphi', false); % VERY SLOW: compute ad and aphi from Zheng and Stramski 2013
  % ECO-BB options
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'bb') & ~contains(lower(cfg.process.instruments2run{i}), {'hyperbb', 'hbb'}))
    cfg.process.calibrate.(cfg.process.instruments2run{i}) = struct('compute_dissolved', true, ...
                                      'TSG_source', cfg.process.TSG_source, ...
                                      'FLOW_source', 'FLOW', ...
                                      'AC_source', cfg.process.AC_source, ...
                                      'CDOM_source', cfg.process.CDOM_source, ... 
                                      'di_method', 'SW_scattering', ... % interpolate constant SW_scattering
                                      'filt_method', 'exponential_fit'); % 25percentil exponential_fit
  % ECO-FL options
  elseif any(contains(lower(cfg.process.instruments2run{i}), {'ws3s'}))
    cfg.process.calibrate.(cfg.process.instruments2run{i}) = struct('compute_dissolved', true, ...
                                      'FLOW_source', 'FLOW', ...
                                      'di_method', 'best_di', ... % best_di interpolate constant SW_scattering
                                      'filt_method', 'exponential_fit'); % 25percentil exponential_fit
  % HyperBB options
  elseif any(contains(lower(cfg.process.instruments2run{i}), {'hyperbb', 'hbb'}))
    cfg.process.calibrate.(cfg.process.instruments2run{i}) = struct('compute_dissolved', false, ...
                                      'TSG_source', cfg.process.TSG_source, ...
                                      'FLOW_source', 'FLOW', ...
                                      'AC_source', cfg.process.AC_source, ...
                                      'CDOM_source', cfg.process.CDOM_source, ... 
                                      'di_method', 'SW_scattering', ... % interpolate constant SW_scattering
                                      'filt_method', 'exponential_fit'); % 25percentil exponential_fit
  % LISST100X options
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'lisst100x'))
    cfg.process.calibrate.(cfg.process.instruments2run{i}) = struct('compute_dissolved', false, ...
                                      'FLOW_source', 'FLOW', ...
                                      'di_method', 'interpolate'); % interpolate constant
  % LISST200X options
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'lisst200x'))
    cfg.process.calibrate.(cfg.process.instruments2run{i}) = struct('compute_dissolved', false, ...
                                      'FLOW_source', 'FLOW', ...
                                      'di_method', 'interpolate'); % interpolate constant
  % LISST-Tau options TODO: add TSG_source and CDOM source and good fDOM inteprolation like ACs
  elseif any(contains(lower(cfg.process.instruments2run{i}), {'lissttau','lisst-tau'}))
    cfg.process.calibrate.(cfg.process.instruments2run{i}) = struct('compute_dissolved', false, ...
                                      'interpolation_method', 'linear', ... % linear CDOM
                                      'CDOM_source', cfg.process.CDOM_source, ...
                                      'FLOW_source', 'FLOW', ...
                                      'di_method', 'normal');
  % SUVF options
  elseif any(contains(lower(cfg.process.instruments2run{i}), {'wscd','suvf'}))
    cfg.process.calibrate.(cfg.process.instruments2run{i}) = struct('compute_dissolved', false);
  end
end

%%% Write %%%
cfg.process.write = struct();
cfg.process.write.mode = 'One day one file'; % 'One file' or 'One day one file'
cfg.process.write.skip = {};

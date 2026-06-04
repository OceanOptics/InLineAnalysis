% TaraPacific Configuration file
% author: Guillaume Bourdin
% created: Aug 16, 2019

cfg = struct('meta', struct(), 'instruments', struct(), 'process', struct());

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                   METADATA (ENTER YOUR METADATA HERE)                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg.meta.investigators = 'Guillaume_Bourdin,Emmanuel_Boss';
cfg.meta.affiliations = 'University_of_Maine';
cfg.meta.emails = 'guillaume.bourdin@maine.edu';
cfg.meta.experiment = 'Tara';
cfg.meta.cruise = 'TaraPacific';
cfg.meta.station = 'NA';
cfg.meta.documents = 'NA';
cfg.meta.calibration_files = 'NA';
cfg.meta.data_type = 'flow_thru';
cfg.meta.data_status = 'preliminary';
cfg.meta.measurement_depth = 1.5;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                INSTRUMENTS (SETUP YOUR INSTRUMENTS HERE)                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
  PATH_ROOT = 'D:\Data\TaraPacific\';
elseif ismac
  PATH_ROOT = '/Volumes/Samsung_T5/Data/TaraPacific/';
end

%%% TSG + GPS %%%
model = 'SBE3845';
SN = '0091';
cfg.instruments.TSG = struct();
cfg.instruments.TSG.model = model;
cfg.instruments.TSG.TSG_source = true;
cfg.instruments.TSG.boat = 'Tara';
cfg.instruments.TSG.logger = 'Matlab'; % TeraTerm Matlab Inlinino
cfg.instruments.TSG.sn = SN;
cfg.instruments.TSG.path = struct('raw',  fullfile(PATH_ROOT, 'raw', 'TSG'),...
                                  'wk',   fullfile(PATH_ROOT, 'wk', 'TSG'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', 'TSG'),...
                                  'prod', fullfile(PATH_ROOT, 'prod'));
cfg.instruments.TSG.view = struct('varname', 't2');
cfg.instruments.TSG.temperature_variable = 't2';
cfg.instruments.TSG.salinity_variable = 's';

%%% FLOW (FlowControl) %%%
SN = 'FlowControl502';
model = 'FTH';
cfg.instruments.FLOW = struct();
cfg.instruments.FLOW.model = model;
cfg.instruments.FLOW.logger = 'FlowControl';
cfg.instruments.FLOW.sn = SN;
cfg.instruments.FLOW.LoadPrevious = true;
cfg.instruments.FLOW.analog1 = '';
cfg.instruments.FLOW.analog2 = '';
cfg.instruments.FLOW.path = struct('raw',  [PATH_ROOT 'raw' filesep 'FlowControl' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'FlowControl' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'FlowControl' filesep]);
cfg.instruments.FLOW.view = struct('varname', 'swt','swt_variable','swt','spd_variable','spd'); % spd1 spd2

%%% AC9 %%%
model = 'AC9';
SN = '245';
cfg.instruments.AC9 = struct();
cfg.instruments.AC9.model = 'AC9';
cfg.instruments.AC9.sn = SN;
cfg.instruments.AC9.logger = 'WetView';
cfg.instruments.AC9.device_file = [PATH_ROOT filesep 'DeviceFiles' filesep 'ac9245_20160816_20161002.dev'];
cfg.instruments.AC9.path = struct('raw',  [PATH_ROOT 'raw' filesep 'AC9' filesep],...
                                  'di',  [PATH_ROOT 'raw' filesep 'AC9' filesep 'DI' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'AC9' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'AC9' filesep]);
cfg.instruments.AC9.view = struct('varname', 'a', 'varcol', 5);

%%% ACS 091 %%% (Aug 20 to ...)
SN = '091';
cfg.instruments.(['ACS' SN]) = struct();
cfg.instruments.(['ACS' SN]).model = 'ACS';
cfg.instruments.(['ACS' SN]).sn = SN;
cfg.instruments.(['ACS' SN]).logger = 'Compass_2.1rc_scheduled_bin';
cfg.instruments.(['ACS' SN]).device_file = [PATH_ROOT filesep 'DeviceFiles' filesep 'acs091_20180530_20180818.dev']; % acs091_20170410_20170904 acs091_20180530_20180818 
cfg.instruments.(['ACS' SN]).path = struct('raw',  [PATH_ROOT 'raw' filesep 'ACS' filesep],...
                                  'di',  [PATH_ROOT 'raw' filesep 'ACS' filesep 'DI' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'ACS' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'ACS' filesep]);
cfg.instruments.(['ACS' SN]).view = struct('varname', 'a', 'varcol', 40);

%%% ACS 111 %%% (Aug 11 to Aug 20)
SN = '111';
cfg.instruments.(['ACS' SN]) = struct();
cfg.instruments.(['ACS' SN]).model = 'ACS';
cfg.instruments.(['ACS' SN]).sn = '111';
cfg.instruments.(['ACS' SN]).ila_prefix = 'ACS';
cfg.instruments.(['ACS' SN]).logger = 'Compass_2.1rc_scheduled_bin';
cfg.instruments.(['ACS' SN]).device_file = [PATH_ROOT filesep 'DeviceFiles' filesep 'acs111_20171212_20180530.dev'];
cfg.instruments.(['ACS' SN]).path = struct('raw',  [PATH_ROOT 'raw' filesep 'ACS' filesep],...
                                  'di',  [PATH_ROOT 'raw' filesep 'ACS' filesep 'DI' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'ACS' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'ACS' filesep]);
cfg.instruments.(['ACS' SN]).view = struct('varname', 'a', 'varcol', 40);

%%% ACS 007 %%% (Aug 11 to Aug 20)
SN = '007';
cfg.instruments.(['ACS' SN]) = struct();
cfg.instruments.(['ACS' SN]).model = 'ACS'; 
cfg.instruments.(['ACS' SN]).sn = SN;
cfg.instruments.(['ACS' SN]).ila_prefix = 'ACS';
cfg.instruments.(['ACS' SN]).logger = 'WetView'; % 'WetView' 'Compass_2.1rc' 'Compass_2.1rc_scheduled' 'Compass_2.1rc_scheduled_bin'
cfg.instruments.(['ACS' SN]).device_file = [PATH_ROOT filesep 'DeviceFiles' filesep 'acs007_20161101_20170220.dev']; % acs007_20160528_20160704 acs007_20161101_20170220
cfg.instruments.(['ACS' SN]).path = struct('raw',  [PATH_ROOT 'raw' filesep 'ACS' filesep],...
                                  'di',  [PATH_ROOT 'raw' filesep 'ACS' filesep 'DI' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'ACS' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'ACS' filesep]);
cfg.instruments.(['ACS' SN]).view = struct('varname', 'a', 'varcol', 40);

%%% ACS 057 %%% (Aug 11 to Aug 20)
SN = '057';
cfg.instruments.(['ACS' SN]) = struct();
cfg.instruments.(['ACS' SN]).model = 'ACS';
cfg.instruments.(['ACS' SN]).sn = SN;
cfg.instruments.(['ACS' SN]).ila_prefix = 'ACS';
cfg.instruments.(['ACS' SN]).logger = 'Compass_2.1rc_scheduled_bin';
cfg.instruments.(['ACS' SN]).device_file = [PATH_ROOT filesep 'DeviceFiles' filesep 'acs057_20160704_20160720.dev'];
cfg.instruments.(['ACS' SN]).path = struct('raw',  [PATH_ROOT 'raw' filesep 'ACS' filesep],...
                                  'di',  [PATH_ROOT 'raw' filesep 'ACS' filesep 'DI' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'ACS' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'ACS' filesep]);
cfg.instruments.(['ACS' SN]).view = struct('varname', 'a', 'varcol', 40);

%%% ACS 279 %%% (Aug 11 to Aug 20)
SN = '279';
cfg.instruments.(['ACS' SN]) = struct();
cfg.instruments.(['ACS' SN]).model = 'ACS';
cfg.instruments.(['ACS' SN]).sn = SN;
cfg.instruments.(['ACS' SN]).ila_prefix = 'ACS';
cfg.instruments.(['ACS' SN]).logger = 'Compass_2.1rc_scheduled_bin';
cfg.instruments.(['ACS' SN]).device_file = [PATH_ROOT filesep 'DeviceFiles' filesep 'acs279_20170902_20171213.dev']; % acs279_20170902_20171213 acs279_20180821_20180920
cfg.instruments.(['ACS' SN]).path = struct('raw',  [PATH_ROOT 'raw' filesep 'ACS' filesep],...
                                  'di',  [PATH_ROOT 'raw' filesep 'ACS' filesep 'DI' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'ACS' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'ACS' filesep]);
cfg.instruments.(['ACS' SN]).view = struct('varname', 'a', 'varcol', 40);

%%% BB3 %%%
SN = '1502';
cfg.instruments.BB3 = struct();
cfg.instruments.BB3.model = 'BB';
cfg.instruments.BB3.sn = SN;
cfg.instruments.BB3.ila_prefix = 'BB3';
cfg.instruments.BB3.logger = 'InlininoBB3';
cfg.instruments.BB3.lambda = [470,532,650];
cfg.instruments.BB3.theta = 120;
cfg.instruments.BB3.slope = [1.066E-05,7.076E-06,3.569E-06];
% cfg.instruments.BB3.slope = [8.407E-06,4.624E-06,4.090E-06];
cfg.instruments.BB3.dark = [50,44,45];
cfg.instruments.BB3.path = struct('raw',  [PATH_ROOT 'raw' filesep 'BB3' filesep],...
                                  'di',  [PATH_ROOT 'raw' filesep 'BB3' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'BB3' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'BB3' filesep]);
cfg.instruments.BB3.view = struct('varname', 'beta', 'varcol', 2);

%%% WSCD %%% (Mai 2016 to Oct 2018)
model = 'WSCD';
SN = '1082P';
cfg.instruments.WSCD1082P = struct();
cfg.instruments.WSCD1082P.model = 'CD';
cfg.instruments.WSCD1082P.CDOM_source = true;
cfg.instruments.WSCD1082P.sn = SN;
cfg.instruments.WSCD1082P.ila_prefix = 'WSCD';
cfg.instruments.WSCD1082P.logger = 'InlininoWSCD';
cfg.instruments.WSCD1082P.slope = 62;
cfg.instruments.WSCD1082P.dark = 0.059;
cfg.instruments.WSCD1082P.path = struct('raw',  [PATH_ROOT 'raw' filesep 'WSCD' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'WSCD' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'WSCD' filesep]);
cfg.instruments.WSCD1082P.view = struct('varname', 'fdom');

%%% PAR %%% (Mai 2016 to Oct 2018)
SN = '50168';
cfg.instruments.PAR = struct();
cfg.instruments.PAR.model = 'PAR';
cfg.instruments.PAR.sn = SN;
cfg.instruments.PAR.logger = 'Inlinino_base';
cfg.instruments.PAR.scale = 6.451E-04; % Volts/(uE/m²sec)
cfg.instruments.PAR.dark = 9.7E-03;
cfg.instruments.PAR.path = struct('raw',  [PATH_ROOT 'raw' filesep 'PAR' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'PAR' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'PAR' filesep]);
cfg.instruments.PAR.view = struct('varname', 'par');

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
  {'flow','tsg','sbe45','sbe3845','nmea','par','qcr','qsp'}));
cfg.process.di.qc = struct('mode', 'ui',... % ui or load
                           'qc_once_for_all', false,... % true = QC all variables | false = QC variables separately);
                           'remove_old', false,... % remove old selection of the same period
                           'remove_when_flow_below', false); % true = remove data when flow <= 0.5 | false = no data removal data depending on flow | number = remove data when flow <= number)
cfg.process.di.bin = struct('bin_size', 30);

%%% Synchronization %%%
cfg.process.sync = struct();
cfg.process.sync.delay = struct();
cfg.process.sync.skip = cfg.process.instruments2run(contains(lower(cfg.process.instruments2run), ...
  {'flow','tsg','sbe45','sbe3845','nmea','par', 'alfa','qcr','qsp'}));
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
  {'flow','tsg','sbe45','sbe3845','nmea','par', 'alfa','qcr','qsp'}));
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
cfg.process.bin.prctile_detection = [2.5, 97.5];
% Ff StepQC with ACS: prctile_average = [2.5, 97.5]; otherwise prctile_average = [5, 75];
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
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = minutes(1); % 1 min for LISST200X
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
% cfg.process.bin.bin_size.LISST1183 = minutes(10);
% cfg.process.bin.bin_size.SUVF6244 = minutes(1);
% cfg.process.bin.bin_size.ALFA = minutes(10);
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
  {'flow','nmea','par','qcr','qsp'}));
cfg.process.qc.specific = struct();
cfg.process.qc.specific.active = true;
cfg.process.qc.specific.run = {cfg.process.qcref.view};

%%% Calibrate %%%
cfg.process.calibrate = struct();
cfg.process.calibrate.skip = cfg.process.instruments2run(contains(lower(cfg.process.instruments2run), ...
  {'flow','tsg','sbe45','sbe3845','nmea','alfa'}));
% cfg.process.min_nb_pts_per_cluster = 200;
% look for TSG and SUVF and automatically turn off CDOM interpolation if not available
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
% if no TSG and CDOM source indicated in cfg: find TSG and CDOM instruments automatically
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

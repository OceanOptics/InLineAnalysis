% Tara Breizh Bloom Configuration file
% author: Nils Haentjens
% created: Mai 27, 2019

cfg = struct('meta', struct(), 'instruments', struct(), 'process', struct());

%%%%%%%%%%%%%%
%% METADATA %%
%%%%%%%%%%%%%%
cfg.meta.investigators = 'Emmanuel_Boss,Nils_Haentjens,Louis_Terrats';
cfg.meta.affiliations = 'University_of_Maine';
cfg.meta.emails = 'emmanuel.boss@maine.edu';
cfg.meta.experiment = 'Breizh_Bloom';
cfg.meta.cruise = 'Tara_Breizh_Bloom';
cfg.meta.station = 'NA';
cfg.meta.documents = 'NA';
cfg.meta.calibration_files = 'NA';
cfg.meta.data_type = 'flow_thru';
cfg.meta.data_status = 'preliminary';
cfg.meta.measurement_depth = 1.5;


%%%%%%%%%%%%%%%%%
%% INSTRUMENTS %%
%%%%%%%%%%%%%%%%%

PATH_ROOT = '/Users/nils/Data/TaraBreizhBloom/';

%%% TSG + GPS %%%
cfg.instruments.TSG = struct();
cfg.instruments.TSG.model = 'TSG';
cfg.instruments.TSG.boat = 'Tara';
cfg.instruments.TSG.path = struct('raw',  [PATH_ROOT 'raw/TSG/'],...
                                  'wk',   [PATH_ROOT 'wk/TSG/'],...
                                  'prod', [PATH_ROOT 'prod/'],...
                                  'ui', [PATH_ROOT 'ui/TSG/']);
cfg.instruments.TSG.view = struct('varname', 't');

%%% FTH (FlowControl) %%%
cfg.instruments.FTH = struct();
cfg.instruments.FTH.model = 'FTH';
cfg.instruments.FTH.logger = 'FlowControl';
cfg.instruments.FTH.path = struct('raw',  [PATH_ROOT 'raw/FlowControl/'],...
                                  'wk',   [PATH_ROOT 'wk/FlowControl/'],...
                                  'prod', [PATH_ROOT 'prod/'],...
                                  'ui', [PATH_ROOT 'ui/FlowControl/']);
cfg.instruments.FTH.view = struct('varname', 'swt');

%%% ACS 301 %%%
cfg.instruments.ACS301 = struct();
cfg.instruments.ACS301.model = 'ACS';
cfg.instruments.ACS301.sn = '301';
cfg.instruments.ACS301.logger = 'Compass_2.1rc_scheduled_bin';
cfg.instruments.ACS301.device_file = [PATH_ROOT 'DeviceFiles/acs301.dev'];
cfg.instruments.ACS301.path = struct('raw',  [PATH_ROOT 'raw/ACS301/'],...
                                    'di',  [PATH_ROOT 'raw/ACS301/DI/'],...
                                    'wk',   [PATH_ROOT 'wk/ACS301/'],...
                                    'prod', [PATH_ROOT 'prod/'],...
                                    'ui', [PATH_ROOT 'ui/ACS301/']);
cfg.instruments.ACS301.view = struct('varname', 'a', 'varcol', 40);

%%% BB3 %%%
cfg.instruments.BB31502 = struct();
cfg.instruments.BB31502.model = 'BB';
cfg.instruments.BB31502.sn = '1502';
cfg.instruments.BB31502.ila_prefix = 'BB3';
cfg.instruments.BB31502.logger = 'InlininoBB3SN';
% Pre-Tara Pacific Calibration from 2016/10/20 ECO BB3-1502
cfg.instruments.BB31502.lambda = [470,532,650];
cfg.instruments.BB31502.theta = 124;
cfg.instruments.BB31502.slope = [1.066e-05, 7.076e-06, 3.569e-06];
cfg.instruments.BB31502.dark = [50, 44, 46];
cfg.instruments.BB31502.path = struct('raw',  [PATH_ROOT 'raw/BB31502/'],...
                                  'di',  [PATH_ROOT 'raw/BB31502/'],...
                                  'wk',   [PATH_ROOT 'wk/BB31502/'],...
                                  'prod', [PATH_ROOT 'prod/'],...
                                  'ui', [PATH_ROOT 'ui/BB31502/']);
cfg.instruments.BB31502.view = struct('varname', 'beta', 'varcol', 2);

%%% WSCD %%%
cfg.instruments.WSCD1082 = struct();
cfg.instruments.WSCD1082.model = 'CD';
cfg.instruments.WSCD1082.sn = '1082';
cfg.instruments.WSCD1082.ila_prefix = 'WSCD1082';
cfg.instruments.WSCD1082.logger = 'InlininoWSCDSN';
cfg.instruments.WSCD1082.slope = NaN;
cfg.instruments.WSCD1082.dark = NaN;
cfg.instruments.WSCD1082.path = struct('raw',  [PATH_ROOT 'raw/WSCD1082/'],...
                                  'wk',   [PATH_ROOT 'wk/WSCD1082/'],...
                                  'prod', [PATH_ROOT 'prod/'],...
                                  'ui', [PATH_ROOT 'ui/WSCD1082/']);
cfg.instruments.WSCD1082.view = struct('varname', 'fdom');

%%%%%%%%%%%%%
%% PROCESS %%
%%%%%%%%%%%%%

%%% General parameters %%%
cfg.process.days2run = datenum(2019,05,27):datenum(2019,06,02);
cfg.process.instruments2run = {'FTH', 'TSG', 'ACS301', 'BB31502', 'WSCD1082'};
cfg.process.write = true;
cfg.process.force_import = false;
cfg.process.parallel = Inf; % 0: disable parallel or Inf: as many thread available
cfg.process.di = struct();
cfg.process.di.skip = {'FTH', 'TSG', 'ACS301', 'BB31502', 'WSCD1082'};
cfg.process.di.qc = struct('mode', 'load');
cfg.process.di.bin = struct('bin_size', 30);

%%% Synchronization %%%
cfg.process.sync = struct();
cfg.process.sync.delay = struct();
cfg.process.sync.delay.FTH = 0;
cfg.process.sync.delay.ACS301 = 66;
cfg.process.sync.delay.BB31502 = 0;
cfg.process.sync.delay.WSCD1082 = 5;
cfg.process.sync.skip = {'TSG'};

%%% QC Reference (Flow Control/FTH) %%%
cfg.process.qcref = struct();
cfg.process.qcref.reference = 'FTH';
cfg.process.qcref.view = 'ACS301';
cfg.process.qcref.mode = 'load'; % load or ui

%%% Split total/filtered %%%
cfg.process.split = struct();
cfg.process.split.reference = 'FTH';
cfg.process.split.buffer = struct();
cfg.process.split.buffer.ACS301 = [180, 30];
cfg.process.split.buffer.BB31502 = [420, 220];
cfg.process.split.buffer.WSCD1082 = [310, 20];
cfg.process.split.skip = {'FTH', 'TSG'};

%%% Binning %%%
cfg.process.bin = struct('bin_size', struct());
cfg.process.bin.prctile_detection = [2.5, 97.5];
cfg.process.bin.prctile_average = [5, 75];
% Bin mode does not affect the outcome of the data but just the way the data is presented to the computer
% cfg.process.bin.mode = 'OneShot'; % Faster for small dataset fiting in the memory of the computer
cfg.process.bin.mode = 'ByDay'; % Slightly slower but can handle a lot more data at once as it will be binned one day at a time
cfg.process.bin.bin_size.ACS301 = 1;
cfg.process.bin.bin_size.BB31502 = 1;
cfg.process.bin.bin_size.WSCD1082 = 1;
cfg.process.bin.bin_size.TSG = 1;
cfg.process.bin.skip = {'FTH'};

%%% Automatically flagging %%%
cfg.process.flag = struct();
cfg.process.flag.skip = {'FTH', 'TSG', 'ACS301', 'BB31502', 'WSCD1082'};
% Default: parameters set to all instruments if not specific parameters set
cfg.process.flag.default = struct();
cfg.process.flag.default.maximum_fudge_factor = 4;
cfg.process.flag.default.variance_fudge_factor = 3;
cfg.process.flag.default.avg_sensitivity = 1;
cfg.process.flag.default.unc1_sensitivity = 1;
cfg.process.flag.default.unc2_sensitivity = 2;
cfg.process.flag.default.smooth_threshold = 60;
cfg.process.flag.default.min_flag_n = 1;
cfg.process.flag.default.filt = struct('smooth_threshold', 2);
  
%%% Manually QC %%%
cfg.process.qc = struct();
cfg.process.qc.mode = 'ui';
cfg.process.qc.global = struct();
cfg.process.qc.global.active = false;
cfg.process.qc.global.view = 'ACS301';
cfg.process.qc.global.apply = {'ACS301', 'BB31502', 'WSCD'};
cfg.process.qc.specific = struct();
cfg.process.qc.specific.active = true;
cfg.process.qc.specific.run = {'TSG', 'ACS301', 'BB31502', 'WSCD1082'};
  
%%% Calibrate %%%
cfg.process.calibrate = struct();
cfg.process.calibrate.ACS301 = struct('compute_dissolved', false,...
                                  'interpolation_method', 'linear', ...
                                  'CDOM_source', 'WSCD1082',... % Pay attention to serial number (could merge WSCD products first)
                                  'FTH_source', 'FTH');
cfg.process.calibrate.BB31502 = struct('compute_dissolved', false,...
                                   'TSG_source', 'TSG',...
                                   'di_method', 'constant');
cfg.process.calibrate.skip = {'FTH', 'TSG', 'WSCD1082'};

%%% Write %%%
cfg.process.write = struct();
cfg.process.write.mode = 'One day one file'; % 'One file' or 'One day one file'
cfg.process.write.skip = {'FTH'};

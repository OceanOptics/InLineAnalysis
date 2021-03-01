% Tara Breizh Bloom Configuration file
% author: Nils Haentjens & Guillaume Bourdin
% created: Mai 27, 2019

cfg = struct('meta', struct(), 'instruments', struct(), 'process', struct());

%%%%%%%%%%%%%%
%% METADATA %%
%%%%%%%%%%%%%%
cfg.meta.investigators = 'Emmanuel_Boss,Nils_Haentjens,Louis_Terrats,Guillaume_Bourdin';
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

PATH_ROOT = 'C:\Users\Gui\Documents\MATLAB\InLineAnalysis\';
PATH_DATA = 'E:\Data\TaraBreizhBloom\';

%%% TSG + GPS %%%
cfg.instruments.TSG = struct();
cfg.instruments.TSG.model = 'TSG';
cfg.instruments.TSG.boat = 'Tara';
cfg.instruments.TSG.path = struct('raw',  [PATH_DATA 'raw' filesep 'TSG' filesep],...
                                  'wk',   [PATH_DATA 'wk' filesep 'TSG' filesep],...
                                  'prod', [PATH_DATA 'prod' filesep],...
                                  'ui', [PATH_DATA 'ui' filesep 'TSG' filesep]);
cfg.instruments.TSG.view = struct('varname', 't');

%%% FTH (FlowControl) %%%
cfg.instruments.FTH = struct();
cfg.instruments.FTH.model = 'FTH';
cfg.instruments.FTH.logger = 'FlowControl';
cfg.instruments.FTH.LoadPrevious = true;
cfg.instruments.FTH.path = struct('raw',  [PATH_DATA 'raw' filesep 'FlowControl' filesep],...
                                  'wk',   [PATH_DATA 'wk' filesep 'FlowControl' filesep],...
                                  'prod', [PATH_DATA 'prod' filesep],...
                                  'ui', [PATH_DATA 'ui' filesep 'FlowControl' filesep]);
cfg.instruments.FTH.view = struct('varname', 'swt');

%%% ACS301 %%%
SN = '301';
cfg.instruments.(['ACS' SN]) = struct();
cfg.instruments.(['ACS' SN]).model = 'ACS';
cfg.instruments.(['ACS' SN]).sn = SN;
cfg.instruments.(['ACS' SN]).logger = 'InlininoACScsv';
cfg.instruments.(['ACS' SN]).device_file = [PATH_DATA filesep 'DeviceFiles' filesep 'acs301_20180806.dev'];
cfg.instruments.(['ACS' SN]).path = struct('raw',  [PATH_DATA 'raw' filesep ['ACS' SN] filesep],...
                                    'di',  [PATH_DATA 'raw' filesep ['ACS' SN] filesep 'DI' filesep],...
                                    'wk',   [PATH_DATA 'wk' filesep ['ACS' SN] filesep],...
                                    'prod', [PATH_DATA 'prod' filesep],...
                                    'ui', [PATH_DATA 'ui' filesep ['ACS' SN] filesep]);
cfg.instruments.ACS57.view = struct('varname', 'a', 'varcol', 40);

%%% BB3 %%%
SN = '1502';
cfg.instruments.(['BB3' SN]) = struct();
cfg.instruments.(['BB3' SN]).model = 'BB';
cfg.instruments.(['BB3' SN]).sn = SN;
cfg.instruments.(['BB3' SN]).ila_prefix = 'BB3';
cfg.instruments.(['BB3' SN]).logger = 'InlininoBB3SN';
% Pre-Tara Pacific Calibration from 2016/10/20 ECO BB3-1502
cfg.instruments.(['BB3' SN]).lambda = [470,532,650];
cfg.instruments.(['BB3' SN]).theta = 124;
cfg.instruments.(['BB3' SN]).slope = [1.066e-05, 7.076e-06, 3.569e-06];
cfg.instruments.(['BB3' SN]).dark = [50, 44, 46];
cfg.instruments.(['BB3' SN]).path = struct('raw',  [PATH_DATA 'raw' filesep ['BB3' SN] filesep],...
                                  'di',  [PATH_DATA 'raw' filesep ['BB3' SN] filesep],...
                                  'wk',   [PATH_DATA 'wk' filesep ['BB3' SN] filesep],...
                                  'prod', [PATH_DATA 'prod' filesep],...
                                  'ui', [PATH_DATA 'ui' filesep ['BB3' SN] filesep]);
cfg.instruments.(['BB3' SN]).view = struct('varname', 'beta', 'varcol', 2);

%%% WSCD %%%
SN = '1082';
cfg.instruments.(['WSCD' SN]) = struct();
cfg.instruments.(['WSCD' SN]).sn = SN;
cfg.instruments.(['WSCD' SN]).model = 'CD';
cfg.instruments.(['WSCD' SN]).ila_prefix = ['WSCD' SN];
cfg.instruments.(['WSCD' SN]).logger = 'InlininoWSCDSN';
cfg.instruments.(['WSCD' SN]).slope = 62;
cfg.instruments.(['WSCD' SN]).dark = 0.059;
cfg.instruments.(['WSCD' SN]).path = struct('raw',  [PATH_DATA 'raw' filesep ['WSCD' SN] filesep],...
                                  'wk',   [PATH_DATA 'wk' filesep ['WSCD' SN] filesep],...
                                  'prod', [PATH_DATA 'prod' filesep],...
                                  'ui', [PATH_DATA 'ui' filesep ['WSCD' SN] filesep]);
cfg.instruments.(['WSCD' SN]).view = struct('varname', 'fdom');

%%%%%%%%%%%%%
%% PROCESS %%
%%%%%%%%%%%%%

%%% General parameters %%%
cfg.process.days2run = datenum(2019,05,18):datenum(2019,06,03);
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
cfg.process.qcref.mode = 'ui'; % load or ui
cfg.process.qcref.MinFiltPeriod = 50; % filter even period in minute

%%% Split total/filtered %%%
cfg.process.split = struct();
cfg.process.split.reference = 'FTH';
cfg.process.split.buffer = struct();
cfg.process.split.buffer.ACS301 = [180, 60];
cfg.process.split.buffer.BB31502 = [400, 220];
cfg.process.split.buffer.WSCD1082 = [310, 20];
cfg.process.split.skip = {'FTH','TSG'};

%%% Binning %%%
cfg.process.bin = struct('bin_size', struct());
cfg.process.bin.prctile_detection = [2.5, 97.5];
cfg.process.bin.prctile_average = [5, 75];
% Bin mode does not affect the outcome of the data but just the way the data is presented to the computer
% cfg.process.bin.mode = 'OneShot'; % Faster for small dataset fiting in the memory of the computer
cfg.process.bin.mode = 'ByDay'; % Slightly slower but can handle a lot more data at once as it will be binned one day at a time
cfg.process.bin.bin_size.FTH = 1;
cfg.process.bin.bin_size.ACS301 = 1;
cfg.process.bin.bin_size.BB31502 = 1;
cfg.process.bin.bin_size.WSCD1082 = 1;
cfg.process.bin.bin_size.TSG = 1;
cfg.process.bin.skip = {};

%%% Automatically flagging %%%
cfg.process.flag = struct();
cfg.process.flag.skip = {'FTH', 'TSG', 'ACS57', 'BB31502', 'WSCD859'};
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
cfg.process.qc.StepQCLim.filtered.a = 3;
cfg.process.qc.StepQCLim.filtered.c = 3;
cfg.process.qc.StepQCLim.total.a = 3;
cfg.process.qc.StepQCLim.total.c = 3;
cfg.process.qc.StepQCLim.filtered.bb = 3;
cfg.process.qc.StepQCLim.total.bb = 3;
cfg.process.qc.Saturation_Threshold_bb = 4000; % (counts)
cfg.process.qc.mode = 'ui';
cfg.process.qc.global = struct();
cfg.process.qc.global.active = false;
cfg.process.qc.global.view = 'ACS301';
cfg.process.qc.global.apply = {'ACS301', 'BB31502', 'WSCD1082'};
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
cfg.process.calibrate.skip = {'FTH', 'TSG'};

%%% Write %%%
cfg.process.write = struct();
cfg.process.write.mode = 'One day one file'; % 'One file' or 'One day one file'
cfg.process.write.skip = {};

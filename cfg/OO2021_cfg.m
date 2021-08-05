% EXPORTS02 Configuration file
% author: Guillaume
% created: April 2021

cfg = struct('meta', struct(), 'instruments', struct(), 'process', struct());

%%%%%%%%%%%%%%
%% METADATA %%
%%%%%%%%%%%%%%
cfg.meta.investigators = 'Guillaume Bourdin, Emmanuel Boss, Collin Roesler, Margaret Estapa, Kenneth Voss, Curtis Mobley, OO2021 Students';
cfg.meta.affiliations = 'University_of_Maine';
cfg.meta.emails = 'emmanuel.boss@maine.edu';
cfg.meta.experiment = 'OO2021';
cfg.meta.cruise = 'OO2021';
cfg.meta.station = 'NA';
cfg.meta.documents = 'NA';
cfg.meta.calibration_files = 'NA';
cfg.meta.data_type = 'flow_thru';
cfg.meta.data_status = 'preliminary';
cfg.meta.measurement_depth = 2;


%%%%%%%%%%%%%%%%%
%% INSTRUMENTS %%
%%%%%%%%%%%%%%%%%

PATH_ROOT = '/Volumes/Samsung_T5/Data/OO2021/';

%%% TSG + GPS %%%
cfg.instruments.TSG = struct();
cfg.instruments.TSG.model = 'TSG';
cfg.instruments.TSG.boat = 'IRA-C';
cfg.instruments.TSG.logger = 'Matlab'; % SBE45TSG matlab_Emmanuel
cfg.instruments.TSG.path = struct('raw',  [PATH_ROOT 'raw' filesep 'TSG' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'TSG' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'TSG' filesep]);
cfg.instruments.TSG.view = struct('varname', 't');

%%% FLOW (FlowControl) %%%
cfg.instruments.FLOW = struct();
cfg.instruments.FLOW.model = 'FTH';
cfg.instruments.FLOW.logger = 'FlowControl';
cfg.instruments.FLOW.LoadPrevious = true;
cfg.instruments.FLOW.path = struct('raw',  [PATH_ROOT 'raw' filesep 'FlowControl' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'FlowControl' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'FlowControl' filesep]);
cfg.instruments.FLOW.view = struct('varname', 'swt');

%%% ACS %%%
SN = '94';
cfg.instruments.(['ACS' SN]) = struct();
cfg.instruments.(['ACS' SN]).di = struct();
cfg.instruments.(['ACS' SN]).di.prefix = ['DIW_ACS' SN '_'];
cfg.instruments.(['ACS' SN]).di.postfix = '';
cfg.instruments.(['ACS' SN]).model = 'ACS';
cfg.instruments.(['ACS' SN]).sn = SN;
cfg.instruments.(['ACS' SN]).logger = 'InlininoACScsv';
cfg.instruments.(['ACS' SN]).device_file = [PATH_ROOT filesep 'DeviceFiles' filesep 'acs094_20190116.dev'];
cfg.instruments.(['ACS' SN]).path = struct('raw',  [PATH_ROOT 'raw' filesep ['ACS' SN] filesep],...
                                    'di',  [PATH_ROOT 'raw' filesep ['ACS' SN] filesep 'DI' filesep],...
                                    'wk',   [PATH_ROOT 'wk' filesep ['ACS' SN] filesep],...
                                    'prod', [PATH_ROOT 'prod' filesep],...
                                    'ui', [PATH_ROOT 'ui' filesep ['ACS' SN] filesep]);
cfg.instruments.(['ACS' SN]).view = struct('varname', 'a', 'varcol', 40);

%%% AC9 %%%
SN = '274';
cfg.instruments.(['AC9' SN]) = struct();
cfg.instruments.(['AC9' SN]).di = struct();
cfg.instruments.(['AC9' SN]).di.prefix = ['DIW_ACS' SN '_'];
cfg.instruments.(['AC9' SN]).di.postfix = '';
cfg.instruments.(['AC9' SN]).model = 'AC9';
cfg.instruments.(['AC9' SN]).sn = SN;
cfg.instruments.(['AC9' SN]).logger = 'WetView';
cfg.instruments.(['AC9' SN]).device_file = [PATH_ROOT filesep 'DeviceFiles' filesep 'ac90274_20210402.dev'];
cfg.instruments.(['AC9' SN]).path = struct('raw',  [PATH_ROOT 'raw' filesep ['AC9' SN] filesep],...
                                    'di',  [PATH_ROOT 'raw' filesep ['AC9' SN] filesep 'DI' filesep],...
                                    'wk',   [PATH_ROOT 'wk' filesep ['AC9' SN] filesep],...
                                    'prod', [PATH_ROOT 'prod' filesep],...
                                    'ui', [PATH_ROOT 'ui' filesep ['AC9' SN] filesep]);
cfg.instruments.(['AC9' SN]).view = struct('varname', 'a', 'varcol', 5);



%%% BB3 %%%
SN = '349';
cfg.instruments.BB3 = struct();
cfg.instruments.BB3.model = 'BB';
cfg.instruments.BB3.sn = SN;
cfg.instruments.BB3.ila_prefix = ['BB3' SN];
cfg.instruments.BB3.logger = 'InlininoBB3SN';
cfg.instruments.BB3.lambda = [470,532,660];
cfg.instruments.BB3.theta = 124;
cfg.instruments.BB3.slope = [8.407E-06,4.624E-06,4.090E-06];
cfg.instruments.BB3.dark = [57,48,45];
cfg.instruments.BB3.path = struct('raw',  [PATH_ROOT 'raw' filesep ['BB3' SN] filesep],...
                                  'di',  [PATH_ROOT 'raw' filesep ['BB3' SN] filesep 'DI' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep ['BB3' SN] filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep ['BB3' SN] filesep]);
cfg.instruments.BB3.view = struct('varname', 'beta', 'varcol', 2);

%%%%%%%%%%%%%
%% PROCESS %%
%%%%%%%%%%%%%

%%% General parameters %%%
cfg.process.days2run = datenum(2021,5,2):datenum(2021,5,5);
cfg.process.instruments2run = {'FLOW', 'TSG', 'ACS94', 'AC9274', 'BB3'};
cfg.process.write = true;
cfg.process.force_import = false;
cfg.process.parallel = Inf; % 0: disable parallel or Inf: as many thread available
cfg.process.di = struct();
cfg.process.di.skip = {'FLOW', 'TSG'};
cfg.process.di.qc = struct('mode', 'ui',... % ui or load
                           'qc_once_for_all', false,... % true = QC all variables | false = QC variables separately)
                           'remove_old', false); % remove old selection of the same period
cfg.process.di.bin = struct('bin_size', 30);

%%% Synchronization %%%
cfg.process.sync = struct();
cfg.process.sync.delay = struct();
cfg.process.sync.delay.FLOW = 0;
cfg.process.sync.delay.ACS94 = 0;
cfg.process.sync.delay.AC9274 = 0;
cfg.process.sync.delay.BB3 = 0;
cfg.process.sync.skip = {'TSG'};

%%% QC Reference (Flow Control/FLOW) %%%
cfg.process.qcref = struct();
cfg.process.qcref.reference = 'FLOW';
cfg.process.qcref.view = 'BB3';
cfg.process.qcref.mode = 'ui'; % load or ui
cfg.process.qcref.remove_old = false; % remove old selection of the same period
cfg.process.qcref.MinFiltPeriod = 50; % filter event period in minute
cfg.process.qcref.szFilt = 10; % filter even length in minute

%%% Split total/filtered %%%
cfg.process.split = struct();
cfg.process.split.reference = 'FLOW';
cfg.process.split.buffer = struct();
cfg.process.split.buffer.ACS94 = [180, 30];
cfg.process.split.buffer.AC9274 = [180, 30];
cfg.process.split.buffer.BB3 = [220, 150];
cfg.process.split.skip = {'FLOW', 'TSG'};

%%% Binning %%%
cfg.process.bin = struct('bin_size', struct());
cfg.process.bin.prctile_detection = [2.5, 97.5];
% if RawAutoQC with ACS or BB: prctile_average = [2.5, 97.5]; otherwise prctile_average = [5, 75];
cfg.process.bin.prctile_average = [2.5, 97.5];
% Bin mode does not affect the outcome of the data but just the way the data is presented to the computer
% cfg.process.bin.mode = 'OneShot'; % Faster for small dataset fiting in the memory of the computer
cfg.process.bin.mode = 'ByDay'; % Slightly slower but can handle a lot more data at once as it will be binned one day at a time
cfg.process.bin.bin_size.FLOW = 1;
cfg.process.bin.bin_size.ACS94 = 1;
cfg.process.bin.bin_size.AC9274 = 1;
cfg.process.bin.bin_size.BB3 = 1;
cfg.process.bin.bin_size.TSG = 1; % TSG does not need to be binned in most cases
cfg.process.bin.skip = {};

%%% Automatically flagging %%%
cfg.process.flag = struct();
cfg.process.flag.skip = cfg.process.instruments2run;
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
cfg.process.qc.global.view = 'ACS94';
cfg.process.qc.global.apply = {'ACS94', 'BB3'};
cfg.process.qc.specific = struct();
cfg.process.qc.specific.active = true;
cfg.process.qc.specific.run = {'TSG', 'ACS94', 'BB3'};
  
%%% Calibrate %%%
cfg.process.calibrate = struct();
cfg.process.calibrate.ACS94 = struct('compute_dissolved', true, ...
                                  'interpolation_method', 'linear', ...
                                  'CDOM_source', [], ... 
                                  'FLOW_source', 'FLOW', ...
                                  'di_method', 'normal', ... % best_di normal
                                  'compute_ad_aphi', false); % VERY SLOW: compute ad and aphi from Zheng and Stramski 2013
cfg.process.calibrate.AC9274 = struct('compute_dissolved', true, ...
                                  'interpolation_method', 'linear', ...
                                  'CDOM_source', [], ... 
                                  'FLOW_source', 'FLOW', ...
                                  'di_method', 'normal', ... % best_di normal
                                  'compute_ad_aphi', false); % VERY SLOW: compute ad and aphi from Zheng and Stramski 2013
cfg.process.calibrate.BB3 = struct('compute_dissolved', true, ...
                                  'TSG_source', 'TSG', ...
                                  'FLOW_source', 'FLOW', ...
                                  'di_method', 'interpolate', ... % interpolate constant SW_scattering
                                  'filt_method', 'exponential_fit'); % 25percentil exponential_fit
cfg.process.calibrate.skip = {'FLOW', 'TSG'};

%%% Write %%%
cfg.process.write = struct();
cfg.process.write.mode = 'One day one file'; % 'One file' or 'One day one file'
cfg.process.write.skip = {};

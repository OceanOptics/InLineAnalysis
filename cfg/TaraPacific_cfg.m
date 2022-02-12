% TaraPacific Configuration file
% author: Guillaume Bourdin
% created: Aug 16, 2019

cfg = struct('meta', struct(), 'instruments', struct(), 'process', struct());

%%%%%%%%%%%%%%
%% METADATA %%
%%%%%%%%%%%%%%
cfg.meta.investigators = 'Emmanuel_Boss,Guillaume_Bourdin';
cfg.meta.affiliations = 'University_of_Maine';
cfg.meta.emails = 'emmanuel.boss@maine.edu';
cfg.meta.experiment = 'Tara';
cfg.meta.cruise = 'TaraPacific';
cfg.meta.station = 'NA';
cfg.meta.documents = 'NA';
cfg.meta.calibration_files = 'NA';
cfg.meta.data_type = 'flow_thru';
cfg.meta.data_status = 'preliminary';
cfg.meta.measurement_depth = 1.5;


%%%%%%%%%%%%%%%%%
%% INSTRUMENTS %%
%%%%%%%%%%%%%%%%%

if ispc
  PATH_ROOT = 'D:\Data\TaraPacific\';
elseif ismac
  PATH_ROOT = '/Volumes/Samsung_T5/Data/TaraPacific/';
end

%%% TSG + GPS %%%
cfg.instruments.TSG = struct();
cfg.instruments.TSG.model = 'TSG';
cfg.instruments.TSG.boat = 'Tara';
cfg.instruments.TSG.logger = 'Matlab'; % TeraTerm Matlab
cfg.instruments.TSG.path = struct('raw',  [PATH_ROOT 'raw' filesep 'TSG' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'TSG' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'TSG' filesep]);
cfg.instruments.TSG.view = struct('varname', 't2');
cfg.instruments.TSG.temperature_variable = 't2';


%%% FLOW (FlowControl) %%%
cfg.instruments.FLOW = struct();
cfg.instruments.FLOW.model = 'FTH';
cfg.instruments.FLOW.logger = 'FlowControl';
cfg.instruments.FLOW.LoadPrevious = true;
cfg.instruments.FLOW.path = struct('raw',  [PATH_ROOT 'raw' filesep 'FlowControl' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'FlowControl' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'FlowControl' filesep]);
cfg.instruments.FLOW.view = struct('varname', 'swt','spd_variable','spd1'); % spd1 spd2

%%% AC9 %%%
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
cfg.instruments.BB3.theta = 124;
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
SN = '1082P';
cfg.instruments.WSCD1082P = struct();
cfg.instruments.WSCD1082P.model = 'CD';
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
cfg.instruments.PAR.logger = 'Inlinino';
cfg.instruments.PAR.scale = 6.451E-04; % Volts/(uE/m²sec)
cfg.instruments.PAR.dark = 9.7E-03;
cfg.instruments.PAR.path = struct('raw',  [PATH_ROOT 'raw' filesep 'PAR' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'PAR' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'PAR' filesep]);
cfg.instruments.PAR.view = struct('varname', 'par');

%% %%%%%%% %%
%  PROCESS  %
%%%%%%%%%%%%%

%%% General parameters %%%
cfg.process.days2run = datenum(2018,5,30,0,0,0):datenum(2018,6,20,0,0,0);
cfg.process.instruments2run = {'FLOW'};
cfg.process.write = true;
cfg.process.force_import = false;
cfg.process.parallel = Inf; % 0: disable parallel or Inf: as many thread available
cfg.process.di = struct();
cfg.process.di.skip = {'FLOW', 'TSG'};
cfg.process.di.qc = struct('mode', 'ui',... % ui or load
                           'qc_once_for_all', false,... % true = QC all variables | false = QC variables separately);
                           'remove_old', false); % remove old selection of the same period
cfg.process.di.bin = struct('bin_size', 30);

%%% Synchronization %%%
cfg.process.sync = struct();
cfg.process.sync.delay = struct();
cfg.process.sync.delay.FLOW = 30;
cfg.process.sync.delay.AC9 = 60; % 
cfg.process.sync.delay.ACS007 = 60; % 
cfg.process.sync.delay.ACS057 = 60; % 
cfg.process.sync.delay.ACS091 = 60; % 280 95 65 82 60 57 52 -10 52 95
cfg.process.sync.delay.ACS111 = 60; % 280 95 65 82 60 57
cfg.process.sync.delay.ACS279 = 60; %70 128 75 62 58 10 60
cfg.process.sync.delay.BB3 = 0; % 55 45 22 20
% cfg.process.sync.delay.LISST = 1;
% cfg.process.sync.delay.WSCD859 = 5;
cfg.process.sync.delay.WSCD1082P = 40;
% cfg.process.sync.delay.ALFA = 15;
cfg.process.sync.skip = {'TSG','PAR'};

%%% QC Reference (Flow Control/FLOW) %%%
cfg.process.qcref = struct();
cfg.process.qcref.reference = 'FLOW';
cfg.process.qcref.view = 'ACS007';
cfg.process.qcref.mode = 'ui'; % load or ui
cfg.process.qcref.remove_old = false; % remove old selection of the same period
cfg.process.qcref.MinFiltPeriod = 50; % filter event period in minute
cfg.process.qcref.szFilt = 10; % filter even length in minute

%%% Split total/filtered %%%
cfg.process.split = struct();
cfg.process.split.reference = 'FLOW';
cfg.process.split.buffer = struct();
cfg.process.split.buffer.AC9 = [180, 60];
cfg.process.split.buffer.ACS279 = [180, 60];
cfg.process.split.buffer.ACS091 = [180, 60];
cfg.process.split.buffer.ACS111 = [180, 60];
cfg.process.split.buffer.ACS007 = [180, 60];
cfg.process.split.buffer.ACS057 = [180, 60];
cfg.process.split.buffer.BB3 = [400, 220];
% cfg.process.split.buffer.LISST = [540, 360];
% cfg.process.split.buffer.WSCD859 = [310, 20];
cfg.process.split.buffer.WSCD1082P = [310, 20];
% cfg.process.split.buffer.ALFA = [180, 30];
% cfg.process.split.skip = {'FLOW', 'TSG','WSCD1082P','PAR','ACS279'};
cfg.process.split.skip = {'FLOW', 'TSG','PAR'};

%%% Binning %%%
cfg.process.bin = struct('bin_size', struct());
cfg.process.bin.prctile_detection = [2.5, 97.5];
% if StepQC with ACS: prctile_average = [2.5, 97.5]; otherwise prctile_average = [5, 75];
cfg.process.bin.prctile_average = [2.5, 97.5];
% Bin mode does not affect the outcome of the data but just the way the data is presented to the computer
% cfg.process.bin.mode = 'OneShot'; % Faster for small dataset fiting in the memory of the computer
cfg.process.bin.mode = 'ByDay'; % Slightly slower but can handle a lot more data at once as it will be binned one day at a time
cfg.process.bin.bin_size.FLOW = 1;
cfg.process.bin.bin_size.AC9 = 1;
cfg.process.bin.bin_size.ACS007 = 1;
cfg.process.bin.bin_size.ACS057 = 1;
cfg.process.bin.bin_size.ACS091 = 1;
cfg.process.bin.bin_size.ACS279 = 1;
cfg.process.bin.bin_size.ACS111 = 1;
cfg.process.bin.bin_size.BB3 = 1;
cfg.process.bin.bin_size.PAR = 1;
% cfg.process.bin.bin_size.LISST = 10;
% cfg.process.bin.bin_size.WSCD859 = 1;
cfg.process.bin.bin_size.WSCD1082P = 1;
cfg.process.bin.bin_size.TSG = 1; % TSG does not need to be binned in most cases
% cfg.process.bin.bin_size.ALFA = 10;
cfg.process.bin.skip = {'TSG'};

%%% Automatically flagging %%%
cfg.process.flag = struct();
cfg.process.flag.skip = {'FLOW', 'TSG', 'AC9', 'ACS007', 'ACS057','ACS091', 'ACS111', 'ACS279', 'BB3', 'LISST', 'WSCD859', 'WSCD1082P', 'ALFA','PAR'};
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
cfg.process.qc.Saturation_Threshold_bb = 4100; % (counts)
  
%%% Manually QC %%%
cfg.process.qc.mode = 'ui';
cfg.process.qc.qc_once_for_all = false;  % true = QC all variables | false = QC variables separately);
cfg.process.qc.remove_old = false; % remove old selection of the same period
cfg.process.qc.global = struct();
cfg.process.qc.global.active = false;
cfg.process.qc.global.view = {'ACS279','BB3','PAR'};
cfg.process.qc.global.apply = {'ACS279','BB3','PAR'};
cfg.process.qc.specific = struct();
cfg.process.qc.specific.active = true;
cfg.process.qc.specific.run = {'ACS279','BB3','PAR'};
  
%%% Calibrate %%%
cfg.process.calibrate = struct();
cfg.process.calibrate.ACS007 = struct('compute_dissolved', false, ...
                                  'interpolation_method', 'CDOM', ...
                                  'CDOM_source', 'WSCD1082P', ... % Pay attention to serial number (must merge WSCD products first)
                                  'FLOW_source', 'FLOW', ...
                                  'di_method', 'normal', ... % best_di normal
                                  'compute_ad_aphi', false); % VERY SLOW: compute ad and aphi from Zheng and Stramski 2013                                
cfg.process.calibrate.AC9 = struct('compute_dissolved', false, ...
                                  'interpolation_method', 'CDOM', ...
                                  'CDOM_source', 'WSCD1082P',... % Pay attention to serial number (must merge WSCD products first)
                                  'FLOW_source', 'FLOW', ...
                                  'di_method', 'normal', ... % best_di normal
                                  'compute_ad_aphi', false); % VERY SLOW: compute ad and aphi from Zheng and Stramski 2013
cfg.process.calibrate.ACS057 = struct('compute_dissolved', false, ...
                                  'interpolation_method', 'CDOM', ...
                                  'CDOM_source', 'WSCD1082P',... % Pay attention to serial number (must merge WSCD products first)
                                  'FLOW_source', 'FLOW', ...
                                  'di_method', 'normal', ... % best_di normal
                                  'compute_ad_aphi', false); % VERY SLOW: compute ad and aphi from Zheng and Stramski 2013
cfg.process.calibrate.ACS091 = struct('compute_dissolved', false, ...
                                  'interpolation_method', 'CDOM', ...
                                  'CDOM_source', 'WSCD1082P',... % Pay attention to serial number (must merge WSCD products first)
                                  'FLOW_source', 'FLOW', ...
                                  'di_method', 'normal', ... % best_di normal
                                  'compute_ad_aphi', false); % VERY SLOW: compute ad and aphi from Zheng and Stramski 2013
cfg.process.calibrate.ACS111 = struct('compute_dissolved', false, ...
                                  'interpolation_method', 'CDOM', ...
                                  'CDOM_source', 'WSCD1082P', ... % Pay attention to serial number (must merge WSCD products first)
                                  'FLOW_source', 'FLOW', ...
                                  'di_method', 'normal', ... % best_di normal
                                  'compute_ad_aphi', false); % VERY SLOW: compute ad and aphi from Zheng and Stramski 2013
cfg.process.calibrate.ACS279 = struct('compute_dissolved', false,...
                                  'interpolation_method', 'CDOM', ...
                                  'CDOM_source', 'WSCD1082P',... % Pay attention to serial number (must merge WSCD products first)
                                  'FLOW_source', 'FLOW', ...
                                  'di_method', 'normal', ... % best_di normal
                                  'compute_ad_aphi', false); % VERY SLOW: compute ad and aphi from Zheng and Stramski 2013
cfg.process.calibrate.BB3 = struct('compute_dissolved', true, ...
                                  'TSG_source', 'TSG', ...
                                  'FLOW_source', 'FLOW', ...
                                  'di_method', 'SW_scattering', ... % interpolate constant SW_scattering
                                  'filt_method', 'exponential_fit'); % 25percentil exponential_fit
cfg.process.calibrate.skip = {'FLOW', 'TSG'};

%%% Write %%%
cfg.process.write = struct();
cfg.process.write.mode = 'One day one file'; % 'One file' or 'One day one file'
cfg.process.write.skip = {'TSG'};

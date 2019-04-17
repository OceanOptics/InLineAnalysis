% NAAMES 3 Configuration file
% author: Nils
% created: Nov 5, 2018

cfg = struct('meta', struct(), 'instruments', struct(), 'process', struct());

%%%%%%%%%%%%%%
%% METADATA %%
%%%%%%%%%%%%%%
cfg.meta.investigators = 'Emmanuel_Boss,Nils_Haentjens';
cfg.meta.affiliations = 'University_of_Maine';
cfg.meta.emails = 'emmanuel.boss@maine.edu';
cfg.meta.experiment = 'NAAMES03';
cfg.meta.cruise = 'AT38';
cfg.meta.station = 'NA';
cfg.meta.documents = 'NA';
cfg.meta.calibration_files = 'NA';
cfg.meta.data_type = 'flow_thru';
cfg.meta.data_status = 'final';
cfg.meta.measurement_depth = 5;


%%%%%%%%%%%%%%%%%
%% INSTRUMENTS %%
%%%%%%%%%%%%%%%%%
if exist('PATH_ROOT', 'var'); fprintf('Using global PATH_ROOT\n');
else; PATH_ROOT = '/Users/nils/Data/NAAMES/NAAMES3/'; end

%%% TSG + GPS %%%
cfg.instruments.TSG = struct();
cfg.instruments.TSG.model = 'TSG';
cfg.instruments.TSG.boat = 'Atlantis';
cfg.instruments.TSG.path = struct('raw',  [PATH_ROOT 'Underway/raw/'],...
                                  'wk',   [PATH_ROOT 'Underway/mat/'],...
                                  'prod', [PATH_ROOT 'Underway/'],...
                                  'ui', [PATH_ROOT 'Underway/ui/']);
cfg.instruments.TSG.view = struct('varname', 't');

%%% FTH (FlowControl) %%%
cfg.instruments.FTH = struct();
cfg.instruments.FTH.model = 'FTH';
cfg.instruments.FTH.logger = 'FlowControl';
cfg.instruments.FTH.path = struct('raw',  [PATH_ROOT 'FlowControl/raw/'],...
                                  'wk',   [PATH_ROOT 'FlowControl/mat/'],...
                                  'prod',   [PATH_ROOT 'FlowControl/'],...
                                  'ui', [PATH_ROOT 'FlowControl/ui/']);
cfg.instruments.FTH.view = struct('varname', 'swt');

%%% ACS 015 %%%
cfg.instruments.ACS = struct();
cfg.instruments.ACS.model = 'ACS';
cfg.instruments.ACS.sn = '015';
cfg.instruments.ACS.logger = 'Compass_2.1rc_scheduled_bin';
cfg.instruments.ACS.device_file = [PATH_ROOT 'ACS/acs015.dev'];
cfg.instruments.ACS.path = struct('raw',  [PATH_ROOT 'ACS/raw/'],...
                                  'di',  [PATH_ROOT 'ACS/raw_di/'],...
                                  'wk',   [PATH_ROOT 'ACS/wk/'],...
                                  'prod', [PATH_ROOT 'ACS/'],...
                                  'ui', [PATH_ROOT 'ACS/ui/']);
cfg.instruments.ACS.view = struct('varname', 'a', 'varcol', 40);

%%% BB3 %%%
cfg.instruments.BB3 = struct();
cfg.instruments.BB3.model = 'BB';
cfg.instruments.BB3.sn = '349';
cfg.instruments.BB3.ila_prefix = 'BB3';
cfg.instruments.BB3.logger = 'Inlinino';
cfg.instruments.BB3.lambda = [470,532,660];
cfg.instruments.BB3.theta = 124;
cfg.instruments.BB3.slope = [8.407E-06,4.624E-06,4.090E-06];
cfg.instruments.BB3.dark = [56,52,45];
cfg.instruments.BB3.path = struct('raw',  [PATH_ROOT 'BB3_WSCD/raw/'],...
                                  'di',  [PATH_ROOT 'BB3_WSCD/raw/'],...
                                  'wk',   [PATH_ROOT 'BB3_WSCD/mat/'],...
                                  'prod', [PATH_ROOT 'BB3_WSCD/'],...
                                  'ui', [PATH_ROOT 'BB3_WSCD/ui/']);
cfg.instruments.BB3.view = struct('varname', 'beta', 'varcol', 2);

%%% LISST %%%
cfg.instruments.LISST = struct();
cfg.instruments.LISST.model = 'LISST';
cfg.instruments.LISST.type = 'B'; 
cfg.instruments.LISST.sn = '1183';
cfg.instruments.LISST.logger = 'TeraTerm';
cfg.instruments.LISST.zsc = [2.203500e+001, 2.568500e+001, 2.503000e+001, 2.986000e+001, 2.842500e+001, 3.283000e+001, 3.077000e+001, 3.659500e+001, 2.978000e+001, 3.552000e+001, 3.198000e+001, 4.216000e+001, 3.916500e+001, 4.662500e+001, 3.974000e+001, 4.454000e+001, 4.403500e+001, 4.604500e+001, 4.430000e+001, 4.510500e+001, 4.719500e+001, 3.850000e+001, 5.373000e+001, 2.664000e+001, 3.180500e+001, 1.655500e+001, 2.205500e+001, 1.554000e+001, 1.422000e+001, 1.123000e+001, 8.780000e+000, 8.555000e+000, 1.515000e+003, 1.167900e+003, 6.410000e+001, 1.055150e+003, 7.700000e+001, 2.116600e+003, 1.807000e+003, 5.476500e+003];
cfg.instruments.LISST.dcal = [1.0000000e+000, 1.0038000e+000, 9.9360000e-001, 1.0027000e+000, 9.9720000e-001, 9.9570000e-001, 9.9030000e-001, 9.9430000e-001, 9.9290000e-001, 9.9000000e-001, 9.9290000e-001, 9.9300000e-001, 9.9150000e-001, 9.9300000e-001, 9.9230000e-001, 9.9090000e-001, 1.1032000e+000, 1.1123000e+000, 1.2430000e+000, 1.1562000e+000, 1.3273000e+000, 1.1999000e+000, 1.0740000e+000, 1.7489000e+000, 1.5382000e+000, 2.5109000e+000, 2.5468000e+000, 3.5504000e+000, 3.9338000e+000, 5.1747342e+000, 7.5143548e+000, 1.2528083e+001];
cfg.instruments.LISST.vcc = 48493;
cfg.instruments.LISST.inversion = 'spherical';
cfg.instruments.LISST.ds = [1.2500,1.4750,1.7405,2.0538,2.4235,2.8597,3.3744,3.9818,4.6986,5.5443,6.5423,7.7199,9.1095,10.7492,12.6841,14.9672,17.6613,20.8403,24.5916,29.0180,34.2413,40.4047,47.6776,56.2595,66.3863,78.3358,92.4362,109.0747,128.7082,151.8757,179.2133,211.4717,249.5366];
cfg.instruments.LISST.theta = [0.082, 0.096, 0.114, 0.134, 0.158, 0.187, 0.221, 0.260, 0.307, 0.362, 0.428, 0.505, 0.596, 0.703, 0.829, 0.979, 1.155, 1.363, 1.609, 1.898, 2.240, 2.643, 3.119, 3.681, 4.344, 5.126, 6.049, 7.138, 8.424, 9.941, 11.73, 13.84];
cfg.instruments.LISST.path = struct('raw',  [PATH_ROOT 'LISST/raw/'],...
                                  'wk',   [PATH_ROOT 'LISST/wk/'],...
                                  'prod', [PATH_ROOT 'LISST/'],...
                                  'ui', [PATH_ROOT 'LISST/ui/']);
cfg.instruments.LISST.view = struct('varname', 'beta', 'varcol', 15);

%%% WSCD %%%
cfg.instruments.WSCD = struct();
cfg.instruments.WSCD.model = 'CD';
cfg.instruments.WSCD.sn = '1299';
cfg.instruments.WSCD.ila_prefix = 'WSCD';
cfg.instruments.WSCD.logger = 'Inlinino';
cfg.instruments.WSCD.path = struct('raw',  [PATH_ROOT 'BB3_WSCD/raw/'],...
                                  'wk',   [PATH_ROOT 'BB3_WSCD/mat/'],...
                                  'prod', [PATH_ROOT 'BB3_WSCD/'],...
                                  'ui', [PATH_ROOT 'BB3_WSCD/ui/']);
cfg.instruments.WSCD.view = struct('varname', 'fdom');

%%%%%%%%%%%%%
%% PROCESS %%
%%%%%%%%%%%%%

%%% General parameters %%%
cfg.process.days2run = datenum('30-Aug-2017'):datenum('22-Sept-2017');
cfg.process.instruments2run = {'FTH', 'TSG', 'ACS', 'BB3', 'LISST', 'WSCD'};
cfg.process.write = true;
cfg.process.force_import = false;
cfg.process.parallel = Inf; % 0: disable parallel or Inf: as many thread available
% cfg.process.di = struct();
% cfg.process.di.skip = {'FTH', 'TSG', 'WSCD'};
% cfg.process.di.qc = struct('mode', 'load');
% cfg.process.di.bin = struct('bin_size', 30);

%%% Synchronization %%%
cfg.process.sync = struct();
cfg.process.sync.delay = struct();
cfg.process.sync.delay.FTH = 30;
cfg.process.sync.delay.ACS = 97;
cfg.process.sync.delay.BB3 = 42;
cfg.process.sync.delay.LISST = 20;
cfg.process.sync.delay.WSCD = 20;
cfg.process.sync.skip = {'TSG'};

%%% Strech %%%
cfg.proces.stretch = struct();
cfg.proces.stretch.delta = struct();
cfg.proces.stretch.delta.BB3 = 130;
cfg.proces.stretch.delta.WSCD = 130;
cfg.process.stretch.skip = {'FTH', 'TSG', 'ACS'};

%%% QC Reference (Flow Control/FTH) %%%
cfg.process.qcref = struct();
cfg.process.qcref.reference = 'FTH';
% cfg.process.qcref.view = 'BB3';
cfg.process.qcref.view = 'ACS';
cfg.process.qcref.mode = 'load'; % load or ui

%%% Split total/filtered %%%
cfg.process.split = struct();
cfg.process.split.reference = 'FTH';
cfg.process.split.buffer = struct();
cfg.process.split.buffer.ACS = [300, 20];
cfg.process.split.buffer.BB3 = [400, 200];
cfg.process.split.buffer.LISST = [400, 200];
cfg.process.split.buffer.WSCD = [210, 0];
cfg.process.split.skip = {'FTH', 'TSG'};

%%% Binning %%%
cfg.process.bin = struct('bin_size', struct());
% Directly read from instrument not from here
% cfg.process.bin.method = '4flag'; % Method to use to flag automatically
% cfg.process.bin.method = 'SB_IN_PRCTL'; % Faster method outputting only requirements for SeaBASS (not compatible with automatic flagging)
% cfg.process.bin.method = 'SB_ALL'; % Faster method outputting only requirements for SeaBASS (not compatible with automatic flagging)
cfg.process.bin.prctile_detection = [2.5, 97.5];
cfg.process.bin.prctile_average = [5, 75];
cfg.process.bin.bin_size.ACS = 1;
cfg.process.bin.bin_size.BB3 = 1;
cfg.process.bin.bin_size.LISST = 10;
cfg.process.bin.bin_size.WSCD = 1;
cfg.process.bin.skip = {'FTH', 'TSG'};

%%% Automatically flagging %%%
cfg.process.flag = struct();
cfg.process.flag.skip = {'FTH', 'TSG', 'ACS', 'BB3', 'LISST', 'WSCD'};
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
% Flag parameters specific to the ACS
cfg.process.flag.ACS = struct();
cfg.process.flag.ACS.abs_uncertainty = 0.004;
cfg.process.flag.ACS.rel_uncertainty = 0.0125;
cfg.process.flag.ACS.min_flag_n = 67;
cfg.process.flag.ACS.primary_varname = 'c';
cfg.process.flag.ACS.filt = struct();
cfg.process.flag.ACS.filt.abs_uncertainty =  0.001;
cfg.process.flag.ACS.filt.rel_uncertainty = 0;
cfg.process.flag.ACS.filt.min_flag_n = 20;
% Flag parameters specific to the BB3
cfg.process.flag.BB3 = struct();
cfg.process.flag.BB3.abs_uncertainty = 0.00001;
cfg.process.flag.BB3.rel_uncertainty = 0.02;
cfg.process.flag.BB3.min_flag_n =  1;
cfg.process.flag.BB3.primary_varname = 'beta';
% Flag parameters specific to the WSCD
cfg.process.flag.WSCD = struct();
cfg.process.flag.WSCD.abs_uncertainty = 0.0007;
cfg.process.flag.WSCD.rel_uncertainty = 0;
cfg.process.flag.WSCD.primary_varname = 'fdom';
  
%%% Manually QC %%%
cfg.process.qc = struct();
cfg.process.qc.mode = 'ui';
cfg.process.qc.global = struct();
cfg.process.qc.global.active = false;
cfg.process.qc.global.view = 'ACS';
cfg.process.qc.global.apply = {'ACS', 'BB3', 'WSCD'};
cfg.process.qc.specific = struct();
cfg.process.qc.specific.active = true;
cfg.process.qc.specific.run = {'TSG', 'ACS', 'BB3', 'LISST', 'WSCD'};
  
%%% Calibrate %%%
cfg.process.calibrate = struct();
cfg.process.calibrate.ACS = struct('compute_dissolved', true,...
                                  'interpolation_method', 'linear', ...
                                  'CDOM_source', 'WSCD',...
                                  'FTH_source', 'FTH');
cfg.process.calibrate.BB3 = struct('compute_dissolved', true,...
                                   'TSG_source', 'TSG',...
                                   'di_method', 'constant');
cfg.process.calibrate.skip = {'FTH', 'TSG', 'WSCD', 'ALFA'};

%%% Write %%%
cfg.process.write = struct();
cfg.process.write.mode = 'One day one file'; % 'One file' or 'One day one file'
cfg.process.write.skip = {'FTH'};

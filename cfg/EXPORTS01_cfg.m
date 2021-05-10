% EXPORTS Configuration file
% author: Nils
% created: Aug 16, 2018

cfg = struct('meta', struct(), 'instruments', struct(), 'process', struct());

%%%%%%%%%%%%%%
%% METADATA %%
%%%%%%%%%%%%%%
cfg.meta.investigators = 'Emmanuel_Boss,Nils_Haentjens';
cfg.meta.affiliations = 'University_of_Maine';
cfg.meta.emails = 'emmanuel.boss@maine.edu';
cfg.meta.experiment = 'EXPORTS';
cfg.meta.cruise = 'EXPORTSNP';
cfg.meta.station = 'NA';
cfg.meta.documents = 'NA';
cfg.meta.calibration_files = 'NA';
cfg.meta.data_type = 'flow_thru';
cfg.meta.data_status = 'preliminary';
cfg.meta.measurement_depth = 5;


%%%%%%%%%%%%%%%%%
%% INSTRUMENTS %%
%%%%%%%%%%%%%%%%%

PATH_ROOT = '/Volumes/Samsung_T5/Data/EXPORTS01/';

%%% TSG + GPS %%%
cfg.instruments.TSG = struct();
cfg.instruments.TSG.model = 'TSG';
cfg.instruments.TSG.boat = 'RRevelle';
cfg.instruments.TSG.logger = 'Matlab'; % TeraTerm Matlab
cfg.instruments.TSG.path = struct('raw',  [PATH_ROOT 'raw' filesep 'TSG' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'TSG' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'TSG' filesep]);
cfg.instruments.TSG.view = struct('varname', 't');

%%% FLOW (FlowControl) %%%
cfg.instruments.FLOW = struct();
cfg.instruments.FLOW.model = 'FTH';
cfg.instruments.FLOW.logger = 'FlowControl';
cfg.instruments.FLOW.path = struct('raw',  [PATH_ROOT 'raw' filesep 'FlowControl' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'FlowControl' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'FlowControl' filesep]);
cfg.instruments.FLOW.view = struct('varname', 'swt');

%%% ACS 301 %%% (Aug 20 to ...)
SN = '301';
cfg.instruments.(['ACS' SN]) = struct();
cfg.instruments.(['ACS' SN]).di = struct();
cfg.instruments.(['ACS' SN]).di.prefix = ['acs_' SN '_'];
cfg.instruments.(['ACS' SN]).di.postfix = '_DI';
cfg.instruments.(['ACS' SN]).model = 'ACS';
cfg.instruments.(['ACS' SN]).sn = SN;
cfg.instruments.(['ACS' SN]).logger = 'Compass_2.1rc_scheduled_bin'; % InlininoACScsv
cfg.instruments.(['ACS' SN]).device_file = [PATH_ROOT 'DeviceFiles' filesep 'acs301.dev'];
cfg.instruments.(['ACS' SN]).path = struct('raw',  [PATH_ROOT 'raw' filesep ['ACS' SN] filesep],...
                                    'di',  [PATH_ROOT 'raw' filesep ['ACS' SN] filesep 'DI' filesep],...
                                    'wk',   [PATH_ROOT 'wk' filesep ['ACS' SN] filesep],...
                                    'prod', [PATH_ROOT 'prod' filesep],...
                                    'ui', [PATH_ROOT 'ui' filesep ['ACS' SN] filesep]);
cfg.instruments.(['ACS' SN]).view = struct('varname', 'a', 'varcol', 40);

%%% ACS 298 %%% (Aug 11 to Aug 20)
SN = '298';
cfg.instruments.(['ACS' SN]) = struct();
cfg.instruments.(['ACS' SN]).di = struct();
cfg.instruments.(['ACS' SN]).di.prefix = ['acs_' SN '_'];
cfg.instruments.(['ACS' SN]).di.postfix = '_DI';
cfg.instruments.(['ACS' SN]).model = 'ACS';
cfg.instruments.(['ACS' SN]).sn = SN;
cfg.instruments.(['ACS' SN]).logger = 'Compass_2.1rc_scheduled_bin'; % InlininoACScsv
cfg.instruments.(['ACS' SN]).device_file = [PATH_ROOT 'DeviceFiles' filesep 'acs298.dev'];cfg.instruments.(['ACS' SN]).path = struct('raw',  [PATH_ROOT 'raw' filesep ['ACS' SN] filesep],...
                                    'di',  [PATH_ROOT 'raw' filesep ['ACS' SN] filesep 'DI' filesep],...
                                    'wk',   [PATH_ROOT 'wk' filesep ['ACS' SN] filesep],...
                                    'prod', [PATH_ROOT 'prod' filesep],...
                                    'ui', [PATH_ROOT 'ui' filesep ['ACS' SN] filesep]);
cfg.instruments.(['ACS' SN]).view = struct('varname', 'a', 'varcol', 40);

%%% BB3 %%%
SN = '349';
cfg.instruments.(['BB3' SN]) = struct();
cfg.instruments.(['BB3' SN]).di = struct();
cfg.instruments.(['BB3' SN]).di.prefix = ['BB3' SN '_'];
cfg.instruments.(['BB3' SN]).di.postfix = 'DI';
cfg.instruments.(['BB3' SN]).model = 'BB';
cfg.instruments.(['BB3' SN]).sn = SN;
cfg.instruments.(['BB3' SN]).ila_prefix = 'Inlinino'; % ['BB3' SN]
cfg.instruments.(['BB3' SN]).logger = 'InlininoBB3'; % InlininoBB3SN
cfg.instruments.(['BB3' SN]).lambda = [470,532,650];
cfg.instruments.(['BB3' SN]).theta = 124;
cfg.instruments.(['BB3' SN]).slope = [8.407E-06,4.624E-06,4.090E-06];
cfg.instruments.(['BB3' SN]).dark = [56,52,45];
cfg.instruments.(['BB3' SN]).path = struct('raw',  [PATH_ROOT 'raw' filesep ['BB3' SN] filesep],...
                                  'di',  [PATH_ROOT 'raw' filesep ['BB3' SN] filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep ['BB3' SN] filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep ['BB3' SN] filesep]);
cfg.instruments.(['BB3' SN]).view = struct('varname', 'beta', 'varcol', 2);

%%% WSCD %%% (Aug 17 to ...)
SN = '1299';
cfg.instruments.(['WSCD' SN]) = struct();
cfg.instruments.(['WSCD' SN]).sn = SN;
cfg.instruments.(['WSCD' SN]).model = 'CD';
cfg.instruments.(['WSCD' SN]).ila_prefix = 'Inlinino'; % ['WSCD' SN]
cfg.instruments.(['WSCD' SN]).logger = 'InlininoWSCD'; % InlininoWSCDSN
cfg.instruments.(['WSCD' SN]).slope = 62;
cfg.instruments.(['WSCD' SN]).dark = 0.059;
cfg.instruments.(['WSCD' SN]).path = struct('raw',  [PATH_ROOT 'raw' filesep ['WSCD' SN] filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep ['WSCD' SN] filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep ['WSCD' SN] filesep]);
cfg.instruments.(['WSCD' SN]).view = struct('varname', 'fdom');

%%% WSCD %%% (Aug 11 to 17)
SN = '859';
cfg.instruments.(['WSCD' SN]) = struct();
cfg.instruments.(['WSCD' SN]).sn = SN;
cfg.instruments.(['WSCD' SN]).model = 'CD';
cfg.instruments.(['WSCD' SN]).ila_prefix = 'Inlinino'; % ['WSCD' SN]
cfg.instruments.(['WSCD' SN]).logger = 'InlininoWSCD'; % InlininoWSCDSN
cfg.instruments.(['WSCD' SN]).slope = 0.0909;
cfg.instruments.(['WSCD' SN]).dark = 44;
cfg.instruments.(['WSCD' SN]).path = struct('raw',  [PATH_ROOT 'raw' filesep ['WSCD' SN] filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep ['WSCD' SN] filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep ['WSCD' SN] filesep]);
cfg.instruments.(['WSCD' SN]).view = struct('varname', 'fdom');


%%% LISST %%%
SN = '1183';
cfg.instruments.LISST = struct();
cfg.instruments.LISST.model = 'LISST';
cfg.instruments.LISST.type = 'B'; 
cfg.instruments.LISST.sn = SN;
cfg.instruments.LISST.logger = 'TeraTerm';
cfg.instruments.LISST.zsc = [2.203500e+001, 2.568500e+001, 2.503000e+001, 2.986000e+001, 2.842500e+001, 3.283000e+001, 3.077000e+001, 3.659500e+001, 2.978000e+001, 3.552000e+001, 3.198000e+001, 4.216000e+001, 3.916500e+001, 4.662500e+001, 3.974000e+001, 4.454000e+001, 4.403500e+001, 4.604500e+001, 4.430000e+001, 4.510500e+001, 4.719500e+001, 3.850000e+001, 5.373000e+001, 2.664000e+001, 3.180500e+001, 1.655500e+001, 2.205500e+001, 1.554000e+001, 1.422000e+001, 1.123000e+001, 8.780000e+000, 8.555000e+000, 1.515000e+003, 1.167900e+003, 6.410000e+001, 1.055150e+003, 7.700000e+001, 2.116600e+003, 1.807000e+003, 5.476500e+003];
% Original dcal file (ring area)
% cfg.instruments.LISST.dcal = [1.0000000e+000, 1.0038000e+000, 9.9360000e-001, 1.0027000e+000, 9.9720000e-001, 9.9570000e-001, 9.9030000e-001, 9.9430000e-001, 9.9290000e-001, 9.9000000e-001, 9.9290000e-001, 9.9300000e-001, 9.9150000e-001, 9.9300000e-001, 9.9230000e-001, 9.9090000e-001, 1.1032000e+000, 1.1123000e+000, 1.2430000e+000, 1.1562000e+000, 1.3273000e+000, 1.1999000e+000, 1.0740000e+000, 1.7489000e+000, 1.5382000e+000, 2.5109000e+000, 2.5468000e+000, 3.5504000e+000, 3.9338000e+000, 5.1747342e+000, 7.5143548e+000, 1.2528083e+001];
% dcal adhoc from email of Wayne Slade Dec 8, 2017
cfg.instruments.LISST.dcal = [ 1.0179083e+00	   9.9213489e-01	   1.0108161e+00	   9.9492883e-01	   1.0043707e+00	   9.9891840e-01	   9.9859055e-01	   1.0042049e+00	   1.0000763e+00	   9.9889997e-01	   1.0009497e+00	   1.0004019e+00	   1.0011130e+00	   1.0004677e+00	   1.0213554e+00	   9.9990262e-01	   1.1115630e+00	   1.1206668e+00	   1.2493699e+00	   1.1643199e+00	   1.3355657e+00	   1.2090892e+00	   1.0781540e+00	   1.7620752e+00	   1.5508563e+00	   2.5304119e+00	   2.5638592e+00	   3.5757212e+00	   3.9631987e+00	   5.0166411e+00	   5.6381118e+00	   8.6881539e+00];
% customize dcal of 1183 on Jan 16, 2018 due to bump in VSF at ring 30
cfg.instruments.LISST.dcal(30) = 2.6; % Good for low values but bad for high values ??
cfg.instruments.LISST.vcc = 48493;
cfg.instruments.LISST.inversion = 'spherical';
cfg.instruments.LISST.ds = [1.2500,1.4750,1.7405,2.0538,2.4235,2.8597,3.3744,3.9818,4.6986,5.5443,6.5423,7.7199,9.1095,10.7492,12.6841,14.9672,17.6613,20.8403,24.5916,29.0180,34.2413,40.4047,47.6776,56.2595,66.3863,78.3358,92.4362,109.0747,128.7082,151.8757,179.2133,211.4717,249.5366];
cfg.instruments.LISST.theta = [0.082, 0.096, 0.114, 0.134, 0.158, 0.187, 0.221, 0.260, 0.307, 0.362, 0.428, 0.505, 0.596, 0.703, 0.829, 0.979, 1.155, 1.363, 1.609, 1.898, 2.240, 2.643, 3.119, 3.681, 4.344, 5.126, 6.049, 7.138, 8.424, 9.941, 11.73, 13.84];
cfg.instruments.LISST.path = struct('raw',  [PATH_ROOT 'raw' filesep ['LISST' SN] filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep ['LISST' SN] filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep ['LISST' SN] filesep]);
cfg.instruments.LISST.view = struct('varname', 'beta', 'varcol', 15);

%%% ALFA %%%
SN = '011';
cfg.instruments.ALFA = struct();
cfg.instruments.ALFA.model = 'ALFA';
cfg.instruments.ALFA.sn = SN;
cfg.instruments.ALFA.logger = 'ALFA_LabView_m';
cfg.instruments.ALFA.path = struct('raw',  [PATH_ROOT 'raw' filesep ['ALFA' SN] filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep ['ALFA' SN] filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep ['ALFA' SN] filesep]);
cfg.instruments.ALFA.view = struct('varname', 'FvFm');

%%%%%%%%%%%%%
%% PROCESS %%
%%%%%%%%%%%%%

%%% General parameters %%%
cfg.process.days2run = datenum(2018,08,11):datenum(2018,09,12);
cfg.process.instruments2run = {'FLOW', 'TSG', 'ACS301', 'BB3349', 'LISST', 'WSCD1299', 'ALFA'};
cfg.process.write = true;
cfg.process.force_import = false;
cfg.process.parallel = Inf; % 0: disable parallel or Inf: as many thread available
cfg.process.di = struct();
cfg.process.di.skip = {'FLOW', 'TSG', 'LISST', 'WSCD859', 'WSCD1299', 'ALFA'};
cfg.process.di.qc = struct('mode', 'load');
cfg.process.di.bin = struct('bin_size', 30);

%%% Synchronization %%%
cfg.process.sync = struct();
cfg.process.sync.delay = struct();
cfg.process.sync.delay.FLOW = 0;
cfg.process.sync.delay.ACS298 = 55;
cfg.process.sync.delay.ACS301 = 66;
cfg.process.sync.delay.BB3349 = 0;
cfg.process.sync.delay.LISST = 1;
cfg.process.sync.delay.WSCD859 = 5;
cfg.process.sync.delay.WSCD1299 = 5;
cfg.process.sync.delay.ALFA = 15;
cfg.process.sync.skip = {'TSG'};

%%% QC Reference (Flow Control/FLOW) %%%
cfg.process.qcref = struct();
cfg.process.qcref.reference = 'FLOW';
cfg.process.qcref.view = 'BB3349';
cfg.process.qcref.mode = 'load'; % load or ui
cfg.process.qcref.MinFiltPeriod = 50; % filter even period in minute

%%% Split total/filtered %%%
cfg.process.split = struct();
cfg.process.split.reference = 'FLOW';
cfg.process.split.buffer = struct();
cfg.process.split.buffer.ACS298 = [180, 60];
cfg.process.split.buffer.ACS301 = [180, 60];
cfg.process.split.buffer.BB3349 = [400, 220];
cfg.process.split.buffer.LISST = [540, 360];
cfg.process.split.buffer.WSCD859 = [310, 20];
cfg.process.split.buffer.WSCD1299 = [310, 20];
cfg.process.split.buffer.ALFA = [180, 30];
cfg.process.split.skip = {'FLOW', 'TSG'};

%%% Binning %%%
cfg.process.bin = struct('bin_size', struct());
cfg.process.bin.prctile_detection = [2.5, 97.5];
% if StepQC with ACS: prctile_average = [2.5, 97.5]; otherwise prctile_average = [5, 75];
cfg.process.bin.prctile_average = [2.5, 97.5];
% Bin mode does not affect the outcome of the data but just the way the data is presented to the computer
% cfg.process.bin.mode = 'OneShot'; % Faster for small dataset fiting in the memory of the computer
cfg.process.bin.mode = 'ByDay'; % Slightly slower but can handle a lot more data at once as it will be binned one day at a time
cfg.process.bin.bin_size.ACS298 = 1;
cfg.process.bin.bin_size.ACS301 = 1;
cfg.process.bin.bin_size.BB3349 = 1;
cfg.process.bin.bin_size.LISST = 10;
cfg.process.bin.bin_size.WSCD859 = 1;
cfg.process.bin.bin_size.WSCD1299 = 1;
cfg.process.bin.bin_size.TSG = 1; % TSG does not need to be binned in most cases
cfg.process.bin.bin_size.ALFA = 10;
cfg.process.bin.bin_size.FLOW = 1;
cfg.process.bin.skip = {};

%%% Automatically flagging %%%
cfg.process.flag = struct();
cfg.process.flag.skip = {'FLOW', 'TSG', 'ACS301', 'ACS298', 'BB3349', 'LISST', 'WSCD859', 'WSCD1299', 'ALFA'};
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
cfg.process.qc.Saturation_Threshold_bb = 4000; % (counts)

%%% Manually QC %%%
cfg.process.qc.mode = 'ui';
cfg.process.qc.global = struct();
cfg.process.qc.global.active = false;
cfg.process.qc.global.view = 'ACS298';
cfg.process.qc.global.apply = {'ACS298', 'ACS301', 'BB3349', 'LISST', 'WSCD859', 'WSCD1299'};
cfg.process.qc.specific = struct();
cfg.process.qc.specific.active = true;
cfg.process.qc.specific.run = {'TSG', 'ACS298', 'ACS301', 'BB3349', 'LISST', 'WSCD859', 'WSCD1299', 'ALFA'};
  
%%% Calibrate %%%
cfg.process.calibrate = struct();
cfg.process.calibrate.ACS298 = struct('compute_dissolved', true,...
                                  'interpolation_method', 'linear', ...
                                  'CDOM_source', 'WSCD859',... % Pay attention to serial number (could merge WSCD products first)
                                  'FLOW_source', 'FLOW');
cfg.process.calibrate.ACS301 = struct('compute_dissolved', true,...
                                  'interpolation_method', 'linear', ...
                                  'CDOM_source', 'WSCD1299',... % Pay attention to serial number (could merge WSCD products first)
                                  'FLOW_source', 'FLOW');
cfg.process.calibrate.BB3349 = struct('compute_dissolved', true,...
                                   'TSG_source', 'TSG',...
                                   'di_method', 'constant');
cfg.process.calibrate.skip = {'FLOW', 'TSG', 'WSCD859', 'WSCD1299', 'ALFA'};

%%% Write %%%
cfg.process.write = struct();
cfg.process.write.mode = 'One day one file'; % 'One file' or 'One day one file'
cfg.process.write.skip = {};

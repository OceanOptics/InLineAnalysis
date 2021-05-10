% EXPORTS Configuration file
% author: Nils
% created: Aug 16, 2018

cfg = struct('meta', struct(), 'instruments', struct(), 'process', struct());

%%%%%%%%%%%%%%
%% METADATA %%
%%%%%%%%%%%%%%
cfg.meta.investigators = 'Emmanuel_Boss, Lee Karp-Boss';
cfg.meta.affiliations = 'University_of_Maine';
cfg.meta.emails = 'emmanuel.boss@maine.edu';
cfg.meta.experiment = 'EXPORTS';
cfg.meta.cruise = 'EXPORTSNA';
cfg.meta.station = 'NA';
cfg.meta.documents = 'NA';
cfg.meta.calibration_files = 'NA';
cfg.meta.data_type = 'flow_thru';
cfg.meta.data_status = 'preliminary';
cfg.meta.measurement_depth = 5;


%%%%%%%%%%%%%%%%%
%% INSTRUMENTS %%
%%%%%%%%%%%%%%%%%

PATH_ROOT = '/Users/emmanuel.boss/Desktop/InLine analysis/Data/EXPORTS02';

%%% TSG + GPS %%%
cfg.instruments.TSG = struct();
cfg.instruments.TSG.model = 'TSG';
cfg.instruments.TSG.boat = 'JCook';
cfg.instruments.TSG.path = struct('raw',  [PATH_ROOT 'raw'],...
                                  'wk',   [PATH_ROOT 'wk'],...
                                  'prod', [PATH_ROOT 'prod'],...
                                  'ui', [PATH_ROOT 'ui']);
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
SN = '91';
cfg.instruments.(['ACS' SN]) = struct();
cfg.instruments.(['ACS' SN]).di = struct();
cfg.instruments.(['ACS' SN]).di.prefix = ['DIW_ACS' SN '_'];
cfg.instruments.(['ACS' SN]).di.postfix = '';
cfg.instruments.(['ACS' SN]).model = 'ACS';
cfg.instruments.(['ACS' SN]).sn = SN;
cfg.instruments.(['ACS' SN]).logger = 'InlininoACScsv';
cfg.instruments.(['ACS' SN]).device_file = [PATH_ROOT filesep 'DeviceFiles' filesep 'acs091_20200217.dev'];
cfg.instruments.(['ACS' SN]).path = struct('raw',  [PATH_ROOT 'raw' filesep ['ACS' SN] filesep],...
                                    'di',  [PATH_ROOT 'raw' filesep ['ACS' SN] filesep 'DI' filesep],...
                                    'wk',   [PATH_ROOT 'wk' filesep ['ACS' SN] filesep],...
                                    'prod', [PATH_ROOT 'prod' filesep],...
                                    'ui', [PATH_ROOT 'ui' filesep ['ACS' SN] filesep]);
cfg.instruments.(['ACS' SN]).view = struct('varname', 'a', 'varcol', 40);

% %%% ACS 298 %%% (Aug 11 to Aug 20)
% cfg.instruments.ACS298 = struct();
% cfg.instruments.ACS298.model = 'ACS';
% cfg.instruments.ACS298.sn = '298';
% cfg.instruments.ACS298.ila_prefix = 'ACS';
% cfg.instruments.ACS298.logger = 'Compass_2.1rc_scheduled_bin';
% cfg.instruments.ACS298.device_file = [PATH_ROOT '../DeviceFiles/ACS298_20171215/acs298.dev'];
% cfg.instruments.ACS298.logger = 'Compass_2.1rc_scheduled';
% cfg.instruments.ACS298.lambda_reference = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
% cfg.instruments.ACS298.lambda_a = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
% cfg.instruments.ACS298.lambda_c = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
% cfg.instruments.ACS298.path = struct('raw',  [PATH_ROOT 'raw/ACS298/'],...
%                                   'di',  [PATH_ROOT 'raw/ACS298/DI/'],...
%                                   'wk',   [PATH_ROOT 'wk/ACS298/'],...
%                                   'prod', [PATH_ROOT 'prod/'],...
%                                   'ui', [PATH_ROOT 'ui/ACS298/']);
% cfg.instruments.ACS298.view = struct('varname', 'a', 'varcol', 40);

%%% BB3 %%%
cfg.instruments.BB3 = struct();
cfg.instruments.BB3.model = 'BB';
cfg.instruments.BB3.sn = '1052';
cfg.instruments.BB3.ila_prefix = 'BB3';
cfg.instruments.BB3.logger = 'InlininoBB3';
cfg.instruments.BB3.lambda = [470,532,660];
cfg.instruments.BB3.theta = 120;
cfg.instruments.BB3.slope = [1.25E-05,1.04E-06,5.49E-06];
cfg.instruments.BB3.dark = [56,52,45];
cfg.instruments.BB3.path = struct('raw',  [PATH_ROOT 'raw/BB31052/'],...
                                  'di',  [PATH_ROOT 'raw/BB31052/'],...
                                  'wk',   [PATH_ROOT 'wk/BB31052/'],...
                                  'prod', [PATH_ROOT 'prod/'],...
                                  'ui', [PATH_ROOT 'ui/BB31052/']);
cfg.instruments.BB3.view = struct('varname', 'beta', 'varcol', 2);

%%% SPCD %%%
cfg.instruments.SPCD = struct();
cfg.instruments.SPCD.model = 'CD';
cfg.instruments.SPCD.sn = '';
cfg.instruments.SPCD.ila_prefix = 'SPCD';
cfg.instruments.SPCD.logger = 'InlininoSPCD';
cfg.instruments.SPCD.slope = NaN;
cfg.instruments.SPCD.dark = NaN;
cfg.instruments.SPCD.path = struct('raw',  [PATH_ROOT 'raw/SPCD/'],...
                                  'wk',   [PATH_ROOT 'wk/SPCD/'],...
                                  'prod', [PATH_ROOT 'prod/'],...
                                  'ui', [PATH_ROOT 'ui/SPCD/']);
cfg.instruments.SPCD.view = struct('varname', 'fdom');

% %%% WSCD %%% (Aug 17 to ...)
% cfg.instruments.WSCD1299 = struct();
% cfg.instruments.WSCD1299.model = 'CD';
% cfg.instruments.WSCD1299.sn = '1299';
% cfg.instruments.WSCD1299.ila_prefix = 'WSCD';
% cfg.instruments.WSCD1299.logger = 'InlininoWSCD';
% cfg.instruments.WSCD1299.slope = NaN;
% cfg.instruments.WSCD1299.dark = NaN;
% cfg.instruments.WSCD1299.path = struct('raw',  [PATH_ROOT 'raw/WSCD1299/'],...
%                                   'wk',   [PATH_ROOT 'wk/WSCD1299/'],...
%                                   'prod', [PATH_ROOT 'prod/'],...
%                                   'ui', [PATH_ROOT 'ui/WSCD1299/']);
% cfg.instruments.WSCD1299.view = struct('varname', 'fdom');
% 
% %%% WSCD %%% (Aug 11 to 17)
% cfg.instruments.WSCD859 = struct();
% cfg.instruments.WSCD859.model = 'CD';
% cfg.instruments.WSCD859.sn = '859';
% cfg.instruments.WSCD859.ila_prefix = 'WSCD';
% cfg.instruments.WSCD859.logger = 'InlininoWSCD';
% cfg.instruments.WSCD859.slope = 0.0909;
% cfg.instruments.WSCD859.dark = 44;
% cfg.instruments.WSCD859.path = struct('raw',  [PATH_ROOT 'raw/WSCD859/'],...
%                                   'wk',   [PATH_ROOT 'wk/WSCD859/'],...
%                                   'prod', [PATH_ROOT 'prod/'],...
%                                   'ui', [PATH_ROOT 'ui/WSCD859/']);
% cfg.instruments.WSCD859.view = struct('varname', 'fdom');

%%% LISST %%%
cfg.instruments.LISST = struct();
cfg.instruments.LISST.model = 'LISST';
cfg.instruments.LISST.type = 'B'; 
cfg.instruments.LISST.sn = '1183';
cfg.instruments.LISST.logger = 'InlininoLISST';
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
cfg.instruments.LISST.path = struct('raw',  [PATH_ROOT 'raw/LISST1183/'],...
                                  'wk',   [PATH_ROOT 'wk/LISST1183/'],...
                                  'prod', [PATH_ROOT 'prod/'],...
                                  'ui', [PATH_ROOT 'ui/LISST1183/']);
cfg.instruments.LISST.view = struct('varname', 'beta', 'varcol', 15);


%%%%%%%%%%%%%
%% PROCESS %%
%%%%%%%%%%%%%

%%% General parameters %%%
cfg.process.days2run = datenum(2021,5,2):datenum(2021,5,5);
cfg.process.instruments2run = {'FLOW', 'TSG', 'ACS91', 'BB3', 'LISST', 'SPCD','Hyperbb'};
cfg.process.write = true;
cfg.process.force_import = false;
cfg.process.parallel = Inf; % 0: disable parallel or Inf: as many thread available
cfg.process.di = struct();
cfg.process.di.skip = {'FLOW', 'TSG', 'LISST', 'SPCD','Hyperbb'};
cfg.process.di.qc = struct('mode', 'load');
cfg.process.di.bin = struct('bin_size', 30);

%%% Synchronization %%%
cfg.process.sync = struct();
cfg.process.sync.delay = struct();
cfg.process.sync.delay.FLOW = 0;
cfg.process.sync.delay.ACS91 = 66;
cfg.process.sync.delay.BB3 = 0;
cfg.process.sync.delay.LISST = 1;
cfg.process.sync.delay.SPCD = 5;
cfg.process.sync.skip = {'TSG'};

%%% QC Reference (Flow Control/FLOW) %%%
cfg.process.qcref = struct();
cfg.process.qcref.reference = 'FLOW';
cfg.process.qcref.view = 'BB3';
cfg.process.qcref.mode = 'load'; % load or ui
cfg.process.qcref.MinFiltPeriod = 50; % filter even period in minute

%%% Split total/filtered %%%
cfg.process.split = struct();
cfg.process.split.reference = 'FLOW';
cfg.process.split.buffer = struct();
cfg.process.split.buffer.ACS91 = [180, 30];
cfg.process.split.buffer.BB3 = [420, 220];
cfg.process.split.buffer.LISST = [540, 360];
cfg.process.split.buffer.SPCD = [310, 20];
cfg.process.split.skip = {'FLOW', 'TSG'};

%%% Binning %%%
cfg.process.bin = struct('bin_size', struct());
cfg.process.bin.prctile_detection = [2.5, 97.5];
% if RawAutoQC with ACS: prctile_average = [2.5, 97.5]; otherwise prctile_average = [5, 75];
cfg.process.bin.prctile_average = [2.5, 97.5];
% Bin mode does not affect the outcome of the data but just the way the data is presented to the computer
% cfg.process.bin.mode = 'OneShot'; % Faster for small dataset fiting in the memory of the computer
cfg.process.bin.mode = 'ByDay'; % Slightly slower but can handle a lot more data at once as it will be binned one day at a time
cfg.process.bin.bin_size.ACS298 = 1;
cfg.process.bin.bin_size.ACS91 = 1;
cfg.process.bin.bin_size.BB3 = 1;
cfg.process.bin.bin_size.LISST = 10;
cfg.process.bin.bin_size.SPCD = 1;
cfg.process.bin.bin_size.TSG = 1; % TSG does not need to be binned in most cases
cfg.process.bin.skip = {};

%%% Automatically flagging %%%
cfg.process.flag = struct();
cfg.process.flag.skip = {'FLOW', 'TSG', 'ACS91',  'BB3', 'LISST', 'SPCD'};
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
cfg.process.qc = struct();
cfg.process.qc.mode = 'ui';
cfg.process.qc.global = struct();
cfg.process.qc.global.active = false;
cfg.process.qc.global.view = 'ACS91';
cfg.process.qc.global.apply = {'ACS91', 'BB3', 'LISST', 'SPCD'};
cfg.process.qc.specific = struct();
cfg.process.qc.specific.active = true;
cfg.process.qc.specific.run = {'TSG', 'ACS91', 'BB3', 'LISST', 'SPCD'};
  
%%% Calibrate %%%
cfg.process.calibrate = struct();
cfg.process.calibrate.ACS298 = struct('compute_dissolved', true,...
                                  'interpolation_method', 'linear', ...
                                  'CDOM_source', 'SPCD',... 
                                  'FLOW_source', 'FLOW');
cfg.process.calibrate.ACS91 = struct('compute_dissolved', true,...
                                  'interpolation_method', 'linear', ...
                                  'CDOM_source', 'SPCD',... 
                                  'FLOW_source', 'FLOW');
cfg.process.calibrate.BB3 = struct('compute_dissolved', true,...
                                   'TSG_source', 'TSG',...
                                   'di_method', 'constant');
cfg.process.calibrate.skip = {'FLOW', 'TSG', 'SPCD'};

%%% Write %%%
cfg.process.write = struct();
cfg.process.write.mode = 'One day one file'; % 'One file' or 'One day one file'
cfg.process.write.skip = {'FLOW'};

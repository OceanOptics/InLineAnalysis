% Tara Pacific Configuration file
% author: Nils
% created: Oct 17, 2018

cfg = struct('meta', struct(), 'instruments', struct(), 'process', struct());

%%%%%%%%%%%%%%
%% METADATA %%
%%%%%%%%%%%%%%
cfg.meta.investigators = 'Emmanuel Boss, Nils Haentjens, Guillaume Bourdin';
cfg.meta.affiliations = 'University of Maine';
cfg.meta.emails = 'emmanuel.boss@maine.edu';
cfg.meta.experiment = 'Tara Pacific';
cfg.meta.cruise = 'Tara Pacific';
cfg.meta.station = 'NA';
cfg.meta.documents = 'NA';
cfg.meta.calibration_files = 'NA';
cfg.meta.data_type = 'flow_thru';
cfg.meta.data_status = 'preliminary';
cfg.meta.measurement_depth = 1.5;


%%%%%%%%%%%%%%%%%
%% INSTRUMENTS %%
%%%%%%%%%%%%%%%%%

PATH_ROOT = '/home/haentjens/InLineAnalysis/Data/TaraPacific/';
% PATH_ROOT = '/Users/nils/Data/TaraPacific/';

%%% TSG + GPS %%%
cfg.instruments.TSG = struct();
cfg.instruments.TSG.model = 'TSG';
cfg.instruments.TSG.boat = 'Tara';
cfg.instruments.TSG.path = struct('raw',  [PATH_ROOT 'TSG/raw/'],...
                                  'wk',   [PATH_ROOT 'TSG/wk/'],...
                                  'ui',   [PATH_ROOT 'TSG/user_input/'],...
                                  'prod', [PATH_ROOT 'TSG/']);
cfg.instruments.TSG.view = struct('varname', 't');

%%% FTH (FlowControl) %%%
cfg.instruments.FTH = struct();
cfg.instruments.FTH.model = 'FTH';
cfg.instruments.FTH.logger = 'FlowControl';
cfg.instruments.FTH.path = struct('raw',  [PATH_ROOT 'FlowControl/raw_dot/'],...
                                  'wk',   [PATH_ROOT 'FlowControl/wk/'],...
                                  'ui',   [PATH_ROOT 'FlowControl/user_input/'],...
                                  'prod', [PATH_ROOT 'FlowControl/']);
cfg.instruments.FTH.view = struct('varname', 'swt');

%%% ACS 301 %%% (Aug 20 to ...)
% cfg.instruments.ACS = struct();
% cfg.instruments.ACS.model = 'ACS';
% cfg.instruments.ACS.sn = '301';
% cfg.instruments.ACS.logger = 'Compass_2.1rc_scheduled';
% cfg.instruments.ACS.lambda_reference = [400.9, 404.6, 408.4, 411.7, 415.5, 419.2, 423.7, 428.1, 432.2, 435.9, 440, 444, 449, 453.4, 457.6, 462, 466.2, 470.8, 475.7, 480.3, 485.2, 489.4, 493.7, 497.7, 502.1, 506.6, 511.4, 515.9, 520.6, 524.9, 529.2, 533.4, 537.5, 541.5, 545.7, 550, 554.3, 558.6, 562.8, 566.9, 570.7, 574.1, 578, 582.2, 586.1, 590.3, 594.5, 599.1, 603.8, 608.1, 612.7, 617.4, 621.7, 626, 630.5, 634.7, 639, 643.3, 647.8, 652.4, 656.8, 661.4, 665.9, 670.3, 674.4, 678.8, 682.8, 687.1, 691, 694.7, 698.4, 702.5, 706.2, 710, 713.7, 717.6, 721.4, 725.1, 728.8, 732.6, 736.2, 739.6];
% cfg.instruments.ACS.lambda_a = [400.9, 404.6, 408.4, 411.7, 415.5, 419.2, 423.7, 428.1, 432.2, 435.9, 440, 444, 449, 453.4, 457.6, 462, 466.2, 470.8, 475.7, 480.3, 485.2, 489.4, 493.7, 497.7, 502.1, 506.6, 511.4, 515.9, 520.6, 524.9, 529.2, 533.4, 537.5, 541.5, 545.7, 550, 554.3, 558.6, 562.8, 566.9, 570.7, 574.1, 578, 582.2, 586.1, 590.3, 594.5, 599.1, 603.8, 608.1, 612.7, 617.4, 621.7, 626, 630.5, 634.7, 639, 643.3, 647.8, 652.4, 656.8, 661.4, 665.9, 670.3, 674.4, 678.8, 682.8, 687.1, 691, 694.7, 698.4, 702.5, 706.2, 710, 713.7, 717.6, 721.4, 725.1, 728.8, 732.6, 736.2, 739.6];
% cfg.instruments.ACS.lambda_c = [400.9, 404.6, 408.4, 411.7, 415.5, 419.2, 423.7, 428.1, 432.2, 435.9, 440, 444, 449, 453.4, 457.6, 462, 466.2, 470.8, 475.7, 480.3, 485.2, 489.4, 493.7, 497.7, 502.1, 506.6, 511.4, 515.9, 520.6, 524.9, 529.2, 533.4, 537.5, 541.5, 545.7, 550, 554.3, 558.6, 562.8, 566.9, 570.7, 574.1, 578, 582.2, 586.1, 590.3, 594.5, 599.1, 603.8, 608.1, 612.7, 617.4, 621.7, 626, 630.5, 634.7, 639, 643.3, 647.8, 652.4, 656.8, 661.4, 665.9, 670.3, 674.4, 678.8, 682.8, 687.1, 691, 694.7, 698.4, 702.5, 706.2, 710, 713.7, 717.6, 721.4, 725.1, 728.8, 732.6, 736.2, 739.6];
% cfg.instruments.ACS.path = struct('raw',  '/Users/nils/Data/EXPORTS/InLine/raw/ACS301/'],...
%                                   'di',  '/Users/nils/Data/EXPORTS/InLine/raw/ACS301/DI/'],...
%                                   'wk',   '/Users/nils/Data/EXPORTS/InLine/wk/ACS301/'],...
%                                   'prod', '/Users/nils/Data/EXPORTS/InLine/prod/ACS301/');
% cfg.instruments.ACS.view = struct('varname', 'a', 'varcol', 40);

%%% ACS 298 %%% (Aug 11 to Aug 20)
% cfg.instruments.ACS298 = struct();
% cfg.instruments.ACS298.model = 'ACS';
% cfg.instruments.ACS298.sn = '298';
% cfg.instruments.ACS298.ila_prefix = 'ACS';
% cfg.instruments.ACS298.logger = 'Compass_2.1rc_scheduled';
% cfg.instruments.ACS298.lambda_reference = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
% cfg.instruments.ACS298.lambda_a = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
% cfg.instruments.ACS298.lambda_c = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
% cfg.instruments.ACS298.path = struct('raw',  '/Users/nils/Data/EXPORTS/InLine/raw/ACS298/'],...
%                                   'di',  '/Users/nils/Data/EXPORTS/InLine/raw/ACS298/DI/'],...
%                                   'wk',   '/Users/nils/Data/EXPORTS/InLine/wk/ACS298/'],...
%                                   'prod', '/Users/nils/Data/EXPORTS/InLine/prod/ACS298/');
% cfg.instruments.ACS298.view = struct('varname', 'a', 'varcol', 40);

%%% BB3 %%%
cfg.instruments.BB3 = struct();
cfg.instruments.BB3.model = 'BB';
cfg.instruments.BB3.sn = '1502';
cfg.instruments.BB3.ila_prefix = 'BB3';
cfg.instruments.BB3.logger = 'InlininoBB3';
% cfg.instruments.BB3.logger = 'TeraTermBB3';
cfg.instruments.BB3.lambda = [470,532,650];
cfg.instruments.BB3.theta = 124;
% Pre-Cruise Calibration from 2016/10/20 ECO BB3-1502
cfg.instruments.BB3.slope = [1.066e-05, 7.076e-06, 3.569e-06];
cfg.instruments.BB3.dark = [50, 44, 46];
cfg.instruments.BB3.path = struct('raw',  [PATH_ROOT 'BB3/raw_inlinino/'],...
                                  'wk',   [PATH_ROOT 'BB3/wk/'],...
                                  'ui',   [PATH_ROOT 'BB3/user_input/'],...
                                  'prod', [PATH_ROOT 'BB3/']);
cfg.instruments.BB3.view = struct('varname', 'beta', 'varcol', 2);

%%% WSCD %%%
cfg.instruments.WSCD = struct();
cfg.instruments.WSCD.model = 'CD';
cfg.instruments.WSCD.sn = '';
cfg.instruments.WSCD.ila_prefix = 'WSCD';
cfg.instruments.WSCD.logger = 'InlininoWSCD';
cfg.instruments.WSCD.slope = NaN;
cfg.instruments.WSCD.dark = NaN;
cfg.instruments.WSCD.path = struct('raw',  [PATH_ROOT 'CDOM/raw/'],...
                                   'wk',   [PATH_ROOT 'CDOM/wk/'],...
                                   'ui',   [PATH_ROOT 'CDOM/user_input/'],...
                                   'prod', [PATH_ROOT 'CDOM/']);
cfg.instruments.WSCD.view = struct('varname', 'fdom');

%%% PAR %%%
cfg.instruments.PAR = struct();
cfg.instruments.PAR.model = 'PAR';
cfg.instruments.PAR.sn = '';
cfg.instruments.PAR.logger = 'Inlinino';
cfg.instruments.PAR.scale = 6.451E-04; % Volts/(uE/m²sec)
cfg.instruments.PAR.path = struct('raw',  [PATH_ROOT 'PAR/raw/'],...
                                   'wk',   [PATH_ROOT 'PAR/wk/'],...
                                   'ui',   [PATH_ROOT 'PAR/user_input/'],...
                                   'prod', [PATH_ROOT 'PAR/']);
cfg.instruments.PAR.view = struct('varname', 'par');


%%%%%%%%%%%%%
%% PROCESS %%
%%%%%%%%%%%%%

%%% General parameters %%%
cfg.process.days2run = datenum(2016,05,28):datenum(2018,10,27);
cfg.process.instruments2run = {'FTH', 'TSG', 'ACS', 'BB3', 'WSCD'};
cfg.process.write = true;
cfg.process.force_import = false;
cfg.process.parallel = 36; % 0: disable parallel or Inf: as many thread available

%%% Synchronization %%%
cfg.process.sync = struct();
cfg.process.sync.delay = struct();
cfg.process.sync.delay.FTH = 0;
cfg.process.sync.delay.ACS = 66;
cfg.process.sync.delay.BB3 = 0;
cfg.process.sync.delay.WSCD = 5;
cfg.process.sync.skip = {'TSG', 'PAR'};

%%% QC Reference (Flow Control/FTH) %%%
cfg.process.qcref = struct();
cfg.process.qcref.reference = 'FTH';
cfg.process.qcref.view = 'BB3';
cfg.process.qcref.mode = 'load'; % load or ui

%%% Split total/filtered %%%
cfg.process.split = struct();
cfg.process.split.reference = 'FTH';
cfg.process.split.buffer = struct();
cfg.process.split.buffer.ACS = [180, 30];
cfg.process.split.buffer.BB3 = [420, 220];
cfg.process.split.buffer.WSCD = [310, 20];
cfg.process.split.skip = {'FTH', 'TSG', 'PAR'};

%%% Binning %%%
cfg.process.bin = struct('bin_size', struct());
% cfg.process.bin.method = '4flag'; % Method to use to flag automatically
cfg.process.bin.method = 'SB_IN_PRCTL'; % Faster method outputting only requirements for SeaBASS (not compatible with automatic flagging)
% cfg.process.bin.method = 'SB_ALL'; % Faster method outputting only requirements for SeaBASS (not compatible with automatic flagging)
cfg.process.bin.prctile_detection = [2.5, 97.5];
cfg.process.bin.prctile_average = [5, 75];
cfg.process.bin.bin_size.ACS = 1;
cfg.process.bin.bin_size.BB3 = 1;
cfg.process.bin.bin_size.WSCD = 1;
cfg.process.bin.bin_size.TSG = 1;
cfg.process.bin.bin_size.PAR = 1; 
cfg.process.bin.skip = {'FTH'};

%%% Automatically flagging %%%
cfg.process.flag = struct();
cfg.process.flag.skip = {'FTH', 'TSG', 'BB3', 'WSCD', 'ACS', 'PAR'};
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
cfg.process.flag.BB3.abs_uncertainty = 2;
cfg.process.flag.BB3.rel_uncertainty = 0.05;
cfg.process.flag.BB3.min_flag_n =  2;
cfg.process.flag.BB3.primary_varname = 'beta';
% Flag parameters specific to the WSCD
cfg.process.flag.WSCD = struct();
cfg.process.flag.WSCD.abs_uncertainty = 0.3;
cfg.process.flag.WSCD.rel_uncertainty = 0.0125;
cfg.process.flag.WSCD.primary_varname = 'fdom';
  
%%% Manually QC %%%
cfg.process.qc = struct();
cfg.process.qc.mode = 'load';
cfg.process.qc.global = struct();
cfg.process.qc.global.active = false;
cfg.process.qc.global.view = 'ACS';
cfg.process.qc.global.apply = {'ACS', 'BB3', 'WSCD'};
cfg.process.qc.specific = struct();
cfg.process.qc.specific.active = true;
cfg.process.qc.specific.run = {'TSG', 'ACS', 'BB3', 'WSCD'};
  
%%% Calibrate %%%
cfg.process.calibrate = struct();
cfg.process.calibrate.ACS = struct('compute_dissolved', false,...
                                  'interpolation_method', 'linear', ...
                                  'CDOM_source', 'WSCD',...
                                  'FTH_source', 'FTH');
cfg.process.calibrate.BB3 = struct('compute_dissolved', false,...
                                   'TSG_source', 'TSG',...
                                   'di_method', 'constant');
cfg.process.calibrate.skip = {'FTH', 'TSG', 'WSCD'};

%%% Write %%%
cfg.process.write = struct();
cfg.process.write.mode = 'One day one file'; % 'One file' or 'One day one file'
cfg.process.write.skip = {'FTH'};

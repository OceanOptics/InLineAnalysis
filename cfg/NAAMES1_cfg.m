% NAAMES 1 Configuration file
% author: Nils
% created: Nov 1, 2018

cfg = struct('meta', struct(), 'instruments', struct(), 'process', struct());

%%%%%%%%%%%%%%
%% METADATA %%
%%%%%%%%%%%%%%
cfg.meta.investigators = 'Emmanuel Boss, Nils Haentjens';
cfg.meta.affiliations = 'University of Maine';
cfg.meta.emails = 'emmanuel.boss@maine.edu';
cfg.meta.experiment = 'NAAMES01';
cfg.meta.cruise = 'AT32';
cfg.meta.station = 'NA';
cfg.meta.documents = 'NA';
cfg.meta.calibration_files = 'NA';
cfg.meta.data_type = 'flow_thru';
cfg.meta.data_status = 'preliminary';
cfg.meta.measurement_depth = 5;


%%%%%%%%%%%%%%%%%
%% INSTRUMENTS %%
%%%%%%%%%%%%%%%%%
if exist('PATH_ROOT', 'var'); fprintf('Using global PATH_ROOT\n');
else; PATH_ROOT = '/Users/nils/Data/NAAMES/NAAMES1/'; end

%%% TSG + GPS %%%
cfg.instruments.TSG = struct();
cfg.instruments.TSG.model = 'TSG';
cfg.instruments.TSG.boat = 'Atlantis';
cfg.instruments.TSG.path = struct('raw',  [PATH_ROOT 'Underway/raw/'],...
                                  'wk',   [PATH_ROOT 'Underway/mat/'],...
                                  'prod', [PATH_ROOT 'Underway/sb/'],...
                                  'ui', [PATH_ROOT 'Underway/ui/']);
cfg.instruments.TSG.view = struct('varname', 't');

%%% FTH (FlowControl) %%%
cfg.instruments.FTH = struct();
cfg.instruments.FTH.model = 'FTH';
cfg.instruments.FTH.logger = 'FlowControl_old';
cfg.instruments.FTH.path = struct('raw',  [PATH_ROOT 'FlowControl/raw/'],...
                                  'wk',   [PATH_ROOT 'FlowControl/mat/'],...
                                  'prod', [PATH_ROOT 'FlowControl/sb/'],...
                                  'ui', [PATH_ROOT 'FlowControl/ui/']);
cfg.instruments.FTH.view = struct('varname', 'swt');

%%% ACS 091 %%%
cfg.instruments.ACS = struct();
cfg.instruments.ACS.model = 'ACS';
cfg.instruments.ACS.sn = '091';
cfg.instruments.ACS.logger = 'Compass_2.1rc_scheduled_bin';
cfg.instruments.ACS.device_file = [PATH_ROOT 'ACS/acs091.dev'];
% cfg.instruments.ACS.lambda_reference = [400.9, 404.6, 408.4, 411.7, 415.5, 419.2, 423.7, 428.1, 432.2, 435.9, 440, 444, 449, 453.4, 457.6, 462, 466.2, 470.8, 475.7, 480.3, 485.2, 489.4, 493.7, 497.7, 502.1, 506.6, 511.4, 515.9, 520.6, 524.9, 529.2, 533.4, 537.5, 541.5, 545.7, 550, 554.3, 558.6, 562.8, 566.9, 570.7, 574.1, 578, 582.2, 586.1, 590.3, 594.5, 599.1, 603.8, 608.1, 612.7, 617.4, 621.7, 626, 630.5, 634.7, 639, 643.3, 647.8, 652.4, 656.8, 661.4, 665.9, 670.3, 674.4, 678.8, 682.8, 687.1, 691, 694.7, 698.4, 702.5, 706.2, 710, 713.7, 717.6, 721.4, 725.1, 728.8, 732.6, 736.2, 739.6];
% cfg.instruments.ACS.lambda_a = [400.9, 404.6, 408.4, 411.7, 415.5, 419.2, 423.7, 428.1, 432.2, 435.9, 440, 444, 449, 453.4, 457.6, 462, 466.2, 470.8, 475.7, 480.3, 485.2, 489.4, 493.7, 497.7, 502.1, 506.6, 511.4, 515.9, 520.6, 524.9, 529.2, 533.4, 537.5, 541.5, 545.7, 550, 554.3, 558.6, 562.8, 566.9, 570.7, 574.1, 578, 582.2, 586.1, 590.3, 594.5, 599.1, 603.8, 608.1, 612.7, 617.4, 621.7, 626, 630.5, 634.7, 639, 643.3, 647.8, 652.4, 656.8, 661.4, 665.9, 670.3, 674.4, 678.8, 682.8, 687.1, 691, 694.7, 698.4, 702.5, 706.2, 710, 713.7, 717.6, 721.4, 725.1, 728.8, 732.6, 736.2, 739.6];
% cfg.instruments.ACS.lambda_c = [400.9, 404.6, 408.4, 411.7, 415.5, 419.2, 423.7, 428.1, 432.2, 435.9, 440, 444, 449, 453.4, 457.6, 462, 466.2, 470.8, 475.7, 480.3, 485.2, 489.4, 493.7, 497.7, 502.1, 506.6, 511.4, 515.9, 520.6, 524.9, 529.2, 533.4, 537.5, 541.5, 545.7, 550, 554.3, 558.6, 562.8, 566.9, 570.7, 574.1, 578, 582.2, 586.1, 590.3, 594.5, 599.1, 603.8, 608.1, 612.7, 617.4, 621.7, 626, 630.5, 634.7, 639, 643.3, 647.8, 652.4, 656.8, 661.4, 665.9, 670.3, 674.4, 678.8, 682.8, 687.1, 691, 694.7, 698.4, 702.5, 706.2, 710, 713.7, 717.6, 721.4, 725.1, 728.8, 732.6, 736.2, 739.6];
cfg.instruments.ACS.path = struct('raw',  [PATH_ROOT 'ACS/raw/'],...
                                  'wk',   [PATH_ROOT 'ACS/wk/'],...
                                  'prod', [PATH_ROOT 'ACS/sb/'],...
                                  'ui', [PATH_ROOT 'ACS/ui/']);
cfg.instruments.ACS.view = struct('varname', 'a', 'varcol', 40);

%%% BB3 %%%
cfg.instruments.BB3 = struct();
cfg.instruments.BB3.model = 'BB';
cfg.instruments.BB3.sn = '349';
cfg.instruments.BB3.ila_prefix = 'BB3';
cfg.instruments.BB3.logger = 'DH4PreProc';
cfg.instruments.BB3.lambda = [470,532,660];
cfg.instruments.BB3.theta = 124;
cfg.instruments.BB3.slope = [1,1,1];
cfg.instruments.BB3.dark = [0,0,0];
cfg.instruments.BB3.path = struct('raw',  [PATH_ROOT 'BB3_WSCD/raw/'],...
                                  'wk',   [PATH_ROOT 'BB3_WSCD/mat/'],...
                                  'prod', [PATH_ROOT 'BB3_WSCD/sb/'],...
                                  'ui', [PATH_ROOT 'BB3_WSCD/ui/']);
cfg.instruments.BB3.view = struct('varname', 'beta', 'varcol', 2);

%%% WSCD %%%
cfg.instruments.WSCD = struct();
cfg.instruments.WSCD.model = 'CD';
cfg.instruments.WSCD.sn = '1299';
cfg.instruments.WSCD.ila_prefix = 'WSCD';
cfg.instruments.WSCD.logger = 'DH4PreProc';
cfg.instruments.WSCD.path = struct('raw',  [PATH_ROOT 'BB3_WSCD/raw/'],...
                                  'wk',   [PATH_ROOT 'BB3_WSCD/mat/'],...
                                  'prod', [PATH_ROOT 'BB3_WSCD/sb/'],...
                                  'ui', [PATH_ROOT 'BB3_WSCD/ui/']);
cfg.instruments.WSCD.view = struct('varname', 'fdom');

%%%%%%%%%%%%%
%% PROCESS %%
%%%%%%%%%%%%%

%%% General parameters %%%
cfg.process.days2run = datenum('06-Nov-2015'):datenum('30-Nov-2015');
cfg.process.instruments2run = {'FTH', 'TSG', 'ACS', 'BB3', 'WSCD'};
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
cfg.process.sync.delay.ACS = 80;
cfg.process.sync.delay.BB3 = -18640;
cfg.process.sync.delay.WSCD = -18600;
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
cfg.process.split.buffer.ACS = [180, 45];
cfg.process.split.buffer.BB3 = [480, 300];
cfg.process.split.buffer.WSCD = [330, 120];
cfg.process.split.skip = {'FTH', 'TSG'};

%%% Automatically flagging %%%
cfg.process.flag = struct();
cfg.process.flag.skip = {'FTH', 'TSG', 'BB3', 'WSCD', 'ACS'};

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
cfg.process.bin.bin_size.WSCD = 1;
cfg.process.bin.bin_size.TSG = 1; % TSG does not need to be binned in most cases
cfg.process.bin.skip = {'FTH'};

%%% Automatically flagging %%%
cfg.process.flag = struct();
cfg.process.flag.skip = {'FTH', 'TSG', 'BB3', 'WSCD', 'ACS'};
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
cfg.process.qc.specific.run = {'TSG', 'ACS', 'BB3', 'WSCD'};
  
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

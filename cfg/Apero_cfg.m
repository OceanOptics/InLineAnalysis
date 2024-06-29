% Tara Europa Configuration file
% author: Guilaume Bourdin & Emmanuel Boss
% created: July 29, 2023

cfg = struct('meta', struct(), 'instruments', struct(), 'process', struct());

%%%%%%%%%%%%%%
%% METADATA %%
%%%%%%%%%%%%%%
cfg.meta.investigators = 'Emmanuel_Boss';
cfg.meta.affiliations = 'University_of_Maine';
cfg.meta.emails = 'emmanuel.boss@maine.edu';
cfg.meta.experiment = 'Apero';
cfg.meta.cruise = 'THALASSA';
cfg.meta.station = 'NA';
cfg.meta.documents = 'NA';
cfg.meta.calibration_files = 'NA';
cfg.meta.data_type = 'flow_thru';
cfg.meta.data_status = 'preliminary';
cfg.meta.measurement_depth = 2;


%%%%%%%%%%%%%%%%%
%% INSTRUMENTS %%
%%%%%%%%%%%%%%%%%

PATH_ROOT = '/Volumes/Samsung_T5/Data/Apero';
% PATH_ROOT = '/Volumes/Data/Apero';

%%% TSG + GPS %%%
model = 'SBE45';
SN = '9999';
cfg.instruments.([model SN]) = struct();
cfg.instruments.([model SN]).model = model;
cfg.instruments.([model SN]).boat = 'Tara';
cfg.instruments.([model SN]).logger = 'Inlinino_base'; % TeraTerm Matlab Inlinino
cfg.instruments.([model SN]).sn = SN;
cfg.instruments.([model SN]).path = struct('raw',  fullfile(PATH_ROOT, 'raw', [model SN]),...
                                  'wk',   fullfile(PATH_ROOT, 'wk', [model SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', [model SN]));
cfg.instruments.([model SN]).view = struct('varname', 't1');
cfg.instruments.([model SN]).temperature_variable = 't1';

%%% FLOW (FlowControl) %%%
model = 'ADU100';
SN = 'B02694';
cfg.instruments.FLOW = struct();
cfg.instruments.FLOW.model = model;
cfg.instruments.FLOW.logger = 'Inlinino_base';
cfg.instruments.FLOW.sn = SN;
cfg.instruments.FLOW.LoadPrevious = true;
cfg.instruments.FLOW.analog1 = '';
cfg.instruments.FLOW.analog2 = 'SUVF6253';
cfg.instruments.FLOW.path = struct('raw',  fullfile(PATH_ROOT, 'raw', [model SN]),...
                                  'wk',   fullfile(PATH_ROOT, 'wk', [model SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', [model SN]));
cfg.instruments.FLOW.view = struct('varname', 'swt','spd_variable','spd1'); % spd1 spd2

%%% ACS348 %%%
SN = '348';
cfg.instruments.(['ACS' SN]) = struct();
cfg.instruments.(['ACS' SN]).di = struct();
cfg.instruments.(['ACS' SN]).di.prefix = ['DIW_ACS' SN];
cfg.instruments.(['ACS' SN]).di.postfix = '';
cfg.instruments.(['ACS' SN]).model = 'ACS';
cfg.instruments.(['ACS' SN]).sn = SN;
cfg.instruments.(['ACS' SN]).logger = 'InlininoACScsv';
cfg.instruments.(['ACS' SN]).device_file = fullfile(PATH_ROOT, 'DeviceFiles', 'ACS-348_20230504.dev');
cfg.instruments.(['ACS' SN]).path = struct('raw',  fullfile(PATH_ROOT, 'raw', ['ACS' SN]),...
                                    'di',  fullfile(PATH_ROOT, 'raw', ['ACS' SN], 'DI'),...
                                    'wk',   fullfile(PATH_ROOT, 'wk', ['ACS' SN]),...
                                    'prod', fullfile(PATH_ROOT, 'prod'),...
                                    'ui', fullfile(PATH_ROOT, 'ui', ['ACS' SN]));
cfg.instruments.(['ACS' SN]).view = struct('varname', 'a', 'varcol', 40);


%%% SUVF %%%
SN = '6253';
cfg.instruments.(['SUVF' SN]) = struct();
cfg.instruments.(['SUVF' SN]).model = 'CD';
cfg.instruments.(['SUVF' SN]).sn = SN;
cfg.instruments.(['SUVF' SN]).slope = 1;
cfg.instruments.(['SUVF' SN]).dark = 0;
cfg.instruments.(['SUVF' SN]).logger = 'Inlinino_base';
cfg.instruments.(['SUVF' SN]).analog_channel = 'Analog2';
if strcmp(cfg.instruments.(['SUVF' SN]).logger, 'InlininoADU100') || contains(cfg.instruments.FLOW.model, 'ADU')
  cfg.instruments.(['SUVF' SN]).ila_prefix = [cfg.instruments.FLOW.model cfg.instruments.FLOW.sn];
  [~, raw_folder_name] = fileparts(cfg.instruments.FLOW.path.raw);
else
  cfg.instruments.(['SUVF' SN]).ila_prefix = ['SUVF' SN];
  raw_folder_name = ['SUVF' SN];
end
cfg.instruments.(['SUVF' SN]).di = struct();
cfg.instruments.(['SUVF' SN]).di.prefix = ['DIW_' cfg.instruments.(['SUVF' SN]).ila_prefix SN];
cfg.instruments.(['SUVF' SN]).di.postfix = '';
cfg.instruments.(['SUVF' SN]).path = struct('raw',  fullfile(PATH_ROOT, 'raw', raw_folder_name),...
                                  'di',  fullfile(PATH_ROOT, 'raw', ['SUVF' SN], 'DI'),...
                                  'wk',   fullfile(PATH_ROOT, 'wk', ['SUVF' SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', ['SUVF' SN]));
cfg.instruments.(['SUVF' SN]).view = struct('varname', 'C2', 'varcol', 1);

%%% BB3 %%%
SN = '1502';
cfg.instruments.(['BB3' SN]) = struct();
cfg.instruments.(['BB3' SN]).di = struct();
cfg.instruments.(['BB3' SN]).di.prefix = ['DIW_BB3' SN];
cfg.instruments.(['BB3' SN]).di.postfix = '';
cfg.instruments.(['BB3' SN]).model = 'BB';
cfg.instruments.(['BB3' SN]).sn = SN;
cfg.instruments.(['BB3' SN]).ila_prefix = ['BB3' SN];
cfg.instruments.(['BB3' SN]).logger = 'Inlinino_base'; % InlininoBB3SN Inlinino_base
cfg.instruments.(['BB3' SN]).lambda = [468,528,655];
cfg.instruments.(['BB3' SN]).theta = 120;
cfg.instruments.(['BB3' SN]).slope = [8.407E-06,4.624E-06,4.090E-06]; % Emmanuel cal 2021
cfg.instruments.(['BB3' SN]).dark = [50, 44, 46];
cfg.instruments.(['BB3' SN]).path = struct('raw',  fullfile(PATH_ROOT, 'raw', ['BB3' SN]),...
                                  'di',  fullfile(PATH_ROOT, 'raw', ['BB3' SN]),...
                                  'wk',   fullfile(PATH_ROOT, 'wk', ['BB3' SN]),...
                                  'prod', fullfile(PATH_ROOT, 'prod'),...
                                  'ui', fullfile(PATH_ROOT, 'ui', ['BB3' SN]));
cfg.instruments.(['BB3' SN]).view = struct('varname', 'beta', 'varcol', 2);



%% %%%%%%% %%
%  PROCESS  %
%%%%%%%%%%%%%

%%% General parameters %%%
cfg.process.days2run = datenum(2023,6,3):datenum(2023,7,11);
% cfg.process.instruments2run = {'FLOW', 'NMEA', 'ACS57', 'BB31502', 'WSCD859', ...
%   'SBE4536073', 'SUVF6244', 'LISST1183', 'HyperBB8005'};
cfg.process.instruments2run = fieldnames(cfg.instruments);
cfg.process.write = true;
cfg.process.force_import = false;
cfg.process.parallel = Inf; % 0: disable parallel or Inf: as many thread available
cfg.process.di = struct();
cfg.process.di.skip = cfg.process.instruments2run(contains(cfg.process.instruments2run, ...
  {'FLOW','TSG','SBE45','SBE3845','NMEA','PAR'}));
cfg.process.di.qc = struct('mode', 'ui',... % ui or load
                           'qc_once_for_all', false,... % true = QC all variables | false = QC variables separately);
                           'remove_old', false); % remove old selection of the same period
cfg.process.di.bin = struct('bin_size', 30);

%%% Synchronization %%%
cfg.process.sync = struct();
cfg.process.sync.delay = struct();
cfg.process.sync.skip = cfg.process.instruments2run(contains(cfg.process.instruments2run, ...
  {'FLOW','TSG','SBE45','SBE3845','ALFA','NMEA', 'PAR'}));
% Set default sync delay.
% To customize sync delay, uncomment section below
for i = 1:size(cfg.process.instruments2run)
  cfg.process.sync.delay.(cfg.process.instruments2run{i}) = 0;
end
% % Manually customize sync delay
% cfg.process.sync.delay.FLOW = 0;
% cfg.process.sync.delay.ACS412 = 0;
% cfg.process.sync.delay.BB3349 = 0;
% cfg.process.sync.delay.WS3S1081 = 0;
% cfg.process.sync.delay.NMEA = 0;
% cfg.process.sync.delay.HyperBB8002 = 0;
% cfg.process.sync.delay.LISST1183 = 0;
% cfg.process.sync.delay.SUVF6244 = 0;

%%% QC Reference (Flow Control/FLOW) %%%
cfg.process.qcref = struct();
cfg.process.qcref.reference = 'FLOW';
cfg.process.qcref.view = cfg.process.instruments2run(find(contains(cfg.process.instruments2run, ...
  {'ACS', 'AC9'}),1, 'first'));
cfg.process.qcref.mode = 'ui'; % load or ui
cfg.process.qcref.remove_old = false; % remove old selection of the same period
cfg.process.qcref.MinFiltPeriod = 50; % filter even period in minute
cfg.process.qcref.szFilt = 10; % filter even length in minute

%%% Split total/filtered %%%
cfg.process.split = struct();
cfg.process.split.reference = 'FLOW';
cfg.process.split.buffer = struct();
cfg.process.split.skip = cfg.process.instruments2run(contains(cfg.process.instruments2run, ...
  {'FLOW','TSG','SBE45','SBE3845','ALFA','NMEA', 'PAR'}));
% Set buffer length depending on instrument type (default).
% To customize buffer length, uncomment section below
for i = 1:size(cfg.process.instruments2run)
  if any(contains(cfg.process.instruments2run{i}, {'ACS', 'AC9'}))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [180, 60]; % [180, 60] for AC meters
  elseif any(contains(cfg.process.instruments2run{i}, {'BB', 'BB'}) & ~contains(cfg.process.instruments2run{i}, {'HyperBB','HBB','hbb'}))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [200, 100]; % [420, 220] for ECO-BB
  elseif any(contains(cfg.process.instruments2run{i}, 'WSCD'))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [200, 100]; % [540, 100] for ECO-fluo
  elseif any(contains(cfg.process.instruments2run{i}, 'WS3S'))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [200, 150]; % [420, 220] for ECO-fluo
  elseif any(contains(cfg.process.instruments2run{i}, 'SUVF'))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [200, 100]; % [240, 100] for Seapoint fluo
  elseif any(contains(cfg.process.instruments2run{i}, {'HyperBB', 'HBB', 'hbb'}))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [200, 100]; % [240, 140] for HyperBB
  elseif any(contains(cfg.process.instruments2run{i}, 'LISST') & ~contains(cfg.process.instruments2run{i}, {'TAU','tau','Tau','200X','200x'}))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [180, 60]; % [540, 360] for LISST
  elseif any(contains(cfg.process.instruments2run{i}, {'LISST200X', 'LISST200x'}))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [180, 60]; % [540, 360] for LISST
  elseif any(contains(cfg.process.instruments2run{i}, {'LISSTTAU','LISSTTau','LISST-TAU'}))
    cfg.process.split.buffer.(cfg.process.instruments2run{i}) = [180, 60]; % [180, 60] for LISST-Tau
  end
end
% % Manually customize buffer length
% cfg.process.split.buffer.ACS57 = [180, 60];
% cfg.process.split.buffer.ACS348 = [180, 60];
% cfg.process.split.buffer.LISSTTau1002G = [180, 60];
% cfg.process.split.buffer.BB31502 = [420, 220];
% cfg.process.split.buffer.WSCD859 = [540, 100];
% cfg.process.split.buffer.SUVF6244 = [240, 100]; % [660, 100]
% cfg.process.split.buffer.HyperBB8005 = [240, 140]; % [540, 340]
% cfg.process.split.buffer.LISST1183 = [540, 360];

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
for i = 1:size(cfg.process.instruments2run)
  if contains(cfg.process.instruments2run{i}, 'FLOW')
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for FLOW
  elseif any(contains(cfg.process.instruments2run{i}, {'GPS', 'NMEA'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for NMEA
  elseif any(contains(cfg.process.instruments2run{i}, {'TSG', 'SBE38', 'SBE45'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for TSG
  elseif any(contains(cfg.process.instruments2run{i}, {'ACS', 'AC9'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for AC meters
  elseif any(contains(cfg.process.instruments2run{i}, {'BB', 'BB'}) & ~contains(cfg.process.instruments2run{i}, {'HyperBB','HBB','hbb'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for ECO-BB
  elseif any(contains(cfg.process.instruments2run{i}, 'WSCD'))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for ECO-fluo
  elseif any(contains(cfg.process.instruments2run{i}, 'WS3S'))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for ECO-fluo
  elseif any(contains(cfg.process.instruments2run{i}, 'SUVF'))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for Seapoint fluo
  elseif any(contains(cfg.process.instruments2run{i}, {'HyperBB', 'HBB', 'hbb'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 5; % 5 min for HyperBB
  elseif any(contains(cfg.process.instruments2run{i}, 'LISST') & ~contains(cfg.process.instruments2run{i}, {'TAU','tau','Tau','200X','200x'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 2; % 2 min for LISST100X
  elseif any(contains(cfg.process.instruments2run{i}, 'LISST200X'))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for LISST200X
  elseif any(contains(cfg.process.instruments2run{i}, {'LISSTTAU','LISSTTau','LISST-TAU'}))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 1; % 1 min for LISST-Tau
  elseif any(contains(cfg.process.instruments2run{i}(1:3), 'ALFA'))
    cfg.process.bin.bin_size.(cfg.process.instruments2run{i}) = 10; % 10 min for ALFA
  end
end
% % Manually customize bin sizes
% cfg.process.bin.bin_size.FLOW = 1;
% cfg.process.bin.bin_size.ACS57 = 1;
% cfg.process.bin.bin_size.ACS348 = 1;
% cfg.process.bin.bin_size.LISSTTau1002G = 1;
% cfg.process.bin.bin_size.BB31502 = 1;
% cfg.process.bin.bin_size.WSCD859 = 1;
% cfg.process.bin.bin_size.SUVF6244 = 1;
% cfg.process.bin.bin_size.SBE38450091 = 1;
% cfg.process.bin.bin_size.NMEA = 1;
% cfg.process.bin.bin_size.HyperBB8005 = 5;
% cfg.process.bin.bin_size.LISST1183 = 10;
% cfg.process.bin.bin_size.SUVF6244 = 1;
% cfg.process.bin.bin_size.ALFA = 10;
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
cfg.process.qc.global = struct();
cfg.process.qc.global.active = false;
cfg.process.qc.global.view = cfg.process.qcref.view;
cfg.process.qc.global.apply = cfg.process.instruments2run(~contains(cfg.process.instruments2run, ...
  {'FLOW','NMEA', 'PAR'}));
cfg.process.qc.specific = struct();
cfg.process.qc.specific.active = true;
cfg.process.qc.specific.run = cfg.process.qcref.view;

%%% Calibrate %%%
cfg.process.calibrate = struct();
cfg.process.calibrate.skip = cfg.process.instruments2run(contains(lower(cfg.process.instruments2run), ...
  {'flow','tsg','sbe45','sbe3845','nmea','alfa'}));
cfg.process.min_nb_pts_per_cluster = 100;
cfg.process.time_weight_for_cluster = 0.9;
% look for TSG and SUVF and automatically turn off CDOM interpolation if not available
cfg.process.TSG_source = '';
cfg.process.CDOM_source = '';
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

% Set calibrate options depending on instrument type (default).
for i = 1:size(cfg.process.instruments2run)
  % AC meter options
  if any(contains(lower(cfg.process.instruments2run{i}), {'acs', 'ac9'}))
    cfg.process.calibrate.(cfg.process.instruments2run{i}) = struct('compute_dissolved', false, ...
                                      'TSG_source', cfg.process.TSG_source, ...
                                      'interpolation_method', 'linear', ... % choose one: linear CDOM
                                      'CDOM_source', cfg.process.CDOM_source, ... %
                                      'FLOW_source', 'FLOW', ...
                                      'di_method', 'best_di', ... % best_di normal
                                      'scattering_correction', 'ZaneveldRottgers_blended', ... % Zaneveld1994_proportional Rottgers2013_semiempirical ZaneveldRottgers_blended
                                      'compute_ad_aphi', false); % VERY SLOW: compute ad and aphi from Zheng and Stramski 2013
  % ECO-BB options
  elseif any(contains(lower(cfg.process.instruments2run{i}), 'bb') & ~contains(lower(cfg.process.instruments2run{i}), {'hyperbb', 'hbb'}))
    cfg.process.calibrate.(cfg.process.instruments2run{i}) = struct('compute_dissolved', true, ...
                                      'TSG_source', cfg.process.TSG_source, ...
                                      'FLOW_source', 'FLOW', ...
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

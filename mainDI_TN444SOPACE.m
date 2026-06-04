% Main InLine Analysis Script
% author: Guillaume Bourdin
% created: May 05, 2021
% clear
% close all
% if feature('IsDebugMode'); dbquit all; end
% 
% % cd('/Users/emmanuel.boss/Desktop/InLine analysis/InLineAnalysis-master/')
% cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master/')
% 
% % Load InLineAnalysis and the configuration
% ila = InLineAnalysis('cfg/TN444SOPACE_cfg.m');

%% Quick cfg update
ila.cfg.days2run = datetime(2025,5,5):datetime(2025,5,12);

%% 'FLOW', 'SBE4536073', 'ACS111', 'HyperBB8005', 'WS3S1081P','SUVF6244','LISST100X1183'
ila.cfg.instruments2run = {'FLOW','ACS111'};
ila.cfg.qcref.view = 'ACS111';
ila.cfg.parallel = 12; % Inf
ila.cfg.calibrate.(ila.cfg.qcref.view).compute_dissolved = true;

%% 1. Read DI 
% To get ag and cg from ACS can run DI day by day
% To get betag from BB3 need to run the entire dataset if select di method constant
ila.cfg.force_import = false;
ila.ReadRawDI();
ila.CheckDataStatus();

%% Load processed data from mat files: 'data' = Raw | 'bin' = Bin | 'qc' = QCed | 'prod' = product
% ila.Read('raw');
ila.Read('bin');
ila.Read('qc');
ila.Read('prod');

%% 1.1. SpectralQC Plot
% check raw spectrums AC, BB, or LISST sensors
ila.SpectralQC('AC',{'raw'});

%% 2. Pass2QC
ila.Pass2QC('dissolved') % Copy data to next level
ila.CheckDataStatus();

%% 2. Automatic QC of raw data for step in ACS spectrum, BB saturated and obvious bad PAR values
% fudge factor for auto QC ACS.
% Varies between ACS: 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.AutoQC_tolerance.dissolved.a = 3;
ila.cfg.qc.AutoQC_tolerance.dissolved.c = 6;
% fudge factor for auto QC BB.
% 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.AutoQC_tolerance.dissolved.bb = 3;
ila.AutoQC('qc');
ila.CheckDataStatus();

%% 2.1. SpectralQC Plot
% check raw spectrums AC, BB, or LISST sensors
ila.SpectralQC('AC',{'qc'});

%% 2.2. Write QC | write only 'part' or 'diw' or 'all'
ila.Write('raw', 'diw')
ila.CheckDataStatus();

%% 3. QC DI
ila.cfg.di.qc.mode = 'ui';
ila.cfg.di.qc.remove_old = false;  % remove old selection of this period
ila.cfg.di.qc.qc_once_for_all = false;  % true = QC all variables | false = QC variables separately)
ila.QCDI();

%% 3.1. Run QC directly on spectra at any level
% ila.SpectralQC inputs:
% 1) 'AC', 'BB', or 'LISST'
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.SpectralQC('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB, BB3, or LISST:  ila.SpectralQC('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'g','ag'})
ila.SpectralQC('AC',{'qc'}, false, {'diw','c'});

%% 3.2. SpectralQC Plot
% check QCed spectrums AC, BB, or LISST sensors
ila.SpectralQC('AC',{'qc'});

%% 4. Write QC | write only 'part' or 'diw' or 'all'
ila.Write('qc', 'diw')
ila.CheckDataStatus();

%% 5. BIN DI
ila.BinDI();

%% 5.1. SpectralQC Plot
% check binned spectrums AC, BB, or LISST sensors
ila.SpectralQC('AC',{'bin'}); % AC or BB

%% 5.2. Run QC directly on spectra at any level
% ila.SpectralQC inputs:
% 1) 'AC', 'BB', or 'LISST'
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.SpectralQC('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB, BB3, or LISST:  ila.SpectralQC('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'g','ag'})
ila.SpectralQC('AC',{'bin'}, false, {'diw','c'});

%% 6. Write bin DI | write only 'part' or 'diw' or 'all'
ila.Write('bin', 'diw')
ila.CheckDataStatus();

%% 7. Calibrate
% update filter event calcualtion method if needed: exponential_fit 25percentil
ila.cfg.calibrate.(ila.cfg.qcref.view).filt_method = 'exponential_fit'; 
% update filter interpolation method if needed: CDOM linear
ila.cfg.calibrate.(ila.cfg.qcref.view).interpolation_method = 'CDOM'; % linear CDOM
% update scattering correction method if needed: Rottgers2013_semiempirical Zaneveld1994_proportional Semiempirical_blended1 Semiempirical_blended2 Semiempirical_blended3 
ila.cfg.calibrate.(ila.cfg.qcref.view).scattering_correction = 'Semiempirical_blended2';

ila.cfg.calibrate.(ila.cfg.qcref.view).compute_dissolved = true;
ila.cfg.calibrate.(ila.cfg.qcref.view).di_method = 'normal'; % best_di normal
ila.Calibrate()
ila.CheckDataStatus();

%% 7.1. Normal and DI prod QC plots
save_figures = false;

%%% AC or BB 3D plots %%%
ila.SpectralQC('AC', {'prod'}, save_figures); % AC or BB

%%% ACS BB3 TSG PAR WSCD final product visualisation %%%
ila.visProd_timeseries()

%% 8. Run QC directly on spectra at any level
% ila.SpectralQC inputs:
% 1) 'AC', 'BB', or 'LISST'
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.SpectralQC('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB, BB3, or LISST:  ila.SpectralQC('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'g','ag'})
ila.SpectralQC('AC',{'prod'}, false, {'g','cg'});

%% 8.1. Load previous qc pick selection at prod level
ila.cfg.qc.mode='load';  % load or ui
ila.cfg.qc.specific.run = {ila.cfg.qcref.view}; % 'FLOW','ACS57','TSG', 'BB31502', 'WSCD859','PAR'
ila.QC();

%% 9. Save products | write only 'part' or 'diw' or 'all'
ila.Write('prod', 'diw')

% % Notify with a song that the job is done
% notif_sound = load('gong'); sound(notif_sound.y, notif_sound.Fs); % handel
% return

%% re-write last version of 'qc' and 'bin' | write only 'part' or 'diw' or 'all'
ila.Write('raw', 'diw')
ila.Write('bin', 'diw')
ila.Write('qc', 'diw')











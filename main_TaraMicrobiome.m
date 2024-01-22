% Main Particulate InLine Analysis Script
% author: Guillaume Bourdin
% created: Jan 05, 2021
clear
close all
closeVar
cd('/PATH/TO/InLineAnalysis-masterDIRECTORY/')

% Load InLineAnalysis and the configuration
ila = InLineAnalysis('cfg/TaraMicrobiome_cfg.m');

% Quick cfg update
%% GPSSC701
% ila.cfg.days2run = datenum(2021,12,4):datenum(2021,12,11);
% ila.cfg.days2run = datenum(2021,12,12):datenum(2021,12,19);
% ila.cfg.days2run = datenum(2021,12,20):datenum(2022,1,8);

%% TSG
% ila.instrument.TSG.logger = 'Matlab';
% ila.cfg.days2run = datenum(2020,12,12):datenum(2021,5,9);
% ila.instrument.TSG.logger = 'TeraTerm';
% ila.cfg.days2run = datenum(2021,1,10):datenum(2021,1,15);

%% ACS57
% ila.cfg.days2run = datenum(2020,12,26):datenum(2021,1,1);
% ila.cfg.days2run = datenum(2021,1,2):datenum(2021,1,7);
% ila.cfg.days2run = datenum(2021,1,8):datenum(2021,1,15);
% ila.cfg.days2run = datenum(2021,1,16):datenum(2021,1,23);
% ila.cfg.days2run = datenum(2021,1,24):datenum(2021,2,5);
% ila.cfg.days2run = datenum(2021,2,17):datenum(2021,2,24);
% ila.cfg.days2run = datenum(2021,2,25):datenum(2021,3,5);
% ila.cfg.days2run = datenum(2021,3,6):datenum(2021,3,13);
% ila.cfg.days2run = datenum(2021,3,14):datenum(2021,3,21);
% ila.cfg.days2run = datenum(2021,3,22):datenum(2021,3,29);
% ila.cfg.days2run = datenum(2021,3,30):datenum(2021,4,6);
% ila.cfg.days2run = datenum(2021,4,7):datenum(2021,4,14);
% ila.cfg.days2run = datenum(2021,4,15):datenum(2021,4,22);
% ila.cfg.days2run = datenum(2021,4,23):datenum(2021,4,30);
% ila.cfg.days2run = datenum(2021,5,1):datenum(2021,5,9);

% ila.cfg.days2run = datenum(2021,8,18):datenum(2021,8,26);
% ila.cfg.days2run = datenum(2021,8,27):datenum(2021,9,10);
ila.cfg.days2run = datenum(2021,9,11):datenum(2021,10,1);
% ila.cfg.days2run = datenum(2021,10,2):datenum(2021,10,18);
% ila.cfg.days2run = datenum(2021,10,19):datenum(2021,11,4);
% ila.cfg.days2run = datenum(2021,11,5):datenum(2021,12,1);
% ila.cfg.days2run = datenum(2021,12,2):datenum(2022,1,1);
% ila.cfg.days2run = datenum(2022,1,2):datenum(2022,1,7);

%% BB31502
% ila.cfg.days2run = datenum(2020,12,26):datenum(2021,1,5);
% ila.cfg.days2run = datenum(2021,1,6):datenum(2021,1,20);
% ila.cfg.days2run = datenum(2021,1,21):datenum(2021,2,5);
% ila.cfg.days2run = datenum(2021,2,17):datenum(2021,3,5);
% ila.cfg.days2run = datenum(2021,3,6):datenum(2021,3,26);
% ila.cfg.days2run = datenum(2021,3,26):datenum(2021,4,14);
% ila.cfg.days2run = datenum(2021,4,15):datenum(2021,5,9);

ila.cfg.days2run = datenum(2020,12,26):datenum(2021,5,9);

%% HyperBB8005
% ila.cfg.days2run = datenum(2021,8,18):datenum(2021,9,8);
% ila.cfg.days2run = datenum(2021,9,9):datenum(2021,10,9);
% ila.cfg.days2run = datenum(2021,10,10):datenum(2021,11,3);
% ila.cfg.days2run = datenum(2021,11,4):datenum(2021,11,27);
% ila.cfg.days2run = datenum(2021,11,28):datenum(2022,1,7);

%% WSCD859
% ila.cfg.days2run = datenum(2020,12,12):datenum(2021,5,9);

%% SUVF
ila.cfg.days2run = datenum(2021,8,18):datenum(2022,1,7);
 
%% LISST
%  ila.cfg.days2run = datenum(2021,5,1):datenum(2021,5,7);

%% LISST-Tau
ila.cfg.days2run = datenum(2022,5,9):datenum(2022,5,9);

%% LISST200X2002
ila.cfg.days2run = datenum(2022,5,1):datenum(2022,5,31);

%% HBB
%  ila.cfg.days2run = datenum(2021,5,1):datenum(2021,5,7);

%% ALPHA
% ila.cfg.days2run = datenum(2020,12,12):datenum(2021,5,9);

%% 'ALPHA','NMEA','FLOW','ACS57','TSG','BB31502','HBB','WSCD859','PAR','SUVF','LISST1183','NMEA','LISSTTau1002G', 'LISST200X2002'
% Easier to process the FLOW with one other instrument at the time.
% Always include the FLOW even if it was already processed
% ila.cfg.instruments2run = {'FLOW'};
ila.cfg.instruments2run = {'FLOW','LISST200X2002'};
ila.cfg.qcref.view = 'LISST200X2002';
ila.cfg.parallel = Inf;
ila.cfg.calibrate.(ila.cfg.qcref.view).compute_dissolved = false;

%% 1. Import | Load raw data
ila.cfg.force_import = false;
ila.ReadRaw();
ila.CheckDataStatus();

%% Or Load data from already processed mat files
ila.Read('raw');
ila.Read('bin');
ila.Read('qc');
ila.Read('prod');

%% 2. Synchronise instruments
% % % Independent of flow rate (for now)
% % % If flow rate varies use the Strech method
% % % Play with delay of synchronisation
% % % TSG is assumed to be set at zero
% % % No noticeable difference was observed between the TSG of EXPORTS and the BB3
% ila.instrument.(ila.cfg.qcref.view).Sync(-90);
% % % Quick visualizzation to sync with TSG
% % fig(30, 'sync TSG');
% % yyaxis('left'); plot(ila.instrument.TSG.data.dt, ila.instrument.TSG.data.t); ylabel('Temperature (^o C)');
% % % yyaxis('right'); plot(ila.instrument.BB3.data.dt, ila.instrument.BB3.data.beta(:,2)); ylabel('\beta (m^{-1} sr^{-1})'); ylim([80 300]);
% % % datetick2_doy();
% % visSync(ila.instrument.BB3.data, ila.instrument.TSG.data.dt, ila.instrument.TSG.data.t, 'Temp (C)');
% visSync(ila.instrument.FLOW.data, ila.instrument.SUVF.data.dt, ila.instrument.SUVF.data.fdom, 'FDOM (counts)');
% visSync(ila.instrument.FLOW.data, ila.instrument.ACS57.data.dt, ila.instrument.ACS57.data.a(:,20), 'a (m^{-1})');
% visSync(ila.instrument.FLOW.data, ila.instrument.ACS57.data.dt, ila.instrument.ACS57.data.c(:,40), 'c (m^{-1})');
% visSync(ila.instrument.FLOW.data, ila.instrument.HBB.data.dt, ila.instrument.HBB.data.beta(:,14), '\beta (counts)');
% visSync(ila.instrument.FLOW.data, ila.instrument.BB31502.data.dt, ila.instrument.BB31502.data.beta(:,1), '\beta (counts)');
% visSync(ila.instrument.FLOW.data, ila.instrument.LISST1183.data.dt, ila.instrument.LISST1183.data.beta(:,10), '\beta (counts)');
% visSync(ila.instrument.FLOW.data, ila.instrument.WSCD859.data.dt, ila.instrument.WSCD859.data.fdom, 'FDOM (counts)');
% visSync(ila.instrument.FLOW.data, ila.instrument.ALFA.data.dt, ila.instrument.ALFA.data.Chlb, 'chlb');yyaxis('left'); ylim([0 2]);
% % % 
% % % % xlim([datenum(2018,08,14,9,55,0) datenum(2018,08,14,11,05,0)]);
% % % % ylim([-0.1 0.2]);
% % % % Once settings are good set them in the configuration file.
% % % % The software is now doing the same with one line of code.
% % ila.Sync()
% % % % ila.instrument.BB31502.Sync(-90);
% % % % ila.instrument.BB31502.Sync(-10);

%% 2. Auto-synchronise: automatic detection of filter events for AC and BB sensors
% ila.cfg.qcref.MinFiltPeriod = 65; % filter even period in minute % ACS: 55 % BB3: 60
% ila.cfg.qcref.szFilt = 10; % filter even length in minute % default = 10
% ila.SplitDetect(ila.cfg.qcref.MinFiltPeriod, ila.cfg.qcref.szFilt);

%% 3. QC Reference
% run with mode ui during first run (it saves your work for the next run)
% run with mode load to load previous QC
% Note: when redoing QC of a given period of time (days2run) the previous
% QC during the same period of time is erased, QC done on other periods of
% time is kept in the json file
ila.cfg.qcref.mode='ui'; % 'ui' or 'load'
ila.cfg.qcref.remove_old = false; % remove old selection of the same period
ila.QCRef();

%% 4. Split fsw and tsw
ila.Split();
ila.CheckDataStatus();

%% 4.1. SpectralQC
% check raw spectrums AC or BB sensors
ila.SpectralQC('AC',{'raw'}); % AC or BB

%% 5. Automatic QC of raw data for step in ACS spectrum, BB saturated and obvious bad PAR & ALFA values
% fudge factor for auto QC ACS.
% Varies between ACS: 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.RawAutoQCLim.filtered.a = 'auto';
ila.cfg.qc.RawAutoQCLim.filtered.c = 'auto';
ila.cfg.qc.RawAutoQCLim.total.a = 'auto';
ila.cfg.qc.RawAutoQCLim.total.c = 'auto';
% define saturation threshold of a and c in uncalibrated m^-1
ila.cfg.qc.AutoQC_Saturation_Threshold.a = 10; % remove any spectra > threshold m^-1 (uncalibrated)
ila.cfg.qc.AutoQC_Saturation_Threshold.c = 30; % remove any spectra > threshold m^-1 (uncalibrated)
% Tolerance factor for auto QC BB
% 0.1 = minimum tolerance and >> 10 = very high tolerance (default = 3)
ila.cfg.qc.AutoQC_tolerance.filtered.bb = 10; % 10
ila.cfg.qc.AutoQC_tolerance.total.bb = 100; % 10
% define saturation threshold of beta in counts
ila.cfg.qc.AutoQC_Saturation_Threshold.bb = 4100; % saturate above 4000 counts
ila.AutoQC('raw');
ila.CheckDataStatus();

%% 5.1. SpectralQC
% check raw spectrums AC or BB sensors
ila.SpectralQC('AC',{'raw'}); % AC or BB

%% 5.2. Run QC directly on spectra at any level
% ila.SpectralQC inputs:
% 1) 'AC or BB'
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.SpectralQC('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB or BB3:  ila.SpectralQC('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'g','ag'})
ila.SpectralQC('AC',{'raw'}, false, {'tsw','all'});

%% 5.3. Loading previous qc pick selection at raw level
ila.cfg.qc.mode='load';  % load or ui
ila.cfg.qc.specific.run = {ila.cfg.qcref.view}; % 'FLOW','ACS57','TSG', 'BB31502', 'WSCD859','PAR'
ila.QC();

%% 5.4. Write clean raw after split for BB3 and HBB | write only 'part' or 'diw' or 'all'
ila.Write('raw', 'part')
ila.CheckDataStatus();

%% 6. Bin
% % Set settings directly in configuration file (no tunning at this step)
ila.cfg.bin.skip = {};
ila.Bin()
ila.CheckDataStatus();

%% 6.1. SpectralQC
% check binned spectrums AC or BB sensors
ila.SpectralQC('AC',{'bin'}); % AC or BB

%% 6.2. Write bin | write only 'part' or 'diw' or 'all'
ila.Write('bin', 'part')
ila.CheckDataStatus();

%% 7. Flag
ila.Flag() % Now deprecated will just copy data to next level
ila.CheckDataStatus();

%% 8. QC Interactive or Loading previous qc selection
%%%%% Settings %%%%%
ila.cfg.qc.mode='ui';  % load or ui
ila.cfg.qc.remove_old = false;  % remove old selection of this period
ila.cfg.qc.qc_once_for_all = true; % true = QC all variables | false = QC variables separately)
ila.cfg.qc.remove_when_flow_below = false; % true = remove data when flow <= 0.5 | false = no data removal data depending on flow | number = remove data when flow <= number)
% Global
ila.cfg.qc.global.view = {ila.cfg.qcref.view};
ila.cfg.qc.global.active = false;
% Specific
ila.cfg.qc.specific.run = {ila.cfg.qcref.view};
%%%%%%%%%%%%%%%%%%%

% QCmap(ila.cfg.days2run); % plot SST maps to help QC in coastal waters
ila.QC();
ila.CheckDataStatus();

%% 8.1. Auto QC at level 'qc': run until it stabilize to 0
ila.RawAutoQC('qc');

%% 8.2. SpectralQC
% check QCed spectrums AC or BB sensors
ila.SpectralQC('AC',{'qc'}); % AC or BB

%% 8.3. Run QC directly on spectra at any level
% ila.SpectralQC inputs:
% 1) 'AC or BB'
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.SpectralQC('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB or BB3:  ila.SpectralQC('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'g','ag'})
ila.SpectralQC('AC',{'qc'}, false, {'tsw','a'});
ila.SpectralQC('BB',{'qc'}, false, {'tsw','all'});

%% 9. QC Switch position
% QC switch position to make sure each filter event is separated by a
% period of total water and eventually move filter events
ila.QCSwitchPosition()

%% 9.1. Write qc | write only 'part' or 'diw' or 'all'
ila.Write('qc', 'part')

%% 10. Calibrate
% ila.cfg.calibrate.skip = {'FLOW', 'TSG', 'ALFA', 'NMEA'};
% update filter event calcualtion method if needed: exponential_fit 25percentil
ila.cfg.calibrate.(ila.cfg.qcref.view).filt_method = 'exponential_fit'; 
% update filter interpolation method if needed: CDOM linear
ila.cfg.calibrate.(ila.cfg.qcref.view).interpolation_method = 'CDOM';
ila.cfg.min_nb_pts_per_cluster = 100;
ila.cfg.time_weight_for_cluster = 0.9;
% update scattering correction method if needed: Rottgers2013_semiempirical Zaneveld1994_proportional
ila.cfg.calibrate.(ila.cfg.qcref.view).scattering_correction = 'Rottgers2013_semiempirical';
ila.Calibrate();
ila.CheckDataStatus()

%% 10.1 Product visualisation plots with option to save
save_figures = false;

%%% AC or BB 3D plots %%%
ila.SpectralQC('BB', {'prod'}, save_figures); % AC or BB

%%% ACS BB3 TSG PAR WSCD SUVF ALFA LISST final product visualisation %%%
ila.visProd_timeseries()

%% 11. Run QC directly on spectra at any level
% ila.SpectralQC inputs:
% 1) 'AC or BB'
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.SpectralQC('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB or BB3:  ila.SpectralQC('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'g','ag'})
ila.SpectralQC('BB',{'prod'}, false, {'p','all'});

%% 11.1. Load previous qc pick selection at prod level
ila.cfg.qc.mode = 'load';  % load or ui
ila.cfg.qc.specific.run = {ila.cfg.qcref.view}; % 'FLOW','ACS57','TSG', 'BB31502', 'WSCD859','PAR'
ila.QC();

%% 12. Save products | write only 'part' or 'diw' or 'all'
ila.Write('prod', 'part')

% % Notify with a song that the job is done
% notif_sound = load('gong'); sound(notif_sound.y, notif_sound.Fs); % handel
% return

%% re-write final version of 'raw', 'qc' and 'bin' | write only 'part' or 'diw' or 'all'
ila.Write('raw', 'part')
ila.Write('bin', 'part')
ila.Write('qc', 'part')











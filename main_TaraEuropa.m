% Main Particulate InLine Anal4sis Script
% author: Guillaume Bourdin
% created: Jun 27, 2023
clear
% close all

% cd('/Volumes/Data/TaraEuropa/InLineAnalysis-master')
cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master')

% Load InLineAnalysis and the configuration
ila = InLineAnalysis('cfg/TaraEuropa_cfg.m');

% Quick cfg update
%% set the date to process
% ila.cfg.days2run = datetime(2023,4,4):datetime(2023,4,19); %datetime(2023,5,10):datetime(2023,5,30);datetime(2023,6,6):datetime(2023,6,30);datetime(2023,7,8):datetime(2023,8,3)
% ila.cfg.days2run = datetime(2023,4,9):datetime(2023,4,19);
% ila.cfg.days2run = datetime(2023,4,20):datetime(2023,5,26);
% ila.cfg.days2run = datetime(2023,5,2):datetime(2023,5,3);
% ila.cfg.days2run = datetime(2023,4,4):datetime(2023,4,6);

% ila.cfg.days2run = datetime(2023,8,24):datetime(2023,8,24);

%%% entire cruise
% ila.cfg.days2run = datetime(2023,4,3):datetime(2024,8,22);

%%% by ACS
% ila.cfg.days2run = datetime(2023,4,3):datetime(2023,11,16);
ila.cfg.days2run = datetime(2024,2,19):datetime(2024,8,22);




% ila.cfg.days2run = datetime(2023,5,1):datetime(2023,5,5);

%% %%%%%%%%%%%%%% PROCESSING CHRONOLOGY RECOMMENDATIONS: %%%%%%%%%%%%%%% %%
% Instruments available: 'NMEA','FLOW','SBE384504970269','SBE384504970286','SUVF6244','ACS3','ACS348','HyperBB8005','LISST100X1183','LISST200X9999','QCR2150A50351'
%%% Run ReadRaw week by week without going further (reloading 'ila' structure between each run to clear memory)
% ila.cfg.instruments2run = {'FLOW','NMEA','SUVF6244','SBE3845TN444','LISST100X1183','HyperBB8005','ACS3','ACS348'};
%%% Process to the end
% ila.cfg.instruments2run = {'FLOW','NMEA','SBE384504970269','SBE384504970286'};
% ila.cfg.instruments2run = {'FLOW','SUVF6244'};
%%% Process each of the following up to qc level (just before Calibrate and save temporay files raw/bin/qc)
% ila.cfg.instruments2run = {'FLOW','ACS3'};
% ila.cfg.instruments2run = {'FLOW','ACS348'};
% ila.cfg.instruments2run = {'FLOW','HyperBB8005'};
% ila.cfg.instruments2run = {'FLOW','LISST100X1183'};
%%% (Optional if DIW runs available) Process ACS DIW up to DIW bin level (just before DIW Calibrate and save temporay files)
% ila.cfg.instruments2run = {'FLOW','ACS3','SUVF6244','SBE384504970269'};
% ila.cfg.instruments2run = {'FLOW','ACS348','SUVF6244','SBE384504970269'};

%%% Once everything is ready to Calibrate, reload ila structure to clear memory
%%% load entire cruise SUVF/TSG prods and ACS qc and run Calibrate on all at once, save particulate prods
ila.cfg.instruments2run = {'SUVF6244','SBE384504970269','SBE384504970286'}; % merge two TSG into one to run calibrate on entire cruise at once
ila.Read('prod');
% ila.cfg.instruments2run = {'FLOW','ACS3'};
% ila.Read('qc');
ila.cfg.instruments2run = {'FLOW','ACS348'};
ila.Read('qc');
%%% run ACS DIW Calibrate (see mainDI code), save dissolved prods
%%% load HyperBB raw and qc and run Calibrate on all at once, save particulate prods
% ila.cfg.instruments2run = {'SUVF6244','SBE384504970269','SBE384504970286','ACS3'}; % merge two TSG into one to run calibrate on entire cruise at once
% ila.Read('prod');
% ila.cfg.instruments2run = {'FLOW','HyperBB8005'};
% ila.Read('raw');
% ila.Read('qc');
%%% load LISST qc and run Calibrate on all at once, save particulate prods
% ila.cfg.instruments2run = {'FLOW','LISST100X1183'};

%%% Process PAR: entire cruise at once until prod
% ila.cfg.instruments2run = {'QCR2150A50351'};

ila.cfg.qcref.view = 'ACS348';
ila.cfg.parallel = 6; % Inf
ila.cfg.calibrate.(ila.cfg.qcref.view).compute_dissolved = false;
ila.cfg.qc.specific.run = {ila.cfg.qcref.view};

%% 1. Import | Load raw data
ila.cfg.force_import = false;
ila.ReadRaw();
ila.CheckDataStatus();

%% Or Load data from already processed mat files if needed
% ila.Read('raw');
% ila.Read('bin');
% ila.Read('qc');
% ila.Read('prod');
% ila.CheckDataStatus();

%% (Optional) Required only when multiple TSG are used and need to be merge to run calibrate on all data of a cruise
% ila.instrument.SBE384504970269.prod.a = [ila.instrument.SBE384504970269.prod.a; ila.instrument.SBE384504970286.prod.a];
% ila.instrument.SBE384504970269.prod.a = sortrows(ila.instrument.SBE384504970269.prod.a, 'dt');
% ila.instrument.SBE384504970286.prod.a = ila.instrument.SBE384504970269.prod.a;

%% 2. (Optional) Synchronise instruments
% % % Independent of flow rate (for now)
% % % If flow rate varies use the Strech method
% % % Play with delay of synchronisation
% % % TSG is assumed to be set at zero
% % % No noticeable difference was observed between the TSG of EXPORTS and the BB3
% % % ila.instrument.FLOW.Sync(seconds(30));
% ila.instrument.TSG.Sync(seconds(0));
% ila.instrument.SUVF.Sync(seconds(0));
% ila.instrument.ACS57.Sync(seconds(0));
% ila.instrument.HBB.Sync(seconds(0));
% ila.instrument.BB31502.Sync(seconds(0));
% ila.instrument.LISST1183.Sync(seconds(0));
% ila.instrument.WSCD859.Sync(seconds(0));
% ila.instrument.ALFA.Sync(seconds(0)); 
% % % Quick visualizzation to sync with TSG
% % fig(30, 'sync TSG');
% % yyaxis('left'); plot(ila.instrument.TSG.data.dt, ila.instrument.TSG.data.t); ylabel('Temperature (^o C)');
% % % yyaxis('right'); plot(ila.instrument.BB3.data.dt, ila.instrument.BB3.data.beta(:,2)); ylabel('\beta (m^{-1} sr^{-1})'); ylim([80 300]);
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
% % % % xlim([datetime(2018,08,14,9,55,0) datetime(2018,08,14,11,05,0)]);
% % % % ylim([-0.1 0.2]);
% % % % Once settings are good set them in the configuration file.
% % % % The software is now doing the same with one line of code.
% % ila.Sync()
% % % % ila.instrument.BB31502.Sync(seconds(-90));
% % % % ila.instrument.BB31502.Sync(seconds(-10));

%% 2. (Optional) Auto-synchronise: automatic detection of filter events for AC and BB sensors
% ila.cfg.qcref.MinFiltPeriod = minutes(60); % filter even period in minute % ACS: 55 % BB3: 60
% ila.cfg.qcref.szFilt = minutes(10); % filter even length in minute % default = 10
% ila.SplitDetect(ila.cfg.qcref.MinFiltPeriod, ila.cfg.qcref.szFilt);

%% 3. QC Reference: Check filter event position
% run with mode ui during first run (it saves your work for the next run)
% run with mode load to load previous QC
% Note: when redoing QC of a given period of time (days2run) the previous
% QC during the same period of time is erased, QC done on other periods of
% time is kept in the json file
ila.cfg.qcref.mode = 'ui'; % 'ui' or 'load'
ila.cfg.qcref.remove_old = false; % clear old selection of the same period
ila.QCRef();

%% 4. Split fsw and tsw
ila.Split();
ila.CheckDataStatus();

%% 4.1. Spectral QC
% check raw spectrums AC | BB | LISST sensors
ila.SpectralQC('AC',{'raw'});

%% 5. Automatic QC of raw data for step in ACS spectrum, spikes in BB and LISST, saturated data, and obvious bad PAR & ALFA values
% Tolerance factor for auto QC ACS.
% Varies between ACS: 0.1 = minimum tolerance and >> 10 = very high tolerance (default = 3)
ila.cfg.qc.AutoQC_tolerance.filtered.a = 'auto'; %
ila.cfg.qc.AutoQC_tolerance.filtered.c = 'auto'; %
ila.cfg.qc.AutoQC_tolerance.total.a = 'auto'; %
ila.cfg.qc.AutoQC_tolerance.total.c = 'auto'; %
% define saturation threshold of a and c in uncalibrated m^-1
ila.cfg.qc.AutoQC_Saturation_Threshold.a = 10; % remove any spectra > threshold m^-1 (uncalibrated)
ila.cfg.qc.AutoQC_Saturation_Threshold.c = 40; % remove any spectra > threshold m^-1 (uncalibrated)
% Tolerance factor for auto QC BB
% 0.1 = minimum tolerance and >> 10 = very high tolerance (default = 3)
ila.cfg.qc.AutoQC_tolerance.filtered.bb = 100; % 10
ila.cfg.qc.AutoQC_tolerance.total.bb = 10; % 10
% define saturation threshold of beta in counts
ila.cfg.qc.AutoQC_Saturation_Threshold.bb = 4100; % saturate above 4100 counts
% Tolerance factor for auto QC LISST
% 0.1 = minimum tolerance and >> 10 = very high tolerance (default = 3)
ila.cfg.qc.AutoQC_tolerance.filtered.lisst = 10; % 10
ila.cfg.qc.AutoQC_tolerance.total.lisst = 10; % 10
ila.AutoQC('raw');
ila.CheckDataStatus();

%% 5.1. Spectral QC
% check raw spectrums AC | BB | LISST sensors
ila.SpectralQC('AC',{'raw'}); % AC or BB

%% 5.2. Run QC directly on spectra at any level
% ila.SpectralQC inputs:
% 1) 'AC' | 'BB' | 'LISST' sensors
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.SpectralQC('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB or BB3:  ila.SpectralQC('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'g','ag'})
ila.SpectralQC('AC',{'raw'}, false, {'fsw','c'});

%% 5.3. (Optional) Loading previous qc pick selection at raw level
ila.cfg.qc.mode='ui';  % load or ui
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

%% 6.1. Spectral QC
% check binned spectrums AC | BB | LISST sensors
ila.SpectralQC('AC',{'bin'});

%% 6.2. Write bin | write only 'part' or 'diw' or 'all'
ila.Write('bin', 'part')
ila.CheckDataStatus();

%% 7. Pass2QC
ila.Pass2QC('particulate') % copy data to next level
ila.CheckDataStatus();

%% 8. QC Interactive or Loading previous qc selection
%%%%% Settings %%%%%
ila.cfg.qc.mode='ui';  % load or ui
ila.cfg.qc.remove_old = false;  % remove old selection of this period
ila.cfg.qc.qc_once_for_all = false; % true = QC all variables | false = QC variables separately)
ila.cfg.qc.remove_when_flow_below = 0; % true = remove data when flow <= 0.5 | false = no data removal data depending on flow | number = remove data when flow <= number)
% Global QC
ila.cfg.qc.global.view = {ila.cfg.qcref.view};
ila.cfg.qc.global.active = false;
% Specific
ila.cfg.qc.specific.run = {ila.cfg.qcref.view};
%%%%%%%%%%%%%%%%%%%

% QCmap(ila.cfg.days2run); % plot SST maps to help QC in coastal waters
ila.QC();
ila.CheckDataStatus();

%% 8.1. (Optional) Auto QC at level 'qc': run until it stabilize to 0
% ila.AutoQC('qc');

%% 8.2. Spectral QC
% check QCed spectrums AC | BB | LISST sensors
ila.SpectralQC('AC',{'qc'});

%% 8.3. Run QC directly on spectra at any level
% ila.SpectralQC inputs:
% 1) 'AC' | 'BB' | 'LISST' sensors
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.SpectralQC('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB or BB3:  ila.SpectralQC('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'g','ag'})
% ila.SpectralQC('AC',{'qc'}, false, {'fsw','a'});
ila.SpectralQC('AC',{'qc'}, false, {'tsw','c'});

%% 9. QC Switch position
% QC switch position to make sure each filter event is separated by a
% period of total water and eventually move filter events
ila.QCSwitchPosition()

%% 9.1. Write qc | write only 'part' or 'diw' or 'all'
ila.Write('qc', 'part')

%% 10. Calibrate
% update filter event calcualtion method if needed: exponential_fit 25percentil
ila.cfg.calibrate.(ila.cfg.qcref.view).filt_method = 'exponential_fit'; 
% update filter interpolation method if needed: CDOM linear
ila.cfg.calibrate.(ila.cfg.qcref.view).interpolation_method = 'CDOM';
% update scattering correction method if needed: Rottgers2013_semiempirical Zaneveld1994_proportional Semiempirical_blended1 Semiempirical_blended2 Semiempirical_blended3 
ila.cfg.calibrate.(ila.cfg.qcref.view).scattering_correction = 'Semiempirical_blended2';
ila.Calibrate();
ila.CheckDataStatus();

%% 10.1 Product visualisation plots with option to save
save_figures = false;

%%% AC or BB 3D plots %%%
ila.SpectralQC('BB', {'prod'}, save_figures); % AC or BB

%%% ACS BB3 TSG PAR WSCD SUVF ALFA LISST final product visualisation %%%
ila.visProd_timeseries()

%% 11. (Optional) Run QC directly on spectra at any level
% ila.SpectralQC inputs:
% 1) 'AC' | 'BB' | 'LISST' sensors
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.SpectralQC('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB or BB3:  ila.SpectralQC('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'g','ag'})
ila.SpectralQC('AC',{'prod'}, false, {'p','cp'});

%% 11.1. (Optional) Load previous qc pick selection at prod level
ila.cfg.qc.mode = 'load';  % load or ui
ila.cfg.qc.specific.run = {ila.cfg.qcref.view}; % 'FLOW','ACS57','TSG', 'BB31502', 'WSCD859','PAR'
ila.QC();

%% 12. Save products | write only 'part' or 'diw' or 'all'
ila.Write('prod', 'part')

%% re-write final version of 'raw', 'qc' and 'bin' | write only 'part' or 'diw' or 'all'
ila.Write('raw', 'part')
ila.Write('bin', 'part')
ila.Write('qc', 'part')





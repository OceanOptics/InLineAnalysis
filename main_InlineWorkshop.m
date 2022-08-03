% Main Particulate InLine Analysis Script
% author: Guillaume Bourdin
% created: Jan 05, 2021
clear
close all

cd('PATH_TO_INLINEANALYSIS-MASTER_FOLDER')

% Load InLineAnalysis and the configuration
ila = InLineAnalysis('cfg/InlineWorkshop_cfg.m');

% Quick cfg update
%% set the date to process
ila.cfg.days2run = datenum(2022,8,1):datenum(2022,8,1);

%%
% instrument_to_clear = 'ACS298';
% ila.instrument.(instrument_to_clear).data = [];
% ila.instrument.(instrument_to_clear).raw.tsw = [];
% ila.instrument.(instrument_to_clear).raw.fsw = [];
% ila.instrument.(instrument_to_clear).raw.bad = [];
% ila.instrument.(instrument_to_clear).bin.tsw = [];
% ila.instrument.(instrument_to_clear).bin.fsw = [];
% ila.instrument.(instrument_to_clear).qc.tsw = [];
% ila.instrument.(instrument_to_clear).qc.fsw = [];
% ila.instrument.(instrument_to_clear).prod.p = [];
% ila.instrument.(instrument_to_clear).prod.QCfailed = [];

%% 'ALPHA','NMEA','FLOW','ACS298','ACS348','TSG','BB31052','HBB','WSCD859','PAR','SUVF6254','LISST1183','NMEA'
ila.cfg.instruments2run = {'FLOW', 'SBE4536073', 'NMEA'};
% ila.cfg.instruments2run = {'SUVF6254'};
ila.cfg.qcref.view = 'SBE4536073';
ila.cfg.parallel = Inf;
ila.cfg.calibrate.(ila.cfg.qcref.view).compute_dissolved = false;

%% 1. Import | Load raw data
ila.cfg.force_import = true;
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
% % % ila.instrument.FLOW.Sync(30);
% ila.instrument.TSG.Sync(0);
% ila.instrument.SUVF.Sync(0);
% ila.instrument.ACS57.Sync(0);
% ila.instrument.HBB.Sync(0);
% ila.instrument.BB31502.Sync(0);
% ila.instrument.LISST1183.Sync(0);
% ila.instrument.WSCD859.Sync(0);
% ila.instrument.ALFA.Sync(0); 
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
% visSync(ila.instrument.FLOW.data, ila.instrument.BB31052.data.dt, ila.instrument.BB31502.data.beta(:,1), '\beta (counts)');
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
% ila.cfg.qcref.MinFiltPeriod = 60; % filter even period in minute % ACS: 55 % BB3: 60
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

%% 4.1. Diagnostic Plot
% check raw spectrums AC or BB sensors
ila.DiagnosticPlot('AC',{'raw'}); % AC or BB

%% 5. Automatic QC of raw data for step in ACS spectrum, BB saturated and obvious bad PAR & ALFA values
% fudge factor for auto QC ACS.
% Varies between ACS: 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.RawAutoQCLim.filtered.a = 2; % 6.5 16
ila.cfg.qc.RawAutoQCLim.filtered.c = 5; % 13
ila.cfg.qc.RawAutoQCLim.total.a = 2; % 4 13
ila.cfg.qc.RawAutoQCLim.total.c = 2; % 15 5
% fudge factor for auto QC BB
% 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.RawAutoQCLim.filtered.bb = 10; % 2
ila.cfg.qc.RawAutoQCLim.total.bb = 15; % 19
% remove saturated periods in BB
ila.cfg.qc.Saturation_Threshold_bb = 4100; % saturate above 4000 counts
ila.RawAutoQC('raw');
ila.CheckDataStatus();

%% 5.1. Diagnostic Plot
% check raw spectrums AC or BB sensors
ila.DiagnosticPlot('AC',{'raw'}); % AC or BB

%% 5.2. Run QC directly on spectra at any level
% ila.DiagnosticPlot inputs:
% 1) 'AC or BB'
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.DiagnosticPlot('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.DiagnosticPlot('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB or BB3:  ila.DiagnosticPlot('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.DiagnosticPlot('AC',{'prod'}, false, {'g','ag'})
ila.DiagnosticPlot('AC',{'raw'}, false, {'fsw','c'});

%% 5.3. Loading previous qc pick selection at raw level
ila.cfg.qc.mode='load';  % load or ui
ila.cfg.qc.specific.run = {ila.cfg.qcref.view}; % 'FLOW','ACS57','TSG', 'BB31502', 'WSCD859','PAR'
ila.QC();

%% 5.4. Write clean raw after split for BB3 and HBB | write only 'part' or 'diw' or 'all'
ila.Write('raw', 'part')
ila.CheckDataStatus();

%% 6. Bin
% % Set settings directly in configuration file (no tunning at this step)
% % run before re-bin only to clear qc tables
% ila.instrument.ACS57.qc.tsw = table(); ila.instrument.ACS57.qc.fsw = table();
ila.cfg.bin.skip = {};
ila.Bin()
ila.CheckDataStatus();

%% 6.1. Diagnostic Plot
% check binned spectrums AC or BB sensors
ila.DiagnosticPlot('BB',{'bin'}); % AC or BB

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
ila.cfg.qc.qc_once_for_all = false; % true = QC all variables | false = QC variables separately)
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

%% 8.2. Diagnostic Plot
% check QCed spectrums AC or BB sensors
ila.DiagnosticPlot('AC',{'qc'}); % AC or BB

%% 8.3. Run QC directly on spectra at any level
% ila.DiagnosticPlot inputs:
% 1) 'AC or BB'
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.DiagnosticPlot('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.DiagnosticPlot('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB or BB3:  ila.DiagnosticPlot('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.DiagnosticPlot('AC',{'prod'}, false, {'g','ag'})
ila.DiagnosticPlot('AC',{'qc'}, false, {'fsw','a'});
ila.DiagnosticPlot('BB',{'qc'}, false, {'tsw','all'});

%% 9. QC Switch position
% QC switch position to make sure each filter event is separated by a
% period of total water and eventually move filter events
ila.QCSwitchPosition()

%% 9.1. Write qc | write only 'part' or 'diw' or 'all'
ila.Write('qc', 'part')

%% 10. Calibrate
% ila.cfg.calibrate.skip = {'FLOW', 'TSG', 'ALFA', 'NMEA'};
ila.cfg.calibrate.(ila.cfg.qcref.view).filt_method = '25percentil'; % exponential_fit 25percentil
ila.Calibrate();
ila.CheckDataStatus();

%% 10.1 Product visualisation plots with option to save
save_figures = true;

%%% AC or BB 3D plots %%%
ila.DiagnosticPlot('AC', {'prod'}, save_figures); % AC or BB

%%% ACS BB3 TSG PAR WSCD SUVF ALFA LISST final product visualisation %%%
ila.visProd_timeseries()

%% 11. Run QC directly on spectra at any level
% ila.DiagnosticPlot inputs:
% 1) 'AC or BB'
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.DiagnosticPlot('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.DiagnosticPlot('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB or BB3:  ila.DiagnosticPlot('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.DiagnosticPlot('AC',{'prod'}, false, {'g','ag'})
ila.DiagnosticPlot('AC',{'prod'}, false, {'p','ap'});

%% 11.1. Load previous qc pick selection at prod level
ila.cfg.qc.mode = 'load';  % load or ui
ila.cfg.qc.specific.run = {ila.cfg.qcref.view}; % 'FLOW','ACS57','TSG', 'BB31502', 'WSCD859','PAR'
ila.QC();

%% 12. Save products | write only 'part' or 'diw' or 'all'
ila.Write('prod', 'part')

%% re-write final version of 'raw', 'qc' and 'bin' | write only 'part' or 'diw' or 'all'
ila.Write('raw', 'part')
ila.Write('bin', 'part')
ila.Write('qc', 'part')





% Main Particulate InLine Analysis Script
% author: Guillaume Bourdin
% created: Jun 27, 2023
clear
close all

cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master/')

% Load InLineAnalysis and the configuration
ila = InLineAnalysis('cfg/OO2023_cfg.m');

% Quick cfg update
%% set the date to process
ila.cfg.days2run = datenum(2023,6,27):datenum(2023,6,27);

%% 'SBE4536073','atlasTSG002','atlasTSG003','atlasTSG004','NMEA','FLOW','ACS412','HyperBB8002','BB3349','SUVF6254','WS3S1081','LISST1183'
ila.cfg.instruments2run = {'FLOW','HyperBB8002'};
ila.cfg.qcref.view = 'HyperBB8002';
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
ila.CheckDataStatus();

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

%% 3. QC Reference: Check filter event position
% run with mode ui during first run (it saves your work for the next run)
% run with mode load to load previous QC
% Note: when redoing QC of a given period of time (days2run) the previous
% QC during the same period of time is erased, QC done on other periods of
% time is kept in the json file
ila.cfg.qcref.mode='ui'; % 'ui' or 'load'
ila.cfg.qcref.remove_old = true; % remove old selection of the same period
ila.QCRef();

%% 4. Split fsw and tsw
ila.Split();
ila.CheckDataStatus();

%% 4.1. Diagnostic Plot
% check raw spectrums AC | BB | LISST sensors
ila.DiagnosticPlot('LISST',{'raw'});

%% 5. Automatic QC of raw data for step in ACS spectrum, BB saturated and obvious bad PAR & ALFA values
% fudge factor for auto QC ACS.
% Varies between ACS: 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.RawAutoQCLim.filtered.a = 1.5; %
ila.cfg.qc.RawAutoQCLim.filtered.c = 5; %
ila.cfg.qc.RawAutoQCLim.total.a = 1.5; %
ila.cfg.qc.RawAutoQCLim.total.c = 3.5; %
% fudge factor for auto QC BB
% 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.RawAutoQCLim.filtered.bb = 2; % 10
ila.cfg.qc.RawAutoQCLim.total.bb = 2; % 10
% remove saturated periods in BB
ila.cfg.qc.Saturation_Threshold_bb = 4100; % saturate above 4000 counts
ila.RawAutoQC('raw');
ila.CheckDataStatus();

%% 5.1. Diagnostic Plot
% check raw spectrums AC | BB | LISST sensors
ila.DiagnosticPlot('LISST',{'raw'}); % AC or BB

%% 5.2. Run QC directly on spectra at any level
% ila.DiagnosticPlot inputs:
% 1) 'AC' | 'BB' | 'LISST' sensors
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.DiagnosticPlot('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.DiagnosticPlot('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB or BB3:  ila.DiagnosticPlot('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.DiagnosticPlot('AC',{'prod'}, false, {'g','ag'})
ila.DiagnosticPlot('LISST',{'raw'}, false, {'fsw','all'});

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

%% 6.1. Diagnostic Plot
% check binned spectrums AC | BB | LISST sensors
ila.DiagnosticPlot('LISST',{'bin'});

%% 6.2. Write bin | write only 'part' or 'diw' or 'all'
ila.Write('bin', 'part')
ila.CheckDataStatus();

%% 7. Flag
ila.Flag() % Now deprecated will just copy data to next level
ila.CheckDataStatus();

%% 8. QC Interactive or Loading previous qc selection
%%%%% Settings %%%%%
ila.cfg.qc.mode='ui';  % load or ui
ila.cfg.qc.remove_old = true;  % remove old selection of this period
ila.cfg.qc.qc_once_for_all = false; % true = QC all variables | false = QC variables separately)
% Global
ila.cfg.qc.global.view = {ila.cfg.qcref.view};
ila.cfg.qc.global.active = true;
% Specific
ila.cfg.qc.specific.run = {ila.cfg.qcref.view};
%%%%%%%%%%%%%%%%%%%

% QCmap(ila.cfg.days2run); % plot SST maps to help QC in coastal waters
ila.QC();
ila.CheckDataStatus();

%% 8.1. Auto QC at level 'qc': run until it stabilize to 0
% ila.RawAutoQC('qc');

%% 8.2. Diagnostic Plot
% check QCed spectrums AC | BB | LISST sensors
ila.DiagnosticPlot('LISST',{'qc'});

%% 8.3. Run QC directly on spectra at any level
% ila.DiagnosticPlot inputs:
% 1) 'AC' | 'BB' | 'LISST' sensors
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.DiagnosticPlot('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.DiagnosticPlot('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB or BB3:  ila.DiagnosticPlot('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.DiagnosticPlot('AC',{'prod'}, false, {'g','ag'})
% ila.DiagnosticPlot('AC',{'qc'}, false, {'fsw','a'});
ila.DiagnosticPlot('LISST',{'qc'}, false, {'fsw','c'});

%% 9. QC Switch position
% QC switch position to make sure each filter event is separated by a
% period of total water and eventually move filter events
ila.QCSwitchPosition()

%% 9.1. Write qc | write only 'part' or 'diw' or 'all'
ila.Write('qc', 'part')

%% 10. Calibrate
% ila.cfg.calibrate.skip = {'FLOW', 'TSG', 'ALFA', 'NMEA'};
ila.cfg.calibrate.(ila.cfg.qcref.view).filt_method = 'exponential_fit'; % exponential_fit 25percentil
ila.Calibrate();
ila.CheckDataStatus();

%% 10.1 Product visualisation plots with option to save
save_figures = true;

%%% AC or BB 3D plots %%%
ila.DiagnosticPlot('LISST', {'prod'}, save_figures); % AC or BB

%%% ACS BB3 TSG PAR WSCD SUVF ALFA LISST final product visualisation %%%
ila.visProd_timeseries()

%% 11. Run QC directly on spectra at any level
% ila.DiagnosticPlot inputs:
% 1) 'AC' | 'BB' | 'LISST' sensors
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.DiagnosticPlot('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.DiagnosticPlot('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB or BB3:  ila.DiagnosticPlot('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.DiagnosticPlot('AC',{'prod'}, false, {'g','ag'})
ila.DiagnosticPlot('LISST',{'prod'}, false, {'p','ap'});

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





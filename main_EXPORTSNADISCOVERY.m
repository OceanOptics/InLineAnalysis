% Main Particulate InLine Anal4sis Script
% author: Guillaume Bourdin
% created: Jun 27, 2023
clear
% close all

cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master')
% cd('C:\All_Work\Ocean\EXPORTS\Data\EXPORTS_2021\Guillaume\InLineAnalysis\InLineAnalysis\')

cruise = 'EXPORTSNADISCOVERY';
% Load InLineAnalysis and the configuration
ila = InLineAnalysis(fullfile('cfg', [cruise '_cfg.m']));

% Quick cfg update
%% set the date to process
% ila.cfg.days2run = datetime(2021,5,2):datetime(2021,5,4);
% ila.cfg.days2run = datetime(2021,5,5):datetime(2021,5,7);
% ila.cfg.days2run = datetime(2021,5,8):datetime(2021,5,10);
% ila.cfg.days2run = datetime(2021,5,11):datetime(2021,5,13);
% ila.cfg.days2run = datetime(2021,5,14):datetime(2021,5,16);
% ila.cfg.days2run = datetime(2021,5,17):datetime(2021,5,19);
% ila.cfg.days2run = datetime(2021,5,20):datetime(2021,5,22);
% ila.cfg.days2run = datetime(2021,5,23):datetime(2021,5,25);
% ila.cfg.days2run = datetime(2021,5,26):datetime(2021,5,28);
% ila.cfg.days2run = datetime(2021,5,29):datetime(2021,5,30);

%%% entire cruise
ila.cfg.days2run = datetime(2021,5,2):datetime(2021,5,30);

%% %%%%%%%%%%%%%% PROCESSING CHRONOLOGY RECOMMENDATIONS: %%%%%%%%%%%%%%% %%
% Instruments available: 'FLOW','WSCD201','SBE3845999','ACS298'
%%% Run ReadRaw week by week without going further (reloading 'ila' structure between each run to clear memory)
% ila.cfg.instruments2run = {'FLOW','WSCD201','SBE3845999','ACS298'};
%%% Process to the end
% ila.cfg.instruments2run = {'FLOW','SBE3845999'};
% ila.cfg.instruments2run = {'FLOW','WSCD201'};
%%% Process each of the following up to qc level (just before Calibrate and save temporay files raw/bin/qc)
% ila.cfg.instruments2run = {'FLOW','ACS298'};
% % ila.cfg.instruments2run = {'FLOW','HyperBB8005'};
% % ila.cfg.instruments2run = {'FLOW','LISST100X1183'};
% %%% (Optional if DIW runs available) Process ACS DIW up to DIW bin level (just before DIW Calibrate and save temporay files)
% % ila.cfg.instruments2run = {'FLOW','ACS348','WSCD201','SBE3845999'};

%%% Once everything is ready to Calibrate, reload ila structure to clear memory
%%% load entire cruise SUVF/TSG prods and ACS qc and run Calibrate on all at once, save particulate prods
ila.cfg.instruments2run = {'WSCD201','SBE3845999'}; % merge two TSG into one to run calibrate on entire cruise at once
ila.Read('prod');
ila.cfg.instruments2run = {'FLOW','ACS298'};
ila.Read('qc');
% %%% run ACS DIW Calibrate (see mainDI code), save dissolved prods
% %%% load HyperBB raw and qc and run Calibrate on all at once, save particulate prods
% % ila.cfg.instruments2run = {'WSCD201','SBE3845999','ACS3'}; % merge two TSG into one to run calibrate on entire cruise at once
% % ila.Read('prod');
% % ila.cfg.instruments2run = {'FLOW','HyperBB8005'};
% % ila.Read('raw');
% % ila.Read('qc');
% %%% load LISST qc and run Calibrate on all at once, save particulate prods
% % ila.cfg.instruments2run = {'FLOW','LISST100X1183'};

ila.cfg.qcref.view = 'ACS298';
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
% ila.instrument.ACS298.Sync(seconds(-320));
ila.instrument.FLOW.Sync(seconds(50));
ila.Sync()

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

%% %%%%%%%%%%%%%%%% Remove border data 16-20 and 30-33 %%%%%%%%%%%%%%%% %%
min_foo = minute(ila.instrument.ACS298.data.dt);
ila.instrument.ACS298.data.a(min_foo >= 16 & min_foo <= 20, :) = NaN;
ila.instrument.ACS298.data.c(min_foo >= 16 & min_foo <= 20, :) = NaN;
ila.instrument.ACS298.data.a(min_foo >= 30 & min_foo <= 33, :) = NaN;
ila.instrument.ACS298.data.c(min_foo >= 30 & min_foo <= 33, :) = NaN;

%% 4. Split fsw and tsw
ila.Split();
ila.CheckDataStatus();

%%
% ila.instrument.ACS298.raw.fsw.a(any(ila.instrument.ACS298.raw.fsw.a < -1, 2), :) = NaN;
% ila.instrument.ACS298.raw.tsw.a(any(ila.instrument.ACS298.raw.tsw.a < -1, 2), :) = NaN;

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
ila.cfg.qc.AutoQC_Saturation_Threshold.a = 1; % remove any spectra > threshold m^-1 (uncalibrated)
ila.cfg.qc.AutoQC_Saturation_Threshold.c = 0.6; % remove any spectra > threshold m^-1 (uncalibrated)
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
ila.SpectralQC('AC',{'raw'}, false, {'tsw','c'});

%% 5.3. (Optional) Loading previous qc pick selection at raw level
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
ila.cfg.qc.remove_when_flow_below = 0.1; % true = remove data when flow <= 0.5 | false = no data removal data depending on flow | number = remove data when flow <= number)
% Global QC
ila.cfg.qc.global.view = {ila.cfg.qcref.view};
ila.cfg.qc.global.active = false;
% Specific
ila.cfg.qc.specific.run = {ila.cfg.qcref.view};
%%%%%%%%%%%%%%%%%%%

% QCmap(ila.cfg.days2run); % plot SST maps to help QC in coastal waters
ila.QC();
ila.CheckDataStatus();

%% %%%%%%%%%%%%%%%%% Remove minutes 16-23 and 30 to 36 %%%%%%%%%%%%%%%%% %%
min_fsw_foo = minute(ila.instrument.ACS298.qc.fsw.dt);
min_tsw_foo = minute(ila.instrument.ACS298.qc.tsw.dt);
% remove data
ila.instrument.ACS298.qc.fsw.a(min_fsw_foo >= 16 & min_fsw_foo <= 23, :) = NaN;
ila.instrument.ACS298.qc.fsw.c(min_fsw_foo >= 16 & min_fsw_foo <= 23, :) = NaN;
ila.instrument.ACS298.qc.fsw.a(min_fsw_foo >= 30 & min_fsw_foo <= 36, :) = NaN;
ila.instrument.ACS298.qc.fsw.c(min_fsw_foo >= 30 & min_fsw_foo <= 36, :) = NaN;
ila.instrument.ACS298.qc.tsw.a(min_tsw_foo >= 16 & min_tsw_foo <= 23, :) = NaN;
ila.instrument.ACS298.qc.tsw.c(min_tsw_foo >= 16 & min_tsw_foo <= 23, :) = NaN;
ila.instrument.ACS298.qc.tsw.a(min_tsw_foo >= 30 & min_tsw_foo <= 36, :) = NaN;
ila.instrument.ACS298.qc.tsw.c(min_tsw_foo >= 30 & min_tsw_foo <= 36, :) = NaN;

ila.instrument.ACS298.qc.tsw.a(min_tsw_foo >= 5 & min_tsw_foo <= 8, :) = NaN;
ila.instrument.ACS298.qc.tsw.c(min_tsw_foo >= 5 & min_tsw_foo <= 8, :) = NaN;

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
ila.SpectralQC('AC',{'qc'}, false, {'fsw','c'});

%% 9. QC Switch position
% QC switch position to make sure each filter event is separated by a
% period of total water and eventually move filter events
ila.QCSwitchPosition()
 
%% %%%%%%%%%%%%%%%%% Remove Split different filter periods %%%%%%%%%%%%%%%%% %%
torm = [datetime(2021,5,4,21,21,0) datetime(2021,5,5,21,46,0)
  datetime(2021,5,4,22,21,0) datetime(2021,5,4,22,46,0)
  datetime(2021,5,5,14,17,0) datetime(2021,5,5,14,47,0)
  datetime(2021,5,6,4,17,0) datetime(2021,5,6,4,47,0)
  datetime(2021,5,6,6,17,0) datetime(2021,5,6,6,27,30)
  datetime(2021,5,6,7,17,0) datetime(2021,5,6,7,47,0)
  datetime(2021,5,6,19,18,0) datetime(2021,5,6,19,44,0)
  datetime(2021,5,6,20,17,0) datetime(2021,5,6,20,46,0)
  datetime(2021,5,6,21,16,0) datetime(2021,5,6,21,47,0)
  datetime(2021,5,6,22,18,0) datetime(2021,5,6,22,47,0)
  datetime(2021,5,7,13,19,0) datetime(2021,5,7,13,47,0)
  datetime(2021,5,7,14,18,0) datetime(2021,5,7,14,47,0)];

for i = 1:size(torm,1)
  id = ila.instrument.ACS298.qc.tsw.dt > torm(i,1) & ila.instrument.ACS298.qc.tsw.dt < torm(i,2);
  ila.instrument.ACS298.qc.tsw(id,:) = [];
end

% every hour after until the end
torm = [datetime(2021,5,9,12,32,0) datetime(2021,5,20,16,50,0)
  datetime(2021,5,21,8,50,0) datetime(2021,5,30,23,59,0)];

min_tsw_foo = minute(ila.instrument.ACS298.qc.tsw.dt);
for i = 1:size(torm,1)
  id = min_tsw_foo >= 15 & min_tsw_foo <= 50 & ila.instrument.ACS298.qc.tsw.dt > torm(i,1) & ila.instrument.ACS298.qc.tsw.dt < torm(i,2);
  ila.instrument.ACS298.qc.tsw(id,:) = [];
end

%% 9.1. Write qc | write only 'part' or 'diw' or 'all'
ila.Write('qc', 'part')

%% 10. Calibrate
% update filter event calcualtion method if needed: exponential_fit 25percentil
ila.cfg.calibrate.(ila.cfg.qcref.view).filt_method = 'exponential_fit'; 
% update filter interpolation method if needed: CDOM linear
ila.cfg.calibrate.(ila.cfg.qcref.view).interpolation_method = 'linear';
% update scattering correction method if needed: Rottgers2013_semiempirical Zaneveld1994_proportional Semiempirical_blended1 Semiempirical_blended2 Semiempirical_blended3 
ila.cfg.calibrate.(ila.cfg.qcref.view).scattering_correction = 'Semiempirical_blended2';
ila.Calibrate();
ila.CheckDataStatus();

%% 10.1 Product visualisation plots with option to save
save_figures = true;

%%% AC or BB 3D plots %%%
ila.SpectralQC('AC', {'prod'}, save_figures); % AC or BB

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
ila.SpectralQC('AC',{'prod'}, false, {'p','all'});

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




% Main Particulate InLine Analysis Script
% author: Guillaume Bourdin
% created: Jan 05, 2021
clear
close all
closeVar
cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master/')

% Load InLineAnalysis and the configuration
ila = InLineAnalysis('cfg/KM24SOPACE_cfg.m');

%% deconsolidate TSG to have
% % load('/Volumes/Samsung_T5/Data/KM24SOPACE/prod/KM24SOPACE_InLine_TSG_20201224_20220917_Product_v20240111.mat')
% % tsg = addvars(tsg, tsg.avg_n, 'NewVariableNames', 't1_avg_n', 'After', 'tcal_avg_sd');
% % tsg = addvars(tsg, tsg.avg_n, 'NewVariableNames', 'c1_avg_n', 'After', 'cond_avg_sd');
% % tsg = addvars(tsg, tsg.avg_n, 'NewVariableNames', 's_avg_n', 'After', 'sss_avg_sd');
% % tsg = addvars(tsg, tsg.avg_n, 'NewVariableNames', 'sv_avg_n', 'After', 'sv_avg_sd');
% % path_tocopy = ('/Volumes/Samsung_T5/Data/KM24SOPACE/wk/SBE38450091');
% % variable_names = ila.instrument.SBE38450091.prod.a.Properties.VariableNames;
% % suffix = 'a';
% % var_torm = [{'lat', 'lon'} variable_names(contains(variable_names, 'sv'))];
% % variable_names(contains(variable_names, 'sv')) = [];
% % deconsolidate(tsg, variable_names, var_torm, path_tocopy, suffix);
% 
% load('/Volumes/Samsung_T5/Data/KM24SOPACE/prod/KM24SOPACE_InLine_TSG_20201224_20220917_Product_v20240111.mat')
% tsg = addvars(tsg, tsg.avg_n, 'NewVariableNames', 'tcal_avg_n', 'After', 'tcal_avg_sd');
% tsg = addvars(tsg, tsg.avg_n, 'NewVariableNames', 'c_avg_n', 'After', 'cond_avg_sd');
% tsg = addvars(tsg, tsg.avg_n, 'NewVariableNames', 's_avg_n', 'After', 'sss_avg_sd');
% tsg = renamevars(tsg, {'cond','cond_avg_sd','sss','sss_avg_sd','sst','sst_avg_sd','avg_n'}, {'c','c_avg_sd','s','s_avg_sd','t2','t2_avg_sd','t2_avg_n'});
% ila.instrument.SBE38450091.qc.tsw = tsg;

% Quick cfg update
%% Whole expedition
% ila.cfg.days2run = datetime(2020,12,12):datetime(2022,9,17);

%% TSG
% ila.instrument.TSG.logger = 'Matlab';
% ila.cfg.days2run = datetime(2020,12,12):datetime(2021,5,9);
% ila.instrument.TSG.logger = 'TeraTerm';
% ila.cfg.days2run = datetime(2021,1,10):datetime(2021,1,15);
% ila.instrument.TSG.logger = 'Inlinino_base';
% ila.cfg.days2run = datetime(2021,1,16):datetime(2022,7,15);
% ila.instrument.TSG.logger = 'Inlinino_base';
% ila.cfg.days2run = datetime(2022,7,15):datetime(2022,10,30);

%% Empty ACs data without losing the rest
instrument_to_clear = 'ACS111'; % ACS91 ACS111
ila.instrument.(instrument_to_clear).data = [];
ila.instrument.(instrument_to_clear).raw.tsw = [];
ila.instrument.(instrument_to_clear).raw.fsw = [];
ila.instrument.(instrument_to_clear).raw.bad = [];
ila.instrument.(instrument_to_clear).bin.tsw = [];
ila.instrument.(instrument_to_clear).bin.fsw = [];
ila.instrument.(instrument_to_clear).qc.tsw = [];
ila.instrument.(instrument_to_clear).qc.fsw = [];
ila.instrument.(instrument_to_clear).prod.p = [];
ila.instrument.(instrument_to_clear).prod.QCfailed = [];

%% pre-processing 
% ila.cfg.days2run = datetime(2024,10,19):datetime(2024,10,20);
% ila.cfg.days2run = datetime(2024,10,24):datetime(2024,10,28);
% ila.cfg.days2run = datetime(2024,10,29):datetime(2024,11,7);
% ila.cfg.days2run = datetime(2024,11,8):datetime(2024,11,16);
% ila.cfg.days2run = datetime(2024,11,17):datetime(2024,11,26);
% ila.cfg.days2run = datetime(2024,11,28):datetime(2024,12,5);
% ila.cfg.days2run = datetime(2024,12,6):datetime(2024,12,16);

%% ACS091
% ila.cfg.days2run = datetime(2024,10,19):datetime(2024,10,20);
% ila.cfg.days2run = datetime(2024,10,24):datetime(2024,10,25);
% ila.cfg.days2run = datetime(2024,10,26):datetime(2024,10,27);
% ila.cfg.days2run = datetime(2024,10,28):datetime(2024,10,29);
% ila.cfg.days2run = datetime(2024,10,30):datetime(2024,10,31);
% ila.cfg.days2run = datetime(2024,11,1):datetime(2024,11,2);
% ila.cfg.days2run = datetime(2024,11,3):datetime(2024,11,4);
% ila.cfg.days2run = datetime(2024,11,5):datetime(2024,11,6);
% ila.cfg.days2run = datetime(2024,11,7):datetime(2024,11,8);
% ila.cfg.days2run = datetime(2024,11,9):datetime(2024,11,10);
% ila.cfg.days2run = datetime(2024,11,11):datetime(2024,11,12);
% ila.cfg.days2run = datetime(2024,11,13):datetime(2024,11,14);
% ila.cfg.days2run = datetime(2024,11,15):datetime(2024,11,16);
% ila.cfg.days2run = datetime(2024,11,17):datetime(2024,11,18);
% ila.cfg.days2run = datetime(2024,11,19):datetime(2024,11,20);
% ila.cfg.days2run = datetime(2024,11,21):datetime(2024,11,22);
% ila.cfg.days2run = datetime(2024,11,23):datetime(2024,11,24);
% ila.cfg.days2run = datetime(2024,11,25):datetime(2024,11,26);

% ila.cfg.days2run = datetime(2024,11,28):datetime(2024,11,29);
% ila.cfg.days2run = datetime(2024,11,30):datetime(2024,11,31);
% ila.cfg.days2run = datetime(2024,12,1):datetime(2024,12,2); % ACS091
% ila.cfg.days2run = datetime(2024,12,2):datetime(2024,12,3); % ACS111
% ila.cfg.days2run = datetime(2024,12,4):datetime(2024,12,5);
% ila.cfg.days2run = datetime(2024,12,6):datetime(2024,12,7);
% ila.cfg.days2run = datetime(2024,12,8):datetime(2024,12,9);
% ila.cfg.days2run = datetime(2024,12,10):datetime(2024,12,11);
% ila.cfg.days2run = datetime(2024,12,12):datetime(2024,12,13);
% ila.cfg.days2run = datetime(2024,12,14):datetime(2024,12,16);

%% calibrate runs
% ila.cfg.days2run = datetime(2024,10,23):datetime(2024,11,07); % ACS091
% ila.cfg.days2run = datetime(2024,11,08):datetime(2024,12,2); % ACS091
% ila.cfg.days2run = datetime(2024,12,2):datetime(2024,12,16); % ACS111

%% entire leg 1
ila.cfg.days2run = datetime(2024,10,23):datetime(2024,11,26);

%% entire leg 2
ila.cfg.days2run = datetime(2024,11,28):datetime(2024,12,16);

%% entire leg with ACS091
ila.cfg.days2run = datetime(2024,10,23):datetime(2024,12,2);

%% entire leg with ACS111
ila.cfg.days2run = datetime(2024,12,2):datetime(2024,12,16);

%% %%%%%%%%%%%%%% PROCESSING CHRONOLOGY RECOMMENDATIONS: %%%%%%%%%%%%%%% %%
% Instruments available: 'NMEA','FLOW','SBE3845KM24','ACS91','ACS111','HyperBB8005','SUVF6244','LISST100X1183','WS3S1081P','ALFA'
%%% Run ReadRaw week by week without going further (reloading 'ila' structure between each run to clear memory)
% ila.cfg.instruments2run = {'FLOW','NMEA','SUVF6244','SBE3845KM24','LISST100X1183','HyperBB8005'};
%%% Process to the end
% ila.cfg.instruments2run = {'FLOW','NMEA','SBE3845KM24'};
% ila.cfg.instruments2run = {'FLOW','SUVF6244'};
%%% Process each of the following up to qc level (just before Calibrate and save temporay files raw/bin/qc)
% ila.cfg.instruments2run = {'FLOW','ACS91'}; % 2 days by 2 days
% ila.cfg.instruments2run = {'FLOW','ACS111'}; % 2 days by 2 days
% ila.cfg.instruments2run = {'FLOW','HyperBB8005'}; % 2 weeks by 2 weeks
% ila.cfg.instruments2run = {'FLOW','LISST100X1183'}; % 2 weeks by 2 weeks
%%% (Optional if DIW runs available) Process ACS DIW up to DIW bin level (just before DIW Calibrate and save temporay files)
% ila.cfg.instruments2run = {'FLOW','ACS91'};
% ila.cfg.instruments2run = {'FLOW','ACS111'};

%%% Once everything is ready to Calibrate, reload ila structure to clear memory
%%% load entire cruise SUVF/TSG prods and ACS qc and run Calibrate on all at once, save particulate prods
% ila.cfg.instruments2run = {'SUVF6244','SBE3845KM24'};
% ila.Read('prod');
% ila.cfg.instruments2run = {'FLOW','ACS91'};
% ila.Read('qc');
% ila.cfg.instruments2run = {'FLOW','ACS111'};
% ila.Read('qc');
%%% run ACS DIW Calibrate, save prod dissolved
%%% run ACS DIW Calibrate (see mainDI code), save dissolved prods
%%% load HyperBB raw and qc and run Calibrate on all at once, save particulate prods
ila.cfg.instruments2run = {'SUVF6244','SBE3845KM24','ACS111','ACS91'}; % merge two TSG into one to run calibrate on entire cruise at once
ila.Read('prod');
ila.cfg.instruments2run = {'FLOW','HyperBB8005'};
ila.Read('raw');
ila.Read('qc');
%%% load LISST qc and run Calibrate, save prod dissolved
% ila.cfg.instruments2run = {'FLOW','LISST100X1183'};

%%% Process PAR: entire cruise at once until prod
% ila.cfg.instruments2run = {'QCR2150A50351'};

ila.cfg.qcref.view = 'HyperBB8005';
ila.cfg.parallel = 6; % Inf
ila.cfg.calibrate.(ila.cfg.qcref.view).compute_dissolved = false;
ila.cfg.qc.specific.run = {ila.cfg.qcref.view};

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

%% 2. Auto-synchronise: automatic detection of filter events for AC and BB sensors
% ila.cfg.qcref.MinFiltPeriod = minutes(60); % filter even period in minute % ACS: 55 % BB3: 60
% ila.cfg.qcref.szFilt = minutes(10); % filter even length in minute % default = 10
% ila.SplitDetect(ila.cfg.qcref.MinFiltPeriod, ila.cfg.qcref.szFilt);

%% 3. QC Reference
% run with mode ui during first run (it saves your work for the next run)
% run with mode load to load previous QC
% Note: when redoing QC of a given period of time (days2run) the previous
% QC during the same period of time is erased, QC done on other periods of
% time is kept in the json file
ila.cfg.qcref.mode='load'; % 'ui' or 'load'
ila.cfg.qcref.remove_old = false; % remove old selection of the same period
ila.QCRef();

%% 4. Split fsw and tsw
ila.Split();
ila.CheckDataStatus();

%% 4.1. SpectralQC Plot
% check raw spectrums AC or BB sensors
ila.SpectralQC('LISST',{'raw'}); % AC, BB, or LISST

%% 5. Automatic QC of raw data for step in ACS spectrum, BB saturated, LISST, and obvious bad PAR & ALFA values
% Tolerance factor for auto QC ACS.
% Varies between ACS: 0.1 = minimum tolerance and >> 10 = very high tolerance (default = 3)
ila.cfg.qc.AutoQC_tolerance.filtered.a = 'auto'; % 2
ila.cfg.qc.AutoQC_tolerance.filtered.c = 'auto'; % 2
ila.cfg.qc.AutoQC_tolerance.total.a = 'auto'; % 2
ila.cfg.qc.AutoQC_tolerance.total.c = 'auto'; % 3
% define saturation threshold of a and c in uncalibrated m^-1
ila.cfg.qc.AutoQC_Saturation_Threshold.a = 10; % remove any spectra > threshold m^-1 (uncalibrated)
ila.cfg.qc.AutoQC_Saturation_Threshold.c = 30; % remove any spectra > threshold m^-1 (uncalibrated)
% Tolerance factor for auto QC BB
% 0.1 = minimum tolerance and >> 10 = very high tolerance (default = 3)
ila.cfg.qc.AutoQC_tolerance.filtered.bb = 10; % 10
ila.cfg.qc.AutoQC_tolerance.total.bb = 20; % 10
% define saturation threshold of beta in counts
ila.cfg.qc.AutoQC_Saturation_Threshold.bb = 4100; % saturate above 4100 counts
% Tolerance factor for auto QC LISST
% 0.1 = minimum tolerance and >> 10 = very high tolerance (default = 3)
ila.cfg.qc.AutoQC_tolerance.filtered.lisst = 10; % 10
ila.cfg.qc.AutoQC_tolerance.total.lisst = 10; % 10
ila.AutoQC('raw');
ila.CheckDataStatus();

%% 5.1. SpectralQC Plot
% check raw spectrums AC, BB, or LISST sensors
ila.SpectralQC('AC',{'raw'}); % AC, BB, or LISST

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
ila.SpectralQC('AC',{'raw'}, false, {'tsw','a'});

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

%% 6.1. SpectralQC Plot
% check binned spectrums AC or BB sensors
ila.SpectralQC('AC',{'bin'}); % AC, BB, or LISST

%% 6.2. Write bin | write only 'part' or 'diw' or 'all'
ila.Write('bin', 'part')
ila.CheckDataStatus();

%% 7. Pass2QC
ila.Pass2QC('particulate') % Copy data to next level
ila.CheckDataStatus();

%% 8. QC Interactive or Loading previous qc selection
%%%%% Settings %%%%%
ila.cfg.qc.mode='ui';  % load or ui
ila.cfg.qc.remove_old = false;  % remove old selection of this period
ila.cfg.qc.qc_once_for_all = false; % true = QC all variables | false = QC variables separately)
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
ila.AutoQC('qc');

%% 8.2. SpectralQC Plot
% check QCed spectrums AC or BB sensors
ila.SpectralQC('AC',{'qc'}); % AC, BB, or LISST

%% 8.3. Run QC directly on spectra at any level
% ila.SpectralQC inputs:
% 1) 'AC', 'BB', or 'LIIST'
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.SpectralQC('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB, BB3, or LISST:  ila.SpectralQC('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'g','ag'})
ila.SpectralQC('LISST',{'qc'}, false, {'tsw','beta'});
ila.SpectralQC('AC',{'qc'}, false, {'fsw','c'});
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
ila.cfg.calibrate.(ila.cfg.qcref.view).interpolation_method = 'CDOM'; % linear CDOM
% update scattering correction method if needed: Rottgers2013_semiempirical Zaneveld1994_proportional Semiempirical_blended1 Semiempirical_blended2 Semiempirical_blended3 
ila.cfg.calibrate.(ila.cfg.qcref.view).scattering_correction = 'Semiempirical_blended2';
ila.cfg.calibrate.(ila.cfg.qcref.view).compute_dissolved = false;
ila.Calibrate();
ila.CheckDataStatus()

% % Compute NAAMES specific chl
% wl = ila.instrument.ACS111.lambda_ref; ACS = ila.instrument.ACS111.prod;
% ap_a=interp1(wl,ACS111.p.ap',[650 676 715],'linear')';
% line_height = (ap_a(:,2)-(39/65*ap_a(:,1)+26/65*ap_a(:,3)));
% p.chl=157*line_height.^1.22;
% ila.instrument.ACS111.prod.p.chl_naames = 95 * line_height .^ 1.06;
% ila.instrument.ACS111.prod.p.chl_naames(real(ila.instrument.ACS111.prod.p.chl_naames) ~= ila.instrument.ACS.prod.p.chl_naames) = NaN;
% fprintf('Done\n');
% 
% % Compute EXPORTS Specific chl
% % Derive Chl (Line heigh at 676 compared to 650 and 715)
% ila_acs = ila.instrument.ACS111.prod.p; ila_acs_wl = ila.instrument.ACS111.lambda_ref;
% ap = interp1(ila_acs_wl, ila_acs.ap', [650 676 715], 'linear')';
% line_height = (ap(:,2)-(39/65*ap(:,1)+26/65*ap(:,3)));
% ila.instrument.ACS111.prod.p.chl_exports = 138.14 * line_height.^1.11; % EXPORTS relation

%% 10.1 Product visualisation plots with option to save
save_figures = true;

%%% AC or BB 3D plots %%%
ila.SpectralQC('BB', {'prod'}, save_figures); % AC or BB

%%% ACS BB3 TSG PAR WSCD SUVF ALFA LISST final product visualisation %%%
ila.visProd_timeseries()

%% 11. Run QC directly on spectra at any level
% ila.SpectralQC inputs:
% 1) 'AC', 'BB', or 'LIIST'
% 2) 'level':  'raw' | 'bin' | 'qc' | 'prod'
% 3) save plot option: boolean
% 4) table and variable to QC as shown in examples below
% Examples:
%     - to QC 'a' of 'tsw' table of 'qc' level of ACs: ila.SpectralQC('AC',{'qc'}, false, {'tsw','a'})
%     - to QC 'cp' of 'p' table of 'prod' level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'p','cp'})
%     - to QC 'beta' of 'fsw' table of 'bin' level of HBB, BB3, or LISST:  ila.SpectralQC('BB',{'bin'}, false, {'fsw','beta'})
%     - to QC 'ag' of 'g' table of prod level of ACs:  ila.SpectralQC('AC',{'prod'}, false, {'g','ag'})
ila.SpectralQC('LISST',{'prod'}, false, {'p','betap'});
ila.SpectralQC('AC',{'prod'}, false, {'p','cp'});
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






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'NMEA'}; % {'FLOW', 'TSG', 'BB31502','PAR', 'WSCD859'}
ila.cfg.days2run = datetime(2022,6,6):datetime(2022,6,13);

% populate ila.instrument
ila.Read('prod');

nmea = ila.instrument.NMEA.prod.a;

AC = ila.instrument.ACS57.prod.p;
latlon_interp = interp1(nmea.dt, [nmea.lat, nmea.lon], AC.dt, 'linear'); % extrap needed for first minute of data

data_AC = table(AC.dt, latlon_interp(:,1), latlon_interp(:,2), ...
  AC.ap, AC.ap_sd, AC.cp, AC.cp_sd, ...
             'VariableNames', {'dt', 'lat', 'lon', 'ap', 'ap_sd', 'cp', 'cp_sd'});
           
data_AC = [data_AC AC(:, 8:end) AC(:, 6:7)];
data_AC.Properties.VariableUnits = {'', 'degrees', 'degrees', '1/m', '1/m', '1/m', '1/m', ...
  '1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m', ...
  'ug/L','unitless','ug/L','unitless','unitless','ug/L','microns','unitless','unitless','unitless','unitless'};
data_AC.Properties.VariableDescriptions = [{''} repmat({'%.4f'}, 1, size(data_AC,2) - 4) repmat({'%i'}, 1, 3)];


ila.visProd_timeseries()

SimpleMap(data_AC.chl_Halh, data_AC(:,1:3), 'Houskeeper [chl] (mg.m^{-3})')
SimpleMap(data_AC.HH_G50, data_AC(:,1:3), 'H&H phytoplankton G50: cross-sectional area (\mum)')
SimpleMap(data_AC.poc, data_AC(:,1:3), '[POC] cp (mg.m^{-3})')
SimpleMap(data_AC.chl_ap676lh, data_AC(:,1:3), 'a_{p676} line height [chl a] (mg.m^{-3})')
SimpleMap(data_AC.gamma, data_AC(:,1:3), 'gamma cp (unitless)')



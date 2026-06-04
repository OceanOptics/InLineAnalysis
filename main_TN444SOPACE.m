% Main Particulate InLine Analysis Script
% author: Guillaume Bourdin
% created: Jan 05, 2021
clear
close all
closeVar
cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master/')

% Load InLineAnalysis and the configuration
ila = InLineAnalysis('cfg/TN444SOPACE_cfg.m');

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
% ila.cfg.days2run = datetime(2025,5,5):datetime(2025,5,12);

%% ACS111
% ila.cfg.days2run = datetime(2025,5,5):datetime(2025,5,6);
% ila.cfg.days2run = datetime(2025,5,7):datetime(2025,5,8);
% ila.cfg.days2run = datetime(2025,5,9):datetime(2025,5,10);
% ila.cfg.days2run = datetime(2025,5,11):datetime(2025,5,12);

% ila.cfg.days2run = datetime(2025,5,27):datetime(2025,5,28);
% ila.cfg.days2run = datetime(2025,5,29):datetime(2025,5,30);
% ila.cfg.days2run = datetime(2025,5,30):datetime(2025,5,31);
% ila.cfg.days2run = datetime(2025,6,1):datetime(2025,6,2);
% ila.cfg.days2run = datetime(2025,6,3):datetime(2025,6,4);
% ila.cfg.days2run = datetime(2025,6,5):datetime(2025,6,6);
% ila.cfg.days2run = datetime(2025,6,7):datetime(2025,6,8);
% ila.cfg.days2run = datetime(2025,6,9):datetime(2025,6,10);
% ila.cfg.days2run = datetime(2025,6,11):datetime(2025,6,12);
% ila.cfg.days2run = datetime(2025,6,13):datetime(2025,6,15);

%% calibrate runs
% ila.cfg.days2run = datetime(2025,5,5):datetime(2025,5,12); % ACS111
% ila.cfg.days2run = datetime(2025,5,27):datetime(2025,6,12); % ACS111

%% entire leg 1
% ila.cfg.days2run = datetime(2025,5,5):datetime(2025,5,12);

%% entire leg 2
% ila.cfg.days2run = datetime(2025,5,27):datetime(2025,6,15);

%% entire cruise
ila.cfg.days2run = datetime(2025,5,5):datetime(2025,6,15);

%% %%%%%%%%%%%%%% PROCESSING CHRONOLOGY RECOMMENDATIONS: %%%%%%%%%%%%%%% %%
% Instruments available: 'NMEA','FLOW','SBE3845TN444','ACS91','ACS111','HyperBB8005','SUVF6266','LISST100X1183','WS3S1081P','ALFA11'
%%% Run ReadRaw week by week without going further (reloading 'ila' structure between each run to clear memory)
% ila.cfg.instruments2run = {'FLOW','NMEA','SUVF6266','SBE3845TN444','LISST100X1183','HyperBB8005','ACS111'};
%%% Process to the end
% ila.cfg.instruments2run = {'FLOW','NMEA','SBE3845TN444'};
% ila.cfg.instruments2run = {'FLOW','SUVF6266'};
%%% Process each of the following up to qc level (just before Calibrate and save temporay files raw/bin/qc)
% ila.cfg.instruments2run = {'FLOW','ACS111'};
% ila.cfg.instruments2run = {'FLOW','HyperBB8005'};
% ila.cfg.instruments2run = {'FLOW','LISST100X1183'};
%%% (Optional if DIW runs available) Process ACS DIW up to DIW bin level (just before DIW Calibrate and save temporay files)
% ila.cfg.instruments2run = {'FLOW','ACS111','SUVF6266','SBE3845TN444'};

%%% Once everything is ready to Calibrate, reload ila structure to clear memory
%%% load entire cruise SUVF/TSG prods and ACS qc and run Calibrate on all at once, save particulate prods
ila.cfg.instruments2run = {'SUVF6266','SBE3845TN444'};
ila.Read('prod');
ila.cfg.instruments2run = {'FLOW','ACS111'};
ila.Read('qc');
%%% run ACS DIW Calibrate (see mainDI code), save dissolved prods
%%% load entire cruise SUVF/TSG/ACS prod
% ila.cfg.instruments2run = {'SUVF6266','SBE3845TN444','ACS111'};
% ila.Read('prod');
%%% load HyperBB raw and qc and run Calibrate on all at once, save particulate prods
% ila.cfg.instruments2run = {'FLOW','HyperBB8005'};
%%% load LISST qc and run Calibrate on all at once, save particulate prods
% ila.cfg.instruments2run = {'FLOW','LISST100X1183'};

ila.cfg.qcref.view = 'ACS111';
ila.cfg.parallel = 6; % Inf
ila.cfg.calibrate.(ila.cfg.qcref.view).compute_dissolved = false;
ila.cfg.qc.specific.run = {ila.cfg.qcref.view};

%%
% figure; scatter(diff(ila.instrument.ACS111.prod.g.cg(:,20)), diff(ila.instrument.ACS111.prod.g.fdom))
% ila.instrument.ACS111.prod.g = table();

%% 1. Import | Load raw data
ila.cfg.force_import = false;
ila.ReadRaw();
ila.CheckDataStatus();

%% Or Load data from already processed mat files
ila.Read('raw');
ila.Read('bin');
ila.Read('qc');
ila.Read('prod');

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

%% 2. (Optional: if no or bad switch data) Auto-synchronise: automatic detection of filter events for AC and BB sensors
% ila.cfg.qcref.MinFiltPeriod = minutes(60); % filter even period in minute % ACS: 55 % BB3: 60
% ila.cfg.qcref.szFilt = minutes(10); % filter even length in minute % default = 10
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

%% 4.1. SpectralQC Plot
% check raw spectrums AC or BB sensors
ila.SpectralQC('LISST',{'raw'}); % AC, BB, or LISST

%%
% ila.instrument.ACS111.raw.fsw.a(any(ila.instrument.ACS111.raw.fsw.a > 2.3,2),:) = NaN;
% ila.instrument.ACS111.raw.fsw.c(any(ila.instrument.ACS111.raw.fsw.c > 3,2),:) = NaN;

%% 5. Automatic QC of raw data for step in ACS spectrum, BB saturated, LISST, and obvious bad PAR & ALFA values
% Tolerance factor for auto QC ACS.
% Varies between ACS: 0.1 = minimum tolerance and >> 10 = very high tolerance (default = 3)
ila.cfg.qc.AutoQC_tolerance.filtered.a = 'auto'; % 2 'auto'
ila.cfg.qc.AutoQC_tolerance.filtered.c = 'auto'; % 2 'auto'
ila.cfg.qc.AutoQC_tolerance.total.a = 'auto'; % 2 auto
ila.cfg.qc.AutoQC_tolerance.total.c = 'auto'; % 3 auto
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
ila.cfg.qc.AutoQC_tolerance.filtered.lisst = 30; % 10
ila.cfg.qc.AutoQC_tolerance.total.lisst = 100; % 10
ila.AutoQC('raw');
ila.CheckDataStatus();

%%


% id650 = abs(650 - ila.instrument.ACS111.lambda_a) == min(abs(650 - ila.instrument.ACS111.lambda_a));
% ila.instrument.ACS111.raw.fsw.a(ila.instrument.ACS111.raw.fsw.a(:,id650) > 0.2,:) = NaN;
% 
% wl2clean = 490;
% id = abs(wl2clean - ila.instrument.ACS111.lambda_a) == min(abs(wl2clean - ila.instrument.ACS111.lambda_a));
% figure(1); histogram(ila.instrument.ACS111.raw.fsw.a(:,id))
% ila.instrument.ACS111.raw.fsw.a(ila.instrument.ACS111.raw.fsw.a(:,id) > 0.76,:) = NaN;


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
ila.SpectralQC('AC',{'raw'}, false, {'fsw','a'});
ila.SpectralQC('LISST',{'raw'}, false, {'tsw','all'});

%% 5.3. Loading previous qc pick selection at raw level
ila.cfg.qc.mode='load';  % load or ui
ila.cfg.qc.specific.run = {ila.cfg.qcref.view}; % 'FLOW','ACS57','TSG', 'BB31502', 'WSCD859','PAR'
ila.QC();

%% 5.4. Write clean raw after split | write only 'part' or 'diw' or 'all'
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
ila.SpectralQC('BB',{'bin'}); % AC, BB, or LISST

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
ila.cfg.qc.remove_when_flow_below = 0.5; % true = remove data when flow <= 0.5 | false = no data removal data depending on flow | number = remove data when flow <= number)
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
ila.SpectralQC('BB',{'qc'}); % AC, BB, or LISST

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
ila.SpectralQC('AC',{'qc'}, false, {'tsw','c'});
ila.SpectralQC('BB',{'qc'}, false, {'tsw','all'});

%% 9. QC Switch position
% QC switch position to make sure each filter event is separated by a
% period of total water and eventually move filter events
ila.QCSwitchPosition()

% %% cut day-long filter event into 30 min pieces (qc level: all instruments except hyperBB)
% idf = find(ila.instrument.FLOW.qc.tsw.dt > datetime(2025,5,10,15,0,0) & ila.instrument.FLOW.qc.tsw.dt < datetime(2025,5,11,12,0,0));
% ila.instrument.FLOW.qc.tsw.(ila.instrument.FLOW.view.swt_variable)(idf(30:30:size(idf,1))) = false;
% 
% %% cut day-long filter event into 30 min pieces (raw level: hyperBB)
% idf = find(ila.instrument.FLOW.raw.tsw.dt > datetime(2025,5,10,15,0,0) & ila.instrument.FLOW.raw.tsw.dt < datetime(2025,5,11,12,0,0));
% ila.instrument.FLOW.raw.tsw.(ila.instrument.FLOW.view.swt_variable)(idf(1800:1800:size(idf,1))) = false;

%% 9.1. Write qc | write only 'part' or 'diw' or 'all'
ila.Write('qc', 'part')

%% 10. Calibrate
% ila.cfg.calibrate.skip = {'FLOW', 'TSG', 'ALFA', 'NMEA'};
% update filter event calcualtion method if needed: exponential_fit 25percentil
ila.cfg.calibrate.(ila.cfg.qcref.view).filt_method = 'exponential_fit';
% update filter interpolation method if needed: CDOM linear
ila.cfg.calibrate.(ila.cfg.qcref.view).interpolation_method = 'linear'; % linear CDOM
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

%%% AC, BB, or LISST 3D plots %%%
ila.SpectralQC('LISST', {'prod'}, save_figures); % AC, BB, or LISST

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
ila.SpectralQC('LISST',{'prod'}, false, {'p','all'});
ila.SpectralQC('AC',{'prod'}, false, {'p','ap'});
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



%% check and correct rust fouling on a and c
wl_a = ila.instrument.ACS111.lambda_a;
wl_c = ila.instrument.ACS111.lambda_c;
st1 = importInlininoACScsv('/Volumes/Extreme SSD/TN444-SOPACE/raw/ACS111/St1_ACS111_20250505_063722.csv');
st1_2 = importInlininoACScsv('/Volumes/Extreme SSD/TN444-SOPACE/raw/ACS111/St1_2_ACS111_20250505_065811.csv');
st2 = importInlininoACScsv('/Volumes/Extreme SSD/TN444-SOPACE/raw/ACS111/St2_ACS111_20250508_075017.csv');

figure; subplot(2,1,1)
plot(wl_a, st1.a')
subplot(2,1,2)
plot(wl_c, st1.c')

figure; subplot(2,1,1)
plot(wl_a, st1_2.a')
subplot(2,1,2)
plot(wl_c, st1_2.c')

figure; subplot(2,1,1)
plot(wl_a, st2.a')
subplot(2,1,2)
plot(wl_c, st2.c')

st1.a(st1.a < prctile(st1.a, 2.5) | st1.a > prctile(st1.a, 97.5)) = NaN;
st1.c(st1.c < prctile(st1.c, 2.5) | st1.c > prctile(st1.c, 97.5)) = NaN;

st1_2.a(st1_2.a < prctile(st1_2.a, 2.5) | st1_2.a > prctile(st1_2.a, 97.5)) = NaN;
st1_2.c(st1_2.c < prctile(st1_2.c, 2.5) | st1_2.c > prctile(st1_2.c, 97.5)) = NaN;

st2.a(st2.a < prctile(st2.a, 2.5) | st2.a > prctile(st2.a, 97.5)) = NaN;
st2.c(st2.c < prctile(st2.c, 2.5) | st2.c > prctile(st2.c, 97.5)) = NaN;

binned_st1 = table();
binned_st1.dt = mean(st1.dt);
binned_st1.a = mean(st1.a, 'omitmissing');
binned_st1.c = mean(st1.c, 'omitmissing');

binned_st1_2 = table();
binned_st1_2.dt = mean(st1_2.dt);
binned_st1_2.a = mean(st1_2.a, 'omitmissing');
binned_st1_2.c = mean(st1_2.c, 'omitmissing');

binned_st2 = table();
binned_st2.dt = mean(st2.dt);
binned_st2.a = mean(st2.a, 'omitmissing');
binned_st2.c = mean(st2.c, 'omitmissing');

figure; subplot(2,1,1); hold on
plot(wl_a, binned_st1.a', 'linewidth', 1.5)
plot(wl_a, binned_st1_2.a', 'linewidth', 1.5)
plot(wl_a, binned_st2.a', 'linewidth', 1.5)
subplot(2,1,2); hold on
plot(wl_c, binned_st1.c', 'linewidth', 1.5)
plot(wl_c, binned_st1_2.c', 'linewidth', 1.5)
plot(wl_c, binned_st2.c', 'linewidth', 1.5)
legend('St1', 'St1_2', 'St2')
set(gca, 'FontSize', 14)

% Find closest inline qc total
id_st1 = abs(ila.instrument.ACS111.qc.tsw.dt - binned_st1.dt) == min(abs(ila.instrument.ACS111.qc.tsw.dt - binned_st1.dt));
id_st1_2 = abs(ila.instrument.ACS111.qc.tsw.dt - binned_st1_2.dt) == min(abs(ila.instrument.ACS111.qc.tsw.dt - binned_st1_2.dt));
id_st2 = abs(ila.instrument.ACS111.qc.tsw.dt - binned_st2.dt) == min(abs(ila.instrument.ACS111.qc.tsw.dt - binned_st2.dt));

figure; subplot(2,1,1); hold on
plot(wl_a, binned_st1.a, 'linewidth', 1.5)
plot(wl_a, ila.instrument.ACS111.qc.tsw.a(id_st1, :), 'linewidth', 1.5)
plot(wl_a, ila.instrument.ACS111.bin.diw.a(1, :), 'linewidth', 1.5)
subplot(2,1,2); hold on
plot(wl_a, binned_st1.c, 'linewidth', 1.5)
plot(wl_a, ila.instrument.ACS111.qc.tsw.c(id_st1, :), 'linewidth', 1.5)
plot(wl_a, ila.instrument.ACS111.bin.diw.c(1, :), 'linewidth', 1.5)
legend('insitu St1', 'underway St1', 'DIW St1')
set(gca, 'FontSize', 14)

figure; subplot(2,1,1); hold on
plot(wl_a, binned_st1_2.a, 'linewidth', 1.5)
plot(wl_a, ila.instrument.ACS111.qc.tsw.a(id_st1_2, :), 'linewidth', 1.5)
plot(wl_a, ila.instrument.ACS111.bin.diw.a(1, :), 'linewidth', 1.5)
subplot(2,1,2); hold on
plot(wl_a, binned_st1_2.c, 'linewidth', 1.5)
plot(wl_a, ila.instrument.ACS111.qc.tsw.c(id_st1_2, :), 'linewidth', 1.5)
plot(wl_a, ila.instrument.ACS111.bin.diw.c(1, :), 'linewidth', 1.5)
legend('insitu St1_2', 'underway St1_2', 'DIW St1_2')
set(gca, 'FontSize', 14)

figure; subplot(2,1,1); hold on
plot(wl_a, binned_st2.a, 'linewidth', 1.5)
plot(wl_a, ila.instrument.ACS111.qc.tsw.a(id_st2, :), 'linewidth', 1.5)
plot(wl_a, ila.instrument.ACS111.bin.diw.a(3, :), 'linewidth', 1.5)
subplot(2,1,2); hold on
plot(wl_a, binned_st2.c, 'linewidth', 1.5)
plot(wl_a, ila.instrument.ACS111.qc.tsw.c(id_st2, :), 'linewidth', 1.5)
plot(wl_a, ila.instrument.ACS111.bin.diw.c(3, :), 'linewidth', 1.5)
legend('insitu St2', 'underway St2', 'DIW St2')
set(gca, 'FontSize', 14)

% plot raw difference closest underway - insitu
figure; subplot(2,1,1); hold on
plot(wl_a, ila.instrument.ACS111.qc.tsw.a(id_st1, :)-binned_st1.a, 'linewidth', 1.5)
plot(wl_a, ila.instrument.ACS111.qc.tsw.a(id_st1_2, :)-binned_st1_2.a, 'linewidth', 1.5)
plot(wl_a, ila.instrument.ACS111.qc.tsw.a(id_st2, :)-binned_st2.a, 'linewidth', 1.5)
subplot(2,1,2); hold on
plot(wl_c, ila.instrument.ACS111.qc.tsw.c(id_st1, :)-binned_st1.c, 'linewidth', 1.5)
plot(wl_c, ila.instrument.ACS111.qc.tsw.c(id_st1_2, :)-binned_st1_2.c, 'linewidth', 1.5)
plot(wl_c, ila.instrument.ACS111.qc.tsw.c(id_st2, :)-binned_st2.c, 'linewidth', 1.5)
legend('St1', 'St1_2', 'St2')
set(gca, 'FontSize', 14)

% T/S correction based Sullivan et al. 2006
psi = table();
psi.wl = [400;402;404;406;408;410;412;414;416;418;420;422;424;426;428;430;432;434;436;438;440;442;444;446;448;450;452;454;456;458;460;462;464;466;468;470;472;474;476;478;480;482;484;486;488;490;492;494;496;498;500;502;504;506;508;510;512;514;516;518;520;522;524;526;528;530;532;534;536;538;540;542;544;546;548;550;552;554;556;558;560;562;564;566;568;570;572;574;576;578;580;582;584;586;588;590;592;594;596;598;600;602;604;606;608;610;612;614;616;618;620;622;624;626;628;630;632;634;636;638;640;642;644;646;648;650;652;654;656;658;660;662;664;666;668;670;672;674;676;678;680;682;684;686;688;690;692;694;696;698;700;702;704;706;708;710;712;714;716;718;720;722;724;726;728;730;732;734;736;738;740;742;744;746;748;750];
psi.psiT = [0.0001;0.0001;0.0001;0.0001;0;0;0;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0002;0.0003;0.0003;0.0004;0.0005;0.0006;0.0006;0.0007;0.0008;0.0009;0.001;0.001;0.001;0.001;0.001;0.0009;0.0009;0.0008;0.0007;0.0006;0.0006;0.0005;0.0004;0.0003;0.0003;0.0002;0.0001;0.0001;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0002;0.0002;0.0002;0.0001;0.0001;0.0001;0;0;-0.0001;-0.0001;-0.0001;-0.0001;-0.0001;-0.0001;0;0;0.0001;0.0002;0.0003;0.0005;0.0007;0.0009;0.0013;0.0017;0.0021;0.0026;0.0032;0.0038;0.0045;0.0054;0.0063;0.0073;0.0083;0.0094;0.0104;0.0113;0.0121;0.0128;0.0133;0.0136;0.0136;0.0133;0.0129;0.0124;0.0116;0.0107];
psi.a_psiS = [3.0e-05;3.0e-05;3.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;3.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;0;0;0;0;0;0;0;0;0;0;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;0;1.0e-05;2.0e-05;3.0e-05;3.0e-05;4.0e-05;5.0e-05;5.0e-05;6.0e-05;6.0e-05;6.0e-05;6.0e-05;6.0e-05;6.0e-05;6.0e-05;5.0e-05;5.0e-05;5.0e-05;5.0e-05;4.0e-05;4.0e-05;4.0e-05;4.0e-05;3.0e-05;3.0e-05;3.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;1.0e-05;1.0e-05;1.0e-05;0;0;-1.0e-05;-2.0e-05;-3.0e-05;-4.0e-05;-6.0e-05;-7.0e-05;-8.0e-05;-9.0e-05;-0.00011;-0.00012;-0.00014;-0.00015;-0.00016;-0.00017;-0.00018;-0.00019;-0.00020;-0.00020;-0.00020;-0.00021;-0.00021;-0.00021;-0.00021;-0.00021;-0.00020;-0.00017;-0.00013;-8.0e-05;-1.0e-05;7.0e-05;0.00016;0.00026;0.00037;0.00046;0.00054;0.00061;0.00067];
psi.c_psiS = [-1.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-3.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-4.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-5.0e-05;-4.0e-05;-3.0e-05;-3.0e-05;-2.0e-05;-1.0e-05;0;1.0e-05;1.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;2.0e-05;1.0e-05;1.0e-05;1.0e-05;1.0e-05;0;0;0;-1.0e-05;-1.0e-05;-1.0e-05;-1.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-2.0e-05;-3.0e-05;-4.0e-05;-5.0e-05;-6.0e-05;-6.0e-05;-8.0e-05;-9.0e-05;-0.00010;-0.00011;-0.00013;-0.00014;-0.00016;-0.00017;-0.00018;-0.00019;-0.00020;-0.00021;-0.00022;-0.00022;-0.00023;-0.00023;-0.00023;-0.00024;-0.00024;-0.00024;-0.00024;-0.00022;-0.00021;-0.00017;-0.00012;-6.0e-05;2.0e-05;0.00012;0.00022;0.00031;0.00041;0.00049;0.00056;0.00062];
psi.sigma_psiT = [0.0002;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0;0;0;0;0.0001;0.0001;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0;0;0;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0002;0.0002;0.0003;0.0003;0.0004;0.0004;0.0004;0.0004;0.0005;0.0005;0.0006;0.0006;0.0007;0.0007;0.0007;0.0006;0.0005;0.0004;0.0003;0.0003;0.0004;0.0005;0.0006;0.0007;0.0008;0.0009];
psi.c_sigma_psiS = [4e-005;4e-005;4e-005;4e-005;4e-005;4e-005;4e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;0;-1e-005;-2e-005;-3e-005;-4e-005;-6e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;2e-005;2e-005;2e-005;3e-005;3e-005;3e-005];
psi.a_sigma_psiS = [3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;3e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;2e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;NaN;NaN;NaN;NaN;NaN;NaN;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;1e-005;2e-005;2e-005;2e-005;3e-005];

% % try removing psiT and psiS (not working, removing too much)
% a_psiS = interp1(psi.wl, psi.a_psiS, wl_a, 'spline');
% c_psiS = interp1(psi.wl, psi.c_psiS, wl_c, 'spline');
% a_psiT = interp1(psi.wl, psi.psiT, wl_a, 'spline');
% c_psiT = interp1(psi.wl, psi.psiT, wl_c, 'spline');
%
% % find closest T/S data 
% id_st1_tsg = abs(ila.instrument.SBE3845TN444.prod.a.dt - ila.instrument.ACS111.qc.tsw.dt(id_st1)) == min(abs(ila.instrument.SBE3845TN444.prod.a.dt - ila.instrument.ACS111.qc.tsw.dt(id_st1)));
% id_st1_2_tsg = abs(ila.instrument.SBE3845TN444.prod.a.dt - ila.instrument.ACS111.qc.tsw.dt(id_st1_2)) == min(abs(ila.instrument.SBE3845TN444.prod.a.dt - ila.instrument.ACS111.qc.tsw.dt(id_st1_2)));
% id_st2_tsg = abs(ila.instrument.SBE3845TN444.prod.a.dt - ila.instrument.ACS111.qc.tsw.dt(id_st2)) == min(abs(ila.instrument.SBE3845TN444.prod.a.dt - ila.instrument.ACS111.qc.tsw.dt(id_st2)));
%
% % ila.instrument.ACS111.qc.tsw.a(id_st1) - a_psiT * ila.instrument.SBE3845TN444.prod.a.sst(id_st1_tsg) - a_psiT * ila.instrument.SBE3845TN444.prod.a.sst(id_st1_tsg)
% figure; hold on
% plot(wl_a, ila.instrument.ACS111.qc.tsw.a(id_st1,:))
% plot(wl_a, ila.instrument.ACS111.qc.tsw.a(id_st1,:) - a_psiT * ila.instrument.SBE3845TN444.prod.a.sst(id_st1_tsg))
% plot(wl_a, ila.instrument.ACS111.qc.tsw.a(id_st1,:) - a_psiT * ila.instrument.SBE3845TN444.prod.a.sst(id_st1_tsg) - a_psiT * ila.instrument.SBE3845TN444.prod.a.sst(id_st1_tsg))
% legend('raw','T corrected', 'T/S corrected')

% Temperature and salinity dependance correction applied on difference underway - insitu
[a_corr_st1, c_corr_st1] = TemperatureAndSalinityDependence(ila.instrument.ACS111.qc.tsw.a(id_st1, :)-binned_st1.a, ila.instrument.ACS111.qc.tsw.c(id_st1, :)-binned_st1.c, wl_a, wl_c, psi, false(size(ila.instrument.ACS111.qc.tsw.dt)));
[a_corr_st1_2, c_corr_st1_2] = TemperatureAndSalinityDependence(ila.instrument.ACS111.qc.tsw.a(id_st1_2, :)-binned_st1_2.a, ila.instrument.ACS111.qc.tsw.c(id_st1_2, :)-binned_st1_2.c, wl_a, wl_c, psi, false(size(ila.instrument.ACS111.qc.tsw.dt)));
[a_corr_st2, c_corr_st2] = TemperatureAndSalinityDependence(ila.instrument.ACS111.qc.tsw.a(id_st2, :)-binned_st2.a, ila.instrument.ACS111.qc.tsw.c(id_st2, :)-binned_st2.c, wl_a, wl_c, psi, false(size(ila.instrument.ACS111.qc.tsw.dt)));

% plot corrected spectra
figure; subplot(2,1,1); hold on
plot(wl_a, a_corr_st1, 'linewidth', 1.5)
plot(wl_a, a_corr_st1_2, 'linewidth', 1.5)
plot(wl_a, a_corr_st2, 'linewidth', 1.5)
ylabel('Brown stuff absorption [m^{-1}]')
set(gca, 'FontSize', 14)
subplot(2,1,2); hold on
plot(wl_c, c_corr_st1, 'linewidth', 1.5)
plot(wl_c, c_corr_st1_2, 'linewidth', 1.5)
plot(wl_c, c_corr_st2, 'linewidth', 1.5)
ylabel('Brown stuff attenuation [m^{-1}]')
legend('St1', 'St1_2', 'St2')
set(gca, 'FontSize', 14)

prc_rust_signal = mean([a_corr_st1; a_corr_st1_2; a_corr_st2]) ./ ila.instrument.ACS111.prod.p.ap .* 100;
figure; hold on
scatter(ila.instrument.ACS111.prod.p.dt, prc_rust_signal(:,10), 20, 'filled')
scatter(ila.instrument.ACS111.prod.p.dt, prc_rust_signal(:,20), 20, 'filled')
scatter(ila.instrument.ACS111.prod.p.dt, prc_rust_signal(:,30), 20, 'filled')
scatter(ila.instrument.ACS111.prod.p.dt, prc_rust_signal(:,40), 20, 'filled')
ylabel('% rust contribution to UDW a_p')
legend([num2str(round(wl_a(10))) 'nm'],[num2str(round(wl_a(20))) 'nm'],[num2str(round(wl_a(30))) 'nm'],[num2str(round(wl_a(40))) 'nm'])

ila.SpectralQC('AC',{'prod'}, false, {'p','ap'});

ila.instrument.ACS111.prod.p.ap = ila.instrument.ACS111.prod.p.ap - mean([a_corr_st1; a_corr_st1_2; a_corr_st2]);

ila.SpectralQC('AC',{'prod'}, false, {'p','ap'});



%% check and correct rust fouling on bb and c
wl_bb = ila.instrument.HyperBB8005.lambda;
load('/Volumes/Extreme SSD/TN444-SOPACE/DeviceFiles/Hbb_Cal_Plaque_20241010_142447.mat')
load('/Volumes/Extreme SSD/TN444-SOPACE/DeviceFiles/Hbb_Cal_Temp_20241008_130132.mat')
st1 = importInlininoHBB('/Volumes/Extreme SSD/TN444-SOPACE/raw/HyperBB8005/St1_HyperBB8005_20250505_063808.csv', cal, cal_temp);
st1_2 = importInlininoHBB('/Volumes/Extreme SSD/TN444-SOPACE/raw/HyperBB8005/St1_2_HyperBB8005_20250505_065720.csv', cal, cal_temp);
st2 = importInlininoHBB('/Volumes/Extreme SSD/TN444-SOPACE/raw/HyperBB8005/St2_HyperBB8005_20250508_075013.csv', cal, cal_temp);

figure;
plot(wl_bb, st1.beta')

figure;
plot(wl_bb, st1_2.beta')

figure;
plot(wl_bb, st2.beta')

st1.beta(st1.beta < prctile(st1.beta, 2.5) | st1.beta > prctile(st1.beta, 97.5)) = NaN;
st1_2.beta(st1_2.beta < prctile(st1_2.beta, 2.5) | st1_2.beta > prctile(st1_2.beta, 97.5)) = NaN;
st2.beta(st2.beta < prctile(st2.beta, 2.5) | st2.beta > prctile(st2.beta, 97.5)) = NaN;

binned_st1 = table();
binned_st1.dt = mean(st1.dt);
binned_st1.beta = mean(st1.beta, 'omitmissing');

binned_st1_2 = table();
binned_st1_2.dt = mean(st1_2.dt);
binned_st1_2.beta = mean(st1_2.beta, 'omitmissing');

binned_st2 = table();
binned_st2.dt = mean(st2.dt);
binned_st2.beta = mean(st2.beta, 'omitmissing');

figure; hold on
plot(wl_bb, binned_st1.beta', 'linewidth', 1.5)
plot(wl_bb, binned_st1_2.beta', 'linewidth', 1.5)
plot(wl_bb, binned_st2.beta', 'linewidth', 1.5)
legend('St1', 'St1_2', 'St2')
set(gca, 'FontSize', 14)

% Find closest inline qc total
id_st1 = abs(ila.instrument.HyperBB8005.qc.tsw.dt - binned_st1.dt) == min(abs(ila.instrument.HyperBB8005.qc.tsw.dt - binned_st1.dt));
id_st1_2 = abs(ila.instrument.HyperBB8005.qc.tsw.dt - binned_st1_2.dt) == min(abs(ila.instrument.HyperBB8005.qc.tsw.dt - binned_st1_2.dt));
id_st2 = abs(ila.instrument.HyperBB8005.qc.tsw.dt - binned_st2.dt) == min(abs(ila.instrument.HyperBB8005.qc.tsw.dt - binned_st2.dt));

figure; hold on
plot(wl_bb, binned_st1.beta, 'linewidth', 1.5)
plot(wl_bb, ila.instrument.HyperBB8005.qc.tsw.beta(id_st1, :), 'linewidth', 1.5)
legend('insitu St1', 'underway St1')
set(gca, 'FontSize', 14)

figure; hold on
plot(wl_bb, binned_st1_2.beta, 'linewidth', 1.5)
plot(wl_bb, ila.instrument.HyperBB8005.qc.tsw.beta(id_st1_2, :), 'linewidth', 1.5)
legend('insitu St1_2', 'underway St1_2')
set(gca, 'FontSize', 14)

figure; hold on
plot(wl_bb, binned_st2.beta, 'linewidth', 1.5)
plot(wl_bb, ila.instrument.HyperBB8005.qc.tsw.beta(id_st2, :), 'linewidth', 1.5)
legend('insitu St2', 'underway St2')
set(gca, 'FontSize', 14)



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



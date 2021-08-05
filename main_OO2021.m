% Main Particulate InLine Analysis Script
% author: Guillaume Bourdin
% created: May 05, 2021
clear
if feature('IsDebugMode'); dbquit all; end

cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master/')

% Load InLineAnalysis and the configuration
ila = InLineAnalysis('cfg/OO2021_cfg.m');

% Quick cfg update
%% date to process
ila.cfg.days2run = datenum(2021,8,4,0,0,0):datenum(2021,8,4,0,0,0);

%%
ila.cfg.instruments2run = {'FLOW','ACS94'}; % 'FLOW', 'TSG', 'BB3', 'ACS94', 'AC9274'
ila.cfg.qcref.view = 'ACS94';
ila.cfg.parallel = Inf;
ila.cfg.calibrate.(ila.cfg.qcref.view).compute_dissolved = false;

%% 1. Import | Load raw data
ila.cfg.force_import = false;
ila.ReadRaw();
ila.CheckDataStatus();

%% 2. Synchronise instruments
% % Independent of flow rate (for now)
% % If flow rate varies use the Strech method
% % Play with delay of synchronisation
% % TSG is assumed to be set at zero
% % ila.instrument.FLOW.Sync(30);
ila.instrument.TSG.Sync(50);
ila.instrument.ACS94.Sync(0);
ila.instrument.AC9274.Sync(0);
ila.instrument.HBB.Sync(0);
ila.instrument.BB3.Sync(0);
% % Quick visualizzation to sync with TSG
% fig(30, 'sync TSG');
% yyaxis('left'); plot(ila.instrument.TSG.data.dt, ila.instrument.TSG.data.t); ylabel('Temperature (^o C)');
% % yyaxis('right'); plot(ila.instrument.BB3.data.dt, ila.instrument.BB3.data.beta(:,2)); ylabel('\beta (m^{-1} sr^{-1})'); ylim([80 300]);
% % datetick2_doy();
% visSync(ila.instrument.BB3.data, ila.instrument.TSG.data.dt, ila.instrument.TSG.data.t, 'Temp (C)');
visSync(ila.instrument.FLOW.data, ila.instrument.ACS94.data.dt, ila.instrument.ACS94.data.a(:,20), 'a (m^{-1})');
visSync(ila.instrument.FLOW.data, ila.instrument.ACS94.data.dt, ila.instrument.ACS94.data.c(:,40), 'c (m^{-1})');
visSync(ila.instrument.FLOW.data, ila.instrument.AC9274.data.dt, ila.instrument.AC9274.data.a(:,5), 'a (m^{-1})');
visSync(ila.instrument.FLOW.data, ila.instrument.AC9274.data.dt, ila.instrument.AC9274.data.c(:,5), 'c (m^{-1})');
visSync(ila.instrument.FLOW.data, ila.instrument.BB3.data.dt, ila.instrument.BB3.data.beta(:,1), '\beta (counts)');
% % % Once settings are good set them in the configuration file.
% % % The software is now doing the same with one line of code.
% ila.Sync()

%% 2. Auto-synchronise: automatic detection of filter events for AC and BB sensors
% ila.cfg.qcref.MinFiltPeriod = 65; % filter even period in minute % ACS: 55 % BB3: 60
% ila.cfg.qcref.szFilt = 12; % filter even length in minute % default = 10
% ila.SplitDetect(ila.cfg.qcref.MinFiltPeriod, ila.cfg.qcref.szFilt);

%% 3. QC Reference
% run with mode ui during first run (it saves your work for the next run)
% run with mode load to load previous QC
% Note: when redoing QC of a given period of time (days2run) the previous
% QC during the same period of time is erased, QC done on other periods of
% time is kept in the json file
ila.cfg.qcref.mode='ui'; % 'ui' or 'load'
ila.cfg.qcref.remove_old = false; % remove old selection of the same period
ila.QCRef

%% 4. Split fsw and tsw
ila.cfg.split.skip = {'FLOW','TSG'};
ila.Split();
ila.CheckDataStatus();

%% 4.1. Diagnostic Plot
% check raw spectrums AC or BB sensors
ila.DiagnosticPlot('AC',{'raw'}); % AC or BB and level to plot

%% 5. Automatic QC of raw data for step in ACS spectrum, BB saturated and obvious bad PAR values
% fudge factor for auto QC ACS.
% Varies between ACS: 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.RawAutoQCLim.filtered.a = 3;
ila.cfg.qc.RawAutoQCLim.filtered.c = 3;
ila.cfg.qc.RawAutoQCLim.total.a = 3;
ila.cfg.qc.RawAutoQCLim.total.c = 3;
% fudge factor for auto QC BB.
% 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.RawAutoQCLim.filtered.bb = 3;
ila.cfg.qc.RawAutoQCLim.total.bb = 3;
% remove saturated periods in BB
ila.cfg.qc.Saturation_Threshold_bb = 4100; % saturate above 4000 counts
ila.RawAutoQC('raw');
ila.CheckDataStatus();

%% 5.1. Diagnostic Plot
% check raw spectrums AC or BB sensors
ila.DiagnosticPlot('AC',{'raw'}); % AC or BB and level to plot

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
ila.cfg.qc.mode='ui';  % load or ui
ila.cfg.qc.specific.run = {ila.cfg.qcref.view};
ila.QC();

%% 5.4. Write raw (keep raw during processing, delete once done)
ila.Write('raw', 'part')
ila.CheckDataStatus();

%% 6. Bin
% % Set settings directly in configuration file (no tunning at this step)
ila.cfg.bin.skip = {};
ila.Bin()
ila.CheckDataStatus();

%% 6.1. Diagnostic Plot
% check binned spectrums AC or BB sensors
ila.DiagnosticPlot('AC',{'bin'}); % AC or BB and level to plot

%% 6.2. Write bin
ila.Write('bin', 'part')
ila.CheckDataStatus();

%% Load processed data from mat files
% ila.Read('raw');
% ila.Read('bin');
% ila.Read('qc');
% ila.Read('prod');

%% 7. Flag
ila.Flag() % Now deprecated will juqst copy data to next level
ila.CheckDataStatus();

%% 8. QC Interactive or Loading previous qc selection
%%%%% Settings %%%%%
ila.cfg.qc.mode='ui';  % load or ui
ila.cfg.qc.remove_old = false;  % remove old selection of this period
ila.cfg.qc.qc_once_for_all = false;  % true = QC all variables | false = QC variables separately)
% Global
ila.cfg.qc.global.view = {ila.cfg.qcref.view};
ila.cfg.qc.global.active = false;
% Specific
ila.cfg.qc.specific.run = {ila.cfg.qcref.view};
%%%%%%%%%%%%%%%%%%%

% Run QC function
ila.QC();
ila.CheckDataStatus();

%% 8.1. Optional auto QC on 'qc' level: run until it stabilize to 0
ila.RawAutoQC('qc');

%% 8.2. Diagnostic Plot
% check QCed spectrums AC or BB sensors
% {'raw','bin','qc','prod'}
ila.DiagnosticPlot('AC',{'qc'}); % AC or BB and level to plot

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
ila.DiagnosticPlot('AC',{'qc'}, false, {'tsw','c'});
ila.DiagnosticPlot('BB',{'qc'}, false, {'fsw','all'});

%% 9. QC Switch position
% QC switch position to make sure each filter event is separated by a
% period of total water and eventually move filter events
ila.QCSwitchPosition()

%% 9.1. Write qc
ila.Write('qc', 'part')
ila.CheckDataStatus();

%% 10. Calibrate and compute products
ila.cfg.calibrate.skip = {'FLOW', 'TSG'};
ila.Calibrate();
ila.CheckDataStatus();

%% 10.1 Product visualisation plots with option to save
save_figures = true;

%%% AC or BB 3D plots %%%
ila.DiagnosticPlot('AC', {'prod'}, save_figures); % AC or BB

%%% ACS BB3 TSG PAR WSCD SPCD ALFA LISST final product visualisation %%%
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
ila.DiagnosticPlot('BB',{'prod'}, false, {'p','all'});

%% 11.1. Load previous qc pick selection at prod level
ila.cfg.qc.mode='load';  % load or ui
ila.cfg.qc.specific.run = {ila.cfg.qcref.view}; % 'FLOW','ACS57','TSG', 'BB31502', 'WSCD859','PAR'
ila.QC();

%% 12. Save products
ila.Write('prod', 'part')

% % Notify with a song that the job is done
% notif_sound = load('gong'); sound(notif_sound.y, notif_sound.Fs); % handel
% return

%% re-write final version of 'qc' and 'bin'
ila.Write('raw', 'part')
ila.Write('bin', 'part')
ila.Write('qc', 'part')













%% 9.1 Save as one file
% tsg = ila.instrument.TSG.prod.a;
% % save([ila.instrument.TSG.path.prod '../TaraPacific_TSG'], 'tsg');
% tsg.dt = datestr(tsg.dt, 'yyyy/mm/dd HH:MM:SS');
% writetable(tsg, [ila.instrument.TSG.path.prod '../TaraPacific_TSG.csv']);

%% 10. Figures !!
%!!! SAVE PREVIOUS WORK BEFORE RUNNING THOSE LINE !!!
% ila = InLineAnalysis('cfg/EXPORTS_cfg.m');
% ila.cfg.instruments2run = {'TSG', 'ACS', 'BB3', 'WSCD'};
% ila.cfg.days2run = datenum(2018,08,15):datenum(2018,09,09);
% ila.Load('prod');
%!!! SAVE PREVIOUS WORK BEFORE RUNNING THOSE LINE !!!

% save_figures = true;

% if save_figures; figure(72); savefig([ila.instrument.ACS.path.prod 'ACS' datestr(ila.cfg.days2run(1), 'yyyymmdd') '_a_p']); end
% if save_figures; figure(73); savefig([ila.instrument.ACS.path.prod 'ACS' datestr(ila.cfg.days2run(1), 'yyyymmdd') '_c_p']); end
% save_fig([ila.instrument.ACS.path.prod 'ACS_stnP_ap'], 1280, 720);

% Display 676 line on 3D plot
% hold('on'); fill3([676 676 676 676], [736961 736937 736937 736961], [0.14 0.14 -0.02 -0.02], [1 1 1 1])

% Check ACS wl registration
% fig(78); hold('on'); plot(wl, ACS111.p.ap(2000:2010, :), 'o-'); plot([676 676], ylim(), 'k'); plot([675 675], ylim(), 'k');

% %% LISST particulate %%%
% % ila.Calibrate();
% 
% % Plot VSF
% sel = ~any(ila.instrument.LISST.prod.p.betap <= 0,2);
% % 3D
% visProd3D(ila.instrument.LISST.theta, ila.instrument.LISST.prod.p.dt(sel), ila.instrument.LISST.prod.p.betap(sel,:), false, 'Log', false, 81);
% set(gca, 'ZScale', 'log', 'XScale', 'log');
% xlabel('\theta (^o)'); zlabel('VSF (m^{-1})');
% % 2D
% visProd2D(ila.instrument.LISST.theta, ila.instrument.LISST.prod.p.dt(sel), ila.instrument.LISST.prod.p.betap(sel,:)); set(gca, 'YScale', 'log');


% Plot PSD (2D)
% sel = 1 <= ila.instrument.LISST.diameters & ila.instrument.LISST.diameters <= 100;
% % sel_bin = ~(any(ila.instrument.LISST.prod.p.PSD < 0,2)); % any(~isnan(ila.instrument.LISST.prod.p.PSD),2);
% visProd2D(ila.instrument.LISST.diameters(sel), ila.instrument.LISST.prod.p.dt, ila.instrument.LISST.prod.p.PSD(:,sel), false);
% set(gca, 'XScale', 'log', 'YScale', 'log');
% xlabel('Diameter (\mum)'); ylabel('PSD (# mL^{-1} \mum^{-1})'); 
% if save_figures; savefig([ila.instrument.LISST.path.prod 'LISST' datestr(ila.cfg.days2run(1), 'yyyymmdd') '_PSD_2D']); end
% save_fig([ila.instrument.LISST.path.prod 'LISST_stnP_PSD_2D'], 1280, 720);
% visProd3D(ila.instrument.LISST.diameters(sel), ila.instrument.LISST.prod.p.dt, ila.instrument.LISST.prod.p.VSD(:,sel), false, 'Log', false, 78);
% set(gca, 'ZScale', 'log', 'XScale', 'log'); %zlim([10^-4 10^4]);
% xlabel('Diameter (\mum)'); zlabel('PSD (# mL^{-1} \mum^{-1})');

% Plot VSD (3D)
% sel = 1 <= ila.instrument.LISST.diameters & ila.instrument.LISST.diameters <= 64;
% visProd3D(ila.instrument.LISST.diameters(sel), ila.instrument.LISST.prod.p.dt, ila.instrument.LISST.prod.p.VSD(:,sel) * 10^-6, false, 'Log', false, 82);
% set(gca, 'ZScale', 'log', 'XScale', 'log'); zlim([10^-4 0.1]);
% xlabel('Diameter (\mum)'); zlabel('VSD (ppm \mum^{-1})');
% if save_figures; savefig([ila.instrument.LISST.path.prod 'LISST' datestr(ila.cfg.days2run(1), 'yyyymmdd') '_PSD_3D']); end
% set(gca, 'FontSize', 14, 'FontName', 'Helvetica Neue');
% save_fig([ila.instrument.LISST.path.prod 'LISST_stnP_VSD'], 1280, 1024);



% %% ALFA %%%
% fig(64); ALFA = ila.instrument.ALFA.prod.a;
% yyaxis('right'); plot(ALFA.dt, ALFA.Chlb, '.-'); ylabel('Chl_b Fl (\mug L^{-1})');
% yyaxis('left'); hold('on'); ylabel('F_{v}/F_{m}');
% plot(ALFA.dt(ALFA.FvFm>0.2), ALFA.FvFm(ALFA.FvFm>0.2), '.-');
% plot(ALFA.dt(ALFA.FvFmG>0.2), ALFA.FvFmG(ALFA.FvFmG>0.2), 'g.-'); 
% datetick2_doy(); %xlabel('Aug 15-17, 2018');
% set(datacursormode(figure(64)),'UpdateFcn',@data_cursor_display_date);
% % if save_figures; savefig([ila.instrument.ALFA.path.prod 'ALFA' datestr(ila.cfg.days2run(1), 'yyyymmdd') '_FvFm']); end
% % save_fig([ila.instrument.ALFA.path.prod 'ALFA_stnP_FvFm'], 1280, 720);

% %% TSG %%%
% fig(65); TSG = ila.instrument.TSG.prod.a;
% yyaxis('right'); plot(TSG.dt, TSG.par, '.-'); ylabel('PAR (\muE/s/m^2)');
% yyaxis('left'); plot(TSG.dt, TSG.t, '.-'); ylabel('Temperature (^oC)');
% datetick2_doy(); set(datacursormode(figure(65)),'UpdateFcn',@data_cursor_display_date);
% % xlabel('Aug 15-17, 2018');
% % set(gca, 'FontSize', 14, 'FontName', 'Helvetica Neue');
% % save_fig([ila.instrument.TSG.path.prod 'Underway_stnP_PAR_T'], 1280, 600);
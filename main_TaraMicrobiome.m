% Main Particulate InLine Analysis Script
% author: Guillaume Bourdin
% created: Jan 05, 2021
clear
close all
cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master/')

% Load InLineAnalysis and the configuration
ila = InLineAnalysis('cfg/TaraMicrobiome_cfg.m');

% Quick cfg update
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
ila.cfg.days2run = datenum(2021,4,7):datenum(2021,4,14);
% ila.cfg.days2run = datenum(2021,4,15):datenum(2021,4,22);
% ila.cfg.days2run = datenum(2021,4,23):datenum(2021,4,30);
% ila.cfg.days2run = datenum(2021,5,1):datenum(2021,5,9);

%% BB31502
% ila.cfg.days2run = datenum(2020,12,26):datenum(2021,1,5);
% ila.cfg.days2run = datenum(2021,1,6):datenum(2021,1,20);
% ila.cfg.days2run = datenum(2021,1,21):datenum(2021,2,5);
% ila.cfg.days2run = datenum(2021,2,17):datenum(2021,3,5);
% ila.cfg.days2run = datenum(2021,3,6):datenum(2021,3,26);
% ila.cfg.days2run = datenum(2021,3,26):datenum(2021,4,14);
% ila.cfg.days2run = datenum(2021,4,15):datenum(2021,5,9);

%% WSCD859
% ila.cfg.days2run = datenum(2020,12,12):datenum(2021,5,9);

%% SPCD
%  ila.cfg.days2run = datenum(2021,5,1):datenum(2021,5,7);
 
%% LISST
%  ila.cfg.days2run = datenum(2021,5,1):datenum(2021,5,7);
 
%% HBB
%  ila.cfg.days2run = datenum(2021,5,1):datenum(2021,5,7);

%% ALPHA
ila.cfg.days2run = datenum(2020,12,12):datenum(2021,5,9);

%%
ila.cfg.instruments2run = {'FLOW','ACS57'}; % 'ALPHA', 'FLOW','ACS57','TSG', 'BB31502', 'WSCD859','PAR','SPCD','ACS91','LISST1183'
ila.cfg.qcref.view = 'ACS57';
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
% % No noticeable difference was observed between the TSG of EXPORTS and the BB3
% % ila.instrument.FLOW.Sync(30);
ila.instrument.TSG.Sync(0);
ila.instrument.SPCD.Sync(0);
ila.instrument.ACS57.Sync(0);
ila.instrument.HBB.Sync(0);
ila.instrument.BB31502.Sync(0);
ila.instrument.LISST1183.Sync(0);
ila.instrument.WSCD859.Sync(0);
ila.instrument.ALFA.Sync(0); 
% % Quick visualizzation to sync with TSG
% fig(30, 'sync TSG');
% yyaxis('left'); plot(ila.instrument.TSG.data.dt, ila.instrument.TSG.data.t); ylabel('Temperature (^o C)');
% % yyaxis('right'); plot(ila.instrument.BB3.data.dt, ila.instrument.BB3.data.beta(:,2)); ylabel('\beta (m^{-1} sr^{-1})'); ylim([80 300]);
% % datetick2_doy();
% visSync(ila.instrument.BB3.data, ila.instrument.TSG.data.dt, ila.instrument.TSG.data.t, 'Temp (C)');
visSync(ila.instrument.FLOW.data, ila.instrument.SPCD.data.dt, ila.instrument.SPCD.data.fdom, 'FDOM (counts)');
visSync(ila.instrument.FLOW.data, ila.instrument.ACS57.data.dt, ila.instrument.ACS57.data.a(:,20), 'a (m^{-1})');
visSync(ila.instrument.FLOW.data, ila.instrument.ACS57.data.dt, ila.instrument.ACS57.data.c(:,40), 'c (m^{-1})');
visSync(ila.instrument.FLOW.data, ila.instrument.HBB.data.dt, ila.instrument.HBB.data.beta(:,14), '\beta (counts)');
visSync(ila.instrument.FLOW.data, ila.instrument.BB31502.data.dt, ila.instrument.BB31502.data.beta(:,1), '\beta (counts)');
visSync(ila.instrument.FLOW.data, ila.instrument.LISST1183.data.dt, ila.instrument.LISST1183.data.beta(:,10), '\beta (counts)');
visSync(ila.instrument.FLOW.data, ila.instrument.WSCD859.data.dt, ila.instrument.WSCD859.data.fdom, 'FDOM (counts)');
visSync(ila.instrument.FLOW.data, ila.instrument.ALFA.data.dt, ila.instrument.ALFA.data.Chlb, 'chlb');yyaxis('left'); ylim([0 2]);
% % 
% % % xlim([datenum(2018,08,14,9,55,0) datenum(2018,08,14,11,05,0)]);
% % % ylim([-0.1 0.2]);
% % % Once settings are good set them in the configuration file.
% % % The software is now doing the same with one line of code.
% ila.Sync()
% % % ila.instrument.BB31502.Sync(-90);
% % % ila.instrument.BB31502.Sync(-10);

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
ila.QCRef();

%% 4. Split fsw and tsw
ila.cfg.split.skip = {'FLOW','TSG','PAR','ALFA'};
ila.Split();
ila.CheckDataStatus();

%% 4.1. Diagnostic Plot
% check raw spectrums AC or BB sensors
ila.DiagnosticPlot('AC',{'raw'}); % AC or BB

%% 5. Automatic QC of raw data for step in ACS spectrum, BB saturated and obvious bad PAR & ALFA values
% fudge factor for auto QC ACS.
% Varies between ACS: 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.RawAutoQCLim.filtered.a = 6;
ila.cfg.qc.RawAutoQCLim.filtered.c = 19;
ila.cfg.qc.RawAutoQCLim.total.a = 4;
ila.cfg.qc.RawAutoQCLim.total.c = 12;
% fudge factor for auto QC BB.
% 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.RawAutoQCLim.filtered.bb = 7;
ila.cfg.qc.RawAutoQCLim.total.bb = 5;
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
ila.DiagnosticPlot('AC',{'raw'}, false, {'tsw','c'});

%% 5.3. Loading previous qc pick selection at raw level
ila.cfg.qc.mode='load';  % load or ui
ila.cfg.qc.specific.run = {ila.cfg.qcref.view}; % 'FLOW','ACS57','TSG', 'BB31502', 'WSCD859','PAR'
ila.QC();

%% 5.4. Write clean raw after split for BB3 and HBB
ila.Write('raw')
ila.CheckDataStatus();

%% 7. Bin
% % Set settings directly in configuration file (no tunning at this step)
% % run before re-bin only to clear qc tables
% ila.instrument.ACS57.qc.tsw = table(); ila.instrument.ACS57.qc.fsw = table();
ila.cfg.bin.skip = {};
ila.Bin()
ila.CheckDataStatus();

%% 6.1. Diagnostic Plot
% check binned spectrums AC or BB sensors
ila.DiagnosticPlot('AC',{'bin'}); % AC or BB

%% 6.2. Write bin
ila.Write('bin')
ila.CheckDataStatus();

%% Load data from mat files
% ila.Read('raw');
% ila.Read('bin');
% ila.Read('qc');
% ila.Read('prod');

%% 7. Flag
ila.Flag() % Now deprecated will just copy data to next level
ila.CheckDataStatus();

%% 8. QC Interactive or Loading previous qc selection
%%%%% Settings %%%%%
ila.cfg.qc.mode='ui';  % load or ui
ila.cfg.qc.qc_once_for_all = true;  % true = QC 'a' and 'c' together | false = QC 'a' and 'c' separately
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
% {'raw','bin','qc','prod'}
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
ila.DiagnosticPlot('AC',{'qc'}, false, {'fsw','c'});
ila.DiagnosticPlot('AC',{'qc'}, false, {'tsw','all'});

%% 9. QC Switch position
% QC switch position to make sure each filter event is separated by a
% period of total water and eventually move filter events
ila.QCSwitchPosition()

%% 9.1. Write qc
ila.Write('qc')
ila.CheckDataStatus();

%% 10. Calibrate
ila.cfg.calibrate.skip = {'FLOW', 'TSG'};
ila.Calibrate();
ila.CheckDataStatus();
% % Calibrate LISST only
% % ila.instrument.LISST.inversion = 'non-spherical';
% % ila.instrument.LISST.Calibrate()
% % Calibrate ACS only
% ila.instrument.ACS111.Calibrate(ila.cfg.calibrate.ACS111.compute_dissolved,...
%                              ila.cfg.calibrate.ACS111.interpolation_method,...
%                              ila.instrument.(ila.cfg.calibrate.ACS111.CDOM_source),...
%                              ila.instrument.(ila.cfg.calibrate.ACS111.FLOW_source))
% 
% % Compute NAAMES specific chl
% wl = ila.instrument.ACS111.lambda_ref; ACS = ila.instrument.ACS111.prod;
% % wl = ila.instrument.ACS298.lambda_ref; ACS = ila.instrument.ACS298.prod;
% ap_a=interp1(wl,ACS111.p.ap',[650 676 715],'linear')';
% line_height = (ap_a(:,2)-(39/65*ap_a(:,1)+26/65*ap_a(:,3)));
% p.chl=157*line_height.^1.22;
% ila.instrument.ACS111.prod.p.chl_naames = 95 * line_height .^ 1.06;
% ila.instrument.ACS111.prod.p.chl_naames(real(ila.instrument.ACS111.prod.p.chl_naames) ~= ila.instrument.ACS.prod.p.chl_naames) = NaN;
% fprintf('Done\n');
% 
% % Compute EXPORTS Specific chl
% % Derive Chl (Line heigh at 676 compared to 650 and 715)
% % ila_acs = ila.instrument.ACS298.prod.p; ila_acs_wl = ila.instrument.ACS298.lambda_ref;
% ila_acs = ila.instrument.ACS111.prod.p; ila_acs_wl = ila.instrument.ACS111.lambda_ref;
% ap = interp1(ila_acs_wl, ila_acs.ap',[650 676 715],'linear')';
% line_height = (ap(:,2)-(39/65*ap(:,1)+26/65*ap(:,3)));
% % ila.instrument.ACS298.prod.p.chl_exports = 138.14 * line_height.^1.11; % EXPORTS relation
% ila.instrument.ACS111.prod.p.chl_exports = 138.14 * line_height.^1.11; % EXPORTS relation

%% 10.1 Product visualisation plots with option to save
save_figures = false;

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
ila.DiagnosticPlot('AC',{'prod'}, false, {'p','ap'});

%% 11.1. Load previous qc pick selection at prod level
ila.cfg.qc.mode='load';  % load or ui
ila.cfg.qc.specific.run = {ila.cfg.qcref.view}; % 'FLOW','ACS57','TSG', 'BB31502', 'WSCD859','PAR'
ila.QC();

%% 12. Save products
ila.Write('prod')

% % Notify with a song that the job is done
% notif_sound = load('gong'); sound(notif_sound.y, notif_sound.Fs); % handel
% return

%% re-write final version of 'qc' and 'bin'
ila.Write('bin')
ila.Write('qc')









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
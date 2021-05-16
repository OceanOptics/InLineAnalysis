% Main InLine Analysis Script
% author: Guillaume Bourdin
% created: Jan 05, 2021
clear
close all
cd('/Users/emmanuel.boss/Desktop/InLine analysis/InLineAnalysis-master/')

% Load InLineAnalysis and the configuration
ila = InLineAnalysis('cfg/EXPORTS02_cfg.m');

% Quick cfg update
%% TSG
 ila.instrument.TSG.logger = 'matlab_Emmanuel'; % SBE45TSG matlab_Emmanuel
 ila.cfg.days2run = datenum(2021,5,1,0,0,0):datenum(2021,5,10,0,0,0);
% ila.instrument.TSG.logger = 'TeraTerm';
% ila.cfg.days2run = datenum(2021,1,10,0,0,0):datenum(2021,1,15,0,0,0);

%% ACS091
ila.cfg.days2run = datenum(2021,5,2,0,0,0):datenum(2021,5,5,0,0,0);
% ila.cfg.days2run = datenum(2021,1,6,0,0,0):datenum(2021,1,20,0,0,0);
% ila.cfg.days2run = datenum(2021,1,20,0,0,0):datenum(2021,2,5,0,0,0);

%% BB31502
 ila.cfg.days2run =datenum(2021,5,2,0,0,0):datenum(2021,5,5,0,0,0);
% ila.cfg.days2run = datenum(2021,1,6,0,0,0):datenum(2021,1,20,0,0,0);
% ila.cfg.days2run = datenum(2021,1,20,0,0,0):datenum(2021,2,5,0,0,0);

%% SeapointCD
 ila.cfg.days2run = datenum(2021,1,2,0,0,0):datenum(2021,5,5,0,0,0);
 
%% HBB
 ila.cfg.days2run = datenum(2021,5,5,0,0,0):datenum(2021,5,10,0,0,0);

%%
ila.cfg.instruments2run = {'FLOW','TSG','HBB'}; % 'FLOW', 'TSG', 'BB31052', 'HBB', 'WSCD859','SPCD'
ila.cfg.qcref.view = 'HBB';
ila.cfg.parallel = Inf;

%% 1.Normal Import | Load raw data
ila.cfg.force_import = false;
ila.ReadRaw();
ila.CheckDataStatus();

%% 1.DI. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read DI 
% To get ag and cg from ACS can run DI day by day
% To get betag from BB3 need to run the entire dataset if select di method constant
% ila.cfg.days2run = ila.cfg.days2run(1)-1:ila.cfg.days2run(end)+1;
ila.cfg.force_import = false;
ila.ReadRawDI();

%% 2.Normal and DI Synchronise instruments
% % Independent of flow rate (for now)
% % If flow rate varies use the Strech method
% % Play with delay of synchronisation
% % TSG is assumed to be set at zero
% % No noticeable difference was observed between the TSG of EXPORTS and the BB3
% % ila.instrument.FLOW.Sync(30);
% ila.instrument.ACS57.Sync(10);
ila.instrument.ACS91.Sync(0);
% % ila.instrument.ACS111.Sync(60);
% % ila.instrument.ACS279.Sync(55);
ila.instrument.HBB.Sync(0);
% ila.instrument.BB3.Sync(0);
% % % ila.instrument.LISST.Sync(1);
% % ila.instrument.WSCD1082P.Sync(40);
% ila.instrument.WSCD859.Sync(40);
% % ila.instrument.TSG.Sync(0);
% % ila.instrument.ALFA.Sync(15); 
% % Quick visualizzation to sync with TSG
% % fig(30, 'sync TSG');
% % yyaxis('left'); plot(ila.instrument.TSG.data.dt, ila.instrument.TSG.data.t); ylabel('Temperature (^o C)');
% % yyaxis('right'); plot(ila.instrument.BB3.data.dt, ila.instrument.BB3.data.beta(:,2)); ylabel('\beta (m^{-1} sr^{-1})'); ylim([80 300]);
% % datetick2_doy();
% visSync(ila.instrument.FLOW.data, ila.instrument.ACS007.data.dt, ila.instrument.ACS007.data.a(:,20), 'a (m^{-1})');
% visSync(ila.instrument.FLOW.data, ila.instrument.ACS007.data.dt, ila.instrument.ACS007.data.c(:,40), 'c (m^{-1})');
visSync(ila.instrument.FLOW.data, ila.instrument.ACS91.data.dt, ila.instrument.ACS91.data.a(:,20), 'a (m^{-1})');
visSync(ila.instrument.FLOW.data, ila.instrument.ACS91.data.dt, ila.instrument.ACS91.data.c(:,40), 'c (m^{-1})');
% % visSync(ila.instrument.FLOW.data, ila.instrument.ACS111.data.dt, ila.instrument.ACS111.data.a(:,20), 'a (m^{-1})');
% % visSync(ila.instrument.FLOW.data, ila.instrument.ACS111.data.dt, ila.instrument.ACS111.data.c(:,40), 'c (m^{-1})');
% % visSync(ila.instrument.FLOW.data, ila.instrument.ACS279.data.dt, ila.instrument.ACS279.data.a(:,20), 'a (m^{-1})');
% % visSync(ila.instrument.FLOW.data, ila.instrument.ACS279.data.dt, ila.instrument.ACS279.data.c(:,40), 'c (m^{-1})');
% % visSync(ila.instrument.FLOW.data, ila.instrument.BB3.data.dt, ila.instrument.BB3.data.beta(:,1), '\beta (counts)');
visSync(ila.instrument.FLOW.data, ila.instrument.HBB.data.dt, ila.instrument.HBB.data.beta(:,14), '\beta (counts)');
% % visSync(ila.instrument.FLOW.data, ila.instrument.LISST.data.dt, ila.instrument.LISST.data.beta(:,10), '\beta (counts)');
% % visSync(ila.instrument.('FLOW').data, ila.instrument.('WSCD').data.dt, ila.instrument.('WSCD').data.fdom, 'FDOM (counts)');
% % visSync(ila.instrument.('FLOW').data, ila.instrument.('WSCD1082P').data.dt, ila.instrument.('WSCD1082P').data.fdom, 'FDOM (counts)');
% % visSync(ila.instrument.('BB3').data, ila.instrument.('TSG').data.dt, ila.instrument.('TSG').data.t, 'Temp (C)');
% % % visSync(ila.instrument.FLOW.data, ila.instrument.ALFA.data.dt, ila.instrument.ALFA.data.Chlb, 'chlb');yyaxis('left'); ylim([0 2]);
% % 
% % % xlim([datenum(2018,08,14,9,55,0) datenum(2018,08,14,11,05,0)]);
% % % ylim([-0.1 0.2]);
% % % Once settings are good set them in the configuration file.
% % % The software is now doing the same with one line of code.
% ila.Sync()
% % % ila.instrument.BB31502.Sync(-90);
% % % ila.instrument.BB31502.Sync(-10);

%% 2.Normal Auto-synchronise: automatic detection of filter events for AC and BB sensors
% % % data = ila.instrument.(ila.cfg.qcref.view).data;
% % % instrument = ila.cfg.qcref.view;
% % % FLOW = ila.instrument.FLOW.data;

ila.cfg.qcref.MinFiltPeriod = 65; % filter even period in minute % ACS: 55 % BB3: 60
ila.cfg.qcref.szFilt = 12; % filter even length in minute % default = 10
ila.SplitDetect(ila.cfg.qcref.MinFiltPeriod, ila.cfg.qcref.szFilt);

%% 3.Normal QC Reference
% run with mode ui during first run (it saves your work for the next run)
% run with mode load to load previous QC
% Note: when redoing QC of a given period of time (days2run) the previous
% QC during the same period of time is erased, QC done on other periods of
% time is kept in the json file
ila.cfg.qcref.mode='ui'; % 'ui' or 'load'
ila.QCRef();

%% 4.Normal Split fsw and tsw
ila.cfg.split.skip = {'FLOW','TSG'};
ila.Split();
ila.CheckDataStatus();

%% Normal and DI Diagnostic Plot
% check raw spectrums AC or BB sensors
ila.DiagnosticPlot('BB',{'raw'}); % AC or BB

%% 5.Normal and DI automatic QC of raw data for step in ACS spectrum, BB saturated and obvious bad PAR values
% fudge factor for auto QC ACS.
% Varies between ACS: 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.RawAutoQCLim.filtered.a = 2; % 6
ila.cfg.qc.RawAutoQCLim.filtered.c = 2; % 15
ila.cfg.qc.RawAutoQCLim.total.a = 2; % 4
ila.cfg.qc.RawAutoQCLim.total.c = 2; % 12
ila.cfg.qc.RawAutoQCLim.dissolved.a = 3; % 4
ila.cfg.qc.RawAutoQCLim.dissolved.c = 3; % 12
% fudge factor for auto QC BB.
% 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.RawAutoQCLim.filtered.bb = 4;
ila.cfg.qc.RawAutoQCLim.total.bb = 4;
% remove saturated periods in BB
ila.cfg.qc.Saturation_Threshold_bb = 4000; % saturate above 4000 counts
ila.RawAutoQC(ila.cfg.qc.RawAutoQCLim, ila.cfg.qc.Saturation_Threshold_bb, 'raw');
% ila.CheckDataStatus();

%% Normal and DI Diagnostic Plot
% check raw spectrums AC or BB sensors
ila.DiagnosticPlot('BB',{'raw'}); % AC or BB

%% 6.Normal Bin
% % Set settings directly in configuration file (no tunning at this step)
% % run before re-bin only to clear qc tables
% ila.instrument.ACS57.qc.tsw = table(); ila.instrument.ACS57.qc.fsw = table();
ila.cfg.bin.skip = {};
ila.Bin()
ila.CheckDataStatus();

%% Normal and DI Diagnostic Plot
% check binned spectrums AC or BB sensors
ila.DiagnosticPlot('BB',{'bin'}); % AC or BB

%% 6.1.Normal Write bin
ila.Write('bin')
ila.CheckDataStatus();

%% Load processed data from mat files: 'data' = Raw | 'bin' = Bin | 'qc' = QCed | 'prod' = product
% ila.Read('data');
ila.Read('bin');
ila.Read('qc');
% ila.Read('prod');

%% 7.Normal Flag
ila.Flag() % Now deprecated will just copy data to next level
ila.CheckDataStatus();

%% 8.DI. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QC DI
ila.cfg.di.qc.mode = 'ui';
ila.QCDI();

%% 8.1.DI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second auto QC DI
ila.RawAutoQC(ila.cfg.qc.RawAutoQCLim, ila.cfg.qc.Saturation_Threshold_bb, 'qc');

%% 8.2.DI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Write QC
ila.Write('qc')
ila.CheckDataStatus();

%% 8.Normal QC
% Interactive or Loading previous qc selection
ila.cfg.qc.mode='ui';  % load or ui
ila.cfg.qc.specific.run = ila.cfg.instruments2run(~contains(ila.cfg.instruments2run, 'FLOW')); % 'FLOW','ACS57','TSG', 'BB31502', 'WSCD859','PAR'
% QCmap(ila.cfg.days2run); % SST & latlon QC
ila.QC();
ila.CheckDataStatus();

%% Normal and DI Diagnostic Plot
% check QCed spectrums AC or BB sensors
% {'raw','bin','qc','prod'}
ila.DiagnosticPlot('BB',{'qc'}); % AC or BB

%% 8.1.Normal Write qc
ila.Write('qc')
ila.CheckDataStatus();

%% 9.DI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BIN DI
% ila.cfg.days2run = ila.cfg.days2run(1)+1:ila.cfg.days2run(end)-1;
ila.BinDI();

%% Normal and DI Diagnostic Plot
% check binned spectrums AC or BB sensors
ila.DiagnosticPlot('AC',{'bin'}); % AC or BB

%% 9.1. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Write bin DI
ila.Write('bin')
ila.CheckDataStatus();

%% 10.Normal and DI Calibrate
% ila.cfg.calibrate.ACS.compute_dissolved = false;
% ila.cfg.calibrate.BB3.compute_dissolved = false;
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

%% 10.1 Normal and DI prod QC plots
% save_figures = true;

%%% ACS 3D %%%
ila.DiagnosticPlot('BB',{'prod'}); % AC or BB

%%% BB 3D %%%
% ila.DiagnosticPlot('BB',{'prod'}); % AC or BB

%%% ACS BB3 TSG PAR WSCD final product visualisation %%%
ila.visProd_timeseries()

%%
% saveGraph('HBB_1minutes_binned_1-4spectre_per_bin', 'fig')
% saveGraph('HBB_1minutes_binned_1-4spectre_per_bin', 'fig')
% saveGraph('HBB_1minutes_binned_1-4spectre_per_bin', 'fig')
%% 11. Normal and DI Save products
ila.Write('prod')

% % Notify with a song that the job is done
% notif_sound = load('gong'); sound(notif_sound.y, notif_sound.Fs); % handel
% return




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
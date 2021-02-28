% Main InLine Analysis Script
% author: Guillaume Bourdin
% created: Jan 05, 2021
cd('C:\Users\Gui\Documents\MATLAB\InLineAnalysis\InLineAnalysis-master\');

% Load InLineAnalysis and the configuration
ila = InLineAnalysis('cfg\TaraMicrobiome_cfg.m');
% Quick Cfg update

%% TSG
% ila.cfg.days2run = datenum(2020,12,12,0,0,0):datenum(2021,1,5,0,0,0);

%% ACS57
% ila.cfg.days2run = datenum(2020,12,26,0,0,0):datenum(2021,1,5,0,0,0);
ila.cfg.days2run = datenum(2021,1,6,0,0,0):datenum(2021,1,20,0,0,0);
% ila.cfg.days2run = datenum(2021,1,20,0,0,0):datenum(2021,2,5,0,0,0);

%% BB31502
% ila.cfg.days2run = datenum(2020,12,12,0,0,0):datenum(2021,1,5,0,0,0);

%% WSCD859
% ila.cfg.days2run = datenum(2020,12,12,0,0,0):datenum(2021,1,5,0,0,0);

%%
ila.cfg.instruments2run = {'FTH','ACS57'}; % 'FTH','ACS57','TSG', 'BB31502', 'WSCD859','PAR'
ila.cfg.qcref.view = 'ACS57';
ila.cfg.parallel = Inf;

%% 1. Import | Load raw data
ila.cfg.force_import = false;
ila.ReadRaw();
ila.CheckDataStatus();

%% 1.1 Read & QC DI
% To get ag and cg from ACS can run DI day by day
% To get betag from BB3 need to run the entire dataset if select di method constant
% ila.cfg.days2run = [ila.cfg.days2run(1)-1:ila.cfg.days2run(end)+1];
% ila.ReadRawDI();
% ila.QCDI();
% ila.BinDI();
% ila.cfg.days2run = [ila.cfg.days2run(1)+1:ila.cfg.days2run(end)-1];

%% 2. Synchronise instruments
% % Independent of flow rate (for now)
% % If flow rate varies use the Strech method
% % Play with delay of synchronisation
% % TSG is assumed to be set at zero
% % No noticeable difference was observed between the TSG of EXPORTS and the BB3
% % ila.instrument.FTH.Sync(30);
% ila.instrument.ACS007.Sync(10);
% % ila.instrument.ACS091.Sync(65);
% % ila.instrument.ACS111.Sync(60);
% % ila.instrument.ACS279.Sync(55);
% % ila.instrument.BB3.Sync(0);
% % % ila.instrument.LISST.Sync(1);
% % ila.instrument.WSCD1082P.Sync(40);
% % ila.instrument.TSG.Sync(0);
% % ila.instrument.ALFA.Sync(15); 
% % Quick visualizzation to sync with TSG
% % fig(30, 'sync TSG');
% % yyaxis('left'); plot(ila.instrument.TSG.data.dt, ila.instrument.TSG.data.t); ylabel('Temperature (^o C)');
% % yyaxis('right'); plot(ila.instrument.BB3.data.dt, ila.instrument.BB3.data.beta(:,2)); ylabel('\beta (m^{-1} sr^{-1})'); ylim([80 300]);
% % datetick2_doy();
% visSync(ila.instrument.FTH.data, ila.instrument.ACS007.data.dt, ila.instrument.ACS007.data.a(:,20), 'a (m^{-1})');
% visSync(ila.instrument.FTH.data, ila.instrument.ACS007.data.dt, ila.instrument.ACS007.data.c(:,40), 'c (m^{-1})');
% % visSync(ila.instrument.FTH.data, ila.instrument.ACS091.data.dt, ila.instrument.ACS091.data.a(:,20), 'a (m^{-1})');
% % visSync(ila.instrument.FTH.data, ila.instrument.ACS091.data.dt, ila.instrument.ACS091.data.c(:,40), 'c (m^{-1})');
% % visSync(ila.instrument.FTH.data, ila.instrument.ACS111.data.dt, ila.instrument.ACS111.data.a(:,20), 'a (m^{-1})');
% % visSync(ila.instrument.FTH.data, ila.instrument.ACS111.data.dt, ila.instrument.ACS111.data.c(:,40), 'c (m^{-1})');
% % visSync(ila.instrument.FTH.data, ila.instrument.ACS279.data.dt, ila.instrument.ACS279.data.a(:,20), 'a (m^{-1})');
% % visSync(ila.instrument.FTH.data, ila.instrument.ACS279.data.dt, ila.instrument.ACS279.data.c(:,40), 'c (m^{-1})');
% % visSync(ila.instrument.('FTH').data, ila.instrument.('BB3').data.dt, ila.instrument.('BB3').data.beta(:,1), '\beta (counts)');
% % visSync(ila.instrument.('FTH').data, ila.instrument.('LISST').data.dt, ila.instrument.('LISST').data.beta(:,10), '\beta (counts)');
% % visSync(ila.instrument.('FTH').data, ila.instrument.('WSCD').data.dt, ila.instrument.('WSCD').data.fdom, 'FDOM (counts)');
% % visSync(ila.instrument.('FTH').data, ila.instrument.('WSCD1082P').data.dt, ila.instrument.('WSCD1082P').data.fdom, 'FDOM (counts)');
% % visSync(ila.instrument.('BB3').data, ila.instrument.('TSG').data.dt, ila.instrument.('TSG').data.t, 'Temp (C)');
% % % visSync(ila.instrument.FTH.data, ila.instrument.ALFA.data.dt, ila.instrument.ALFA.data.Chlb, 'chlb');yyaxis('left'); ylim([0 2]);
% % 
% % % xlim([datenum(2018,08,14,9,55,0) datenum(2018,08,14,11,05,0)]);
% % % ylim([-0.1 0.2]);
% % % Once settings are good set them in the configuration file.
% % % The software is now doing the same with one line of code.
% % ila.Sync()
% % % ila.instrument.BB31502.Sync(-90);
% % % ila.instrument.BB31502.Sync(-10);

%% Automatic detection of filter events for AC and BB sensors
% data = ila.instrument.(ila.cfg.qcref.view).data;
% instrument = ila.cfg.qcref.view;
% FTH = ila.instrument.FTH.data;
% MinFiltPeriod=ila.cfg.qcref.MinFiltPeriod;

ila.cfg.qcref.MinFiltPeriod = 65; % filter even period in minute % ACS: 55 % BB3: 60
ila.SplitDetect(ila.cfg.qcref.MinFiltPeriod);

%% 3. QC Reference
% run with mode ui during first run (it saves your work for the next run)
% run with mode load to load previous QC
% Note: when redoing QC of a given period of time (days2run) the previous
% QC during the same period of time is erased, QC done on other periods of
% time is kept in the json file
ila.cfg.qcref.mode='ui'; % 'ui' or 'load'
ila.QCRef();

%% 4. Split fsw and tsw
ila.cfg.split.skip = {'FTH','TSG'};
ila.Split();
ila.CheckDataStatus();

%% Diagnostic Plot
% check raw spectrums AC or BB sensors
ila.DiagnosticPlot('AC',{'raw'}); % AC or BB

%% automatic QC of raw data for step in ACS spectrum, BB saturated and obvious bad PAR values
% fudge factor for auto QC ACS. Varies between ACS must be >= 3 (default = 3 = maximum filtration)
ila.cfg.qc.StepQCLim.filtered.a = 6;
ila.cfg.qc.StepQCLim.filtered.c = 15;
ila.cfg.qc.StepQCLim.total.a = 4;
ila.cfg.qc.StepQCLim.total.c = 12;
% fudge factor for auto QC BB. Must be >= 3 (default = 3 = maximum filtration)
ila.cfg.qc.StepQCLim.bb = 3;
% remove saturated periods
ila.cfg.qc.Saturation_Threshold_bb = 4000;
ila.StepQC(ila.cfg.qc.StepQCLim, ila.cfg.qc.Saturation_Threshold_bb);
ila.CheckDataStatus();

%% Diagnostic Plot
% check raw spectrums AC or BB sensors
ila.DiagnosticPlot('AC',{'raw'}); % AC or BB

%% 5. Bin
% % Set settings directly in configuration file (no tunning at this step)
% % run before re-bin only to clear qc tables
ila.instrument.ACS57.qc.tsw = table(); ila.instrument.ACS57.qc.fsw = table();
ila.cfg.bin.skip = {};
ila.Bin()

%% Diagnostic Plot
% check binned spectrums AC or BB sensors
ila.DiagnosticPlot('AC',{'bin'}); % AC or BB

%% Write bin
ila.Write('bin')
ila.CheckDataStatus();

%% 5.1 Read Raw and Bin data
% ila.Read('raw');
% ila.Read('bin');
% ila.Read('qc');
% ila.Read('prod');

% ila.instrument.ACS091.ReadDeviceFile()

%% 6. Flag
ila.Flag() % Now deprecated will just copy data to next level

%% 7. QC
% Interactive or Loading previous qc selection
ila.cfg.qc.mode='ui';  % load or ui
ila.cfg.qc.specific.run = {'ACS57'}; % 'FTH','ACS57','TSG', 'BB31502', 'WSCD859','PAR'
% QCmap(ila.cfg.days2run); % SST & latlon QC
ila.QC();
ila.CheckDataStatus();

%% Diagnostic Plot
% check QCed spectrums AC or BB sensors
% {'raw','bin','qc','prod'}
ila.DiagnosticPlot('AC',{'qc'}); % AC or BB

%% Write qc
ila.Write('qc')
ila.CheckDataStatus();

%% 8. Calibrate
% ila.cfg.calibrate.ACS.compute_dissolved = false;
% ila.cfg.calibrate.BB3.compute_dissolved = false;

ila.Calibrate();
ila.CheckDataStatus();
% % Calibrate LISST only
% % ila.instrument.LISST.inversion = 'non-spherical';
% % ila.instrument.LISST.Calibrate()
% % Calibrate ACS only
% ila.instrument.ACS111.Calibrate(ila.cfg.calibrate.ACS111.compute_dissolved,...
%                              ila.cfg.calibrate.ACS111.interpolation_method,...
%                              ila.instrument.(ila.cfg.calibrate.ACS111.CDOM_source),...
%                              ila.instrument.(ila.cfg.calibrate.ACS111.FTH_source))
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

%% 3D QC plots
% save_figures = true;

%%% ACS %%%
ila.DiagnosticPlot('AC',{'prod'}); % AC or BB

%%% BB %%%
ila.DiagnosticPlot('BB',{'prod'}); % AC or BB

%%% PAR %%%
% figure(77);
% scatter(datetime(ila.instrument.PAR.prod.a.dt,'ConvertFrom','datenum'), ila.instrument.PAR.prod.a.par, 4, 'filled'); ylabel('PAR (\muE.m^-^2.s^-^1)');

%% WSCD %%%
figure(78);
scatter(datetime(ila.instrument.WSCD859.prod.pd.dt,'ConvertFrom','datenum'), ila.instrument.WSCD859.prod.pd.fdom, 4, 'filled'); ylabel('fdom ppb');

%%% TSG %%%
figure(79);
yyaxis('left')
scatter(datetime(ila.instrument.TSG.prod.a.dt,'ConvertFrom','datenum'), ila.instrument.TSG.prod.a.t, 4, 'filled'); ylabel('TSG T (DegC)');
yyaxis('right')
scatter(datetime(ila.instrument.TSG.prod.a.dt,'ConvertFrom','datenum'), ila.instrument.TSG.prod.a.s, 4, 'filled'); ylabel('TSG S (PSU)');

%% 9. Save products
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

%%% ACS Chl %%%
% fig(79);
figure(79); hold('on');
yyaxis('left'); plot(ila.instrument.ACS57.prod.p.dt, ila.instrument.ACS57.prod.p.chl, '.-'); ylabel('Chl (\mug L^{-1})');
yyaxis('right'); plot(ila.instrument.ACS57.prod.p.dt, ila.instrument.ACS57.prod.p.gamma, '.-'); ylabel('\gamma');
datetick2_doy(); set(datacursormode(figure(79)),'UpdateFcn',@data_cursor_display_date);
set(gca, 'FontSize', 14, 'FontName', 'Helvetica Neue');
% if save_figures; savefig([ila.instrument.ACS.path.prod 'ACS' datestr(ila.cfg.days2run(1), 'yyyymmdd') '_chl_gamma']); end
% save_fig([ila.instrument.ACS.path.prod 'ACS_stnP_chl_gamma'], 1280, 600);

% Check ACS wl registration
fig(78); hold('on'); plot(wl, ACS111.p.ap(2000:2010, :), 'o-'); plot([676 676], ylim(), 'k'); plot([675 675], ylim(), 'k');



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

% Main InLine Analysis Script
% author: Nils Haentjens
% created: Oct 22, 2018

% Load InLineAnalysis and the configuration
ila = InLineAnalysis('cfg/default_cfg.m');
% Quick Cfg update
ila.cfg.days2run = datenum(2018,08,20):datenum(2018,09,12);
% ila.cfg.instruments2run = {'FTH', 'TSG', 'BB3', 'LISST', 'WSCD1299', 'ALFA'};
ila.cfg.instruments2run = {'FTH', 'ACS301'};

%% 1. Import | Load raw data
ila.cfg.force_import = false;
ila.ReadRaw();

%% 1.1 Read & QC DI
% To get ag and cg from ACS can run DI day by day
% To get betag from BB3 need to run the entire dataset if select di method constant
ila.cfg.days2run = [ila.cfg.days2run(1)-1:ila.cfg.days2run(end)+1];
ila.ReadRawDI();
ila.QCDI();
ila.BinDI();
ila.cfg.days2run = [ila.cfg.days2run(1)+1:ila.cfg.days2run(end)-1];

%% 2. Synchronise instruments
% Independent of flow rate (for now)
% If flow rate varies use the Strech method
% Play with delay of synchronisation
% TSG is assumed to be set at zero
% No noticeable difference was observed between the TSG of EXPORTS and the BB3
% ila.instrument.FTH.Sync(30);
% ila.instrument.ACS.Sync(97);
% ila.instrument.ACS298.Sync(93);
%ila.instrument.BB3.Sync(-90);
% ila.instrument.LISST.Sync(1);
% ila.instrument.WSCD.Sync(5);
% ila.instrument.ALFA.Sync(15); 
% Quick visualizzation to sync with TSG
% fig(30, 'sync TSG');
% yyaxis('left'); plot(ila.instrument.TSG.data.dt, ila.instrument.TSG.data.t); ylabel('Temperature (^o C)');
% yyaxis('right'); plot(ila.instrument.BB3.data.dt, ila.instrument.BB3.data.beta(:,2)); ylabel('\beta (m^{-1} sr^{-1})'); ylim([80 300]);
% datetick2_doy();
% visSync(ila.instrument.FTH.data, ila.instrument.ACS.data.dt, ila.instrument.ACS.data.a(:,20), 'a (m^{-1})');
% visSync(ila.instrument.FTH.data, ila.instrument.ACS298.data.dt, ila.instrument.ACS298.data.c(:,40), 'c (m^{-1})');
% visSync(ila.instrument.('FTH').data, ila.instrument.('BB3').data.dt, ila.instrument.('BB3').data.beta(:,1), '\beta (counts)');
% visSync(ila.instrument.('FTH').data, ila.instrument.('LISST').data.dt, ila.instrument.('LISST').data.beta(:,10), '\beta (counts)');
% visSync(ila.instrument.('FTH').data, ila.instrument.('WSCD').data.dt, ila.instrument.('WSCD').data.fdom, 'FDOM (counts)');
% visSync(ila.instrument.('FTH').data, ila.instrument.('WSCD1299').data.dt, ila.instrument.('WSCD1299').data.fdom, 'FDOM (counts)');
% visSync(ila.instrument.FTH.data, ila.instrument.ALFA.data.dt, ila.instrument.ALFA.data.Chlb, 'chlb');yyaxis('left'); ylim([0 2]);

% xlim([datenum(2018,08,14,9,55,0) datenum(2018,08,14,11,05,0)]);
% ylim([-0.1 0.2]);
% Once settings are good set them in the configuration file.
% The software is now doing the same with one line of code.
ila.Sync()
%ila.instrument.BB3.Sync(-90);
% ila.instrument.BB3.Sync(-10);

%% 3. QC Reference
% run with mode ui during first run (it saves your work for the next run)
% run with mode load to load previous QC
% Note: when redoing QC of a given period of time (days2run) the previous
% QC during the same period of time is erased, QC done on other periods of
% time is kept in the json file
ila.cfg.qcref.mode='load'; % 'ui' or 'load'
ila.QCRef();

%% 4. Split fsw and tsw
% ila.instrument.ACS.Split(ila.instrument.FTH, [120, 30])
% i=20; visSplit(ila.instrument.FTH.data,...
%               ila.instrument.ACS.raw.tsw.dt, ila.instrument.ACS.raw.tsw.a(:,i),...
%               ila.instrument.ACS.raw.fsw.dt, ila.instrument.ACS.raw.fsw.a(:,i),...
%               ila.instrument.ACS.raw.bad.dt, ila.instrument.ACS.raw.bad.a(:,i),...
%               'a (m^{-1})'); % ylim([-0.02 0.04]);
% i=20; visSplit(ila.instrument.FTH.data,...
%               ila.instrument.ACS301.raw.tsw.dt, ila.instrument.ACS301.raw.tsw.a(:,i),...
%               ila.instrument.ACS301.raw.fsw.dt, ila.instrument.ACS301.raw.fsw.a(:,i),...
%               ila.instrument.ACS301.raw.bad.dt, ila.instrument.ACS301.raw.bad.a(:,i),...
%               'a (m^{-1})'); % ylim([-0.02 0.04]);
% ila.instrument.BB3.Split(ila.instrument.FTH, [420, 220])
% i=1; visSplit(ila.instrument.FTH.data,...
%               ila.instrument.BB3.raw.tsw.dt, ila.instrument.BB3.raw.tsw.beta(:,i),...
%               ila.instrument.BB3.raw.fsw.dt, ila.instrument.BB3.raw.fsw.beta(:,i),...
%               ila.instrument.BB3.raw.bad.dt, ila.instrument.BB3.raw.bad.beta(:,i),...
%               '\beta (counts)');
% ila.instrument.LISST.Split(ila.instrument.FTH, [540, 360])
% i=15; visSplit(ila.instrument.FTH.data,...
%               ila.instrument.LISST.raw.tsw.dt, ila.instrument.LISST.raw.tsw.beta(:,i),...
%               ila.instrument.LISST.raw.fsw.dt, ila.instrument.LISST.raw.fsw.beta(:,i),...
%               ila.instrument.LISST.raw.bad.dt, ila.instrument.LISST.raw.bad.beta(:,i),...
%               '\beta (counts)');
% ila.instrument.WSCD.Split(ila.instrument.FTH, [310, 20])
% visSplit(ila.instrument.FTH.data,...
%         ila.instrument.WSCD.raw.tsw.dt, ila.instrument.WSCD.raw.tsw.fdom,...
%         [], [],...
%         ila.instrument.WSCD.raw.bad.dt, ila.instrument.WSCD.raw.bad.fdom,...
%         'FDOM (counts)');
% ila.instrument.ALFA.Split(ila.instrument.FTH, [120, 30])
% visSplit(ila.instrument.FTH.data,...
%         ila.instrument.ALFA.raw.tsw.dt, ila.instrument.ALFA.raw.tsw.Chlb,...
%         ila.instrument.ALFA.raw.fsw.dt, ila.instrument.ALFA.raw.fsw.Chlb,...
%         ila.instrument.ALFA.raw.bad.dt, ila.instrument.ALFA.raw.bad.Chlb,...
%         'Chlb (counts)'); yyaxis('left'); ylim([0 2]);
ila.Split();

% i=20; visSplit(ila.instrument.FTH.data,...
%               ila.instrument.ACS.raw.tsw.dt, ila.instrument.ACS.raw.tsw.a(:,i),...
%               ila.instrument.ACS.raw.fsw.dt, ila.instrument.ACS.raw.fsw.a(:,i),...
%               ila.instrument.ACS.raw.bad.dt, ila.instrument.ACS.raw.bad.a(:,i),...
%               'a (m^{-1})'); % ylim([-0.02 0.04]);
% visSplit(ila.instrument.FTH.data,...
%               ila.instrument.BB3.raw.tsw.dt, ila.instrument.BB3.raw.tsw.beta(:,2),...
%               ila.instrument.BB3.raw.fsw.dt, ila.instrument.BB3.raw.fsw.beta(:,2),...
%               ila.instrument.BB3.raw.bad.dt, ila.instrument.BB3.raw.bad.beta(:,2),...
%               '\beta (counts)');
% return
%% 5. Bin
% Set settings directly in configuration file (no tunning at this step)
ila.Bin()
% ila.Write('bin')
% return

%% 5.1 Read Raw and Bin data
% ila.instrument.ACS.ReadDeviceFile()
% ila.Read('raw');
% ila.Read('bin');

%% 6. Flag
% Visualize flag parameter before applying them
% Global params
% params.maximum_fudge_factor = 4;
% params.variance_fudge_factor = 3;
% params.avg_sensitivity = 1;
% params.unc1_sensitivity = 1;
% params.unc2_sensitivity = 2;
% params.smooth_threshold = 60;
% ACS
% params.abs_uncertainty = 0.004;
% params.rel_uncertainty = 0.0125;
% params.min_flag_n = floor(size(ila.instrument.ACS.lambda_ref,2)*0.8);
% params.primary_varname = 'c';
% visFlagParams(params, ila.instrument.ACS.bin.tsw, 20); ylim([0 0.02]);
% params_fsw = params;
% params_fsw.abs_uncertainty = 0.001;
% params_fsw.smooth_threshold = 2;
% visFlagParams(params_fsw, ila.instrument.ACS.bin.fsw, 30);
% % BB3
% params.abs_uncertainty = 2;
% params.rel_uncertainty = 0.05;
% params.min_flag_n = 2;
% params.primary_varname = 'beta';
% visFlagParams(params, ila.instrument.BB3.bin.tsw, 1);
% params_fsw = params;
% params_fsw.smooth_threshold = 2;
% visFlagParams(params_fsw, ila.instrument.BB3.bin.fsw, 'beta', 1);
% % LISST
% params.abs_uncertainty = 1;
% params.rel_uncertainty = 0.05;
% params.min_flag_n = 27;
% params.smooth_threshold = 6;
% params.primary_varname = 'beta';
% visFlagParams(params, ila.instrument.LISST.bin.tsw, 15);
% params_fsw = params;
% params_fsw.smooth_threshold = 2; 
% visFlagParams(params, ila.instrument.LISST.bin.fsw, 'beta', 15);
% % WSCD
% params.abs_uncertainty=0.3;
% params.rel_uncertainty=0.0125;
% visFlagParams(params, ila.instrument.WSCD.bin.tsw, 'fdom', 1);
% ALFA
% params.abs_uncertainty=0;
% params.rel_uncertainty=0.25;
% params.primary_varname = 'FvFm';
% visFlagParams(params, ila.instrument.ALFA.bin.tsw, 1); ylim([0 10]);

% Update flag parameters live
% ila.cfg.flag.BB3.tot.rel_uncertainty = 0.05;
% ila.cfg.flag.LISST.tot.rel_uncertainty = 0.05;
% cfg.process.flag.ALFA.rel_uncertainty = 0.25;
ila.Flag() % Now deprecated will just copy data to next level
 
% visFlag(ila.instrument.ACS.raw.tsw, ila.instrument.ACS.raw.fsw,...
%         ila.instrument.ACS.qc.tsw, ila.instrument.ACS.suspect.tsw,...
%         ila.instrument.ACS.qc.fsw, ila.instrument.ACS.suspect.fsw,...
%         'c', 67); %ylim([-0.1 0.3]);
% visFlag(ila.instrument.ALFA.raw.tsw, ila.instrument.ALFA.raw.fsw,...
%         ila.instrument.ALFA.qc.tsw, ila.instrument.ALFA.suspect.tsw,...
%         ila.instrument.ALFA.qc.fsw, ila.instrument.ALFA.suspect.fsw,...
%         'FvFmG', 1); ylim([-0.1 10]); % FvFm % Chlb

% return
%% 7. QC
% Interactive or Loading previous qc selection
ila.cfg.qc.mode='load';  % load or ui
% cfg.process.qc.specific.run = {'ACS298', 'ACS301', 'BB3', 'LISST', 'WSCD859', 'WSCD1299', 'ALFA'};
% ila.cfg.qc.specific.run = {'TSG', 'ACS', 'BB3', 'LISST', 'WSCD', 'ALFA'};
% ila.instrument.TSG.view.varname = 's';
% QC a few
% ila.cfg.qc.specific.run = {'BB3'};
% re-QC only LISST
% ila.cfg.qc.specific.run = {'LISST'};
% re-QC only ACS
% ila.cfg.qc.specific.run = {'ACS'};
% ila.instrument.ACS.view.varname = 'c';
% QC TSG only
% ila.cfg.qc.global.active = false;
% ila.cfg.qc.specific.run = {'TSG'};
% QC WSCD
% ila.cfg.qc.global.active = false;
% ila.cfg.qc.specific.run = {'WSCD'};
% QC ALFA
% ila.cfg.qc.specific.run = {'ALFA'};
% Run QC

% fig(51);
% scatter(ila.instrument.TSG.data.dt, ila.instrument.TSG.data.t, 'filled');
% ylabel('t (^{o}C)');
% datetick2_doy();
% visSync(ila.instrument.FTH.data, ila.instrument.TSG.data.dt, ila.instrument.TSG.data.s, 's (PSU)');
% visSync(ila.instrument.FTH.data, ila.instrument.TSG.data.dt, ila.instrument.TSG.data.fchl, 'fchl (counts)');
ila.QC();

ila.Write('qc')
% ila.Read('qc')

%% 8. Calibrate
% ila.cfg.calibrate.ACS.compute_dissolved = true;
% ila.cfg.calibrate.BB3.compute_dissolved = true;
ila.Calibrate();
% Calibrate LISST only
% ila.instrument.LISST.inversion = 'non-spherical';
% ila.instrument.LISST.Calibrate()
% Calibrate ACS only
% ila.instrument.ACS.Calibrate(ila.cfg.calibrate.ACS.compute_dissolved,...
%                              ila.cfg.calibrate.ACS.interpolation_method,...
%                              ila.instrument.(ila.cfg.calibrate.ACS.CDOM_source),...
%                              ila.instrument.(ila.cfg.calibrate.ACS.FTH_source))

% Compute NAAMES specific chl
% wl = ila.instrument.ACS.lambda_ref; ACS = ila.instrument.ACS.prod;
% wl = ila.instrument.ACS298.lambda_ref; ACS = ila.instrument.ACS298.prod;
% ap_a=interp1(wl,ACS.p.ap',[650 676 715],'linear')';
% line_height = (ap_a(:,2)-(39/65*ap_a(:,1)+26/65*ap_a(:,3)));
% p.chl=157*line_height.^1.22;
% ila.instrument.ACS.prod.p.chl_naames = 95 * line_height .^ 1.06;
% ila.instrument.ACS.prod.p.chl_naames(real(ila.instrument.ACS.prod.p.chl_naames) ~= ila.instrument.ACS.prod.p.chl_naames) = NaN;
% fprintf('Done\n');

% Compute EXPORTS Specific chl
% Derive Chl (Line heigh at 676 compared to 650 and 715)
% % ila_acs = ila.instrument.ACS298.prod.p; ila_acs_wl = ila.instrument.ACS298.lambda_ref;
% ila_acs = ila.instrument.ACS301.prod.p; ila_acs_wl = ila.instrument.ACS301.lambda_ref;
% ap = interp1(ila_acs_wl, ila_acs.ap',[650 676 715],'linear')';
% line_height = (ap(:,2)-(39/65*ap(:,1)+26/65*ap(:,3)));
% % ila.instrument.ACS298.prod.p.chl_exports = 138.14 * line_height.^1.11; % EXPORTS relation
% ila.instrument.ACS301.prod.p.chl_exports = 138.14 * line_height.^1.11; % EXPORTS relation

%% 9. Save products
% ila.Write('prod')

% Notify with a song that the job is done
notif_sound = load('gong'); sound(notif_sound.y, notif_sound.Fs); % handel
return

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



%%% ACS %%%
% wl = ila.instrument.ACS298.lambda_ref; ACS = ila.instrument.ACS298.prod; ACS_DIW = ila.instrument.ACS298.bin.diw;
wl = ila.instrument.ACS301.lambda_ref; ACS = ila.instrument.ACS301.prod; ACS_DIW = ila.instrument.ACS301.bin.diw;
% wl = ila.instrument.ACS.lambda_ref; ACS = ila.instrument.ACS.prod; ACS_DIW = ila.instrument.ACS.bin.diw;
% visProd3D(wl, ACS.p.dt, ACS.p.ap_sd./ACS.p.ap_n, false); zlabel('ste(a_p) (m^{-1})');
visProd3D(wl, ACS.p.dt, ACS.p.ap, false, 'Wavelength', false, 72); zlabel('a_p (m^{-1})'); %, 'Wavelength', true
visProd3D(wl, ACS.p.dt, ACS.p.cp, false, 'Wavelength', false, 73); zlabel('c_p (m^{-1})');
visProd3D(wl, ACS.g.dt, ACS.g.ag, false, 'Wavelength', false, 74); zlabel('a_g (m^{-1})');
visProd3D(wl, ACS.g.dt, ACS.g.cg, false, 'Wavelength', false, 75); zlabel('c_g (m^{-1})');
visProd3D(wl,ACS_DIW.dt,ACS_DIW.a, false, 'Wavelength', false, 76); view(0,0); zlabel('a_0 (m^-1)');
visProd3D(wl,ACS_DIW.dt,ACS_DIW.c, false, 'Wavelength', false, 77); view(0,0); zlabel('c_0 (m^-1)');
xlabel('\lambda (nm)');
% ylabel(datestr(ila.cfg.days2run(1)));

% if save_figures; figure(72); savefig([ila.instrument.ACS.path.prod 'ACS' datestr(ila.cfg.days2run(1), 'yyyymmdd') '_a_p']); end
% if save_figures; figure(73); savefig([ila.instrument.ACS.path.prod 'ACS' datestr(ila.cfg.days2run(1), 'yyyymmdd') '_c_p']); end
% save_fig([ila.instrument.ACS.path.prod 'ACS_stnP_ap'], 1280, 720);

% Display 676 line on 3D plot
% hold('on'); fill3([676 676 676 676], [736961 736937 736937 736961], [0.14 0.14 -0.02 -0.02], [1 1 1 1])

%%% ACS Chl %%%
% fig(79);
figure(79); hold('on');
yyaxis('left'); plot(ACS.p.dt, ACS.p.chl, '.-'); ylabel('Chl (\mug L^{-1})');
yyaxis('right'); plot(ACS.p.dt, ACS.p.gamma, '.-'); ylabel('\gamma');
% yyaxis('right'); plot(ila.instrument.TSG.prod.a.dt, ila.instrument.TSG.prod.a.fchl, '.-'); ylabel('fchl (counts)');
datetick2_doy(); set(datacursormode(figure(79)),'UpdateFcn',@data_cursor_display_date);
set(gca, 'FontSize', 14, 'FontName', 'Helvetica Neue');
% if save_figures; savefig([ila.instrument.ACS.path.prod 'ACS' datestr(ila.cfg.days2run(1), 'yyyymmdd') '_chl_gamma']); end
% save_fig([ila.instrument.ACS.path.prod 'ACS_stnP_chl_gamma'], 1280, 600);

% Check ACS wl registration
fig(78); hold('on'); plot(wl, ACS.p.ap(2000:2010, :), 'o-'); plot([676 676], ylim(), 'k'); plot([675 675], ylim(), 'k');

%% BB3 particulate %%%
BB3 = ila.instrument.BB3.prod.p;
% visProd3D(ila.instrument.BB3.lambda, ila.instrument.BB3.prod.p.dt, ila.instrument.BB3.prod.p.bbp, false); view(90,0);
% visProd3D(ila.instrument.BB3.lambda, ila.instrument.BB3.prod.g.dt, ila.instrument.BB3.prod.g.betag, false); view(90,0);
% visProd3D(ila.instrument.BB3.lambda, ila.instrument.BB3.bin.diw.dt, ila.instrument.BB3.bin.diw.beta, false); view(90,0);
fig(80, 'beta_p'); hold('on'); CS = lines(5); csi = [1,5,2];
for i=1:3; plot(ila.instrument.BB3.prod.p.dt, ila.instrument.BB3.prod.p.betap(:,i), '.-', 'Color', CS(csi(i),:)); end
datetick2_doy(); set(datacursormode(figure(80)),'UpdateFcn',@data_cursor_display_date);
fig(81, 'beta_g'); hold('on'); CS = lines(5); csi = [1,5,2];
for i=1:3; plot(ila.instrument.BB3.prod.g.dt, ila.instrument.BB3.prod.g.betag(:,i), '.-', 'Color', CS(csi(i),:)); end
plot(xlim(), [0 0], 'k');
datetick2_doy(); set(datacursormode(figure(81)),'UpdateFcn',@data_cursor_display_date);
% fig(82, 'beta_0'); hold('on'); CS = lines(5); csi = [1,5,2];
% for i=1:3; plot(ila.instrument.BB3.bin.diw.dt, ila.instrument.BB3.bin.diw.beta(:,i), '.-', 'Color', CS(csi(i),:)); end
% datetick2_doy(); set(datacursormode(figure(82)),'UpdateFcn',@data_cursor_display_date);
% pause(); 
% if save_figures; figure(80); savefig([ila.instrument.BB3.path.prod 'BB3' datestr(ila.cfg.days2run(1), 'yyyymmdd') '_bbp']); end

%% WSCD %%%
WSCD = ila.instrument.WSCD1299.prod.a;
fig(63);
plot(WSCD.dt, WSCD.fdom, '.-'); datetick2_doy();
ylabel('FDOM (counts)');
% if save_figures; savefig([ila.instrument.WSCD.path.prod 'WSCD' datestr(ila.cfg.days2run(1), 'yyyymmdd') '_FDOM']); end

%%% TSG %%%
fig(65); TSG = ila.instrument.TSG.prod.a;
yyaxis('right'); plot(TSG.dt, TSG.par, '.-'); ylabel('PAR (\muE/s/m^2)');
yyaxis('left'); plot(TSG.dt, TSG.t, '.-'); ylabel('Temperature (^oC)');
datetick2_doy(); set(datacursormode(figure(65)),'UpdateFcn',@data_cursor_display_date);
% xlabel('Aug 15-17, 2018');
% set(gca, 'FontSize', 14, 'FontName', 'Helvetica Neue');
% save_fig([ila.instrument.TSG.path.prod 'Underway_stnP_PAR_T'], 1280, 600);

%% LISST particulate %%%
% ila.Calibrate();

% Plot VSF
sel = ~any(ila.instrument.LISST.prod.p.betap <= 0,2);
% 3D
visProd3D(ila.instrument.LISST.theta, ila.instrument.LISST.prod.p.dt(sel), ila.instrument.LISST.prod.p.betap(sel,:), false, 'Log', false, 81);
set(gca, 'ZScale', 'log', 'XScale', 'log');
xlabel('\theta (^o)'); zlabel('VSF (m^{-1})');
% 2D
visProd2D(ila.instrument.LISST.theta, ila.instrument.LISST.prod.p.dt(sel), ila.instrument.LISST.prod.p.betap(sel,:)); set(gca, 'YScale', 'log');


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



%% ALFA %%%
fig(64); ALFA = ila.instrument.ALFA.prod.a;
yyaxis('right'); plot(ALFA.dt, ALFA.Chlb, '.-'); ylabel('Chl_b Fl (\mug L^{-1})');
yyaxis('left'); hold('on'); ylabel('F_{v}/F_{m}');
plot(ALFA.dt(ALFA.FvFm>0.2), ALFA.FvFm(ALFA.FvFm>0.2), '.-');
plot(ALFA.dt(ALFA.FvFmG>0.2), ALFA.FvFmG(ALFA.FvFmG>0.2), 'g.-'); 
datetick2_doy(); %xlabel('Aug 15-17, 2018');
set(datacursormode(figure(64)),'UpdateFcn',@data_cursor_display_date);
% if save_figures; savefig([ila.instrument.ALFA.path.prod 'ALFA' datestr(ila.cfg.days2run(1), 'yyyymmdd') '_FvFm']); end
% save_fig([ila.instrument.ALFA.path.prod 'ALFA_stnP_FvFm'], 1280, 720);

%% TSG %%%
fig(65); TSG = ila.instrument.TSG.prod.a;
yyaxis('right'); plot(TSG.dt, TSG.par, '.-'); ylabel('PAR (\muE/s/m^2)');
yyaxis('left'); plot(TSG.dt, TSG.t, '.-'); ylabel('Temperature (^oC)');
datetick2_doy(); set(datacursormode(figure(65)),'UpdateFcn',@data_cursor_display_date);
% xlabel('Aug 15-17, 2018');
% set(gca, 'FontSize', 14, 'FontName', 'Helvetica Neue');
% save_fig([ila.instrument.TSG.path.prod 'Underway_stnP_PAR_T'], 1280, 600);

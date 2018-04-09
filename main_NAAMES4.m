% In-Line Analysis main script
% run script entirely or step by step 
% author: nils
% created: Oct 12, 2017

% NAAMES 4; start on March 20, 2018 (737139)
%           end on April... 2018 (

% Load InLineAnalysis
ila = InLineAnalysis('cfg/NAAMES4_cfg.json');
% ila.cfg.days2run = datenum(2018,03,20):datenum(2018,04,07);
ila.cfg.days2run = datenum(2018,04,06):datenum(2018,04,07);%:datenum(2018,04,01);
% ila.cfg.instruments2run = {'ACS', 'TSG', 'FTH'};
% ila.cfg.instruments2run = {'LISST', 'TSG', 'FTH'};
% ila.cfg.instruments2run = {'TSG'};

% TODO: ACS on 2018-03-30 -> too much data removed in the evening
% TODO: Redo 2018-03-03
% TODO: QC Salinity

%% 1. Import | Load raw data
ila.cfg.force_import = true;
ila.Read();

%% 1.1 Load & Bin DI
ila.cfg.days2run = [ila.cfg.days2run(1)-1:ila.cfg.days2run(end)+1];
ila.ReadDI();
ila.cfg.qc.mode = 'load';
ila.QCDI();
ila.BinDI();
ila.cfg.days2run = [ila.cfg.days2run(1)+1:ila.cfg.days2run(end)-1];

%% 2. Synchronise instruments
% Independent of flow rate (for now)
% If flow rate varies then a new method needs to be implemented
% Play with delay of add for synchronisation
% TSG is assumed to be set at zero
% ila.instrument.FTH.Sync(30);
% ila.instrument.ACS.Sync(30+60);
% ila.instrument.BB3.Sync(30);
% ila.instrument.LISST.Sync(30);
% ila.instrument.WSCD.Sync(30);
% Quick visualization
% visSync(ila.instrument.FTH.data, ila.instrument.ACS.data.dt, ila.instrument.ACS.data.a(:,20), 'a (m^{-1})');
% visSync(ila.instrument.FTH.data, ila.instrument.ACS.data.dt, ila.instrument.ACS.data.c(:,40), 'c (m^{-1})');
% visSync(ila.instrument.('FTH').data, ila.instrument.('BB3').data.dt, ila.instrument.('BB3').data.beta(:,1), '\beta (counts)');
% visSync(ila.instrument.('FTH').data, ila.instrument.('LISST').data.dt, ila.instrument.('LISST').data.beta(:,10), '\beta (counts)');
% visSync(ila.instrument.('FTH').data, ila.instrument.('WSCD').data.dt, ila.instrument.('WSCD').data.fdom, '\beta (counts)');
% xlim([datenum(2018,03,23,19,55,0) datenum(2018,03,23,22,05,0)]);
% ylim([-0.1 0.2]);
% Once settings are good set them in the configuration file.
% The software is now doing the same with one line of code.
ila.Sync()

% visSync(ila.instrument.FTH.data, ila.instrument.ACS.data.dt, ila.instrument.ACS.data.a(:,20), 'a (m^{-1})');

%% 3. QC Reference
% run with mode ui during first run (it saves your work for the next run)
% run with mode load to load previous QC
ila.cfg.qcref.mode='ui'; % 'ui' or 'load'
ila.QCRef();

%% 4. Split fsw and tsw
% ila.instrument.ACS.Split(ila.instrument.FTH, [100, 30])
% i=20; visSplit(ila.instrument.FTH.data,...
%               ila.instrument.ACS.raw.tsw.dt, ila.instrument.ACS.raw.tsw.a(:,i),...
%               ila.instrument.ACS.raw.fsw.dt, ila.instrument.ACS.raw.fsw.a(:,i),...
%               ila.instrument.ACS.raw.bad.dt, ila.instrument.ACS.raw.bad.a(:,i),...
%               'a (m^{-1})'); % ylim([-0.02 0.04]);
% xlim([datenum(2018,03,23,19,55,0) datenum(2018,03,23,22,05,0)]);
% ila.instrument.BB3.Split(ila.instrument.FTH, [420, 220])
% i=1; visSplit(ila.instrument.FTH.data,...
%               ila.instrument.BB3.raw.tsw.dt, ila.instrument.BB3.raw.tsw.beta(:,i),...
%               ila.instrument.BB3.raw.fsw.dt, ila.instrument.BB3.raw.fsw.beta(:,i),...
%               ila.instrument.BB3.raw.bad.dt, ila.instrument.BB3.raw.bad.beta(:,i),...
%               '\beta (counts)');
% ila.instrument.LISST.Split(ila.instrument.FTH, [470, 300])
% i=15; visSplit(ila.instrument.FTH.data,...
%               ila.instrument.LISST.raw.tsw.dt, ila.instrument.LISST.raw.tsw.beta(:,i),...
%               ila.instrument.LISST.raw.fsw.dt, ila.instrument.LISST.raw.fsw.beta(:,i),...
%               ila.instrument.LISST.raw.bad.dt, ila.instrument.LISST.raw.bad.beta(:,i),...
%               '\beta (counts)');
% ila.instrument.WSCD.Split(ila.instrument.FTH, [10, 10])
% visSplit(ila.instrument.FTH.data,...
%         ila.instrument.WSCD.raw.tsw.dt, ila.instrument.WSCD.raw.tsw.fdom,...
%         [], [],...
%         ila.instrument.WSCD.raw.bad.dt, ila.instrument.WSCD.raw.bad.fdom,...
%         'FDOM (counts)');
ila.Split();

% i=20; visSplit(ila.instrument.FTH.data,...
%               ila.instrument.ACS.raw.tsw.dt, ila.instrument.ACS.raw.tsw.a(:,i),...
%               ila.instrument.ACS.raw.fsw.dt, ila.instrument.ACS.raw.fsw.a(:,i),...
%               ila.instrument.ACS.raw.bad.dt, ila.instrument.ACS.raw.bad.a(:,i),...
%               'a (m^{-1})'); % ylim([-0.02 0.04]);

%% 5. Bin
% Set settings directly in configuration file (no tunning at this step)
ila.Bin()

%% 6. Flag
% Visualize flag parameter before applying them
% % Global params
% params.maximum_fudge_factor = 4;
% params.variance_fudge_factor = 3;
% params.avg_sensitivity = 1;
% params.unc1_sensitivity = 1;
% params.unc2_sensitivity = 2;
% params.smooth_threshold = 60;
% % ACS
% params.abs_uncertainty = 0.004;
% params.rel_uncertainty = 0.0125;
% params.min_flag_n = floor(size(ila.instrument.ACS.lambda_ref,2)*0.8);
% params.primary_varname = 'c';
% visFlagParams(params, ila.instrument.ACS.bin.tsw, 20); ylim([0 0.02]);
% params_fsw = params;
% params_fsw.abs_uncertainty = 0.001;
% params_fsw.smooth_threshold = 2;
% visFlagParams(params_fsw, ila.instrument.ACS.bin.fsw, 'c', 30);
% BB3
% params.abs_uncertainty = 2;
% params.rel_uncertainty = 0.025;
% params.min_flag_n = 2;
% visFlagParams(params, ila.instrument.BB3.bin.fsw, 'beta', 1);
% params_fsw = params;
% params_fsw.smooth_threshold = 2;
% visFlagParams(params_fsw, ila.instrument.BB3.bin.fsw, 'beta', 1);
% LISST
% params.abs_uncertainty = 1;
% params.rel_uncertainty = 0.005;
% params.min_flag_n = 27;
% params.smooth_threshold = 6;
% % visFlagParams(params, ila.instrument.LISST.bin.tsw, 'beta', 15);
% params_fsw = params;
% params_fsw.smooth_threshold = 2; 
% visFlagParams(params, ila.instrument.LISST.bin.fsw, 'beta', 15);
% WSCD
% params.abs_uncertainty=0.3;
% params.rel_uncertainty=0.0125;
% visFlagParams(params, ila.instrument.WSCD.bin.tsw, 'fdom', 1);

% Flag data
% ila.cfg.flag.ACS.tot = params;
ila.Flag()

% visFlag(ila.instrument.ACS.raw.tsw, ila.instrument.ACS.raw.fsw,...
%         ila.instrument.ACS.qc.tsw, ila.instrument.ACS.suspect.tsw,...
%         ila.instrument.ACS.qc.fsw, ila.instrument.ACS.suspect.fsw,...
%         'a', 67); %ylim([-0.1 0.3]);
% return
%% 7. QC
% Interactive or Loading previous qc selection
ila.cfg.qc.mode='ui';  % load or ui
% ila.instrument.TSG.view.varname = 's';
% QC a few
% ila.cfg.qc.global.active = false;
% ila.cfg.qc.specific.run = {'BB3', 'LISST', 'WSCD'};
% re-QC only LISST
% ila.cfg.qc.global.active = false;
% ila.cfg.qc.specific.run = {'LISST'};
% re-QC only ACS
% ila.cfg.qc.global.active = false;z
% ila.cfg.qc.specific.run = {'ACS'};
% ila.instrument.ACS.view.varname = 'a';
% QC TSG only
% ila.cfg.qc.global.active = false;
% ila.cfg.qc.specific.run = {'TSG'};
% QC WSCD
% ila.cfg.qc.global.active = false;
% ila.cfg.qc.specific.run = {'WSCD'};
% Run QC

% visSync(ila.instrument.FTH.data, ila.instrument.TSG.data.dt, ila.instrument.TSG.data.t, 't (^{o}C)');
% visSync(ila.instrument.FTH.data, ila.instrument.TSG.data.dt, ila.instrument.TSG.data.s, 's (PSU)');
% visSync(ila.instrument.FTH.data, ila.instrument.TSG.data.dt, ila.instrument.TSG.data.fchl, 'fchl (counts)');
ila.QC();

%% 8. Calibrate
ila.cfg.calibrate.ACS.compute_dissolved = true;
ila.cfg.calibrate.BB3.compute_dissolved = true;
ila.Calibrate()
% Calibrate LISST only
% ila.instrument.LISST.inversion = 'non-spherical';
% ila.instrument.LISST.Calibrate()
% Calibrate ACS only
% ila.instrument.ACS.Calibrate(ila.cfg.calibrate.ACS.compute_dissolved,...
%                              ila.cfg.calibrate.ACS.interpolation_method,...
%                              ila.instrument.(ila.cfg.calibrate.ACS.CDOM_source),...
%                              ila.instrument.(ila.cfg.calibrate.ACS.FTH_source))

% TODO reprocess ACS: ap looked like crap, need to re-sync ap
%% ACS
wl = ila.instrument.ACS.lambda_ref; ACS = ila.instrument.ACS.prod;
% visProd3D(wl, ACS.p.dt, ACS.p.ap_sd./ACS.p.ap_n, false); zlabel('ste(a_p) (m^{-1})');
visProd3D(wl, ACS.p.dt, ACS.p.ap, false); zlabel('a_p (m^{-1})'); %, 'Wavelength', true
% visProd3D(wl, ACS.p.dt, ACS.p.cp, false); zlabel('c_p (m^{-1})');
% visProd3D(wl, ACS.g.dt, ACS.g.ag, false); zlabel('a_g (m^{-1})');
xlabel('\lambda (nm)');
ylabel(datestr(ila.cfg.days2run(1)));
%% ACS Chl
fig(63);
yyaxis('left'); plot(ACS.p.dt, ACS.p.chl, '.-'); ylabel('Chl (\mug L^{-1})');
% yyaxis('right'); plot(ACS.p.dt, ACS.p.gamma, '.-'); ylabel('\gamma');
yyaxis('right'); plot(ila.instrument.TSG.prod.a.dt, ila.instrument.TSG.prod.a.fchl, '.-'); ylabel('fchl (counts)');
datetick2_doy(); 
%% BB3 particulate
visProd3D(ila.instrument.BB3.lambda, ila.instrument.BB3.prod.p.dt, ila.instrument.BB3.prod.p.bbp, false); view(90,0);
% visProd3D(ila.instrument.BB3.lambda, ila.instrument.BB3.prod.g.dt, ila.instrument.BB3.prod.g.betag, false); view(90,0);
% pause(); 
%% LISST particulate
sel = 1 <= ila.instrument.LISST.diameters & ila.instrument.LISST.diameters <= 100;
% sel_bin = ~(any(ila.instrument.LISST.prod.p.PSD < 0,2)); % any(~isnan(ila.instrument.LISST.prod.p.PSD),2);
visProd2D(ila.instrument.LISST.diameters(sel), ila.instrument.LISST.prod.p.dt, ila.instrument.LISST.prod.p.PSD(:,sel), false);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Diameter (\mum)'); ylabel('PSD (# mL^{-1} \mum^{-1})'); 
visProd3D(ila.instrument.LISST.diameters(sel), ila.instrument.LISST.prod.p.dt, ila.instrument.LISST.prod.p.PSD(:,sel), false, 'Log');
set(gca, 'ZScale', 'log', 'XScale', 'log');
xlabel('Diameter (\mum)'); zlabel('PSD (# mL^{-1} \mum^{-1})');
%% WSCD
fig(63);
plot(ila.instrument.WSCD.prod.a.dt, ila.instrument.WSCD.prod.a.fdom, '.-'); datetick2_doy();
ylabel('FDOM (counts)');
%% 10. Save products
% ila.cfg.write.skip = ["BB3", "LISST", "WSCD", "TSG", "FTH"];
ila.cfg.write.skip = ["FTH"];
ila.Write()
% ila.SeaBASS()
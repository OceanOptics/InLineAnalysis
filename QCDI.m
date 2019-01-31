% Process DI runs
% author: Nils
% created: April 7, 2018


% QC & BIN step order are switched for DI runs

ila = InLineAnalysis('cfg/default_cfg.m');
ila.cfg.days2run = datenum(2018,08,11):datenum(2018,09,12);
% ila.cfg.days2run = datenum(2018,09,12);
% ila.cfg.instruments2run = {'BB3'};
ila.instrument.ACS = ila.instrument.ACS301;
ila.cfg.instruments2run = {'ACS'};

%% 1. Read
ila.cfg.force_import = false;
ila.ReadRawDI();

%% 2. QC
%%% ACS %%%
ila.instrument.ACS.view.varcol = [3, 40, 50, 66];
ila.instrument.ACS.view.varname = {'a', 'c'};
%%% BB3 %%%
% ila.instrument.BB3.view.varcol = 1:3;
% QC
ila.cfg.di.qc.mode = 'load'; % load | ui | skip
ila.QCDI();

% 2.1 visualize data
%%% ACS %%%
visProd3D(ila.instrument.ACS.lambda_ref,...
          ila.instrument.ACS.qc.diw.dt,...
          ila.instrument.ACS.qc.diw.c, false, 'Wavelength', false, 73); view(0,0); zlabel('c (m^{-1})');
visProd3D(ila.instrument.ACS.lambda_ref,...
          ila.instrument.ACS.qc.diw.dt,...
          ila.instrument.ACS.qc.diw.a, false, 'Wavelength', false, 74); view(0,0); zlabel('a (m^{-1})');
%%% BB3 %%%
% visProd3D(ila.instrument.BB3.lambda,...
%           ila.instrument.BB3.qc.diw.dt,...
%           ila.instrument.BB3.qc.diw.beta); view(50,45);
%% 3. Bin
ila.BinDI();

%% 3.1 Check QC (one day at a time)
%%% ACS %%%
ACS = ila.instrument.ACS; wl = ila.instrument.ACS.lambda_ref;
fig(53, 'c'); hold('on'); 
plot(wl, ACS.qc.diw.c, 'k.');
plot(wl, ACS.bin.diw.c, 'LineWidth', 3);
fig(54, 'a'); hold('on');
plot(wl, ACS.qc.diw.a, 'k.');
plot(wl, ACS.bin.diw.a, 'LineWidth', 3);
%%% BB3 %%%
% BB3 = ila.instrument.BB3; wl = ila.instrument.BB3.lambda;
% fig(53, 'beta'); hold('on');
% plot(wl, BB3.qc.diw.beta, 'k.');
% plot(wl, BB3.bin.diw.beta, 'LineWidth', 3);

% return
%% 3.2 Super check BIN (multiple days at a time)
%%% ACS %%%
visProd3D(ila.instrument.ACS.lambda_ref,...
          ila.instrument.ACS.bin.diw.dt,...
          ila.instrument.ACS.bin.diw.c, false, 'Wavelength', false, 75); view(0,0); zlabel('c (m^{-1})'); datetick('y', 'dd-mmm');
visProd3D(ila.instrument.ACS.lambda_ref,...
          ila.instrument.ACS.bin.diw.dt,...
          ila.instrument.ACS.bin.diw.a, false, 'Wavelength', false, 76); view(0,0);  zlabel('a (m^{-1})'); datetick('y', 'dd-mmm');
%%% BB3 %%%
% fig(73); hold('on'); CS = lines(5); csi = [1,5,2];
% for i=1:3; plot(ila.instrument.BB3.bin.diw.dt, ila.instrument.BB3.bin.diw.beta(:,i), '.-', 'Color', CS(csi(i),:)); end
% datetick2_doy(); set(datacursormode(figure(73)),'UpdateFcn',@data_cursor_display_date);

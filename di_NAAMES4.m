% Process DI runs for NAAMES 4
% author: Nils
% created: April 7, 2018


% QC & BIN step order are switched for DI runs

ila = InLineAnalysis('cfg/NAAMES4_cfg.json');
ila.cfg.days2run = datenum(2018,03,20):datenum(2018,04,08);
% ila.cfg.days2run = datenum(2018,4,8);
ila.cfg.instruments2run = {'ACS'};

%% 1. Read
ila.cfg.force_import = false;
ila.ReadDI();

%% 2. QC
ila.cfg.qc.mode = 'load'; % load | ui | skip
ila.QCDI();

%% 3. Bin
ila.BinDI();

%% 3.1 Visualize
% visFlag(ila.instrument.ACS.qc.diw, [],...
%         ila.instrument.ACS.bin.diw, [], [], [],...
%         'a', 67);
visFlag(ila.instrument.BB3.qc.diw, [],...
        ila.instrument.BB3.bin.diw, [], [], [],...
        'beta', 2);
% visFlag(ila.instrument.LISST.qc.diw, [],...
%         ila.instrument.LISST.bin.diw, [], [], [],...
%         'beta', 15);
% visFlag(ila.instrument.WSCD.qc.diw, [],...
%         ila.instrument.WSCD.bin.diw, [], [], [],...
%         'fdom', 1);

%% 4. Save
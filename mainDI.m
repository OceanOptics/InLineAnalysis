% Main InLine Analysis Script
% author: Guillaume Bourdin
% created: May 05, 2021
% clear
% close all
% if feature('IsDebugMode'); dbquit all; end
% 
% % cd('/Users/emmanuel.boss/Desktop/InLine analysis/InLineAnalysis-master/')
% cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master/')
% 
% % Load InLineAnalysis and the configuration
% ila = InLineAnalysis('cfg/EXPORTS02_cfg.m');

% Quick cfg update
%% update dates
ila.cfg.days2run = datenum(2021,5,1,0,0,0):datenum(2021,5,7,0,0,0);

%%
ila.cfg.instruments2run = {'FLOW','ACS57'}; % 'FLOW', 'TSG', 'BB3', 'HBB', 'WSCD','SUVF','ACS91','LISST'
ila.cfg.qcref.view = 'ACS57';
ila.cfg.parallel = Inf;
ila.cfg.calibrate.(ila.cfg.qcref.view).compute_dissolved = true;

%% 1. Read DI 
% To get ag and cg from ACS can run DI day by day
% To get betag from BB3 need to run the entire dataset if select di method constant
ila.cfg.force_import = false;
ila.ReadRawDI();
ila.CheckDataStatus();

%% Load processed data from mat files: 'data' = Raw | 'bin' = Bin | 'qc' = QCed | 'prod' = product
% ila.Read('data');
ila.Read('bin');
ila.Read('qc');
% ila.Read('prod');

%% 1.1. Diagnostic Plot
% check raw spectrums AC or BB sensors
ila.DiagnosticPlot('AC',{'raw'}); % AC or BB

%% 2. Automatic QC of raw data for step in ACS spectrum, BB saturated and obvious bad PAR values
% fudge factor for auto QC ACS.
% Varies between ACS: 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.RawAutoQCLim.dissolved.a = 3;
ila.cfg.qc.RawAutoQCLim.dissolved.c = 3;
% fudge factor for auto QC BB.
% 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.RawAutoQCLim.dissolved.bb = 3;
ila.RawAutoQC('raw');
ila.CheckDataStatus();

%% 2.1. Diagnostic Plot
% check raw spectrums AC or BB sensors
ila.DiagnosticPlot('AC',{'raw'}); % AC or BB

%% 2.2. Run QC directly on spectra at any level
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
ila.DiagnosticPlot('AC',{'raw'}, false, {'diw','all'});

%% 3. QC DI
ila.cfg.di.qc.mode = 'ui';
ila.cfg.di.qc.remove_old = false;  % remove old selection of this period
ila.cfg.di.qc.qc_once_for_all = false;  % true = QC all variables | false = QC variables separately)
ila.QCDI();

%% 4. Second auto QC DI
ila.RawAutoQC('qc');

%% 4.1. Diagnostic Plot
% check QCed spectrums AC or BB sensors
% {'raw','bin','qc','prod'}
ila.DiagnosticPlot('AC',{'qc'}); % AC or BB

%% 4.2. Run QC directly on spectra at any level
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
ila.DiagnosticPlot('AC',{'qc'}, false, {'diw','beta'});

%% 5. Write QC | write only 'part' or 'diw' or 'all'
ila.Write('qc', 'diw')
ila.CheckDataStatus();

%% 6. BIN DI
ila.BinDI();

%% 6.1. Diagnostic Plot
% check binned spectrums AC or BB sensors
ila.DiagnosticPlot('AC',{'bin'}); % AC or BB

%% 7. Write bin DI | write only 'part' or 'diw' or 'all'
ila.Write('bin', 'diw')
ila.CheckDataStatus();

%% 8. Calibrate
ila.cfg.calibrate.skip = {'FLOW', 'TSG'};
ila.Calibrate();
ila.CheckDataStatus();

%% 8.1. Normal and DI prod QC plots
save_figures = false;

%%% AC or BB 3D plots %%%
ila.DiagnosticPlot('AC', {'prod'}, save_figures); % AC or BB

%%% ACS BB3 TSG PAR WSCD final product visualisation %%%
ila.visProd_timeseries()

%% 9. Run QC directly on spectra at any level
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
ila.DiagnosticPlot('AC',{'prod'}, false, {'g','ag'});

%% 9.1. Load previous qc pick selection at prod level
ila.cfg.qc.mode='load';  % load or ui
ila.cfg.qc.specific.run = {ila.cfg.qcref.view}; % 'FLOW','ACS57','TSG', 'BB31502', 'WSCD859','PAR'
ila.QC();

%% 10. Save products | write only 'part' or 'diw' or 'all'
ila.Write('prod', 'diw')

% % Notify with a song that the job is done
% notif_sound = load('gong'); sound(notif_sound.y, notif_sound.Fs); % handel
% return

%% re-write last version of 'qc' and 'bin' | write only 'part' or 'diw' or 'all'
ila.Write('bin', 'diw')
ila.Write('qc', 'diw')











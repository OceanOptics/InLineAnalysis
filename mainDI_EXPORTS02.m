% Main InLine Analysis Script
% author: Guillaume Bourdin
% created: Jan 05, 2021
% clear
% close all
% if feature('IsDebugMode'); dbquit all; end
% 
% % cd('/Users/emmanuel.boss/Desktop/InLine analysis/InLineAnalysis-master/')
% cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master/')
% 
% % Load InLineAnalysis and the configuration
% ila = InLineAnalysis('cfg/EXPORTS02_cfg.m');

%% Quick cfg update
% ACS91
ila.cfg.days2run = datenum(2021,5,1,0,0,0):datenum(2021,5,7,0,0,0);
% ila.cfg.days2run = datenum(2021,1,6,0,0,0):datenum(2021,1,20,0,0,0);
% ila.cfg.days2run = datenum(2021,1,20,0,0,0):datenum(2021,2,5,0,0,0);

%% BB3
 ila.cfg.days2run = datenum(2021,5,1,0,0,0):datenum(2021,5,7,0,0,0);
% ila.cfg.days2run = datenum(2021,1,6,0,0,0):datenum(2021,1,20,0,0,0);
% ila.cfg.days2run = datenum(2021,1,20,0,0,0):datenum(2021,2,5,0,0,0);

%% SeapointCD
 ila.cfg.days2run = datenum(2021,5,1,0,0,0):datenum(2021,5,7,0,0,0);
 
%% LISST
 ila.cfg.days2run = datenum(2021,5,1,0,0,0):datenum(2021,5,7,0,0,0);
 
%% HBB
 ila.cfg.days2run = datenum(2021,5,1,0,0,0):datenum(2021,5,7,0,0,0);

%%
ila.cfg.instruments2run = {'FLOW','SPCD'}; % 'FLOW', 'TSG', 'BB3', 'HBB', 'WSCD','SPCD','ACS91','LISST'
ila.cfg.qcref.view = 'SPCD';
ila.cfg.parallel = Inf;
ila.cfg.calibrate.(ila.cfg.qcref.view).compute_dissolved = true;

%% 1. Read DI 
% To get ag and cg from ACS can run DI day by day
% To get betag from BB3 need to run the entire dataset if select di method constant
% ila.cfg.days2run = ila.cfg.days2run(1)-1:ila.cfg.days2run(end)+1;
ila.cfg.force_import = false;
ila.ReadRawDI();
ila.CheckDataStatus();

%% 1.1. Diagnostic Plot
% check raw spectrums AC or BB sensors
ila.DiagnosticPlot('BB',{'raw'}); % AC or BB

%% 2. Automatic QC of raw data for step in ACS spectrum, BB saturated and obvious bad PAR values
% fudge factor for auto QC ACS.
% Varies between ACS: 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.RawAutoQCLim.dissolved.a = 3;
ila.cfg.qc.RawAutoQCLim.dissolved.c = 3;
% fudge factor for auto QC BB.
% 0.1 = maximum filtration and >> 10 = very small filtration (default = 3)
ila.cfg.qc.RawAutoQCLim.dissolved.bb = 3;
% remove saturated periods in BB
ila.cfg.qc.Saturation_Threshold_bb = 4000; % saturate above 4000 counts
ila.RawAutoQC(ila.cfg.qc.RawAutoQCLim, ila.cfg.qc.Saturation_Threshold_bb, 'raw');
ila.CheckDataStatus();

%% 2.1. Diagnostic Plot
% check raw spectrums AC or BB sensors
ila.DiagnosticPlot('BB',{'raw'}); % AC or BB

%% Load processed data from mat files: 'data' = Raw | 'bin' = Bin | 'qc' = QCed | 'prod' = product
% ila.Read('data');
ila.Read('bin');
ila.Read('qc');
% ila.Read('prod');

%% 3. QC DI
ila.cfg.di.qc.mode = 'ui';
ila.cfg.di.qc.qc_once_for_AandC = false;  % true = QC 'a' and 'c' together | false = QC 'a' and 'c' separately
ila.QCDI();

%% 4. Second auto QC DI
ila.RawAutoQC(ila.cfg.qc.RawAutoQCLim, ila.cfg.qc.Saturation_Threshold_bb, 'qc');

%% 4.1. Diagnostic Plot
% check QCed spectrums AC or BB sensors
% {'raw','bin','qc','prod'}
ila.DiagnosticPlot('BB',{'qc'}); % AC or BB

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
ila.DiagnosticPlot('BB',{'qc'}, false, {'diw','beta'});

%% 5. Write QC
ila.Write('qc')
ila.CheckDataStatus();

%% 6. BIN DI
% ila.cfg.days2run = ila.cfg.days2run(1)+1:ila.cfg.days2run(end)-1;
ila.BinDI();

%% 6.1. Diagnostic Plot
% check binned spectrums AC or BB sensors
ila.DiagnosticPlot('BB',{'bin'}); % AC or BB

%% 7. Write bin DI
ila.Write('bin')
ila.CheckDataStatus();

%% 8. Calibrate
ila.cfg.calibrate.BB3.filt_method = '25percentil'; % 25percentil exponential_fit
ila.cfg.calibrate.skip = {'FLOW', 'TSG'};
ila.Calibrate();
ila.CheckDataStatus();

%% 8.1. Normal and DI prod QC plots
% save_figures = true;

%%% ACS 3D %%%
ila.DiagnosticPlot('BB',{'prod'}); % AC or BB

%%% BB 3D %%%
% ila.DiagnosticPlot('BB',{'prod'}); % AC or BB

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
ila.cfg.qc.specific.run = ila.cfg.instruments2run(~contains(ila.cfg.instruments2run, 'FLOW')); % 'FLOW','ACS57','TSG', 'BB31502', 'WSCD859','PAR'
ila.QC();

%% 10. Save products
ila.Write('prod')

% % Notify with a song that the job is done
% notif_sound = load('gong'); sound(notif_sound.y, notif_sound.Fs); % handel
% return

%% re-write last version of 'qc' and 'bin'
ila.Write('bin')
ila.Write('qc')











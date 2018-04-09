% In-Line Analysis main script
% run script entirely or step by step 
% author: nils
% created: Oct 12, 2017

% Load InLineAnalysis
ila = InLineAnalysis('cfg/NAAMES3_cfg.json');

%% 1. Import | Load raw data
ila.Read();

%% 2. Synchronise instruments
% Independent of flow rate (for now)
% If flow rate varies then a new method needs to be implemented
% Play with delay of add for synchronisation
% TSG is assumed to be set at zero
% ila.instrument.FTH.Sync(30);
% ila.instrument.ACS.Sync(30+67);
% ila.instrument.BB3.Sync(30+12);
% ila.instrument.LISST.Sync(30+18);
% ila.instrument.WSCD.Sync(30+18);
% Quick visualization
% visSync(ila.instrument.FTH.data, ila.instrument.ACS.data.dt, ila.instrument.ACS.data.a(:,40), 'a (m^{-1})');
% visSync(ila.instrument.('FTH').data, ila.instrument.('BB3').data.dt, ila.instrument.('BB3').data.beta(:,1), '\beta (m^{-1} sr^{-1})');
% xlim([datenum(2017,09,14,1,30,0) datenum(2017,09,14,3,30,0)]);

% Once settings are good set them in the configuration file.
% The software is now doing the same with one line of code.
ila.Sync()

%% 3. Split filtered and tswal periods
% 3,1 Manually QC the reference
ila.QCRef();

% 3.2 Split and remove buffer periods
% Same as previous step, let's find the right settings first.
% ila.instrument.BB3.Split(ila.instrument.FTH.data, [400, 200])
% Quick visualization
% i=1; visSplit(ila.instrument.FTH.data,...
%               ila.instrument.BB3.raw.tsw.dt, ila.instrument.BB3.raw.tsw.beta(:,i),...
%               ila.instrument.BB3.raw.fsw.dt, ila.instrument.BB3.raw.fsw.beta(:,i),...
%               ila.instrument.BB3.raw.bad.dt, ila.instrument.BB3.raw.bad.beta(:,i),...
%               '\beta (counts)');
% i=20; visSplit(ila.instrument.FTH.data,...
%               ila.instrument.ACS.raw.tsw.dt, ila.instrument.ACS.raw.tsw.a(:,i),...
%               ila.instrument.ACS.raw.fsw.dt, ila.instrument.ACS.raw.fsw.a(:,i),...
%               ila.instrument.ACS.raw.bad.dt, ila.instrument.ACS.raw.bad.a(:,i),...
%               'a (m^{-1})'); ylim([-0.1 1]);

% All set ? Copy your settings in the cfg file
% Comment the previous lines and run this single line of code instead
%   need to upadte configuration of object)
ila.Split()

%% 4. Bin data
% No need to tunne the settings here, you know what you want right, it's
% set in the configuration file so let's just run it !
% ...
% You can take a break, your computer is working for you.
% This step takes a few minutes for large dataset. (It's the longest of the process)
ila.Bin()

% Let's check that it work, never trust black boxes if you take this code as it
% fh=visFlag(ila.instrument.BB3.raw.tsw, ila.instrument.BB3.raw.fsw,...
%         ila.instrument.BB3.bin.tsw, ila.instrument.BB3.bin.tsw([],:),...
%         ila.instrument.BB3.bin.fsw, ila.instrument.BB3.bin.fsw([],:),...
%         'beta', 1);

%% 5. Flag data
% 5.1 Visualize automatic flags
params.maximum_fudge_factor = 4;
params.variance_fudge_factor = 3;
params.avg_sensitivity = 1;
params.unc1_sensitivity = 1;
params.unc2_sensitivity = 2;
params.smooth_threshold = 60;
params.abs_uncertainty = 2; % 2;      
params.rel_uncertainty = 0.025; % 0.02
params_fsw = params;
params_fsw.smooth_threshold = 2;
visFlagParams(params, ila.instrument.BB3.bin.tsw, 'beta', 1);
% 5.2 Flag
ila.Flag()


%% 6. Quality Check data

ila.QC()


%% 6. Calibrate, Correct, Adjust, and Compute Products
% Load ACS DI calibrations 
acs_di = load([ila.instrument.ACS.path.raw '../NAAMES03_ACS_DIW_cals.mat']);
acs_di = table(acs_di.data(:,1), acs_di.data(:,2:89), acs_di.data(:, 90:177),...
               acs_di.data(:,178:265), acs_di.data(:,266:353),...
               'VariableNames', {'dt', 'a', 'c', 'a_avg_sd', 'c_avg_sd'});
acs_di(median(acs_di.a_avg_sd,2) > ila.cfg.flag.ACS.tot.abs_uncertainty,:) = []; % Remove bad calibrations
ila.instrument.acs.qc.diw = sortrows(acs_di);
% Load BB3 DI calibrations
ila.instrument.acs.qc.diw = table(datenum(2017,09,20,12,17,00), [84.0 87.0 69.0], [1.9 3.4 3.1], 'VariableNames', {'dt','beta', 'beta_avg_sd'});
% Process
% TODO give access to other object when calibrate
%   ex: ACS needs WSCD and FTH
%       BB3 needs TSG
ila.Calibrate()

%% 7. Preapare SeaBASS files
% ila.WriteSeaBASS()

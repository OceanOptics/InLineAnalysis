% Process NAAMES 1 - BB3
% In-Line Analysis main script
% run script entirely or step by step 
% Nils Haentjens
% created: Feb 19, 2018

% Load InLineAnalysis
ila = InLineAnalysis('cfg/NAAMES1_cfg.json');

%% 1. Import | Load raw data
ila.Read();

%% 2. Synchronise instruments
% Set delay of FTH (assume same delay as for NAAMES 3)
% ila.instrument.FTH.Sync(30);

% Synchronize DH4 (~5 hours delay)
% ila.instrument.BB3.Sync(-18640);
% ila.instrument.WSCD.Sync(-18600);

% Stretch DH4 data by few minutes (DH4 clock drift in time)
% ila.instrument.BB3.Stretch(130);
% ila.instrument.WSCD.Stretch(130);

% visSync(ila.instrument.FTH.data, ila.instrument.BB3.data.dt, ila.instrument.BB3.data.beta(:,1), '\beta (counts)');
% visSync(ila.instrument.FTH.data, ila.instrument.WSCD.data.dt, ila.instrument.WSCD.data.fdom(:,1), 'CDOM (counts)');
% xlim([datenum(2015,11,25,05,00,00) datenum(2015,11,25,16,00,00)]); % Time where switch system was off
% xlim([datenum(2015,11,7,07,00,00) datenum(2015,11,7,9,00,00)]); % Beginning sync
% xlim([datenum(2015,11,29,22,00,00) datenum(2015,11,29,23,30,00)]); % End sync

% Once settings are good set them in the configuration file.
% The software is now doing the same with one line of code
ila.Sync()
ila.Stretch()

%% 3. Split filtered and total periods
% 3.1 Manually QC the reference
% ila.cfg.qcref.mode='ui'; % if uncomment this line it shows the user interface
% (otherwise load data previously QC)
ila.QCRef();
% 3.2 Split and remove buffer periods
ila.Split();
% Quick visualization
% i=1; visSplit(ila.instrument.FTH.data,...
%               ila.instrument.BB3.raw.tsw.dt, ila.instrument.BB3.raw.tsw.beta(:,i),...
%               ila.instrument.BB3.raw.fsw.dt, ila.instrument.BB3.raw.fsw.beta(:,i),...
%               ila.instrument.BB3.raw.bad.dt, ila.instrument.BB3.raw.bad.beta(:,i),...
%               '\beta (counts)');

%% 4. Bin data
ila.Bin();
            
%% 5. Flag data
% Wow, this section is the opposite of the previous one, there is so many
% parameters to tune, make sure you understand what each parameter does.
% The best way to understand how it works is to look at the code.
% 5.1 Visualize automatic flags
% params.maximum_fudge_factor = 4;
% params.variance_fudge_factor = 3;
% params.avg_sensitivity = 1;
% params.unc1_sensitivity = 1;
% params.unc2_sensitivity = 2;
% params.smooth_threshold = 60;
% params.abs_uncertainty = 0.00001; % BB3
% params.rel_uncertainty = 0.02; % 0.025 % BB3
% % params.abs_uncertainty = 0.0007; % WSCD
% % params.rel_uncertainty = 0; % WSCD
% 
% params_fsw = params;
% params_fsw.smooth_threshold = 2;
% visFlagParams(params, ila.instrument.BB3.bin.tsw, 'beta', 1);
% ylim([0 10^-4]);
%% 5.2 Flag
% Now that you're happy with your tunning let's run apply the flagging
% ila.cfg.flag.ACS.filt.abs_uncertainty = 0.01; % 736829:736832
ila.cfg.flag.BB3.tot.abs_uncertainty = 0.0001;
ila.cfg.flag.BB3.filt.abs_uncertainty = 0.0001;
ila.Flag()

%% 6. Quality Check data
% This is either interactive or loading your previous work.
% Filter out everything that looks bad or suspect.
ila.cfg.qc.mode='load';
ila.QC();

%% 7. Calibrate, Correct, Adjust, and Compute Products
% The fun starts now ! If you did a good job up to now, this section will run
% smoothly and give you the desired products.
% 6.1 Process
ila.Calibrate()

%% Check
fig(2);
plot(ila.instrument.BB3.prod.p.dt, ila.instrument.BB3.prod.p.bbp(:,1));
datetick();

%% 8. Write data
ila.cfg.write.skip = ['FTH'];
ila.cfg.write.mode = "One file";
ila.Write();

%% 9. Merge data in one file
% Load data if start from here
ila = InLineAnalysis('cfg/NAAMES1_cfg.json');
ila.LoadProducts();
% Rename data for ease of use
tsg = ila.instrument.TSG.prod.a;
bb3 = ila.instrument.BB3.prod.p;
wscd = ila.instrument.WSCD.prod.a;

% Interp TSG on BB3
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], bb3.dt, 'linear');
tsg_bb3 = table(bb3.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4),...
                        bb3.betap, bb3.betap_sd, bb3.betap_n, bb3.bbp);
tsg_bb3.Properties.VariableNames = {'dt', 'lat', 'lon', 't', 's', 'betap', 'betap_sd', 'betap_n', 'bbp'};
tsg_bb3.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', '1/m/sr', '1/m/sr', 'none', '1/m'};
tsg_bb3.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.3d', '%.3d', '%d', '%.3d'};

% Interp TSG on WSCD
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], wscd.dt, 'linear');
tsg_wscd = table(wscd.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4),...
                        wscd.fdom, wscd.fdom_avg_sd, wscd.fdom_avg_n);
tsg_wscd.Properties.VariableNames = {'dt', 'lat', 'lon', 't', 's', 'fdom', 'fdom_sd', 'fdom_n'};
tsg_wscd.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'V', 'V', 'none'};
tsg_wscd.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.4f', '%.4f', '%d'};

% Save data to Matlab format
if true
  data = tsg;
  save([ila.instrument.TSG.path.prod 'NAAMES1_InLine_TSG.mat'], 'data');
  data = tsg_bb3; wl = ila.instrument.BB3.lambda;
  save([ila.instrument.BB3.path.prod 'NAAMES1_InLine_BB3.mat'], 'data', 'wl');
  data = tsg_wscd;
  save([ila.instrument.WSCD.path.prod 'NAAMES1_InLine_WSCD.mat'], 'data');
end

%% 10. Export to SeaBASS
% BB3
ila.meta.documents = 'NAAMES1_BB3_Processing.pdf';
ila.meta.calibration_files = 'NAAMES1_BB3_Processing.pdf';
exportSeaBASS([ila.instrument.BB3.path.prod 'NAAMES1_InLine_BB3.sb'], ila.meta, tsg_bb3, {string(ila.instrument.BB3.lambda), string(ila.instrument.BB3.lambda), '', string(ila.instrument.BB3.lambda)});
% WSCD
ila.meta.documents = 'NAAMES1_WSCD_Processing.pdf';
ila.meta.calibration_files = 'NAAMES1_WSCD_Processing.pdf';
exportSeaBASS([ila.instrument.WSCD.path.prod 'NAAMES1_InLine_WSCD.sb'], ila.meta, tsg_wscd, {'', '', ''});



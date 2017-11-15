% In-Line Analysis main script
% run script entirely or step by step 
% author: nils
% created: Oct 12, 2017

% Load InLineAnalysis
ila = InLineAnalysis('cfg/PEACETIME_cfg.json');

%% 1. Import | Load raw data
ila.Read();

%% 2. Synchronise instruments
% Independent of flow rate (for now)
% If flow rate varies then a new method needs to be implemented
% Play with delay of add for synchronisation
% TSG is assumed to be set at zero
% ila.instrument.FTH.Sync(60);
% ila.instrument.ACS.Sync(60+67);
% ila.instrument.BB.Sync(60+9);
% ila.instrument.FL.Sync(60+9);
% ila.instrument.CD.Sync(60+9);
% Quick visualization
% visSync(ila.instrument.FTH.data, ila.instrument.ACS.data.dt, ila.instrument.ACS.data.a(:,20), 'a (m^{-1})');
% visSync(ila.instrument.FTH.data, ila.instrument.BB.data.dt, ila.instrument.BB.data.beta(:,1), '\beta (counts)');
% visSync(ila.instrument.FTH.data, ila.instrument.CD.data.dt, ila.instrument.CD.data.fdom(:,1), 'CDOM (counts)');
% Sync with TSG
% fh = fig(30);
% yyaxis('right');
% sel = ila.instrument.TSG.data.s > 0;
% plot(ila.instrument.TSG.data.dt(sel), ila.instrument.TSG.data.t(sel), 'k');
% ylabel('Temperature (^{o}C)');
% ax = gca; ax.YAxis(2).Color = 'k';
% yyaxis('left'); hold('on');
% plot(ila.instrument.ACS.data.dt, ila.instrument.ACS.data.a(:,end) - ila.instrument.ACS.data.a(:,end-15), '.');
% ylabel('a(end) - a(end-15) (m^{-1})');
% datetick2_doy();
% return

% Once settings are good set them in the configuration file.
% The software is now doing the same with one line of code.
ila.Sync()

%% 3. Split filtered and tswal periods
% 3,1 Manually QC the reference
ila.QCRef();
% 3.2 Split and remove buffer periods
% Same as previous step, let's find the right settings first.
% ila.instrument.BB.Split(ila.instrument.FTH, [480, 240])
% Quick visualization
% i=1; visSplit(ila.instrument.FTH.data,...
%               ila.instrument.BB.raw.tsw.dt, ila.instrument.BB.raw.tsw.beta(:,i),...
%               ila.instrument.BB.raw.fsw.dt, ila.instrument.BB.raw.fsw.beta(:,i),...
%               ila.instrument.BB.raw.bad.dt, ila.instrument.BB.raw.bad.beta(:,i),...
%               '\beta (counts)');
% i=20; visSplit(ila.instrument.FTH.data,...
%               ila.instrument.ACS.raw.tsw.dt, ila.instrument.ACS.raw.tsw.a(:,i),...
%               ila.instrument.ACS.raw.fsw.dt, ila.instrument.ACS.raw.fsw.a(:,i),...
%               ila.instrument.ACS.raw.bad.dt, ila.instrument.ACS.raw.bad.a(:,i),...
%               'a (m^{-1})'); ylim([-0.1 1]);

% All set ? Copy your settings in the cfg file
% Comment the previous lines and run this single line of code instead
%   do not forget to upadte the configuration of the object
ila.Split()

%% 4. Bin data
% No need to tunne the settings here, you know what you want right, it's
% set in the configuration file so let's just run it !
% ...
% You can have a break, your computer is working for you.
% This step takes a few minutes for large dataset. (It's the longest of the process)
ila.Bin()

% Let's check that it work, never trust black boxes if you take this code as it
% fh=visFlag(ila.instrument.BB3.raw.tsw, ila.instrument.BB3.raw.fsw,...
%         ila.instrument.BB3.bin.tsw, ila.instrument.BB3.bin.tsw([],:),...
%         ila.instrument.BB3.bin.fsw, ila.instrument.BB3.bin.fsw([],:),...
%         'beta', 1);

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
% params.abs_uncertainty = 0.01; % 2;      
% params.rel_uncertainty = 0.02; % 0.02
% params_fsw = params;
% params_fsw.smooth_threshold = 2;
% visFlagParams(params, ila.instrument.ACS.bin.fsw, 'a', 1);

%% 5.2 Flag
% Now that you're happy with your tunning let's run apply the flagging
% ila.cfg.flag.ACS.filt.abs_uncertainty = 0.01; % 736829:736832
ila.Flag()

%% 6. Quality Check data
% This is either interactive or loading your previous work.
% Filter out everything that looks bad or suspect.
% ila.QC()
% ila.cfg.qc.mode='load';
ila.QC();

% 6. Calibrate, Correct, Adjust, and Compute Products
% The fun starts now ! If you did a good job up to now, this section will run
% smoothly and give you the desired products.
% 6.1 Process
ila.Calibrate()

%% 6.2 Checkout the result !
% ACS
wl = ila.instrument.ACS.lambda_ref;
% ACS = ila.instrument.ACS.prod;
ACS.p = data;
% visProd3D(wl, ACS.p.dt, ACS.p.ap_sd./ACS.p.ap_n, true);
visProd3D(wl, ACS.p.dt, ACS.p.ap, true);
zlabel('a_p (m^{-1})');
xlabel('\lambda (nm)');
ylabel(datestr(ila.cfg.days2run(1)));
% ylim([ila.cfg.days2run(1) ila.cfg.days2run(1)+1]);

%% 7. Save data
% 7.1 Write data in matlab format
ila.Write()

% 7.2 Write data in SeaBASS format
% ila.WriteSeaBASS()

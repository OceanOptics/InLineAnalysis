% Polish script is intended to merge the instrument as desired to export
% for SeaBASS or other formats
% author: Nils Haentjens
% created: Oct 25, 2017


% Load InLineAnalysis
ila = InLineAnalysis('cfg/PEACETIME_cfg.json');
ila.cfg.days2run = 736827:736856;

% 1. Load data products made with main_PEACETIME
% ACS was saved with mode One day one file
for d=ila.cfg.days2run
  ila.instrument.ACS.LoadProducts(['ACS' datestr(d, 'yyyymmdd')], d);
end
% Other instrument were saved with mode One file
ila.cfg.write.skip = ['ACS','FTH'];
ila.cfg.write.mode = 'One file';
ila.LoadProducts();
% Otherway to load data
% ila.instrument.FL.LoadProducts('FLALL', ila.cfg.days2run);
% ila.instrument.BB.LoadProducts('BBALL', ila.cfg.days2run);
% ila.instrument.CD.LoadProducts('CDALL', ila.cfg.days2run);
% ila.instrument.TSG.LoadProducts('TSGALL', ila.cfg.days2run);
% Create shortcuts to access instruments
acs = ila.instrument.ACS.prod.p;
tsg = ila.instrument.TSG.prod.a;
fl = ila.instrument.FL.prod.p; bb = ila.instrument.BB.prod.p; cd = ila.instrument.CD.prod.pd; 

% Clear NaN values
fl(isnan(fl.fchl), :) = [];
bb(isnan(bb.betap), :) = [];
cd(isnan(cd.fdom), :) = [];

% Clear FLBBCD data
% FL and CD data is acceptable before 2017,05,15,13,20,0 but spiky
% sel = fl.dt <= datenum(2017, 06,10,6,0,0);
sel = datenum(2017,05,15,13,20,0) < fl.dt & fl.dt <= datenum(2017, 06,10,6,0,0);
fl(~sel,:) = [];
sel = datenum(2017,05,15,13,20,0) < cd.dt & cd.dt <= datenum(2017, 06,10,6,0,0);
%       ~(cd.fdom > 1.44 & cd.dt < datenum(2017,05,15,13,20,0));
cd(~sel,:) = [];
sel = datenum(2017,05,15,13,20,0) < bb.dt & bb.dt <= datenum(2017, 06,10,6,0,0);
bb(~sel,:) = [];

% Post-processing Cleanning of ACS spectrums
wl=ila.instrument.ACS.lambda_ref;
sel_bad = any(acs.ap(:,wl < 430) < 0,2)...
          | any(acs.ap(:,:) < -0.0015,2)...
          | std(acs.ap(:,wl<430)')' > 6 * 10^-3;
acs(sel_bad,:) = [];
% Cut start and end with bubbles (start might be acceptable but chla values are unexpected)
sel = datenum(2017,05,15,13,20,0) < acs.dt & acs.dt <= datenum(2017, 06,10,6,0,0);
acs(~sel,:) = [];

%% 2. Merge instruments (get them ready for SeaBASS)
% 2.1 Merge TSG, FL, and BB
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], bb.dt, 'linear', 'extrap'); % extrap needed for first minute of data
fl_interp = interp1(fl.dt, [fl.fchl, fl.fchl_sd, fl.fchl_n], bb.dt, 'linear', 'extrap'); % extrap needed for first minute of data
tsg_flbb = table(bb.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4),...
                        fl_interp(:,1), fl_interp(:,2), fl_interp(:,3),...
                        bb.betap, bb.betap_sd, bb.betap_n, bb.bbp);
tsg_flbb.Properties.VariableNames = {'dt', 'lat', 'lon', 't', 's', 'fchl', 'fchl_sd', 'fchl_n', 'betap', 'betap_sd', 'betap_n', 'bbp'};
tsg_flbb.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'ug/L', 'ug/L', 'none', '1/m/sr', '1/m/sr', 'none', '1/m'};
tsg_flbb.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.4f', '%.4f', '%d', '%.3d', '%.3d', '%d', '%.3d'};

% 2.2 Merge TSG and CD
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], cd.dt, 'linear', 'extrap'); % extrap needed for first minute of data
tsg_cd = table(cd.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4),...
                        cd.fdom, cd.fdom_sd, cd.fdom_n);
tsg_cd.Properties.VariableNames = {'dt', 'lat', 'lon', 't', 's', 'fdom', 'fdom_sd', 'fdom_n'};
tsg_cd.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'ppb', 'ppb', 'none'};
tsg_cd.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.4f', '%.4f', '%d'};

% 2.3 Merge TSG and ACS
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], acs.dt, 'linear', 'extrap'); % extrap needed for first minute of data
tsg_acs = table(acs.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4),...
                      acs.ap, acs.ap_sd, acs.ap_n, acs.cp, acs.cp_sd, acs.cp_n,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'ap', 'ap_sd', 'ap_n', 'cp', 'cp_sd', 'cp_n'});
tsg_acs.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', '1/m', '1/m', 'none', '1/m', '1/m', 'none'};
tsg_acs.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.4f', '%.4f', '%d', '%.4f', '%.4f', '%d'};


%% 3. Compute ACS products
fprintf('Computing POC, chl, and gamma... ');
wl = ila.instrument.ACS.lambda_ref;
% 3.1 Derive POC
cp660 = interp1(wl,tsg_acs.cp',660,'linear')';
poc = cp660.*380;
% 3.2 Derive Chl
% Line heigh at 676 compared to 650 and 715
ap_a=interp1(wl,tsg_acs.ap',[650 676 715],'linear'); ap_a=ap_a';
line_height = (ap_a(:,2)-(39/65*ap_a(:,1)+26/65*ap_a(:,3)));
% Get chl
chl=157*line_height.^1.22;
chl(real(chl) ~= chl) = NaN;
% 3.3 Derive Gamma (size parameter)
% remove NaN value first
[cp0,gamma,fiterr] = FitSpectra_HM2(wl(:,1:end-2),tsg_acs.cp(:,1:end-2));
fprintf('Done\n');
% 3.4 Make table
tsg_acs_prod = table(tsg_acs.dt, tsg_acs.lat, tsg_acs.lon, tsg_acs.t, tsg_acs.s, poc, chl, gamma);
tsg_acs_prod.Properties.VariableNames = {'dt', 'lat', 'lon', 't', 's', 'POC', 'Chl', 'gamma'};
tsg_acs_prod.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'ug/L', 'ug/L', 'none'};
tsg_acs_prod.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.1f', '%.4f', '%.4f'};

%% 4. Save as Matlab
data = tsg;
save([ila.instrument.TSG.path.prod 'PEACETIME_InLine_TSG.mat'], 'data');
data = tsg_flbb;
save([ila.instrument.BB.path.prod 'PEACETIME_InLine_FLBB.mat'], 'data');
data = tsg_cd;
save([ila.instrument.CD.path.prod 'PEACETIME_InLine_CD.mat'], 'data');
data = tsg_acs;
save([ila.instrument.ACS.path.prod 'PEACETIME_InLine_ACS.mat'], 'data');
data = tsg_acs_prod;
save([ila.instrument.ACS.path.prod 'PEACETIME_InLine_ACS_prod.mat'], 'data');

%% 5. Export to SeaBASS
% 5.1 Set SeaBASS standards
% FLBB
tsg_flbb.fchl_se = tsg_flbb.fchl_sd ./ tsg_flbb.fchl_n;
tsg_flbb.betap_se = tsg_flbb.betap_sd ./ tsg_flbb.betap_n;
tsg_flbb(:,8) = []; tsg_flbb(:,10) = []; % delete _n
tsg_flbb = [tsg_flbb(:,1:7) tsg_flbb(:,11) tsg_flbb(:,8:9) tsg_flbb(:,12) tsg_flbb(:,10)];% reorder column
tsg_flbb.Properties.VariableNames = {'dt', 'lat', 'lon', 't', 's', 'Chl', 'Chl_sd', 'Chl_se', 'VSFp_124ang', 'VSFp_124ang_sd', 'VSFp_124ang_se', 'bbp'};
tsg_flbb.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'ug/L', 'ug/L', 'ug/L', '1/m/sr', '1/m/sr', '1/m/sr', '1/m'};
tsg_flbb.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.4f', '%.4f', '%.2d', '%.3d', '%.2d', '%.2d', '%.3d'};
% CD (issue need to convert to Volts ??)
tsg_cd.Properties.VariableNames = {'dt', 'lat', 'lon', 't', 's', 'cdmf', 'cdmf_sd', 'bincount'};
tsg_cd.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'ppb', 'ppb', 'none'};
% ACS
tsg_acs = tsg_acs(:,[1:7 9:11]);
tsg_acs.Properties.VariableNames{end} = 'bincount';
% ACS_prod

%% 5.2 Export to SeaBASS
% flbb
ila.meta.documents = 'PEACETIME_FLBBCD_Processing.pdf';
ila.meta.calibration_files = 'PEACETIME_FLBBCD_Processing.pdf';
exportSeaBASS([ila.instrument.BB.path.prod 'PEACETIME_InLine_FLBB.sb'], ila.meta, tsg_flbb, {'', '', '', string(ila.instrument.BB.lambda), string(ila.instrument.BB.lambda), string(ila.instrument.BB.lambda), ''});
% cd
ila.meta.documents = 'PEACETIME_FLBBCD_Processing.pdf';
ila.meta.calibration_files = 'PEACETIME_FLBBCD_Processing.pdf';
exportSeaBASS([ila.instrument.CD.path.prod 'PEACETIME_InLine_CD.sb'], ila.meta, tsg_cd, {'', '', ''});
% acs
ila.meta.documents = 'PEACETIME_ACS_Processing.pdf';
ila.meta.calibration_files = 'acs111.dev';
exportSeaBASS([ila.instrument.ACS.path.prod 'PEACETIME_InLine_ACS.sb'], ila.meta, tsg_acs, {string(ila.instrument.ACS.lambda_ref), string(ila.instrument.ACS.lambda_ref), string(ila.instrument.ACS.lambda_ref), string(ila.instrument.ACS.lambda_ref), ''});
% acs_prod
ila.meta.documents = 'PEACETIME_ACS_Products.pdf';
ila.meta.calibration_files = 'acs111.dev';
exportSeaBASS([ila.instrument.ACS.path.prod 'PEACETIME_InLine_ACS_prod.sb'], ila.meta, tsg_acs_prod, {'', '', ''});
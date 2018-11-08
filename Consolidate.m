% Consolidate daily processing, merge variables, and export to SeaBASS
% author: Nils
% created: Aug 17, 2018

%% Import data
ila = InLineAnalysis('cfg/default_cfg.m');
% ACS 298: Aug 11 to Aug 20
% ACS 301: Aug 20 to Sept ...
% WSCD 859: Aug 11 to 17
% WSCD 1299: Aug 17 to ...

% Update cfg
% ila.cfg.days2run = datenum(2018,08,11):datenum(2018,09,12);
% ila.cfg.instruments2run = {'TSG', 'ACS', 'ACS298', 'BB3', 'LISST', 'WSCD', 'WSCD859', 'ALFA'};
ila.cfg.instruments2run = {'TSG', 'BB3'};
ila.cfg.write.mode = 'One day one file';
% Load products
ila.instrument.ACS298.ReadDeviceFile()
ila.Read('prod');
% Merge Products coming from multiple instruments
% ila.MergeProducts('WSCD', 'WSCD859'); % Works only for instruments having variable of similar width (not the ACS)

%% Merge (GPS+TSG+<products>)
fprintf('Consolidate... ');
tsg = ila.instrument.TSG.prod.a;
% tsg = table(tsg.dt, tsg.lat, tsg.lon, tsg.t, tsg.s, tsg.o2, tsg.o2_sat, tsg.fchl, tsg.par,...
%             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'o2', 'o2_sat', 'fchl', 'par'});
% tsg.Properties.VariableUnits = {'', 'degN', 'degE', 'degC', 'PSU', 'mL/L', 'mL/L', 'ug/L', 'uE/s/m^2'};
% tsg.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f', '%.3f', '%.2f'};
tsg = table(tsg.dt, tsg.lat, tsg.lon, tsg.t, tsg.s, tsg.fchl,...
            'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'fchl'});
tsg.Properties.VariableUnits = {'', 'degN', 'degE', 'degC', 'PSU', 'ug/L'};
tsg.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f','%.3f'};
% ACS Particulate
ila_acs = ila.instrument.ACS298.prod.p;
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], ila_acs.dt, 'linear', 'extrap'); % extrap needed for first minute of data
acs = table(ila_acs.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4),...
                      ila_acs.ap, ila_acs.ap_sd, ila_acs.ap_n, ila_acs.cp, ila_acs.cp_sd, ila_acs.cp_n,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'ap', 'ap_sd', 'ap_n', 'cp', 'cp_sd', 'cp_n'});
acs.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', '1/m', '1/m', 'none', '1/m', '1/m', 'none'};
acs.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.4f', '%.4f', '%d', '%.4f', '%.4f', '%d'};
% ACS Dissolved
ila_acsg = ila.instrument.ACS298.prod.g;
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], ila_acsg.dt, 'linear', 'extrap'); % extrap needed for first minute of data
acs_g = table(ila_acsg.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4),...
                      ila_acsg.ag, ila_acsg.ag_sd, ila_acsg.ag_n, ila_acsg.cg, ila_acsg.cg_sd, ila_acsg.cg_n,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'ag', 'ag_sd', 'ag_n', 'cg', 'cg_sd', 'cg_n'});
acs_g.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', '1/m', '1/m', 'none', '1/m', '1/m', 'none'};
acs_g.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.4f', '%.4f', '%d', '%.4f', '%.4f', '%d'};
% ACS 298
% ila_acs298 = ila.instrument.ACS298.prod.p;
% tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], ila_acs298.dt, 'linear', 'extrap'); % extrap needed for first minute of data
% acs298 = table(ila_acs298.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4),...
%                       ila_acs298.ap, ila_acs298.ap_sd, ila_acs298.ap_n, ila_acs298.cp, ila_acs298.cp_sd, ila_acs298.cp_n,...
%              'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'ap', 'ap_sd', 'ap_n', 'cp', 'cp_sd', 'cp_n'});
% acs298.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', '1/m', '1/m', 'none', '1/m', '1/m', 'none'};
% acs298.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.4f', '%.4f', '%d', '%.4f', '%.4f', '%d'};
% ACS merged prod
% ila_acsm = [ila_acs298.dt, ila_acs298.chl, ila_acs298.poc, ila_acs298.gamma;...
%            ila_acs.dt, ila_acs.chl, ila_acs.poc, ila_acs.gamma];
% tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], ila_acsm(:,1), 'linear', 'extrap'); % extrap needed for first minute of data
% acs_prod = table(ila_acsm(:,1), tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4),...
%                       ila_acsm(:,2), ila_acsm(:,3), ila_acsm(:,4),...
%              'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'Chl', 'POC', 'gamma'});
% acs_prod.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'ug/L', 'ug/L', 'unitless'};
% acs_prod.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.4f', '%.2f', '%.2f'};
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], ila_acs.dt, 'linear', 'extrap'); % extrap needed for first minute of data
acs_prod = table(ila_acs.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4),...
                      ila_acs.chl, ila_acs.poc, ila_acs.gamma, ila_acs.chl_naames,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'Chl', 'POC', 'gamma', 'Chl_NAAMES'});
acs_prod.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'ug/L', 'ug/L', 'unitless', 'ug/L'};
acs_prod.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.4f', '%.2f', '%.2f', '%.4f'};
% BB3
ila_bb3 = ila.instrument.BB3.prod.p;
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], ila_bb3.dt, 'linear', 'extrap'); % extrap needed for first minute of data
bb3 = table(ila_bb3.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4), ila_bb3.betap, ila_bb3.betap_sd, ila_bb3.bbp, ila_bb3.betap_n,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'VSF_124ang', 'VSF_124ang_sd','bbp', 'bincount'});
bb3.Properties.VariableUnits = {'', 'degN', 'degE', 'degC', 'PSU', '1/m/sr', '1/m/sr', '1/m', 'none'};
bb3.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.3d', '%.3d', '%.3d', '%d'};
% LISST particulate
lisst = ila.instrument.LISST.prod.p;
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], lisst.dt, 'linear', 'extrap'); % extrap needed for first minute of data
lisst = table(lisst.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4), lisst.betap, lisst.betap_sd, lisst.cp, lisst.VD, lisst.VD_sd, lisst.PSD, lisst.VSD,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'betap', 'betap_sd', 'cp', 'VD', 'VD_sd', 'PSD', 'VSD'});
lisst.Properties.VariableUnits = {'', 'degN', 'degE', 'degC', 'PSU', 'counts', 'counts', '1/m', 'uL/L', 'uL/L', '#/um^3/um', '#/um'};
lisst.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f'};
% WSCD particulate
wscd = ila.instrument.WSCD.prod.a;
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s, tsg.fchl], wscd.dt, 'linear', 'extrap'); % extrap needed for first minute of data
fl = table(wscd.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4), tsg_interp(:,5), wscd.fdom, wscd.fdom_avg_sd, wscd.fdom_avg_n,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'fchl', 'cdmf', 'cdmf_sd', 'bincount'});
fl.Properties.VariableUnits = {'', 'degN', 'degE', 'degC', 'PSU', 'counts', 'counts', 'counts', 'none'};
fl.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.1f', '%.1f', '%.1f', '%d'};
% ALFA particulate
alfa = ila.instrument.ALFA.prod.a;
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s, tsg.fchl, tsg.par], alfa.dt, 'linear', 'extrap'); % extrap needed for first minute of data
alfa = table(alfa.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4), tsg_interp(:,5), tsg_interp(:,6),...
             alfa.Chlb, alfa.CFRb, alfa.WLCFb, alfa.CDOMRb,...
             alfa.R613Rb, alfa.R625Rb, alfa.R642Rb, alfa.R662Rb,...
             alfa.Chlg, alfa.CFRg, alfa.WLCFg, alfa.PE1Rg, alfa.PE2Rg, alfa.PE3Rg,...
             alfa.R642Rg, alfa.R662Rg, alfa.PE1CFg, alfa.PE2CFg, alfa.PE3CFg,...
             alfa.PE12Rg, alfa.PE12CFg, alfa.WLPE12g,...
             alfa.FvFm, alfa.FvFmC, alfa.FvFmG, alfa.FvFmCG, alfa.Chlb_avg_n,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'fchl_ws3s', 'par',...
                               'Chlb', 'CFRb', 'WLCFb', 'CDOMRb', 'R613Rb',...
                               'R625Rb', 'R642Rb', 'R662Rb', 'Chlg', 'CFRg', 'WLCFg',...
                               'PE1Rg', 'PE2Rg', 'PE3Rg', 'R642Rg', 'R662Rg', 'PE1CFg',...
                               'PE2CFg', 'PE3CFg', 'PE12Rg', 'PE12CFg', 'WLPE12g',...
                               'FvFm', 'FvFmC', 'FvFmG', 'FvFmCG', 'bincount'});
alfa.Properties.VariableUnits = {'', 'degN', 'degE', 'degC', 'PSU', 'ug/L', 'uE/s/m^2', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', 'none'};
alfa.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.3f', '%.2f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%d'};
fprintf('Done\n');

%% Export to matlab format
fprintf('Export to mat... ');
save([ila.instrument.TSG.path.prod 'EXPORTS_TSG'], 'tsg');
wl = ila.instrument.ACS.lambda_ref;
save([ila.instrument.ACS.path.prod 'EXPORTS_ACS301_Particulate'], 'acs', 'wl');
wl = ila.instrument.ACS298.lambda_ref; acs301 = acs; acs = acs298;
save([ila.instrument.ACS298.path.prod 'NAAMES4_InLine_ACS_Particulate'], 'acs', 'wl');
save([ila.instrument.ACS298.path.prod 'NAAMES4_InLine_ACS_Dissolved'], 'acs_g', 'wl');
save([ila.instrument.ACS298.path.prod 'NAAMES4_InLine_ACS_Products'], 'acs_prod');
wl = ila.instrument.BB3.lambda;
save([ila.instrument.BB3.path.prod 'EXPORTS_BB3_Particulate'], 'bb3', 'wl');
diameter = ila.instrument.LISST.diameters;
save([ila.instrument.LISST.path.prod 'EXPORTS_LISST_Particulate'], 'lisst', 'diameter');
save([ila.instrument.WSCD.path.prod 'EXPORTS_WSCD'], 'fl');
save([ila.instrument.ALFA.path.prod 'EXPORTS_ALFA'], 'alfa');
fprintf('Done\n');

clear('acs'); % remove variable acs from here to avoid confusion
%% Export to SeaBASS format
% Reformat
% ACS (can only keep 1 bincount)
acs = acs(:,[1:7 9:11]);
acs.Properties.VariableNames{end} = 'bincount';
acs_g = acs_g(:,[1:7 9:11]);
acs_g.Properties.VariableNames{end} = 'bincount';

% Keep only chl specific (NAAMES) instead of chl_a global
acs_prod.Chl = acs_prod.Chl_NAAMES;
acs_prod.Chl_NAAMES = [];


acs298 = acs298(:,[1:7 9:11]);
acs298.Properties.VariableNames{end} = 'bincount';
%% Export
fprintf('Export to sb... ');
wl = ila.instrument.ACS298.lambda_ref;
% acs Particulate
ila.meta.documents = 'NAAMES4_ACS_Processing.pdf';
ila.meta.calibration_files = 'acs298.dev';
exportSeaBASS([ila.instrument.ACS298.path.prod 'NAAMES4_InLine_ACS_Particulate.sb'], ila.meta, acs, {string(wl), string(wl), string(wl), string(wl), ''});
% acs Dissovled
ila.meta.documents = 'NAAMES4_ACS_Processing.pdf';
ila.meta.calibration_files = 'acs298.dev';
exportSeaBASS([ila.instrument.ACS298.path.prod 'NAAMES4_InLine_ACS_Dissolved.sb'], ila.meta, acs_g, {string(wl), string(wl), string(wl), string(wl), ''});
%
% ila.meta.documents = 'EXPORTS_ACS298_Processing.pdf';
% ila.meta.calibration_files = 'acs298.dev';
% exportSeaBASS([ila.instrument.ACS298.path.prod '../EXPORTS_ACS298.sb'], ila.meta, acs298, {string(ila.instrument.ACS298.lambda_ref), string(ila.instrument.ACS298.lambda_ref), string(ila.instrument.ACS298.lambda_ref), string(ila.instrument.ACS298.lambda_ref), ''});
% acs_prod
ila.meta.documents = 'NAAMES4_ACS_Processing.pdf';
ila.meta.calibration_files = 'acs298.dev';
exportSeaBASS([ila.instrument.ACS298.path.prod 'NAAMES4_InLine_ACS_Products.sb'], ila.meta, acs_prod, {'', '', ''});
%% bb3
ila.meta.documents = 'NAAMES4_BB3_Processing.pdf';
ila.meta.calibration_files = 'NAAMES4_BB3_Processing.pdf';
exportSeaBASS([ila.instrument.BB3.path.prod 'NAAMES4_InLine_BB3.sb'], ila.meta, bb3, {string(wl), string(wl), string(wl), ''});
% wscd
ila.meta.documents = 'EXPORTS_WSCD_Processing.pdf';
ila.meta.calibration_files = 'EXPORTS_WSCD_Processing.pdf';
exportSeaBASS([ila.instrument.WSCD.path.prod '../EXPORTS_WSCD.sb'], ila.meta, fl, {'', '', '', ''});
fprintf('Done\n');
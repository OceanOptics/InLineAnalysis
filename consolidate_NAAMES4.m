% NAAMES IV Consolidate files
% author: Nils
% created: March 31, 2018

% Load InLineAnalysis
ila = InLineAnalysis('cfg/NAAMES4_cfg.json');
% Update cfg
ila.cfg.days2run = datenum(2018,03,20):datenum(2018,04,13);
% ila.cfg.days2run = datenum(2018,04,2);
ila.cfg.instruments2run = ["TSG", "ACS", "BB3", "LISST", "WSCD"];

% Load processed data
% ila.cfg.write.mode = 'One file';
ila.cfg.write.mode = 'One day one file';
ila.LoadProducts();

%% Consolidate GPS+TSG+Other instrument
fprintf('Consolidate... ');
tsg = ila.instrument.TSG.prod.a;
% ACS
ila_acs = ila.instrument.ACS.prod.p;
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], ila_acs.dt, 'linear', 'extrap'); % extrap needed for first minute of data
acs = table(ila_acs.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4),...
                      ila_acs.ap, ila_acs.ap_sd, ila_acs.ap_n, ila_acs.cp, ila_acs.cp_sd, ila_acs.cp_n,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'ap', 'ap_sd', 'ap_n', 'cp', 'cp_sd', 'cp_n'});
acs.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', '1/m', '1/m', 'none', '1/m', '1/m', 'none'};
acs.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.4f', '%.4f', '%d', '%.4f', '%.4f', '%d'};
% ACS prod
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], ila_acs.dt, 'linear', 'extrap'); % extrap needed for first minute of data
acs_prod = table(ila_acs.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4),...
                      ila_acs.chl, ila_acs.poc, ila_acs.gamma,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'chl', 'poc', 'gamma'});
acs_prod.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'ug/L', 'ug/L', 'unitless'};
acs_prod.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.4f', '%.2f', '%.2f'};
% BB3
ila_bb3 = ila.instrument.BB3.prod.p;
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], ila_bb3.dt, 'linear', 'extrap'); % extrap needed for first minute of data
bb3 = table(ila_bb3.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4), ila_bb3.betap, ila_bb3.betap_sd, ila_bb3.bbp, ila_bb3.betap_n,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'VSF_124ang', 'VSF_124ang_sd','bbp', 'bincount'});
bb3.Properties.VariableUnits = {'', 'degN', 'degE', 'degC', 'PSU', '1/sr/m', '1/sr/m', '1/m', 'none'};
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
fprintf('Done\n');

%% Save to matlab format
save([ila.instrument.TSG.path.prod '../NAAMES4_InLine_TSG'], 'tsg');
wl = ila.instrument.ACS.lambda_ref;
save([ila.instrument.ACS.path.prod '../NAAMES4_InLine_ACS'], 'acs', 'wl');
save([ila.instrument.ACS.path.prod '../NAAMES4_InLine_ACS_Products'], 'acs_prod');
wl = ila.instrument.BB3.lambda;
save([ila.instrument.BB3.path.prod '../NAAMES4_InLine_BB3'], 'bb3', 'wl');
diameter = ila.instrument.LISST.diameters;
save([ila.instrument.LISST.path.prod '../NAAMES4_InLine_LISST_Particulate'], 'lisst', 'diameter');
save([ila.instrument.WSCD.path.prod '../NAAMES4_InLine_WSCD_Particulate'], 'fl');

% Save to csv
% writetable(acs_prod, [ila.instrument.ACS.path.prod '../NAAMES4_InLine_ACS_Products.csv']);
return
%% Check result
sel = 1 <= diameter & diameter <= 100;
visProd2D(diameter(sel), lisst.dt, lisst.PSD(:,sel), false);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Diameter (\mum)'); ylabel('PSD (# mL^{-1} \mum^{-1})'); 
visProd3D(diameter(sel), lisst.dt, lisst.PSD(:,sel), true, 'Log'); %, 'Log'
set(gca, 'ZScale', 'log', 'XScale', 'log');
xlabel('Diameter (\mum)'); zlabel('PSD (# mL^{-1} \mum^{-1})');

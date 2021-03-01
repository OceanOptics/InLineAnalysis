% Consolidate daily processing, merge variables, and export to SeaBASS
% author: Guillaume Bourdin
% created: March 1, 2021

%% Import data
cd 'C:\Users\Gui\Documents\MATLAB\InLineAnalysis\InLineAnalysis-master'
cruise = 'TaraMicrobiome';
ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
path_dev = [strrep(ila.instrument.FTH.path.prod, ...
  ['prod' filesep], '') 'DeviceFiles'];

% create Graph folder if it doesn't exist
if ~isfolder([ila.instrument.FTH.path.prod 'Graphs'])
  mkdir([ila.instrument.FTH.path.prod 'Graphs'])
end

% whenever TSG is processed load with
load([ila.instrument.FTH.path.prod cruise '_InLine_TSG_prod.mat'])

% load lat lon vector
% load([[ila.instrument.FTH.path.prod cruise '_LatLon.mat'])
% latlon=table(nav_data.dt_utc, nav_data.lat, nav_data.lon,'VariableNames',{'dt','lat','lon'});
% % [tsg2, index] = unique(tsg.dt); 
% % tsg = tsg(index,:);
% % yi = interp1(tsg2, y(index), tsg.dt);

% Update cfg
ila.cfg.write.mode = 'One day one file';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TSG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila.cfg.instruments2run = {'TSG'}; % {'FTH', 'TSG', 'BB31502','PAR', 'WSCD1082P'}
ila.cfg.write.skip = {}; % {'FTH', 'TSG', 'BB31502','PAR', 'WSCD1082P'}
ila.cfg.days2run = datenum(2020,12,24):datenum(2021,02,05);

% populate ila.instrument
ila.Read('prod');

% extract TSG data from obj
tsg_temp = ila.instrument.(cell2mat(ila.cfg.instruments2run)).prod.a;
tsg_temp.dt = datetime(tsg_temp.dt,'ConvertFrom','datenum');

% TSG Products
tsg = table(tsg_temp.dt, tsg_temp.lat, tsg_temp.lon, tsg_temp.t, ...
    tsg_temp.t_avg_sd, tsg_temp.s, tsg_temp.s_avg_sd, tsg_temp.s_avg_n,...
             'VariableNames', {'dt', 'lat', 'lon', 'sst', 'sst_sd', 'sss', 'sss_sd', 'bincount'});
tsg.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'degreesC', 'PSU', 'PSU', 'none'};
tsg.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f'};

[~,b] = sort(tsg.dt); % sort by date
tsg = tsg(b,:);

figure()
yyaxis('left')
scatter(tsg.dt, tsg.sst, 5, 'filled');
ylabel([tsg.Properties.VariableNames{strcmp(tsg.Properties.VariableNames, 'sst')} ' (' ...
    tsg.Properties.VariableUnits{strcmp(tsg.Properties.VariableNames, 'sst')} ')'])
yyaxis('right')
scatter(tsg.dt, tsg.sss, 5, 'filled'); 
ylabel([tsg.Properties.VariableNames{strcmp(tsg.Properties.VariableNames, 'sss')} ' (' ...
    tsg.Properties.VariableUnits{strcmp(tsg.Properties.VariableNames, 'sss')} ')'])
saveGraph([ila.instrument.TSG.path.prod 'Graphs' filesep cruise '_tsg'], 'jpg')
close figure 1

% save TSG prod
fprintf('Export to mat and csv... ');
save([ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
    cruise '_InLine_TSG_prod'], 'tsg');
writetable(tsg, [ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
    cruise '_InLine_TSG_prod.csv']);
fprintf('Done\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_leg = {...
    datenum(2020,12,24):datenum(2021,02,05);...
%     datenum(2016,06,30):datenum(2016,07,14);...
    };

list_dev = {...
    [path_dev filesep 'acs057_20200129.dev'];...
%     [path_dev filesep 'acs301_20160704_20160720.dev'];...
    };

list_instru = {...
  'ACS57';...
%   'ACS301';...
  };

data_AC = struct('particulate', [], 'product', []);
acs = [];

for i=1:size(list_instru,1)
ila.cfg.instruments2run = list_instru(i);
ila.cfg.days2run = list_leg{i};

% populate ila.instrument
ila.Read('prod');

[ila.instrument.(list_instru{i}).lambda_c, ...
  ila.instrument.(list_instru{i}).lambda_a] = importACSDeviceFile(list_dev{i});
if contains(list_dev(i),{'ac9', 'AC9'})
    ila.instrument.(list_instru{i}).lambda_a = [412 440 488 510 532 555 650 676 715];
    ila.instrument.(list_instru{i}).lambda_c = [412 440 488 510 532 555 650 676 715];
end

% interpolate SST / SSS / LatLon
ref = [list_instru{i} '_' datestr(ila.cfg.days2run(1),'yyyymmdd') '_' datestr(ila.cfg.days2run(end),'yyyymmdd')];
AC = ila.instrument.(list_instru{i}).prod.p;
AC.dt = datetime(AC.dt,'ConvertFrom','datenum');
% latlon_interp = interp1(latlon.dt, [latlon.lat, latlon.lon], AC.dt, 'linear', 'extrap'); % extrap needed for first minute of data
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.sst, tsg.sss], AC.dt, 'linear', 'extrap'); % extrap needed for first minute of data

% ACS Particulate for SeaBASS
data_AC.particulate.(ref) = table(AC.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4), AC.ap, AC.ap_sd, AC.cp, AC.cp_sd, AC.cp_n,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'ap', 'ap_sd', 'cp', 'cp_sd', 'bincount'});
data_AC.particulate.(ref).Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', '1/m', '1/m', '1/m', '1/m', 'none'};
data_AC.particulate.(ref).Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%d'};

% ACS Products
data_AC.product.(ref) = table(AC.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4), AC.poc, AC.chl, AC.gamma,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'POC_cp','Chl_lineheight','cp_gamma'});
data_AC.product.(ref).Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'ug/L', 'ug/L', 'unitless'};
data_AC.product.(ref).Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f', '%.2f'};

%Remove NaN records
data_AC.product.(ref).POC_cp(data_AC.product.(ref).POC_cp < 0) = NaN;
data_AC.product.(ref).Chl_lineheight(data_AC.product.(ref).Chl_lineheight < 0) = NaN;
data_AC.product.(ref).cp_gamma(data_AC.product.(ref).cp_gamma < 0) = NaN;
data_AC.product.(ref)(all(isnan(table2array(data_AC.product.(ref)(:,6:8))),2),:)=[];

visProd3D(ila.instrument.(list_instru{i}).lambda_a, data_AC.particulate.(ref).dt, ...
  data_AC.particulate.(ref).ap, false, 'Wavelength', false, 72);
zlabel('a_p (m^{-1})'); xlabel('{\lambda} (nm)'); ylabel('Time');
saveGraph([ila.instrument.(list_instru{i}).path.prod 'Graphs' filesep cruise '_' (ref) '_particulate_ap'], 'fig')
saveGraph([ila.instrument.(list_instru{i}).path.prod 'Graphs' filesep cruise '_' (ref) '_particulate_ap'], 'jpg')
close figure 72
visProd3D(ila.instrument.(list_instru{i}).lambda_c, data_AC.particulate.(ref).dt, ...
  data_AC.particulate.(ref).cp, false, 'Wavelength', false, 73);
zlabel('c_p (m^{-1})'); xlabel('{\lambda} (nm)'); ylabel('Time');
saveGraph([ila.instrument.(list_instru{i}).path.prod 'Graphs' filesep cruise '_' (ref) '_particulate_cp'], 'fig')
saveGraph([ila.instrument.(list_instru{i}).path.prod 'Graphs' filesep cruise '_' (ref) '_particulate_cp'], 'jpg')
close figure 73

% pause(2)
% % prompt approuval
% dlg_title = 'Validation';
% num_lines = 1;
% prompt = 'Do you validate the time series? ';
% defaultans = {'ok'};
% answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
% if ~strcmp(answer(1),{'ok','yes','OK','YES','Ok','Yes'})
%     error('operation aborted by user')
% end
% close all

% export particulate to SeaBASS format
ila.meta.documents = [cruise '_ACS_ProcessingReport.pdf'];
[~, calfile] = fileparts(list_dev{i});
ila.meta.calibration_files = [calfile '.dev'];
exportSeaBASS([ila.instrument.(list_instru{i}).path.prod cruise '_InLine_' ref '_Particulate.sb'],...
    ila.meta,...
    data_AC.particulate.(ref),...
    {string(ila.instrument.(list_instru{i}).lambda_a),...
    string(ila.instrument.(list_instru{i}).lambda_a),...
    string(ila.instrument.(list_instru{i}).lambda_c),...
    string(ila.instrument.(list_instru{i}).lambda_c),''});
sprintf('%s_InLine_%s_Particulate.sb saved', cruise, ref)

% export product to SeaBASS format
exportSeaBASS([ila.instrument.(list_instru{i}).path.prod cruise '_InLine_' ref '_Product.sb'],...
    ila.meta,...
    data_AC.product.(ref),...
    {'', '', ''});
sprintf('%s_InLine_%s_Product.sb saved', cruise, ref)

% ACS merged prod
acs = [acs; data_AC.product.(ref)];
ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
end

% sort AC and save prod file .mat
[~,b] = sort(acs.dt); % sort by date
acs = acs(b,:);

acs.Properties.VariableNames = {'dt', 'lat', 'lon', 'sst', 'sss', 'POC', 'chl', 'gamma'};

% load TaraPacific_InLine_ACS_prod.mat
figure('WindowState', 'maximized')
yyaxis('left')
scatter(acs.dt, acs.POC,10,[0 0 1],'filled'); ylabel('[POC] (mg.m^{-3})'); hold on;
yyaxis('right')
scatter(acs.dt, acs.chl,10,[0 1 0],'filled');
scatter(acs.dt, acs.gamma,10,[1 0 0],'filled'); ylabel('[chl a] (mg.m^{-3}) and gamma (unitless)');
legend('[POC](mg.m^{-3})','[chl a](mg.m^{-3})','gamma (unitless)')
saveGraph([ila.instrument.(list_instru{i}).path.prod 'Graphs' filesep cruise '_ACS_prod_timeseries'], 'fig')
saveGraph([ila.instrument.(list_instru{i}).path.prod 'Graphs' filesep cruise '_ACS_prod_timeseries'], 'jpg')
close figure 1

SimpleMap(acs.POC, acs(:,1:3), '[POC](mg.m^{-3})')
saveGraph([ila.instrument.TSG.path.prod 'Graphs' filesep cruise '_ACS_POC_map'], 'jpg')
close figure 1

SimpleMap(acs.chl, acs(:,1:3), '[chl a](mg.m^{-3})')
saveGraph([ila.instrument.TSG.path.prod 'Graphs' filesep cruise '_ACS_chl_map'], 'jpg')
close figure 1

SimpleMap(acs.gamma, acs(:,1:3), 'gamma (unitless)')
saveGraph([ila.instrument.TSG.path.prod 'Graphs' filesep cruise '_ACS_gamma_map'], 'jpg')
close figure 1

% save AC prod
fprintf('Export to mat and csv... ');
save([ila.instrument.(list_instru{i}).path.prod ...
    cruise '_InLine_ACS_prod'], 'acs');
writetable(acs, [ila.instrument.(list_instru{i}).path.prod ...
    cruise '_InLine_ACS_prod.csv']);
fprintf('Done\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WSCD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'WSCD859'}; % {'FTH', 'TSG', 'BB31502','PAR', 'WSCD859'}
ila.cfg.days2run = datenum(2020,12,24):datenum(2021,02,05);

% populate ila.instrument
ila.Read('prod');

% interpolate SST / SSS / LatLon
wscd_temp = ila.instrument.(cell2mat(ila.cfg.instruments2run)).prod.pd;
wscd_temp.dt = datetime(wscd_temp.dt,'ConvertFrom','datenum');
% latlon_interp = interp1(latlon.dt, [latlon.lat, latlon.lon], fdom_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.sst, tsg.sss], wscd_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data

% FDOM Products
wscd = table(wscd_temp.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4), wscd_temp.fdom, wscd_temp.fdom_sd, wscd_temp.fdom_n,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'fdom','fdom_sd','bincount'});
wscd.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'ppb', 'ppb', 'none'};
wscd.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f', '%.2f'};

[~,b] = sort(wscd.dt); % sort by date
wscd = wscd(b,:);

figure()
scatter(wscd.dt, wscd.fdom, 5, 'filled'); ylabel('fdom ppb')
saveGraph([ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
  'Graphs' filesep cruise '_FDOM_timeseries'], 'jpg')
close figure 1

SimpleMap(wscd.fdom, wscd(:,1:3), 'WSCD fdom ppb')
saveGraph([ila.instrument.TSG.path.prod 'Graphs' filesep cruise '_WSCD_fdom_map'], 'jpg')
close figure 1

% % export product to SeaBASS format
% ila.meta.documents = [cruise '_WSCD_ProcessingReport_V2.pdf';
% ila.meta.calibration_files = cell2mat(list_dev(i));
% exportSeaBASS([ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod cruise '_InLine_' cell2mat(ila.cfg.instruments2run) '_Product.sb'],...
%     ila.meta,...
%     wscd,...
%     {'', '', ''});
% sprintf('%s_InLine_%s_Product.sb saved', cruise, cell2mat(ila.cfg.instruments2run))

wscd.Properties.VariableNames = {'dt', 'lat', 'lon', 'sst', 'sss','fdom','fdom_sd','fdom_n'};

% save WSCD prod
fprintf('Export to mat and csv... ');
save([ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
    cruise '_InLine_WSCD_prod'], 'wscd');
writetable(wscd, [ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
    cruise '_InLine_WSCD_prod.csv']);
fprintf('Done\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'PAR'}; % {'FTH', 'TSG', 'BB3','PAR', 'WSCD1082P'}
ila.cfg.days2run = datenum(2020,12,24):datenum(2021,02,05);

% populate ila.instrument
ila.Read('prod');

% interpolate SST / SSS / LatLon
par_temp = ila.instrument.(cell2mat(ila.cfg.instruments2run)).prod.a;
par_temp.dt = datetime(par_temp.dt,'ConvertFrom','datenum');
latlon_interp = interp1(latlon.dt, [latlon.lat, latlon.lon], par_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data
tsg_interp = interp1(tsg.dt, [tsg.sst, tsg.sss_adj], par_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data

% par Products in uE/cm^2/s for SeaBASS
par = table(par_temp.dt, latlon_interp(:,1), latlon_interp(:,2), tsg_interp(:,1), tsg_interp(:,2), par_temp.par./10000, par_temp.par_sd./10000, par_temp.par_n,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'par','par_sd','bincount'});
par.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'uE/cm^2/s', 'uE/cm^2/s', 'none'};
par.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f', '%.2f'};

[~,b] = sort(par.dt); % sort by date
par = par(b,:);

figure()
scatter(par.dt, par.par, 5, 'filled'); ylabel('PAR (\muE.cm^{-2}.s^{-1})')
saveGraph([ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
  'Graphs' filesep cruise '_PAR_timeseries'], 'jpg')
close figure 1

% export product to SeaBASS format
ila.meta.documents = [cruise '_PAR_ProcessingReport_V2.pdf'];
ila.meta.calibration_files = 'PAR-50168_CalSheet.pdf';
exportSeaBASS([ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod cruise '_InLine_' cell2mat(ila.cfg.instruments2run) '_Product.sb'],...
    ila.meta,...
    par,...
    {'', '', ''});
sprintf('%s_InLine_%s_Product.sb saved', cruise, cell2mat(ila.cfg.instruments2run))

par.Properties.VariableNames = {'dt', 'lat', 'lon', 'sst', 'sss', 'par', 'par_sd', 'par_n'};

% convert to uE/m^2/s for mat file
par.par = par.par.*10000;
par.par_sd = par.par_sd.*10000;
par.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'uE/m^2/s', 'uE/m^2/s', 'none'};

% save PAR prod
fprintf('Export to mat and csv... ');
save([ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
    cruise '_InLine_PAR_prod'], 'par');
writetable(par, [ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
    cruise '_InLine_PAR_prod.csv']);
fprintf('Done\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BB3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'BB31502'}; % {'FTH', 'TSG', 'BB31502','PAR', 'WSCD1082'}
ila.cfg.days2run = datenum(2020,12,24):datenum(2021,02,05);

% populate ila.instrument
ila.Read('prod');

% interpolate SST / SSS / LatLon
lambda = ila.instrument.BB31502.lambda;
bb3_temp = ila.instrument.(cell2mat(ila.cfg.instruments2run)).prod.p;
bb3_temp.dt = datetime(bb3_temp.dt,'ConvertFrom','datenum');
% latlon_interp = interp1(latlon.dt, [latlon.lat, latlon.lon], bb3_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.sst, tsg.sss], bb3_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data

% bb3 Products
bb3 = table(bb3_temp.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4), bb3_temp.betap, bb3_temp.bbp, bb3_temp.betap_sd, bb3_temp.betap_n,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'VSF_124ang','bbp','VSF_124ang_sd','bincount'});
bb3.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', '1/m/sr', '1/m', '1/m/sr', 'none'};
bb3.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f'};

[~,b] = sort(bb3.dt); % sort by date
bb3 = bb3(b,:);
bb3.bbp(bb3.bbp<0)=NaN;
bb3.VSF_124ang_sd(bb3.VSF_124ang<0)=NaN;
bb3.VSF_124ang(bb3.VSF_124ang<0)=NaN;
bb3(all(isnan(bb3.bbp),2),:)=[];

% delete red channel
% bb3.VSF_124ang(:,3) = [];
% bb3.bbp(:,3) = [];
% bb3.VSF_124ang_sd(:,3) = [];
% bb3  = bb3(any(~isnan(bb3.bbp),2),:);
% load [cruise '_InLine_BB3_particulate.mat]

% plot time series
figure()
h = errorbar(datenum(bb3.dt), bb3.VSF_124ang(:,2), bb3.VSF_124ang_sd(:,2), 'vertical',...
    'LineStyle', 'none','LineWidth', 2, 'Color', [0.9 0.9 0.9]);
% % % Set transparency level (0:1)
% % alpha = 0.15;   
% % % Set transparency (undocumented)
% % set([h.Bar, h.Line], 'LineWidth', 5, 'ColorType', 'truecoloralpha', ...
% %     'ColorData', [h.Line.ColorData(1:3); 255*alpha])
% % set(h.Bar, 'ColorType', 'truecoloralpha', 'ColorData', [h.Bar.ColorData(1:3); 255*alpha])
h.CapSize = 0;
hold on
scatter(datenum(bb3.dt), bb3.VSF_124ang(:,2), 10, [0.3 0.8 0.3], 'filled');
ylabel('VSF_124ang (532 nm) [m^-^1]')
datetick2_doy()
saveGraph([ila.instrument.TSG.path.prod 'Graphs' filesep cruise '_VSF_124ang_timeseries'], 'jpg')
close figure 1

SimpleMap(bb3.bbp(:,2), bb3(:,1:3), 'bbp (532 nm) [m^-^1]')
saveGraph([ila.instrument.TSG.path.prod 'Graphs' filesep cruise '_bb3_bbp_map'], 'jpg')
close figure 1

% export product to SeaBASS format
ila.meta.documents = [cruise '_BB3_ProcessingReport_V2.pdf'];
ila.meta.calibration_files = 'BB3-1502_(470-532-650nm)_CharSheet.pdf';
exportSeaBASS([ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod cruise '_InLine_' cell2mat(ila.cfg.instruments2run) '_Particulate.sb'],...
    ila.meta,...
    bb3,...
    {string(lambda(1:2)),string(lambda(1:2)),string(lambda(1:2)),''});
sprintf('%s_InLine_%s_Particulate.sb saved', cruise, cell2mat(ila.cfg.instruments2run))

bb3.Properties.VariableNames = {'dt', 'lat', 'lon', 'sst', 'sss', 'VSF124', 'bbp', 'VSF124_sd','bincount'};
bb3.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'm^-1.sr^-1', 'm^-1', 'm^-1.sr^-1', 'none'};
bb3.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f'};

[bb3.poc, ~, ~, bb3.cphyto, ~, ~,] = estimate_poc_cphyto(bb3.bbp,[470 532 650], 'soccom');

% plot compare POC from different lambda
figure()
scatter(bb3.dt, bb3.poc(:,1), 5, [0.3 0.3 0.8], 'filled'); hold on
scatter(bb3.dt, bb3.poc(:,2), 5, [0.3 0.8 0.3], 'filled');
scatter(bb3.dt, bb3.poc(:,3), 5, [0.8 0.3 0.3], 'filled'); ylabel('POC (mg.m^{-3})')
saveGraph([ila.instrument.TSG.path.prod 'Graphs' filesep cruise '_POC_timeseries'], 'jpg')
close figure 1

figure()
scatter(bb3.dt, bb3.cphyto(:,1), 5, [0.3 0.3 0.8], 'filled'); hold on
scatter(bb3.dt, bb3.cphyto(:,2), 5, [0.3 0.8 0.3], 'filled');
scatter(bb3.dt, bb3.cphyto(:,3), 5, [0.8 0.3 0.3], 'filled'); ylabel('Cphyto (mg.m^{-3})')
saveGraph([ila.instrument.TSG.path.prod 'Graphs' filesep cruise '_CPhyto_timeseries'], 'jpg')
close figure 1

% save BB3 prod
fprintf('Export to mat and csv... ');
save([ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
    cruise '_InLine_BB3_particulate'], 'bb3', 'lambda');
writetable(bb3,[ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
    cruise '_InLine_BB3_particulate.csv']);
fprintf('Done\n');


%% Merge (GPS+TSG+<products>)
% fprintf('Consolidate... ');

% ACS Dissolved
% ila_acsg = ila.instrument.ACS298.prod.g;
% tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], ila_acsg.dt, 'linear', 'extrap'); % extrap needed for first minute of data
% acs_g = table(ila_acsg.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4),...
%                       ila_acsg.ag, ila_acsg.ag_sd, ila_acsg.ag_n, ila_acsg.cg, ila_acsg.cg_sd, ila_acsg.cg_n,...
%              'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'ag', 'ag_sd', 'ag_n', 'cg', 'cg_sd', 'cg_n'});
% acs_g.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', '1/m', '1/m', 'none', '1/m', '1/m', 'none'};
% acs_g.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.4f', '%.4f', '%d', '%.4f', '%.4f', '%d'};

% LISST particulate
% lisst = ila.instrument.LISST.prod.p;
% tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s], lisst.dt, 'linear', 'extrap'); % extrap needed for first minute of data
% lisst = table(lisst.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4), lisst.betap, lisst.betap_sd, lisst.cp, lisst.VD, lisst.VD_sd, lisst.PSD, lisst.VSD,...
%              'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'betap', 'betap_sd', 'cp', 'VD', 'VD_sd', 'PSD', 'VSD'});
% lisst.Properties.VariableUnits = {'', 'degN', 'degE', 'degC', 'PSU', '1/m/sr', '1/m/sr', '1/m', 'uL/L', 'uL/L', '#/um^3/um', '#/um'};
% lisst.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f'};
% WSCD particulate

% % ALFA particulate
% alfa = ila.instrument.ALFA.prod.a;
% tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.t, tsg.s, tsg.fchl, tsg.par], alfa.dt, 'linear', 'extrap'); % extrap needed for first minute of data
% % All parameters from ALFA
% % alfa = table(alfa.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4), tsg_interp(:,5), tsg_interp(:,6),...
% %              alfa.Chlb, alfa.CFRb, alfa.WLCFb, alfa.CDOMRb,...
% %              alfa.R613Rb, alfa.R625Rb, alfa.R642Rb, alfa.R662Rb,...
% %              alfa.Chlg, alfa.CFRg, alfa.WLCFg, alfa.PE1Rg, alfa.PE2Rg, alfa.PE3Rg,...
% %              alfa.R642Rg, alfa.R662Rg, alfa.PE1CFg, alfa.PE2CFg, alfa.PE3CFg,...
% %              alfa.PE12Rg, alfa.PE12CFg, alfa.WLPE12g,...
% %              alfa.FvFm, alfa.FvFmC, alfa.FvFmG, alfa.FvFmCG, alfa.Chlb_avg_n,...
% %              'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'fchl_ws3s', 'par',...
% %                                'Chlb', 'CFRb', 'WLCFb', 'CDOMRb', 'R613Rb',...
% %                                'R625Rb', 'R642Rb', 'R662Rb', 'Chlg', 'CFRg', 'WLCFg',...
% %                                'PE1Rg', 'PE2Rg', 'PE3Rg', 'R642Rg', 'R662Rg', 'PE1CFg',...
% %                                'PE2CFg', 'PE3CFg', 'PE12Rg', 'PE12CFg', 'WLPE12g',...
% %                                'FvFm', 'FvFmC', 'FvFmG', 'FvFmCG', 'bincount'});
% % alfa.Properties.VariableUnits = {'', 'degN', 'degE', 'degC', 'PSU', 'ug/L', 'uE/s/m^2', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', '?', 'none'};
% % alfa.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.3f', '%.2f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%d'};
% % Same parameters as Mike
% alfa = table(alfa.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4), tsg_interp(:,5), tsg_interp(:,6) * 0.0001,... % convert PAR from uE/m^2/s to uE/cm^2/s for SeaBASS
%              alfa.Chlb, alfa.Chlg, ...
%              alfa.FvFmC, alfa.FvFmCG, alfa.Chlb_avg_n,...
%              'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'Chl_stimf_ex460', 'par',...
%                                'Chl_stimf_ex405', 'Chl_stimf_ex514',...
%                                'Fv_Fm_ex405', 'Fv_Fm_ex514',... % 'Sigma_PSII'??, 'F0''Fm',
%                                'bincount'});
% alfa.Properties.VariableUnits = {'', 'degN', 'degE', 'degC', 'PSU', 'mg/m^3', 'uE/cm^2/s', 'mg/m^3', 'mg/m^3', 'unitless', 'unitless','none'};
% alfa.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f', '%.3f', '%.2f', '%.4f', '%.4f', '%.4f', '%.4f', '%d'};
% fprintf('Done\n');

%% Additional QC on ALFA
% alfa(alfa.Fv_Fm_ex405 < 0,:) = [];

%% Additional QC for LISST
% lisst.betap(lisst.betap < 0) = NaN;

%% Export to matlab format
% fprintf('Export to mat... ');
% diameter = ila.instrument.LISST.diameters;
% theta = ila.instrument.LISST.theta;
% save([ila.instrument.LISST.path.prod 'EXPORTS_InLine_LISST_Particulate'], 'lisst', 'diameter', 'theta');
% save([ila.instrument.ALFA.path.prod 'EXPORTS_InLine_ALFA'], 'alfa');
% fprintf('Done\n');


%% LISST
% VSF Angles
% LISST_VSF_ANGLES = ila.instrument.LISST.theta;
% LISST_VSF_ANGLES_STR = strings(size(LISST_VSF_ANGLES));
% for i=1:length(LISST_VSF_ANGLES)
%   LISST_VSF_ANGLES_STR(i) = string(sprintf('_%2.3fang', LISST_VSF_ANGLES(i)));
% end
% LISST_PSD_DIAMETERS = ila.instrument.LISST.diameters;
% LISST_PSD_DIAMETERS_STR = strings(size(LISST_VSF_ANGLES));
% for i=1:length(LISST_VSF_ANGLES)
%   LISST_PSD_DIAMETERS_STR(i) = string(sprintf('_%.3fsize', LISST_PSD_DIAMETERS(i)));
% end
% lisst = removevars(lisst, {'PSD', 'VSD'});
% lisst.Properties.VariableNames{strcmp(lisst.Properties.VariableNames, 'betap')} = 'VSF670';
% lisst.Properties.VariableNames{strcmp(lisst.Properties.VariableNames, 'betap_sd')} = 'VSF670_sd';
% lisst.Properties.VariableNames{strcmp(lisst.Properties.VariableNames, 'VD')} = 'PSD';
% lisst.Properties.VariableNames{strcmp(lisst.Properties.VariableNames, 'VD_sd')} = 'PSD_sd';

%% LISST
% ila.meta.documents = 'EXPORTS_InLine_LISST_Processing.pdf';
% ila.meta.calibration_files = 'EXPORTS_InLine_LISST_Processing.pdf';
% exportSeaBASS([ila.instrument.LISST.path.prod 'EXPORTS_InLine_LISST.sb'], ila.meta, lisst, {LISST_VSF_ANGLES_STR, LISST_VSF_ANGLES_STR, '', LISST_PSD_DIAMETERS_STR, LISST_PSD_DIAMETERS_STR});
%% ALFA
% ila.meta.documents = 'EXPORTS_InLine_ALFA_Processing.pdf';
% ila.meta.calibration_files = 'EXPORTS_InLine_ALFA_Processing.pdf';
% exportSeaBASS([ila.instrument.ALFA.path.prod 'EXPORTS_InLine_ALFA.sb'], ila.meta, alfa, {'', '', '', '', '', '', ''});
% fprintf('Done\n');
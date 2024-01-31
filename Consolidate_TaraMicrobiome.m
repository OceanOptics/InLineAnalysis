% Consolidate daily processing, merge variables, and export to SeaBASS
% author: Guillaume Bourdin
% created: March 1, 2021

%% Import data
cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master/')

cruise = 'TaraMicrobiome';
% Load InLineAnalysis and the configuration
ila = InLineAnalysis(['cfg' filesep cruise '_cfg.m']);

path_dev = strrep(ila.instrument.FLOW.path.prod, ...
  'prod', 'DeviceFiles');

% create Graph folder if it doesn't exist
if ~isfolder(fullfile(ila.instrument.FLOW.path.prod, 'plots'))
  mkdir(fullfile(ila.instrument.FLOW.path.prod, 'plots'))
end

% whenever TSG is processed load with
load(fullfile(ila.instrument.FLOW.path.prod, ['TaraChile&' cruise '_InLine_TSG_20201224_20220917_Product_v20240110.mat']))

% Update cfg
ila.cfg.write.mode = 'One day one file';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'NMEA'}; % {'FLOW', 'TSG', 'BB31502','PAR', 'WSCD859'}
ila.cfg.qcref.view = 'NMEA';
ila.cfg.days2run = datenum(2021,8,18):datenum(2022,9,17);

% populate ila.instrument
ila.Read('prod');

nmea = ila.instrument.NMEA.prod.a;
nmea.dt = datetime(nmea.dt, 'ConvertFrom', 'datenum');

%% GPS Chile
load('/Volumes/Samsung_T5/Data/TaraMicrobiome/prod/TaraChile_latlon.mat')

%% GPS meteo
load('/Volumes/Samsung_T5/Data/TaraMicrobiome/prod/TaraChile_meteo_prod.mat')

%% delete data in ZEE
to_remove = [datenum(2022, 7, 3, 2, 30, 0) datenum(2022, 7, 3, 5, 56, 0);
  datenum(2022, 7, 3, 13, 37, 0) datenum(2022, 7, 3, 16, 48, 0);
  datenum(2022, 7, 3, 19, 45, 0) datenum(2022, 7, 4, 2, 25, 0)];

for i = 1:size(to_remove, 1)
  % raw
  if ~isempty(ila.instrument.(ila.cfg.qcref.view).raw.tsw)
    ila.instrument.(ila.cfg.qcref.view).raw.tsw(ila.instrument.(ila.cfg.qcref.view).raw.tsw.dt > to_remove(i, 1) & ...
      ila.instrument.(ila.cfg.qcref.view).raw.tsw.dt < to_remove(i, 2), :) = [];
  end
  % bin
  if ~isempty(ila.instrument.(ila.cfg.qcref.view).bin.tsw)
    ila.instrument.(ila.cfg.qcref.view).bin.tsw(ila.instrument.(ila.cfg.qcref.view).bin.tsw.dt > to_remove(i, 1) & ...
      ila.instrument.(ila.cfg.qcref.view).bin.tsw.dt < to_remove(i, 2), :) = [];
  end
  % qc
  if ~isempty(ila.instrument.(ila.cfg.qcref.view).qc.tsw)
    ila.instrument.(ila.cfg.qcref.view).qc.tsw(ila.instrument.(ila.cfg.qcref.view).qc.tsw.dt > to_remove(i, 1) & ...
      ila.instrument.(ila.cfg.qcref.view).qc.tsw.dt < to_remove(i, 2), :) = [];
  end
  % prod
  fname = fieldnames(ila.instrument.(ila.cfg.qcref.view).prod);
  fname = fname{~strcmp(fname, 'QCfailed') & ~strcmp(fname, 'g') & ~strcmp(fname, 'FiltStat')};
  if ~isempty(ila.instrument.(ila.cfg.qcref.view).prod.(fname))
    ila.instrument.(ila.cfg.qcref.view).prod.(fname)(ila.instrument.(ila.cfg.qcref.view).prod.(fname).dt > to_remove(i, 1) & ...
      ila.instrument.(ila.cfg.qcref.view).prod.(fname).dt < to_remove(i, 2), :) = [];
  end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TSG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'SBE38450091'}; % {'FLOW', 'SBE38450091', 'BB31502','PAR', 'WSCD859'}
ila.cfg.days2run = datenum(2020,12,24):datenum(2022,9,17);

% populate ila.instrument
ila.Read('prod');

% interpolate SST / SSS / LatLon
tsg_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.a;
tsg_temp.dt = datetime(tsg_temp.dt,'ConvertFrom','datenum');

% build TSG table
if any(strcmp(tsg_temp.Properties.VariableNames, 'lat'))
  tsg = tsg_temp;
  tsg = removevars(tsg,{'lat_avg_sd','lat_avg_n','lon_avg_sd','lon_avg_n'});
  load(fullfile(ila.instrument.FLOW.path.prod, 'TaraChile_latlon'))
else
  % merge lat, lon
  replace_consecutive_nan = 72*60; % 72h
  tsg = merge_timeseries(tsg_temp, nmea, {'lat', 'lon'});
  tsg = merge_timeseries(tsg, latlon, {'lat', 'lon'});
  tsg = merge_timeseries(tsg, meteo_ftp, {'lat', 'lon'}, '', replace_consecutive_nan);
  % Remove NaN
  if any(any(isnan([tsg.lat tsg.lon]),2))
    warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
      sum(any(isnan([tsg.lat tsg.lon]),2)), replace_consecutive_nan/60)
    tsg(any(isnan([tsg.lat tsg.lon]),2), :) = [];
  end
end
if any(strcmp(tsg.Properties.VariableNames, 'tcal'))
  tsg = removevars(tsg,{'tcal_avg_n','c_avg_n','s_avg_n'});
  tsg = renamevars(tsg,{'t_avg_n','t','t_avg_sd','s','s_avg_sd','c','c_avg_sd'},...
    {'avg_n','sst','sst_avg_sd','sss','sss_avg_sd','cond','cond_avg_sd'});
else
  if any(strcmp(tsg.Properties.VariableNames, 'sv'))
    tsg = removevars(tsg,'sv_avg_n');
  end
  tsg = removevars(tsg,{'t1_avg_n','c1_avg_n','s_avg_n'});
  tsg = renamevars(tsg,{'t2_avg_n','t1','t1_avg_sd','t2','t2_avg_sd','s','s_avg_sd','c1','c1_avg_sd'},...
    {'avg_n','tcal','tcal_avg_sd','sst','sst_avg_sd','sss','sss_avg_sd','cond','cond_avg_sd'});
end
% add units and precision
if any(strcmp(tsg.Properties.VariableNames, 'sv'))
  tsg.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'degreesC', ...
    'S/m', 'S/m', 'PSU', 'PSU', 'm/sec', 'm/sec', 'degreesC', 'degreesC', 'none'};
  tsg.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', ...
    '%.4f', '%.4f', '%.4f', '%.4f', '%.2f', '%.2f', '%.4f', '%.4f', '%.2f'};
else
  tsg.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'degreesC', ...
    'S/m', 'S/m', 'PSU', 'PSU', 'degreesC', 'degreesC', 'none'};
  tsg.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', ...
    '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f'};
end

% sort by date
tsg = sortrows(tsg, 'dt');

ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_TSG_timeseries']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_TSG_timeseries']), 'fig')
close figure 30

SimpleMap(tsg.sst, tsg(:,1:3), 'TSG SST [°C]')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_TSG_SST_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_TSG_SST_map']), 'fig')
close figure 1
SimpleMap(tsg.sss, tsg(:,1:3), 'TSG SSS [PSU]')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_TSG_SSS_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_TSG_SSS_map']), 'fig')
close figure 1

filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, 'TSG', ...
  datetime(min(tsg.dt), 'Format', 'yyyyMMdd'), datetime(max(tsg.dt), 'Format', 'yyyyMMdd'), ...
  datetime('today', 'Format', 'yyyyMMdd'));

% % export product to SeaBASS format
% ila.meta.documents = [cruise '_TSG_ProcessingReport_V2.pdf';
% ila.meta.calibration_files = cell2mat(list_dev(i));
% exportSeaBASS(fullfile(ila.instrument.FLOW.path.prod, filename),...
%     ila.meta,...
%     tsg,...
%     {'', '', ''});
% sprintf('%s_InLine_%s_Product.sb saved', cruise, cell2mat(ila.cfg.instruments2run))

% save TSG prod
fprintf('Export to mat and csv... ');
filename = strrep(filename, '.sb', '');
save(fullfile(ila.instrument.FLOW.path.prod, filename), 'tsg');
writetable(tsg, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
fprintf('Done\n');

%% T/S scatter plot for GAYOSO leg
stations = readtable('/Users/gui/Documents/Maine/mobilisation/Microbiome/Gayoso/Rosette/Station_UTC.txt');
stations.Properties.VariableNames = {'st', 'dt_start', 'dt_end', 'lat', 'lon'};

% lagrangian #1 A061-A063
% lagrangian #2 A064-A067
% create table
id_lagrang_tot = tsg.dt >= stations.dt_start(strcmp(stations.st,'A061')) & ...
  tsg.dt <= stations.dt_end(strcmp(stations.st,'A067'));
id_lagrang = id_lagrang_tot;
id_repo_tot = tsg.dt >= (stations.dt_start(strcmp(stations.st,'A059')) - days(1)) & ...
  tsg.dt <= (stations.dt_end(strcmp(stations.st,'A070')) + days(1));
id_repo = id_repo_tot;
stations.id = cell(size(stations,1), 1);
for i = 1:size(stations, 1)
  stations.id{i} = tsg.dt >= stations.dt_start(i) & tsg.dt <= stations.dt_end(i);
  id_lagrang = id_lagrang - stations.id{i};
  id_repo = id_repo - stations.id{i};
end
id_lagrang(id_lagrang < 0) = 0;
id_lagrang = logical(id_lagrang);
id_repo(id_repo < 0) = 0;
id_repo = logical(id_repo);

% Whole leg
wholeleg = stations([1:3 10:12],:);
wholeleg.st{3} = 'Lagrangian #1';
wholeleg.st{4} = 'Lagrangian #2';
wholeleg.dt_start(3) = stations.dt_start(3);
wholeleg.dt_end(3) = stations.dt_end(5);
wholeleg.dt_start(4) = stations.dt_start(6);
wholeleg.dt_end(4) = stations.dt_end(9);
wholeleg.id{3} = tsg.dt >= wholeleg.dt_start(3) & tsg.dt <= wholeleg.dt_end(3);
wholeleg.id{4} = tsg.dt >= wholeleg.dt_start(4) & tsg.dt <= wholeleg.dt_end(4);

% Plot
col = brewermap(6, 'Paired');
TS_density(tsg.sst(id_repo_tot), tsg.sss(id_repo_tot)); hold on
lg = [];
lg(1) = scatter(tsg.sss(id_repo), tsg.sst(id_repo), 40, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.8);
for i = 1:size(wholeleg, 1)
  lg(i+1) = scatter(tsg.sss(wholeleg.id{i}), tsg.sst(wholeleg.id{i}), ...
    40, col(i,:), 'filled', 'MarkerFaceAlpha', 0.5);
end
[~,hlg] = legend(lg, cat(1, {'Repositioning'}, wholeleg.st), 'FontSize', 14, 'Location', 'North' );
% remove alpha on legend
PatchInLegend = findobj(hlg, 'type', 'patch');
set(PatchInLegend, 'FaceAlpha', 1);
ylabel('Temperature [°C]', 'FontSize', 14)
xlabel('Salinity [PSU]', 'FontSize', 14)
set(gca, 'FontSize', 14)
saveGraph(fullfile(ila.instrument.FLOW.path.prod, ...
  'plots', [cruise '_gayoso_TSplot']), 'jpg')
close figure 1


% Zoom on lagrangian
% Plot
st_lagrang = stations(3:9,:);
col = brewermap(10, 'PuOr');
col(4:6,:) = [];
TS_density(tsg.sst(id_lagrang_tot), tsg.sss(id_lagrang_tot)); hold on
lg(1) = scatter(tsg.sss(id_lagrang), tsg.sst(id_lagrang), 40, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.8);
for i = 1:size(st_lagrang, 1)
  lg(i+1) = scatter(tsg.sss(st_lagrang.id{i}), tsg.sst(st_lagrang.id{i}), ...
    40, col(i,:), 'filled', 'MarkerFaceAlpha', 0.5);
end
[~,hlg] = legend(lg, cat(1, {'Repositioning'}, st_lagrang.st), 'FontSize', 14, 'Location', 'North' );
% remove alpha on legend
PatchInLegend = findobj(hlg, 'type', 'patch');
set(PatchInLegend, 'FaceAlpha', 1);
ylabel('Temperature [°C]', 'FontSize', 14)
xlabel('Salinity [PSU]', 'FontSize', 14)
set(gca, 'FontSize', 14)
saveGraph(fullfile(ila.instrument.FLOW.path.prod, ...
  'plots', [cruise '_gayoso_zoomTSplot']), 'jpg')
close figure 1

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WSCD859 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'WSCD859'}; % {'FLOW', 'TSG', 'BB31502','PAR', 'WSCD859'}
ila.cfg.days2run = datenum(2020,12,24):datenum(2021,5,9);

% populate ila.instrument
ila.Read('prod');

wscd_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.pd;
wscd_temp.dt = datetime(wscd_temp.dt,'ConvertFrom','datenum');
% merge lat, lon, sst, sss
replace_consecutive_nan = 72*60; % 72h
wscd = merge_timeseries(wscd_temp, tsg, {'lat', 'lon', 'sst', 'sss'});
wscd = merge_timeseries(wscd, nmea, {'lat', 'lon'});
wscd = merge_timeseries(wscd, latlon, {'lat', 'lon'});
wscd = merge_timeseries(wscd, meteo_ftp, {'lat', 'lon'}, '', replace_consecutive_nan);

% Remove NaN
if any(any(isnan([wscd.lat wscd.lon]),2))
  warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
    sum(any(isnan([wscd.lat wscd.lon]),2)), replace_consecutive_nan/60)
  wscd(any(isnan([wscd.lat wscd.lon]),2), :) = [];
end

% add units and precision
wscd = renamevars(wscd, {'sst','sss','fdom_n'} ,{'t','s','bincount'});
wscd.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'v_uncalibrated', 'v_uncalibrated', 'none'};
wscd.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f', '%.2f'};

% sort by date
wscd = sortrows(wscd, 'dt');

ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_FDOM_timeseries']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_FDOM_timeseries']), 'fig')
close figure 50

SimpleMap(wscd.fdom, wscd(:,1:3), 'WSCD fdom [v uncalibrated]')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_WSCD_fdom_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_WSCD_fdom_map']), 'fig')
close figure 1

filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, ila.cfg.instruments2run{:}, ...
  datetime(min(wscd.dt), 'Format', 'yyyyMMdd'), datetime(max(wscd.dt), 'Format', 'yyyyMMdd'), ...
  datetime('today', 'Format', 'yyyyMMdd'));

% % export product to SeaBASS format
% ila.meta.documents = [cruise '_WSCD_ProcessingReport_V2.pdf';
% ila.meta.calibration_files = cell2mat(list_dev(i));
% exportSeaBASS(fullfile(ila.instrument.FLOW.path.prod, filename),...
%     ila.meta,...
%     wscd,...
%     {'', '', ''});
% sprintf('%s_InLine_%s_Product.sb saved', cruise, cell2mat(ila.cfg.instruments2run))

wscd = renamevars(wscd, {'t','s','bincount'}, {'sst','sss','fdom_n'});

% save WSCD prod
fprintf('Export to mat and csv... ');
filename = strrep(filename, '.sb', '');
save(fullfile(ila.instrument.FLOW.path.prod, filename), 'wscd');
writetable(wscd, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
fprintf('Done\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUVF6244 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'SUVF6244'}; % {'FLOW', 'TSG', 'BB31502','PAR', 'WSCD859'}
ila.cfg.days2run = datenum(2021,08,15):datenum(2022,9,17);

% populate ila.instrument
ila.Read('prod');

% interpolate SST / SSS / LatLon
suvf_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.pd;
suvf_temp.dt = datetime(suvf_temp.dt,'ConvertFrom','datenum');
% merge lat, lon, sst, sss
replace_consecutive_nan = 1*60; % 1h
suvf = merge_timeseries(suvf_temp, tsg, {'lat', 'lon', 'sst', 'sss'});
suvf = merge_timeseries(suvf, nmea, {'lat', 'lon'});
suvf = merge_timeseries(suvf, latlon, {'lat', 'lon'});
suvf = merge_timeseries(suvf, meteo_ftp, {'lat', 'lon'}, '', replace_consecutive_nan);

% Remove NaN
if any(any(isnan([suvf.lat suvf.lon]),2))
  warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
    sum(any(isnan([suvf.lat suvf.lon]),2)), replace_consecutive_nan/60)
  suvf(any(isnan([suvf.lat suvf.lon]),2), :) = [];
end

% add units and precision
suvf = renamevars(suvf, {'sst','sss','fdom_n'}, {'t','s','bincount'});
suvf.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'v_uncalibrated', 'v_uncalibrated', 'none'};
suvf.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f', '%.2f'};

% sort by date
suvf = sortrows(suvf, 'dt');

ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_SUVF_fdom_timeseries']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_SUVF_fdom_timeseries']), 'fig')
close figure 50

SimpleMap(suvf.fdom, suvf(:,1:3), 'SUVF fdom [v uncalibrated]')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_SUVF_fdom_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_SUVF_fdom_map']), 'fig')
close figure 1

filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, ila.cfg.instruments2run{:}, ...
  datetime(min(suvf.dt), 'Format', 'yyyyMMdd'), datetime(max(suvf.dt), 'Format', 'yyyyMMdd'), datetime('today', 'Format', 'yyyyMMdd'));

% export product to SeaBASS format
ila.meta.documents = [cruise '_SUVF_ProcessingReport_V2.pdf'];
ila.meta.calibration_files = 'SUVF6244_CharSheet.pdf';
exportSeaBASS(fullfile(ila.instrument.FLOW.path.prod, filename),...
    ila.meta,...
    suvf,...
    {'', '', ''});
sprintf('%s_InLine_%s_Product.sb saved', cruise, cell2mat(ila.cfg.instruments2run))

suvf = renamevars(suvf, {'t','s','bincount'}, {'sst','sss','fdom_n'});

% save SUVF prod
fprintf('Export to mat and csv... ');
filename = strrep(filename, '.sb', '');
save(fullfile(ila.instrument.FLOW.path.prod, filename), 'suvf');
writetable(suvf, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
fprintf('Done\n');

%% merge both fDOM tables
wscd.fdom_sensor = repmat({'WSCD859'}, size(wscd.dt));
suvf.fdom_sensor = repmat({'SUVF6244'}, size(suvf.dt));
fdom = [wscd; suvf];

filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s', cruise, 'fdom_WSCD859_SUVF6244_merged', ...
  datetime(min(fdom.dt), 'Format', 'yyyyMMdd'), datetime(max(fdom.dt), 'Format', 'yyyyMMdd'), datetime('today', 'Format', 'yyyyMMdd'));

% save SUVF prod
fprintf('Export to mat and csv... ');
save(fullfile(ila.instrument.FLOW.path.prod, filename), 'fdom');
writetable(fdom, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
fprintf('Done\n');

load(fullfile(ila.instrument.FLOW.path.prod, 'TaraMicrobiome_InLine_fdom_WSCD859_SUVF6244_merged_20201226_20220917_Product_v20240111.mat'), 'fdom')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'PAR'}; % {'FLOW', 'TSG', 'BB3','PAR', 'WSCD1082P'}
ila.cfg.days2run = datenum(2020,12,24):datenum(2021,5,9);

% populate ila.instrument
ila.Read('prod');

% interpolate SST / SSS / LatLon
par_temp = ila.instrument.(cell2mat(ila.cfg.instruments2run)).prod.a;
par_temp.dt = datetime(par_temp.dt,'ConvertFrom','datenum');
% Convert to uE/cm^2/s for SeaBASS
par_temp.par = par_temp.par./10000;
par_temp.par_sd = par_temp.par_sd./10000;

% merge lat, lon, sst, sss
replace_consecutive_nan = 1*60; % 1h
par = merge_timeseries(par_temp, tsg, {'lat', 'lon', 'sst', 'sss'});
par = merge_timeseries(par, nmea, {'lat', 'lon'});
par = merge_timeseries(par, latlon, {'lat', 'lon'});
par = merge_timeseries(par, meteo_ftp, {'lat', 'lon'}, '', replace_consecutive_nan);

% Remove NaN
par(any(isnan([par.lat par.lon]),2), :) = [];
if any(isnan([par.lat par.lon]),2)
  par(any(isnan([par.lat par.lon]),2), :) = [];
  warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
    sum(any(isnan([suvf.lat suvf.lon]),2)), replace_consecutive_nan/60)
end

% add units and precision
par = renamevars(par, {'sst','sss','par_n'} ,{'t','s','bincount'});
par.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'uE/cm^2/s', 'uE/cm^2/s', 'none'};
par.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f', '%.2f'};

% sort by date
par = sortrows(par, 'dt');

ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod, ...
  'plots', [cruise '_PAR_timeseries']), 'jpg')
close figure 1

filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, ila.cfg.instruments2run{:}, ...
  datetime(min(par.dt), 'Format', 'yyyyMMdd'), datetime(max(par.dt), 'Format', 'yyyyMMdd'), ...
  datetime('today', 'Format', 'yyyyMMdd'));

% export product to SeaBASS format
ila.meta.documents = [cruise '_PAR_ProcessingReport_V2.pdf'];
ila.meta.calibration_files = 'PAR-50168_CalSheet.pdf';
exportSeaBASS(fullfile(ila.instrument.FLOW.path.prod, filename),...
    ila.meta,...
    par,...
    {'', '', ''});
sprintf('%s_InLine_%s_Product.sb saved', cruise, cell2mat(ila.cfg.instruments2run))

par = renamevars(par, {'t','s','bincount'}, {'sst','sss','par_n'});

% convert to uE/m^2/s for mat file
par.par = par.par.*10000;
par.par_sd = par.par_sd.*10000;
par.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'uE/m^2/s', 'uE/m^2/s', 'none'};

% save PAR prod
fprintf('Export to mat and csv... ');
filename = strrep(filename, '.sb', '');
save(fullfile(ila.instrument.FLOW.path.prod, filename), 'par');
writetable(par, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
fprintf('Done\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALFA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'ALFA011'}; % {'FLOW', 'TSG', 'BB3','PAR', 'WSCD1082P'}
ila.cfg.days2run = datenum(2020,12,24):datenum(2021,5,9);

% populate ila.instrument
ila.Read('prod');

% interpolate SST / SSS / LatLon
alfa_temp = ila.instrument.(cell2mat(ila.cfg.instruments2run)).prod.pd;
alfa_temp.dt = datetime(alfa_temp.dt,'ConvertFrom','datenum');
% merge lat, lon, sst, sss
replace_consecutive_nan = 72*60; % 72h
alfa = merge_timeseries(alfa_temp, tsg, {'lat', 'lon', 'sst', 'sss'});
alfa = merge_timeseries(alfa, nmea, {'lat', 'lon'});
alfa = merge_timeseries(alfa, latlon, {'lat', 'lon'});
alfa = merge_timeseries(alfa, meteo_ftp, {'lat', 'lon'}, '', replace_consecutive_nan);

% Remove NaN
alfa(any(isnan([alfa.lat alfa.lon]),2), :) = [];
if any(isnan([alfa.lat alfa.lon]),2)
  alfa(any(isnan([alfa.lat alfa.lon]),2), :) = [];
  warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
    sum(any(isnan([alfa.lat alfa.lon]),2)), replace_consecutive_nan/60)
end

% add units and precision
% alfa = renamevars(alfa, {'sst','sss','par_n'}, {'t','s','bincount'});
alfa.Properties.VariableUnits = repmat({''}, 1, size(alfa, 2));
alfa.Properties.VariableDescriptions = repmat({''}, 1, size(alfa, 2));

% sort by date
alfa = sortrows(alfa, 'dt');

ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod, ...
  'plots', [cruise '_ALFA_timeseries']), 'jpg')
close figure 96

filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, ila.cfg.instruments2run{:}, ...
  datetime(min(alfa.dt), 'Format', 'yyyyMMdd'), datetime(max(alfa.dt), 'Format', 'yyyyMMdd'), ...
  datetime('today', 'Format', 'yyyyMMdd'));

% % export product to SeaBASS format
% ila.meta.documents = [cruise '_' ila.cfg.instruments2run{:} '_ProcessingReport_V2.pdf'];
% ila.meta.calibration_files = [ila.cfg.instruments2run{:} '_CalSheet.pdf'];
% exportSeaBASS(fullfile(ila.instrument.FLOW.path.prod, filename),...
%     ila.meta,...
%     alfa,...
%     {'', '', ''});
% sprintf('%s_InLine_%s_Product.sb saved', cruise, cell2mat(ila.cfg.instruments2run))

% alfa = renamevars(alfa, {'t','s','bincount'}, {'sst','sss','par_n'});

% % convert to uE/m^2/s for mat file
% alfa.par = alfa.par.*10000;
% alfa.par_sd = alfa.par_sd.*10000;
% alfa.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'uE/m^2/s', 'uE/m^2/s', 'none'};

% save ALFA prod
fprintf('Export to mat and csv... ');
filename = strrep(filename, '.sb', '');
save(fullfile(ila.instrument.FLOW.path.prod, filename), 'alfa');
writetable(alfa, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
fprintf('Done\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BB3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'BB31502'}; % {'FLOW', 'TSG', 'BB31502','PAR', 'WSCD1082'}
ila.cfg.days2run = datenum(2020,12,24):datenum(2021,5,9);

% populate ila.instrument
ila.Read('prod');

% interpolate SST / SSS / LatLon
bb3_lambda = ila.instrument.(ila.cfg.instruments2run{:}).lambda;
bb3_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.p;
bb3_temp.dt = datetime(bb3_temp.dt,'ConvertFrom','datenum');
% merge lat, lon, sst, sss
replace_consecutive_nan = 72*60; % 72h
bb3_temp = merge_timeseries(bb3_temp, tsg, {'lat', 'lon', 'sst', 'sss'});
bb3_temp = merge_timeseries(bb3_temp, nmea, {'lat', 'lon'});
bb3_temp = merge_timeseries(bb3_temp, latlon, {'lat', 'lon'});
bb3_temp = merge_timeseries(bb3_temp, meteo_ftp, {'lat', 'lon'}, '', replace_consecutive_nan);

% Remove NaN
bb3_temp(any(isnan([bb3_temp.lat bb3_temp.lon]),2), :) = [];
if any(isnan([bb3_temp.lat bb3_temp.lon]),2)
  bb3_temp(any(isnan([bb3_temp.lat bb3_temp.lon]),2), :) = [];
  warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
    sum(any(isnan([bb3_temp.lat bb3_temp.lon]),2)), replace_consecutive_nan/60)
end

% add units and precision
bb3 = bb3_temp(:, ~contains(bb3_temp.Properties.VariableNames, {'gamma_bbp', 'poc', 'cphyto'}));
bb3.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', '1/m/sr', '1/m', '1/m/sr', 'none'};
bb3.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f'};
bb3 = renamevars(bb3, {'sst', 'sss', 'betap','betap_sd','betap_n'}, {'t','s','VSF_124ang','VSF_124ang_sd','bincount'});

% sort by date
bb3 = sortrows(bb3, 'dt'); 
bb3.bbp(bb3.bbp<0)=NaN;
bb3.VSF_124ang_sd(bb3.VSF_124ang<0)=NaN;
bb3.VSF_124ang(bb3.VSF_124ang<0)=NaN;
bb3(all(isnan(bb3.bbp),2),:)=[];

% remove dobious values
bb3(bb3.dt > datetime(2021,2,17,0,0,0) & bb3.dt < datetime(2021,2,23,18,30,0), :) = [];
bb3(bb3.dt > datetime(2021,3,23,21,25,0) & bb3.dt < datetime(2021,3,25,0,0,0), :) = [];
bb3(bb3.dt > datetime(2021,3,25,3,14,0) & bb3.dt < datetime(2021,3,31,0,0,0), :) = [];
% bb3(any(isnan([bb3.lat bb3.lon]),2), :) = [];

%%% BB 3D plots %%%
save_figures = true;
ila.SpectralQC('BB', {'prod'}, save_figures); % AC or BB
close all

% plot time series
ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_BB3_BBparticulate_timeseries']), 'fig')
close figure 86
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_BB3_POCparticulate_timeseries']), 'fig')
close figure 85
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_BB3_BBdissolved_timeseries']), 'fig')
close figure 87

SimpleMap(bb3.bbp(:,2), bb3(:,1:3), 'bbp (532 nm) [m^-^1]')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_bb3_bbp_map']), 'jpg')
close figure 1

filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, ila.cfg.instruments2run{:}, ...
  datetime(min(bb3.dt), 'Format', 'yyyyMMdd'), datetime(max(bb3.dt), 'Format', 'yyyyMMdd'), ...
  datetime('today', 'Format', 'yyyyMMdd'));

% export product to SeaBASS format
ila.meta.documents = sprintf('%s_BB3_ProcessingReport_v%s.pdf', cruise, datetime('today', 'Format', 'yyyyMMdd'));
ila.meta.calibration_files = 'BB3-1502_(470-532-650nm)_CalSheet_20220917.pdf';
exportSeaBASS(fullfile(ila.instrument.FLOW.path.prod, filename),...
    ila.meta,...
    bb3,...
    {string(bb3_lambda),string(bb3_lambda),string(bb3_lambda),''});
fprintf('%s_InLine_%s_Particulate_v%s.sb saved\n', cruise, ila.cfg.instruments2run{:}, ...
  datetime('today', 'Format', 'yyyyMMdd'))

bb3_temp.Properties.VariableNames = {'dt','lat','lon','sst','sss','VSF124','bbp','VSF124_sd','bincount','gamma_bbp','poc','cphyto_bbp'};
bb3_temp.Properties.VariableUnits = {'','degrees','degrees','degreesC','PSU','m^-1.sr^-1','m^-1','m^-1.sr^-1','none','unitless','ug/L','ug/L'};
bb3_temp.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.2f','%.4f','%.2f','%.2f'};

bb3 = bb3_temp;

% save BB3 prod
fprintf('Export to mat and csv... ');
filename = strrep(filename, '.sb', '');
save(fullfile(ila.instrument.FLOW.path.prod, filename), 'bb3', 'bb3_lambda');
writetable(bb3, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
fprintf('Done\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HBB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'HyperBB8005'}; % {'FLOW', 'TSG', 'BB31502','PAR', 'WSCD1082', 'HBB'}
ila.cfg.days2run = datenum(2021,8,16):datenum(2022,9,17);

% populate ila.instrument
ila.Read('prod');

% interpolate SST / SSS / LatLon
hbb_lambda = ila.instrument.(ila.cfg.instruments2run{:}).lambda;
hbb_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.p;
hbb_temp.dt = datetime(hbb_temp.dt,'ConvertFrom','datenum');
% merge lat, lon, sst, sss
replace_consecutive_nan = 72*60; % 72h
hbb_temp = merge_timeseries(hbb_temp, tsg, {'lat', 'lon', 'sst', 'sss'});
hbb_temp = merge_timeseries(hbb_temp, nmea, {'lat', 'lon'});
hbb_temp = merge_timeseries(hbb_temp, latlon, {'lat', 'lon'});
hbb_temp = merge_timeseries(hbb_temp, meteo_ftp, {'lat', 'lon'}, '', replace_consecutive_nan);

% Remove NaN
if any(any(isnan([hbb_temp.lat hbb_temp.lon]),2))
  warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
    sum(any(isnan([hbb_temp.lat hbb_temp.lon]),2)), replace_consecutive_nan/60)
  hbb_temp(any(isnan([hbb_temp.lat hbb_temp.lon]),2), :) = [];
end

% add units and precision
hbb = hbb_temp(:, ~contains(hbb_temp.Properties.VariableNames, {'gamma_bbp', 'poc', 'cphyto'}));
hbb.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', '1/m/sr', '1/m', '1/m/sr', 'none'};
hbb.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f'};
hbb = renamevars(hbb, {'sst', 'sss', 'betap','betap_sd','betap_n'}, {'t','s','VSF_124ang','VSF_124ang_sd','bincount'});

% sort by date
hbb = sortrows(hbb, 'dt'); 
hbb.bbp(hbb.bbp<0)=NaN;
hbb.VSF_124ang_sd(hbb.VSF_124ang<0)=NaN;
hbb.VSF_124ang(hbb.VSF_124ang<0)=NaN;
hbb(all(isnan(hbb.bbp),2),:)=[];

%%% HBB 3D plots %%%
save_figures = true;
ila.SpectralQC('BB', {'prod'}, save_figures); % AC or BB
close all

% plot time series
ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_HBB_gammabbp_bbp550']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_HBB_gammabbp_bbp550']), 'fig')
close figure 11
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_HBB_BBparticulate_timeseries']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_HBB_BBparticulate_timeseries']), 'fig')
close figure 21
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_HBB_POCparticulate_timeseries']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_HBB_POCparticulate_timeseries']), 'fig')
close figure 20
% saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_HBB_BBdissolved_timeseries']), 'jpg')
% saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_HBB_BBdissolved_timeseries']), 'fig')
% close figure 25

SimpleMap(hbb.bbp(:,hbb_lambda == 530), hbb(:,1:3), 'bbp (530 nm) [m^-^1]')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_hbb_bbp_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_hbb_bbp_map']), 'fig')
close figure 1

filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, ila.cfg.instruments2run{:}, ...
  datetime(min(hbb.dt), 'Format', 'yyyyMMdd'), datetime(max(hbb.dt), 'Format', 'yyyyMMdd'), ...
  datetime('today', 'Format', 'yyyyMMdd'));

% export product to SeaBASS format
ila.meta.documents = sprintf('%s_HBB_ProcessingReport_v%s.pdf', cruise, datetime('today', 'Format', 'yyyyMMdd'));
ila.meta.calibration_files = 'HBB8005_CharSheet.pdf';
exportSeaBASS(fullfile(ila.instrument.FLOW.path.prod, filename),...
    ila.meta,...
    hbb,...
    {string(hbb_lambda)',string(hbb_lambda)',string(hbb_lambda)',''});
fprintf('%s_InLine_%s_Particulate_v%s.sb saved\n', cruise, ila.cfg.instruments2run{:}, ...
  datetime('today', 'Format', 'yyyyMMdd'))

hbb_temp.Properties.VariableNames = {'dt','lat','lon','sst','sss','VSF124','bbp','VSF124_sd','bincount','gamma_bbp','poc','cphyto_bbp'};
hbb_temp.Properties.VariableUnits = {'','degrees','degrees','degreesC','PSU','m^-1.sr^-1','m^-1','m^-1.sr^-1','none','unitless','ug/L','ug/L'};
hbb_temp.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.2f','%.4f','%.2f','%.2f'};

hbb = hbb_temp;

% save HBB prod
fprintf('Export to mat and csv... ');
filename = strrep(filename, '.sb', '');
save(fullfile(ila.instrument.FLOW.path.prod, filename), 'hbb', 'hbb_lambda');
writetable(hbb, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
fprintf('Done\n');


%% diel cycle gamma_bbp
foo_hbb = table();
foo_dt = [min(hbb.dt) max(hbb.dt)];
foo_hbb.dt = dateshift((foo_dt(1):minutes(1):foo_dt(2))','start','minute');
id_nan = isnan(hbb.gamma_bbp);
% create new fdom table interpolating for gaps < 60 min
foo_hbb.gamma_bbp = interp1(hbb.dt, hbb.gamma_bbp, foo_hbb.dt, 'linear');
foo_hbb.gamma_bbp = fillmissing(foo_hbb.gamma_bbp, 'nearest');
% smooth WSCD cdom data removing frequencies lower than 1800min == 30h
d = designfilt('highpassiir', 'FilterOrder', 1, 'PassbandFrequency', 1/(60*1800), 'SampleRate', 1/60);
foo_hbb.gamma_bbp = filtfilt(d, foo_hbb.gamma_bbp);
% smooth WSCD cdom data removing frequencies higher than 1080min == 18h
d = designfilt('lowpassiir', 'FilterOrder', 1, 'PassbandFrequency', 1/(60*1080), 'SampleRate', 1/60);
foo_hbb.gamma_bbp = filtfilt(d, foo_hbb.gamma_bbp);
hbb.gamma_bbp_18_30h_cycle = interp1(foo_hbb.dt, foo_hbb.gamma_bbp, hbb.dt, 'linear');
hbb.gamma_bbp_18_30h_cycle(id_nan) = NaN;

figure; hold on;
scatter(hbb.dt, hbb.gamma_bbp, 20, 'filled')
scatter(hbb.dt, hbb.gamma_bbp_18_30h_cycle, 20, 'filled')
close figure 1

% save HBB prod
fprintf('Export to mat and csv... ');
filename = strrep(filename, 'HyperBB8005', 'HyperBB8005_with_smoothed_gamma_bbp');
save(fullfile(ila.instrument.FLOW.path.prod, filename), 'hbb', 'hbb_lambda');
writetable(hbb, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
fprintf('Done\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_leg = {...
    datenum(2020,12,24):datenum(2021,02,05);...
    datenum(2021,02,17):datenum(2021,3,17);...
    datenum(2021,03,20):datenum(2021,4,7);...
    datenum(2021,04,11):datenum(2021,4,20);...
    datenum(2021,04,26):datenum(2021,5,9);...
    datenum(2021,04,26):datenum(2021,5,9);...
    datenum(2021,8,18):datenum(2021,9,18);...
    datenum(2021,9,24):datenum(2021,10,9);...
    datenum(2021,10,17):datenum(2021,11,3);...
    datenum(2021,11,11):datenum(2021,11,25);...
    datenum(2021,12,5):datenum(2021,12,27);...
    datenum(2022,1,1):datenum(2022,2,27);...
    datenum(2022,3,5):datenum(2022,4,23);...
    datenum(2022,5,1):datenum(2022,6,1);...
    datenum(2022,6,6):datenum(2022,7,9);...
    datenum(2022,7,16):datenum(2022,8,5);...
    datenum(2022,8,6):datenum(2022,9,17);
    };

list_dev = [repmat({fullfile(path_dev, 'acs057_20200129.dev')}, size(list_leg, 1)-1, 1);...
  {fullfile(path_dev, 'acs348_20220606.dev')}];

list_instru = [repmat({'ACS57'}, size(list_leg, 1)-1, 1); {'ACS348'}];

data_AC = struct('particulate', [], 'product', []);
acs = [];

for i=1:size(list_instru,1)
  cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master/')
  ila = InLineAnalysis(['cfg' filesep cruise '_cfg.m']);
  ila.cfg.instruments2run = list_instru(i);
  ila.cfg.days2run = list_leg{i};

  % populate ila.instrument
  ila.Read('prod');

  % get wavelength
  [ila.instrument.(list_instru{i}).lambda_c, ...
    ila.instrument.(list_instru{i}).lambda_a] = importACSDeviceFile(list_dev{i});
  
  % set section name
  ref = sprintf('%s_%s_%s', list_instru{i}, datetime(ila.cfg.days2run(1), 'ConvertFrom','datenum', 'Format', 'yyyyMMdd'), ...
    datetime(ila.cfg.days2run(end), 'ConvertFrom', 'datenum', 'Format', 'yyyyMMdd'));
  
  % remove flagged data
  flag = read_flagbit(ila.instrument.(list_instru{i}).prod.p.flag_bit, 'AC');
  % flag_info = set_flagbit('ACS');

  % flag.HH_G50_flag(flag.HH_G50_flag & ila.instrument.(list_instru{i}).prod.p.HH_G50 > 0 & ila.instrument.(list_instru{i}).prod.p.HH_G50 <= 500) = false;
  % ila.instrument.(list_instru{i}).prod.p.flag_bit = set_flagbit(flag);
  
  % remove flagged products
  ila.instrument.(list_instru{i}).prod.p.poc(flag.poc_flag) = NaN;
  ila.instrument.(list_instru{i}).prod.p.chl_ap676lh(flag.chl_ap676lh_flag) = NaN;
  ila.instrument.(list_instru{i}).prod.p.gamma(flag.cp_bubbles) = NaN;
  ila.instrument.(list_instru{i}).prod.p.poc(flag.cp_bubbles) = NaN;
  ila.instrument.(list_instru{i}).prod.p.chl_Halh(flag.chl_Halh_flag) = NaN;
  ila.instrument.(list_instru{i}).prod.p.chl_Halh(flag.cp_bubbles) = NaN;
  ila.instrument.(list_instru{i}).prod.p.HH_mphi(flag.cp_bubbles) = NaN;
  ila.instrument.(list_instru{i}).prod.p.HH_G50(flag.HH_G50_flag) = NaN;
  ila.instrument.(list_instru{i}).prod.p.HH_G50(flag.cp_bubbles) = NaN;
  % ila.instrument.(list_instru{i}).prod.p.chl_Halh(flag.chlratio_flag) = NaN;
  
  % % remove suspicious products
  % ila.instrument.(list_instru{i}).prod.p.gamma(flag.gamma_suspicious) = NaN;
  % % ila.instrument.(list_instru{i}).prod.p.gamma(acs_prod.gamma > 2) = NaN;
  % ila.instrument.(list_instru{i}).prod.p.poc(flag.ap_bubbles) = NaN;
  % ila.instrument.(list_instru{i}).prod.p.poc(flag.poc_suspicious) = NaN;
  % ila.instrument.(list_instru{i}).prod.p.chl_ap676lh(flag.chl_ap676lh_suspicious) = NaN;
  % ila.instrument.(list_instru{i}).prod.p.chl_Halh(flag.chl_Halh_suspicious) = NaN;
  % ila.instrument.(list_instru{i}).prod.p.HH_mphi(flag.HH_G50_mphi_suspicious) = NaN;
  % ila.instrument.(list_instru{i}).prod.p.HH_G50(flag.HH_G50_mphi_suspicious) = NaN;

  ila.visProd_timeseries()
  saveGraph(fullfile(ila.instrument.(list_instru{i}).path.prod, 'plots', [cruise '_' ref '_ACS_prod_regressions']), 'jpg')
  saveGraph(fullfile(ila.instrument.(list_instru{i}).path.prod, 'plots', [cruise '_' ref '_ACS_prod_regressions']), 'fig')
  close figure 11
  saveGraph(fullfile(ila.instrument.(list_instru{i}).path.prod, 'plots', [cruise '_' ref '_ACS_prod_timeseries']), 'jpg')
  saveGraph(fullfile(ila.instrument.(list_instru{i}).path.prod, 'plots', [cruise '_' ref '_ACS_prod_timeseries']), 'fig')
  close figure 10

  % extract AC data from obj
  AC = ila.instrument.(list_instru{i}).prod.p;
  AC.dt = datetime(AC.dt,'ConvertFrom','datenum');

  % build AC table: merge lat, lon, sst, sss, and fdom
  replace_consecutive_nan = 1*60; % 1h
  AC = merge_timeseries(AC, fdom, {'fdom'}, '', replace_consecutive_nan);
  AC = merge_timeseries(AC, tsg, {'lat', 'lon', 'sst', 'sss'});
  AC = merge_timeseries(AC, nmea, {'lat', 'lon'});
  AC = merge_timeseries(AC, latlon, {'lat', 'lon'});
  replace_consecutive_nan = 72*60; % 72h
  AC = merge_timeseries(AC, meteo_ftp, {'lat', 'lon'}, '', replace_consecutive_nan);
  % rename variables for SeaBASS
  AC = renamevars(AC, {'sst','sss'}, {'t','s'});
  
  % Remove NaN
  if any(any(isnan([AC.lat AC.lon]),2))
    warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
      sum(any(isnan([AC.lat AC.lon]),2)), replace_consecutive_nan/60)
    AC(any(isnan([AC.lat AC.lon]),2), :) = [];
  end

  % split into particulate table
  id_particulate = logical(sum(categorical(AC.Properties.VariableNames) == {'dt','lat','lon','t','s','fdom','ap','ap_sd','cp','cp_sd','ap_n'}'));
  data_AC.particulate.(ref) = AC(:, id_particulate);
  % add variable units and description
  data_AC.particulate.(ref).Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU','ppb', '1/m', '1/m', '1/m', '1/m', 'none'};
  data_AC.particulate.(ref).Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.0f'};
  % rename variables for SeaBASS
  data_AC.particulate.(ref) = renamevars(data_AC.particulate.(ref), {'fdom','ap_n'}, {'cdomf','bincount'});

  % split into product table
  id_product = logical(sum(categorical(AC.Properties.VariableNames) == {'ap', 'ap_sd', 'cp', 'cp_sd'}'));
  data_AC.product.(ref) = AC(:, ~id_product);
  % add variable units and description
  data_AC.product.(ref).Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU','ppb', ...
              '1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m', ...
              'ug/L','unitless','ug/L','unitless','unitless','ug/L','microns','unitless','unitless','unitless','unitless'};
  data_AC.product.(ref).Properties.VariableDescriptions = [{''} repmat({'%.4f'}, 1, size(data_AC.product.(ref),2) - 4) ...
    repmat({'%i'}, 1, 3)];
  % rename variables for SeaBASS
  data_AC.product.(ref) = renamevars(data_AC.product.(ref), {'poc','gamma','chl_ap676lh', 'fdom'}, ...
    {'POC_cp','cp_gamma','Chl_lineheight', 'cdomf'});

  %%% AC 3D plots %%%
  save_figures = true;
  ila.SpectralQC('AC', {'prod'}, save_figures); % AC or BB
  close all
  
  filename = sprintf('%s_InLine_%s_Particulate_v%s.sb', cruise, ref, datetime('today', 'Format', 'yyyyMMdd'));
  % export particulate to SeaBASS format
  ila.meta.documents = sprintf('%s_ACS_ProcessingReport_v%s.pdf', cruise, datetime('today', 'Format', 'yyyyMMdd'));
  [~, calfile] = fileparts(list_dev{i});
  ila.meta.calibration_files = [calfile '.dev'];
  exportSeaBASS(fullfile(ila.instrument.(list_instru{i}).path.prod, filename),...
      ila.meta,...
      data_AC.particulate.(ref),...
      {'', string(ila.instrument.(list_instru{i}).lambda_a),...
      string(ila.instrument.(list_instru{i}).lambda_c),...
      string(ila.instrument.(list_instru{i}).lambda_a),...
      string(ila.instrument.(list_instru{i}).lambda_c),''});
  fprintf('%s saved\n', filename)
  
  lambda.a = ila.instrument.(list_instru{i}).lambda_a;
  lambda.c = ila.instrument.(list_instru{i}).lambda_c;
  acs_part = data_AC.particulate.(ref);
  
  filename = strrep(filename, '.sb', '');
  % rename variables for .mat
  acs_part = renamevars(acs_part, {'t','s','cdomf'}, {'sst','sss','fdom'});
  save(fullfile(ila.instrument.(list_instru{i}).path.prod, filename), 'acs_part', 'lambda');

  % ACS merged prod
  acs = [acs; data_AC.product.(ref)];
end

% sort by date
acs_prod = sortrows(acs, 'dt');

if any(strcmp(acs_prod.Properties.VariableNames, 'chl_ap676lh'))
  acs_prod = renamevars(acs_prod, 'chl_ap676lh', 'Chl_lineheight');
end
if any(strcmp(acs.Properties.VariableNames, 'poc'))
  acs_prod = renamevars(acs_prod, 'poc', 'POC_cp');
end
if any(strcmp(acs.Properties.VariableNames, 'gamma'))
  acs_prod = renamevars(acs_prod, 'gamma', 'cp_gamma');
end

ila.instrument.(list_instru{1}).prod.p = acs_prod;

% figure()
% subplot(2,3,1); histogram(acs_prod.POC_cp); xlabel('[POC] cp (mg.m^{-3})')
% subplot(2,3,2); histogram(acs_prod.Chl_lineheight); xlabel('a_{p676} line height [chl a] (mg.m^{-3})')
% subplot(2,3,3); histogram(acs_prod.cp_gamma); xlabel('gamma cp (mg.m^{-3})')
% subplot(2,3,4); histogram(acs_prod.chl_Halh); xlabel('Houskeeper [chl] (mg.m^{-3})')
% subplot(2,3,5); histogram(acs_prod.HH_mphi); xlabel('H&H phytoplankton slope size distribution')
% subplot(2,3,6); histogram(acs_prod.HH_G50); xlabel('H&H phytoplankton G50: cross-sectional area (\mum)')

ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.(list_instru{i}).path.prod, 'plots', [cruise '_ACS_prod_regressions']), 'jpg')
saveGraph(fullfile(ila.instrument.(list_instru{i}).path.prod, 'plots', [cruise '_ACS_prod_regressions']), 'fig')
close figure 11
saveGraph(fullfile(ila.instrument.(list_instru{i}).path.prod, 'plots', [cruise '_ACS_prod_timeseries']), 'jpg')
saveGraph(fullfile(ila.instrument.(list_instru{i}).path.prod, 'plots', [cruise '_ACS_prod_timeseries']), 'fig')
close figure 10

SimpleMap(acs_prod.chl_Halh, acs_prod(:,1:3), 'Houskeeper [chl] (mg.m^{-3})')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_ACS_chl_Houskeeper_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_ACS_chl_Houskeeper_map']), 'fig')
close figure 1

SimpleMap(acs_prod.HH_G50, acs_prod(:,1:3), 'H&H phytoplankton G50: cross-sectional area (\mum)')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_ACS_H&H_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_ACS_H&H_map']), 'fig')
close figure 1

SimpleMap(acs_prod.POC_cp, acs_prod(:,1:3), '[POC] cp (mg.m^{-3})')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_ACS_POC_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_ACS_POC_map']), 'fig')
close figure 1

SimpleMap(acs_prod.Chl_lineheight, acs_prod(:,1:3), 'a_{p676} line height [chl a] (mg.m^{-3})')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_ACS_chl_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_ACS_chl_map']), 'fig')
close figure 1

SimpleMap(acs_prod.cp_gamma, acs_prod(:,1:3), 'gamma cp (unitless)')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_ACS_gamma_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_ACS_gamma_map']), 'fig')
close figure 1

% acs_prod.chl_Halh = [];
% acs_prod.HH_G50 = [];

filename = sprintf('%s_InLine_ACS_%s_%s_Products_full_v%s', cruise, ...
  datetime(min(acs_prod.dt), 'Format', 'yyyyMMdd'),...
  datetime(max(acs_prod.dt), 'Format', 'yyyyMMdd'),...
  datetime('today', 'Format', 'yyyyMMdd'));

% rename variables for .mat
acs_prod = renamevars(acs_prod, {'t','s','cdomf'}, {'sst','sss','fdom'});

% save AC prod
fprintf('Export to mat and csv... ');
save(fullfile(ila.instrument.(list_instru{i}).path.prod, filename), 'acs_prod');
writetable(acs_prod, fullfile(ila.instrument.(list_instru{i}).path.prod, [filename '.csv']));
fprintf('Done\n');

% rename variables for SeaBASS
acs_prod = renamevars(acs_prod, {'sst','sss','fdom'}, {'t','s','cdomf'});

% keep only old variables
acs_prod = acs_prod(:, contains(acs_prod.Properties.VariableNames, {'dt', 'lat', 'lon', ...
  't','s','cdomf','Chl_lineheight','POC_cp','cp_gamma','ap_n','cp_n','flag_bit'}) & ...
  ~contains(acs_prod.Properties.VariableNames, {'agaus'}));

% export product to SeaBASS format
filename = [filename '.sb'];
filename = strrep(filename, 'Products_full', 'Products');
exportSeaBASS(fullfile(ila.instrument.(list_instru{i}).path.prod, filename), ila.meta, acs_prod);
fprintf('%s saved\n', filename)

% rename variables for .mat
acs_prod = renamevars(acs_prod, {'t','s','cdomf'}, {'sst','sss','fdom'});

% save AC prod
fprintf('Export to mat and csv... ');
filename = strrep(filename, '.sb', '');
save(fullfile(ila.instrument.(list_instru{i}).path.prod, filename), 'acs_prod');
writetable(acs_prod, fullfile(ila.instrument.(list_instru{i}).path.prod, [filename '.csv']));
fprintf('Done\n');


%%
acs_part(~ismember(acs_part.dt, acs.dt), :) = [];
bp_div_chl = (acs_part.cp(:,40) - acs_part.ap(:,40)) ./ acs.Chl_lineheight;

SimpleMap(bp_div_chl, acs(:,1:3), 'b_p / [chl a]')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_ACS_bp_div_chl_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_ACS_bp_div_chl_map']), 'fig')
close figure 1

%%
acs_part(~ismember(acs_part.dt, acs.dt), :) = [];

SimpleMap(acs_part.cp(:,40) - acs_part.ap(:,40), acs(:,1:3), 'b_p [m^{-1}')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_ACS_bp_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_ACS_bp_map']), 'fig')
close figure 1

%% bbp/bp graph to validate BB3 
bp_interp = interp1(bp_bb3(:,1), bp_bb3(:,2:end), ...
  datenum(bb3.dt), 'linear', 'extrap'); % extrap needed for first minute of data

% get spectrumRGB
C = reshape(spectrumRGB(bb3_lambda), max(size(bb3_lambda)),  3);

figure(); hold on
scatter(bb3.dt, bb3.bbp./bp_interp, 15, C, 'filled')
xlim([min(bb3.dt) max(bb3.dt)])
ylabel('{\itbb_p}/{\itb_p}')
legend(cellfun(@(x) [x 'nm'], cellstr(num2str(bb3_lambda')), 'un', 0), 'FontSize', 14)
set(gca, 'FontSize', 13)
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_BB3_bbp_over_bp']), 'jpg')
close figure 1

%% bbp/bp graph to validate HBB data
bp_interp = interp1(bp_hbb(:,1), bp_hbb(:,2:end), ...
  datenum(hbb.dt), 'linear', 'extrap'); % extrap needed for first minute of data

% get spectrumRGB
C = reshape(spectrumRGB(hbb_lambda), max(size(hbb_lambda)),  3);

figure(); hold on
scatter(hbb.dt, hbb.bbp./bp_interp, 15, C, 'filled')
xlim([min(hbb.dt) max(hbb.dt)])
ylabel('{\itbb_p}/{\itb_p}')
legend(cellfun(@(x) [x 'nm'], cellstr(num2str(hbb_lambda')), 'un', 0), 'FontSize', 14)
set(gca, 'FontSize', 13)
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_HBB_bbp_over_bp']), 'jpg')
close figure 1

%% GAYOSO map centered on stations
stations = readtable('/Users/gui/Documents/Maine/mobilisation/Microbiome/Gayoso/Rosette/Station_UTC.txt');
stations.Properties.VariableNames = {'st', 'dt_start', 'dt_end', 'lat', 'lon'};

% lagrangian #1 A061-A063
% lagrangian #2 A064-A067
% create table
sz_tbl = size(stations,1);
stations.st{sz_tbl+1} = 'Lagrangian #1';
stations.dt_start(sz_tbl+1) = stations.dt_start(3);
stations.dt_end(sz_tbl+1) = stations.dt_end(5);
stations.lat(sz_tbl+1) = stations.lat(4);
stations.lon(sz_tbl+1) = stations.lon(4);

stations.st{sz_tbl+2}  = 'Lagrangian #2';
stations.dt_start(sz_tbl+2) = stations.dt_start(6);
stations.dt_end(sz_tbl+2) = stations.dt_end(9);
stations.lat(sz_tbl+2) = stations.lat(7);
stations.lon(sz_tbl+2) = stations.lon(7);

% create table
stations.id = cell(size(stations,1), 1);
stations.toplot = cell(size(stations,1), 1);
for i = 1:size(stations, 1)
  stations.id{i} = tsg.dt >= stations.dt_start(i) & tsg.dt <= stations.dt_end(i);
  if strcmp(stations.st{i}, 'Lagrangian #1')
    stations.toplot{i} = {'A061','A062','A063'};
  elseif strcmp(stations.st{i}, 'Lagrangian #2')
    stations.toplot{i} = {'A064','A065','A066','A067'};
  else
    stations.toplot{i} = stations.st{i};
  end
end

stations.bbp_text_color = {'white','black','white','black','white','white','black',...
  'black','black','black','black','white','white','black'}';
stations.chl_text_color = {'white','white','white','white','white','white','white',...
  'white','white','black','black','black','white','white'}';
stations.chlbbp_text_color = {'black','black','black','black','black','black','black',...
  'black','black','black','black','black','black','black'}';

path_sat = '/Users/gui/Documents/Maine/Data/TaraPacific/remote_sensing/time_series/level3_timeseries/gayoso';
list_sat = {dir(fullfile(path_sat, 'gayoso*mat')).name}';
sat_start = datetime(cellfun(@(s) s{2}, cellfun(@(c) strsplit(c, '_'), list_sat, 'un', 0), 'un', 0), ...
  'InputFormat', 'yyyyMMdd-HHmm');
sat_end = datetime(cellfun(@(s) s{3}, cellfun(@(c) strsplit(c, '_'), list_sat, 'un', 0), 'un', 0), ...
  'InputFormat', 'yyyyMMdd-HHmm');

% merge chl into HBB table
hbb.chl = interp1(acs.dt, acs.Chl_lineheight, hbb.dt,'linear');
hbb.gamma = interp1(acs.dt, acs.cp_gamma, hbb.dt,'linear');
hbb.POC = interp1(acs.dt, acs.POC_cp, hbb.dt,'linear');
hbb.HH_G50 = interp1(acs.dt, acs.HH_G50, hbb.dt,'linear');
hbb(hbb.dt > datetime(2022,1,1),:) = [];
load(fullfile(path_sat, list_sat{1}))
hbb = hbb(hbb.lat >= min(min(lat))+0.25 & ...
  hbb.lat <= max(max(lat))-0.25 & ...
  hbb.lon >= min(min(lon))+0.25 & ...
  hbb.lon <= max(max(lon))-0.25, :);

for i = 1:size(stations,1)
  dif_sat = abs(stations.dt_start(i) - sat_start) + abs(stations.dt_end(i) - sat_end);
  sel_sat = list_sat{find(dif_sat == min(dif_sat), 1, 'first')};
  load(fullfile(path_sat, sel_sat))
  
  id_subset = lat < stations.lat(i) + 1 & lat > stations.lat(i) - 1 & ...
    lon < stations.lon(i) + 1 & lon > stations.lon(i) - 1;
  
  bbp550(~id_subset) = NaN;
  chl_poly(~id_subset) = NaN;
  lat(~id_subset) = NaN;
  lon(~id_subset) = NaN;
  
  sel_hbb = hbb(hbb.lat >= min(min(lat(id_subset))) & ...
    hbb.lat <= max(max(lat(id_subset))) & ...
    hbb.lon >= min(min(lon(id_subset))) & ...
    hbb.lon <= max(max(lon(id_subset))), :);
  
  % map bbp + hbb
  [~, cb] = QuickMap(bbp550, lat, lon, 'morgenstemning', [0 0.012], 'log'); % morgenstemning 0.00005 max(bbp550(:))
  cb.Label.String = 'bbp550 [m^{-1}]';
  hold on
  m_scatter(sel_hbb.lon, sel_hbb.lat, 100, sel_hbb.bbp(:, hbb_lambda == 530), ...
    'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
  m_scatter(sel_hbb.lon, sel_hbb.lat, 100, sel_hbb.bbp(:, hbb_lambda == 530), 'filled');
  if contains(stations.st{i}, 'Lagrangian')
    for j = 1:size(stations.toplot{i},2)
      topl = strcmp(stations.toplot{i}{j}, stations.st);
      m_scatter(stations.lon(topl), stations.lat(topl), 200, 'x', ...
        stations.bbp_text_color{i},'LineWidth', 3)
      m_text(stations.lon(topl)+0.2, stations.lat(topl), stations.st{topl}, 'Color',...
        stations.bbp_text_color{i}, 'FontSize', 15)
    end
  else
    m_scatter(stations.lon(i), stations.lat(i), 200, 'x', ...
      stations.bbp_text_color{i},'LineWidth', 3)
    m_text(stations.lon(i)+0.2, stations.lat(i), stations.st{i}, 'Color',...
      stations.bbp_text_color{i}, 'FontSize', 15)
  end
  title(sprintf('Station %s %s', stations.st{i}, datetime(stations.dt_start(i),'yyyy-mm-dd')), 'FontSize', 13)
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', ...
    sprintf('%s_Gayoso_bbp_map_station_%s_%s', cruise, stations.st{i}, ...
    datetime(stations.dt_start(i),'yyyymmdd'))), 'jpg')
  pause(5)
  close figure 1
  
  % map chl + acs chl
  [~, cb] = QuickMap(chl_poly, lat, lon, 'morgenstemning', [0 7], 'log'); % morgenstemning 0.00005 max(bbp550(:))
  cb.Label.String = '[Chl] [mg.m^{-3}]';
  hold on
  m_scatter(sel_hbb.lon, sel_hbb.lat, 100, sel_hbb.chl, ...
    'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
  m_scatter(sel_hbb.lon, sel_hbb.lat, 100, sel_hbb.chl, 'filled');
  if contains(stations.st{i}, 'Lagrangian')
    for j = 1:size(stations.toplot{i},2)
      topl = strcmp(stations.toplot{i}{j}, stations.st);
      m_scatter(stations.lon(topl), stations.lat(topl), 200, 'x', ...
        stations.chl_text_color{i},'LineWidth', 3)
      m_text(stations.lon(topl)+0.2, stations.lat(topl), stations.st{topl}, 'Color', ...
        stations.chl_text_color{i}, 'FontSize', 15)
    end
  else
    m_scatter(stations.lon(i), stations.lat(i), 200, 'x', ...
      stations.chl_text_color{i},'LineWidth', 3)
    m_text(stations.lon(i)+0.2, stations.lat(i), stations.st{i}, 'Color', ...
      stations.chl_text_color{i}, 'FontSize', 15)
  end
  title(sprintf('Station %s %s', stations.st{i}, datetime(stations.dt_start(i),'yyyy-mm-dd')), 'FontSize', 13)
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', ...
    sprintf('%s_Gayoso_chl_map_station_%s_%s', cruise, stations.st{i}, ...
    datetime(stations.dt_start(i),'yyyymmdd'))), 'jpg')
  pause(5)
  close figure 1
  
  % map bbp/chl + inline bbp/chl
  [~, cb] = QuickMap(bbp550./chl_poly, lat, lon, 'morgenstemning', [0.001 0.015], 'log'); % morgenstemning 0.00005 max(bbp550(:))
  cb.Label.String = 'bbp550/[Chl]';
  hold on
  m_scatter(sel_hbb.lon, sel_hbb.lat, 100, sel_hbb.bbp(:, hbb_lambda == 530)./sel_hbb.chl, ...
    'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
  m_scatter(sel_hbb.lon, sel_hbb.lat, 100, sel_hbb.bbp(:, hbb_lambda == 530)./sel_hbb.chl, 'filled');
  if contains(stations.st{i}, 'Lagrangian')
    for j = 1:size(stations.toplot{i},2)
      topl = strcmp(stations.toplot{i}{j}, stations.st);
      m_scatter(stations.lon(topl), stations.lat(topl), 200, 'x', ...
        stations.chlbbp_text_color{i},'LineWidth', 3)
      m_text(stations.lon(topl)+0.2, stations.lat(topl), stations.st{topl}, 'Color', ...
        stations.chlbbp_text_color{i}, 'FontSize', 15)
    end
  else
    m_scatter(stations.lon(i), stations.lat(i), 200, 'x', ...
      stations.chlbbp_text_color{i},'LineWidth', 3)
    m_text(stations.lon(i)+0.2, stations.lat(i), stations.st{i}, 'Color',...
      stations.chlbbp_text_color{i}, 'FontSize', 15)
  end
  title(sprintf('Station %s %s', stations.st{i}, datetime(stations.dt_start(i),'yyyy-mm-dd')), 'FontSize', 13)
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', ...
    sprintf('%s_Gayoso_bbpchl_map_station_%s_%s', cruise, stations.st{i}, ...
    datetime(stations.dt_start(i),'yyyymmdd'))), 'jpg')
  pause(5)
  close figure 1
end

%% Whole month Gayoso
satmonth = list_sat(sat_start > min(stations.dt_start) & sat_end < max(stations.dt_end));
satmonth = satmonth([3 7 11 15 21 24]);
load(fullfile(path_sat, satmonth{1}))
merged_chl = NaN(size(chl_poly,1), size(chl_poly,2), size(satmonth,1));
merged_bbp = NaN(size(chl_poly,1), size(chl_poly,2), size(satmonth,1));
merged_chl(:,:,1) = chl_poly;
merged_bbp(:,:,1) = bbp550;
for i = 2:size(satmonth,1)
  load(fullfile(path_sat, satmonth{i}))
  merged_chl(:,:,i) = chl_poly;
  merged_bbp(:,:,i) = bbp550;
end
month_chl = mean(merged_chl, 3, 'omitnan');
month_bbp = mean(merged_bbp, 3, 'omitnan');

% delete bad data point inline
hbb(hbb.bbp(:, hbb_lambda == 530)./hbb.chl > 0.02, :) = [];

% change text color
stations.bbp_text_color = {'white','white','white','black','black','white','white',...
  'white','white','black','black','white','white','black'}';
stations.chl_text_color = {'white','white','white','black','black','white','white',...
  'white','white','black','white','white','white','white'}';
stations.chlbbp_text_color = {'white','white','black','white','white','white','white',...
  'white','black','white','white','white','black','black'}';

stations.bbp_mark_color = {'white','white','white','white','white','white','white',...
  'white','white','black','black','white','white','black'}';
stations.chl_mark_color = {'white','white','white','white','white','white','white',...
  'white','white','black','white','white','white','white'}';
stations.chlbbp_mark_color = {'white','white','white','white','white','white','white',...
  'white','white','white','white','white','black','black'}';
stations.text_xpos = [1   1  -0.4 -2.5 -2.5  1   1  1  0   1 1 1 1 1]';
stations.text_ypos = [0 -0.5 -0.5   0   0.3 -0.4 0 0.3 0.5 0 0 0 0 0]';

% map bbp + hbb
[~, cb] = QuickMap(month_bbp, lat, lon, 'morgenstemning', [0 0.015], 'log'); % morgenstemning 0.00005 max(bbp550(:))
cb.Label.String = 'bbp550 [m^{-1}]';
hold on
% plot Inline
m_scatter(hbb.lon, hbb.lat, 100, hbb.bbp(:, hbb_lambda == 530), ...
  'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
m_scatter(hbb.lon, hbb.lat, 100, hbb.bbp(:, hbb_lambda == 530), ...
  'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
m_scatter(hbb.lon, hbb.lat, 100, hbb.bbp(:, hbb_lambda == 530), 'filled');
for j = 1:size(stations,1)-2
  m_scatter(stations.lon(j), stations.lat(j), 200, 'x', ...
    stations.bbp_mark_color{j},'LineWidth', 1)
  m_text(stations.lon(j)+stations.text_xpos(j), stations.lat(j)+stations.text_ypos(j), stations.st{j}, 'Color', ...
    stations.bbp_text_color{j}, 'FontSize', 15)
end
title(sprintf('Gayoso leg, %i days satellite average: %s - %s', ...
  round(days(max(stations.dt_start) - min(stations.dt_start))), ...
  datetime(stations.dt_start(1),'yyyy-mm-dd'), ...
  datetime(stations.dt_end(12),'yyyy-mm-dd')), 'FontSize', 13)
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', ...
  sprintf('%s_Gayoso_bbp_map_wholeleg', cruise)), 'jpg')
pause(1)
close figure 1

% map chl + acs chl
[~, cb] = QuickMap(month_chl, lat, lon, 'morgenstemning', [0 7], 'log'); % morgenstemning 0.00005 max(bbp550(:))
cb.Label.String = '[Chl] [mg.m^{-3}]';
hold on
% plot Inline
m_scatter(hbb.lon, hbb.lat, 100, hbb.chl, ...
  'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
m_scatter(hbb.lon, hbb.lat, 100, hbb.chl, ...
  'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
m_scatter(hbb.lon, hbb.lat, 100, hbb.chl, 'filled');
for j = 1:size(stations,1)-2
  m_scatter(stations.lon(j), stations.lat(j), 200, 'x', ...
    stations.chl_mark_color{j},'LineWidth', 1)
  m_text(stations.lon(j)+stations.text_xpos(j), stations.lat(j)+stations.text_ypos(j), stations.st{j}, 'Color', ...
    stations.chl_text_color{j}, 'FontSize', 15)
end
title(sprintf('Gayoso leg, %i days satellite average: %s - %s', ...
  round(days(max(stations.dt_start) - min(stations.dt_start))), ...
  datetime(stations.dt_start(1),'yyyy-mm-dd'), ...
  datetime(stations.dt_end(12),'yyyy-mm-dd')), 'FontSize', 13)
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', ...
  sprintf('%s_Gayoso_chl_map_wholeleg', cruise)), 'jpg')
pause(1)
close figure 1

% map bbp/chl + inline bbp/chl
month_bbp_month_chl = month_bbp./month_chl;

% month_bbp_month_chl(month_bbp < 0.003) = 0;
[~, cb] = QuickMap(month_bbp_month_chl, lat, lon, 'morgenstemning', [0.001 0.015], 'log'); % morgenstemning 0.00005 max(bbp550(:))
cb.Label.String = 'bbp550/[Chl]';
hold on
% plot Inline
m_scatter(hbb.lon, hbb.lat, 100, hbb.bbp(:, hbb_lambda == 530)./hbb.chl, ...
  'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
m_scatter(hbb.lon, hbb.lat, 100, hbb.bbp(:, hbb_lambda == 530)./hbb.chl, ...
  'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
m_scatter(hbb.lon, hbb.lat, 100, hbb.bbp(:, hbb_lambda == 530)./hbb.chl, 'filled');
for j = 1:size(stations,1)-2
  m_scatter(stations.lon(j), stations.lat(j), 200, 'x', ...
    stations.chlbbp_mark_color{j},'LineWidth', 1)
  m_text(stations.lon(j)+stations.text_xpos(j), stations.lat(j)+stations.text_ypos(j), stations.st{j}, 'Color', ...
    stations.chlbbp_text_color{j}, 'FontSize', 15)
end
title(sprintf('Gayoso leg, %i days satellite average: %s - %s', ...
  round(days(max(stations.dt_start) - min(stations.dt_start))), ...
  datetime(stations.dt_start(1),'yyyy-mm-dd'), ...
  datetime(stations.dt_end(12),'yyyy-mm-dd')), 'FontSize', 13)
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', ...
  sprintf('%s_Gayoso_bbpchl_map_wholeleg', cruise)), 'svg')
pause(1)
close figure 1




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
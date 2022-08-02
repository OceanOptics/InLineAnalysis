% Consolidate daily processing, merge variables, and export to SeaBASS
% author: Guillaume Bourdin
% created: March 1, 2021

%% Import data
if ispc
  cd('C:\Users\Gui\Documents\MATLAB\InLineAnalysis\InLineAnalysis-master\');
elseif ismac
  cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master/')
end
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
load(fullfile(ila.instrument.TSG.path.prod, ['TaraChile&' cruise '_InLine_TSG_prod.mat']))

% Update cfg
ila.cfg.write.mode = 'One day one file';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'NMEA'}; % {'FLOW', 'TSG', 'BB31502','PAR', 'WSCD859'}
ila.cfg.days2run = datenum(2021,8,18):datenum(2022,7,15);

% populate ila.instrument
ila.Read('prod');

nmea = ila.instrument.NMEA.prod.a;
nmea.dt = datetime(nmea.dt, 'ConvertFrom', 'datenum');

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

%% TSG chilean leg and create lat/lon file Chile leg
% load lat lon vector
nav_data = readtable(fullfile(ila.instrument.FLOW.path.prod, 'TaraChile_LatLon_boat.csv'));
nav_data.Properties.VariableNames = {'dt','lat','lon'};
nav_data.dt = datetime(nav_data.dt, 'InputFormat', 'yyyy/MM/dd HH:mm');
tsg_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.a;
tsg_temp.dt = datetime(tsg_temp.dt,'ConvertFrom','datenum');

% [tsg2, index] = unique(tsg.dt); 
% tsg = tsg(index,:);
% yi = interp1(tsg2, y(index), tsg.dt);

% check for duplicats in flow data and delete
[~, L, ~] = unique(nav_data.dt,'first');
indexToDump = not(ismember(1:numel(nav_data.dt),L));
if sum(indexToDump) > 0
  fprintf('Warning: %i identical dates in instrument data => deleted\n', sum(indexToDump))
  nav_data(indexToDump, :) = [];
end

tsg_latlon_missing = interp1(nav_data.dt, [nav_data.lat nav_data.lon], tsg_temp.dt(isnan(tsg_temp.lat)));
tsg_temp.lat(isnan(tsg_temp.lat)) = tsg_latlon_missing(:, 1);
tsg_temp.lon(isnan(tsg_temp.lon)) = tsg_latlon_missing(:, 2);

latlon = [tsg_temp(:,1) tsg_temp(:,2) tsg_temp(:,5)];

% save latlon prod
fprintf('Export to mat and csv... ');
save(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    'TaraChile_latlon'), 'latlon');
writetable(latlon, fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    'TaraChile_latlon.csv'));
fprintf('Done\n');

% clean TSG table
tsg = tsg_temp;
tsg = removevars(tsg,{'lat_avg_sd','lat_avg_n','lon_avg_sd','lon_avg_n'});
  
if any(strcmp(tsg.Properties.VariableNames, 'tcal'))
  tsg = removevars(tsg,{'tcal_avg_n','c_avg_n','s_avg_n'});
  tsg = renamevars(tsg,'t_avg_n','avg_n');
  tsg = renamevars(tsg,'t','sst');
  tsg = renamevars(tsg,'t_avg_sd','sst_avg_sd');
  tsg = renamevars(tsg,'s','sss');
  tsg = renamevars(tsg,'s_avg_sd','sss_avg_sd');
  tsg = renamevars(tsg,'c','cond');
  tsg = renamevars(tsg,'c_avg_sd','cond_avg_sd');
else
  tsg = removevars(tsg,{'t1_avg_n','c1_avg_n','s_avg_n','sv_avg_n'});
  tsg = renamevars(tsg,'t2_avg_n','avg_n');
  tsg = renamevars(tsg,'t1','tcal');
  tsg = renamevars(tsg,'t1_avg_sd','tcal_avg_sd');
  tsg = renamevars(tsg,'t2','sst');
  tsg = renamevars(tsg,'t2_avg_sd','sst_avg_sd');
  tsg = renamevars(tsg,'s','sss');
  tsg = renamevars(tsg,'s_avg_sd','sss_avg_sd');
  tsg = renamevars(tsg,'c1','cond');
  tsg = renamevars(tsg,'c1_avg_sd','cond_avg_sd');
end
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
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
  'plots', 'TaraChile_TSG_timeseries'), 'jpg')
close figure 96

SimpleMap(tsg.sst, tsg(:,1:3), 'TSG SST [°C]')
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', 'TaraChile_TSG_SST_map'), 'jpg')
close figure 1
SimpleMap(tsg.sss, tsg(:,1:3), 'TSG SSS [PSU]')
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', 'TaraChile_TSG_SSS_map'), 'jpg')
close figure 1

% save TSG prod
fprintf('Export to mat and csv... ');
filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
  sprintf('%s_InLine_%s_%s_%s_Product_v%s', 'TaraChile', 'TSG', ...
  datestr(min(tsg.dt), 'yyyymmdd'), datestr(max(tsg.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));
save(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    filename), 'tsg');
writetable(tsg, fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    [filename '.csv']));
fprintf('Done\n');

% tsg = renamevars(tsg,'sst','t');
% tsg = renamevars(tsg,'sss','s');

% filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
%   sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', 'TaraChile', 'TSG', ...
%   datestr(min(tsg.dt), 'yyyymmdd'), datestr(max(tsg.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));

% % export product to SeaBASS format
% ila.meta.documents = 'TaraChile_TSG_ProcessingReport_V2.pdf';
% ila.meta.calibration_files = cell2mat(list_dev(i));
% exportSeaBASS(filename,...
%     ila.meta,...
%     tsg,...
%     {'', '', ''});
% sprintf('%s_InLine_%s_Product.sb saved', 'TaraChile', cell2mat(ila.cfg.instruments2run))

tsg_chile = tsg;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TSG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'TSG'}; % {'FLOW', 'TSG', 'BB31502','PAR', 'WSCD859'}
ila.cfg.days2run = datenum(2021,8,16):datenum(2022,7,9);

% populate ila.instrument
ila.Read('prod');

% interpolate SST / SSS / LatLon
tsg_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.a;
tsg_temp.dt = datetime(tsg_temp.dt,'ConvertFrom','datenum');

% build TSG table
if any(strcmp(tsg_temp.Properties.VariableNames, 'lat'))
  tsg = tsg_temp;
  tsg = removevars(tsg,{'lat_avg_sd','lat_avg_n','lon_avg_sd','lon_avg_n'});
  load(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    'TaraChile_latlon'))
else
  latlon_interp = interp1(nmea.dt, [nmea.lat, nmea.lon], tsg_temp.dt, 'linear', 'extrap');
  tsg = addvars(tsg_temp, latlon_interp(:,1), latlon_interp(:,2), 'After', 'dt',...
    'NewVariableNames',{'lat','lon'});
end
if any(strcmp(tsg.Properties.VariableNames, 'tcal'))
  tsg = removevars(tsg,{'tcal_avg_n','c_avg_n','s_avg_n'});
  tsg = renamevars(tsg,'t_avg_n','avg_n');
  tsg = renamevars(tsg,'t','sst');
  tsg = renamevars(tsg,'t_avg_sd','sst_avg_sd');
  tsg = renamevars(tsg,'s','sss');
  tsg = renamevars(tsg,'s_avg_sd','sss_avg_sd');
  tsg = renamevars(tsg,'c','cond');
  tsg = renamevars(tsg,'c_avg_sd','cond_avg_sd');
else
  tsg = removevars(tsg,{'t1_avg_n','c1_avg_n','s_avg_n','sv_avg_n'});
  tsg = renamevars(tsg,'t2_avg_n','avg_n');
  tsg = renamevars(tsg,'t1','tcal');
  tsg = renamevars(tsg,'t1_avg_sd','tcal_avg_sd');
  tsg = renamevars(tsg,'t2','sst');
  tsg = renamevars(tsg,'t2_avg_sd','sst_avg_sd');
  tsg = renamevars(tsg,'s','sss');
  tsg = renamevars(tsg,'s_avg_sd','sss_avg_sd');
  tsg = renamevars(tsg,'c1','cond');
  tsg = renamevars(tsg,'c1_avg_sd','cond_avg_sd');
end
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
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
  'plots', [cruise '_TSG_timeseries']), 'jpg')
close figure 96

SimpleMap(tsg.sst, tsg(:,1:3), 'TSG SST [°C]')
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', [cruise '_TSG_SST_map']), 'jpg')
close figure 1
SimpleMap(tsg.sss, tsg(:,1:3), 'TSG SSS [PSU]')
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', [cruise '_TSG_SSS_map']), 'jpg')
close figure 1

% filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
%   sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, 'TSG', ...
%   datestr(min(tsg.dt), 'yyyymmdd'), datestr(max(tsg.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));

% % export product to SeaBASS format
% ila.meta.documents = [cruise '_TSG_ProcessingReport_V2.pdf';
% ila.meta.calibration_files = cell2mat(list_dev(i));
% exportSeaBASS(filename,...
%     ila.meta,...
%     tsg,...
%     {'', '', ''});
% sprintf('%s_InLine_%s_Product.sb saved', cruise, cell2mat(ila.cfg.instruments2run))

% save TSG prod
fprintf('Export to mat and csv... ');
filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
  sprintf('%s_InLine_%s_%s_%s_Product_v%s', cruise, 'TSG', ...
  datestr(min(tsg.dt), 'yyyymmdd'), datestr(max(tsg.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));
save(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    filename), 'tsg');
writetable(tsg, fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    [filename '.csv']));
fprintf('Done\n');

%% merge Chile and Microbiome TSG
tsg_chile = addvars(tsg_chile, NaN(size(tsg_chile, 1), 1), 'After', 'sss_avg_sd', 'NewVariableNames', 'sv');
tsg_chile = addvars(tsg_chile, NaN(size(tsg_chile, 1), 1), 'After', 'sv', 'NewVariableNames', 'sv_avg_sd');
tsg = [tsg_chile; tsg];

SimpleMap(tsg.sst, tsg(:,1:3), 'TSG SST [°C]')
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
  'plots', ['TaraChile&' cruise '_TSG_SST_map']), 'jpg')
close figure 1
SimpleMap(tsg.sss, tsg(:,1:3), 'TSG SSS [PSU]')
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
  'plots', ['TaraChile&' cruise '_TSG_SSS_map']), 'jpg')
close figure 1

fprintf('Export to mat and csv... ');
filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
  sprintf('%s_InLine_%s_%s_%s_Product_v%s', ['TaraChile&' cruise], 'TSG', ...
  datestr(min(tsg.dt), 'yyyymmdd'), datestr(max(tsg.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));
save(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    filename), 'tsg');
writetable(tsg, fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    [filename '.csv']));
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
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
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
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
  'plots', [cruise '_gayoso_zoomTSplot']), 'jpg')
close figure 1

%% T/S scatter plot for Namibian leg
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
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
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
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
  'plots', [cruise '_gayoso_zoomTSplot']), 'jpg')
close figure 1

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WSCD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'WSCD859'}; % {'FLOW', 'TSG', 'BB31502','PAR', 'WSCD859'}
ila.cfg.days2run = datenum(2020,12,24):datenum(2021,5,9);

% populate ila.instrument
ila.Read('prod');

% interpolate SST / SSS / LatLon
wscd_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.pd;
wscd_temp.dt = datetime(wscd_temp.dt,'ConvertFrom','datenum');
% latlon_interp = interp1(latlon.dt, [latlon.lat, latlon.lon], fdom_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.sst, tsg.sss], wscd_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data

% FDOM Products
wscd = table(wscd_temp.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4), wscd_temp.fdom, wscd_temp.fdom_sd, wscd_temp.fdom_n,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'fdom','fdom_sd','bincount'});
wscd.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'ppb', 'ppb', 'none'};
wscd.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f', '%.2f'};

% sort by date
wscd = sortrows(wscd, 'dt');

ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
  'plots', [cruise '_FDOM_timeseries']), 'jpg')
close figure 94

SimpleMap(wscd.fdom, wscd(:,1:3), 'WSCD fdom ppb')
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', [cruise '_WSCD_fdom_map']), 'jpg')
close figure 1

% filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
%   sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, ila.cfg.instruments2run{:}, ...
%   datestr(min(wscd.dt), 'yyyymmdd'), datestr(max(wscd.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));

% % export product to SeaBASS format
% ila.meta.documents = [cruise '_WSCD_ProcessingReport_V2.pdf';
% ila.meta.calibration_files = cell2mat(list_dev(i));
% exportSeaBASS(filename,...
%     ila.meta,...
%     wscd,...
%     {'', '', ''});
% sprintf('%s_InLine_%s_Product.sb saved', cruise, cell2mat(ila.cfg.instruments2run))

wscd.Properties.VariableNames = {'dt', 'lat', 'lon', 'sst', 'sss','fdom','fdom_sd','fdom_n'};

% save WSCD prod
fprintf('Export to mat and csv... ');
filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
  sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, ila.cfg.instruments2run{:}, ...
  datestr(min(wscd.dt), 'yyyymmdd'), datestr(max(wscd.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));
save(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    filename), 'wscd');
writetable(wscd, fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    [filename '.csv']));
fprintf('Done\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUVF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'SUVF'}; % {'FLOW', 'TSG', 'BB31502','PAR', 'WSCD859'}
% ila.cfg.days2run = datenum(2021,08,15):datenum(2022,01,7);

% populate ila.instrument
ila.Read('prod');

% interpolate SST / SSS / LatLon
suvf_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.pd;
suvf_temp.dt = datetime(suvf_temp.dt,'ConvertFrom','datenum');
latlon_interp = interp1(nmea.dt, [nmea.lat, nmea.lon], suvf_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data
tsg_interp = interp1(tsg.dt, [tsg.sst, tsg.sss], suvf_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data

% FDOM Products
suvf = table(suvf_temp.dt, latlon_interp(:,1), latlon_interp(:,2), tsg_interp(:,1), tsg_interp(:,2), suvf_temp.fdom, suvf_temp.fdom_sd, suvf_temp.fdom_n,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'fdom','fdom_sd','bincount'});
suvf.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'ppb', 'ppb', 'none'};
suvf.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f', '%.2f'};

% sort by date
suvf = sortrows(suvf, 'dt');

ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
  'plots', [cruise '_FDOM_timeseries']), 'jpg')
close figure 94

SimpleMap(suvf.fdom, suvf(:,1:3), 'SUVF fdom ppb')
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
  'plots', [cruise '_SUVF_fdom_map']), 'jpg')
close figure 1

filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
  sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, ila.cfg.instruments2run{:}, ...
  datestr(min(suvf.dt), 'yyyymmdd'), datestr(max(suvf.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));

% export product to SeaBASS format
ila.meta.documents = [cruise '_SUVF_ProcessingReport_V2.pdf'];
ila.meta.calibration_files = 'SUVF6244_CharSheet.pdf';
exportSeaBASS(filename,...
    ila.meta,...
    suvf,...
    {'', '', ''});
sprintf('%s_InLine_%s_Product.sb saved', cruise, cell2mat(ila.cfg.instruments2run))

suvf.Properties.VariableNames = {'dt', 'lat', 'lon', 'sst', 'sss','fdom','fdom_sd','fdom_n'};

% save SUVF prod
filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
  sprintf('%s_InLine_%s_%s_%s_Product_v%s', cruise, ila.cfg.instruments2run{:}, ...
  datestr(min(suvf.dt), 'yyyymmdd'), datestr(max(suvf.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));
fprintf('Export to mat and csv... ');
save(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    filename), 'suvf');
writetable(suvf, fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    [filename '.csv']));
fprintf('Done\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'PAR'}; % {'FLOW', 'TSG', 'BB3','PAR', 'WSCD1082P'}
ila.cfg.days2run = datenum(2020,12,24):datenum(2021,5,9);

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

ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod, ...
  'plots', [cruise '_PAR_timeseries']), 'jpg')
close figure 1

filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
  sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, ila.cfg.instruments2run{:}, ...
  datestr(min(par.dt), 'yyyymmdd'), datestr(max(par.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));

% export product to SeaBASS format
ila.meta.documents = [cruise '_PAR_ProcessingReport_V2.pdf'];
ila.meta.calibration_files = 'PAR-50168_CalSheet.pdf';
exportSeaBASS(filename,...
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
filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
  sprintf('%s_InLine_%s_%s_%s_Product_v%s', cruise, ila.cfg.instruments2run{:}, ...
  datestr(min(par.dt), 'yyyymmdd'), datestr(max(par.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));
save(fullfile(ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod, ...
    filename), 'par');
writetable(par, fullfile(ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod, ...
    [filename '.csv']));
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
% latlon_interp = interp1(latlon.dt, [latlon.lat, latlon.lon], bb3_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.sst, tsg.sss], bb3_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data

% bb3 Products
bb3 = table(bb3_temp.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4), bb3_temp.betap, bb3_temp.bbp, bb3_temp.betap_sd, bb3_temp.betap_n,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'VSF_124ang','bbp','VSF_124ang_sd','bincount'});
bb3.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', '1/m/sr', '1/m', '1/m/sr', 'none'};
bb3.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f'};

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
ila.DiagnosticPlot('BB', {'prod'}, save_figures); % AC or BB
close all

% plot time series
ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_BB3_BBparticulate_timeseries']), 'jpg')
close figure 86
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_BB3_POCparticulate_timeseries']), 'jpg')
close figure 85
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_BB3_BBdissolved_timeseries']), 'jpg')
close figure 87

SimpleMap(bb3.bbp(:,2), bb3(:,1:3), 'bbp (532 nm) [m^-^1]')
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_bb3_bbp_map']), 'jpg')
close figure 1

filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
  sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, ila.cfg.instruments2run{:}, ...
  datestr(min(bb3.dt), 'yyyymmdd'), datestr(max(bb3.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));

% export product to SeaBASS format
ila.meta.documents = sprintf('%s_BB3_ProcessingReport_v%s.pdf', cruise, datestr(now, 'yyyymmdd'));
ila.meta.calibration_files = 'BB3-1502_(470-532-650nm)_CharSheet.pdf';
exportSeaBASS(filename,...
    ila.meta,...
    bb3,...
    {string(bb3_lambda),string(bb3_lambda),string(bb3_lambda),''});
fprintf('%s_InLine_%s_Particulate_v%s.sb saved\n', cruise, ila.cfg.instruments2run{:}, ...
  datestr(now, 'yyyymmdd'))

bb3.Properties.VariableNames = {'dt', 'lat', 'lon', 'sst', 'sss', 'VSF124', 'bbp', 'VSF124_sd','bincount'};
bb3.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'm^-1.sr^-1', 'm^-1', 'm^-1.sr^-1', 'none'};
bb3.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f'};

% save BB3 prod
fprintf('Export to mat and csv... ');
filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
  sprintf('%s_InLine_%s_%s_%s_Product_v%s', cruise, ila.cfg.instruments2run{:}, ...
  datestr(min(bb3.dt), 'yyyymmdd'), datestr(max(bb3.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));
save(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    filename), 'bb3', 'bb3_lambda');
writetable(bb3, fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    [filename '.csv']));
fprintf('Done\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HBB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'HBB'}; % {'FLOW', 'TSG', 'BB31502','PAR', 'WSCD1082', 'HBB'}
ila.cfg.days2run = datenum(2021,8,16):datenum(2022,7,9);

% populate ila.instrument
ila.Read('prod');

% interpolate SST / SSS / LatLon
hbb_lambda = ila.instrument.(ila.cfg.instruments2run{:}).lambda;
hbb_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.p;
hbb_temp.dt = datetime(hbb_temp.dt,'ConvertFrom','datenum');
latlon_interp = interp1(nmea.dt, [nmea.lat, nmea.lon], hbb_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data
tsg_interp = interp1(tsg.dt, [tsg.sst, tsg.sss], hbb_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data

% bb3 Products
hbb = table(hbb_temp.dt, latlon_interp(:,1), latlon_interp(:,2), tsg_interp(:,1), tsg_interp(:,2), ...
  hbb_temp.betap, hbb_temp.bbp, hbb_temp.betap_sd, hbb_temp.betap_n, 'VariableNames', ...
  {'dt', 'lat', 'lon', 't', 's', 'VSF_124ang','bbp','VSF_124ang_sd','bincount'});
hbb.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', '1/m/sr', '1/m', '1/m/sr', 'none'};
hbb.Properties.VariableDescriptions = {'', '%.4f', '%.4f','%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f'};

% sort by date
hbb = sortrows(hbb, 'dt'); 
hbb.bbp(hbb.bbp<0)=NaN;
hbb.VSF_124ang_sd(hbb.VSF_124ang<0)=NaN;
hbb.VSF_124ang(hbb.VSF_124ang<0)=NaN;
hbb(all(isnan(hbb.bbp),2),:)=[];

%%% BB 3D plots %%%
save_figures = true;
ila.DiagnosticPlot('BB', {'prod'}, save_figures); % AC or BB
close all

% plot time series
ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_HBB_BBparticulate_timeseries']), 'jpg')
close figure 86
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_HBB_POCparticulate_timeseries']), 'jpg')
close figure 85
% saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_HBB_BBdissolved_timeseries']), 'jpg')
% close figure 87

SimpleMap(hbb.bbp(:,hbb_lambda == 530), hbb(:,1:3), 'bbp (530 nm) [m^-^1]')
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_hbb_bbp_map']), 'jpg')
close figure 1

filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
  sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, ila.cfg.instruments2run{:}, ...
  datestr(min(hbb.dt), 'yyyymmdd'), datestr(max(hbb.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));

% export product to SeaBASS format
ila.meta.documents = sprintf('%s_HBB_ProcessingReport_v%s.pdf', cruise, datestr(now, 'yyyymmdd'));
ila.meta.documents = [cruise '_HBB_ProcessingReport.pdf'];
ila.meta.calibration_files = 'HBB8005_CharSheet.pdf';
exportSeaBASS(filename,...
    ila.meta,...
    hbb,...
    {string(hbb_lambda),string(hbb_lambda),string(hbb_lambda),''});
fprintf('%s_InLine_%s_Particulate_v%s.sb saved\n', cruise, ila.cfg.instruments2run{:}, ...
  datestr(now, 'yyyymmdd'))

hbb.Properties.VariableNames = {'dt', 'lat', 'lon', 'sst', 'sss', 'VSF124', 'bbp', 'VSF124_sd','bincount'};
hbb.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'm^-1.sr^-1', 'm^-1', 'm^-1.sr^-1', 'none'};
hbb.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f'};

% save HBB prod
fprintf('Export to mat and csv... ');
filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
  sprintf('%s_InLine_%s_%s_%s_Product_v%s', cruise, ila.cfg.instruments2run{:}, ...
  datestr(min(hbb.dt), 'yyyymmdd'), datestr(max(hbb.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));
save(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    filename), 'hbb', 'hbb_lambda');
writetable(hbb, fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, ...
    [filename '.csv']));
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
    };

list_dev = cellstr(repmat(fullfile(path_dev, 'acs057_20200129.dev'), size(list_leg)));
list_instru = cellstr(repmat('ACS57', size(list_leg)));

data_AC = struct('particulate', [], 'product', []);
acs = [];
if exist('bb3', 'var')
  bp_bb3 = [];
end
if exist('hbb', 'var')
  bp_hbb = [];
end

for i=1:size(list_instru,1)
  if ispc
    cd('C:\Users\Gui\Documents\MATLAB\InLineAnalysis\InLineAnalysis-master\');
  elseif ismac
    cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master/')
  end
  ila = InLineAnalysis(['cfg' filesep cruise '_cfg.m']);
  ila.cfg.instruments2run = list_instru(i);
  ila.cfg.days2run = list_leg{i};

  % populate ila.instrument
  ila.Read('prod');

  [ila.instrument.(list_instru{i}).lambda_c, ...
    ila.instrument.(list_instru{i}).lambda_a] = importACSDeviceFile(list_dev{i});

  % interpolate SST / SSS / LatLon
  ref = [list_instru{i} '_' datestr(ila.cfg.days2run(1),'yyyymmdd') '_' datestr(ila.cfg.days2run(end),'yyyymmdd')];
  AC = ila.instrument.(list_instru{i}).prod.p;
  AC.dt = datetime(AC.dt,'ConvertFrom','datenum');
  latlon_interp = interp1(nmea.dt, [nmea.lat, nmea.lon], AC.dt, 'linear'); % extrap needed for first minute of data
  if all(isnan(latlon_interp(:)))
    latlon_interp = interp1(tsg.dt, [tsg.lat, tsg.lon], AC.dt, 'linear'); % extrap needed for first minute of data
  end
%   latlon_interp = interp1(latlon.dt, [latlon.lat, latlon.lon], AC.dt, 'linear', 'extrap'); % extrap needed for first minute of data
%   tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.sst, tsg.sss], AC.dt, 'linear', 'extrap'); % extrap needed for first minute of data
  tsg_interp = interp1(tsg.dt, [tsg.sst, tsg.sss], AC.dt, 'linear', 'extrap'); % extrap needed for first minute of data

  % ACS Particulate for SeaBASS
%   data_AC.particulate.(ref) = table(AC.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4), AC.ap, AC.ap_sd, AC.cp, AC.cp_sd, AC.cp_n,...
%                'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'ap', 'ap_sd', 'cp', 'cp_sd', 'bincount'});
%   data_AC.particulate.(ref).Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', '1/m', '1/m', '1/m', '1/m', 'none'};
%   data_AC.particulate.(ref).Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.0f'};

  data_AC.particulate.(ref) = table(AC.dt, latlon_interp(:,1), latlon_interp(:,2), ...
    tsg_interp(:,1), tsg_interp(:,2), AC.ap, AC.ap_sd, AC.cp, AC.cp_sd, AC.cp_n,...
               'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'ap', 'ap_sd', 'cp', 'cp_sd', 'bincount'});
  data_AC.particulate.(ref).Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', '1/m', '1/m', '1/m', '1/m', 'none'};
  data_AC.particulate.(ref).Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.0f'};

  % ACS Products
%   data_AC.product.(ref) = table(AC.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4),...
%                'VariableNames', {'dt', 'lat', 'lon', 't', 's'});
%   data_AC.product.(ref) = [data_AC.product.(ref) AC(:, 8:end)];
%   data_AC.product.(ref).Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', ...
%               '1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m', ...
%               'ug/L','unitless','ug/L','unitless','unitless','ug/L','microns','??'};

  data_AC.product.(ref) = table(AC.dt, latlon_interp(:,1), latlon_interp(:,2), tsg_interp(:,1), tsg_interp(:,2),...
               'VariableNames', {'dt', 'lat', 'lon', 't', 's'});
  data_AC.product.(ref) = [data_AC.product.(ref) AC(:, 8:end) AC(:, 6:7)];
  data_AC.product.(ref).Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', ...
              '1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m', ...
              'ug/L','unitless','ug/L','unitless','unitless','ug/L','microns','unitless','unitless','unitless','unitless'};
  data_AC.product.(ref).Properties.VariableDescriptions = [{''} repmat({'%.4f'}, 1, size(data_AC.product.(ref),2) - 4) ...
    repmat({'%i'}, 1, 3)];

  data_AC.product.(ref).Properties.VariableNames{...
    strcmp(data_AC.product.(ref).Properties.VariableNames, 'poc')} = 'POC_cp';
  data_AC.product.(ref).Properties.VariableNames{...
    strcmp(data_AC.product.(ref).Properties.VariableNames, 'gamma')} = 'cp_gamma';
  data_AC.product.(ref).Properties.VariableNames{...
    strcmp(data_AC.product.(ref).Properties.VariableNames, 'chl_ap676lh')} = 'Chl_lineheight';

  % Remove NaN and aberrant
  data_AC.particulate.(ref)(any(isnan([data_AC.particulate.(ref).lat data_AC.particulate.(ref).lon]),2), :) = [];
  data_AC.product.(ref).POC_cp(data_AC.product.(ref).POC_cp < 0) = NaN;
  data_AC.product.(ref).Chl_lineheight(data_AC.product.(ref).Chl_lineheight < 0) = NaN;
%   data_AC.product.(ref).cp_gamma(data_AC.product.(ref).cp_gamma < 0) = NaN;
  data_AC.product.(ref)(all(isnan(table2array(data_AC.product.(ref)(:,6:8))),2),:)=[];
  data_AC.product.(ref)(any(isnan([data_AC.product.(ref).lat data_AC.product.(ref).lon]),2), :) = [];
  
  % interpolate over BB3 wavelength to QC/QA with bbp/bp ratio
  bp = data_AC.particulate.(ref).cp - data_AC.particulate.(ref).ap;
  try
    if exist('bp_bb3', 'var')
      bp_bb3 = [bp_bb3; datenum(data_AC.particulate.(ref).dt) interp1(ila.instrument.(list_instru{i}).lambda_a, ...
        bp',bb3_lambda,'linear')'];
    end
  catch
    if exist('bp_hbb', 'var')
      bp_hbb = [bp_hbb; datenum(data_AC.particulate.(ref).dt) interp1(ila.instrument.(list_instru{i}).lambda_a, ...
        bp',hbb_lambda,'linear')'];
    end
  end
  
  %%% AC 3D plots %%%
  save_figures = false;
  ila.DiagnosticPlot('AC', {'prod'}, save_figures); % AC or BB
  close all
  
  filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
    sprintf('%s_InLine_%s_Particulate_v%s.sb', cruise, ref, datestr(now, 'yyyymmdd')));
  % export particulate to SeaBASS format
  ila.meta.documents = sprintf('%s_ACS_ProcessingReport_v%s.pdf', cruise, datestr(now, 'yyyymmdd'));
  [~, calfile] = fileparts(list_dev{i});
  ila.meta.calibration_files = [calfile '.dev'];
  exportSeaBASS(filename,...
      ila.meta,...
      data_AC.particulate.(ref),...
      {string(ila.instrument.(list_instru{i}).lambda_a),...
      string(ila.instrument.(list_instru{i}).lambda_a),...
      string(ila.instrument.(list_instru{i}).lambda_c),...
      string(ila.instrument.(list_instru{i}).lambda_c),''});
  fprintf('%s saved\n', filename)
  
  lambda.a = ila.instrument.(list_instru{i}).lambda_a;
  lambda.c = ila.instrument.(list_instru{i}).lambda_c;
  acs_part = data_AC.particulate.(ref);
  
  save(fullfile(ila.instrument.(list_instru{i}).path.prod, ...
    sprintf('%s_InLine_%s_Particulate_v%s', cruise, ref, datestr(now, 'yyyymmdd'))), ...
    'acs_part', 'lambda');

  % ACS merged prod
  acs = [acs; data_AC.product.(ref)];
end

% sort AC and save prod file .mat
[~,b] = sort(acs.dt); % sort by date
acs_prod = acs(b,:);

if any(strcmp(acs_prod.Properties.VariableNames, 'chl_ap676lh'))
  acs_prod = renamevars(acs_prod, 'chl_ap676lh', 'Chl_lineheight');
end
if any(strcmp(acs.Properties.VariableNames, 'poc'))
  acs_prod = renamevars(acs_prod, 'poc', 'POC_cp');
end
if any(strcmp(acs.Properties.VariableNames, 'gamma'))
  acs_prod = renamevars(acs_prod, 'gamma', 'cp_gamma');
end

flag = read_flagbit(acs_prod.flag_bit, 'AC');

% remove flagged products
acs_prod.POC_cp(flag.poc_flag) = NaN;
acs_prod.Chl_lineheight(flag.chl_ap676lh_flag) = NaN;
acs_prod.cp_gamma(flag.gamma_flag) = NaN;
acs_prod.cp_gamma(flag.noisy600_650) = NaN;
acs_prod.chl_Halh(flag.chl_Halh_flag) = NaN;
acs_prod.HH_mphi(flag.HH_mphi_flag) = NaN;
acs_prod.HH_G50(flag.HH_G50_flag) = NaN;
acs_prod.chl_Halh(flag.chlratio_flag) = NaN;

% % remove suspicious products
% acs_prod.cp_gamma(flag.gamma_suspicious) = NaN;
% % acs_prod.cp_gamma(acs_prod.cp_gamma > 2) = NaN;
% acs_prod.POC_cp(flag.ap_bubbles) = NaN;
% acs_prod.POC_cp(flag.poc_suspicious) = NaN;
% acs_prod.Chl_lineheight(flag.chl_ap676lh_suspicious) = NaN;
% acs_prod.chl_Halh(flag.chl_Halh_suspicious) = NaN;
% acs_prod.HH_mphi(flag.HH_G50_mphi_suspicious) = NaN;
% acs_prod.HH_G50(flag.HH_G50_mphi_suspicious) = NaN;

figure()
subplot(2,3,1); histogram(acs_prod.POC_cp); xlabel('[POC] cp (mg.m^{-3})')
subplot(2,3,2); histogram(acs_prod.Chl_lineheight); xlabel('a_{p676} line height [chl a] (mg.m^{-3})')
subplot(2,3,3); histogram(acs_prod.cp_gamma); xlabel('gamma cp (mg.m^{-3})')
subplot(2,3,4); histogram(acs_prod.chl_Halh); xlabel('Houskeeper [chl] (mg.m^{-3})')
subplot(2,3,5); histogram(acs_prod.HH_mphi); xlabel('H&H phytoplankton slope size distribution')
subplot(2,3,6); histogram(acs_prod.HH_G50); xlabel('H&H phytoplankton G50: cross-sectional area (\mum)')

ila.instrument.(list_instru{i}).prod.p = acs_prod;

ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.(list_instru{i}).path.prod, 'plots', [cruise '_ACS_prod_timeseries']), 'jpg')
close figure 78
saveGraph(fullfile(ila.instrument.(list_instru{i}).path.prod, 'plots', [cruise '_ACS_prod_regressions']), 'jpg')
close figure 77

SimpleMap(acs_prod.chl_Halh, acs_prod(:,1:3), 'Houskeeper [chl] (mg.m^{-3})')
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_ACS_chl_Houskeeper_map']), 'jpg')
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_ACS_chl_Houskeeper_map']), 'jpg')
close figure 1

SimpleMap(acs_prod.HH_G50, acs_prod(:,1:3), 'H&H phytoplankton G50: cross-sectional area (\mum)')
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_ACS_H&H_map']), 'jpg')
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_ACS_H&H_map']), 'fig')
close figure 1

SimpleMap(acs_prod.POC_cp, acs_prod(:,1:3), '[POC] cp (mg.m^{-3})')
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_ACS_POC_map']), 'jpg')
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_ACS_POC_map']), 'fig')
close figure 1

SimpleMap(acs_prod.Chl_lineheight, acs_prod(:,1:3), 'a_{p676} line height [chl a] (mg.m^{-3})')
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_ACS_chl_map']), 'jpg')
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_ACS_chl_map']), 'fig')
close figure 1

SimpleMap(acs_prod.cp_gamma, acs_prod(:,1:3), 'gamma cp (unitless)')
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_ACS_gamma_map']), 'jpg')
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_ACS_gamma_map']), 'fig')
close figure 1

acs_prod.chl_Halh = [];
acs_prod.HH_G50 = [];

% save AC prod
fprintf('Export to mat and csv... ');
save(fullfile(ila.instrument.(list_instru{i}).path.prod, [cruise '_InLine_ACS_prod_full']), 'acs_prod');
writetable(acs_prod, fullfile(ila.instrument.(list_instru{i}).path.prod, [cruise '_InLine_ACS_prod_full.csv']));
fprintf('Done\n');

% keep only old variables
acs_prod = acs_prod(:, contains(acs_prod.Properties.VariableNames, {'dt', 'lat', 'lon', ...
  't','s','Chl_lineheight','POC_cp','cp_gamma','ap_n','cp_n','flag_bit'}) & ...
  ~contains(acs_prod.Properties.VariableNames, {'agaus'}));

% save AC prod
fprintf('Export to mat and csv... ');
filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
  sprintf('%s_InLine_ACS_%s_%s_Product_v%s', cruise, datestr(min(acs_prod.dt), 'yyyymmdd'), ...
  datestr(max(acs_prod.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));
save(fullfile(ila.instrument.(list_instru{i}).path.prod, filename), 'acs_prod');
writetable(acs_prod, fullfile(ila.instrument.(list_instru{i}).path.prod, [filename '.csv']));
fprintf('Done\n');

% export product to SeaBASS format
filename = fullfile(ila.instrument.(list_instru{i}).path.prod, ...
  sprintf('%s_InLine_ACS_%s_%s_Product_v%s.sb', cruise, datestr(min(acs_prod.dt), 'yyyymmdd'), ...
  datestr(max(acs_prod.dt), 'yyyymmdd'), datestr(now, 'yyyymmdd')));
exportSeaBASS(filename,...
    ila.meta,...
    acs_prod);
fprintf('%s saved\n', filename)

%%
acs_part(~ismember(acs_part.dt, acs.dt), :) = [];
bp_div_chl = (acs_part.cp(:,40) - acs_part.ap(:,40)) ./ acs.Chl_lineheight;

SimpleMap(bp_div_chl, acs(:,1:3), 'b_p / [chl a]')
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_ACS_bp_div_chl_map']), 'jpg')
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_ACS_bp_div_chl_map']), 'fig')
close figure 1

%%
acs_part(~ismember(acs_part.dt, acs.dt), :) = [];

SimpleMap(acs_part.cp(:,40) - acs_part.ap(:,40), acs(:,1:3), 'b_p [m^{-1}')
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_ACS_bp_map']), 'jpg')
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_ACS_bp_map']), 'fig')
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
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_BB3_bbp_over_bp']), 'jpg')
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
saveGraph(fullfile(ila.instrument.TSG.path.prod, 'plots', [cruise '_HBB_bbp_over_bp']), 'jpg')
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
  title(sprintf('Station %s %s', stations.st{i}, datestr(stations.dt_start(i),'yyyy-mm-dd')), 'FontSize', 13)
  saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', ...
    sprintf('%s_Gayoso_bbp_map_station_%s_%s', cruise, stations.st{i}, ...
    datestr(stations.dt_start(i),'yyyymmdd'))), 'jpg')
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
  title(sprintf('Station %s %s', stations.st{i}, datestr(stations.dt_start(i),'yyyy-mm-dd')), 'FontSize', 13)
  saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', ...
    sprintf('%s_Gayoso_chl_map_station_%s_%s', cruise, stations.st{i}, ...
    datestr(stations.dt_start(i),'yyyymmdd'))), 'jpg')
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
  title(sprintf('Station %s %s', stations.st{i}, datestr(stations.dt_start(i),'yyyy-mm-dd')), 'FontSize', 13)
  saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', ...
    sprintf('%s_Gayoso_bbpchl_map_station_%s_%s', cruise, stations.st{i}, ...
    datestr(stations.dt_start(i),'yyyymmdd'))), 'jpg')
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
  datestr(stations.dt_start(1),'yyyy-mm-dd'), ...
  datestr(stations.dt_end(12),'yyyy-mm-dd')), 'FontSize', 13)
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', ...
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
  datestr(stations.dt_start(1),'yyyy-mm-dd'), ...
  datestr(stations.dt_end(12),'yyyy-mm-dd')), 'FontSize', 13)
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', ...
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
  datestr(stations.dt_start(1),'yyyy-mm-dd'), ...
  datestr(stations.dt_end(12),'yyyy-mm-dd')), 'FontSize', 13)
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', ...
  sprintf('%s_Gayoso_bbpchl_map_wholeleg', cruise)), 'svg')
pause(1)
close figure 1


%% leg Namibia-Congo map centered on stations
stations = readtable('/Volumes/Samsung_T5/Data/TaraMicrobiome/prod/stations_Namibia-Congo.csv');

% lagrangian #1 A061-A063
% lagrangian #2 A064-A067
% create table
sz_tbl = size(stations,1);
stations.dt_start = datetime(stations.dt_start, 'InputFormat', 'yyyy/MM/dd HH:mm');
stations.dt_end = datetime(stations.dt_end, 'InputFormat', 'yyyy/MM/dd HH:mm');
stations.station{sz_tbl+1} = 'Diurnal experiment day 1';
stations.dt_start(sz_tbl+1) = stations.dt_start(13);
stations.dt_end(sz_tbl+1) = stations.dt_end(16);
stations.lat(sz_tbl+1) = stations.lat(15);
stations.lon(sz_tbl+1) = stations.lon(15);

stations.station{sz_tbl+2}  = 'Diurnal experiment day 2';
stations.dt_start(sz_tbl+2) = stations.dt_start(17);
stations.dt_end(sz_tbl+2) = stations.dt_end(20);
stations.lat(sz_tbl+2) = stations.lat(19);
stations.lon(sz_tbl+2) = stations.lon(19);

stations.station{sz_tbl+3}  = 'Diurnal experiment day 3';
stations.dt_start(sz_tbl+3) = stations.dt_start(21);
stations.dt_end(sz_tbl+3) = stations.dt_end(24);
stations.lat(sz_tbl+3) = stations.lat(23);
stations.lon(sz_tbl+3) = stations.lon(23);

% create table
stations.id = cell(size(stations,1), 1);
stations.toplot = cell(size(stations,1), 1);
for i = 1:size(stations, 1)
  stations.id{i} = tsg.dt >= stations.dt_start(i) & tsg.dt <= stations.dt_end(i);
  if strcmp(stations.station{i}, 'Diurnal experiment day 1')
    stations.toplot{i} = {'D134-A','D134-B','D134-C','D134-D'};
  elseif strcmp(stations.station{i}, 'Diurnal experiment day 2')
    stations.toplot{i} = {'D134-E','D134-F','D134-G','D134-H'};
  elseif strcmp(stations.station{i}, 'Diurnal experiment day 3')
    stations.toplot{i} = {'D134-I','D134-J','D134-K','D134-L'};
  else
    stations.toplot{i} = stations.station{i};
  end
end

stations.bbp_text_color = {'white','black','white','black','white','white','black',...
  'black','black','black','black','white','white','black','black','black','black',...
  'black','black','black','black','black','black','black','black','black','black'}';
stations.chl_text_color = {'white','white','white','white','white','white','white',...
  'white','white','black','black','black','white','white','black','black','black',...
  'black','black','black','black','black','black','black','black','black','black'}';
stations.chlbbp_text_color = {'black','black','black','black','black','black','black',...
  'black','black','black','black','black','black','black','black','black','black',...
  'black','black','black','black','black','black','black','black','black','black'}';

% path_sat = '/Users/gui/Documents/Maine/Data/TaraPacific/remote_sensing/time_series/level3_timeseries/gayoso';
% list_sat = {dir(fullfile(path_sat, 'gayoso*mat')).name}';
% sat_start = datetime(cellfun(@(s) s{2}, cellfun(@(c) strsplit(c, '_'), list_sat, 'un', 0), 'un', 0), ...
%   'InputFormat', 'yyyyMMdd-HHmm');
% sat_end = datetime(cellfun(@(s) s{3}, cellfun(@(c) strsplit(c, '_'), list_sat, 'un', 0), 'un', 0), ...
%   'InputFormat', 'yyyyMMdd-HHmm');

% remove duplicates
[~, L, ~] = unique(hbb.dt,'first');
indexToDump = not(ismember(1:numel(hbb.dt), L));
if sum(indexToDump) > 0
  fprintf('Warning: %i identical dates in HBB data => deleted\n', sum(indexToDump))
  hbb(indexToDump, :) = [];
end

% remove duplicates
[~, L, ~] = unique(acs.dt,'first');
indexToDump = not(ismember(1:numel(acs.dt), L));
if sum(indexToDump) > 0
  fprintf('Warning: %i identical dates in HBB data => deleted\n', sum(indexToDump))
  acs(indexToDump, :) = [];
end

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
  title(sprintf('Station %s %s', stations.st{i}, datestr(stations.dt_start(i),'yyyy-mm-dd')), 'FontSize', 13)
  saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', ...
    sprintf('%s_Gayoso_bbp_map_station_%s_%s', cruise, stations.st{i}, ...
    datestr(stations.dt_start(i),'yyyymmdd'))), 'jpg')
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
  title(sprintf('Station %s %s', stations.st{i}, datestr(stations.dt_start(i),'yyyy-mm-dd')), 'FontSize', 13)
  saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', ...
    sprintf('%s_Gayoso_chl_map_station_%s_%s', cruise, stations.st{i}, ...
    datestr(stations.dt_start(i),'yyyymmdd'))), 'jpg')
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
  title(sprintf('Station %s %s', stations.st{i}, datestr(stations.dt_start(i),'yyyy-mm-dd')), 'FontSize', 13)
  saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', ...
    sprintf('%s_Gayoso_bbpchl_map_station_%s_%s', cruise, stations.st{i}, ...
    datestr(stations.dt_start(i),'yyyymmdd'))), 'jpg')
  pause(5)
  close figure 1
end

%% Whole leg Namibia-Congo
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
  datestr(stations.dt_start(1),'yyyy-mm-dd'), ...
  datestr(stations.dt_end(12),'yyyy-mm-dd')), 'FontSize', 13)
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', ...
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
  datestr(stations.dt_start(1),'yyyy-mm-dd'), ...
  datestr(stations.dt_end(12),'yyyy-mm-dd')), 'FontSize', 13)
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', ...
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
  datestr(stations.dt_start(1),'yyyy-mm-dd'), ...
  datestr(stations.dt_end(12),'yyyy-mm-dd')), 'FontSize', 13)
saveGraph(fullfile(ila.instrument.(ila.cfg.instruments2run{:}).path.prod, 'plots', ...
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
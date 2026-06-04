% Consolidate daily processing, merge variables, and export to SeaBASS
% author: Guillaume Bourdin
% created: Jan 30, 2024

%% Import data
cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master')

cruise = 'TaraEuropa';
% Load InLineAnalysis and the configuration
ila = InLineAnalysis(['cfg' filesep cruise '_cfg.m']);

path_dev = strrep(ila.instrument.FLOW.path.prod, ...
  'prod', 'DeviceFiles');

% create Graph folder if it doesn't exist
if ~isfolder(fullfile(ila.instrument.FLOW.path.prod, 'plots'))
  mkdir(fullfile(ila.instrument.FLOW.path.prod, 'plots'))
end

% whenever TSG is processed load with
%load([ila.instrument.FLOW.path.prod cruise '_InLine_TSG_prod.mat'])
% load([ila.instrument.SPCD.path.prod cruise '_InLine_SPCD_prod.mat'])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load lat lon vector from mergedGPScaptain
% load(fullfile(ila.instrument.FLOW.path.prod, 'TaraEuropa_merged_GPScaptain_v20240131.mat'))
% 
% ila.cfg.days2run = datetime(2023,4,3):datetime(2024,08,22);
% ila.instrument.NMEA.model = 'GPSSC701Tara'; % GPSSC701Tara GPS32Tara GPSCOMPASSAT
% ila.cfg.instruments2run = {'NMEA'};
% ila.Read('prod');
% GPSSC701Tara = ila.instrument.NMEA.prod.a;
% 
% latlon = table();
% latlon.dt = (min([merged_GPScaptain.dt; GPSSC701Tara.dt]):minutes(1):max([merged_GPScaptain.dt; GPSSC701Tara.dt]))';
% latlon = merge_timeseries(latlon, GPSSC701Tara, {'lat', 'lon'});
% latlon = merge_timeseries(latlon, merged_GPScaptain, {'lat', 'lon'});
% latlon = round_timestamp(latlon);
% 
% figure; hold on
% scatter(latlon.dt, latlon.lat, [], 'b')
% scatter(latlon.dt(isnan(latlon.lat) & ~isnan(merged_GPScaptain.lat)), merged_GPScaptain.lat(isnan(latlon.lat) & ~isnan(merged_GPScaptain.lat)), [], 'r')
% ylabel('Latitude')
% legend('GPSSC701', 'missing GPSSC701 filled with GPS captain')
% saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_GPSrecovery']), 'fig')
% close figure 1
% 
%  save(fullfile(ila.instrument.FLOW.path.prod, 'TaraEuropa_merged_GPSall_v20240906.mat'), 'latlon')
 load(fullfile(ila.instrument.FLOW.path.prod, 'TaraEuropa_merged_GPSall_v20240906.mat'))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TSG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ila.cfg.days2run = datetime(2023,4,3):datetime(2023,11,16);
% ila.cfg.days2run = datetime(2024,2,19):datetime(2024,8,22);
ila.cfg.days2run = datetime(2023,4,3):datetime(2024,8,22);
ila.cfg.instruments2run = {'SBE384504970269','SBE384504970286'};

% populate ila.instrument
ila.Read('prod');

% extract TSG data from obj
tsg_temp = [ila.instrument.('SBE384504970269').prod.a; ila.instrument.('SBE384504970286').prod.a]; % both years
% tsg_temp = ila.instrument.('SBE384504970269').prod.a;
tsg_temp = round_timestamp(tsg_temp);

% build TSG table: merge lat, lon
replace_consecutive_nan = 3*60; % 3h
tsg = merge_timeseries(tsg_temp, latlon, {'lat', 'lon'}, '', replace_consecutive_nan);
% Remove NaN
if any(any(isnan([tsg.lat tsg.lon]),2))
  warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
    sum(any(isnan([tsg.lat tsg.lon]),2)), replace_consecutive_nan/60)
  tsg(any(isnan([tsg.lat tsg.lon]),2), :) = [];
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

%%
load(fullfile(ila.instrument.FLOW.path.prod, 'TaraEuropa_InLine_TSG_20230403_20240822_Product_v20250710.mat'))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WSCD859 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ila.cfg.instruments2run = {'WSCD859'};
% ila.cfg.days2run = datetime(2020,12,24):datetime(2021,5,9);
% 
% % populate ila.instrument
% ila.Read('prod');
% 
% wscd_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.pd;
% % merge lat, lon, sst, sss
% replace_consecutive_nan = 72*60; % 72h
% % merge lat, lon, sst, sss
% replace_consecutive_nan = 72*60; % 72h
% wscd = merge_timeseries(wscd_temp, tsg, {'lat', 'lon', 'sst', 'sss'});
% wscd = merge_timeseries(wscd, nmea, {'lat', 'lon'});
% wscd = merge_timeseries(wscd, latlon, {'lat', 'lon'});
% wscd = merge_timeseries(wscd, meteo_ftp, {'lat', 'lon'}, '', replace_consecutive_nan);
% 
% % Remove NaN
% wscd(any(isnan([wscd.lat wscd.lon]),2), :) = [];
% if any(isnan([wscd.lat wscd.lon]),2)
%   wscd(any(isnan([wscd.lat wscd.lon]),2), :) = [];
%   warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
%     sum(any(isnan([wscd.lat wscd.lon]),2)), replace_consecutive_nan/60)
% end
% 
% % add units and precision
% wscd = renamevars(wscd, {'sst','sss','fdom_n'} ,{'t','s','bincount'});
% wscd.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'v_uncalibrated', 'v_uncalibrated', 'none'};
% wscd.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f', '%.2f'};
% 
% % sort by date
% wscd = sortrows(wscd, 'dt');
% 
% ila.visProd_timeseries()
% saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_FDOM_timeseries']), 'jpg')
% close figure 94
% 
% SimpleMap(wscd.fdom, wscd(:,1:3), 'WSCD fdom [v uncalibrated]')
% saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_WSCD_fdom_map']), 'jpg')
% close figure 1
% 
% filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, ila.cfg.instruments2run{:}, ...
%   datetime(min(wscd.dt), 'Format', 'yyyyMMdd'), datetime(max(wscd.dt), 'Format', 'yyyyMMdd'), ...
%   datetime('today', 'Format', 'yyyyMMdd'));
% 
% % % export product to SeaBASS format
% % ila.meta.documents = [cruise '_WSCD_ProcessingReport_V2.pdf';
% % ila.meta.calibration_files = cell2mat(list_dev(i));
% % exportSeaBASS(fullfile(ila.instrument.FLOW.path.prod, filename),...
% %     ila.meta,...
% %     wscd,...
% %     {'', '', ''});
% % sprintf('%s_InLine_%s_Product.sb saved', cruise, cell2mat(ila.cfg.instruments2run))
% 
% wscd = renamevars(wscd, {'t','s','bincount'}, {'sst','sss','fdom_n'});
% 
% % save WSCD prod
% fprintf('Export to mat and csv... ');
% filename = strrep(filename, '.sb', '');
% save(fullfile(ila.instrument.FLOW.path.prod, filename), 'wscd');
% writetable(wscd, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
% fprintf('Done\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ila.cfg.instruments2run = {'PAR'};
% ila.cfg.days2run = datetime(2020,12,24):datetime(2021,5,9);
% 
% % populate ila.instrument
% ila.Read('prod');
% 
% % interpolate SST / SSS / LatLon
% par_temp = ila.instrument.(cell2mat(ila.cfg.instruments2run)).prod.a;
% latlon_interp = interp1(latlon.dt, [latlon.lat, latlon.lon], par_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data
% tsg_interp = interp1(tsg.dt, [tsg.sst, tsg.sss_adj], par_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data
% 
% % par Products in uE/cm^2/s for SeaBASS
% par = table(par_temp.dt, latlon_interp(:,1), latlon_interp(:,2), tsg_interp(:,1), tsg_interp(:,2), par_temp.par./10000, par_temp.par_sd./10000, par_temp.par_n,...
%              'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'par','par_sd','bincount'});
% par.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'uE/cm^2/s', 'uE/cm^2/s', 'none'};
% par.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f', '%.2f'};
% 
% [~,b] = sort(par.dt); % sort by date
% par = par(b,:);
% 
% ila.visProd_timeseries()
% saveGraph([ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
%   'Graphs' filesep cruise '_PAR_timeseries'], 'jpg')
% close figure 1
%
% filename = [ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod cruise ...
%   '_InLine_' cell2mat(ila.cfg.instruments2run) '_Product_v' datestr(now, 'yyyymmdd') '.sb'];
%
% % export product to SeaBASS format
% ila.meta.documents = [cruise '_PAR_ProcessingReport_V2.pdf'];
% ila.meta.calibration_files = 'PAR-50168_CalSheet.pdf';
% exportSeaBASS(filename,...
%     ila.meta,...
%     par,...
%     {'', '', ''});
% sprintf('%s_InLine_%s_Product.sb saved', cruise, cell2mat(ila.cfg.instruments2run))
% 
% par.Properties.VariableNames = {'dt', 'lat', 'lon', 'sst', 'sss', 'par', 'par_sd', 'par_n'};
% 
% % convert to uE/m^2/s for mat file
% par.par = par.par.*10000;
% par.par_sd = par.par_sd.*10000;
% par.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'uE/m^2/s', 'uE/m^2/s', 'none'};
% 
% % save PAR prod
% fprintf('Export to mat and csv... ');
% save([ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
%     cruise '_InLine_PAR_prod'], 'par');
% writetable(par, [ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
%     cruise '_InLine_PAR_prod.csv']);
% fprintf('Done\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUVF6244 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila.cfg.instruments2run = {'SUVF6244'};
%ila.cfg.days2run = datetime(2023,4,3):datetime(2023,11,16);
ila.cfg.days2run = datetime(2023,4,3):datetime(2024,8,22);

% populate ila.instrument
ila.Read('prod');

% extract SUVF data from obj
suvf_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.pd;
suvf_temp = round_timestamp(suvf_temp);

% build suvf table: merge lat, lon, sst, sss
replace_consecutive_nan = 3*60; % 3h
suvf = merge_timeseries(suvf_temp, tsg, {'lat', 'lon', 'sst', 'sss'});
suvf = merge_timeseries(suvf, latlon, {'lat', 'lon'}, '', replace_consecutive_nan);

% Remove NaN
if any(isnan([suvf.lat suvf.lon]),2)
  suvf(any(isnan([suvf.lat suvf.lon]),2), :) = [];
  warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
    sum(any(isnan([suvf.lat suvf.lon]),2)), replace_consecutive_nan/60)
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

%%
load(fullfile(ila.instrument.FLOW.path.prod, 'TaraEuropa_InLine_SUVF6244_20230404_20240822_Product_v20250710.mat'))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HBB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila.cfg.instruments2run = {'HyperBB8005'};
ila.cfg.days2run = datetime(2023,4,3):datetime(2023,11,16);

% populate ila.instrument
ila.Read('prod');

% extract HyperBB data from obj
hbb_lambda = ila.instrument.(ila.cfg.instruments2run{:}).lambda;
hbb_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.p;
hbb_temp = round_timestamp(hbb_temp);

% build hbb table: merge lat, lon, sst, sss
replace_consecutive_nan = 3*60; % 3h
hbb_temp = merge_timeseries(hbb_temp, tsg, {'lat', 'lon', 'sst', 'sss'});
hbb_temp = merge_timeseries(hbb_temp, latlon, {'lat', 'lon'}, '', replace_consecutive_nan);

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
hbb = renamevars(hbb, {'sst', 'sss', 'betap','betap_se','betap_n'}, {'t','s','VSF_124ang','VSF_124ang_sd','bincount'});

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
close figure 26
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_HBB_BBparticulate_timeseries']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_HBB_BBparticulate_timeseries']), 'fig')
close figure 24
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_HBB_POCparticulate_timeseries']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_HBB_POCparticulate_timeseries']), 'fig')
close figure 23
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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_leg = {...
    datetime(2023,4,4):datetime(2023,4,26);...
    datetime(2023,5,3):datetime(2023,5,30);...
    datetime(2023,6,6):datetime(2023,6,30);...
    datetime(2023,7,8):datetime(2023,7,27);...
    datetime(2023,7,31):datetime(2023,8,17);...
    datetime(2023,8,19):datetime(2023,9,8);...
    datetime(2023,9,15):datetime(2023,10,11);...
    datetime(2023,10,18):datetime(2023,11,16);
    datetime(2024,2,19):datetime(2024,3,24);...
    datetime(2024,3,27):datetime(2024,4,4);...
    datetime(2024,4,7):datetime(2024,4,28);...
    datetime(2024,4,30):datetime(2024,5,28);...
    datetime(2024,5,31):datetime(2024,6,19);...
    datetime(2024,6,22):datetime(2024,7,10);...
    datetime(2024,7,17):datetime(2024,7,30);...
    datetime(2024,8,2):datetime(2024,8,22);
    };
% list_leg = {...
%     datetime(2024,2,19):datetime(2024,3,24);...
%     datetime(2024,3,27):datetime(2024,4,4);...
%     datetime(2024,4,7):datetime(2024,4,28);...
%     datetime(2024,4,30):datetime(2024,5,28);...
%     datetime(2024,5,31):datetime(2024,6,19);...
%     datetime(2024,6,22):datetime(2024,7,10);...
%     datetime(2024,7,17):datetime(2024,7,30);...
%     datetime(2024,8,2):datetime(2024,8,22);
%     };

list_dev = {'acs003_20220602.dev';
    'acs003_20220602.dev';
    'acs003_20220602.dev';
    'acs003_20220602.dev';
    'acs003_20220602.dev';
    'acs003_20220602.dev';
    'acs003_20220602.dev';
    'acs003_20220602.dev';
    'ACS348_20231013.dev';
    'ACS348_20231013.dev';
    'ACS348_20231013.dev';
    'ACS348_20231013.dev';
    'ACS348_20231013.dev';
    'ACS348_20231013.dev';
    'ACS348_20231013.dev';
    'ACS348_20231013.dev';
    };


% list_dev = repmat({fullfile(path_dev, 'ACS348_20231013.dev')}, size(list_leg, 1), 1);
% list_dev = [repmat({fullfile(path_dev, 'acs003_20220602.dev')}, size(list_leg, 1)-1, 1);...
%   {fullfile(path_dev, )}];

list_instru = {'ACS3';
    'ACS3';
    'ACS3';
    'ACS3';
    'ACS3';
    'ACS3';
    'ACS3';
    'ACS3';
    'ACS348';
    'ACS348';
    'ACS348';
    'ACS348';
    'ACS348';
    'ACS348';
    'ACS348';
    'ACS348';
    };

% list_instru = repmat({'ACS348'}, size(list_leg, 1), 1);
% list_instru = [repmat({'ACS3'}, size(list_leg, 1)-1, 1); {'ACS348'}];

data_AC = struct('particulate', [], 'product', []);
acs = [];

for i=1:size(list_instru,1)
  cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master')
  % cd('/Volumes/Data2/InLineAnalysis-master')
  ila = InLineAnalysis(['cfg' filesep cruise '_cfg.m']);
  ila.cfg.instruments2run = list_instru(i);
  ila.cfg.days2run = list_leg{i};

  % populate ila.instrument
  ila.Read('prod');
  
  % get wavelength
  [ila.instrument.(list_instru{i}).lambda_c, ...
    ila.instrument.(list_instru{i}).lambda_a] = importACSDeviceFile(fullfile(path_dev, list_dev{i}));

  % set section name
  ref = sprintf('%s_%s_%s', list_instru{i}, datetime(ila.cfg.days2run(1), 'Format', 'yyyyMMdd'), ...
    datetime(ila.cfg.days2run(end), 'Format', 'yyyyMMdd'));
  
  % remove flagged data
  flag = read_flagbit(ila.instrument.(list_instru{i}).prod.p.flag_bit, 'AC');
  % flag_info = InlineFlagInfo('ACS');

  % flag.HH_G50_flag(flag.HH_G50_flag & ila.instrument.(list_instru{i}).prod.p.HH_G50 > 0 & ila.instrument.(list_instru{i}).prod.p.HH_G50 <= 500) = false;
  % ila.instrument.(list_instru{i}).prod.p.flag_bit = set_flagbit(flag);
  
  % remove flagged products
  ila.instrument.(list_instru{i}).prod.p.poc(flag.poc_flag) = NaN;
  ila.instrument.(list_instru{i}).prod.p.chl_ap676lh(flag.chl_ap676lh_flag) = NaN;
  ila.instrument.(list_instru{i}).prod.p.gamma(flag.cp_step) = NaN;
  ila.instrument.(list_instru{i}).prod.p.poc(flag.cp_step) = NaN;
  ila.instrument.(list_instru{i}).prod.p.chl_Halh(flag.chl_Halh_flag) = NaN;
  ila.instrument.(list_instru{i}).prod.p.chl_Halh(flag.cp_step) = NaN;
  ila.instrument.(list_instru{i}).prod.p.HH_mphi(flag.cp_step) = NaN;
  ila.instrument.(list_instru{i}).prod.p.HH_G50(flag.HH_G50_flag) = NaN;
  ila.instrument.(list_instru{i}).prod.p.HH_G50(flag.cp_step) = NaN;
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
  close figure 14
  saveGraph(fullfile(ila.instrument.(list_instru{i}).path.prod, 'plots', [cruise '_' ref '_ACS_prod_timeseries']), 'jpg')
  saveGraph(fullfile(ila.instrument.(list_instru{i}).path.prod, 'plots', [cruise '_' ref '_ACS_prod_timeseries']), 'fig')
  close figure 10

  % extract AC data from obj
  AC = ila.instrument.(list_instru{i}).prod.p;

  % build AC table: merge lat, lon, sst, sss, and fdom
  replace_consecutive_nan = 3*60; % 3h
  AC = merge_timeseries(AC, suvf, {'fdom'}, '', replace_consecutive_nan);
  AC = merge_timeseries(AC, tsg, {'lat', 'lon', 'sst', 'sss'});
  AC = merge_timeseries(AC, latlon, {'lat', 'lon'}, '', replace_consecutive_nan);

  % rename variables for SeaBASS
  AC = renamevars(AC, {'sst','sss'}, {'t','s'});
  
  % Remove NaN
  if any(any(isnan([AC.lat AC.lon]),2))
    warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
      sum(any(isnan([AC.lat AC.lon]),2)), replace_consecutive_nan/60)
    AC(any(isnan([AC.lat AC.lon]),2), :) = [];
  end

  % split into particulate table
  id_particulate = logical(sum(categorical(AC.Properties.VariableNames) == {'dt','lat','lon','t','s','fdom','ap','ap_se','cp','cp_se','ap_n'}'));
  data_AC.particulate.(ref) = AC(:, id_particulate);
  % add variable units and description
  data_AC.particulate.(ref).Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU','ppb', '1/m', '1/m', '1/m', '1/m', 'none'};
  data_AC.particulate.(ref).Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.0f'};
  % rename variables for SeaBASS
  data_AC.particulate.(ref) = renamevars(data_AC.particulate.(ref), {'fdom','ap_n'}, {'cdomf','bincount'});

  % split into product table
  id_product = logical(sum(categorical(AC.Properties.VariableNames) == {'ap', 'ap_se', 'cp', 'cp_se'}'));
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
close figure 14
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

foo = acs_prod.cp_gamma > prctile(acs_prod.cp_gamma, 1) & acs_prod.cp_gamma < prctile(acs_prod.cp_gamma, 99);
SimpleMap(acs_prod.cp_gamma(foo), acs_prod(foo,1:3), 'gamma cp (unitless)')
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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QCR2150A50351 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila.cfg.instruments2run = {'QCR2150A50351'};
%ila.cfg.days2run = datetime(2023,4,3):datetime(2023,11,16);
ila.cfg.days2run = datetime(2023,4,3):datetime(2024,8,22);

% populate ila.instrument
ila.Read('prod');

% extract par data from obj
par_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.a;
par_temp = round_timestamp(par_temp);

% build par table: merge lat, lon
replace_consecutive_nan = 0*60; % 0h
par = merge_timeseries(par_temp, latlon, {'lat', 'lon'}, '', replace_consecutive_nan);

% Remove NaN
if any(any(isnan([par.lat par.lon]),2))
  par(any(isnan([par.lat par.lon]),2), :) = [];
  warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
    sum(any(isnan([par.lat par.lon]),2)), replace_consecutive_nan/60)
end

% add units and precision
par = renamevars(par, {'par_n'}, {'bincount'});
par.Properties.VariableUnits = {'', 'degrees', 'degrees', 'uE.m^{-2}.s^{-1}', 'uE.m^{-2}.s^{-1}', 'none'};
par.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f'};

% sort by date
par = sortrows(par, 'dt');

ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_PAR_timeseries']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_PAR_timeseries']), 'fig')
close figure 40

% par.par = par.par.*10^-6; % change from uE.m^-2.s^-1 to E.m^-2.s^-1
% par.par_sd = par.par_sd.*10^-6; % change from uE.m^-2.s^-1 to E.m^-2.s^-1

SimpleMap(par.par, par(:,1:3), 'PAR [uE.m^{-2}.s^{-1}]')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_PAR_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_PAR_map']), 'fig')
close figure 1

filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, ila.cfg.instruments2run{:}, ...
  datetime(min(par.dt), 'Format', 'yyyyMMdd'), datetime(max(par.dt), 'Format', 'yyyyMMdd'), datetime('today', 'Format', 'yyyyMMdd'));

% export product to SeaBASS format
ila.meta.documents = [cruise '_QCR2150A50351_ProcessingReport_V2.pdf'];
ila.meta.calibration_files = 'QCR2150A50351_CharSheet.pdf';
exportSeaBASS(fullfile(ila.instrument.FLOW.path.prod, filename),...
    ila.meta,...
    par,...
    {'', '', ''});
sprintf('%s_InLine_%s_Product.sb saved', cruise, cell2mat(ila.cfg.instruments2run))

par = renamevars(par, {'bincount'}, {'par_n'});

% save PAR prod
fprintf('Export to mat and csv... ');
filename = strrep(filename, '.sb', '');
save(fullfile(ila.instrument.FLOW.path.prod, filename), 'par');
writetable(par, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
fprintf('Done\n');







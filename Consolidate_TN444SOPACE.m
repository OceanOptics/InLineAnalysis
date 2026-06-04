% Consolidate daily processing, merge variables, and export to SeaBASS
% author: Guillaume Bourdin
% created: Feb 12, 2025

%% Import data
cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master')

official_name = 'PVST_SOPACE_TN444';
cruise = 'TN444SOPACE';
% Load InLineAnalysis and the configuration
ila = InLineAnalysis(['cfg' filesep cruise '_cfg.m']);

path_dev = strrep(ila.instrument.FLOW.path.prod, 'prod', 'DeviceFiles');

% create Graph folder if it doesn't exist
if ~isfolder(fullfile(ila.instrument.FLOW.path.prod, 'plots'))
  mkdir(fullfile(ila.instrument.FLOW.path.prod, 'plots'))
end

%% 'NMEA','FLOW','SBE3845TN444','ACS91','ACS111','HyperBB8005','SUVF6266','LISST100X1183','WS3S1081P','ALFA11'
if ~isfile(fullfile(ila.instrument.FLOW.path.prod, [official_name '_InLine_GPS_20241024_20241216_Product_v20250829.mat']))
  ila.cfg.days2run = datetime(2025,5,5):datetime(2025,6,15);
  ila.cfg.instruments2run = {'NMEA'};
  
  % populate ila.instrument
  ila.Read('prod');
  
  % extract data from ila
  latlon = ila.instrument.('NMEA').prod.a;
  % latlon = removevars(latlon, {'dt_instrument','dt_instrument_avg_n','lat_avg_sd','lat_avg_n','lon_avg_sd','lon_avg_n',...
  %   'mag_track','mag_track_avg_sd','mag_track_avg_n','altitude','altitude_avg_sd','altitude_avg_n','gps_qual','gps_qual_avg_sd','gps_qual_avg_n',...
  %   'num_sats','num_sats_avg_sd','num_sats_avg_n','horizontal_dil','horizontal_dil_avg_sd','horizontal_dil_avg_n',...
  %   'true_course','true_course_avg_sd','true_course_avg_n'}); % CSV import
  var2rm = latlon.Properties.VariableNames(contains(latlon.Properties.VariableNames, {'gps_dt','_avg_n','_avg_sd','wind','atm','spd','heading','depth','water','air','magnetic'}));
  latlon = removevars(latlon, var2rm); % RAW import
  latlon = round_timestamp(latlon);

  SimpleMap(NaN(size(latlon.dt)), latlon(:,1:3), 'TN444-SOPACE track')
  close Figure 1

  filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s.mat', official_name, 'GPS', ...
    datetime(min(latlon.dt), 'Format', 'yyyyMMdd'), datetime(max(latlon.dt), 'Format', 'yyyyMMdd'), ...
    datetime('today', 'Format', 'yyyyMMdd'));
  save(fullfile(ila.instrument.FLOW.path.prod, filename), 'latlon')
  writetable(latlon, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
else
  load(fullfile(ila.instrument.FLOW.path.prod, [official_name '_InLine_GPS_20250505_20250615_Product_v20250829.mat']))
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TSG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfile(fullfile(ila.instrument.FLOW.path.prod, [official_name '_InLine_TSG_20250505_20250615_Product_v20250829.mat']))
  ila.cfg.days2run = datetime(2025,5,5):datetime(2025,6,15);
  ila.cfg.instruments2run = {'SBE3845TN444'};
  
  % populate ila.instrument
  ila.Read('prod');
  
  % extract TSG data from ila
  tsg_temp = round_timestamp(ila.instrument.('SBE3845TN444').prod.a);
  
  % build TSG table: merge lat, lon
  replace_consecutive_nan = 3*60; % 3h
  tsg = merge_timeseries(tsg_temp, latlon, {'lat', 'lon'}, '', replace_consecutive_nan);
  % Remove NaN
  if any(any(isnan([tsg.lat tsg.lon]),2))
    warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
      sum(any(isnan([tsg.lat tsg.lon]),2)), replace_consecutive_nan/60)
    tsg(any(isnan([tsg.lat tsg.lon]),2), :) = [];
  end
  tsg = removevars(tsg, {'dt_instrument','dt_instrument_avg_n','cond_avg_n','sss_avg_n','sst_avg_n'});
  tsg = renamevars(tsg,{'sv_avg_n'},{'avg_n'});
  tsg.Properties.VariableUnits = {'','degrees','degrees','degreesC','degreesC','S/m','S/m','PSU','PSU','m/sec','m/sec','none'};
  tsg.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.2f','%.2f','%i'};
  
  % sort by date
  tsg = sortrows(tsg, 'dt');
  
  ila.visProd_timeseries()
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_TSG_timeseries']), 'jpg')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_TSG_timeseries']), 'fig')
  close figure 30
  
  SimpleMap(tsg.sst, tsg(:,1:3), 'TSG SST [°C]')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_TSG_SST_map']), 'jpg')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_TSG_SST_map']), 'fig')
  close figure 1
  SimpleMap(tsg.sss, tsg(:,1:3), 'TSG SSS [PSU]')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_TSG_SSS_map']), 'jpg')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_TSG_SSS_map']), 'fig')
  close figure 1
  
  filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', official_name, 'TSG', ...
    datetime(min(tsg.dt), 'Format', 'yyyyMMdd'), datetime(max(tsg.dt), 'Format', 'yyyyMMdd'), ...
    datetime('today', 'Format', 'yyyyMMdd'));
  
  % % export product to SeaBASS format
  % ila.meta.documents = [official_name '_TSG_ProcessingReport_V2.pdf';
  % ila.meta.calibration_files = cell2mat(list_dev(i));
  % exportSeaBASS(fullfile(ila.instrument.FLOW.path.prod, filename),...
  %     ila.meta,...
  %     tsg,...
  %     {'', '', ''});
  % sprintf('%s_InLine_%s_Product.sb saved', official_name, cell2mat(ila.cfg.instruments2run))
  
  % save TSG prod
  fprintf('Export to mat and csv... ');
  filename = strrep(filename, '.sb', '');
  save(fullfile(ila.instrument.FLOW.path.prod, filename), 'tsg');
  writetable(tsg, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
  fprintf('Done\n');
else
  load(fullfile(ila.instrument.FLOW.path.prod, [official_name '_InLine_TSG_20250505_20250615_Product_v20250829.mat']))
end

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
% if any(any(isnan([wscd.lat wscd.lon]),2))
%   warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
%     sum(any(isnan([wscd.lat wscd.lon]),2)), replace_consecutive_nan/60)
%   wscd(any(isnan([wscd.lat wscd.lon]),2), :) = [];
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
% saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_FDOM_timeseries']), 'jpg')
% close figure 94
% 
% SimpleMap(wscd.fdom, wscd(:,1:3), 'WSCD fdom [v uncalibrated]')
% saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_WSCD_fdom_map']), 'jpg')
% close figure 1
% 
% filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', official_name, ila.cfg.instruments2run{:}, ...
%   datetime(min(wscd.dt), 'Format', 'yyyyMMdd'), datetime(max(wscd.dt), 'Format', 'yyyyMMdd'), ...
%   datetime('today', 'Format', 'yyyyMMdd'));
% 
% % % export product to SeaBASS format
% % ila.meta.documents = [official_name '_WSCD_ProcessingReport_V2.pdf';
% % ila.meta.calibration_files = cell2mat(list_dev(i));
% % exportSeaBASS(fullfile(ila.instrument.FLOW.path.prod, filename),...
% %     ila.meta,...
% %     wscd,...
% %     {'', '', ''});
% % sprintf('%s_InLine_%s_Product.sb saved', official_name, cell2mat(ila.cfg.instruments2run))
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
%   'Graphs' filesep official_name '_PAR_timeseries'], 'jpg')
% close figure 1
%
% filename = [ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod official_name ...
%   '_InLine_' cell2mat(ila.cfg.instruments2run) '_Product_v' datestr(now, 'yyyymmdd') '.sb'];
%
% % export product to SeaBASS format
% ila.meta.documents = [official_name '_PAR_ProcessingReport_V2.pdf'];
% ila.meta.calibration_files = 'PAR-50168_CalSheet.pdf';
% exportSeaBASS(filename,...
%     ila.meta,...
%     par,...
%     {'', '', ''});
% sprintf('%s_InLine_%s_Product.sb saved', official_name, cell2mat(ila.cfg.instruments2run))
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
%     official_name '_InLine_PAR_prod'], 'par');
% writetable(par, [ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
%     official_name '_InLine_PAR_prod.csv']);
% fprintf('Done\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUVF6266 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila.cfg.instruments2run = {'SUVF6266'};
if ~isfile(fullfile(ila.instrument.FLOW.path.prod, [official_name '_InLine_' ila.cfg.instruments2run{:} '_20250505_20250615_Product_v20250829.mat']))
  ila.cfg.days2run = datetime(2025,5,5):datetime(2025,6,15);
  
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
  if any(any(isnan([suvf.lat suvf.lon]),2))
    warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
      sum(any(isnan([suvf.lat suvf.lon]),2)), replace_consecutive_nan/60)
    suvf(any(isnan([suvf.lat suvf.lon]),2), :) = [];
  end
  
  % add units and precision
  suvf = renamevars(suvf, {'sst','sss','fdom_n'}, {'t','s','bincount'});
  suvf.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'v_uncalibrated', 'v_uncalibrated', 'none'};
  suvf.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f'};
  
  % sort by date
  suvf = sortrows(suvf, 'dt');
  
  ila.visProd_timeseries()
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_SUVF_fdom_timeseries']), 'jpg')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_SUVF_fdom_timeseries']), 'fig')
  close figure 50
  
  SimpleMap(suvf.fdom, suvf(:,1:3), 'SUVF fdom [v uncalibrated]')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_SUVF_fdom_map']), 'jpg')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_SUVF_fdom_map']), 'fig')
  close figure 1
  
  filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', official_name, ila.cfg.instruments2run{:}, ...
    datetime(min(suvf.dt), 'Format', 'yyyyMMdd'), datetime(max(suvf.dt), 'Format', 'yyyyMMdd'), datetime('today', 'Format', 'yyyyMMdd'));
  
  % export product to SeaBASS format
  ila.meta.documents = [official_name '_SUVF_ProcessingReport_V2.pdf'];
  ila.meta.calibration_files = [ila.cfg.instruments2run{:} '_CharSheet.pdf'];
  exportSeaBASS(fullfile(ila.instrument.FLOW.path.prod, filename),...
      ila.meta,...
      suvf,...
      {'', '', ''});
  sprintf('%s_InLine_%s_Product.sb saved', official_name, cell2mat(ila.cfg.instruments2run))
  
  suvf = renamevars(suvf, {'t','s','bincount'}, {'sst','sss','fdom_n'});
  
  % save SUVF prod
  fprintf('Export to mat and csv... ');
  filename = strrep(filename, '.sb', '');
  save(fullfile(ila.instrument.FLOW.path.prod, filename), 'suvf');
  writetable(suvf, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
  fprintf('Done\n');
else
  load(fullfile(ila.instrument.FLOW.path.prod, [official_name '_InLine_' ila.cfg.instruments2run{:} '_20250505_20250615_Product_v20250829.mat']))
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HBB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_leg = table();
list_leg.dt = [...
  datetime(2025,5,5) datetime(2025,6,15);
  ];
list_leg.experiment = {'PVST_SOPACE'};
list_leg.cruise = {'TN444'};

ila.cfg.instruments2run = {'HyperBB8005'};
hbb = [];
for i=1:size(list_leg, 1)
  ila.cfg.days2run = list_leg.dt(i, 1):list_leg.dt(i, 2);
  
  % populate ila.instrument
  ila.Read('prod');
  
  % update metadata for each leg
  ila.meta.experiment = list_leg.experiment{i};
  ila.meta.cruise = list_leg.cruise{i};

  % % remove flagged data
  % flag = read_flagbit(ila.instrument.(list_leg.instru{:}).prod.p.flag_bit, 'BB');

  % extract HyperBB data from obj
  hbb_lambda = ila.instrument.(ila.cfg.instruments2run{:}).lambda;
  hbb_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.p;
  hbb_temp = round_timestamp(hbb_temp);
  
  % build hbb table: merge lat, lon, sst, sss
  replace_consecutive_nan = 3*60/5; % 3h (5min binning)
  hbb_temp = merge_timeseries(hbb_temp, suvf, {'fdom'});
  hbb_temp = merge_timeseries(hbb_temp, tsg, {'lat', 'lon', 'sst', 'sss'});
  hbb_temp = merge_timeseries(hbb_temp, latlon, {'lat', 'lon'}, '', replace_consecutive_nan);
  
  % Remove NaN
  if any(any(isnan([hbb_temp.lat hbb_temp.lon]),2))
    warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
      sum(any(isnan([hbb_temp.lat hbb_temp.lon]),2)), replace_consecutive_nan/60*5); % (5min binning)
    hbb_temp(any(isnan([hbb_temp.lat hbb_temp.lon]),2), :) = [];
  end

  % add units and precision
  hbb_temp.Properties.VariableNames = {'dt','lat','lon','sst','sss','fdom','VSF124','bbp','VSF124_sd','bincount','gamma_bbp','poc','cphyto_bbp','flag_bit'};
  hbb_temp.Properties.VariableUnits = {'','degrees','degrees','degreesC','PSU','v_uncalibrated','1/m/sr','1/m','1/m/sr','none','unitless','ug/L','ug/L','unitless'};
  hbb_temp.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%i','%.4f','%.2f','%.2f','%i'};

  % add units and precision
  hbb4seabass = hbb_temp(:, ~contains(hbb_temp.Properties.VariableNames, {'fdom', 'gamma_bbp', 'poc', 'cphyto','flag_bit'}));
  hbb4seabass = renamevars(hbb4seabass, {'sst', 'sss', 'VSF124','VSF124_sd','bincount'}, {'t','s','VSF_124ang','VSF_124ang_sd','bincount'});
  
  % sort by date
  hbb4seabass = sortrows(hbb4seabass, 'dt'); 
  hbb4seabass.bbp(hbb4seabass.bbp<0)=NaN;
  hbb4seabass.VSF_124ang_sd(hbb4seabass.VSF_124ang<0)=NaN;
  hbb4seabass.VSF_124ang(hbb4seabass.VSF_124ang<0)=NaN;
  hbb4seabass(all(isnan(hbb4seabass.bbp),2),:)=[];
  
  %%% HBB 3D plots %%%
  save_figures = false;
  ila.SpectralQC('BB', {'prod'}, save_figures); % AC or BB
  close all
  
  % plot time series
  ila.visProd_timeseries()
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_HBB_gammabbp_bbp550']), 'jpg')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_HBB_gammabbp_bbp550']), 'fig')
  close figure 26
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_HBB_BBparticulate_timeseries']), 'jpg')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_HBB_BBparticulate_timeseries']), 'fig')
  close figure 24
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_HBB_POCparticulate_timeseries']), 'jpg')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_HBB_POCparticulate_timeseries']), 'fig')
  close figure 23
  % saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_HBB_BBdissolved_timeseries']), 'jpg')
  % saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_HBB_BBdissolved_timeseries']), 'fig')
  % close figure 25
  
  SimpleMap(hbb4seabass.bbp(:,hbb_lambda == 530), hbb4seabass(:,1:3), 'bbp (530 nm) [m^-^1]')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_hbb_bbp_map']), 'jpg')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_hbb_bbp_map']), 'fig')
  close figure 1
  
  filename_part = sprintf('%s_%s_InLine_%s_%s_%s_Particulate_v%s.sb', list_leg.experiment{i}, list_leg.cruise{i}, ila.cfg.instruments2run{:}, ...
    datetime(min(hbb4seabass.dt), 'Format', 'yyyyMMdd'), datetime(max(hbb4seabass.dt), 'Format', 'yyyyMMdd'), ...
    datetime('today', 'Format', 'yyyyMMdd'));
  
  % export product to SeaBASS format
  ila.meta.documents = sprintf('%s_HBB_ProcessingReport_v%s.pdf', official_name, datetime(2025,2,20, 'Format', 'yyyyMMdd'));
  ila.meta.calibration_files = 'Hbb_Cal_Plaque_20241010_142447.mat,Hbb_Cal_Temp_20241008_130132.mat';
  exportSeaBASS(fullfile(ila.instrument.FLOW.path.prod, filename_part),...
      ila.meta,...
      hbb4seabass,...
      {string(hbb_lambda)',string(hbb_lambda)',string(hbb_lambda)',''});
  fprintf('%s saved\n', filename_part)
  
  hbb = [hbb; hbb_temp];
end

% sort by date
hbb = sortrows(hbb, 'dt');
ila.cfg.days2run = list_leg.dt(1, 1):list_leg.dt(end, 2);
ila.instrument.(ila.cfg.instruments2run{:}).prod.p = hbb;
ila.instrument.(ila.cfg.instruments2run{:}).prod.p = renamevars(ila.instrument.(ila.cfg.instruments2run{:}).prod.p, {'VSF124', 'VSF124_sd'}, {'betap','betap_sd'});

% plot time series
ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_HBB_gammabbp_bbp550']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_HBB_gammabbp_bbp550']), 'fig')
close figure 26
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_HBB_BBparticulate_timeseries']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_HBB_BBparticulate_timeseries']), 'fig')
close figure 24
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_HBB_POCparticulate_timeseries']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_HBB_POCparticulate_timeseries']), 'fig')
close figure 23
% saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_HBB_BBdissolved_timeseries']), 'jpg')
% saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_HBB_BBdissolved_timeseries']), 'fig')
% close figure 25

SimpleMap(hbb.bbp(:,hbb_lambda == 530), hbb(:,1:3), 'bbp (530 nm) [m^-^1]')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_hbb_bbp_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_hbb_bbp_map']), 'fig')
close figure 1

filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', official_name, ila.cfg.instruments2run{:}, ...
  datetime(min(hbb.dt), 'Format', 'yyyyMMdd'), datetime(max(hbb.dt), 'Format', 'yyyyMMdd'), ...
  datetime('today', 'Format', 'yyyyMMdd'));
  
% save HBB prod
fprintf('Export to mat and csv... ');
filename = strrep(filename, '.sb', '');
save(fullfile(ila.instrument.FLOW.path.prod, filename), 'hbb', 'hbb_lambda');
writetable(hbb, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
fprintf('Done\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_leg = table();
list_leg.dt = [...
  datetime(2025,5,5) datetime(2025,6,15);
  ];

list_leg.experiment = {'PVST_SOPACE'};
list_leg.cruise = {'TN444'};
list_leg.devfile = {'ACS111_20231006.dev'};
list_leg.instru = {'ACS111'};

data_AC = struct('particulate', [], 'product', []);
acs_prod = [];
acs_prod_diw = [];

for i=1:size(list_leg, 1)
  % cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master')
  ila = InLineAnalysis(['cfg' filesep cruise '_cfg.m']);
  ila.cfg.instruments2run = list_leg.instru(i);
  ila.cfg.days2run = list_leg.dt(i, 1):list_leg.dt(i, 2);

  % populate ila.instrument
  ila.Read('prod');
  % delete dissolved data => rubbish in this case
  ila.instrument.ACS111.prod.g = table();

  % update metadata for each leg
  ila.meta.experiment = list_leg.experiment{i};
  ila.meta.cruise = list_leg.cruise{i};
  
  % get wavelength
  [ila.instrument.(list_leg.instru{i}).lambda_c, ...
    ila.instrument.(list_leg.instru{i}).lambda_a] = importACSDeviceFile(fullfile(path_dev, list_leg.devfile{i}));

  % set section name
  ref = sprintf('%s_%s_%s', list_leg.instru{i}, datetime(ila.cfg.days2run(1), 'Format', 'yyyyMMdd'), ...
    datetime(ila.cfg.days2run(end), 'Format', 'yyyyMMdd'));
  
  % remove flagged data
  flag = read_flagbit(ila.instrument.(list_leg.instru{i}).prod.p.flag_bit, 'AC');
  % flag_info = InlineFlagInfo('ACS');

  % flag.HH_G50_flag(flag.HH_G50_flag & ila.instrument.(list_leg.instru{i}).prod.p.HH_G50 > 0 & ila.instrument.(list_leg.instru{i}).prod.p.HH_G50 <= 500) = false;
  % ila.instrument.(list_leg.instru{i}).prod.p.flag_bit = set_flagbit(flag);
  
  % remove flagged products
  ila.instrument.(list_leg.instru{i}).prod.p.poc(flag.poc_flag) = NaN;
  ila.instrument.(list_leg.instru{i}).prod.p.chl_ap676lh(flag.chl_ap676lh_flag) = NaN;
  ila.instrument.(list_leg.instru{i}).prod.p.gamma(flag.cp_step) = NaN;
  ila.instrument.(list_leg.instru{i}).prod.p.poc(flag.cp_step) = NaN;
  ila.instrument.(list_leg.instru{i}).prod.p.chl_Halh(flag.chl_Halh_flag) = NaN;
  ila.instrument.(list_leg.instru{i}).prod.p.chl_Halh(flag.cp_step) = NaN;
  ila.instrument.(list_leg.instru{i}).prod.p.HH_mphi(flag.cp_step) = NaN;
  ila.instrument.(list_leg.instru{i}).prod.p.HH_G50(flag.HH_G50_flag) = NaN;
  ila.instrument.(list_leg.instru{i}).prod.p.HH_G50(flag.cp_step) = NaN;
  % ila.instrument.(list_leg.instru{i}).prod.p.chl_Halh(flag.chlratio_flag) = NaN;
  
  % % remove suspicious products
  % ila.instrument.(list_leg.instru{i}).prod.p.gamma(flag.gamma_suspicious) = NaN;
  % % ila.instrument.(list_leg.instru{i}).prod.p.gamma(acs_prod.gamma > 2) = NaN;
  % ila.instrument.(list_leg.instru{i}).prod.p.poc(flag.ap_step) = NaN;
  % ila.instrument.(list_leg.instru{i}).prod.p.poc(flag.poc_suspicious) = NaN;
  % ila.instrument.(list_leg.instru{i}).prod.p.chl_ap676lh(flag.chl_ap676lh_suspicious) = NaN;
  % ila.instrument.(list_leg.instru{i}).prod.p.chl_Halh(flag.chl_Halh_suspicious) = NaN;
  % ila.instrument.(list_leg.instru{i}).prod.p.HH_mphi(flag.HH_G50_mphi_suspicious) = NaN;
  % ila.instrument.(list_leg.instru{i}).prod.p.HH_G50(flag.HH_G50_mphi_suspicious) = NaN;
  
  if isfield(ila.instrument.(list_leg.instru{i}).prod, 'g')
    if ~isempty(ila.instrument.(list_leg.instru{i}).prod.g)
      % keep only ag/cg for good lambda: between 410 and 580 nm
      sel_a = ila.instrument.(list_leg.instru{i}).lambda_a >= 410 & ila.instrument.(list_leg.instru{i}).lambda_a < 580;
      sel_c = ila.instrument.(list_leg.instru{i}).lambda_c >= 410 & ila.instrument.(list_leg.instru{i}).lambda_c < 580;
      % FIX ISSUE WITH AG/CG and FDOM in the SECOND PART OF THE LEG
      % recompute exponential fit on ag
      [ila.instrument.(list_leg.instru{i}).prod.g.y_intercp_fit_ag, ...
        ila.instrument.(list_leg.instru{i}).prod.g.base_fit_ag, ~, ~, ...
        ila.instrument.(list_leg.instru{i}).prod.g.RMSE_fit_ag] = FitExp(ila.instrument.(list_leg.instru{i}).lambda_a(sel_a), ...
        ila.instrument.(list_leg.instru{i}).prod.g.ag(:, sel_a), ila.instrument.(list_leg.instru{i}).prod.g.ag_sd(:, sel_a));
      % add fit flag
      ila.instrument.(list_leg.instru{i}).prod.g.ag_fitflag = false(size(ila.instrument.(list_leg.instru{i}).prod.g, 1), 1);
      ila.instrument.ACS91.prod.g.ag_fitflag(ila.instrument.(list_leg.instru{i}).prod.g.RMSE_fit_ag > 0.0025) = true;
      % recompute exponential fit on cg
      [ila.instrument.(list_leg.instru{i}).prod.g.y_intercp_fit_cg, ila.instrument.(list_leg.instru{i}).prod.g.base_fit_cg, ~, ~, ...
        ila.instrument.(list_leg.instru{i}).prod.g.RMSE_fit_cg] = FitExp(ila.instrument.(list_leg.instru{i}).lambda_c(sel_c), ...
        ila.instrument.(list_leg.instru{i}).prod.g.cg(:, sel_c), ila.instrument.(list_leg.instru{i}).prod.g.cg_sd(:, sel_c));
      % add fit flag
      ila.instrument.(list_leg.instru{i}).prod.g.cg_fitflag = false(size(ila.instrument.(list_leg.instru{i}).prod.g, 1), 1);
      ila.instrument.(list_leg.instru{i}).prod.g.cg_fitflag(ila.instrument.(list_leg.instru{i}).prod.g.RMSE_fit_cg > 0.0025) = true;
      % re-merge cdom
      replace_consecutive_nan = 3*60; % 3h
      ila.instrument.(list_leg.instru{i}).prod.g.fdom = NaN(size(ila.instrument.(list_leg.instru{i}).prod.g.dt));
      ila.instrument.(list_leg.instru{i}).prod.g = merge_timeseries(ila.instrument.(list_leg.instru{i}).prod.g, suvf, {'fdom'}, '', replace_consecutive_nan);
      
      ila.Write('prod', 'diw')
    end
  end

  ila.visProd_timeseries()
  saveGraph(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_' ref '_ACS_prod_regressions']), 'jpg')
  saveGraph(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_' ref '_ACS_prod_regressions']), 'fig')
  close figure 14
  saveGraph(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_' ref '_ACS_prod_part_timeseries']), 'jpg')
  saveGraph(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_' ref '_ACS_prod_part_timeseries']), 'fig')
  close figure 10
  if isfield(ila.instrument.(list_leg.instru{i}).prod, 'g')
    if ~isempty(ila.instrument.(list_leg.instru{i}).prod.g)
      saveGraph(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_' ref '_ACS_prod_diw_timeseries']), 'jpg')
      saveGraph(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_' ref '_ACS_prod_diw_timeseries']), 'fig')
      close figure 16
    end
  end

  %%% AC 3D plots %%%
  save_figures = true;
  ila.SpectralQC('AC', {'prod'}, save_figures); % AC or BB
  close all

  % extract ACp data from obj
  ACp = ila.instrument.(list_leg.instru{i}).prod.p;
  ACp = sortrows(ACp, 'dt');
  % build AC table: merge lat, lon, sst, sss, and fdom
  replace_consecutive_nan = 3*60; % 3h
  ACp = merge_timeseries(ACp, suvf, {'fdom'}, '', replace_consecutive_nan);
  ACp = merge_timeseries(ACp, tsg, {'lat', 'lon', 'sst', 'sss'});
  ACp = merge_timeseries(ACp, latlon, {'lat', 'lon'}, '', replace_consecutive_nan);
  % rename variables for SeaBASS
  ACp = renamevars(ACp, {'sst','sss'}, {'t','s'});
  % Remove NaN
  if any(any(isnan([ACp.lat ACp.lon]),2))
    warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
      sum(any(isnan([ACp.lat ACp.lon]),2)), replace_consecutive_nan/60)
    ACp(any(isnan([ACp.lat ACp.lon]),2), :) = [];
  end
  % add variable units and description
  ACp.Properties.VariableUnits = {'','degrees','degrees','degreesC','PSU','ppb','1/m','1/m','1/m','1/m', ...
    'none','none','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m', ...
     'ug/L','unitless','ug/L','unitless','unitless','ug/L','microns','unitless','unitless'};
  ACp.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f', ...
    '%i','%i','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f', ...
    '%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%i'};

  % split into particulate table
  id_ancil = logical(sum(categorical(ACp.Properties.VariableNames) == {'dt','lat','lon','t','s','fdom','ap_n'}'));
  id_particulate = logical(sum(categorical(ACp.Properties.VariableNames) == {'ap','ap_se','cp','cp_se'}'));
  data_AC.particulate.(ref) = ACp(:, id_ancil | id_particulate);
  % rename variables for SeaBASS
  data_AC.particulate.(ref) = renamevars(data_AC.particulate.(ref), {'fdom','ap_n'}, {'cdomf','bincount'});

  % export particulate to SeaBASS format
  filename_part = sprintf('%s_%s_InLine_%s_Particulate_v%s.sb', list_leg.experiment{i}, list_leg.cruise{i}, ref, ...
    datetime('today', 'Format', 'yyyyMMdd'));
  ila.meta.documents = sprintf('%s_ACS_ProcessingReport_v%s.pdf', official_name, datetime(2025,2,20, 'Format', 'yyyyMMdd'));
  ila.meta.calibration_files = list_leg.devfile{i};
  exportSeaBASS(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, filename_part),...
      ila.meta,...
      data_AC.particulate.(ref),...
      {'', string(ila.instrument.(list_leg.instru{i}).lambda_a),...
      string(ila.instrument.(list_leg.instru{i}).lambda_c),...
      string(ila.instrument.(list_leg.instru{i}).lambda_a),...
      string(ila.instrument.(list_leg.instru{i}).lambda_c),''});
  fprintf('%s saved\n', filename_part)

  % save in .mat file
  lambda.a = ila.instrument.(list_leg.instru{i}).lambda_a;
  lambda.c = ila.instrument.(list_leg.instru{i}).lambda_c;
  acs_part = data_AC.particulate.(ref);
  % rename variables for .mat
  acs_part = renamevars(acs_part, {'t','s','cdomf'}, {'sst','sss','fdom'});
  save(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, strrep(filename_part, '.sb', '')), 'acs_part', 'lambda');
  
  % split into product table
  data_AC.product.(ref) = ACp(:, id_ancil| ~id_particulate);
  % rename variables for SeaBASS
  data_AC.product.(ref) = renamevars(data_AC.product.(ref), {'poc','gamma','chl_ap676lh', 'fdom'}, ...
    {'POC_cp','cp_gamma','Chl_lineheight', 'cdomf'});

  % export products to SeaBASS format
  filename_prod = strrep(filename_part, '_Particulate_', '_ProductsFull_');
  ila.meta.documents = sprintf('%s_ACS_ProcessingReport_v%s.pdf', official_name, datetime(2025,2,20, 'Format', 'yyyyMMdd'));
  exportSeaBASS(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, filename_prod), ila.meta, data_AC.product.(ref));
  fprintf('%s saved\n', filename_prod)
  
  % ACS merged prod
  acs_prod = [acs_prod; data_AC.product.(ref)];

  % keep only old variables
  acs = data_AC.product.(ref)(:, contains(data_AC.product.(ref).Properties.VariableNames, {'dt', 'lat', 'lon', ...
    't','s','cdomf','Chl_lineheight','POC_cp','cp_gamma','ap_n','cp_n','flag_bit'}) & ...
    ~contains(data_AC.product.(ref).Properties.VariableNames, {'agaus'}));
  
  % export product to SeaBASS format
  filename_prod = strrep(filename_prod, '_ProductsFull_', '_ProductsLegacy_');
  exportSeaBASS(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, filename_prod), ila.meta, acs);
  fprintf('%s saved\n', filename_prod)

  if isfield(ila.instrument.(list_leg.instru{i}).prod, 'g')
    if ~isempty(ila.instrument.(list_leg.instru{i}).prod.g)

      % % rebuild ag from exponential fit parameters
      % ag_rebuilt = NaN(size(ACg.ag));
      % cg_rebuilt = NaN(size(ACg.cg));
      % % define exponential function
      % expfun = @(p, xd) p(1) * exp(p(2) * (xd - 440));
      % for t = 1:size(ACg, 1)
      %   ag_rebuilt(t, :) = expfun([ACg.y_intercp_fit_ag(t) ACg.base_fit_ag(t)], ila.instrument.(list_leg.instru{i}).lambda_a);
      %   cg_rebuilt(t, :) = expfun([ACg.y_intercp_fit_cg(t) ACg.base_fit_cg(t)], ila.instrument.(list_leg.instru{i}).lambda_c);
      % end
      % % merge ag and rebuilt ag
      % ila.instrument.(list_leg.instru{i}).prod.g.ag(isnan(ila.instrument.(list_leg.instru{i}).prod.g.ag)) = ag_rebuilt(isnan(ila.instrument.(list_leg.instru{i}).prod.g.ag));
      % % merge cg and rebuilt cg
      % ila.instrument.(list_leg.instru{i}).prod.g.cg(isnan(ila.instrument.(list_leg.instru{i}).prod.g.cg)) = cg_rebuilt(isnan(ila.instrument.(list_leg.instru{i}).prod.g.cg));

      % extract ACg data from obj
      ACg = ila.instrument.(list_leg.instru{i}).prod.g;
      ACg = sortrows(ACg, 'dt');
      % build AC table: merge lat, lon, sst, sss, and fdom
      replace_consecutive_nan = 3*60; % 3h
      ACg = merge_timeseries(ACg, tsg, {'lat', 'lon', 'sst', 'sss'});
      ACg = merge_timeseries(ACg, latlon, {'lat', 'lon'}, '', replace_consecutive_nan);
      % rename variables for SeaBASS
      ACg = renamevars(ACg, {'sst','sss'}, {'t','s'});
      % Remove NaN
      if any(any(isnan([ACg.lat ACg.lon]),2))
        warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
          sum(any(isnan([ACg.lat ACg.lon]),2)), replace_consecutive_nan/60)
        ACg(any(isnan([ACg.lat ACg.lon]),2), :) = [];
      end
      % keep only good lambda: between 410 and 580 nm
      ACg.ag = ACg.ag(:, sel_a);
      ACg.ag_sd = ACg.ag_sd(:, sel_c);
      ACg.cg = ACg.cg(:, sel_a);
      ACg.cg_sd = ACg.cg_sd(:, sel_c);
      % move fdom variable after s
      ACg = movevars(ACg, 'fdom', 'After', 's');
      % add variable units and description
      ACg.Properties.VariableUnits = {'','degrees','degrees','degreesC','PSU','ppb','1/m','1/m','1/m','1/m','unitless','unitless',...
        'unitless','unitless','unitless','unitless','1/m','unitless','unitless','none','1/m','unitless','unitless'};
      ACg.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%i','%i',...
        '%.4f','%.4f','%.4f','%.4f','%.4f','%i','%.4f','%.4f','%.4f','%i','%i'};
    
      % ACS merged prod DIW
      id_diw_prod = logical(all(categorical(ACg.Properties.VariableNames) ~= {'ag','ag_sd','cg','cg_sd'}'));
      acs_prod_diw = [acs_prod_diw; ACg(:, id_diw_prod)];
      % split into dissolved table
      id_dissolved = logical(sum(categorical(ACg.Properties.VariableNames) == {'dt','lat','lon','t','s','fdom','ag','ag_sd','cg','cg_sd','ag_n'}'));
      data_AC.dissolved.(ref) = ACg(:, id_dissolved);
      % rename variables for SeaBASS
      data_AC.dissolved.(ref) = renamevars(data_AC.dissolved.(ref), {'fdom','ag_n'}, {'cdomf','bincount'});
    
      % export dissolved to SeaBASS format
      filename_diw = strrep(filename_part, '_Particulate_', '_Dissolved_');
      ila.meta.documents = sprintf('%s_ACS_ProcessingReport_v%s.pdf', official_name, datetime(2025,2,20, 'Format', 'yyyyMMdd'));
      ila.meta.calibration_files = list_leg.devfile{i};
      exportSeaBASS(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, filename_diw),...
          ila.meta,...
          data_AC.dissolved.(ref),...
          {'', string(ila.instrument.(list_leg.instru{i}).lambda_a(sel_a)),...
          string(ila.instrument.(list_leg.instru{i}).lambda_c(sel_c)),...
          string(ila.instrument.(list_leg.instru{i}).lambda_a(sel_a)),...
          string(ila.instrument.(list_leg.instru{i}).lambda_c(sel_c)),''});
      fprintf('%s saved\n', filename_diw)

      % save in .mat file
      lambda_diw.a = ila.instrument.(list_leg.instru{i}).lambda_a(sel_a);
      lambda_diw.c = ila.instrument.(list_leg.instru{i}).lambda_c(sel_c);
      acs_diw = data_AC.dissolved.(ref);

      % rename variables for .mat
      acs_diw = renamevars(acs_diw, {'t','s','cdomf'}, {'sst','sss','fdom'});
      save(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, strrep(filename_diw, '.sb', '')), 'acs_diw', 'lambda_diw');
    end
  end
end

% sort by date
acs_prod = sortrows(acs_prod, 'dt');
if isfield(ila.instrument.(list_leg.instru{i}).prod, 'g')
  if ~isempty(ila.instrument.(list_leg.instru{i}).prod.g)
    acs_prod_diw = sortrows(acs_prod_diw, 'dt');
  end
end

if any(strcmp(acs_prod.Properties.VariableNames, 'chl_ap676lh'))
  acs_prod = renamevars(acs_prod, 'chl_ap676lh', 'Chl_lineheight');
end
if any(strcmp(acs.Properties.VariableNames, 'poc'))
  acs_prod = renamevars(acs_prod, 'poc', 'POC_cp');
end
if any(strcmp(acs.Properties.VariableNames, 'gamma'))
  acs_prod = renamevars(acs_prod, 'gamma', 'cp_gamma');
end

ila.cfg.days2run = list_leg.dt(1, 1):list_leg.dt(end, 2);
ila.instrument.(list_leg.instru{i}).prod.p = acs_prod;
if isfield(ila.instrument.(list_leg.instru{i}).prod, 'g')
  if ~isempty(ila.instrument.(list_leg.instru{i}).prod.g)
    ila.instrument.(list_leg.instru{i}).prod.g = acs_prod_diw;
  end
end

% figure()
% subplot(2,3,1); histogram(acs_prod.POC_cp); xlabel('[POC] cp (mg.m^{-3})')
% subplot(2,3,2); histogram(acs_prod.Chl_lineheight); xlabel('a_{p676} line height [chl a] (mg.m^{-3})')
% subplot(2,3,3); histogram(acs_prod.cp_gamma); xlabel('gamma cp (mg.m^{-3})')
% subplot(2,3,4); histogram(acs_prod.chl_Halh); xlabel('Houskeeper [chl] (mg.m^{-3})')
% subplot(2,3,5); histogram(acs_prod.HH_mphi); xlabel('H&H phytoplankton slope size distribution')
% subplot(2,3,6); histogram(acs_prod.HH_G50); xlabel('H&H phytoplankton G50: cross-sectional area (\mum)')

ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, 'plots', [official_name '_ACS_prod_regressions']), 'jpg')
saveGraph(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, 'plots', [official_name '_ACS_prod_regressions']), 'fig')
close figure 14
saveGraph(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, 'plots', [official_name '_ACS_prod_timeseries']), 'jpg')
saveGraph(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, 'plots', [official_name '_ACS_prod_timeseries']), 'fig')
close figure 10
if isfield(ila.instrument.(list_leg.instru{i}).prod, 'g')
  if ~isempty(ila.instrument.(list_leg.instru{i}).prod.g)
    saveGraph(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, 'plots', [official_name '_' ref '_ACS_prod_diw_timeseries']), 'jpg')
    saveGraph(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, 'plots', [official_name '_' ref '_ACS_prod_diw_timeseries']), 'fig')
    close figure 16
  end
end

SimpleMap(acs_prod.chl_Halh, acs_prod(:,1:3), 'Houskeeper [chl] (mg.m^{-3})')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_ACS_chl_Houskeeper_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_ACS_chl_Houskeeper_map']), 'fig')
close figure 1

SimpleMap(acs_prod.HH_G50, acs_prod(:,1:3), 'H&H phytoplankton G50: cross-sectional area (\mum)')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_ACS_H&H_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_ACS_H&H_map']), 'fig')
close figure 1

SimpleMap(acs_prod.POC_cp, acs_prod(:,1:3), '[POC] cp (mg.m^{-3})')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_ACS_POC_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_ACS_POC_map']), 'fig')
close figure 1

SimpleMap(acs_prod.Chl_lineheight, acs_prod(:,1:3), 'a_{p676} line height [chl a] (mg.m^{-3})')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_ACS_chl_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_ACS_chl_map']), 'fig')
close figure 1

foo = acs_prod.cp_gamma > prctile(acs_prod.cp_gamma, 1) & acs_prod.cp_gamma < prctile(acs_prod.cp_gamma, 99);
SimpleMap(acs_prod.cp_gamma(foo), acs_prod(foo,1:3), 'gamma cp (unitless)')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_ACS_gamma_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_ACS_gamma_map']), 'fig')
close figure 1

% acs_prod.chl_Halh = [];
% acs_prod.HH_G50 = [];

filename = sprintf('%s_InLine_ACS_%s_%s_ProductsFull_v%s', official_name, ...
  datetime(min(acs_prod.dt), 'Format', 'yyyyMMdd'),...
  datetime(max(acs_prod.dt), 'Format', 'yyyyMMdd'),...
  datetime('today', 'Format', 'yyyyMMdd'));

% rename variables for .mat
acs_prod = renamevars(acs_prod, {'t','s','cdomf'}, {'sst','sss','fdom'});

% save AC prod
fprintf('Export to mat and csv... ');
save(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, filename), 'acs_prod');
writetable(acs_prod, fullfile(ila.instrument.(list_leg.instru{i}).path.prod, [filename '.csv']));
fprintf('Done\n');

% rename variables for SeaBASS
acs_prod = renamevars(acs_prod, {'sst','sss','fdom'}, {'t','s','cdomf'});

% keep only old variables
acs_prod = acs_prod(:, contains(acs_prod.Properties.VariableNames, {'dt', 'lat', 'lon', ...
  't','s','cdomf','Chl_lineheight','POC_cp','cp_gamma','ap_n','cp_n','flag_bit'}) & ...
  ~contains(acs_prod.Properties.VariableNames, {'agaus'}));

% export product to SeaBASS format
filename = [filename '.sb'];
filename = strrep(filename, '_ProductsFull_', '_ProductsLegacy_');
exportSeaBASS(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, filename), ila.meta, acs_prod);
fprintf('%s saved\n', filename)

% rename variables for .mat
acs_prod = renamevars(acs_prod, {'t','s','cdomf'}, {'sst','sss','fdom'});

% save AC prod
fprintf('Export to mat and csv... ');
filename = strrep(filename, '.sb', '');
save(fullfile(ila.instrument.(list_leg.instru{i}).path.prod, filename), 'acs_prod');
writetable(acs_prod, fullfile(ila.instrument.(list_leg.instru{i}).path.prod, [filename '.csv']));
fprintf('Done\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LISST100X1183 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_leg = table();
list_leg.dt = [...
  datetime(2025,5,5) datetime(2025,6,15);
  ];
list_leg.experiment = {'PVST_SOPACE';};
list_leg.cruise = {'TN444'};

ila.cfg.instruments2run = {'LISST100X1183'};
lisst = [];
for i=1:size(list_leg, 1)
  ila.cfg.days2run = list_leg.dt(i, 1):list_leg.dt(i, 2);
  
  % populate ila.instrument
  ila.Read('prod');

  % update metadata for each leg
  ila.meta.experiment = list_leg.experiment{i};
  ila.meta.cruise = list_leg.cruise{i};

  % % remove flagged data
  % flag = read_flagbit(ila.instrument.(list_leg.instru{:}).prod.p.flag_bit, 'BB');

  % extract HyperBB data from obj
  lisst_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.p;
  lisst_temp = round_timestamp(lisst_temp);
  
  % build hbb table: merge lat, lon, sst, sss
  replace_consecutive_nan = 3*60; % 3h
  lisst_temp = merge_timeseries(lisst_temp, tsg, {'lat', 'lon', 'sst', 'sss'});
  lisst_temp = merge_timeseries(lisst_temp, latlon, {'lat', 'lon'}, '', replace_consecutive_nan);
  
  % Remove NaN
  if any(any(isnan([lisst_temp.lat lisst_temp.lon]),2))
    warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
      sum(any(isnan([lisst_temp.lat lisst_temp.lon]),2)), replace_consecutive_nan/60)
    lisst_temp(any(isnan([lisst_temp.lat lisst_temp.lon]),2), :) = [];
  end

  % add units and precision
  lisst_temp.Properties.VariableNames = {'dt','lat','lon','sst','sss','cp670','VSF670','VSF670_sd','VD','VD_sd','PSD','VSD','bincount'};
  lisst_temp.Properties.VariableUnits = {'','degrees','degrees','degreesC','PSU','1/m','1/m/sr','1/m/sr','uL/L','uL/L','nb/mL/micron','ppm/m','none'};
  lisst_temp.Properties.VariableDescriptions = {'','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%.4f','%i'};
  lisst_parameters = ila.instrument.(ila.cfg.instruments2run{:}).prod.p.Properties.CustomProperties;
  lisst_temp.Properties.CustomProperties = lisst_parameters;

  % add units and precision
  lisst4seabass = lisst_temp(:, ~contains(lisst_temp.Properties.VariableNames, {'PSD','VSD'}));
  lisst4seabass = renamevars(lisst4seabass, {'sst', 'sss','VD','VD_sd'}, {'t','s','PSD','PSD_sd'});
  
  % sort by date
  lisst4seabass = sortrows(lisst4seabass, 'dt'); 
  lisst4seabass(lisst4seabass.cp670 < 0, :) = [];
  lisst4seabass(any(lisst4seabass.VSF670 < 0, 2), :) = [];
  lisst4seabass(all(isnan(lisst4seabass.VSF670), 2), :) = [];
  
  %%% HBB 3D plots %%%
  save_figures = true;
  ila.SpectralQC('LISST', {'prod'}, save_figures); % AC or BB
  close all
  
  % plot time series
  ila.visProd_timeseries()
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_LISST100X_VSD']), 'jpg')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_LISST100X_VSD']), 'fig')
  close figure 103
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_LISST100X_PSD']), 'jpg')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_LISST100X_PSD']), 'fig')
  close figure 102
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_LISST100X_VSF']), 'jpg')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_LISST100X_VSF']), 'fig')
  close figure 101
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_LISST100X_cp']), 'jpg')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_LISST100X_cp']), 'fig')
  close figure 80
  
  SimpleMap(lisst4seabass.cp670, lisst4seabass(:,1:3), 'c_p (670 nm) [m^{-1}]')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_LISST100X_cp_map']), 'jpg')
  saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [list_leg.experiment{i} '_' list_leg.cruise{i} '_LISST100X_cp_map']), 'fig')
  close figure 1
  
  filename_part = sprintf('%s_%s_InLine_%s_%s_%s_Particulate_v%s.sb', list_leg.experiment{i}, list_leg.cruise{i}, ila.cfg.instruments2run{:}, ...
    datetime(min(lisst4seabass.dt), 'Format', 'yyyyMMdd'), datetime(max(lisst4seabass.dt), 'Format', 'yyyyMMdd'), ...
    datetime('today', 'Format', 'yyyyMMdd'));
  
  % export product to SeaBASS format
  ila.meta.documents = sprintf('%s_LISST100X_ProcessingReport_v%s.pdf', official_name, datetime(2025,3,1, 'Format', 'yyyyMMdd'));
  ila.meta.calibration_files = 'LISST100X1183_20230215_factory_zsc.asc,LISST100X1183_20230215_InstrumentData.txt,LISST100X1183_20230215_Lisst.ini,LISST100X1183_20230215_Ringarea.asc';
  exportSeaBASS(fullfile(ila.instrument.FLOW.path.prod, filename_part),...
      ila.meta,...
      lisst4seabass,...
      {'',string(strrep(cellstr(num2str(ila.instrument.(ila.cfg.instruments2run{:}).theta(:), '_%.3fang')),' ',''))', ...
      string(strrep(cellstr(num2str(ila.instrument.(ila.cfg.instruments2run{:}).theta(:), '_%.3fang')),' ',''))',...
      string(strrep(cellstr(num2str(ila.instrument.(ila.cfg.instruments2run{:}).diameters(:), '_%.3fsize')),' ',''))',...
      string(strrep(cellstr(num2str(ila.instrument.(ila.cfg.instruments2run{:}).diameters(:), '_%.3fsize')),' ',''))', ''});
  fprintf('%s saved\n', filename_part)
  
  lisst = [lisst; lisst_temp];
end

lisst.Properties.CustomProperties = lisst_parameters;

% sort by date
lisst = sortrows(lisst, 'dt');
ila.cfg.days2run = list_leg.dt(1, 1):list_leg.dt(end, 2);
ila.instrument.(ila.cfg.instruments2run{:}).prod.p = renamevars(lisst, {'cp670','VSF670','VSF670_sd'},{'cp','betap','betap_sd'});

% plot time series
ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_LISST100X_VSD']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_LISST100X_VSD']), 'fig')
close figure 103
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_LISST100X_PSD']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_LISST100X_PSD']), 'fig')
close figure 102
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_LISST100X_VSF']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_LISST100X_VSF']), 'fig')
close figure 101
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_LISST100X_cp']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_LISST100X_cp']), 'fig')
close figure 80

SimpleMap(lisst.cp670, lisst(:,1:3), 'c_p (670 nm) [m^{-1}]')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_LISST100X_cp_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [official_name '_LISST100X_cp_map']), 'fig')
close figure 1

filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', official_name, ila.cfg.instruments2run{:}, ...
  datetime(min(lisst.dt), 'Format', 'yyyyMMdd'), datetime(max(lisst.dt), 'Format', 'yyyyMMdd'), ...
  datetime('today', 'Format', 'yyyyMMdd'));
  
% save LISST100X prod
fprintf('Export to mat and csv... ');
filename = strrep(filename, '.sb', '');
save(fullfile(ila.instrument.FLOW.path.prod, filename), 'lisst', 'lisst_parameters');
writetable(lisst, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
fprintf('Done\n');














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
if ~isfolder([ila.instrument.FLOW.path.prod 'Graphs'])
  mkdir([ila.instrument.FLOW.path.prod 'Graphs'])
end

% whenever TSG is processed load with
load([ila.instrument.TSG.path.prod cruise '_InLine_TSG_prod.mat'])

% load lat lon vector
% load([[ila.instrument.FLOW.path.prod cruise '_LatLon.mat'])
% latlon=table(nav_data.dt_utc, nav_data.lat, nav_data.lon,'VariableNames',{'dt','lat','lon'});
% % [tsg2, index] = unique(tsg.dt); 
% % tsg = tsg(index,:);
% % yi = interp1(tsg2, y(index), tsg.dt);

% Update cfg
ila.cfg.write.mode = 'One day one file';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TSG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila.cfg.instruments2run = {'TSG'}; % {'FLOW', 'TSG', 'BB31502','PAR', 'WSCD1082P'}
ila.cfg.write.skip = {}; % {'FLOW', 'TSG', 'BB31502','PAR', 'WSCD1082P'}
ila.cfg.days2run = datenum(2020,12,24):datenum(2021,05,09);

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

ila.visProd_timeseries()
saveGraph([ila.instrument.TSG.path.prod 'Graphs' filesep cruise '_tsg'], 'jpg')
close figure 90

SimpleMap(tsg.sst, tsg(:,1:3), 'SST (°C)')
saveGraph([ila.instrument.TSG.path.prod 'Graphs' filesep cruise '_TSG_SST_map'], 'jpg')
close figure 1

SimpleMap(tsg.sss, tsg(:,1:3), 'SSS (PSU)')
saveGraph([ila.instrument.TSG.path.prod 'Graphs' filesep cruise '_TSG_SSS_map'], 'jpg')
close figure 1

% save TSG prod
fprintf('Export to mat and csv... ');
save([ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
    cruise '_InLine_TSG_prod'], 'tsg');
writetable(tsg, [ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
    cruise '_InLine_TSG_prod.csv']);
fprintf('Done\n');

% % add missing lat/lon from meteo file
load([ila.instrument.FLOW.path.prod cruise '_meteo_prod.mat'])
load([ila.instrument.FLOW.path.prod cruise '_InLine_TSG_prod.mat'])
meteo_merged_tsg = MergeTimeSeries(tsg.dt, {meteo_ftp}, {'meteo_ftp'});
tsg.lat(isnan(tsg.lat)) = meteo_merged_tsg.lat(isnan(tsg.lat));
tsg.lon(isnan(tsg.lon)) = meteo_merged_tsg.lon(isnan(tsg.lon));
% save TSG prod
fprintf('Export to mat and csv... ');
save([ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
    cruise '_InLine_TSG_prod'], 'tsg');
writetable(tsg, [ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod ...
    cruise '_InLine_TSG_prod.csv']);
fprintf('Done\n');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WSCD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
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
saveGraph([ila.instrument.(ila.cfg.instruments2run{:}).path.prod ...
  'Graphs' filesep cruise '_FDOM_timeseries'], 'jpg')
close figure 94

SimpleMap(wscd.fdom, wscd(:,1:3), 'WSCD fdom ppb')
saveGraph([ila.instrument.(ila.cfg.instruments2run{:}).path.prod 'plots' filesep cruise '_WSCD_fdom_map'], 'jpg')
close figure 1

% filename = [ila.instrument.(cell2mat(ila.cfg.instruments2run)).path.prod cruise ...
%   '_InLine_' cell2mat(ila.cfg.instruments2run) '_Product_v' datestr(now, 'yyyymmdd') '.sb'];

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
save([ila.instrument.(ila.cfg.instruments2run{:}).path.prod ...
    cruise '_InLine_WSCD_prod'], 'wscd');
writetable(wscd, [ila.instrument.(ila.cfg.instruments2run{:}).path.prod ...
    cruise '_InLine_WSCD_prod.csv']);
fprintf('Done\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
% ila.cfg.instruments2run = {'PAR'}; % {'FLOW', 'TSG', 'BB3','PAR', 'WSCD1082P'}
% ila.cfg.days2run = datenum(2020,12,24):datenum(2021,5,9);
% 
% % populate ila.instrument
% ila.Read('prod');
% 
% % interpolate SST / SSS / LatLon
% par_temp = ila.instrument.(cell2mat(ila.cfg.instruments2run)).prod.a;
% par_temp.dt = datetime(par_temp.dt,'ConvertFrom','datenum');
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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BB3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
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

%%% BB 3D plots %%%
save_figures = true;
ila.DiagnosticPlot('BB', {'prod'}, save_figures); % AC or BB
close all

% plot time series
ila.visProd_timeseries()
saveGraph([ila.instrument.TSG.path.prod 'plots' filesep cruise '_BB3_BBparticulate_timeseries'], 'jpg')
close figure 86
saveGraph([ila.instrument.TSG.path.prod 'plots' filesep cruise '_BB3_POCparticulate_timeseries'], 'jpg')
close figure 85
saveGraph([ila.instrument.TSG.path.prod 'plots' filesep cruise '_BB3_BBdissolved_timeseries'], 'jpg')
close figure 87

SimpleMap(bb3.bbp(:,2), bb3(:,1:3), 'bbp (532 nm) [m^-^1]')
saveGraph([ila.instrument.TSG.path.prod 'plots' filesep cruise '_bb3_bbp_map'], 'jpg')
close figure 1

filename = [ila.instrument.(ila.cfg.instruments2run{:}).path.prod cruise ...
  '_InLine_' ila.cfg.instruments2run{:} '_Particulate_v' datestr(now, 'yyyymmdd') '.sb'];

% export product to SeaBASS format
ila.meta.documents = [cruise '_BB3_ProcessingReport.pdf'];
ila.meta.calibration_files = 'BB3-1502_(470-532-650nm)_CharSheet.pdf';
exportSeaBASS(filename,...
    ila.meta,...
    bb3,...
    {string(bb3_lambda(1:2)),string(bb3_lambda(1:2)),string(bb3_lambda(1:2)),''});
sprintf('%s_InLine_%s_Particulate.sb saved', cruise, ila.cfg.instruments2run{:})

bb3.Properties.VariableNames = {'dt', 'lat', 'lon', 'sst', 'sss', 'VSF124', 'bbp', 'VSF124_sd','bincount'};
bb3.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'm^-1.sr^-1', 'm^-1', 'm^-1.sr^-1', 'none'};
bb3.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f'};

% save BB3 prod
fprintf('Export to mat and csv... ');
save([ila.instrument.(ila.cfg.instruments2run{:}).path.prod ...
    cruise '_InLine_BB3_particulate'], 'bb3', 'bb3_lambda');
writetable(bb3,[ila.instrument.(ila.cfg.instruments2run{:}).path.prod ...
    cruise '_InLine_BB3_particulate.csv']);
fprintf('Done\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HBB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila = InLineAnalysis(['cfg/' cruise '_cfg.m']);
ila.cfg.instruments2run = {'HBB'}; % {'FLOW', 'TSG', 'BB31502','PAR', 'WSCD1082', 'HBB'}
ila.cfg.days2run = datenum(2020,12,24):datenum(2021,5,9);

% populate ila.instrument
ila.Read('prod');

% interpolate SST / SSS / LatLon
hbb_lambda = ila.instrument.(ila.cfg.instruments2run{:}).lambda;
hbb_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.p;
hbb_temp.dt = datetime(hbb_temp.dt,'ConvertFrom','datenum');
% latlon_interp = interp1(latlon.dt, [latlon.lat, latlon.lon], bb3_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data
tsg_interp = interp1(tsg.dt, [tsg.lat, tsg.lon, tsg.sst, tsg.sss], hbb_temp.dt, 'linear', 'extrap'); % extrap needed for first minute of data

% bb3 Products
hbb = table(hbb_temp.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4), hbb_temp.betap, bb3_temp.bbp, bb3_temp.betap_sd, bb3_temp.betap_n,...
             'VariableNames', {'dt', 'lat', 'lon', 't', 's', 'VSF_124ang','bbp','VSF_124ang_sd','bincount'});
hbb.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', '1/m/sr', '1/m', '1/m/sr', 'none'};
hbb.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f'};

% sort by date
hbb = sortrows(hbb, 'dt'); 
hbb.bbp(hbb.bbp<0)=NaN;
hbb.VSF_124ang_sd(hbb.VSF_124ang<0)=NaN;
hbb.VSF_124ang(hbb.VSF_124ang<0)=NaN;
hbb(all(isnan(hbb.bbp),2),:)=[];

% remove dobious values
hbb(hbb.dt > datetime(2021,2,17,0,0,0) & hbb.dt < datetime(2021,2,23,18,30,0), :) = [];
hbb(hbb.dt > datetime(2021,3,23,21,25,0) & hbb.dt < datetime(2021,3,25,0,0,0), :) = [];
hbb(hbb.dt > datetime(2021,3,25,3,14,0) & hbb.dt < datetime(2021,3,31,0,0,0), :) = [];

%%% BB 3D plots %%%
save_figures = true;
ila.DiagnosticPlot('BB', {'prod'}, save_figures); % AC or BB
close all

% plot time series
ila.visProd_timeseries()
saveGraph([ila.instrument.TSG.path.prod 'plots' filesep cruise '_HBB_BBparticulate_timeseries'], 'jpg')
close figure 86
saveGraph([ila.instrument.TSG.path.prod 'plots' filesep cruise '_HBB_POCparticulate_timeseries'], 'jpg')
close figure 85
saveGraph([ila.instrument.TSG.path.prod 'plots' filesep cruise '_HBB_BBdissolved_timeseries'], 'jpg')
close figure 87

SimpleMap(hbb.bbp(:,hbb_lambda == 530), hbb(:,1:3), 'bbp (530 nm) [m^-^1]')
saveGraph([ila.instrument.TSG.path.prod 'plots' filesep cruise '_hbb_bbp_map'], 'jpg')
close figure 1

filename = [ila.instrument.(ila.cfg.instruments2run{:}).path.prod cruise ...
  '_InLine_' ila.cfg.instruments2run{:} '_Particulate_v' datestr(now, 'yyyymmdd') '.sb'];

% export product to SeaBASS format
ila.meta.documents = [cruise '_HBB_ProcessingReport.pdf'];
ila.meta.calibration_files = 'HBB8004_CharSheet.pdf';
exportSeaBASS(filename,...
    ila.meta,...
    hbb,...
    {string(lambda(1:2)),string(lambda(1:2)),string(lambda(1:2)),''});
sprintf('%s_InLine_%s_Particulate.sb saved', cruise, ila.cfg.instruments2run{:})

hbb.Properties.VariableNames = {'dt', 'lat', 'lon', 'sst', 'sss', 'VSF124', 'bbp', 'VSF124_sd','bincount'};
hbb.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', 'm^-1.sr^-1', 'm^-1', 'm^-1.sr^-1', 'none'};
hbb.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.2f'};

% save HBB prod
fprintf('Export to mat and csv... ');
save([ila.instrument.(ila.cfg.instruments2run{:}).path.prod ...
    cruise '_InLine_HBB_particulate'], 'hbb', 'lambda');
writetable(hbb,[ila.instrument.(ila.cfg.instruments2run{:}).path.prod ...
    cruise '_InLine_HBB_particulate.csv']);
fprintf('Done\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_leg = {...
    datenum(2020,12,24):datenum(2021,02,05);...
    datenum(2021,02,17):datenum(2021,3,17);...
    datenum(2021,03,20):datenum(2021,4,7);...
    datenum(2021,04,11):datenum(2021,4,20);...
    datenum(2021,04,26):datenum(2021,5,9);...
    };

list_dev = {...
    [path_dev filesep 'acs057_20200129.dev'];...
    [path_dev filesep 'acs057_20200129.dev'];...
    [path_dev filesep 'acs057_20200129.dev'];...
    [path_dev filesep 'acs057_20200129.dev'];...
    [path_dev filesep 'acs057_20200129.dev'];...
%     [path_dev filesep 'acs301_20160704_20160720.dev'];...
    };

list_instru = {...
  'ACS57';...
  'ACS57';...
  'ACS57';...
  'ACS57';...
  'ACS57';...
%   'ACS301';...
  };

data_AC = struct('particulate', [], 'product', []);
acs = [];
if exist('bb3', 'var')
  bb3_b = [];
end
if exist('hbb', 'var')
  hbb_b = [];
end

for i=1:size(list_instru,1)
  ila = InLineAnalysis(['cfg' filesep cruise '_cfg.m']);
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
  data_AC.product.(ref) = table(AC.dt, tsg_interp(:,1), tsg_interp(:,2), tsg_interp(:,3), tsg_interp(:,4),...
               'VariableNames', {'dt', 'lat', 'lon', 't', 's'});
  data_AC.product.(ref) = [data_AC.product.(ref) AC(:, 8:end)];
  data_AC.product.(ref).Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU', ...
              '1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m','1/m', ...
              'ug/L','unitless','ug/L','unitless','unitless','ug/L','microns','??'};
  data_AC.product.(ref).Properties.VariableDescriptions = [{''}, repmat({'%.4f'}, 1, size(data_AC.product.(ref),2) - 1)];

  data_AC.product.(ref).Properties.VariableNames{...
    strcmp(data_AC.product.(ref).Properties.VariableNames, 'poc')} = 'POC_cp';
  data_AC.product.(ref).Properties.VariableNames{...
    strcmp(data_AC.product.(ref).Properties.VariableNames, 'gamma')} = 'cp_gamma';
  data_AC.product.(ref).Properties.VariableNames{...
    strcmp(data_AC.product.(ref).Properties.VariableNames, 'chl_ap676lh')} = 'Chl_lineheight';

  % Remove NaN and aberrant 
  data_AC.product.(ref).POC_cp(data_AC.product.(ref).POC_cp < 0) = NaN;
  data_AC.product.(ref).Chl_lineheight(data_AC.product.(ref).Chl_lineheight < 0) = NaN;
  data_AC.product.(ref).cp_gamma(data_AC.product.(ref).cp_gamma < 0) = NaN;
  data_AC.product.(ref)(all(isnan(table2array(data_AC.product.(ref)(:,6:8))),2),:)=[];
  
  % interpolate over BB3 wavelength to QC/QA with bbp/bp ratio
  bp = data_AC.particulate.(ref).cp - data_AC.particulate.(ref).ap;
  if exist('bb3_b', 'var')
    bb3_b = [bb3_b; datenum(data_AC.particulate.(ref).dt) interp1(ila.instrument.(list_instru{i}).lambda_a, ...
      bp',bb3_lambda,'linear')'];
  end
  if exist('hbb_b', 'var')
    hbb_b = [hbb_b; datenum(data_AC.particulate.(ref).dt) interp1(ila.instrument.(list_instru{i}).lambda_a, ...
      bp',hbb_lambda,'linear')'];
  end
  
  %%% AC 3D plots %%%
  save_figures = true;
  ila.DiagnosticPlot('AC', {'prod'}, save_figures); % AC or BB
  close all
  
  filename = [ila.instrument.(list_instru{i}).path.prod cruise '_InLine_' ref ...
    '_Particulate_v' datestr(now, 'yyyymmdd') '.sb'];
  % export particulate to SeaBASS format
  ila.meta.documents = [cruise '_ACS_ProcessingReport.pdf'];
  [~, calfile] = fileparts(list_dev{i});
  ila.meta.calibration_files = [calfile '.dev'];
  exportSeaBASS(filename,...
      ila.meta,...
      data_AC.particulate.(ref),...
      {string(ila.instrument.(list_instru{i}).lambda_a),...
      string(ila.instrument.(list_instru{i}).lambda_a),...
      string(ila.instrument.(list_instru{i}).lambda_c),...
      string(ila.instrument.(list_instru{i}).lambda_c),''});
  sprintf('%s_InLine_%s_Particulate.sb saved', cruise, ref)

  % % export product to SeaBASS format
  % exportSeaBASS([ila.instrument.(list_instru{i}).path.prod cruise '_InLine_' ref '_Product.sb'],...
  %     ila.meta,...
  %     data_AC.product.(ref),...
  %     {'', '', ''});
  % sprintf('%s_InLine_%s_Product.sb saved', cruise, ref)

  % ACS merged prod
  acs = [acs; data_AC.product.(ref)];
end

% sort AC and save prod file .mat
[~,b] = sort(acs.dt); % sort by date
acs = acs(b,:);

ila.visProd_timeseries()
saveGraph([ila.instrument.(list_instru{i}).path.prod 'plots' filesep cruise '_ACS_prod_timeseries'], 'jpg')
close figure 78
saveGraph([ila.instrument.(list_instru{i}).path.prod 'plots' filesep cruise '_ACS_prod_regressions'], 'jpg')
close figure 77

SimpleMap(acs.chl_Halh, acs(:,1:3), 'Houskeeper [chl] (mg.m^{-3})')
saveGraph([ila.instrument.TSG.path.prod 'plots' filesep cruise '_ACS_chl_Houskeeper_map'], 'jpg')
close figure 1

SimpleMap(acs.HH_G50, acs(:,1:3), 'H&H phytoplankton G50: cross-sectional area (\mum)')
saveGraph([ila.instrument.TSG.path.prod 'plots' filesep cruise '_ACS_H&H_map'], 'jpg')
close figure 1

SimpleMap(acs.POC_cp, acs(:,1:3), '[POC] cp (mg.m^{-3})')
saveGraph([ila.instrument.TSG.path.prod 'plots' filesep cruise '_ACS_POC_map'], 'jpg')
close figure 1

SimpleMap(acs.Chl_lineheight, acs(:,1:3), 'a_{p676} line height [chl a] (mg.m^{-3})')
saveGraph([ila.instrument.TSG.path.prod 'plots' filesep cruise '_ACS_chl_map'], 'jpg')
close figure 1

SimpleMap(acs.cp_gamma, acs(:,1:3), 'gamma cp (unitless)')
saveGraph([ila.instrument.TSG.path.prod 'plots' filesep cruise '_ACS_gamma_map'], 'jpg')
close figure 1

% save AC prod
fprintf('Export to mat and csv... ');
save([ila.instrument.(list_instru{i}).path.prod cruise '_InLine_ACS_prod'], 'acs');
writetable(acs, [ila.instrument.(list_instru{i}).path.prod cruise '_InLine_ACS_prod.csv']);
fprintf('Done\n');

%% bbp/bp graph to validate BB3 
bp_interp = interp1(bb3_b(:,1), bb3_b(:,2:end), ...
  datenum(bb3.dt), 'linear', 'extrap'); % extrap needed for first minute of data

% get spectrumRGB
C = reshape(spectrumRGB(bb3_lambda), max(size(bb3_lambda)),  3);

figure(); hold on
scatter(bb3.dt, bb3.bbp./bp_interp, 15, C, 'filled')
xlim([min(bb3.dt) max(bb3.dt)])
ylabel('{\itbb_p}/{\itb_p}')
legend(cellfun(@(x) [x 'nm'], cellstr(num2str(bb3_lambda')), 'un', 0), 'FontSize', 14)
set(gca, 'FontSize', 13)
saveGraph([ila.instrument.TSG.path.prod 'plots' filesep cruise '_BB3_bbp_over_bp'], 'jpg')
close figure 1

%% bbp/bp graph to validate HBB data
bp_interp = interp1(hbb_b(:,1), hbb_b(:,2:end), ...
  datenum(hbb.dt), 'linear', 'extrap'); % extrap needed for first minute of data

% get spectrumRGB
C = reshape(spectrumRGB(hbb_lambda), max(size(hbb_lambda)),  3);

figure(); hold on
scatter(hbb.dt, hbb.bbp./bp_interp, 15, C, 'filled')
xlim([min(hbb.dt) max(hbb.dt)])
ylabel('{\itbb_p}/{\itb_p}')
legend(cellfun(@(x) [x 'nm'], cellstr(num2str(hbb_lambda')), 'un', 0), 'FontSize', 14)
set(gca, 'FontSize', 13)
saveGraph([ila.instrument.TSG.path.prod 'plots' filesep cruise '_HBB_bbp_over_bp'], 'jpg')
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
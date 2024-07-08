% Consolidate daily processing, merge variables, and export to SeaBASS
% author: Guillaume Bourdin
% created: May 27, 2024

%% Import data

% cd('/Volumes/Data2/TaraEuropa/InLineAnalysis-master')
cd('/Users/gui/Documents/MATLAB/InLineAnalysis/InLineAnalysis-master')

cruise = 'Apero';
% Load InLineAnalysis and the configuration
ila = InLineAnalysis(['cfg' filesep cruise '_cfg.m']);

path_dev = strrep(ila.instrument.FLOW.path.prod, ...
  'prod', 'DeviceFiles');

% create Graph folder if it doesn't exist
if ~isfolder([ila.instrument.FLOW.path.prod '/plots'])
  mkdir([ila.instrument.FLOW.path.prod '/plots'])
end

% whenever TSG and fdom are processed, load with:
% load([ila.instrument.FLOW.path.prod cruise '_InLine_TSG_prod.mat'])
% load([ila.instrument.SUVF6253.path.prod cruise '_InLine_SUVF6253_prod.mat'])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TSG & GPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila.cfg.days2run = datenum(2023,6,4):datenum(2023,7,11);
ila.cfg.instruments2run = {'SBE459999'};

% populate ila.instrument
ila.Read('prod');

% extract TSG data from obj
tsg_temp = ila.instrument.('SBE459999').prod.a;
tsg_temp.dt = datetime(tsg_temp.dt,'ConvertFrom','datenum');
tsg_temp = round_timestamp(tsg_temp);

% build TSG table: linearly interpolate missing lat, lon when missing data < replace_consecutive_nan
replace_consecutive_nan = 3*60; % 3h
tsg = merge_timeseries(tsg_temp, tsg_temp, {'lat', 'lon'}, '', replace_consecutive_nan);
% Remove NaN
if any(any(isnan([tsg.lat tsg.lon]),2))
  warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
    sum(any(isnan([tsg.lat tsg.lon]),2)), replace_consecutive_nan/60)
  tsg(any(isnan([tsg.lat tsg.lon]),2), :) = [];
end
tsg(:, contains(tsg.Properties.VariableNames, {'_sd', '_n'})) = [];
tsg = movevars(tsg, {'lat', 'lon'}, 'After','dt');
tsg = renamevars(tsg,{'t1','s'},{'sst','sss'});

% add units and precision
tsg.Properties.VariableUnits = {'', 'degrees', 'degrees', 'degreesC', 'PSU'};
tsg.Properties.VariableDescriptions = {'', '%.4f', '%.4f', '%.4f', '%.4f'};

% sort by date
tsg = sortrows(tsg, 'dt');
ila.instrument.('SBE459999').prod.a = renamevars(ila.instrument.('SBE459999').prod.a, {'t1'}, {'sst'});

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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUVF6253 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila.cfg.days2run = datenum(2023,6,4):datenum(2023,7,11);
ila.cfg.instruments2run = {'SUVF6253'};

% populate ila.instrument
ila.Read('prod');

% extract SUVF data from obj
suvf_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.a;
suvf_temp.dt = datetime(suvf_temp.dt,'ConvertFrom','datenum');
suvf_temp = round_timestamp(suvf_temp);

% build suvf table: merge lat, lon, sst, sss
replace_consecutive_nan = 3*60; % 3h
suvf = merge_timeseries(suvf_temp, tsg, {'lat', 'lon', 'sst', 'sss'}, '', replace_consecutive_nan);

% Remove NaN
if any(isnan([suvf.lat suvf.lon]),2)
  suvf(any(isnan([suvf.lat suvf.lon]),2), :) = [];
  warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
    sum(any(isnan([suvf.lat suvf.lon]),2)), replace_consecutive_nan/60)
end

% add units and precision
suvf = renamevars(suvf, {'sst','sss','fdom_avg_n'}, {'t','s','bincount'});
suvf = removevars(suvf, {'swt','swt_avg_sd','swt_avg_n','spd1','spd1_avg_sd','spd1_avg_n'});
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

% % export product to SeaBASS format
% ila.meta.documents = [cruise '_SUVF_ProcessingReport_V2.pdf'];
% ila.meta.calibration_files = 'SUVF6244_CharSheet.pdf';
% exportSeaBASS(fullfile(ila.instrument.FLOW.path.prod, filename),...
%     ila.meta,...
%     suvf,...
%     {'', '', ''});
% sprintf('%s_InLine_%s_Product.sb saved', cruise, cell2mat(ila.cfg.instruments2run))

suvf = renamevars(suvf, {'t','s','bincount'}, {'sst','sss','fdom_n'});

% save SUVF prod
fprintf('Export to mat and csv... ');
filename = strrep(filename, '.sb', '');
save(fullfile(ila.instrument.FLOW.path.prod, filename), 'suvf');
writetable(suvf, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
fprintf('Done\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BB31502 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ila.cfg.days2run = datenum(2023,6,4):datenum(2023,7,11);
ila.cfg.instruments2run = {'BB31502'};

% populate ila.instrument
ila.Read('prod');

% extract HyperBB data from obj
bb3_lambda = ila.instrument.(ila.cfg.instruments2run{:}).lambda;
bb3_temp = ila.instrument.(ila.cfg.instruments2run{:}).prod.p;
bb3_temp.dt = datetime(bb3_temp.dt,'ConvertFrom','datenum');
bb3_temp = round_timestamp(bb3_temp);

% build hbb table: merge lat, lon, sst, sss
replace_consecutive_nan = 3*60; % 3h
bb3_temp = merge_timeseries(bb3_temp, tsg, {'lat', 'lon', 'sst', 'sss'}, '', replace_consecutive_nan);

% Remove NaN
if any(any(isnan([bb3_temp.lat bb3_temp.lon]),2))
  warning('%i row with missing lat/lon (longer than %ih consecutive): deleted', ...
    sum(any(isnan([bb3_temp.lat bb3_temp.lon]),2)), replace_consecutive_nan/60)
  bb3_temp(any(isnan([bb3_temp.lat bb3_temp.lon]),2), :) = [];
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

%%% BB3 3D plots %%%
save_figures = true;
ila.SpectralQC('BB', {'prod'}, save_figures); % AC or BB
close all

% plot time series
ila.visProd_timeseries()
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_BB3_gammabbp_bbp550']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_BB3_gammabbp_bbp550']), 'fig')
close figure 11
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_BB3_BBparticulate_timeseries']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_BB3_BBparticulate_timeseries']), 'fig')
close figure 21
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_BB3_POCparticulate_timeseries']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_BB3_POCparticulate_timeseries']), 'fig')
close figure 20
% saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_BB3_BBdissolved_timeseries']), 'jpg')
% saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_BB3_BBdissolved_timeseries']), 'fig')
% close figure 25

lambda_id = abs(bb3_lambda - 530) == min(abs(bb3_lambda - 530));
SimpleMap(bb3.bbp(:,lambda_id), bb3(:,1:3), ['bbp (' num2str(bb3_lambda(lambda_id)) 'nm) [m^-^1]'] )
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_bb3_bbp_map']), 'jpg')
saveGraph(fullfile(ila.instrument.FLOW.path.prod, 'plots', [cruise '_bb3_bbp_map']), 'fig')
close figure 1

filename = sprintf('%s_InLine_%s_%s_%s_Product_v%s.sb', cruise, ila.cfg.instruments2run{:}, ...
  datetime(min(bb3.dt), 'Format', 'yyyyMMdd'), datetime(max(bb3.dt), 'Format', 'yyyyMMdd'), ...
  datetime('today', 'Format', 'yyyyMMdd'));

% export product to SeaBASS format
ila.meta.documents = sprintf('%s_%s_ProcessingReport_v%s.pdf', cruise, ...
  ila.cfg.instruments2run, datetime('today', 'Format', 'yyyyMMdd'));
ila.meta.calibration_files = [ila.cfg.instruments2run '_CharSheet.pdf'];
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

% save HBB prod
fprintf('Export to mat and csv... ');
filename = strrep(filename, '.sb', '');
save(fullfile(ila.instrument.FLOW.path.prod, filename), 'bb3', 'bb3_lambda');
writetable(bb3, fullfile(ila.instrument.FLOW.path.prod, [filename '.csv']));
fprintf('Done\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_leg = {...
    datenum(2023,6,4):datenum(2023,7,11);
    };

list_dev = repmat({fullfile(path_dev, 'ACS-348_20230504.dev')}, size(list_leg, 1), 1);
list_instru = repmat({'ACS348'}, size(list_leg, 1), 1);

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
  % flag_info = InlineFlagInfo('ACS');

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
  replace_consecutive_nan = 3*60; % 3h
  AC = merge_timeseries(AC, suvf, {'fdom'}, '', replace_consecutive_nan);
  AC = merge_timeseries(AC, tsg, {'lat', 'lon', 'sst', 'sss'}, '', replace_consecutive_nan);

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









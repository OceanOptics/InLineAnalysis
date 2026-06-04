function [ data ] = importInlininoNMEA( filename, verbose )
  % importInlininoNMEA Import NMEA data from csv files
  % Author: Guillaume Bourdin
  % Date: August 1st, 2022
  %
  % Input:
  %   - filename: <char> filename including full path
  %   - verbose (optional)
  % 
  % Example: [ data, lambda] = importInlininoNMEA( filename, verbose )
  %%
  if nargin < 2; verbose = false; end
  if verbose
    foo = strsplit(filename, '/');
    fprintf('Importing %s ... ', foo{end});
  end
  
  % % Get date from filename
  % s = strsplit(filename, '_');
  % dt_ref = s{end-1};
  
  % Open file
  fid=fopen(filename);
  if fid==-1
    error('Unable to open file: %s', filename);
  end
  
  % Get header
  hd = strip(strsplit(fgetl(fid), ','));
  if any(cellfun('isempty', hd))
    emptyvarname = strcat({'Var'}, cellstr(num2str((1:size(hd, 2))')));
    if any(ismember(emptyvarname, hd))
      duplic = find(ismember(emptyvarname, hd));
      for i = 1:size(duplic,1)
        emptyvarname{duplic(i)} = ['Var', num2str(str2double(emptyvarname{duplic(i)}(end))+100)];
      end
    end
    hd(cellfun('isempty', hd)) = emptyvarname(cellfun('isempty', hd));
  end
  if endsWith(filename, '.log')
    raw_format = true;
  elseif strcmp(hd{2}, 'packet')
    raw_format = true;
  else
    raw_format = false;
  end
  
  if raw_format
    % read text file
    t = string();
    foo = fgetl(fid);
    while ischar(foo)
      t{end+1,:} = foo;
      foo = fgetl(fid);
    end
    fclose(fid);
    % t = readlines(filename);
    % parse RMC or GGA sentences: Recommended Minimum Navigation Information
    t_rmc = t(contains(t,'RMC,'));
    t_gga = t(contains(t,'GGA,'));
    if ~isempty(t_rmc)
      data = parseRMC(t_rmc);
      data = round_timestamp(data, seconds(1));
    elseif ~isempty(t_gga)
      data = parseGGA(t_gga);
      data = round_timestamp(data, seconds(1));
    else
      warning('No RMC or GGA sentences found, datetime, latitude, and longitude data missing, check NMEA sentences')
    end
    % parse MWD sentences: True Wind
    t_mwd = t(contains(t,'MWD,'));
    if ~isempty(t_mwd)
      data_mwd = parseMWD(t_mwd);
      data_mwd = round_timestamp(data_mwd, seconds(1));
      % merge mwd
      if exist('data','var')
        data.true_wind_dir_mwd = interp1(data_mwd.dt, data_mwd.true_wind_dir, data.dt, 'nearest');
        data.true_wind_spd_mwd = interp1(data_mwd.dt, data_mwd.true_wind_spd, data.dt, 'nearest');
      else
        data = data_mwd;
        data = renamevars(data, {'true_wind_dir','true_wind_spd'},{'true_wind_dir_mwd','true_wind_spd_mwd'});
      end
    else
      if exist('data','var')
        data.true_wind_dir_mwd = NaN(size(data.dt));
        data.true_wind_spd_mwd = NaN(size(data.dt));
      end
    end
    % parse MWV sentences: Apparent Wind
    t_mwv = t(contains(t,'MWV,'));
    if ~isempty(t_mwv)
      data_mwv = parseMWV(t_mwv);
      data_mwv = round_timestamp(data_mwv, seconds(1));
      % merge mwv
      if exist('data','var')
        data.true_wind_dir_mwv = interp1(data_mwv.dt, data_mwv.true_wind_dir, data.dt, 'nearest');
        data.true_wind_spd_mwv = interp1(data_mwv.dt, data_mwv.true_wind_spd, data.dt, 'nearest');
        data.apparent_wind_dir = interp1(data_mwv.dt, data_mwv.apparent_wind_dir, data.dt, 'nearest');
        data.apparent_wind_spd = interp1(data_mwv.dt, data_mwv.apparent_wind_spd, data.dt, 'nearest');
      else
        data = data_mwv;
        data = renamevars(data, {'true_wind_dir','true_wind_spd'},{'true_wind_dir_mwv','true_wind_spd_mwv'});
      end
    else
      if exist('data','var')
        data.true_wind_dir_mwv = NaN(size(data.dt));
        data.true_wind_spd_mwv = NaN(size(data.dt));
        data.apparent_wind_dir = NaN(size(data.dt));
        data.apparent_wind_spd = NaN(size(data.dt));
      end
    end
    % parse MMB sentences: Atmospheric Pressure
    t_mmb = t(contains(t,'MMB,'));
    if ~isempty(t_mmb)
      data_mmb = parseMMB(t_mmb);
      data_mmb = round_timestamp(data_mmb, seconds(1));
      % merge mmb
      if exist('data','var')
        data.atm_press = interp1(data_mmb.dt, data_mmb.atm_press, data.dt, 'nearest');
      else
        data = data_mmb;
      end
    else
      if exist('data','var')
        data.atm_press = NaN(size(data.dt));
      end
    end
    % parse VPW sentences: Speed - Measured Parallel to Wind
    t_vpw = t(contains(t,'VPW,'));
    if ~isempty(t_vpw)
      data_vpw = parseVPW(t_vpw);
      data_vpw = round_timestamp(data_vpw, seconds(1));
      % merge vpw
      if exist('data','var')
        data.spd_wind_relative = interp1(data_vpw.dt, data_vpw.spd_wind_relative, data.dt, 'nearest');
      else
        data = data_vpw;
      end
    else
      if exist('data','var')
        data.spd_wind_relative = NaN(size(data.dt));
      end
    end
    % parse VHW sentences: Water speed and heading
    t_vhw = t(contains(t,'VHW,'));
    if ~isempty(t_vhw)
      data_vhw = parseVHW(t_vhw);
      data_vhw = round_timestamp(data_vhw, seconds(1));
      % merge vhw
      if exist('data','var')
        data.heading_true_vhw = interp1(data_vhw.dt, data_vhw.heading_true, data.dt, 'nearest');
        data.heading_magnetic_vhw = interp1(data_vhw.dt, data_vhw.heading_magnetic, data.dt, 'nearest');
        data.spd_water_relative = interp1(data_vhw.dt, data_vhw.spd_water_relative, data.dt, 'nearest');
      else
        data = data_vhw;
        data = renamevars(data, {'heading_true','heading_magnetic'},{'heading_true_vhw','heading_magnetic_vhw'});
      end
    else
      if exist('data','var')
        data.heading_true_vhw = NaN(size(data.dt));
        data.heading_magnetic_vhw = NaN(size(data.dt));
        data.spd_water_relative = NaN(size(data.dt));
      end
    end
    % parse DBT sentences: Depth below transducer
    t_dbt = t(contains(t,'DBT,'));
    if ~isempty(t_dbt)
      data_dbt = parseDBT(t_dbt);
      data_dbt = round_timestamp(data_dbt, seconds(1));
      % merge dbt
      if exist('data','var')
        data.depth_m = interp1(data_dbt.dt, data_dbt.depth_m, data.dt, 'nearest');
      else
        data = data_dbt;
      end
    else
      if exist('data','var')
        data.depth_m = NaN(size(data.dt));
      end
    end
    % parse MTW sentences: Mean Temperature of Water
    t_mtw = t(contains(t,'MTW,'));
    if ~isempty(t_mtw)
      data_mtw = parseMTW(t_mtw);
      data_mtw = round_timestamp(data_mtw, seconds(1));
      % merge mtw
      if exist('data','var')
        data.water_temp = interp1(data_mtw.dt, data_mtw.water_temp, data.dt, 'nearest');
      else
        data = data_mtw;
      end
    else
      if exist('data','var')
        data.water_temp = NaN(size(data.dt));
      end
    end
    % parse MTA sentences: Air Temperature
    t_mta = t(contains(t,'MTA,'));
    if ~isempty(t_mta)
      data_mta = parseMTA(t_mta);
      data_mta = round_timestamp(data_mta, seconds(1));
      % merge mta
      if exist('data','var')
        data.air_temp = interp1(data_mta.dt, data_mta.air_temp, data.dt, 'nearest');
      else
        data = data_mta;
      end
    else
      if exist('data','var')
        data.air_temp = NaN(size(data.dt));
      end
    end
    % parse VTG sentences: Course over ground
    t_vtg = t(contains(t,'VTG,'));
    if ~isempty(t_vtg)
      data_vtg = parseVTG(t_vtg);
      data_vtg = round_timestamp(data_vtg, seconds(1));
      % merge mta
      if exist('data','var')
        data.cog = interp1(data_vtg.dt, data_vtg.cog, data.dt, 'nearest');
        data.magnetic_cog = interp1(data_vtg.dt, data_vtg.magnetic_cog, data.dt, 'nearest');
        data.sog = interp1(data_vtg.dt, data_vtg.sog_kt, data.dt, 'nearest');
      else
        data = data_vtg;
      end
    else
      if exist('data','var')
        data.cog = NaN(size(data.dt));
        data.magnetic_cog = NaN(size(data.dt));
        data.sog = NaN(size(data.dt));
      end
    end
    clear t
    if ~any(strcmp(data.Properties.VariableNames,'dt_instrument'))
      data = addvars(data, NaT(size(data.dt)),'NewVariableNames','dt_instrument','After','dt');
    end
    if ~any(strcmp(data.Properties.VariableNames,'lat'))
      data = addvars(data, NaN(size(data.dt)),'NewVariableNames','lat','After','dt_instrument');
    end
    if ~any(strcmp(data.Properties.VariableNames,'lon'))
      data = addvars(data, NaN(size(data.dt)),'NewVariableNames','lon','After','lat');
    end
    % % remove empty rows
    % data(all(isnat(data.dt) & isnat(data.dt_instrument) & isnan(data.lat) & isnan(data.lon), 2), :) = [];
    % data(all(isnat(data.dt_instrument) & isnan(data.lat) & isnan(data.lon), 2), :) = [];
    % data.dt = datenum(data.dt);
    % data.dt_instrument = datenum(data.dt_instrument);
  else
    % get units skipping empty lines (bug in old Inlinino)
    unit = fgetl(fid);
    while isempty(unit)
      unit = fgetl(fid);
    end
    % get units and lambda
    unit = strip(strsplit(strrep(unit, ',,', ', ,'), ','));
    first_row = strip(strsplit(fgetl(fid), ','));
    pars = repmat({''}, size(hd));
    % rename variables
    for i = 1:size(hd, 2)
      if strcmpi(hd{i}, 'time')
        hd{i} = 'dt';
        pars{i} = '%s';
      elseif strcmp(hd{i}, 'datetime')
        hd{i} = 'dt_instrument';
        pars{i} = '%s';
      elseif strcmpi(hd{i}, 'latitude')
        hd{i} = 'lat';
        pars{i} = '%f';
      elseif strcmpi(hd{i}, 'longitude')
        hd{i} = 'lon';
        pars{i} = '%f';
      elseif all(contains(lower(hd{i}), 'latitude') & ~any(strcmpi(hd(i~=1:size(hd,2)), 'latitude')) & ~any(strcmpi(hd(i~=1:size(hd,2)), 'lat')))
        hd{i} = 'lat';
        pars{i} = '%f';
      elseif all(contains(lower(hd{i}), 'longitude') & ~any(strcmpi(hd(i~=1:size(hd,2)), 'longitude')) & ~any(strcmpi(hd(i~=1:size(hd,2)), 'lon')))
        hd{i} = 'lon';
        pars{i} = '%f';
      elseif ~isnan(str2double(first_row{i}))
        pars{i} = '%f';
      else
        pars{i} = '%s';
      end
      hd{i} = strrep(strrep(hd{i}, '(', ''), ')', '');
    end
    % Set parser
    parser = strjoin(pars, '');
    
    % Read data
    t = textscan(fid, parser, 'delimiter',',');
    % Close file
    fclose(fid);
    
    data = table();
    % extract all variables
    sz_dt = size(t{strcmp(hd, 'dt')},1);
    for i = 1:size(hd, 2)
      if sz_dt ~= size(t{i},1) && size(data,1) > 0
        warning('Missing data on last %i row of file: deleted\n', sz_dt-size(t{i},1))
        data = data(1:size(t{i}, 1), :);
      end
      if contains(hd{i}, 'dt') && any(cell2mat(cellfun(@(c) contains(c, 'T'), t{i}, 'un', 0))) && ...
          any(cell2mat(cellfun(@(c) contains(c, 'Z'), t{i}, 'un', 0)))
        data.(hd{i}) = datetime(t{i}, 'inputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
      elseif contains(hd{i}, 'dt_instrument') && any(cell2mat(cellfun(@(c) contains(c, '/'), t{i}, 'un', 0))) && ...
          any(cell2mat(cellfun(@(c) contains(c, ':'), t{i}, 'un', 0)))
        data.(hd{i}) = datetime(t{i}, 'InputFormat', 'yyyy/MM/dd HH:mm:ss''+00:00');
      elseif contains(hd{i}, 'dt') && any(cell2mat(cellfun(@(c) contains(c, '/'), t{i}, 'un', 0))) && ...
          any(cell2mat(cellfun(@(c) contains(c, ':'), t{i}, 'un', 0)))
        data.(hd{i}) = datetime(t{i}, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
      elseif contains(hd{i}, 'dt_instrument') && any(cell2mat(cellfun(@(c) contains(c, '-'), t{i}, 'un', 0))) && ...
          any(cell2mat(cellfun(@(c) contains(c, ':'), t{i}, 'un', 0)))
        data.(hd{i}) = datetime(t{i}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss''+00:00');
      elseif contains(hd{i}, 'dt') && any(cell2mat(cellfun(@(c) contains(c, '-'), t{i}, 'un', 0))) && ...
          any(cell2mat(cellfun(@(c) contains(c, ':'), t{i}, 'un', 0)))
        data.(hd{i}) = datetime(t{i}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
      elseif contains(hd{i}, 'dt') && ~all(strcmp(t{i}, 'nan'))
        data.(hd{i}) = datetime(cellfun(@(x) [dt_ref x], t{i}, 'UniformOutput', false), 'InputFormat', 'yyyyMMddHH:mm:ss.SSS');
      elseif contains(hd{i}, 'dt')
        data.(hd{i}) = NaT(size(t{i}));
      elseif isnumeric(t{i})
        data.(hd{i}) = t{i};
      elseif all(strcmpi(t{i}, 'nan')) || any(~isnan(str2double(t{i})))
        data.(hd{i}) = str2double(t{i});
      elseif all(contains(lower(t{i}), {'true','false'})) % contains(hd{i}, 'swt')
        t{i}(contains(lower(t{i}), 'true')) = {'1'};
        t{i}(contains(lower(t{i}), 'false')) = {'0'};
        data.(hd{i}) = logical(str2num(cell2mat(t{i})));
      else
        data.(hd{i}) = t{i};
      end
    end
    % remove empty rows
    % data(all(isnat(data.dt) & isnat(data.dt_instrument) & isnan(data.lat) & isnan(data.lon), 2), :) = [];
    data(all(isnat(data.dt_instrument) & isnan(data.lat) & isnan(data.lon), 2), :) = [];

    % add units
    data.Properties.VariableUnits = unit;
  end
  
  % Remove last line if it's past midnight (bug in old Inlinino)
  if ~isempty(data) && size(data,1) > 1
    if data.dt(end-1) > data.dt(end)
      data(end,:) = [];
    end
  end
  
  if verbose; fprintf('Done\n'); end
end

function data = parseRMC(t)
  ndel = count(t, "*")+count(t, ",");
  if any(ndel ~= median(ndel)); t(ndel ~= median(ndel)) = []; end
  foo = split(t, {',', '*'});
  if size(foo, 2) == 15
    data = cell2table(cellstr(foo), 'VariableNames', ...
      {'dt','msg_ID','time','status','lat','lat_SN','lon','lon_WE','sog','cog',...
      'date','magnetic_var','magnetic_var_WE','nav_status','check_sum'});
  elseif size(foo, 2) == 14
    data = cell2table(cellstr(foo), 'VariableNames', ...
      {'dt','msg_ID','time','status','lat','lat_SN','lon','lon_WE','sog','cog',...
      'date','magnetic_var','magnetic_var_WE','check_sum'});
  else
    error('wrong number of variable names, check importInlininoNMEA - parseRMC')
  end
  data(cellfun('isempty', data.lat) | cellfun('isempty', data.lon), :) = [];
  data.dt = datetime(data.dt, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
  if contains(data.time, '.')
    data = addvars(data, datetime(strcat(data.date, data.time), 'InputFormat', 'ddMMyyHHmmss.SS'), ...
      'NewVariableNames', 'dt_instrument', 'After','dt');
  else
    data = addvars(data, datetime(strcat(data.date, data.time), 'InputFormat', 'ddMMyyHHmmss'), ...
      'NewVariableNames', 'dt_instrument', 'After','dt');
  end
  data.lat = str2double(cellfun(@(c) c(1:2), data.lat, 'un', 0)) + str2double(cellfun(@(c) c(3:end), data.lat, 'un', 0))/100; % should probably be /60 not /100 if format ddmm.mmmm
  data.lat(strcmp(data.lat_SN, 'S')) = -data.lat(strcmp(data.lat_SN, 'S'));
  data.lon = str2double(cellfun(@(c) c(1:3), data.lon, 'un', 0)) + str2double(cellfun(@(c) c(4:end), data.lon, 'un', 0))/100; % should probably be /60 not /100 if format ddmm.mmmm
  data.lon(strcmp(data.lon_WE, 'W')) = -data.lon(strcmp(data.lon_WE, 'W'));
  data.sog = str2double(data.sog);
  data.cog = str2double(data.cog);
  data = addvars(data, deg2rad(data.cog), 'NewVariableNames', 'magnetic_cog', 'After','cog');
  magnetic_var = deg2rad(str2double(data.magnetic_var));
  data.magnetic_cog(strcmp(data.magnetic_var_WE, 'E')) = data.magnetic_cog(strcmp(data.magnetic_var_WE, 'E')) - magnetic_var(strcmp(data.magnetic_var_WE, 'E'));
  data.magnetic_cog(strcmp(data.magnetic_var_WE, 'W')) = data.magnetic_cog(strcmp(data.magnetic_var_WE, 'W')) + magnetic_var(strcmp(data.magnetic_var_WE, 'W'));
  data.magnetic_cog = rad2deg(data.magnetic_cog);
  data = removevars(data, {'msg_ID','status','lat_SN','lon_WE','date','time',...
    'magnetic_var' ,'magnetic_var_WE','check_sum'});
  if any(strcmp(data.Properties.VariableNames, 'nav_status'))
    data = removevars(data, {'nav_status'});
  end
end

function data = parseGGA(t)
  ndel = count(t, "*")+count(t, ",");
  if any(ndel ~= median(ndel)); t(ndel ~= median(ndel)) = []; end
  foo = split(t, {',', '*'});
  if size(foo, 2) == 17
    data = cell2table(cellstr(foo), 'VariableNames', ...
      {'dt','msg_ID','time','lat','lat_SN','lon','lon_WE','nb_sat','h_precision','whateveritis', ...
      'altitude','altitude_unit','geoid_dist','geoid_dist_unit','DGPS_age','id_DGPS','check_sum'});
  else
    error('wrong number of variable names, check importInlininoNMEA - parseGGA')
  end
  data(cellfun('isempty', data.lat) | cellfun('isempty', data.lon), :) = [];
  data.dt = datetime(data.dt, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
  data = addvars(data, datetime(strcat(string(dateshift(data.dt,'start','day')), data.time), 'InputFormat', 'dd-MMM-yyyyHHmmss.SSS'), ...
    'NewVariableNames', 'dt_instrument', 'After','dt');
  data.lat = str2double(cellfun(@(c) c(1:2), data.lat, 'un', 0)) + str2double(cellfun(@(c) c(3:end), data.lat, 'un', 0))/60;
  data.lat(strcmp(data.lat_SN, 'S')) = -data.lat(strcmp(data.lat_SN, 'S'));
  data.lon = str2double(cellfun(@(c) c(1:3), data.lon, 'un', 0)) + str2double(cellfun(@(c) c(4:end), data.lon, 'un', 0))/60;
  data.lon(strcmp(data.lon_WE, 'W')) = -data.lon(strcmp(data.lon_WE, 'W'));
  data = removevars(data, {'msg_ID','time','lat_SN','lon_WE','nb_sat','h_precision',...
    'whateveritis','altitude','altitude_unit','geoid_dist','geoid_dist_unit','DGPS_age','id_DGPS','check_sum'});

end

function data = parseMWD(t)
  ndel = count(t, "*")+count(t, ",");
  if any(ndel ~= median(ndel)); t(ndel ~= median(ndel)) = []; end
  data = cell2table(cellstr(split(t, {',', '*'})), 'VariableNames', ...
    {'dt','msg_ID','wind_dir1','ref1','wind_dir2','ref2','wind_spd1','wind_spd1_unit','wind_spd2','wind_spd2_unit','check_sum'});
  data.dt = datetime(data.dt, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
  data.wind_dir1 = str2double(data.wind_dir1);
  data.wind_dir2 = str2double(data.wind_dir2);
  data.wind_spd1 = str2double(data.wind_spd1);
  data.wind_spd2 = str2double(data.wind_spd2);
  data.true_wind_dir = NaN(size(data.dt));
  data.true_wind_dir(strcmp(data.ref1, 'T')) = data.wind_dir1(strcmp(data.ref1, 'T'));
  data.true_wind_dir(strcmp(data.ref2, 'T')) = data.wind_dir2(strcmp(data.ref2, 'T'));
  data.true_wind_spd = NaN(size(data.dt));
  data.true_wind_spd(strcmp(data.wind_spd1_unit, 'N')) = data.wind_spd1(strcmp(data.wind_spd1_unit, 'N'));
  data.true_wind_spd(strcmp(data.wind_spd2_unit, 'N')) = data.wind_spd2(strcmp(data.wind_spd2_unit, 'N'));
  data = removevars(data, {'msg_ID','wind_dir1','wind_dir2','wind_spd1',...
    'wind_spd1_unit','wind_spd2','wind_spd2_unit','ref1','ref2','check_sum'});
end

function data = parseMWV(t)
  ndel = count(t, "*")+count(t, ",");
  if any(ndel ~= median(ndel)); t(ndel ~= median(ndel)) = []; end
  data = cell2table(cellstr(split(t, {',', '*'})), 'VariableNames', ...
    {'dt','msg_ID','wind_angle','ref','wind_spd','wind_spd_unit','status','check_sum'});
  data.dt = datetime(data.dt, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
  data.wind_angle = str2double(data.wind_angle);
  data.wind_spd = str2double(data.wind_spd);
  % check if ref == T, get true wind data
  data.true_wind_dir = NaN(size(data.dt));
  data.true_wind_spd = NaN(size(data.dt));
  data.true_wind_dir(strcmp(data.ref, 'T')) = data.wind_angle(strcmp(data.ref, 'T'));
  data.true_wind_spd(strcmp(data.ref, 'T')) = data.wind_spd(strcmp(data.ref, 'T'));
  % check if ref == R, get apparent wind data
  data.apparent_wind_dir = NaN(size(data.dt));
  data.apparent_wind_spd = NaN(size(data.dt));
  data.apparent_wind_dir(strcmp(data.ref, 'R')) = data.wind_angle(strcmp(data.ref, 'R'));
  data.apparent_wind_spd(strcmp(data.ref, 'R')) = data.wind_spd(strcmp(data.ref, 'R'));
  data = removevars(data, {'msg_ID','wind_angle','ref','wind_spd','wind_spd_unit','status','check_sum'});
end

function data = parseMMB(t)
  ndel = count(t, "*")+count(t, ",");
  if any(ndel ~= median(ndel)); t(ndel ~= median(ndel)) = []; end
  data = cell2table(cellstr(split(t, {',', '*'})), 'VariableNames', ...
    {'dt','msg_ID','atm_press1','atm_press_unit1','atm_press2','atm_press_unit2','check_sum'});
  data.dt = datetime(data.dt, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
  data.atm_press1 = str2double(data.atm_press1);
  data.atm_press2 = str2double(data.atm_press2);
  % get atm pressure in bar
  data.atm_press = NaN(size(data.dt));
  data.atm_press(strcmp(data.atm_press_unit1, 'B')) = data.atm_press1(strcmp(data.atm_press_unit1, 'B'));
  data.atm_press(strcmp(data.atm_press_unit2, 'B')) = data.atm_press2(strcmp(data.atm_press_unit2, 'B'));
  data = removevars(data, {'msg_ID','atm_press1','atm_press_unit1','atm_press2','atm_press_unit2','check_sum'});
end

function data = parseVPW(t)
  ndel = count(t, "*")+count(t, ",");
  if any(ndel ~= median(ndel)); t(ndel ~= median(ndel)) = []; end
  data = cell2table(cellstr(split(t, {',', '*'})), 'VariableNames', ...
    {'dt','msg_ID','spd_wind_relative1','spd_unit1','spd_wind_relative2','spd_unit2','check_sum'});
  data.dt = datetime(data.dt, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
  data.spd_wind_relative1 = str2double(data.spd_wind_relative1);
  data.spd_wind_relative2 = str2double(data.spd_wind_relative2);
  % get speed - measured parallel to wind in knots
  data.spd_wind_relative = NaN(size(data.dt));
  data.spd_wind_relative(strcmp(data.spd_unit1, 'N')) = data.spd_wind_relative1(strcmp(data.spd_unit1, 'N'));
  data.spd_wind_relative(strcmp(data.spd_unit2, 'N')) = data.spd_wind_relative2(strcmp(data.spd_unit2, 'N'));
  data = removevars(data, {'msg_ID','spd_wind_relative1','spd_unit1','spd_wind_relative2','spd_unit2','check_sum'});
end

function data = parseVHW(t)
  ndel = count(t, "*")+count(t, ",");
  if any(ndel ~= median(ndel)); t(ndel ~= median(ndel)) = []; end
  data = cell2table(cellstr(split(t, {',', '*'})), 'VariableNames', ...
    {'dt','msg_ID','heading_true','ref_T','heading_magnetic','ref_M','spd_water_relative1',...
    'spd_unit1','spd_water_relative2','spd_unit2','check_sum'});
  data.dt = datetime(data.dt, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
  data.heading_true = str2double(data.heading_true);
  data.heading_magnetic = str2double(data.heading_magnetic);
  data.spd_water_relative1 = str2double(data.spd_water_relative1);
  data.spd_water_relative2 = str2double(data.spd_water_relative2);
  % get speed - measured parallel to wind in knots
  data.spd_water_relative = NaN(size(data.dt));
  data.spd_water_relative(strcmp(data.spd_unit1, 'N')) = data.spd_water_relative1(strcmp(data.spd_unit1, 'N'));
  data.spd_water_relative(strcmp(data.spd_unit2, 'N')) = data.spd_water_relative2(strcmp(data.spd_unit2, 'N'));
  data = removevars(data, {'msg_ID','ref_T','ref_M','spd_water_relative1','spd_unit1','spd_water_relative2','spd_unit2','check_sum'});
end

function data = parseDBT(t)
  ndel = count(t, "*")+count(t, ",");
  if any(ndel ~= median(ndel)); t(ndel ~= median(ndel)) = []; end
  data = cell2table(cellstr(split(t, {',', '*'})), 'VariableNames', ...
    {'dt','msg_ID','depth_f','unit_f','depth_m','unit_m','depth_F','unit_F','check_sum'});
  data.dt = datetime(data.dt, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
  data.depth_m = str2double(data.depth_m);
  data = removevars(data, {'msg_ID','depth_f','unit_f','unit_m','depth_F','unit_F','check_sum'});
end

function data = parseMTW(t)
  ndel = count(t, "*")+count(t, ",");
  if any(ndel ~= median(ndel)); t(ndel ~= median(ndel)) = []; end
  data = cell2table(cellstr(split(t, {',', '*'})), 'VariableNames', ...
    {'dt','msg_ID','water_temp','unit','check_sum'});
  data.dt = datetime(data.dt, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
  data.water_temp = str2double(data.water_temp);
  data = removevars(data, {'msg_ID','unit','check_sum'});
end

function data = parseMTA(t)
  ndel = count(t, "*")+count(t, ",");
  if any(ndel ~= median(ndel)); t(ndel ~= median(ndel)) = []; end
  data = cell2table(cellstr(split(t, {',', '*'})), 'VariableNames', ...
    {'dt','msg_ID','air_temp','unit','check_sum'});
  data.dt = datetime(data.dt, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
  data.air_temp = str2double(data.air_temp);
  data = removevars(data, {'msg_ID','unit','check_sum'});
end

function data = parseVTG(t)
  ndel = count(t, "*")+count(t, ",");
  if any(ndel ~= median(ndel)); t(ndel ~= median(ndel)) = []; end
  data = cell2table(cellstr(split(t, {',', '*'})), 'VariableNames', ...
    {'dt','msg_ID','cog','unit_cog','magnetic_cog','unit_magnetic_cog',...
    'sog_kt','sog_kt_unit','sog_kmh','sog_kmh_unit','mode','check_sum'});
  data.dt = datetime(data.dt, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
  data.cog = str2double(data.cog);
  data.magnetic_cog = str2double(data.magnetic_cog);
  data.sog_kt = str2double(data.sog_kt);
  data.sog_kmh = str2double(data.sog_kmh);
  data = removevars(data, {'msg_ID','unit_cog','unit_magnetic_cog','sog_kt_unit','sog_kmh_unit','mode','check_sum'});
end
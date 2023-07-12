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
  if endsWith(filename, '.log')
    raw_format = true;
  elseif strcmp(hd{2}, 'packet')
    raw_format = true;
  else
    raw_format = false;
  end
  
  if raw_format
    t = fileread(filename);
    tab_sz = count(t,'RMC,');
    % preallocate table
    data = table(NaN(tab_sz, 1), NaN(tab_sz, 1), NaN(tab_sz, 1), NaN(tab_sz, 1), ...
      'VariableNames', {'dt', 'gps_dt', 'lat', 'lon'});
    % data = table(NaT(tab_sz, 1), NaT(tab_sz, 1), NaN(tab_sz, 1), NaN(tab_sz, 1), ...
    %   'VariableNames', {'dt', 'gps_dt', 'lat', 'lon'});
    for i = progress(1:size(data, 1))
      lin = fgetl(fid);
      if contains(lin, 'RMC,')
        [data.dt(i), data.gps_dt(i), data.lat(i), data.lon(i)] = parseGPRMC(lin);
      end
    end
    data(all(isnan(data.dt) & isnan(data.gps_dt) & isnan(data.lat) & isnan(data.lon), 2), :) = [];
    % data(all(isnat(data.dt) & isnat(data.gps_dt) & isnan(data.lat) & isnan(data.lon), 2), :) = [];
  else
    hd{strcmp(hd, 'time')} = 'dt';
    hd{strcmp(hd, 'datetime')} = 'gps_dt';
    hd{strcmp(hd, 'latitude')} = 'lat';
    hd{strcmp(hd, 'longitude')} = 'lon';
    % get units skipping empty lines (bug in old Inlinino)
    unit = fgetl(fid);
    while isempty(unit)
        unit = fgetl(fid);
    end
    % get units and lambda
    unit = strip(strsplit(unit, ','));
  
    % Set parser
    parser = ['%s%s' repmat('%f', 1, size(hd,2)-2)];
  
    % Read data
    t = textscan(fid, parser, 'delimiter',',');
    % Close file
    fclose(fid);
  
    % remove GPS time
    unit = unit(~strcmp(hd, 'gps_dt'));
    t = t(~strcmp(hd, 'gps_dt'));
    hd = hd(~strcmp(hd, 'gps_dt'));
  
    % Build table
    dat = [];
    for i = 1:size(hd, 2)
      if strcmp(hd{i}, 'dt')
        dat = [dat datenum(t{i}, 'yyyy/mm/dd HH:MM:SS.FFF')];
      elseif contains(hd{i}, 'swt')
        foo2 = t{i};
        foo2(contains(foo2, 'True')) = {'1'};
        foo2(contains(foo2, 'False')) = {'0'};
        dat = [dat logical(cell2mat(foo2))];
      else
        dat = [dat t{i}];
      end
    end
    data = array2table(dat, 'VariableNames', hd);
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

function [dt, gps_dt, lat, lon] = parseGPRMC(sentence)
  foo = strsplit(sentence, ',');
  if size(foo, 2) >= 11
    if contains(foo{1}, {'[',']'})
      foo_dt = strsplit(foo{1}, ' $');
      dt = datetime(foo_dt{1}, 'InputFormat', '[yyyy-MM-dd HH:mm:ss.SSS]');
      gps_dt = datetime([foo{10} foo{2}], 'InputFormat', 'ddMMyyHHmmss.SS');
      lat = str2double(foo{4}(1:2)) + str2double(foo{4}(3:end))/100;
      if strcmp(foo{5}, 'S'); lat = -lat; end
      lon = str2double(foo{6}(1:3)) + str2double(foo{6}(4:end))/100;
      if strcmp(foo{7}, 'W'); lon = -lon; end
    else
      if startsWith(foo{1}, '$')
        gps_dt = datetime([foo{10} foo{2}], 'InputFormat', 'ddMMyyHHmmss.SS');
        dt = gps_dt;
        lat = str2double(foo{4}(1:2)) + str2double(foo{4}(3:end))/100;
        if strcmp(foo{5}, 'S'); lat = -lat; end
        lon = str2double(foo{6}(1:3)) + str2double(foo{6}(4:end))/100;
        if strcmp(foo{7}, 'W'); lon = -lon; end
      else
        dt = datetime(foo{1}, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
        gps_dt = datetime([foo{11} foo{3}], 'InputFormat', 'ddMMyyHHmmss.SSS');
        lat = str2double(foo{5}(1:2)) + str2double(foo{5}(3:end))/100;
        if strcmp(foo{6}, 'S'); lat = -lat; end
        lon = str2double(foo{7}(1:3)) + str2double(foo{7}(4:end))/100;
        if strcmp(foo{8}, 'W'); lon = -lon; end
      end
    end
    dt = datenum(dt);
    gps_dt = datenum(gps_dt);
  else
    dt = NaN;
    gps_dt = NaN;
    lat = NaN;
    lon = NaN;
  end
end

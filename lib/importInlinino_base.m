function [ data, lambda ] = importInlinino_base( filename, verbose )
  % importInlininoNMEA Import data from basic Inlinino csv files
  % Author: Guillaume Bourdin
  % Date: June 30st, 2023
  %
  % Input:
  %   - filename: <char> filename including full path
  %   - verbose (optional)
  % 
  % Example: [ data, lambda] = importInlinino_base( filename, verbose )
  %%
  if nargin < 2; verbose = false; end
  if verbose
    foo = strsplit(filename, '/');
    fprintf('Importing %s ... ', foo{end});
  end
  
  % Get date from filename
  s = strsplit(filename, '_');
  dt_ref = s{end-1};
  
  % Open file
  fid=fopen(filename);
  if fid==-1
    error('Unable to open file: %s', filename);
  end
  
  % Get header
  hd = strip(strsplit(fgetl(fid), ','));
  if strcmp(filename(end-3:end), '.csv')
    % get units skipping empty lines (bug in old Inlinino)
    unit = fgetl(fid);
    while isempty(unit)
      unit = fgetl(fid);
    end
    % get units and lambda
    unit = strip(strsplit(strrep(unit, ',,', ', ,'), ','));
    first_row = fgetl(fid);
    if ~isnumeric(first_row) & first_row ~= -1
      first_row = strip(strsplit(first_row, ','));
      pars = repmat({''}, size(hd));
      % rename variables
      for i = 1:size(hd, 2)
        if strcmp(hd{i}, 'time')
          hd{i} = 'dt';
          pars{i} = '%s';
        elseif strcmp(hd{i}, 'datetime')
          hd{i} = 'dt_instrument';
          pars{i} = '%s';
        elseif strcmpi(hd{i}, 'switch') || contains(lower(hd{i}), 'pump')
          hd{i} = 'swt1';
          pars{i} = '%s';
        elseif strcmpi(hd{i}, 'switch(0)') || contains(lower(hd{i}), 'pump(0)')
          hd{i} = 'swt1';
          pars{i} = '%s';
        elseif strcmpi(hd{i}, 'switch(1)') || contains(lower(hd{i}), 'pump(1)')
          hd{i} = 'swt2';
          pars{i} = '%s';
        elseif strcmpi(hd{i}, 'flow(0)')
          hd{i} = 'spd1';
          pars{i} = '%f';
        elseif strcmpi(hd{i}, 'flow(1)')
          hd{i} = 'spd2';
          pars{i} = '%f';
        elseif strcmpi(hd{i}, 'flow(2)')
          hd{i} = 'spd3';
          pars{i} = '%f';
        elseif strcmpi(hd{i}, 'flow(3)')
          hd{i} = 'spd4';
          pars{i} = '%f';
        elseif ~isnan(str2double(first_row{i}))
          pars{i} = '%f';
        else
          pars{i} = '%s';
        end
        hd{i} = strrep(strrep(hd{i}, '(', ''), ')', '');
      end
      if size(unit,2) < size(hd,2)
        unit = [unit repmat({''},1,size(hd,2)-size(unit,2))];
      end
    else
      warning('Empty file: %s - ignored\n', filename)
      data = [];
      lambda = [];
      return
    end
  elseif strcmp(filename(end-3:end), '.raw')
    % hd = strip(strsplit(fgetl(fid), ','));
    % unit = cell(size(hd));
    unit = strip(strsplit(fgetl(fid), ','));
    first_row = fgetl(fid);
    if ~isnumeric(first_row) & first_row ~= -1
      first_row = strip(strsplit(first_row, ','));
      if size(first_row, 2) > size(hd,2)
        cols = first_row;
        nb_missing_hd = size(first_row, 2) - size(hd,2);
        hd = [hd strcat(repmat({'Var'},nb_missing_hd,1), cellstr(num2str((1:nb_missing_hd)')))'];
      else
        cols = hd;
      end
      pars = repmat({''}, size(cols));
      % rename variables
      for i = 1:size(cols, 2)
        if (contains(cols{i}, '/') && contains(cols{i}, ':')) || strcmpi(hd{i}, 'time')
          hd{i} = 'dt';
          pars{i} = '%s';
        elseif contains(cols{i}, '-') && contains(cols{i}, ':')
          hd{i} = 'dt_instrument';
          pars{i} = '%s';
        elseif contains(cols{i}, '-')
          hd{i} = 'date_instrument';
          pars{i} = '%s';
        elseif contains(cols{i}, '/')
          hd{i} = 'date_instrument';
          pars{i} = '%s';
        elseif contains(cols{i}, ':')
          hd{i} = 'time_instrument';
          pars{i} = '%s';
        elseif contains(lower(cols{i}), 'true') && contains(lower(cols{i}), 'false') || contains(lower(cols{i}), 'switch')
          hd{i} = 'swt';
          pars{i} = '%s';
        elseif any(strcmpi(cols, 'swt')) && i == 3 && ~isnan(str2double(first_row{i}))
          hd{i} = 'spd1';
          pars{i} = '%f';
          unit{i} = 'L/min';
        elseif any(strcmpi(cols, 'swt')) && i == 4 && ~isnan(str2double(first_row{i}))
          hd{i} = 'analog_gain';
          pars{i} = '%f';
        elseif any(strcmpi(cols, 'swt')) && i == 5 && ~isnan(str2double(first_row{i}))
          hd{i} = 'Analog2';
          pars{i} = '%f';
          unit{i} = 'V';
        elseif ~isnan(str2double(first_row{i}))
          if size(hd,2) < i
            hd{i} = ['Var' num2str(i)];
          end
          pars{i} = '%f';
        else
          if size(hd,2) < i
            hd{i} = ['Var' num2str(i)];
          end
          pars{i} = '%s';
        end
        hd{i} = strrep(strrep(hd{i}, '(', ''), ')', '');
      end
    else
      warning('Empty file: %s - ignored\n', filename)
      data = [];
      lambda = [];
      return
    end
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
    if contains(hd{i}, 'dt') && any(cell2mat(cellfun(@(c) contains(c, '/'), t{i}, 'un', 0))) && ...
        any(cell2mat(cellfun(@(c) contains(c, ':'), t{i}, 'un', 0)))
      data.(hd{i}) = datetime(t{i}, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
    elseif contains(hd{i}, 'dt') && any(cell2mat(cellfun(@(c) contains(c, '-'), t{i}, 'un', 0))) && ...
        any(cell2mat(cellfun(@(c) contains(c, ':'), t{i}, 'un', 0)))
      data.(hd{i}) = datetime(t{i}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    elseif contains(hd{i}, 'dt')
      data.(hd{i}) = datetime(cellfun(@(x) [dt_ref x], t{i}, 'UniformOutput', false), 'InputFormat', 'yyyyMMddHH:mm:ss.SSS');
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
  % reformat dt_instrument if date and time instruments are separated
  if any(strcmp(hd, 'date_instrument')) && any(strcmp(hd, 'time_instrument'))
    if contains(data.date_instrument, '-')
      data = addvars(data, datetime(strcat(data.date_instrument, repmat({' '}, size(data.dt)), ...
        data.time_instrument), 'InputFormat','dd-MM-yyyy HH:mm:ss'), 'NewVariableNames', 'dt_instrument', 'After', 'dt');
    elseif contains(data.date_instrument, '/')
      data = addvars(data, datetime(strcat(data.date_instrument, repmat({' '}, size(data.dt)), ...
        data.time_instrument), 'InputFormat','dd/MM/yyyy HH:mm:ss'), 'NewVariableNames', 'dt_instrument', 'After', 'dt');
    end
    data.date_instrument = [];
    data.time_instrument = [];
    unit(strcmp(hd,'date_instrument')) = strcat(unit(strcmp(hd,'date_instrument')), {' '}, unit(strcmp(hd,'time_instrument')));
    unit(strcmp(hd,'time_instrument')) = [];
    hd = strrep(hd, 'date_instrument','dt_instrument');
    hd(strcmp(hd,'time_instrument')) = [];
  end
  
  % detect if lambda to extract
  hd_tx = cell(size(hd));
  hd_digit = cell(size(hd));
  for i = 1:size(hd, 2)
    hd_tx{i} = hd{i}(isstrprop(hd{i},'alpha'));
    hd_digit{i} = hd{i}(isstrprop(hd{i},'digit'));
  end
  % [~, d] = unique(hd_tx, 'first');
  % var_lambda = hd_tx(not(ismember(1:numel(hd_tx), d)));
  [var_lambda, d] = unique(hd_tx, 'stable');
  if size(var_lambda,2) < size(hd_tx,2) & sum(str2double(hd_digit) > 300 & str2double(hd_digit) < 1200); merge_lambda_bool = true; else; merge_lambda_bool = false; end
  var_lambda = var_lambda(logical(sum(categorical(var_lambda) == ...
    {'beta','betap','bb','bbp','a','ap','c','cp','b','bp'}', 1)));
  % find wavelength and extract
  if merge_lambda_bool
    lambda = str2double(hd_digit(contains(hd, unique(var_lambda))));
  end
  
  % add units
  data.Properties.VariableUnits = unit;
  
  % merge into one column when spectral data
  if merge_lambda_bool
    merged_lambda = cat(2, table2array(data(:, contains(hd, unique(hd_tx(not(ismember(1:numel(hd_tx), d))))))));
    data(:, contains(hd, unique(hd_tx(not(ismember(1:numel(hd_tx), d)))))) = [];
    data.(cell2mat(unique(hd_tx(not(ismember(1:numel(hd_tx), d)))))) = merged_lambda;
    % data.Properties.VariableUnits(end) = unique(unit(not(ismember(1:numel(hd_tx), d))));
    data.Properties.VariableUnits = unit(d);
  end

  % Remove last line if it's past midnight (bug in old Inlinino)
  if ~isempty(data) && size(data,1) > 1
    if data.dt(end-1) > data.dt(end)
      data(end,:) = [];
    end
  end
  
  if verbose; fprintf('Done\n'); end
end


% % Build table starting by datetime
% id_dt = find(strcmp(hd, 'dt'));
% for i = 1:size(id_dt, 2)
% % if size(id_dt, 2) >= 1
%   if contains(t{id_dt(i)}(1), '/')
%     dt = datetime(t{id_dt}, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
%     % dt = datenum(t{i}, 'yyyy/mm/dd HH:MM:SS.FFF');
%   else
%     dt = datetime(cellfun(@(x) [dt_ref x], t{id_dt}, 'UniformOutput', false), 'InputFormat', 'yyyyMMddHH:mm:ss.SSS');
%     % dt = datenum(cellfun(@(x) [dt_ref x], t{1}, 'UniformOutput', false), 'yyyymmddHH:MM:SS.FFF');
%   end
% else
%   error('No datetime variable found in file')
% end
% % extract all other variables
% t = t(~strcmp(hd, 'dt'));
% hd = hd(~strcmp(hd, 'dt'));
% dat = [];
% for i = 1:size(hd, 2)
%   if size(dat,1) ~= size(t{i},1) && size(dat,1) > 0
%     warning('Missing data on last %i row of file: deleted\n', size(dat,1)-size(t{i},1))
%     dt = dt(1:size(t{i}, 1), :);
%     dat = dat(1:size(t{i}, 1), :);
%   end
%   if contains(hd{i}, 'swt')
%     foo2 = t{i};
%     foo2(contains(foo2, 'True')) = {'1'};
%     foo2(contains(foo2, 'False')) = {'0'};
%     dat = [dat logical(str2num(cell2mat(foo2)))];
%   else
%     dat = [dat t{i}];
%   end
% end
% 
% % detect if lambda to extract
% hd_tx = cell(size(hd));
% hd_digit = cell(size(hd));
% for i = 1:size(hd, 2)
%   hd_tx{i} = hd{i}(isstrprop(hd{i},'alpha'));
%   hd_digit{i} = hd{i}(isstrprop(hd{i},'digit'));
% end
% [~, d] = unique(hd_tx, 'first');
% var_lambda = hd_tx(not(ismember(1:numel(hd_tx), d)));
% var_lambda = var_lambda(logical(sum(categorical(var_lambda) == ...
%   {'beta','betap','bb','bbp','a','ap','c','cp'}', 1)));
% % find wavelength and extract
% if ~isempty(var_lambda)
%   lambda = str2double(hd_digit(contains(hd, unique())));
% end
% 
% data = array2table(dat, 'VariableNames', hd);
% data = addvars(data, dt, 'NewVariableNames', 'dt', 'Before', hd{1});
% data.Properties.VariableUnits = unit;
% 
% % merge into one column when spectral data
% if ~isempty(var_lambda)
%   merged_lambda = cat(2, dat(:, contains(hd, unique(hd_tx(not(ismember(1:numel(hd_tx), d)))))));
%   data(:, contains(hd, unique(hd_tx(not(ismember(1:numel(hd_tx), d)))))) = [];
%   data.(cell2mat(unique(hd_tx(not(ismember(1:numel(hd_tx), d)))))) = merged_lambda;
%   data.Properties.VariableUnits(end) = unique(unit(not(ismember(1:numel(hd_tx), d))));
% end
% 
% % Remove last line if it's past midnight (bug in old Inlinino)
% if ~isempty(data) && size(data,1) > 1
%   if data.dt(end-1) > data.dt(end)
%     data(end,:) = [];
%   end
% end
% 
% if verbose; fprintf('Done\n'); end

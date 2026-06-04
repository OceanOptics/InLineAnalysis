function cal_data = importLISST100XDeviceFile(filename)
  % import data from LISST100X device files
  fid=fopen(filename); % Open file
  if fid==-1
    error('Unable to open file: %s', filename);
  end
  % Get header to know how many columns
  foo = fgetl(fid);
  foo_split = strip(strsplit(foo, {',', ' ', '\t'}));
  fclose(fid); % Close file
  id_str = isnan(cell2mat(cellfun(@(c) str2double(c), foo_split, 'un', 0)));
  parser = repmat({''}, 1, size(foo_split, 2)); % Init parser
  parser(id_str) = {'%s'}; % Set string parser
  parser(~id_str) = {'%.f'}; % Set numeric parser
  parser = strjoin(parser, '');
  fid=fopen(filename); % reopen file
  if contains(foo, ',')
    t = textscan(fid, parser, 'Delimiter', ','); % Read data
  else
    t = textscan(fid, parser); % Read data
  end
  fclose(fid); % Close file
  if ~contains(parser, '%s')
    cal_data = cell2mat(t);
  else
    cal_data = t;
  end
  cal_data = cal_data(:)';
end
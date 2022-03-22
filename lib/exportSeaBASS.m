function exportSeaBASS(filename, meta, data, subfields)
% EXPORTSEABASS export data to SeaBASS format
[file_path,file_name,file_ext] = fileparts(filename);
if isempty(file_ext); file_ext = '.sb'; end

if isfile([file_path filesep file_name file_ext])
  delete([file_path filesep file_name file_ext])
end

% Create file
fid = fopen([file_path filesep file_name file_ext], 'W+', 'n', 'US-ASCII');

% Check fields
i_dt = strcmp(data.Properties.VariableNames, 'dt');
i_lat = strcmp(data.Properties.VariableNames, 'lat');
i_lon = strcmp(data.Properties.VariableNames, 'lon');
i_t = strcmp(data.Properties.VariableNames, 't');
i_s = strcmp(data.Properties.VariableNames, 's');
if sum(i_dt) ~= 1 || sum(i_lat) ~= 1 || sum(i_lon) ~= 1 || sum(i_t) ~= 1 || sum(i_s) ~= 1
  error('Missing one or more required field(s) is missing');
end
% Find fields to export
i_specific_fields = find(~(i_dt | i_lat | i_lon | i_t | i_s));
% Check subfileds (field with multiple columns)
if nargin < 4
%   subfields = {[]};
  subfields = repmat({[]}, 1, size(i_specific_fields,2));
else
  if size(i_specific_fields,2) ~= size(subfields,2)
    error('Missing subfields information');
  end
end
% Generate fieldnames
k = 1;
specific_fields = {}; specific_units = {};
for i=1:size(i_specific_fields,2)
  if isempty(subfields{i})
    specific_fields{k} = data.Properties.VariableNames{i_specific_fields(i)};
    specific_units{k} = data.Properties.VariableUnits{i_specific_fields(i)};
    k=k+1;
  elseif strfind(data.Properties.VariableNames{i_specific_fields(i)}, '_')
    foo = strsplit(data.Properties.VariableNames{i_specific_fields(i)}, '_');
    for j=1:size(subfields{i},2)
      specific_fields{k} = sprintf('%s%s%s', foo{1}, subfields{i}(j), sprintf('_%s', foo{2:end}));
      specific_units{k} = data.Properties.VariableUnits{i_specific_fields(i)};
      k=k+1;
    end
  else
    for j=1:size(subfields{i},2)
      specific_fields{k} = sprintf('%s%s', data.Properties.VariableNames{i_specific_fields(i)}, subfields{i}(j));
      specific_units{k} = data.Properties.VariableUnits{i_specific_fields(i)};
      k=k+1;
    end
  end
end

% Replace bad values
% Assume column 1-3 are dt, lat and lon
for i=4:size(data,2)
  % Replace Missing values
  data.(data.Properties.VariableNames{i})(isnan(table2array(data(:,i)))) = -9999;
  % Replace Out of range values
%   data.(data.Properties.VariableNames{i})(table2array(data(:,i)) < 0 & table2array(data(:,i)) ~= -9999) = -8888;
%   data.(data.Properties.VariableNames{i})(table2array(data(:,i)) == Inf) = -7777;
end

% Make header
fprintf(fid,'/begin_header\n');
fprintf(fid,'/investigators=%s\n', meta.investigators);
fprintf(fid,'/affiliations=%s\n', meta.affiliations);
fprintf(fid,'/contact=%s\n', meta.emails);
fprintf(fid,'/experiment=%s\n', meta.experiment);
fprintf(fid,'/cruise=%s\n', meta.cruise);
fprintf(fid,'/station=%s\n', meta.station);
fprintf(fid,'/data_file_name=%s\n', [file_name file_ext]);
fprintf(fid,'/documents=%s\n', meta.documents);
fprintf(fid,'/calibration_files=%s\n', meta.calibration_files);
fprintf(fid,'/data_type=%s\n',meta.data_type);
fprintf(fid,'/data_status=%s\n', meta.data_status);
% Date & Time
dt_min = min(data.dt);
dt_max = max(data.dt);
fprintf(fid,'/start_date=%s\n', datestr(dt_min,'yyyymmdd'));
fprintf(fid,'/end_date=%s\n', datestr(dt_max,'yyyymmdd'));
fprintf(fid,'/start_time=%s[GMT]\n', datestr(dt_min,'HH:MM:SS'));
fprintf(fid,'/end_time=%s[GMT]\n', datestr(dt_max,'HH:MM:SS'));
% Location
fprintf(fid,'/north_latitude=%.3f[DEG]\n', max(data.lat));
fprintf(fid,'/south_latitude=%.3f[DEG]\n', min(data.lat));
fprintf(fid,'/east_longitude=%.3f[DEG]\n', max(data.lon));
fprintf(fid,'/west_longitude=%.3f[DEG]\n', min(data.lon));
% More metadata
fprintf(fid,'/water_depth=NA\n');
fprintf(fid,'/measurement_depth=%.1f\n', meta.measurement_depth); % not allowed if depth is in the data
fprintf(fid,'/missing=-9999\n');
% fprintf(fid,'/below_detection_limit=-8888\n');
fprintf(fid,'/delimiter=comma\n'); % tab, space, or comma
% Fields
core_fields_sb = {'date', 'time', 'lat', 'lon', 'Wt', 'sal'};
core_units_sb = {'yyyymmdd', 'hh:mm:ss', 'degrees', 'degrees', 'degreesC', 'PSU'};
fields = [core_fields_sb(:)' specific_fields(:)'];
units = [core_units_sb(:)' specific_units(:)'];
foo = sprintf('%s,', fields{:}); fprintf(fid,'/fields=%s\n', foo(1:end-1));
foo = sprintf('%s,', units{:});fprintf(fid,'/units=%s\n', foo(1:end-1));
fprintf(fid,'/end_header\n');
fclose(fid);

dat = table(datestr(data.dt, 'yyyymmdd'), datestr(data.dt, 'HH:MM:SS'), round(data.lat, 4), ...
  round(data.lon, 4), round(data.t, 4), round(data.s, 4), 'VariableNames', {'Date', 'Time', ...
  'lat', 'lon', 't', 's'});
var_precision = data.Properties.VariableDescriptions(i_specific_fields);
int = find(contains(var_precision, {'%d','%i'}));
if any(int)
  for j = 1:size(int,2)
    var_precision{int(j)} = '%.0d';
  end
end
var_precision = str2num(cell2mat(regexprep(var_precision','\D','')));
for j=progress(1:size(i_specific_fields,2))
  var = data.Properties.VariableNames{i_specific_fields(j)};
  dat = [dat array2table(round(data.(var), var_precision(j)), 'VariableNames', ...
    cellfun(@(x) [var '_' x], cellstr(num2str((1:size(data.(var), 2))')), 'un', 0))];
end
writetable(dat,[file_path filesep file_name file_ext],'Encoding','US-ASCII', ...
  'Filetype', 'text','WriteMode','Append',...
    'WriteVariableNames',false,'WriteRowNames',false)
  
% % set separator
% com = repmat(',', size(data, 1), 1);
% % build date, time, lat, lon
% dat = [datestr(data.dt, 'yyyymmdd') com datestr(data.dt, 'HH:MM:SS') com char(num2str(round(data.lat, 4))), ...
%   com char(num2str(round(data.lon, 4))) com char(num2str(round(data.t, 4))) com char(num2str(round(data.s, 4)))];
% dat = cellstr(dat);
% % populate with instrument specific variables
% for j=progress(1:size(i_specific_fields,2))
%   dat = [dat join(compose('%g', (round(table2array(data(:,i_specific_fields(j))), 4))), ',', 2)];
% end
% dat = join(dat, ',');
% dat = [dat, repmat({'\n'}, size(data, 1), 1)];
% dat = char(join(dat,''));


% DEPRECATED (works well but very slow)
% Make content
% Pre-format dt
% date_str = datestr(data.dt, 'yyyymmdd');
% time_str = datestr(data.dt, 'HH:MM:SS');
% for i=progress(1:size(data.dt,1))
%   fprintf('%s,%s,%.4f,%.4f,%.4f,%.4f', date_str(i,:), time_str(i,:), data.lat(i), data.lon(i), data.t(i), data.s(i));
%   fprintf(fid,'%s,%s,%.4f,%.4f,%.4f,%.4f', date_str(i,:), time_str(i,:), data.lat(i), data.lon(i), data.t(i), data.s(i));
%   for j=1:size(i_specific_fields,2)
%     fprintf(fid,[',' data.Properties.VariableDescriptions{i_specific_fields(j)}], table2array(data(i,i_specific_fields(j))));
%   end
%   fprintf(fid,'\n');
% end
% fclose(fid);
end
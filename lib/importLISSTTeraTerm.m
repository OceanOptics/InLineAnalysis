function [ data ] = importLISSTTeraTerm( filename, verbose )
%IMPORTINLININO Import LISST data logged and Timestamped with TeraTerm during NAAMES 3

if nargin < 2; verbose = false; end
if verbose
  foo = strsplit(filename, '/');
  fprintf('Importing %s ... ', foo{end});
end

% Init data arrays
num_line = 1800; % 3.4 Mb files is 89745, one set of data is 49
dt = NaN(num_line, 1);
beta = NaN(num_line, 32);
laser_transmission = NaN(num_line, 1);
battery = NaN(num_line,1);
external_instrument = NaN(num_line,1);
laser_reference = NaN(num_line,1);
depth = NaN(num_line,1);
temperature = NaN(num_line,1);
lisst_dt = NaN(num_line,1);

% Open file
fid=fopen(filename);
if fid==-1
  error('Unable to open file: %s', filename);
end

% Read file line by line
flag_sample = false;
i=1; % sample number
j=1; % ring value | aux parameter index
while ~feof(fid)
  % Get line
  l = fgetl(fid);
  
  % Skip empty lines
  if isempty(l)
    continue
  end
  
  % Start sample
  if l(end) == '{'
    flag_sample = true;
    j = 1;
    % Check date format based on first caracter of date
    if ~isnan(str2double(l(2)))
      dt(i) = datenum(l(2:20),'yyyy-mm-dd HH:MM:SS');
    else
      dt(i) = datenum([l(26:29) l(6:20)], 'yyyymmm dd HH:MM:SS');
    end
    continue
  end
  
  % End sample
  if l(end) == '}'
    flag_sample = false;
    i=i+1;
    continue
  end
  
  % Read sample parameter by parameter
  if flag_sample
    foo = str2double(l(end-5:end));
    if j <= 32
      beta(i,j) = foo;
    elseif j == 33 % 
      laser_transmission(i) = foo; % mW % Transmission
    elseif j == 34
      battery(i) = foo * 0.01; % volts
    elseif j == 35
      external_instrument(i) = foo; % volts
    elseif j == 36
      laser_reference(i) = foo; % mw
    elseif j == 37
      depth(i) = foo * 0.01; % meters
    elseif j == 38
      temperature(i) = foo * 0.01; % deg C
    elseif j == 39
      lisst_dt(i) = floor(foo/100) + mod(foo,100)/24;
    elseif j == 40
      lisst_dt(i) = lisst_dt(i) + floor(foo/100)/24/60 + mod(foo,100)/24/3600;
    else
      fprintf('Ignoring extra parameters in sample %d', i);
    end
    j=j+1;
  end
end

fclose(fid);

% remove last sample if incomplete
%   need to be before table in case size of arrays are different
if flag_sample
  if size(dt, 1) >= i; dt(i) = []; end
  if size(beta, 1) >= i; beta(i,:) = []; end
  if size(laser_transmission, 1) >= i; laser_transmission(i) = []; end
  if size(battery, 1) >= i; battery(i) = []; end
  if size(external_instrument, 1) >= i; external_instrument(i) = []; end
  if size(laser_reference, 1) >= i; laser_reference(i) = []; end
  if size(depth, 1) >= i; depth(i) = []; end
  if size(temperature, 1) >= i; temperature(i) = []; end
  if size(lisst_dt, 1) >= i; lisst_dt(i) = []; end
%   end; end; end; end; end; end; end; end; end
end

% build table
data = table(dt, beta, laser_transmission, battery, external_instrument,...
             laser_reference, depth, temperature, lisst_dt);

% remove empty lines
data(isnan(data.dt), :) = [];

if verbose; fprintf('Done\n'); end

end



function [lambda_c, lambda_a, param] = importACSDeviceFile(filename)
% Get info out of ACS device file
param = struct();
if contains(filename,{'ac9','AC9'})
  lambda_c = [412 440 488 510 532 555 650 676 715];
  lambda_a = lambda_c;
  param.T_cal = [];
  param.I_cal = [];
else
  fid=fopen(filename);
  param.n_wl = NaN;
  lambda_c = [];
  lambda_a = [];
  param.n_tbin_cal = NaN;
  param.tbin_cal = [];
  i = 1;
  next_tbin = false;
  while ~feof(fid)
    l = fgetl(fid);
    if contains(lower(l), 'tcal')
      foo = strsplit(l, ' ');
      param.T_cal = str2double(foo{2});
      param.I_cal = str2double(foo{5});
    elseif contains(l, 'output wavelengths')
      % pre-allocate arrays with number of wavelength
      foo = strsplit(l, ';');
      param.n_wl = str2double(foo{1});
      lambda_c = NaN(1, param.n_wl);
      lambda_a = NaN(1, param.n_wl);
    elseif contains(l, 'number of temperature bins')
      foo = strsplit(l, ';');
      param.n_tbin_cal = str2double(foo{1});
      next_tbin = true;
    elseif next_tbin
      foo = strsplit(l, ';');
      foo0 = strsplit(foo{1}, '\t');
      param.tbin_cal = str2double(foo0(~cellfun('isempty', foo0)));
      next_tbin = false;
    elseif l(1) == 'C'
      % Get wavelength at each line
      foo = strsplit(l, '\t\t');
      bar = strsplit(foo{1}, '\t');
      lambda_c(i) = str2double(bar{1}(2:end));
      lambda_a(i) = str2double(bar{2}(2:end));
      i = i + 1;
    end
  end
  fclose(fid);
end

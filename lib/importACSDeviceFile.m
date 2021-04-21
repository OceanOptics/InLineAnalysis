function [lambda_c, lambda_a, T_cal, n] = importACSDeviceFile(filename)
% Get wavelength out of ACS device file
if contains(filename,{'ac9','AC9'})
  lambda_c = [412 440 488 510 532 555 650 676 715];
  lambda_a = lambda_c;
  T_cal = 25.50628;
else
  fid=fopen(filename);
  lambda_c = [];
  lambda_a = [];
  n = NaN; i = 1;
  while ~feof(fid)
    l = fgetl(fid);
    if contains(l, 'tcal')
      foo = strsplit(l, ' ');
      T_cal = str2double(foo{2});
    elseif contains(l, 'number of temperature bins')
      % Init arrays with number of wavelength
      foo = strsplit(l, ';');
      n = str2double(foo{1});
      lambda_c = NaN(1,n);
      lambda_a = NaN(1,n);
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

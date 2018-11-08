function [lambda_c, lambda_a, n] = importACSDeviceFile(filename, verbose)
% Get wavelength out of ACS device file

fid=fopen(filename);
lambda_c = [];
lambda_a = [];
n = NaN; i = 1;
while ~feof(fid)
    l = fgetl(fid);
    if contains(l, 'number of temperature bins')
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

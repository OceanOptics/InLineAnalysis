function [data_out, bad, Nbad] = AutoQC (instrument, data_in, lambda, tolerance, bb_dark, saturation_threshold, DI)
% author: Guillaume Bourdin
% created: Nov 13, 2019
%
% AutoQC: Auto QC ACs/AC9/BB spectra
%
% INPUT:
%   - data_in: raw data in <Nx3 Table> time series of data that must contains:
%         - dt <Nx1 datenum> date and time precise to the second
%         for each other field the folowwing field must exist:
%         - a spectra <NxM double> in column
%         - c spectra <NxM double> in column
%     or for BB data
%         - bb spectra <NxM double> in column
%   - lambda: structure with three fields:
%       - a: wavelenghts for absorption <1xM double> 
%       - c: wavelenghts for attenuation <1xM double>
%       - bb: wavelenghts of BB sensor <1xM double>
%   - tolerance: structure with three fields:
%       - a: threshold for automatic QC of absorption spectra (varies between ACS, must be >= 3)
%           (default = 3: Step = 3x > median(diff))
%       - c: threshold for automatic QC of attenuation spectra (varies between ACS, must be >= 3)
%           (default = 3: Step = 3x > median(diff))
%       - bb: threshold for automatic QC of beta spectra (must be >= 3)
%           (default = 3: diff(beta(lambda(i))) > 3 * median(beta(lambda(~i))
%   - bb_dark: beta dark measurements <1xM double>
%   - saturation_threshold: structure with three fields:
%       - a: wavelenghts for absorption <1xM double> 
%       - c: wavelenghts for attenuation <1xM double>
%       - bb: wavelenghts of BB sensor <1xM double>
% 
% threshold for automatic QC of beta spectra saturated (default = 4100)
%           (bb sensor saturation = 4130 counts)
%   - DI: boolean DI or not DI
%
% OUTPUT:
%   - data_out: cleaned data out <Nx3 Table> time series of data that must contains:
%         - dt <Nx1 datenum> date and time precise to the second
%         for each other field the folowwing field must exist:
%         - a spectra <NxM double> in column
%         - c spectra <NxM double> in column
%     or for BB data
%         - bb spectra <NxM double> in column
%   - bad: QC failed data table
%   - Nbad: structure with three fields:
%       - a: <double> percentage of absorption spectra deleted
%       - c: <double> percentage of attenuation spectra deleted
%       - bb: <1xM double> percentage of beta vales deleted

if nargin < 3
  error('input wavelenghts missing');
end
if nargin < 2
  error('input data in missing');
end
if nargin < 4
  tolerance.a = 3;
  tolerance.c = 3;
  warning('tolerance missing, set to default (3)');
end
if nargin < 7
  DI = false;
end

if contains(instrument, 'AC')
  if any(tolerance.a < 3 | tolerance.c < 3)
    warning('QC threshold might be too low, data might be lost');
  end
  if nargin < 6
    saturation_threshold.a = 50;
    warning('a threshold missing, set to default (50 m^{-1} uncalibrated)');
    saturation_threshold.c = 50;
    warning('c threshold missing, set to default (50 m^{-1} uncalibrated)');
  end
  
  data_out = data_in;
  Nbad.a = sum(any(data_in.a >= saturation_threshold.a, 2));
  Nbad.c = sum(any(data_in.c >= saturation_threshold.c, 2));
  bad = data_in(any(data_in.a >= saturation_threshold.a, 2) | any(data_in.c >= saturation_threshold.c, 2), :);
  data_out.a(any(data_in.a >= saturation_threshold.a, 2), :) = NaN;
  data_out.c(any(data_in.c >= saturation_threshold.c, 2), :) = NaN;

  % normalize data to homogenize auto QC
  datanorm = data_in;
  datanorm.a = datanorm.a - min(datanorm.a(:)) + 2;
  datanorm.c = datanorm.c - min(datanorm.c(:)) + 2;
  
  wl_a = lambda.a(lambda.a > 500 & lambda.a < 600);
  wl_c = lambda.c(lambda.c > 500 & lambda.c < 600);

  diff_a = [diff(datanorm.a(:, lambda.a > 500 & lambda.a < 600),[],2) NaN(size(datanorm,1),1)];
  diff_c = [diff(datanorm.c(:, lambda.c > 500 & lambda.c < 600),[],2) NaN(size(datanorm,1),1)];
  
  bad_a = max(abs(diff_a(:,wl_a > 560 & wl_a < 600)),[],2)...
    > tolerance.a*median(abs(diff_a(:,wl_a > 500 & wl_a < 550)),2);
  bad_c = max(abs(diff_c(:,wl_c > 560 & wl_c < 600)),[],2)...
    > tolerance.c*median(abs(diff_c(:,wl_c > 500 & wl_c < 550)),2);
  if DI 
    % segment database per DI event
    foo = find(diff(datetime(datanorm.dt, 'ConvertFrom', 'datenum')) > hours(0.5)) + 1;
    if isempty(foo)
      foo = [1 size(datanorm, 1)];
    else
      foo = [[1; foo], [foo-1; foo(end) + foo(1)]];
      foo(end, end) = size(datanorm, 1);
    end
    for i = 1:size(foo,1)
      seg = datanorm(foo(i, 1):foo(i, 2), :);
      % normalize each segments independently
      seg.a = seg.a - min(seg.a(:)) + 2;
      seg.c = seg.c - min(seg.c(:)) + 2;
      
      % remove spectra if any derivative over time between 500 and 660 nm is
      % 10x (a) and 5x (c) higher than average derivation over time AND any value between
      % 475 and 700 nm is greater than 99.7 percentile
      deriv = [NaN(1, size(seg.a,2)); diff(seg.a + min(seg.a(:)) + 500)];
      perc_a = prctile(seg.a, 99.7);
      segbad_a = false(size(seg.a, 1), 1);
      segbad_a(any(abs(deriv(:,lambda.a > 500 & lambda.a < 660)) > 10 * median(abs(deriv(:, lambda.a > 500 & lambda.a < 660)),1, 'omitnan'), 2) & ...
        any(seg.a(:,lambda.a > 500 & lambda.a < 660) > perc_a(lambda.a > 500 & lambda.a < 660), 2)) = true;

      deriv = [NaN(1, size(seg.c,2)); diff(seg.c + min(seg.c(:)) + 500)];
      perc_c = prctile(seg.c, 99.5);
      segbad_c = false(size(seg.c, 1), 1);
      segbad_c(any(abs(deriv(:,lambda.c > 500 & lambda.c < 660)) > 5 * median(abs(deriv(:, lambda.c > 500 & lambda.c < 660)),1, 'omitnan'), 2) & ...
        any(seg.c(:,lambda.c > 475 & lambda.c < 700) > perc_c(lambda.c > 475 & lambda.c < 700), 2)) = true;
      
      bad_a(foo(i, 1):foo(i, 2)) = bad_a(foo(i, 1):foo(i, 2)) + segbad_a;
      bad_c(foo(i, 1):foo(i, 2)) = bad_c(foo(i, 1):foo(i, 2)) + segbad_c;
    end
    bad_a = bad_a > 0;
    bad_c = bad_c > 0;
  end
  Nbad.a = (Nbad.a + sum(bad_a)) / size(data_in,1)*100;
  Nbad.c = (Nbad.c + sum(bad_c)) / size(data_in,1)*100;
  bad = [bad; data_in(bad_a | bad_c, :)];
  data_out.a(bad_a,:) = NaN;
  data_out.c(bad_c,:) = NaN;
  data_out(bad_a & bad_c,:) = [];
    
elseif contains(instrument, 'BB3')
  if nargin < 6
    saturation_threshold.bb = 4100;
    warning('beta threshold missing, set to default (4100 counts)');
  end

  if nargin < 5; error('beta dark missing (4th argument)'); end
  data_out = data_in;

  Nbad.bb = NaN(1, size(lambda.bb,2));
  for ii = 1:size(lambda.bb,2)
    Nbad.bb(ii) = sum(data_in.beta(:,ii) >= saturation_threshold.bb);
    data_out.beta(data_in.beta(:,ii) >= saturation_threshold.bb, ii) = NaN;
  end

  if any(tolerance.bb < 3); warning('QC threshold might be too low, data might be lost'); end
  bad_bb = false(size(data_in,1), size(lambda.bb, 2));
  for ii = 1:size(lambda.bb, 2)
    other = 1:size(lambda.bb, 2); other(ii)=[];
    bad_bb (:,ii) = data_in.beta(:,ii) - bb_dark(ii) > tolerance.bb * ...
      (median(data_in.beta(:,other),2, 'omitnan') - median(bb_dark(other),2, 'omitnan'));
%     Nbad.bb(ii) = Nbad.bb(ii) + sum(bad_bb(:,ii));
  end
  
  if DI 
    % segment database per DI event
    foo = find(diff(datetime(data_in.dt, 'ConvertFrom', 'datenum')) > hours(0.5)) + 1;
    if isempty(foo)
      foo = [1 size(data_in, 1)];
    else
      foo = [[1; foo], [foo-1; foo(end) + foo(1)]];
      foo(end, end) = size(data_in, 1);
    end
    for i = 1:size(foo,1)
      seg = data_in(foo(i, 1):foo(i, 2), :);
      % normalize each segments independently
      seg.beta = seg.beta - min(seg.beta(:)) + 2;
      
      % remove spectra if any derivative over time is
      % 5 times higher than average derivative over time AND
      % any value is greater than 99.7 percentile
      deriv = [NaN(1, size(seg.beta,2)); diff(seg.beta + min(seg.beta(:)) + 500)];
      perc = prctile(seg.beta, 95);
      segbad = false(size(seg.beta));
      segbad(abs(deriv) > tolerance.bb * median(abs(deriv),1, 'omitnan') & ...
        seg.beta > perc) = true;
      
      bad_bb(foo(i, 1):foo(i, 2), :) = bad_bb(foo(i, 1):foo(i, 2), :) + segbad;
    end
    bad_bb = bad_bb > 0;
  end
  
  Nbad.bb = sum(bad_bb) / size(data_in.beta,1) * 100;
  bad = data_in(bad_bb, :);
  data_out.beta(bad_bb) = NaN;
  data_out(all(isnan(data_out.beta), 2),:) = [];
  
elseif any(contains(instrument, {'HBB', 'HyperBB'}))
  % delete empty rows
  data_in(all(isnan(data_in.beta),2), :) = [];
  data_out = data_in;
  dup = data_in;
  dup.beta = dup.beta + median(dup.beta(:, end), 'omitnan');
  dup.beta = fillmissing(dup.beta, 'movmedian', 30);
  bad_bb = false(size(dup.beta));
  
  % detecting spike
  movdiff = movmedian(dup.beta, 50);
  bad_bb(dup.beta > movdiff + movdiff*tolerance.bb / 100) = true;
  
%   % detecting spike in dimenssion 1
%   diff_bb = [NaN(size(dup,1),1) diff(dup.beta,[],2)];
%   bad_bb = diff_bb > tolerance.bb * prctile(diff_bb, 95, 2);
%   
%   % detecting spike in dimenssion 2
%   diff_bb = [NaN(1, size(dup.beta, 2)); diff(dup.beta,[],1)];
%   movdiff = movmean(dup.beta, 10);
%   bad_bb(diff_bb > movdiff) = true;
%   
%   % old method
%   bad_bb(diff_bb > tolerance.bb * abs(median(movdiff, 1, 'omitnan'))) = true;
  
  data_out.beta(bad_bb) = NaN;
  
%   visProd3D(lambda.bb, dup.dt, dup.beta, ...
%     false, 'Wavelength', false, 71); zlabel('beta (m^{-1})'); %, 'Wavelength', true
%   visProd3D(lambda.bb, dup.dt, movdiff, ...
%     false, 'Wavelength', false, 72); zlabel('beta (m^{-1})'); %, 'Wavelength', true
%   visProd3D(lambda.bb, data.dt, diff_bb, ...
%     false, 'Wavelength', false, 73); zlabel('beta (m^{-1})'); %, 'Wavelength', true
%   visProd3D(lambda.bb, data.dt, abs(diff_bb), ...
%     false, 'Wavelength', false, 73); zlabel('beta (m^{-1})'); %, 'Wavelength', true
%   
%   visProd3D(lambda.bb, data.dt, data.beta, ...
%     false, 'Wavelength', false, 72); zlabel('beta (m^{-1})'); %, 'Wavelength', true

  if any(tolerance.bb < 3)
    warning('QC threshold might be too low, data might be lost');
  end
%   bad_bb = bad_bb_up + bad_bb_chl + bad_bb_down + bad_bb_time;
  bad_bb = bad_bb > 0;
  bad = data_in(bad_bb, :);
  data_out.beta(bad_bb) = NaN;
  data_out(sum(isnan(data_out.beta), 2) > size(data_out.beta, 2) / 3, :) = [];
  % count only bad that were not already NaN (interpolation)
  bad_bb(isnan(data_in.beta)) = false;
  Nbad.bb = sum(bad_bb) / size(data_in.beta,1) * 100;
end

fnames = fieldnames(Nbad);
for i = 1:size(fnames,1)
  if isempty(Nbad.(fnames{i}))
    Nbad.(fnames{i}) = 0;
  end
end

% visProd3D(lambda.bb, datanorm.dt, datanorm.beta, ...
%   false, 'Wavelength', false, 72); zlabel('beta (m^{-1})'); %, 'Wavelength', true
% visProd3D(lambda.bb, datanorm.dt, diff_bb, ...
%   false, 'Wavelength', false, 73); zlabel('beta (m^{-1})'); %, 'Wavelength', true
% visProd2D(lambda.bb, datanorm.dt(1), median(diff_bb), ...
%   false, 75); zlabel('beta (m^{-1})'); %, 'Wavelength', true


% visProd3D(lambda.bb, datanorm.dt, bad_bb_time, ...
%   false, 'Wavelength', false, 78); zlabel('beta (m^{-1})'); %, 'Wavelength', true


% visProd3D(lambda_ref, dataini.dt(100000:101000,:), dataini.a(100000:101000,:), false, 'Wavelength', false, 72); zlabel('a_p (m^{-1})'); %, 'Wavelength', true
% visProd3D(lambda_ref, dataini.dt(100000:101000,:), dataini.c(100000:101000,:), false, 'Wavelength', false, 73); zlabel('c_p (m^{-1})');
% visProd3D(lambda_ref, data.dt(100000:121000,:), data.a(100000:121000,:), false, 'Wavelength', false, 74); zlabel('a_p (m^{-1})'); %, 'Wavelength', true
% visProd3D(lambda_ref, data.dt(100000:121000,:), data.c(100000:121000,:), false, 'Wavelength', false, 75); zlabel('c_p (m^{-1})');

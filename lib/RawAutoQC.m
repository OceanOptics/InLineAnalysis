function [data, Nbad] = RawAutoQC (instrument, data, lambda, fudge_factor, bb_dark, bb_threshold, DI)
% author: Guillaume Bourdin
% created: Nov 13, 2019
%
% RawAutoQC: Auto QC ACs/AC9/BB spectrum and PAR data
%
% INPUT:
%   - data_in: ACS data <Nx3 Table> time series of data that must contains:
%         - dt <Nx1 datenum> date and time precise to the second
%         for each other field the folowwing field must exist:
%         - a spectrum <NxM double> in column
%         - c spectrum <NxM double> in column
%     or for BB data
%         - bb spectrum <NxM double> in column
%   - lambda: structure with three fields:
%       - a: wavelenghts for absorption <1xM double> 
%       - c: wavelenghts for attenuation <1xM double>
%       - bb: wavelenghts of BB sensor <1xM double>
%   - fudge_factor: structure with three fields:
%       - a: threshold for automatic QC of absorption spectrum (varies between ACS, must be >= 3)
%           (default = 3: Step = 3x > mean(diff))
%       - c: threshold for automatic QC of attenuation spectrum (varies between ACS, must be >= 3)
%           (default = 3: Step = 3x > mean(diff))
%       - bb: threshold for automatic QC of beta spectrum (must be >= 3)
%           (default = 3: Step = beta(lambda(i)) > 3 * mean(beta(lambda(~i))
%   - bb_dark: beta dark measurements <1xM double>
%   - bb_threshold: threshold for automatic QC of beta spectrum saturated (default = 4100)
%           (bb sensor saturation = 4130 counts)
%   - DI: boolean DI or not DI
%
% OUTPUT:
%   - data_in: quality filtered ACS data <Nx3 Table> time series of data that must contains:
%         - dt <Nx1 datenum> date and time precise to the second
%         for each other field the folowwing field must exist:
%         - a spectrum <NxM double> in column
%         - c spectrum <NxM double> in column
%     or for BB data
%         - bb spectrum <NxM double> in column
%   - Nbad: structure with three fields:
%       - a: <double> percentage of absorption spectrum deleted
%       - c: <double> percentage of attenuation spectrum deleted
%       - bb: <1xM double> percentage of beta vales deleted

if nargin < 3
  error('input wavelenghts missing');
end
if nargin < 2
  error('input data in missing');
end
if nargin < 4
  fudge_factor.a = 3;
  fudge_factor.c = 3;
  warning('fudge_factor missing, set to default (3)');
end
if nargin < 7
  DI = false;
end

if contains(instrument, 'AC')
  if any(fudge_factor.a < 3 | fudge_factor.c < 3)
    warning('QC threshold might be too low, data might be lost');
  end
  % normalize data to homogenize auto QC
  datanorm = data;
  datanorm.a = datanorm.a - min(datanorm.a(:)) + 2;
  datanorm.c = datanorm.c - min(datanorm.c(:)) + 2;
  
  wl_a = lambda.a(lambda.a > 500 & lambda.a < 600);
  wl_c = lambda.c(lambda.c > 500 & lambda.c < 600);

  diff_a = [diff(datanorm.a(:,lambda.a > 500 & lambda.a < 600),[],2) NaN(size(datanorm,1),1)];
  diff_c = [diff(datanorm.c(:,lambda.c > 500 & lambda.c < 600),[],2) NaN(size(datanorm,1),1)];

  bad_a = max(abs(diff_a(:,wl_a > 560 & wl_a < 600)),[],2)...
    > fudge_factor.a*mean(abs(diff_a(:,wl_a > 500 & wl_a < 550)),2);
  bad_c = max(abs(diff_c(:,wl_c > 560 & wl_c < 600)),[],2)...
    > fudge_factor.c*mean(abs(diff_c(:,wl_c > 500 & wl_c < 550)),2);
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
      
      % remove spectrum if any derivation over time between 500 and 660 nm is
      % 10x (a) and 5x (c) higher than average derivation over time AND any value between
      % 475 and 700 nm is greater than 99.7 percentile
      deriv = [NaN(1, size(seg.a,2)); diff(seg.a + min(seg.a(:)) + 500)];
      perc_a = prctile(seg.a, 99.7);
      segbad_a = false(size(seg.a, 1), 1);
      segbad_a(any(abs(deriv(:,lambda.a > 500 & lambda.a < 660)) > 10 * mean(abs(deriv(:, lambda.a > 500 & lambda.a < 660)),1, 'omitnan'), 2) & ...
        any(seg.a(:,lambda.a > 500 & lambda.a < 660) > perc_a(lambda.a > 500 & lambda.a < 660), 2)) = true;

      deriv = [NaN(1, size(seg.c,2)); diff(seg.c + min(seg.c(:)) + 500)];
      perc_c = prctile(seg.c, 99.5);
      segbad_c = false(size(seg.c, 1), 1);
      segbad_c(any(abs(deriv(:,lambda.c > 500 & lambda.c < 660)) > 5 * mean(abs(deriv(:, lambda.c > 500 & lambda.c < 660)),1, 'omitnan'), 2) & ...
        any(seg.c(:,lambda.c > 475 & lambda.c < 700) > perc_c(lambda.c > 475 & lambda.c < 700), 2)) = true;
      
      bad_a(foo(i, 1):foo(i, 2)) = bad_a(foo(i, 1):foo(i, 2)) + segbad_a;
      bad_c(foo(i, 1):foo(i, 2)) = bad_c(foo(i, 1):foo(i, 2)) + segbad_c;
    end
    bad_a = bad_a > 0;
    bad_c = bad_c > 0;
  end

  % Nbad_tot = sum(bad_a & bad_c)/size(data,1)*100;
  Nbad.a = sum(bad_a)/size(data,1)*100;%-Nbad_tot;
  Nbad.c = sum(bad_c)/size(data,1)*100;%-Nbad_tot;
  data.a(bad_a,:) = NaN;
  data.c(bad_c,:) = NaN;
  data(bad_a & bad_c,:)=[];
    
elseif contains(instrument, 'BB3')
  if nargin < 6
    bb_threshold = 4100;
    warning('beta threshold missing, set to default (4100 counts)');
  end

  if nargin < 5
    error('beta dark missing (4th argument)');
  end
  Nbad.bb = NaN(1, size(lambda.bb,2));
  for ii = 1:size(lambda.bb,2)
    Nbad.bb(ii) = sum(data.beta(:,ii) >= bb_threshold);
    data.beta(data.beta(:,ii) >= bb_threshold, ii) = NaN;
  end

  if any(fudge_factor.bb < 3)
    warning('QC threshold might be too low, data might be lost');
  end
  bad_bb = NaN(size(data,1), size(lambda.bb, 2));
  for ii = 1:size(lambda.bb, 2)
    other = 1:size(lambda.bb, 2); other(ii)=[];
    bad_bb (:,ii) = data.beta(:,ii) - bb_dark(ii) > fudge_factor.bb * ...
      (mean(data.beta(:,other),2, 'omitnan') - mean(bb_dark(other),2, 'omitnan'));
    Nbad.bb(ii) = Nbad.bb(ii) + sum(bad_bb(:,ii));
  end
  data.beta(logical(bad_bb)) = NaN;
  Nbad.bb = Nbad.bb / size(data.beta(:,ii),1) * 100;
  
elseif contains(instrument, 'HBB')

  % normalize data to homogenize auto QC
  datanorm = data;
  datanorm.beta = fillmissing(datanorm.beta, 'linear', 2);
  datanorm.beta = datanorm.beta - min(datanorm.beta(:)) + 2;
  
  diff_bb = [NaN(size(datanorm,1),1) diff(datanorm.beta,[],2)];
  diff_bb(:,1) = median(diff_bb,2, 'omitnan');
  diff_neg = diff_bb;
  diff_neg(diff_neg > 0) = NaN;
  bad_bb_up = diff_bb > abs(median(diff_neg, 2, 'omitnan'));
%   bad_bb = diff_bb > fudge_factor.bb * median(abs(diff_bb),2, 'omitnan');
  bad_bb_down = diff_bb < fudge_factor.bb * median(diff_neg,2, 'omitnan');

  if any(fudge_factor.bb < 3)
    warning('QC threshold might be too low, data might be lost');
  end
  bad_bb = bad_bb_up + bad_bb_down;
  bad_bb = bad_bb > 0;
  data.beta(bad_bb) = NaN;
  data(sum(isnan(data.beta), 2) > size(data.beta, 2) /3, :) = [];
  Nbad.bb = sum(bad_bb) / size(data.beta,1) * 100;
end

% visProd3D(lambda.bb, datanorm.dt, datanorm.beta, ...
%   false, 'Wavelength', false, 72); zlabel('beta (m^{-1})'); %, 'Wavelength', true


% visProd3D(lambda_ref, dataini.dt(100000:101000,:), dataini.a(100000:101000,:), false, 'Wavelength', false, 72); zlabel('a_p (m^{-1})'); %, 'Wavelength', true
% visProd3D(lambda_ref, dataini.dt(100000:101000,:), dataini.c(100000:101000,:), false, 'Wavelength', false, 73); zlabel('c_p (m^{-1})');
% visProd3D(lambda_ref, data.dt(100000:121000,:), data.a(100000:121000,:), false, 'Wavelength', false, 74); zlabel('a_p (m^{-1})'); %, 'Wavelength', true
% visProd3D(lambda_ref, data.dt(100000:121000,:), data.c(100000:121000,:), false, 'Wavelength', false, 75); zlabel('c_p (m^{-1})');

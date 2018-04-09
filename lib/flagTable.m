function [ dflag, stats ] = flagTable( dbin, params )
%CLEANTABLE clean data based on statistics of the dataset
%   
% INPUT:
%   - data_in <N Table> time series of data that must contains:
%         dt <Nx1 datenum> date and time precise to the second
%         for each other field the folowwing field must exist:
%             _p5 <NxM double>  5th percentile
%             _p95 <NxM double> 95th percentile
%             _std <NxM double> standard deviation
%
% OUTPUT:
%   - data_out <N Table> time series without all the flagged data
%   - flags <Nx1 Table> list of flags matching data_in timestamp
%       must contain fields:
%         dt <Nx1 datenum>
%         flag <Nx1 int32>
%           0: No QC
%           1: No data
%           2: Good
%           4: Questionnable Manual
%           8: Bad Manual
%           16: 2x > median(TimeSeries)
%           32: 95th percentile - 5th percentile > X counts

%TODO
% Check coefficient of variance in

% Set constants
if nargin < 2
  params.maximum_fudge_factor = 2;
  params.variance_fudge_factor = 2;
  params.abs_variance_coef = 0.005; % value for ACS
  params.rel_variance_coef = 0.05;
end

% Get dt field
idt = strcmp(dbin.Properties.VariableNames,'dt');
if ~any(idt) || sum(idt) > 1
  error('dt field cannot be find');
end
% Fields to exclude
in = strcmp(dbin.Properties.VariableNames,'n');
if ~any(idt)
  error('n field cannot be find, must run bin before cleanTable');
end
istats = ~cellfun(@isempty, strfind(dbin.Properties.VariableNames,'_'));
% Get other fields
lvar = dbin.Properties.VariableNames(~(idt|in|istats));

% Compute standard stats for time series
% assume most of the time series is good data
for j=1:size(lvar,2)
  % Maximum value threshol
  stats.(lvar{j}).med_median = nanmedian(dbin.(lvar{j}));
  stats.(lvar{j}).maximum_threshold = stats.(lvar{j}).med_median * params.maximum_fudge_factor;
  % Variance threshold
  stats.(lvar{j}).med_variance = nanmedian(dbin.([lvar{j} '_dtc_var']));
  stats.(lvar{j}).variance_threshold = stats.(lvar{j}).med_variance * params.variance_fudge_factor;
%   stats.(lvar{j}).uncertainty_threshold = params.abs_uncertainty + params.rel_uncertainty .* dbin.(lvar{j});
  stats.(lvar{j}).uncertainty_threshold = max(params.abs_uncertainty, params.rel_uncertainty .* dbin.(lvar{j}));
  flag_smooth_error = false;
  while params.smooth_threshold > 1
    try
      stats.(lvar{j}).smoothed_uncertainty_threshold = filtfilt(ones(params.smooth_threshold,1), params.smooth_threshold, stats.(lvar{j}).uncertainty_threshold);
      break;
    catch
      params.smooth_threshold = params.smooth_threshold - 1;
      flag_smooth_error = true;
    end
  end
  if flag_smooth_error; warning('smooth_threshold was recuiced to: %d', params.smooth_threshold); end
end

% Copy input to output
dflag = table();
% Flag data
for j=1:size(lvar,2)
  % Init flag for variable
  flags = zeros(size(dbin,1),1, 'uint32');
  % Test 16=2^4: Values too large
  flags = flags + uint32(16 * (dbin.(lvar{j}) > stats.(lvar{j}).maximum_threshold));
  % Test 32=2^5: Variance too large
  flags = flags + uint32(32 * (dbin.([lvar{j} '_dtc_var']) > stats.(lvar{j}).variance_threshold));
  % Test 64=2^6: Outisde of uncertainty bracket (run in binTable)
%   flags = flags +uint32(64);
  % Test 128=2^7: |mean - median| (Wendy's criteria)
  flags = flags + uint32(128 * (abs(dbin.([lvar{j} '_dtc_mn']) - dbin.([lvar{j} '_dtc_md'])) > 1 / params.avg_sensitivity .*...
          stats.(lvar{j}).smoothed_uncertainty_threshold));
  % Test 256=2^8: Uncertainty is too large (based on variance)
  flags = flags + uint32(256 * (dbin.([lvar{j} '_dtc_unc']) > 1 / params.unc1_sensitivity .*...
          stats.(lvar{j}).smoothed_uncertainty_threshold));
  % Test 512=2^9: Standard error is too large
  flags = flags + uint32(512 * (dbin.([lvar{j} '_dtc_se']) > 1 / params.unc2_sensitivity .*...
          stats.(lvar{j}).smoothed_uncertainty_threshold));
  
  % Test 2 =2^1: Good data
  flags = flags + uint32(2 * (flags == 0));
  % Write flag in data
  dflag.([lvar{j} '_flag']) = flags;
end

end

function read_flag_minimal(flag)
  % flag should be uint32()
  fprintf('%s\n', mat2str(bitget(flag,32:-1:1)));
end

function read_flag_verbose(flag)
  flag_name{1} = 'No Data';
  flag_name{2} = 'Good';
  flag_name{3} = 'Suspect Manual';
  flag_name{4} = 'Bad Manual';
  flag_name{5} = 'Value too large';
  flag_name{6} = 'Uncertainty too large';
  flag_name{7} = 'Outside percentile range';
  flag_name{8} = 'Wendy''s test';
  foo = boolean(bitget(flag,1:32));
  flag_bits = find(1 == foo);
  
  if isempty(flag_bits)
    fprintf('0\t No QC\n');
  else
    for bit=flag_bits
      fprintf('%2d\t%2d\t%s\n', bit, 2^(bit-1), flag_name{bit});
    end
  end
end
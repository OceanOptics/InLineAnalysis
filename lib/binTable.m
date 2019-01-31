function [ bin, gflag ] = binTable( raw, bin_size, mode, prctile_dtc, prctile_avg, dt_discontinus, parallel_flag, verbose )
%BINTABLE compute median,, 5 and 95 percentiles, standard deviation, and number
%  observation for each bin. The size of each bin is defined by bin_size
% bin_size is in seconds
% raw must contain a dt column
% parallel: 0 -> False | Inf: -> Max number of thread running
% mode: 
%     4flag: Method to use to flag automatically
%              prctile (low & high),
%              mean, median, standard and deviation withtin percentile interval, 
%              n, variance, uncertainty for median (1), and uncertainty for mean (2)
%     SB_IN_PRCTL:  Faster method outputting only requirements for SeaBASS (not compatible with automatic flagging)
%     SB_ALL:  Faster method outputting only requirements for SeaBASS (not compatible with automatic flagging)
% dt_discontinuous: recommended to use on sparse dataset (typically DI runs)
%                   false (default), much slower

%
% OUTPUT:
%   bin <Table> binned by minute with statics based on mode
%   gflag <Table> only available for 4flag mode
%             synchronise with raw.dt and NOT bin.dt

% Bench Mark (2.8 GHz Intel Core i7, OSX: 10.12.6, Matlab r2017a)
% >>tic; for i=1:10; binTable(bb3_tot, 60/3600/24, 'mnmd', [2.5 97.5], 8, true); end; toc
% Elapsed time is 42.037515 seconds.
% >> tic; for i=1:10; binTable(bb3_tot, 60/3600/24, 'prctile', [2.5 97.5], 8, true); end; toc
% Elapsed time is 42.298939 seconds.
% >> tic; for i=1:10; binTable(bb3_tot, 60/3600/24, 'full', [2.5 97.5], 8, true); end; toc
% Elapsed time is 40.625535 seconds.
% >> tic; for i=1:10; binTable(bb3_tot, 60/3600/24, 'full', [2.5 97.5], 4, true); end; toc
% Elapsed time is 41.553749 seconds.
% >> tic; for i=1:10; binTable(bb3_tot, 60/3600/24, 'full', [2.5 97.5], 0, true); end; toc
% Elapsed time is 66.181466 seconds.
% ==> Mode does not matter for speed
% => similar speed for 4 and 8 Threads

% Set default parameters
if nargin < 2; bin_size = 1/60/24; end
if nargin < 3; mode = '4flag'; end
if nargin < 4; prctile_dtc = [2.5 97.5]; end
if nargin < 5; prctile_avg = [5 75]; end
if nargin < 6; dt_discontinus = false; end
if nargin < 7; parallel_flag = 0; end
if nargin < 8; verbose = false; end

% Get dt field
idt = strcmp(raw.Properties.VariableNames,'dt');
if ~any(idt) || sum(idt) > 1
  error('dt field cannot be find');
end
% Get other fields
lvar = raw.Properties.VariableNames(~idt);


if verbose; fprintf('Binning %s ... ', inputname(1)); end

% Get start dt
start_dt = min(raw.dt);
start_dt = start_dt - second(start_dt)/3600/24; % Start at 0 second

% Init bin array with NaN values
if dt_discontinus  
  dt = [];
  i = 1; n = size(raw.dt, 1);
  while i <= n
    dti = raw.dt(i);
    sel = find(dti <= raw.dt & raw.dt < dti + bin_size);
    if isempty(sel)
      error('Issue with discontinuous time code');
    end
    dt(end+1,1) = mean(raw.dt(sel));
    i = sel(end) + 1;
  end
else
  dt(:,1) = start_dt:bin_size:max(raw.dt); % max range
end
bin = table(dt);
gflag = table(); % Used only in clean 1 mode;
% Compute stats for each bin
% It takes some time
switch mode
  case '4flag'
    % Compute upper and lower percentile (2.5 97.5 percentile by defaul)
    % Compute mean and median of data within percentile interval
    for j=1:size(lvar,2)
      % Init Var In
      raw_dt = raw.dt;
      raw_var = raw.(lvar{j});
      % Init var Out for detection (used for flagging data)
      dtc_pl = NaN(size(dt,1),size(raw_var,2));
      dtc_ph = NaN(size(dt,1),size(raw_var,2));
      dtc_md = NaN(size(dt,1),size(raw_var,2));
      dtc_mn = NaN(size(dt,1),size(raw_var,2));
      dtc_sd = NaN(size(dt,1),size(raw_var,2));
%       dtc_n  = NaN(size(dt,1),size(raw_var,2));
      dtc_n  = NaN(size(dt));
      % Init car out for averaging (used in rest of program)
      avg_pl = NaN(size(dt,1),size(raw_var,2));
      avg_ph = NaN(size(dt,1),size(raw_var,2));
      avg_md = NaN(size(dt,1),size(raw_var,2));
      avg_mn = NaN(size(dt,1),size(raw_var,2));
      avg_sd = NaN(size(dt,1),size(raw_var,2));
%       avg_n  = NaN(size(dt,1),size(raw_var,2));
      avg_n  = NaN(size(dt));
      parfor (i=1:size(dt,1), parallel_flag)
%       for i=1:size(dt,1)
        % Select minute to bin
        sel = dt(i) - bin_size/2 <= raw_dt & raw_dt < dt(i) + bin_size/2;
        % Unselect raw with NaN values
        foo = ~isnan(raw_var); % step required for 2d variables
        sel = sel & all(foo,2);
        % Check number of samples
        if sum(sel) > 0
          raw_var_sel = raw_var(sel,:);
          % Compute percentiles
          foo = prctile(raw_var_sel,[prctile_dtc, prctile_avg],1);
          dtc_pl(i,:) = foo(1,:); % low percentile
          dtc_ph(i,:) = foo(2,:); % high percentile
          avg_pl(i,:) = foo(3,:); % low percentile
          avg_ph(i,:) = foo(4,:); % high percentile
          
          % Select data within percentile interval
          %     ideal way, each column of the variable is independent
          %     does not run in parallel
          %     if one column bad, the entire bin will be flagged anyway
%           dtc_sel = dtc_pl(i,:) <= raw_var_sel & raw_var_sel <= dtc_ph(i,:);
%           avg_sel = avg_pl(i,:) <= raw_var_sel & raw_var_sel <= avg_ph(i,:);
%           dtc_n(i,:) = sum(dtc_sel);
%           avg_n(i,:) = sum(avg_sel);
%           for k=1:size(raw_var,2) % for each column in 2d array (might be a way to optimize that section)
%             dtc_md(i, k) = median(raw_var_sel(dtc_sel(:,k)));
%             dtc_mn(i, k) = mean(raw_var_sel(dtc_sel(:,k)));
%             dtc_sd(i, k) = std(raw_var_sel(dtc_sel(:,k)));
%             avg_md(i, k) = median(raw_var_sel(avg_sel(:,k)));
%             avg_mn(i, k) = mean(raw_var_sel(avg_sel(:,k)));
%             avg_sd(i, k) = std(raw_var_sel(avg_sel(:,k)));
%           end
          % Consider all variables as one
          dtc_sel = any(dtc_pl(i,:) <= raw_var_sel & raw_var_sel <= dtc_ph(i,:),2);
          avg_sel = any(avg_pl(i,:) <= raw_var_sel & raw_var_sel <= avg_ph(i,:),2);
          dtc_n(i) = sum(dtc_sel);
          avg_n(i) = sum(avg_sel);
          dtc_md(i,:) = median(raw_var_sel(dtc_sel,:));
          dtc_mn(i,:) = mean(raw_var_sel(dtc_sel,:));
          dtc_sd(i,:) = std(raw_var_sel(dtc_sel,:));
          avg_md(i,:) = median(raw_var_sel(avg_sel,:));
          avg_mn(i,:) = mean(raw_var_sel(avg_sel,:));
          avg_sd(i,:) = std(raw_var_sel(avg_sel,:));
        end
      end
      % Detection stats
%       bin.([lvar{j} '_dtc_pl']) = dtc_pl; % low percentile
%       bin.([lvar{j} '_dtc_ph']) = dtc_ph; % high percentile
      bin.([lvar{j} '_dtc_md']) =  dtc_md;% median
      bin.([lvar{j} '_dtc_mn']) = dtc_mn; % mean
%       bin.([lvar{j} '_dtc_sd']) = dtc_sd; % standard deviation
%       bin.([lvar{j} '_dtc_n']) = dtc_n;   % number of samples
      bin.([lvar{j} '_dtc_var']) = (dtc_ph - dtc_pl) ./ 2; % variance
      bin.([lvar{j} '_dtc_unc']) = bin.([lvar{j} '_dtc_var']) ./ sqrt(dtc_n); % uncertainty median
      bin.([lvar{j} '_dtc_se']) = dtc_sd ./ sqrt(dtc_n);  % uncertainty mean
      % Average Stats
      bin.([lvar{j} '_avg_pl']) = avg_pl; % low percentile
      bin.([lvar{j} '_avg_ph']) = avg_ph; % high percentile
      bin.([lvar{j} '_avg_md']) =  avg_md;% median
      bin.(lvar{j}) = avg_mn; % mean
      bin.([lvar{j} '_avg_sd']) = avg_sd; % standard deviation
      bin.([lvar{j} '_avg_n']) = avg_n;   % number of samples
%       bin.([lvar{j} '_avg_var']) = (bin.([lvar{j} '_avg_ph']) - bin.([lvar{j} '_avg_pl'])) ./ 2; % variance
%       bin.([lvar{j} '_avg_unc']) = bin.([lvar{j} '_avg_var']) ./ sqrt(bin.([lvar{j} '_avg_n'])); % uncertainty median 
      bin.([lvar{j} '_avg_se']) = bin.([lvar{j} '_avg_sd']) ./ sqrt(bin.([lvar{j} '_avg_n']));   % uncertainty mean   (coefficient of error | standard error)
    end
    % Remove empty lines
    sel = all(isnan(bin.(lvar{1})),2);
    for j=2:size(lvar,2)
      sel = sel & all(isnan(bin.(lvar{j})),2);
    end
    bin(sel,:) = [];
  case 'SB_IN_PRCTL'
    % Compute mean, standard deviation, and number of point
    %   for points within averaging percentile within the minute
    %   products recommended by SeaBASS
    for j=1:size(lvar,2)
      % Init Var In
      raw_dt = raw.dt;
      raw_var = raw.(lvar{j});

      % Init for averaging
      avg_pl = NaN(size(dt,1),size(raw_var,2));
      avg_ph = NaN(size(dt,1),size(raw_var,2));
      avg_mn = NaN(size(dt,1),size(raw_var,2));
      avg_sd = NaN(size(dt,1),size(raw_var,2));
      avg_n  = NaN(size(dt));
      
      % Init Parpool (needed for proper display of completion)
      %parpool(parallel_flag)
      
      % Init Display Progress
%       if verbose
%         fprintf('\t Completion: ');
%         showTimeToCompletion;
%         p = parfor_progress( size(dt,1) / 100);
%       end
%       startTime=tic;
      parfor (i=1:size(dt,1), parallel_flag)
%       for i=1:size(dt,1)
%         fprintf('\t%s %s\n', lvar{j}, datestr(dt(i))); 
        % Select minute to bin
        sel = dt(i) - bin_size/2 <= raw_dt & raw_dt < dt(i) + bin_size/2;
        % Unselect raw with NaN values
        foo = ~isnan(raw_var); % step required for 2d variables
        sel = sel & all(foo,2);
        % Check number of samples
        if sum(sel) > 0
          raw_var_sel = raw_var(sel,:);
          % Compute percentiles
          foo = prctile(raw_var_sel,[prctile_avg],1);
          avg_pl(i,:) = foo(1,:); % low percentile
          avg_ph(i,:) = foo(2,:); % high percentile
          
          % Consider all variables as one
          avg_sel = any(avg_pl(i,:) <= raw_var_sel & raw_var_sel <= avg_ph(i,:),2);
          avg_n(i) = sum(avg_sel);
          avg_mn(i,:) = mean(raw_var_sel(avg_sel,:));
          avg_sd(i,:) = std(raw_var_sel(avg_sel,:));
        end
        
        % Display Progress
%         if verbose
%             if ~mod(i, 100)
%               p = parfor_progress;
%               showTimeToCompletion( p/100, [], [], startTime );
%             end
%         end
      end
      % Average Stats
%       bin.([lvar{j} '_avg_pl']) = avg_pl; % low percentile
%       bin.([lvar{j} '_avg_ph']) = avg_ph; % high percentile
      bin.(lvar{j}) = avg_mn; % mean
      bin.([lvar{j} '_avg_sd']) = avg_sd; % standard deviation
      bin.([lvar{j} '_avg_n']) = avg_n;
    end
    % Remove empty lines
    sel = all(isnan(bin.(lvar{1})),2);
    for j=2:size(lvar,2)
      sel = sel & all(isnan(bin.(lvar{j})),2);
    end
    bin(sel,:) = [];
  case 'SB_ALL'
    % Compute mean, standard deviation, and number of point
    %   for all points within the minute
    %   products recommended by SeaBASS
    for j=1:size(lvar,2)
      % Init Var In
      raw_dt = raw.dt;
      raw_var = raw.(lvar{j});
      % Init var Out
      avg_mn = NaN(size(dt,1),size(raw_var,2));
      avg_sd = NaN(size(dt,1),size(raw_var,2));
      avg_n  = NaN(size(dt));
      
      % Init Display Progress
%       if verbose
%         fprintf('\t Completion: ');
%         showTimeToCompletion;
%         p = parfor_progress( size(dt,1) );
%       end
%       startTime=tic;
      parfor (i=1:size(dt,1), parallel_flag)
%       for i=1:size(dt,1)
        % Select minute to bin
        sel = dt(i) - bin_size/2 <= raw_dt & raw_dt < dt(i) + bin_size/2;
        % Unselect raw with NaN values
        foo = ~isnan(raw_var); % step required for 2d variables
        sel = sel & all(foo,2);
        % Check number of samples
        avg_n(i) = sum(sel);
        if sum(sel) > 0
%           avg_md(i,:) = median(raw_var(sel,:));
          avg_mn(i,:) = mean(raw_var(sel,:));
          avg_sd(i,:) = std(raw_var(sel,:));
        end
        % Display Progress
%         if verbose
%           p = parfor_progress;
%           showTimeToCompletion( p/100, [], [], startTime );
%         end
      end
      % Average Stats
      bin.(lvar{j}) = avg_mn; % mean
      bin.([lvar{j} '_avg_sd']) = avg_sd; % standard deviation
      bin.([lvar{j} '_avg_n']) = avg_n;   % number of samples
    end
    % Remove empty lines
    sel = all(isnan(bin.(lvar{1})),2);
    for j=2:size(lvar,2)
      sel = sel & all(isnan(bin.(lvar{j})),2);
    end
    bin(sel,:) = [];
%   case 'classic'
%     prctile_avg = [50 prctile_avg];
%     for j=1:size(lvar,2)
%       % Init Var In
%       raw_dt = raw.dt;
%       raw_var = raw.(lvar{j});
%       % Init Var Out
%       md = NaN(size(dt,1),size(raw.(lvar{j}),2)); % Median
%       pl = NaN(size(dt,1),size(raw.(lvar{j}),2));
%       ph = NaN(size(dt,1),size(raw.(lvar{j}),2));
%       sd = NaN(size(dt,1),size(raw.(lvar{j}),2));
%       n = NaN(size(dt,1),1);
%       parfor (i=1:size(dt,1), parallel_flag)
% %       for i=1:size(dt,1)
%         % Select minute to bin
%         sel = dt(i) - bin_size/2 <= raw_dt & raw_dt < dt(i) + bin_size/2;
%         % Unselect raw with NaN values
%         foo = ~isnan(raw_var); % step required for 2d variables
%         sel = sel & foo(:,1);
%         % Check number of samples
%         n(i) = sum(sel);
%         if n(i) > 0
%           foo = prctile(raw_var(sel,:),prctile_avg,1);
%           md(i,:) = foo(1,:); % median
%           pl(i,:) = foo(2,:); % low percentile
%           ph(i,:) = foo(3,:); % high percentile
%           sd(i,:) = std(raw_var(sel,:),0,1);
%         end
%       end
%       bin.(lvar{j}) =  md;
%       bin.([lvar{j} '_pl']) = pl;
%       bin.([lvar{j} '_ph']) = ph;
%       bin.([lvar{j} '_sd']) = sd;
%       bin.([lvar{j} '_n']) = n;
%     end
%     % Remove empty lines
%     sel = bin.([lvar{1} '_n']) == 0;
%     for j=2:size(lvar,2)
%       sel = sel & bin.([lvar{j} '_n']) == 0;
%     end
%     bin(sel,:) = [];
%   case 'average'
%     for j=1:size(lvar,2)
%       % Init Var In
%       raw_dt = raw.dt;
%       raw_var = raw.(lvar{j});
%       % Init Var Out
%       md = NaN(size(dt,1),size(raw.(lvar{j}),2)); % Median
%       mn = NaN(size(dt,1),size(raw.(lvar{j}),2)); % Mean
%       parfor (i=1:size(dt,1), parallel_flag)
% %       for i=1:size(dt,1)
%         % Select minute to bin
%         sel = dt(i) - bin_size/2 <= raw_dt & raw_dt < dt(i) + bin_size/2;
%         % Unselect raw with NaN values
%         foo = ~isnan(raw_var); % step required for 2d variables
%         sel = sel & foo(:,1);
%         % Check number of samples
%         n = sum(sel);
%         if n > 0
%           md(i,:) = median(raw_var(sel,:));
%           mn(i,:) = mean(raw_var(sel,:));
%         end
%       end
%       bin.(lvar{j}) =  md; % median
%       bin.([lvar{j} '_mn']) = mn; % mean
%     end
%     % Remove empty lines
%     sel = isnan(bin.(lvar{1})(1,:));
%     for j=2:size(lvar,2)
%       sel = sel & isnan(bin.(lvar{j})(1,:));
%     end
%     bin(sel,:) = [];
  otherwise
    error('Unknown mode of binning data\n');
end

if verbose; fprintf('Done\n'); end

end



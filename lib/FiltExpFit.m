function [filt_avg, filt_stat] = FiltExpFit(filt_avg, filt, filt_st, filt_end)
% Author: Guillaume bourdin
% Date: June 25, 2021
%
% Based on method in: Dall'Olmo G., R. J. W. Brewin, F.Nencioli, E. Organelli, 
%     I. Lefering, D. McKee, R. Rottgers, C. Mitchell, E. Boss, A. Bricaud, 
%     and G. Tilstone , 2017. Determination of the absorption coefficient of 
%     chromophoric dissolved organic matter from underway spectrophotometry. 
%     Optics Express, 25, 24, https://doi.org/10.1364/OE.25.0A1079
%
% Exponential fit to filter event
%%
[nspect, nlamda] = size(filt_avg.beta);

%setting options for fmisearch
opts = optimset('fminsearch');
opts = optimset(opts, 'MaxIter', 5000, 'Display', 'none');
opts = optimset(opts, 'MaxFunEvals', 1000);   % usually 100*number of params
opts = optimset(opts, 'TolFun', 1e-6);

filt_stat = table(filt_avg.dt, 'VariableNames', {'dt'});
filt_stat.slope = NaN(nspect, size(filt.beta, 2));

% figure()
% scatter(datetime(filt.dt, 'ConvertFrom', 'datenum'), filt.beta(:,1), 10, 'filled')
% vline(filt_st, '-g')
% vline(filt_end, '-r')

% compute exp fit at each wavelength and each filter event
for i = progress(1:nspect)
%   select filter event
  sel_filt = filt_st(i) <= filt.dt & filt.dt <= filt_end(i);
  if sum(sel_filt) > 3
    foo = filt(sel_filt,:);
    if max(foo.dt) - min(foo.dt) > 0.0007 * 4
      avg_beta = median(foo.beta,2,'omitnan');
      foo(foo.dt < max(foo.dt(avg_beta == max(avg_beta))),:) = [];
      if max(foo.dt) - min(foo.dt) > 0.0007 * 4
        % remove short peaks and smooth signal
        j = 1;
        deriv = diff(foo.beta);
        deriv_neg = deriv;
        deriv_neg(deriv > 0) = NaN;
        lim_neg = median(deriv_neg, 'omitnan') * 3;
        foo.beta([deriv; zeros(1, size(foo.beta, 2))] < lim_neg & ...
          [zeros(1, size(foo.beta, 2)); deriv] > 0) = NaN;
        % remove duplicates
        [~, L, ~] = unique(foo.dt,'first');
        indexToDump = not(ismember(1:numel(foo.dt), L));
        foo(indexToDump, :) = [];
        while j < 5 && any(isnan(foo.beta(:)))
          foo.beta = fillmissing(foo.beta,'linear','SamplePoints', foo.dt);
%           figure();
%           plot(foo.dt, foo.beta)
          deriv = diff(foo.beta);
          deriv_neg = deriv;
          deriv_neg(deriv > 0) = NaN;
          lim_neg = median(deriv_neg, 'omitnan') * 1.5;
          foo.beta([deriv; zeros(1, size(foo.beta, 2))] < lim_neg & ...
            [zeros(1, size(foo.beta, 2)); deriv] > 0) = NaN;
          j = j + 1;
        end
  %       foo.beta([diff(foo.beta); zeros(1, size(foo.beta, 2))] < -10 | ...
  %         [zeros(1, size(foo.beta, 2)); diff(foo.beta)] > 10) = NaN;
  %       foo([diff(foo.beta); 0] > 10^-3, :) = [];
        beta_temp = NaN(1, size(foo.beta, 2));
        beta_avg_sd_temp = NaN(1, size(foo.beta, 2));
        slope_temp = NaN(1, size(foo.beta, 2));
        fval_temp = NaN(1, size(foo.beta, 2));
        exitflag_temp = false(1, size(foo.beta, 2));
        parfor j = 1:nlamda
          foo_wl = foo.beta(:,j);
          if all(isfinite(foo_wl)) && all(~isnan(foo_wl))
            % define exponential function
            expfun = @(p, xd) p(1) * exp(p(2) * (xd - min(foo.dt))) + p(3);
            % define x0
            x0 = [max(foo_wl) - min(foo_wl) -830 min(foo_wl)];
            % define weight
            weig = 1 - 0.5 * (1:size(foo, 1))' / size(foo, 1);
            % define error function: sum_err/std
            errfun = @(p) sum(abs(expfun(p, foo.dt) - foo_wl) ./ weig); 
            % run the minimizer
            [pfit, FVAL, EXITFLAG] = fminsearch(errfun, x0, opts);

  %           figure(j);  hold on;
  % %           plot(foo.dt, foo_wl)
  %           scatter(foo.dt, foo_wl, 10, 'filled')
  %           vline(filt_st(i), '-g')
  %           vline(filt_end(i), '-r')
  %           plot(foo.dt, expfun(pfit, foo.dt), 'r-');

            % populate table and propagate error
            beta_temp(j) = pfit(3);
            beta_avg_sd_temp(j) = sum(abs(expfun(pfit, foo.dt) - foo_wl) / ...
              sum(~isnan(foo_wl)));
            slope_temp(j) = pfit(2);
            fval_temp(j) = FVAL;
            exitflag_temp(j) = EXITFLAG;
  %           filt_avg.beta(i,j) = pfit(3);
  %           filt_avg.beta_avg_sd(i,j) = sum(abs(expfun(pfit, foo.dt) - foo_wl) / ...
  %             sum(~isnan(foo_wl)));
  %           filt_stat.slope(i,j) = pfit(2);
  %           filt_stat.fval(i,j) = FVAL;
  %           filt_stat.exitflag(i,j) = EXITFLAG;
          end
        end
        filt_avg.beta(i,:) = beta_temp;
        filt_avg.beta_avg_sd(i,:) = beta_avg_sd_temp;
        filt_stat.slope(i,:) = slope_temp;
        filt_stat.fval(i,:) = fval_temp;
        filt_stat.exitflag(i,:) = exitflag_temp;
        filt_avg.beta_avg_n(i) = median(sum(~isnan(foo.beta)),2, 'omitnan');
      end
    end
  end
end
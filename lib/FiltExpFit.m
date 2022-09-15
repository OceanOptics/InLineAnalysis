function [filt_avg, filt_stat] = FiltExpFit(filt_avg, filt_good, filt_bad, filt_st, filt_end)
% Author: Guillaume bourdin
% Date: June 25, 2021
%
% Based on method in: Dall’Olmo, G., Westberry, T.K., Behrenfeld, M.J., Boss, 
%                     E., Slade, W.H., 2009. Direct contribution of phytoplankton-sized particles 
%                     to optical backscattering in the open ocean. Biogeosciences Discuss 6, 291–340. 
%                     https://doi.org/10.5194/bgd-6-291-2009
%
% Exponential fit to filter event
%%
[nspect, nlamda] = size(filt_avg.beta);

% figure(1);  hold on;
% scatter(filt.dt, filt.beta, 10, 'filled')
% hold off

%setting options for fmisearch
opts = optimset('fminsearch');
opts = optimset(opts, 'MaxIter', 5000, 'Display', 'none');
opts = optimset(opts, 'MaxFunEvals', 1000);   % usually 100*number of params
opts = optimset(opts, 'TolFun', 1e-6);

filt = sortrows([filt_good; filt_bad], 'dt');

filt_stat = table(filt_avg.dt, 'VariableNames', {'dt'});
filt_stat.slope = NaN(nspect, size(filt.beta, 2));
% compute exp fit at each wavelength and each filter event
for i = progress(1:nspect)
%   select filter event
  sel_filt = filt_st(i) <= filt.dt & filt.dt <= filt_end(i);
  sel_filt_good = filt_st(i) <= filt_good.dt & filt_good.dt <= filt_end(i);
  if sum(sel_filt) > 3
    foo = filt(sel_filt,:);
    if max(foo.dt) - min(foo.dt) > 0.0007 * 2
      avg_beta = median(foo.beta,2,'omitnan');
%       cut_tail = foo.dt > max(foo.dt(avg_beta == min(avg_beta)));
%       foo(cut_tail,:) = [];
%       avg_beta(cut_tail,:) = [];
      foo(foo.dt < max(foo.dt(avg_beta == max(avg_beta))),:) = [];
      if max(foo.dt) - min(foo.dt) > 0.0007 * 2
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
        % remove spikes
        while j < 5 && any(isnan(foo.beta(:)))
          foo.beta = fillmissing(foo.beta,'linear','SamplePoints', foo.dt);
          deriv = diff(foo.beta);
          deriv_neg = deriv;
          deriv_neg(deriv > 0) = NaN;
          lim_neg = median(deriv_neg, 'omitnan') * 1.5;
          foo.beta([deriv; zeros(1, size(foo.beta, 2))] < lim_neg & ...
            [zeros(1, size(foo.beta, 2)); deriv] > 0) = NaN;
          j = j + 1;
        end
        foo.beta = fillmissing(foo.beta,'linear','SamplePoints', foo.dt);
        beta_temp = NaN(1, size(foo.beta, 2));
        beta_avg_sd_temp = NaN(1, size(foo.beta, 2));
        slope_temp = NaN(1, size(foo.beta, 2));
        fval_temp = NaN(1, size(foo.beta, 2));
        exitflag_temp = false(1, size(foo.beta, 2));
        parfor j = 1:nlamda
          foo_wl = foo.beta(:,j);
          if all(isfinite(foo_wl)) && all(~isnan(foo_wl))
            expfun = @(p, xd) p(1) * exp(p(2) * (xd - min(foo.dt))) + p(3); % define exponential function
            x0 = [max(foo_wl) - min(foo_wl) -830 min(foo_wl)]; % define x0
%             weig = 1 - 0.5 * (1:size(foo, 1))' / size(foo, 1); % define weight
            weig = (1:size(foo, 1))' / size(foo, 1); % define weight
            errfun = @(p) sum(abs(expfun(p, foo.dt) - foo_wl) ./ weig); % define error function: sum_err/std
            [pfit, FVAL, EXITFLAG] = fminsearch(errfun, x0, opts); % run the minimizer

%             f1 = figure(1);  hold on;
%   %           plot(foo.dt, foo_wl)
%             sc = scatter(foo.dt, foo_wl, 10, 'filled');
%             vline(filt_st(i), '-g')
%             vline(filt_end(i), '-r')
%             plot(foo.dt, expfun(pfit, foo.dt), 'Color', sc.CData);

            % populate table and propagate error
            beta_temp(j) = pfit(3);
            beta_avg_sd_temp(j) = sum(abs(expfun(pfit, foo.dt) - foo_wl) / ...
              sum(~isnan(foo_wl)));
            slope_temp(j) = pfit(2);
            fval_temp(j) = FVAL;
            exitflag_temp(j) = EXITFLAG;
          end
        end          
        % remove fit when constant > 25% percentile
        if sum(sel_filt_good) > 3
          foo_good = filt_good(sel_filt_good,:);
          exitflag_temp(beta_temp > prctile(foo_good.beta, 25)) = false;
        else
          exitflag_temp(beta_temp > prctile(foo.beta, 25)) = false;
        end
        % remove fit when constant < prctile(foo.beta, 25)/2
        exitflag_temp(beta_temp < prctile(foo.beta, 25)/2) = false;
        % eliminate bad fit when error > 20%
        exitflag_temp(beta_avg_sd_temp./beta_temp > 0.2) = false;
%         clf(f1)
        
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
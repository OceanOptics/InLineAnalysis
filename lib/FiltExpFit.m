function [filt_avg, filt_stat] = FiltExpFit(var, filt_avg, filt_good, filt_bad, filt_st, filt_end)
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
[nspect, nlamda] = size(filt_avg.(var));

% figure(1);  hold on;
% scatter(filt.dt, filt.(var), 10, 'filled')
% hold off

%setting options for fmisearch
opts = optimset('fminsearch');
opts = optimset(opts, 'MaxIter', 5000, 'Display', 'none');
opts = optimset(opts, 'MaxFunEvals', 1000);   % usually 100*number of params
opts = optimset(opts, 'TolFun', 1e-6);

filt = sortrows([filt_good; filt_bad], 'dt');

filt_stat = table(filt_avg.dt, 'VariableNames', {'dt'});
filt_stat.slope = NaN(nspect, size(filt.(var), 2));
filt_stat.slope = NaN(nspect, size(filt.(var), 2));
% compute exp fit at each wavelength and each filter event
for i = progress(1:nspect)
%   select filter event
  sel_filt = filt_st(i) <= filt.dt & filt.dt <= filt_end(i);
  sel_filt_good = filt_st(i) <= filt_good.dt & filt_good.dt <= filt_end(i);
  if sum(sel_filt) > 3
    filt_subset = filt(sel_filt,:);
    if max(filt_subset.dt) - min(filt_subset.dt) > minutes(2)
      avg_var = median(filt_subset.(var),2,'omitnan');
%       cut_tail = foo.dt > max(foo.dt(avg_var == min(avg_var)));
%       foo(cut_tail,:) = [];
%       avg_var(cut_tail,:) = [];
      filt_subset(filt_subset.dt < max(filt_subset.dt(avg_var == max(avg_var))),:) = [];
      if max(filt_subset.dt) - min(filt_subset.dt) > minutes(2)
        % remove short peaks and smooth signal
        j = 1;
        deriv = diff(filt_subset.(var));
        deriv_neg = deriv;
        deriv_neg(deriv > 0) = NaN;
        lim_neg = median(deriv_neg, 'omitnan') * 3;
        filt_subset.(var)([deriv; zeros(1, size(filt_subset.(var), 2))] < lim_neg & ...
          [zeros(1, size(filt_subset.(var), 2)); deriv] > 0) = NaN;
        % remove duplicates
        [~, L, ~] = unique(filt_subset.dt,'first');
        indexToDump = not(ismember(1:numel(filt_subset.dt), L));
        filt_subset(indexToDump, :) = [];
        % remove spikes
        while j < 5 && any(isnan(filt_subset.(var)(:)))
          filt_subset.(var) = fillmissing(filt_subset.(var),'linear','SamplePoints', filt_subset.dt);
          deriv = diff(filt_subset.(var));
          deriv_neg = deriv;
          deriv_neg(deriv > 0) = NaN;
          lim_neg = median(deriv_neg, 'omitnan') * 1.5;
          filt_subset.(var)([deriv; zeros(1, size(filt_subset.(var), 2))] < lim_neg & ...
            [zeros(1, size(filt_subset.(var), 2)); deriv] > 0) = NaN;
          j = j + 1;
        end
        filt_subset.(var) = fillmissing(filt_subset.(var),'linear','SamplePoints', filt_subset.dt);
        filt_subset.dt = juliandate(filt_subset.dt);
        var_temp = NaN(1, size(filt_subset.(var), 2));
        var_avg_sd_temp = NaN(1, size(filt_subset.(var), 2));
        slope_temp = NaN(1, size(filt_subset.(var), 2));
        fval_temp = NaN(1, size(filt_subset.(var), 2));
        exitflag_temp = false(1, size(filt_subset.(var), 2));
        parfor j = 1:nlamda
          data2fit = filt_subset.(var)(:,j);
          if all(isfinite(data2fit)) && all(~isnan(data2fit))
            expfun = @(p, xd) p(1) * exp(p(2) * (xd - min(filt_subset.dt))) + p(3); % define exponential function
            x0 = [max(data2fit) - min(data2fit) -830 min(data2fit)]; % define x0
            % % weig = 1 - 0.5 * (1:size(foo, 1))' / size(foo, 1); % define weight
            % weig = (1:size(filt_subset, 1))' / size(filt_subset, 1); % define weight
            weig = ones(size(filt_subset, 1), 1); % define weight
            errfun = @(p) sum(abs(expfun(p, filt_subset.dt) - data2fit) .* weig); % define error function: sum_err/std
            [pfit, FVAL, EXITFLAG] = fminsearch(errfun, x0, opts); % run the minimizer

  %           figure(1);  hold on;
  % %           plot(foo.dt, foo_wl)
  %           sc = scatter(datetime(filt_subset.dt, 'ConvertFrom', 'juliandate'), data2fit, 10, 'filled');
  %           % vline(filt_st(i), '-g')
  %           % vline(filt_end(i), '-r')
  %           plot(datetime(filt_subset.dt, 'ConvertFrom', 'juliandate'), expfun(pfit, filt_subset.dt), 'Color', sc.CData);
  %           visProd3D(430:10:700, filt_subset.dt, filt_subset.beta, false, 'Wavelength', false, 68);
  %           pause (0.5)
  %           clf

            % populate table and propagate error
            var_temp(j) = pfit(3);
            var_avg_sd_temp(j) = sum(abs(expfun(pfit, filt_subset.dt) - data2fit) / sum(~isnan(data2fit)));
            slope_temp(j) = pfit(2);
            fval_temp(j) = FVAL;
            exitflag_temp(j) = EXITFLAG;
          end
        end

        % lambda = 430:10:700;
        % f2 = figure(2);
        % subplot(4,1,1); scatter(lambda(exitflag_temp), var_temp(exitflag_temp))
        % hold on; scatter(lambda(~exitflag_temp), var_temp(~exitflag_temp))
        % subplot(4,1,2); scatter(lambda(exitflag_temp), var_avg_sd_temp(exitflag_temp))
        % hold on; scatter(lambda(~exitflag_temp), var_avg_sd_temp(~exitflag_temp))
        % subplot(4,1,3); scatter(lambda(exitflag_temp), slope_temp(exitflag_temp))
        % hold on; scatter(lambda(~exitflag_temp), slope_temp(~exitflag_temp))
        % subplot(4,1,4); scatter(lambda(exitflag_temp), fval_temp(exitflag_temp))
        % hold on; scatter(lambda(~exitflag_temp), fval_temp(~exitflag_temp))

        % remove fit when fval > 25% percentile
        if sum(sel_filt_good) > 3
          foo_good = filt_good(sel_filt_good,:);
          exitflag_temp(var_temp > median(foo_good.(var))) = false;
        else
          exitflag_temp = false(size(exitflag_temp));
        end
        % remove fit when variance < prctile(filt_subset.(var), 25)/2
        exitflag_temp(var_temp < prctile(filt_subset.(var), 25)/2) = false;
        % eliminate bad fit when error > 20%
        exitflag_temp(var_avg_sd_temp./var_temp > 0.2) = false;
        % eliminate bad fit when any slope > 0
        exitflag_temp(slope_temp > 0) = false;
%         clf(f1)
        
        % f3 = figure(3);
        % subplot(4,1,1); scatter(lambda(exitflag_temp), var_temp(exitflag_temp))
        % hold on; scatter(lambda(~exitflag_temp), var_temp(~exitflag_temp))
        % subplot(4,1,2); scatter(lambda(exitflag_temp), var_avg_sd_temp(exitflag_temp))
        % hold on; scatter(lambda(~exitflag_temp), var_avg_sd_temp(~exitflag_temp))
        % subplot(4,1,3); scatter(lambda(exitflag_temp), slope_temp(exitflag_temp))
        % hold on; scatter(lambda(~exitflag_temp), slope_temp(~exitflag_temp))
        % subplot(4,1,4); scatter(lambda(exitflag_temp), fval_temp(exitflag_temp))
        % hold on; scatter(lambda(~exitflag_temp), fval_temp(~exitflag_temp))
        % pause(1)
        % clf(f2); clf(f3);

        filt_avg.(var)(i,:) = var_temp;
        filt_avg.([var '_avg_sd'])(i,:) = var_avg_sd_temp;
        filt_stat.slope(i,:) = slope_temp;
        filt_stat.fval(i,:) = fval_temp;
        filt_stat.exitflag(i,:) = exitflag_temp;
        filt_avg.([var '_avg_n'])(i) = median(sum(~isnan(filt_subset.(var))),2, 'omitnan');
      end
    end
  end
end
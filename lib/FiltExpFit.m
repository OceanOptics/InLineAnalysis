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
opts = optimset(opts, 'MaxIter', 5000);
opts = optimset(opts, 'MaxFunEvals', 1000);   % usually 100*number of params
opts = optimset(opts, 'TolFun', 1e-6);

filt_stat = table(filt_avg.dt, 'VariableNames', {'dt'});
filt_stat.slope = NaN(nspect, size(filt.beta, 2));
% compute exp fit at each wavelength and each filter event
for i = progress(1:nspect)
  % select filter event
  sel_filt = filt_st(i) <= filt.dt & filt.dt <= filt_end(i);
  if sum(sel_filt) > 240
    foo = filt(sel_filt,:);
    avg_beta = median(foo.beta,2,'omitnan');
    foo(foo.dt < max(foo.dt(avg_beta == max(avg_beta))),:) = [];
    % remove short peaks and smooth signal
    foo.beta([diff(foo.beta); zeros(1, size(foo.beta, 2))] < -10 | ...
      [zeros(1, size(foo.beta, 2)); diff(foo.beta)] > 10) = NaN;
    foo.beta = fillmissing(foo.beta,'linear','SamplePoints', foo.dt);
%     foo([diff(foo.beta); 0] > 10^-3, :) = [];
    for j = 1:nlamda
      if all(isfinite(foo.beta(:,j))) && all(~isnan(foo.beta(:,j)))
        % define exponential function
        expfun = @(p, xd) p(1) * exp(p(2) * (xd - min(foo.dt))) + p(3);
        % define x0
        x0 = [max(foo.beta(:,j)) - min(foo.beta(:,j)) -830 min(foo.beta(:,j))];
        % define weight
        weig = 1 - 0.5 * (1:size(foo, 1))' / size(foo, 1);
        % define error function: sum_err/std
        errfun = @(p) sum(abs(expfun(p, foo.dt) - foo.beta(:,j)) ./ weig); 
        % run the minimizer
        [pfit, FVAL, EXITFLAG] = fminsearch(errfun, x0, opts);

%         figure(j);  hold on;
%         plot(foo.dt, foo.beta)
%         vline(filt_st(i), '-b')
%         vline(filt_end(i), '-r')
%         plot(foo.dt, expfun(pfit, foo.dt), 'r-');
        
        % populate table and propagate error
        filt_avg.beta(i,j) = pfit(3);
        filt_avg.beta_avg_sd(i,j) = sum(abs(expfun(pfit, foo.dt) - foo.beta(:,j)) / sum(~isnan(foo.beta(:,j))));
        filt_stat.slope(i,j) = pfit(2);
        filt_stat.fval(i,j) = FVAL;
        filt_stat.exitflag(i,j) = EXITFLAG;
      end
    end
    filt_avg.beta_avg_n(i) = median(sum(~isnan(foo.beta)),2, 'omitnan');
  end
end
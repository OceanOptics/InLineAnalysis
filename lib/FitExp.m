function [y_intercp, base, SSE, MSE, RMSE] = FitExp(lambda, meas, meas_sd)
% Author: Guillaume bourdin
% Date: May 19, 2021
%
% Exponential fit to ag and cg
%%
lambda = lambda(:)';
[nspect, nlambd] = size(meas);

if length(lambda)~=nlambd, error('Invalid lambda'); end

% x0 = [mean(meas(:, abs(lambda-440)==min(abs(lambda-440)))) -0.015];
x0 = [meas(:, abs(lambda-440)==min(abs(lambda-440))) repmat(-0.015, nspect, 1)];

y_intercp = NaN(nspect,1);
base = NaN(nspect,1);
SSE = NaN(nspect,1);
MSE = NaN(nspect,1);
RMSE = NaN(nspect,1);

%setting options for fmisearch
opts = optimset('fminsearch');
opts = optimset(opts, 'MaxIter', 500000);
opts = optimset(opts, 'MaxFunEvals', 100000);   % usually 100*number of params
opts = optimset(opts, 'TolFun', 1e-9);

parfor i = 1:nspect
  if all(isfinite(meas(i,:))) && ~any(isnan(meas(i, :)), 2)
    % define exponential function
    expfun = @(p, xd) p(1) * exp(p(2) * (xd - 440));
    % define sum/std error
    errfun = @(p) sum(abs((expfun(p, lambda) - meas(i, :)) ./ meas_sd(i,:)));
    %run the minimizer
    pfit = fminsearch(errfun, x0(i,:), opts);
    y_intercp(i) = pfit(1);
    base(i) = pfit(2);
    % evaluate goodness
    SSE(i) = sum((expfun(pfit, lambda) - meas(i, :)) .^2);
    MSE(i) = SSE(i) / (size(lambda, 2) - size(x0, 2));
    RMSE(i) = sqrt(MSE(i)); 
%     figure(i)
%     plot(lambda, meas(i, :),'bo');  hold on;   %plot your raw data
%     plot(lambda, expfun(pfit, lambda), 'r-');  %plot the fit data
  end
end

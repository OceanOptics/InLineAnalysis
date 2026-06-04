function [y_intercp, base, SSE, MSE, RMSE, nRMSE] = FitExp(lambda, data, data_sd)
  % Author: Guillaume bourdin
  % Date: May 19, 2021
  %
  % Exponential fit to ag and cg
  %%
  lambda = lambda(:)';
  [n, nlambd] = size(data);
  
  if length(lambda)~=nlambd, error('Invalid lambda'); end
  if nargin < 3; data_sd = ones(size(data)); elseif isempty(data_sd); data_sd = ones(size(data)); end
  
  % x0 = [mean(data(:, abs(lambda-440)==min(abs(lambda-440)))) -0.015];
  x0 = [data(:, abs(lambda-440)==min(abs(lambda-440))) repmat(-0.015, n, 1)];
  
  y_intercp = NaN(n,1);
  base = NaN(n,1);
  SSE = NaN(n,1);
  MSE = NaN(n,1);
  RMSE = NaN(n,1);
  nRMSE = NaN(n,1);
  
  % setting options for fmisearch
  opts = optimset('fminsearch');
  opts = optimset(opts, 'MaxIter', 500000, 'Display', 'none');
  opts = optimset(opts, 'MaxFunEvals', 100000);   % usually 100*number of params
  opts = optimset(opts, 'TolFun', 1e-9);
  
  % define exponential function
  expfun = @(p, xd) p(1) * exp(p(2) * (xd - 440));

  if n > 5000
    parfor i = 1:n
      if all(isfinite(data(i,:))) && ~any(isnan(data(i, :)), 2)
        % define sum/std error
        errfun = @(p) sum(abs((expfun(p, lambda) - data(i, :)) ./ data_sd(i,:)));
        % run the minimizer
        pfit = fminsearch(errfun, x0(i,:), opts);
        y_intercp(i) = pfit(1);
        base(i) = pfit(2);
        % evaluate fit
        SSE(i) = sum((expfun(pfit, lambda) - data(i, :)) .^2);
        MSE(i) = SSE(i) / (size(lambda, 2) - size(x0, 2));
        RMSE(i) = sqrt(MSE(i));
        nRMSE(i) = RMSE(i)/mean(data(i, :))*100;
      end
    end
  else
    for i = progress(1:n)
      if all(isfinite(data(i,:))) && ~any(isnan(data(i, :)), 2)
        % define sum/std error
        errfun = @(p) sum(abs((expfun(p, lambda) - data(i, :)) ./ data_sd(i,:)));
        % run the minimizer
        pfit = fminsearch(errfun, x0(i,:), opts);
        y_intercp(i) = pfit(1);
        base(i) = pfit(2);
        % evaluate fit
        SSE(i) = sum((expfun(pfit, lambda) - data(i, :)) .^2);
        MSE(i) = SSE(i) / (size(lambda, 2) - size(x0, 2));
        RMSE(i) = sqrt(MSE(i));
        nRMSE(i) = RMSE(i)/mean(data(i, :))*100;
      end
    end
  end


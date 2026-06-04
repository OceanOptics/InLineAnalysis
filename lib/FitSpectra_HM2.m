function [meas0, gamma, fiterr] = FitSpectra_HM2(lambda, meas)
  % Fit power law to data
  lambda = lambda(:)';
  [nspect, nlambd] = size(meas);
  
  if length(lambda)~=nlambd, error('Invalid lambda'); end
  
  meas0 = NaN(nspect,1);
  gamma = NaN(nspect,1);
  fiterr = NaN(nspect,1);
  
  %setting options for fmisearch
  opts = optimset('fminsearch');
  opts = optimset(opts,'MaxIter',4000, 'Display', 'none');
  opts = optimset(opts,'MaxFunEvals',2000);   % usually 100*number of params
  opts = optimset(opts,'TolFun',1e-9);
  %opts = optimset('LevenbergMarquardt','on');
  if nspect > 30000
    parfor k = 1:nspect
      if all(isfinite(meas(k,:))) && ~any(isnan(meas(k,:)), 2)
        % guess for paramters (data at lambda0, beamc slope)
        x0 = [1.0 -0.8];
        % minimization routine a la Nelder Mead
        [x, fiterr(k)] = fminsearch(@least_square, x0, opts, meas(k,:), lambda);    
        meas0(k) = x(1);
        gamma(k) = x(2);
      end
    end
  else
    for k = progress(1:nspect)
      if all(isfinite(meas(k,:))) && ~any(isnan(meas(k,:)), 2)
        % guess for paramters (data at lambda0, beamc slope)
        x0 = [1.0 -0.8];
        % minimization routine a la Nelder Mead
        [x, fiterr(k)] = fminsearch(@least_square, x0, opts, meas(k,:), lambda);
        meas0(k) = x(1);
        gamma(k) = x(2);
      end
    end
  end
return


function y = least_square(x0, spec, lambda)
  y = sum(((spec - x0(1) .* (532 ./ lambda) .^ x0(2))) .^ 2);
return
% Fit power law to data
function [meas0, gamma, fiterr] = FitSpectra_HM2(lambda, meas)
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

parfor k = 1:nspect
  if all(isfinite(meas(k,:))) && ~any(isnan(meas(k,:)), 2)
    % guess for paramters (data at lambda0, beamc slope)
    x0 = [1.0 -0.8];

    % minimization routine a la Nelder Mead
    [x, fiterr(k)] = fminsearch(@least_square, x0, opts, meas(k,:), lambda); %lambda,lambda0

    meas0(k) = x(1);
    gamma(k) = x(2);

    %disp(['  Fitting power-law model: k=' ...
        %num2str(k) ' meas0=' num2str(x(1)) ' gamma=' num2str(x(2)) ' err=' num2str(fiterr(k))]);
  end
end
return 

% function err = power_law(x, meas, lambda, lambda0);
% 
% 	yhat = x(1)*(lambda./lambda0).^(-x(2));
% 
%     err = sum(abs(meas-yhat));
% 
% return

function y = least_square(x0, spec, lambda)
y = sum(((spec - x0(1) .* (532 ./ lambda) .^ x0(2))) .^ 2);
return
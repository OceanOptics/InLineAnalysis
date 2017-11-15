
% Fit power law to beamc data
function [cp0, gamma, fiterr] = FitSpectra_HM2(lambda, cpmeas)

lambda = lambda(:)';

[nspect, nlambd] = size(cpmeas);

if length(lambda)~=nlambd, error('Invalid lambda'); end

cp0 = repmat(nan, nspect,1);
gamma = repmat(nan, nspect,1);
fiterr = repmat(nan, nspect,1);

%setting options for fmisearch
opts = optimset('fminsearch');      
opts = optimset(opts,'MaxIter',4000); 
opts = optimset(opts,'MaxFunEvals',2000);   % usually 100*number of params
opts = optimset(opts,'TolFun',1e-9);
%opts = optimset('LevenbergMarquardt','on');

for k = 1:nspect

    if all(isfinite(cpmeas(k,:)))
        % guess for paramters (beamc at lambda0, beamc slope)
        x0 = [1.0 0.8];

        % minimization routine a la Nelder Mead
        [x, fiterr(k)] = fminsearch(@least_square, x0, opts, cpmeas(k,:), lambda); %lambda,lambda0

        cp0(k) = x(1);
        gamma(k) = x(2);
		
        %disp(['  Fitting power-law model: k=' ...
            %num2str(k) ' cp0=' num2str(x(1)) ' gamma=' num2str(x(2)) ' err=' num2str(fiterr(k))]);
			
    end
    
end

return 


% function err = power_law(x, cpmeas, lambda, lambda0);
% 
% 	cphat = x(1)*(lambda./lambda0).^(-x(2));
% 
%     err = sum(abs(cpmeas-cphat));
% 
% return

function y = least_square(x0,cp_spec,lambda)
y=sum(((cp_spec-x0(1).*(532./lambda).^x0(2))).^2);
return
    
function [ad, aph] = partition_ap_vec(ap_i, lambda, prcntl)
% 
% Description: This code partitions the spectral total particulate absorption coefficient, ap(lambda) into
% phytoplankton, aph(lambda), and non-algal, ad(lambda), components.
% 
% Input:  ap_i, spectrum of total particulate absorption coefficient,  must be a column vector
%         lamda, wavelengths of ap_i, must be a column vector, and must be hyperspectral from 400-700 nm
%         prcntl, percentiles of feasible solutions, can be a scalor or vector
% Output: ad, non-algal (detrital) particulate absorption spectrum
%         aph, phytoplankton absorption spectrum
% 
% If no feasible solutions are found, the output values are NaNs.
% 
% Author:  Guangming Zheng, gzheng@ucsd.edu
% Version: 0.0, released date
% 
% An example to use this code:
% 
% [ad_scm, aph_scm] = partition_ap(ap_input, lambda_input, [10 50 90]);
%             , where lambda_input is a column vector, [400:700]'
%                     ap_input is the total non-water absorption coefficient at the same wavelengths as lambda_input
%                     [10 50 90] specifies the output values as the 10th, 50th (median, i.e., optimal solution), and 90th percentiles of all feasible solutions.
%                                If you only want the optimal solution, just ignore this input option.
% 
% Reference: 
% Zheng, G., and D. Stramski (2013), A model based on stacked-constraints
% approach for partitioning the light absorption coefficient of seawater 
% into phytoplankton and non-phytoplankton components, J. Geophys. Res. Oceans, 
% 118, 2155?2174, doi:10.1002/jgrc.20115.
% 
% Zheng, G., and D. Stramski (2013), A model for partitioning the light
% absorption coefficient of suspended marine particles into phytoplankton 
% and nonalgal components, J. Geophys. Res. Oceans, 118, 2977?2991, 
% doi:10.1002/jgrc.20206.
%
%%
    m = length(ap_i);
    [fsbl_Ad,fsbl_Sd] = scm_ap(ap_i,lambda); % solve for all feasible solutions of Ad and Sd
    n = length(fsbl_Ad);
    adtemp  = repmat(fsbl_Ad,1,m) .* exp(-fsbl_Sd*lambda'); % all feasible solution of ad
    aphtemp = repmat(ap_i',n,1) - adtemp; % all feasible solution of aph
    ad      = prctile(adtemp,prcntl);  % percentiles of all feasible solutions of ad
    aph     = prctile(aphtemp,prcntl); % percentiles of all feasible solutions of aph
end


function [Ad,Sd] = scm_ap(ap_i, lambda)
% Input:  ap_i, spectrum of total particulate absorption coefficient
%         lamda, wavelengths of ap_i
% Output: Ad, amplitude of spectral non-algal particulate absorption coefficient, ad
%         Sd, slope of ad
%         Ad and Sd are vectors representing values for all feasible solutions
%
%% find the index of wavelengths involved in all constraints
idx400 = find(lambda==400);
idx412 = lambda==412;
idx420 = find(lambda==420);
idx430 = find(lambda==430);
idx443 = lambda==443;
idx450 = find(lambda==450);
idx467 = lambda==467;
idx490 = lambda==490;
idx500 = find(lambda==500);
idx510 = lambda==510;
idx550 = find(lambda==550);
idx555 = lambda==555;
idx630 = lambda==630;
idx650 = lambda==650;
idx670 = lambda==670;
idx700 = lambda==700;

% create vectors representing all possible values of aph ratios
x  = .01:.01:1;   m=numel(x);
y  = (.01:.01:1)';   n=numel(y);

% upper bound of constraint #5
ubrph12 = .38 + ap_i(idx467, :) ./ ap_i(idx412, :);

% upper bound of constraint #7
ubrph45 = 2.3 * ap_i(idx510, :) ./ ap_i(idx555, :) - 1.3;

% calculate speculative solutions and identify all feasible solutions
Axy = NaN(m, n, size(ap_i, 2));
Sxy = NaN(m, n, size(ap_i, 2));
flag = false(m, n, size(ap_i, 2));
aps = ap_i(idx443 | idx412 | idx490 | idx510, :);
for i = progress(1:m)
  xi = x(i);
%   [Ad_xy, Sd_xy] = solver4s_ap(aps, [443 412 490 510], xi, y);
  for j = 1:n
%     calculate Ad and Sd from x and y using a subroutine "solver4s_ap.m"
    yi = y(j);
    [Ad_xy, Sd_xy] = solver4s_ap(aps, [443 412 490 510], xi, yi);
    Axy(i, j, :) = reshape(Ad_xy, 1, 1, size(ap_i, 2));
    Sxy(i, j, :) = reshape(Sd_xy, 1, 1, size(ap_i, 2));
    
%     Ad_xy1 = NaN(1, size(Ad_xy, 2));
%     Sd_xy1 = NaN(1, size(Sd_xy, 2));
%     for k = 1:size(aps, 2)
%       [Ad_xy1(k), Sd_xy1(k)] = solver4s_ap_1(aps(:,k), [443 412 490 510], xi, yi);
%     end
%     if nansum(Ad_xy - Ad_xy1) ~= 0 || any(isnan(Ad_xy1(~isnan(Ad_xy)))) || any(isnan(Ad_xy(~isnan(Ad_xy1))))
%       fprintf('Difference= %i | i= %i | j= %i\n', sum([nansum(Ad_xy1 - Ad_xy1) ...
%         nansum(isnan(Ad_xy1(~isnan(Ad_xy)))) nansum(isnan(Ad_xy(~isnan(Ad_xy1))))]), i, j)
%     end

    if all(isnan(Ad_xy)) || all(isnan(Sd_xy))
      break
    end 
    
    % calculate adg and aph using derived Ad and Sd
    ad_temp = Ad_xy .* exp( -Sd_xy .* lambda) ;
    aph_temp = ap_i - ad_temp ;
    
    % Constraint #4 and #5
    rph12 = aph_temp(idx467,:) ./ aph_temp(idx412,:);
    flag1 = false(size(rph12));
    flag1(all(rph12 < 1.54 & rph12 < ubrph12 & rph12 > .74)) = true;
    if any(flag1)
      break
    end
    
    % Constraint #6 and #7
    rph45 = aph_temp(idx510,:) ./ aph_temp(idx555,:);
    flag2 = false(size(rph45));
    flag2(all(rph45 < 10 & rph45 < ubrph45 & rph45 > 1.3)) = true;
    if any(flag2)
      break
    end
    
    % Constraint #8
    b2r = aph_temp(idx443,:) ./ aph_temp(idx670,:);
    flag3 = false(size(b2r));
    flag3(all(b2r > 1.4 & b2r < 9.1)) = true;
    if any(flag3)
      break
    end
    
    % Constraint #9
    x0  = ap_i(idx412,:) ./ ap_i(idx443,:);
    y0   = ad_temp(idx412,:) ./ ap_i(idx412,:);
    ylb = x0 - .91; 
    yub = x0 - .45;
    flag4 = false(size(x0));
    flag4(all(y0 > yub & y0 < ylb)) = true;
    if any(flag4)
      break
    end
    
    % Constraint #10
    naph  = aph_temp ./ (aph_temp(idx443,:));
    sph45 = (naph(idx510,:) - naph(idx555,:))/(555-510);
    flag5 = false(size(sph45));
    flag5(all(sph45 > 8.7e-3 & sph45 < 3e-3)) = true;
    if any(flag5)
      break
    end
    
    % Constraint #11 #12 #13
    g2rsum  = sum(aph_temp(find(idx550):find(idx630),:)) ./ sum(aph_temp(find(idx650):find(idx700),:)); % Constraint parameter #11
    b2gsum = sum(aph_temp(idx450:find(idx490),:)) ./ sum(aph_temp(find(idx510):idx550,:)); % Constraint parameter #12
    
    % calculate slopes of naph used in constraint #13
    [rmax,irmax] = max(naph(idx450:idx500,:)); % max naph 450-500 nm
    [rmin,irmin] = min(naph(idx500:idx550,:)); % min naph 500-550 nm
    bgslp = (rmax - rmin) ./ (51-irmax+irmin);  % slope of naph to the right-hand side of absorption peak
    [lmax,ilmax] = max(naph(idx430:idx450,:)); % max naph 430-450 nm
    [lmin,ilmin] = min(naph(idx400:idx420,:)); % min naph 400-420 nm
    bvslp = (lmax - lmin)./ (30+ilmax-ilmin);  % slope of naph to the left-hand side of absorption peak
    
    difslp = bgslp - bvslp; % Constraint parameter #13
    
    flag6 = false(size(g2rsum));
    flag6(all(g2rsum < 1.28 & b2gsum > 2 & difslp> -.0047)) = true;
    
    % merge all flags
    flag_merged = false(size(flag1));
    flag_merged(flag1 & flag2 & flag3 & flag4 & flag5 & flag6 & ... 
      all(aph_temp > 0) & all(ad_temp > 0)) = true;
    if any(flag_merged)
      break
    end

%     flag(i, j, reshape(flag_merged, 1, 1, size(ap_i, 2))) = true;
    flag(i, j, :) =  reshape(flag_merged, 1, 1, size(ap_i, 2));
  end
end
% output
Ad = Axy(flag);
Sd = Sxy(flag);
end


function [Ad, Sd] = solver4s_ap(ap, lam, x1, y1)
% 
% This subroutine solves the equation, fres = 0.
% 
% Input:
%     ap, total particulate absorption coefficient at 4 wavelengths
%     lam, wavelengths of ap
%     x1, ratio of aph(lam(2))/aph(lam(1))
%     y1, ratio of aph(lam(4))/aph(lam(3))
% Output:
%     Ad, amplitude of ad,
%     Sd, slope of ad,         where ad = Ad exp(-Sd lambda)
%%
ss = 0.005:.0001:.018; % Constraint #3

% calculate values of the residual function, fres
fres = ((ap(2, :) - x1 .* ap(1, :)) .* (exp(-ss * lam(4)) - y1 * exp(-ss * lam(3)))' - ...
    (exp(-ss * lam(2)) - x1 * exp(-ss * lam(1)))' * (ap(4, :) - y1 .* ap(3, :)))';

% testspec = i;
% fres1 = ( ap(2,testspec) - x1*ap(1,testspec) ) * ( exp(-ss*lam(4)) - y1*exp(-ss*lam(3)) ) - ...
%   ( exp(-ss*lam(2)) - x1*exp(-ss*lam(1)) ) * ( ap(4,testspec) - y1*ap(3,testspec) ); 
% idx = find(diff(fres1>0));
% isempty(idx)
% Sd1 = interp1([fres1(idx) fres1(idx+1)], [ss(idx) ss(idx+1)], 0)
% Ad1 = ( ap(2,testspec) - x1*ap(1,testspec) ) / ( exp(-lam(2)*Sd1) - x1*exp(-lam(1)*Sd1) )

% locate the index of ss after which the "fres" changes sign
foo = diff(fres > 0, [], 2); % e.g. idx = 3, if fres = [1 2 3 -1 -2].

% calculate Sd
Sd = NaN(1, size(ap, 2));
for i = 1:size(foo, 1)
  if sum(foo(i, :) ~= 0) > 0
    try
      Sd(i) = interp1([fres(i, [foo(i, :) ~= 0 false]) fres(i, [false foo(i, :) ~= 0])], ...
        [ss([foo(i, :) ~= 0 false]) ss([false foo(i, :) ~= 0])], 0);
    catch
      fprintf("Warning: solver4s_ap's interpolation crashed for spectrum %i, replaced by NaN\n", i)
    end
  end
end
% calculate Ad from Sd
Ad = (ap(2, :) - x1 * ap(1, :)) ./ (exp(-lam(2) * Sd) - x1 * exp(-lam(1) * Sd)); 
end


function [Ad, Sd] = solver4s_ap_1(ap, lam, x1, y1)
% 
% This subroutine solves the equation, fres = 0.
% 
% Input:
%     ap, total particulate absorption coefficient at 4 wavelengths
%     lam, wavelengths of ap
%     x1, ratio of aph(lam(2))/aph(lam(1))
%     y1, ratio of aph(lam(4))/aph(lam(3))
% Output:
%     Ad, amplitude of ad,
%     Sd, slope of ad,         where ad = Ad exp(-Sd lambda)
% 
% 
    ss   = 0.005:.0001:.018; % Constraint #3
% 
% cacluate values of the residual function, fres
    fres = ( ap(2) - x1*ap(1) ) * ( exp(-ss*lam(4)) - y1*exp(-ss*lam(3)) ) - ...
        ( exp(-ss*lam(2)) - x1*exp(-ss*lam(1)) ) * ( ap(4) - y1*ap(3) ); 
%
% locate the index of ss after which the "fres" changes sign
    idx = find(diff(fres>0)); % e.g. idx = 3, if fres = [1 2 3 -1 -2].
% 
% calcualte the value of ss where fres = 0
    if isempty(idx)
        Sd = NaN;  % no solutions found, if there is no sign change in fres
        Ad = NaN;
    else
        Sd = interp1([fres(idx) fres(idx+1)], [ss(idx) ss(idx+1)], 0); % interpolation between the two points to find where fres=0
        Ad = ( ap(2) - x1*ap(1) ) / ( exp(-lam(2)*Sd) - x1*exp(-lam(1)*Sd) ); % calculate Ad from Sd
    end
end

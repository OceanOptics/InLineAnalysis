function [ad, aph] = partition_ap(ap_i, lambda, prcntl)
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
% 
    m = length(ap_i);
% 
    [fsbl_Ad,fsbl_Sd] = scm_ap(ap_i,lambda); % solve for all feasible solutions of Ad and Sd
    n = length(fsbl_Ad);
    adtemp  = repmat(fsbl_Ad,1,m) .* exp(-fsbl_Sd*lambda'); % all feasible solution of ad
    aphtemp = repmat(ap_i',n,1) - adtemp; % all feasible solution of aph
%     
    ad      = prctile(adtemp,prcntl);  % percentiles of all feasible solutions of ad
    aph     = prctile(aphtemp,prcntl); % percentiles of all feasible solutions of aph
end


function [Ad,Sd] = scm_ap(ap_i, lambda)
%
% Input:  ap_i, spectrum of total particulate absorption coefficient
%         lamda, wavelengths of ap_i
% Output: Ad, amplitude of spectral non-algal particulate absorption coefficient, ad
%         Sd, slope of ad
%         Ad and Sd are vectors representing values for all feasible solutions
% 
% 
% 
% find the index of wavelengths involved in all constraints
    idx400 = find(lambda==400);
    idx412 = find(lambda==412);
    idx420 = find(lambda==420);
    idx430 = find(lambda==430);
    idx443 = find(lambda==443);
    idx450 = find(lambda==450);
    idx467 = find(lambda==467);
    idx490 = find(lambda==490);
    idx500 = find(lambda==500);
    idx510 = find(lambda==510);
    idx550 = find(lambda==550);
    idx555 = find(lambda==555);
    idx630 = find(lambda==630);
    idx650 = find(lambda==650);
    idx670 = lambda==670;
    idx700 = find(lambda==700);
% 
% create vectors representing all possible values of aph ratios
    x  = .01:.01:1;   m=numel(x);
    y  = .01:.01:1;   n=numel(y);
% 
% define empty matrices to store speculative solutions of Ad and Sd derived from x and y
    Axy  = NaN(m,n);
    Sxy  = NaN(m,n);
%
% define flag matrix for index of feasible solutions
    flag = zeros(m,n);
% 
% upper bound of constraint #5
    ubrph12 = .38 + ap_i(idx467)/ap_i(idx412) ;
% 
% upper bound of constraint #7
    ubrph45 = 2.3*ap_i(idx510)/ap_i(idx555) - 1.3 ;
%
% calculate speculative solutions and identify all feasible solutions
for i=1:m
    for j=1:n
        % calculate Ad and Sd from x and y using a subroutine "solver4s_ap.m"
            [Axy(i,j),Sxy(i,j)] = solver4s_ap(ap_i([idx443 idx412 idx490 idx510]),[443 412 490 510],x(i),y(j));
        % calculate adg and aph using derived Ad and Sd
            ad_temp = Axy(i,j)*exp( -Sxy(i,j)*lambda ) ;
            aph_temp = ap_i - ad_temp ;
        % 
        % Constraint #4 and #5
            rph12 = aph_temp(idx467)./aph_temp(idx412) ;
            if rph12<1.54 && rph12<ubrph12 && rph12>.74
              flag1  = 1;
            else
              flag1 = 0;
            end
        % 
        % Constraint #6 and #7
            rph45 = aph_temp(idx510)/aph_temp(idx555); 
            if rph45<10 && rph45<ubrph45 && rph45>1.3
              flag2  = 1;
            else
              flag2 = 0;
            end
        % 
        % Constraint #8
            b2r = aph_temp(idx443)/aph_temp(idx670); 
            if b2r>1.4 && b2r<9.1
              flag3  = 1;
            else
              flag3 = 0;
            end
        % 
        % Constraint #9
            x0  = ap_i(idx412)/ap_i(idx443) ;
            y0   = ad_temp(idx412)/ap_i(idx412) ;
            ylb = x0 - .91; 
            yub = x0 - .45;
            if y0<yub && y0>ylb
              flag4  = 1;
            else
              flag4 = 0;
            end
        % 
        % Constraint #10
            naph  = aph_temp/(aph_temp(idx443));
            sph45 = (naph(idx510)-naph(idx555))/(555-510);
            if sph45<8.7e-3 && sph45>3e-3
              flag5  = 1;
            else
              flag5 = 0;
            end
        %
        % Constraint #11 #12 #13
            g2rsum  = sum(aph_temp(idx550:idx630))/sum(aph_temp(idx650:idx700)); % Constraint parameter #11
            %
            b2gsum = sum(aph_temp(idx450:idx490))/sum(aph_temp(idx510:idx550)); % Constraint parameter #12
            %
            % calculate slopes of naph used in constraint #13
                [rmax,irmax] = max(naph(idx450:idx500)); % max naph 450-500 nm
                [rmin,irmin] = min(naph(idx500:idx550)); % min naph 500-550 nm
                bgslp = (rmax - rmin)/(51-irmax+irmin);  % slope of naph to the right-hand side of absorption peak
                [lmax,ilmax] = max(naph(idx430:idx450)); % max naph 430-450 nm
                [lmin,ilmin] = min(naph(idx400:idx420)); % min naph 400-420 nm
                bvslp = (lmax - lmin)/(30+ilmax-ilmin);  % slope of naph to the left-hand side of absorption peak
            %
            difslp = bgslp - bvslp; % Constraint parameter #13
            %      
            if g2rsum<1.28 && b2gsum>2 && difslp>-.0047
              flag6  = 1;
            else
              flag6 = 0;
            end
        % 
        % Check if all constraints are concurrently satisfied
            if flag1==1 && flag2==1 && flag3==1 && ...
               flag4==1 && flag5==1 && flag6==1 && ... 
               all(aph_temp>0) && all(ad_temp>0)  % non-negative requirement
                flag(i,j) = 1;
            end
    end
end
% 
% output
Ad = Axy(flag==1);
Sd = Sxy(flag==1);
% 
end


function [Ad, Sd] = solver4s_ap(ap,lam,x1,y1)
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
      if size(idx, 2) == 2
        if abs(idx(1) - idx(2)) == 1
          Sd = NaN;  % no solutions found, if there is no sign change in fres
          Ad = NaN;
        else
          Sd = interp1([fres(idx) fres(idx+1)], [ss(idx) ss(idx+1)], 0); % interpolation between the two points to find where fres=0
          Ad = ( ap(2) - x1*ap(1) ) / ( exp(-lam(2)*Sd) - x1*exp(-lam(1)*Sd) ); % calculate Ad from Sd
        end
      else
        Sd = interp1([fres(idx) fres(idx+1)], [ss(idx) ss(idx+1)], 0); % interpolation between the two points to find where fres=0
        Ad = ( ap(2) - x1*ap(1) ) / ( exp(-lam(2)*Sd) - x1*exp(-lam(1)*Sd) ); % calculate Ad from Sd
      end
    end
end
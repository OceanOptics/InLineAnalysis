function [a_sw,bb_sw] = get_water_iops(wave,T,S)

%function to obtain seawater absorption and backscattering spectra

%pure water absorption from Mason et al 2016 for 250-550, Pope and Frye for
%550-730 nm, and Smith and Baker for 730-800 nm
%salt water backscattering from Zhang et al 2009
%corrected for in situ temperature and salinity conditions Sullivan et al. 2006

wl_water1 = [250   260   270   280   290   300   302   304   306   308   310   312   314   316   318   320   322   324   326   328   330 ...
    332   334   336   338   340   342   344   346   348   350   352   354   356   358   360   362   364   366   368   370   372 ...
    374   376   378   380   382   384   386   388   390   392   394   396   398   400   402   404   406   408   410   412   414 ...
    416   418   420   422   424   426   428   430   432   434   436   438   440   442   444   446   448   450   452   454   456 ...
    458   460   462   464   466   468   470   472   474   476   478   480   482   484   486   488   490   492   494   496   498 ...
    500   502   504   506   508   510   512   514   516   518   520   522   524   526   528   530   532   534   536   538   540 ...
    542   544   546   548   550 ];
wl_water2 = [552.5   555 557.5	560	562.5	565	567.5	570	572.5	575	577.5	580	582.5	585	587.5	590	592.5	595	597.5 ...
    600	602.5	605	607.5	610	612.5	615	617.5	620	622.5	625	627.5	630	632.5	635	637.5	640	642.5 ...
    645	647.5	650	652.5	655	657.5	660	662.5	665	667.5	670	672.5	675	677.5	680	682.5	685	687.5 ...
    690	692.5	695	697.5	700	702.5	705	707.5	710	712.5	715	717.5	720	722.5	725	727.5	730  732.5 ...
    735.0  737.5  740.0  742.5  745.0  747.5  750.0  752.5  755.0  757.5  760.0  762.5  765.0  767.5  770.0  772.5 ...
    775.0  777.5 780.0  782.5  785.0  787.5  790.0  792.5  795.0  797.5  800.0];

wl_water = [wl_water1 wl_water2];

aw1 = [58.7100   51.5000   43.5700   22.3000    9.3900    4.6700    4.3300    3.5600    3.1300    2.7500    2.3600    2.0500 ...
    1.8500    1.7600    1.6300    1.4700    1.3300    1.2400    1.1800    1.1200    1.0700    1.0100    0.9900    0.9500 ...
    0.9100    0.8500    0.8200    0.8100    0.8200    0.8400    0.8900    0.9400    0.9700    0.9800    0.9900    1.0600 ...
    1.1500    1.2000    1.2100    1.2200    1.2400    1.2700    1.2900    1.3300    1.3700    1.4300    1.4700    1.5100 ...
    1.5500    1.6200    1.7000    1.7500    1.8500    1.9600    2.0800    2.2200    2.3700    2.4800    2.5700    2.5900 ...
    2.6600    2.7100    2.8000    2.8800    3.0000    3.1200    3.2200    3.3100    3.4400    3.5800    3.7600    3.9500 ...
    4.1700    4.4200    4.8000    5.2200    5.7400    6.2600    6.9100    7.5100    8.0800    8.4200    8.6300    8.7700 ...
    8.9300    9.0900    9.3300    9.5500    9.7900    9.9900   10.3000   10.6500   11.0000   11.3800   11.7700   12.1400 ...
    12.5400   12.9400   13.3600   13.9100   14.6000   15.4500   16.4800   17.7400   19.2600   20.7300   22.4200   24.2400 ...
    26.6800   29.7100   33.0000   35.6900   37.3800   38.2100   38.7800   39.1700   39.6200   40.1700   40.8800   41.6200 ...
    42.4200   43.3000   44.3600   45.4100   46.4500   47.5400   48.8200   50.4000   52.2400   54.2500   56.2900];

aw2 = [0.0593 0.0596 0.0606 0.0619 0.064...
    0.0642 0.0672 0.0695 0.0733 0.0772 0.0836 0.0896 0.0989 0.11 0.122 0.1351 0.1516 0.1672 0.1925 0.2224 0.247 0.2577...
    0.2629 0.2644 0.2665 0.2678 0.2707 0.2755 0.281 0.2834 0.2904 0.2916 0.2995 0.3012 0.3077 0.3108 0.322 0.325 0.335...
    0.34 0.358 0.371 0.393 0.41 0.424 0.429 0.436 0.439 0.448 0.448 0.461 0.465 0.478 0.486 0.502 0.516 0.538 0.559...
    0.592 0.624	0.663 0.704 0.756 0.827 0.914 1.007 1.119 1.231 1.356 1.489 1.678 1.7845 1.9333 2.0822 2.2311 2.3800...
    2.4025 2.4250 2.4475 2.4700 2.4900 2.5100 2.5300 2.5500 2.5400 2.5300 2.5200 2.5100 2.4725 2.4350 2.3975 2.3600...
    2.3100 2.2600 2.2100 2.1600 2.1375 2.1150 2.0925 2.0700];

a_water = [aw1*1e-3 aw2]; % 1e-3 comes from Mason et al. 2016, Table 2

a_pw = interp1(wl_water,a_water,wave,'linear','extrap');

%temp and salinity correction for water absorption (need to know at what T it was measured):
S(isnan(S)) = 35;
T(isnan(T)) = 22;
T_pope=22.0;

% % use salt water scattering from Zhang et al 2009
% [betasw124, bb_sw, beta90sw, theta] = betasw124_ZHH2009(wave, S, T);

% use salt water scattering from:
%   L. Hu, X. Zhang, and M. J. Perry, "Light scattering by pure seawater:
%   Effect of pressure," Deep-Sea Res. I 146, 103-109 (2019).
[~, bb_sw] = betasw_HZP2019(wave,90,T,S,1);

% use the temp and salinity corrections from Sullivan et al. 2006
[psiT,psiS] = tempsal_corr(wave);
%temperature and salinity corrections:
a_sw = a_pw(:)'+psiT(:)'.*(T-T_pope)+psiS(:)'.*S;

end

function [psiT,psiS] = tempsal_corr(wavel)
tmp = [400.0000    0.0001         0
    402.0000    0.0001         0
    404.0000    0.0001         0
    406.0000    0.0001         0
    408.0000         0    0.0000
    410.0000         0    0.0000
    412.0000         0    0.0000
    414.0000    0.0001    0.0000
    416.0000    0.0001    0.0000
    418.0000         0    0.0000
    420.0000         0    0.0000
    422.0000         0         0
    424.0000         0         0
    426.0000         0         0
    428.0000         0         0
    430.0000         0         0
    432.0000         0         0
    434.0000         0         0
    436.0000         0         0
    438.0000         0         0
    440.0000         0         0
    442.0000         0         0
    444.0000         0         0
    446.0000         0         0
    448.0000         0         0
    450.0000         0         0
    452.0000         0         0
    454.0000         0         0
    456.0000         0         0
    458.0000         0         0
    460.0000         0         0
    462.0000         0         0
    464.0000         0         0
    466.0000         0         0
    468.0000         0         0
    470.0000         0         0
    472.0000         0         0
    474.0000         0         0
    476.0000         0         0
    478.0000         0         0
    480.0000         0   -0.0000
    482.0000         0   -0.0000
    484.0000         0   -0.0000
    486.0000         0   -0.0000
    488.0000         0   -0.0000
    490.0000         0   -0.0000
    492.0000         0   -0.0000
    494.0000         0   -0.0000
    496.0000         0   -0.0000
    498.0000         0   -0.0000
    500.0000         0   -0.0000
    502.0000         0   -0.0000
    504.0000         0   -0.0000
    506.0000         0   -0.0000
    508.0000    0.0001   -0.0000
    510.0000    0.0001   -0.0000
    512.0000    0.0001   -0.0000
    514.0000    0.0001   -0.0000
    516.0000    0.0001         0
    518.0000    0.0001         0
    520.0000    0.0001         0
    522.0000    0.0001         0
    524.0000    0.0001         0
    526.0000    0.0001         0
    528.0000         0         0
    530.0000         0         0
    532.0000         0         0
    534.0000         0         0
    536.0000         0         0
    538.0000         0         0
    540.0000         0         0
    542.0000         0   -0.0000
    544.0000         0   -0.0000
    546.0000         0         0
    548.0000         0         0
    550.0000         0         0
    552.0000         0         0
    554.0000         0         0
    556.0000         0         0
    558.0000         0         0
    560.0000         0   -0.0000
    562.0000         0   -0.0000
    564.0000         0   -0.0000
    566.0000         0   -0.0000
    568.0000         0   -0.0000
    570.0000         0   -0.0000
    572.0000         0   -0.0000
    574.0000         0   -0.0000
    576.0000    0.0001   -0.0000
    578.0000    0.0001   -0.0000
    580.0000    0.0002   -0.0000
    582.0000    0.0002   -0.0000
    584.0000    0.0003   -0.0000
    586.0000    0.0004   -0.0000
    588.0000    0.0004   -0.0000
    590.0000    0.0005   -0.0000
    592.0000    0.0006   -0.0000
    594.0000    0.0008   -0.0000
    596.0000    0.0009   -0.0000
    598.0000    0.0010   -0.0000
    600.0000    0.0011   -0.0000
    602.0000    0.0011         0
    604.0000    0.0011    0.0000
    606.0000    0.0011    0.0000
    608.0000    0.0011    0.0000
    610.0000    0.0010    0.0000
    612.0000    0.0009    0.0001
    614.0000    0.0008    0.0001
    616.0000    0.0007    0.0001
    618.0000    0.0007    0.0001
    620.0000    0.0006    0.0001
    622.0000    0.0005    0.0001
    624.0000    0.0004    0.0001
    626.0000    0.0003    0.0001
    628.0000    0.0002    0.0001
    630.0000    0.0002    0.0001
    632.0000    0.0001    0.0001
    634.0000         0    0.0001
    636.0000         0    0.0000
    638.0000         0    0.0000
    640.0000   -0.0001    0.0000
    642.0000   -0.0001    0.0000
    644.0000   -0.0001    0.0000
    646.0000   -0.0001    0.0000
    648.0000   -0.0001    0.0000
    650.0000         0    0.0000
    652.0000         0    0.0000
    654.0000    0.0001    0.0000
    656.0000    0.0001    0.0000
    658.0000    0.0002    0.0000
    660.0000    0.0002    0.0000
    662.0000    0.0002    0.0000
    664.0000    0.0002    0.0000
    666.0000    0.0002    0.0000
    668.0000    0.0001    0.0000
    670.0000    0.0001    0.0000
    672.0000         0    0.0000
    674.0000         0    0.0000
    676.0000   -0.0001         0
    678.0000   -0.0001   -0.0000
    680.0000   -0.0001   -0.0000
    682.0000   -0.0002   -0.0000
    684.0000   -0.0002   -0.0000
    686.0000   -0.0002   -0.0001
    688.0000   -0.0001   -0.0001
    690.0000   -0.0001   -0.0001
    692.0000   -0.0001   -0.0001
    694.0000         0   -0.0001
    696.0000    0.0001   -0.0001
    698.0000    0.0002   -0.0001
    700.0000    0.0004   -0.0002
    702.0000    0.0006   -0.0002
    704.0000    0.0009   -0.0002
    706.0000    0.0012   -0.0002
    708.0000    0.0017   -0.0002
    710.0000    0.0022   -0.0002
    712.0000    0.0027   -0.0002
    714.0000    0.0033   -0.0002
    716.0000    0.0040   -0.0002
    718.0000    0.0049   -0.0002
    720.0000    0.0060   -0.0002
    722.0000    0.0071   -0.0003
    724.0000    0.0084   -0.0003
    726.0000    0.0096   -0.0003
    728.0000    0.0108   -0.0002
    730.0000    0.0120   -0.0002
    732.0000    0.0130   -0.0002
    734.0000    0.0139   -0.0001
    736.0000    0.0146    0.0000
    738.0000    0.0150    0.0001
    740.0000    0.0150    0.0003
    742.0000    0.0147    0.0004
    744.0000    0.0142    0.0005
    746.0000    0.0138    0.0007
    748.0000    0.0131    0.0008
    750.0000    0.0124    0.0009];

psiT = interp1(tmp(:,1),tmp(:,2), wavel,'linear',0); % both this line and the next changed from 'spline' on 14 Dec 2015
psiS = interp1(tmp(:,1),tmp(:,3), wavel,'linear',0);

end




% function [betasw124, bsw, beta90sw, theta] = betasw124_ZHH2009(lambda, S, Tc, delta)
% %
% % function [betasw124, bsw, beta90sw, theta]= betasw124_ZHH2009(lambda,S,Tc,delta)
% %
% %
% % Scattering by pure seawater: Effect of salinity
% % Xiaodong Zhang, Lianbo Hu, and Ming-Xia He, Optics Express, 2009, accepted
% % lambda (nm): wavelength
% % Tc: temperauter in degree Celsius, must be a scalar
% % S: salinity, must be scalar
% % delta: depolarization ratio, if not provided, default = 0.039 will be
% % used.
% % betasw: volume scattering at angles defined by theta. Its size is [x y],
% % where x is the number of angles (x = length(theta)) and y is the number
% % of wavelengths in lambda (y = length(lambda))
% % beta90sw: volume scattering at 90 degree. Its size is [1 y]
% % bw: total scattering coefficient. Its size is [1 y]
% % for backscattering coefficients, divide total scattering by 2
% %
% % Xiaodong Zhang, March 10, 2009
% %
% %
% % MODIFIED on 17/05/2011 to be able to process bbp profiles with coincident T and sal profiles
% % MODIFIED on 05 Apr 2013 to use 124 degs instead of 117 degs
% %
% 
% % make sure that S and Tc are column vectors
% S = S(:);
% Tc = Tc(:);
% 
% % values of the constants
% Na = 6.0221417930e23 ;   %  Avogadro's constant
% Kbz = 1.3806503e-23 ;    %  Boltzmann constant
% Tk = Tc+273.15 ;         %  Absolute tempearture
% M0 = 18e-3;              %  Molecular weigth of water in kg/mol
% 
% %error(nargchk(3, 4, nargin));
% if nargin == 3
%     delta = 0.039; % Farinato and Roswell (1976)
% end
% 
% theta=[0:.01:180];
% 
% %if ~isscalar(Tc) || ~isscalar (S)
% %    error('Both Tc and S need to be scalar variable');
% %end
% 
% lambda = lambda(:)'; % a row variable
% rad = theta(:)*pi/180; % angle in radian as a colum variable
% 
% % nsw: absolute refractive index of seawater
% % dnds: partial derivative of seawater refractive index w.r.t. salinity
% [nsw dnds] = RInw(lambda,Tc,S);
% 
% % isothermal compressibility is from Lepple & Millero (1971,Deep
% % Sea-Research), pages 10-11
% % The error ~ +/-0.004e-6 bar^-1
% IsoComp = BetaT(Tc,S);
% 
% % density of water and seawater,unit is Kg/m^3, from UNESCO,38,1981
% density_sw = rhou_sw(Tc, S);
% 
% % water activity data of seawater is from Millero and Leung (1976,American
% % Journal of Science,276,1035-1077). Table 19 was reproduced using
% % Eq.(14,22,23,88,107) then were fitted to polynominal equation.
% % dlnawds is partial derivative of natural logarithm of water activity
% % w.r.t.salinity
% dlnawds = dlnasw_ds(Tc, S);
% 
% % density derivative of refractive index from PMH model
% DFRI = PMH(nsw);  %% PMH model
% 
% % volume scattering at 90 degree due to the density fluctuation
% beta_df = pi*pi./2.*((lambda*1e-9).^(-4)).*Kbz.*Tk.*IsoComp.*DFRI.^2.*(6+6*delta)./(6-7*delta);
% 
% % volume scattering at 90 degree due to the concentration fluctuation
% flu_con = S.*M0.*dnds.^2./density_sw./(-dlnawds)./Na;
% beta_cf = 2*pi*pi*((lambda*1e-9).^(-4)).*nsw.^2.*(flu_con)*(6+6*delta)/(6-7*delta);
% 
% % total volume scattering at 90 degree
% beta90sw = beta_df+beta_cf;
% bsw = 8*pi/3.*beta90sw.*(2+delta)./(1+delta);
% 
% rad124 = find (rad2deg(rad)>=124, 1);
% 
% betasw124 = NaN(size(beta90sw));
% %for i=1:length(beta90sw)
% %    betasw124(i,:) = beta90sw(i) * (  1+((cos(rad(rad124))).^2).*(1-delta)./(1+delta)  );
% %end
% betasw124 = beta90sw * (  1+((cos(rad(rad124))).^2).*(1-delta)./(1+delta)  );
% 
% end

% function [nsw dnswds]= RInw(lambda,Tc,S)
% % refractive index of air is from Ciddor (1996,Applied Optics)
% n_air = 1.0+(5792105.0./(238.0185-1./(lambda/1e3).^2)+167917.0./(57.362-1./(lambda/1e3).^2))/1e8;
% 
% % refractive index of seawater is from Quan and Fry (1994, Applied Optics)
% n0 = 1.31405; n1 = 1.779e-4 ; n2 = -1.05e-6 ; n3 = 1.6e-8 ; n4 = -2.02e-6 ;
% n5 = 15.868; n6 = 0.01155;  n7 = -0.00423;  n8 = -4382 ; n9 = 1.1455e6;
% 
% nsw = n0+(n1+n2*Tc+n3*Tc.^2).*S+n4*Tc.^2+(n5+n6*S+n7*Tc)./lambda+n8./lambda.^2+n9./lambda.^3; % pure seawater
% nsw = nsw.*n_air;
% dnswds = (n1+n2*Tc+n3*Tc.^2+n6./lambda).*n_air;
% 
% end

% function IsoComp = BetaT(Tc, S)
% % pure water secant bulk Millero (1980, Deep-sea Research)
% kw = 19652.21+148.4206*Tc-2.327105*Tc.^2+1.360477e-2*Tc.^3-5.155288e-5*Tc.^4;
% Btw_cal = 1./kw;
% 
% % isothermal compressibility from Kell sound measurement in pure water
% % Btw = (50.88630+0.717582*Tc+0.7819867e-3*Tc.^2+31.62214e-6*Tc.^3-0.1323594e-6*Tc.^4+0.634575e-9*Tc.^5)./(1+21.65928e-3*Tc)*1e-6;
% 
% % seawater secant bulk
% a0 = 54.6746-0.603459*Tc+1.09987e-2*Tc.^2-6.167e-5*Tc.^3;
% b0 = 7.944e-2+1.6483e-2*Tc-5.3009e-4*Tc.^2;
% 
% Ks =kw + a0.*S + b0.*S.^1.5;
% 
% % calculate seawater isothermal compressibility from the secant bulk
% IsoComp = 1./Ks*1e-5; % unit is pa
% end

% function density_sw = rhou_sw(Tc, S)
% 
% % density of water and seawater,unit is Kg/m^3, from UNESCO,38,1981
% a0 = 8.24493e-1;  a1 = -4.0899e-3; a2 = 7.6438e-5; a3 = -8.2467e-7; a4 = 5.3875e-9;
% a5 = -5.72466e-3; a6 = 1.0227e-4;  a7 = -1.6546e-6; a8 = 4.8314e-4;
% b0 = 999.842594; b1 = 6.793952e-2; b2 = -9.09529e-3; b3 = 1.001685e-4;
% b4 = -1.120083e-6; b5 = 6.536332e-9;
% 
% % density for pure water
% density_w = b0 + b1.*Tc + b2.*Tc.^2 + b3.*Tc.^3 + b4.*Tc.^4 + b5.*Tc.^5;
% % density for pure seawater
% density_sw = density_w +((a0 + a1.*Tc + a2.*Tc.^2 + a3.*Tc.^3 + a4.*Tc.^4).*S + (a5 + a6.*Tc + a7.*Tc.^2).*S.^1.5 + a8.*S.^2);
% 
% end

% function dlnawds = dlnasw_ds(Tc, S)
% % water activity data of seawater is from Millero and Leung (1976,American
% % Journal of Science,276,1035-1077). Table 19 was reproduced using
% % Eqs.(14,22,23,88,107) then were fitted to polynominal equation.
% % dlnawds is partial derivative of natural logarithm of water activity
% % w.r.t.salinity
% % lnaw = (-1.64555e-6-1.34779e-7*Tc+1.85392e-9*Tc.^2-1.40702e-11*Tc.^3)+......
% %            (-5.58651e-4+2.40452e-7*Tc-3.12165e-9*Tc.^2+2.40808e-11*Tc.^3).*S+......
% %            (1.79613e-5-9.9422e-8*Tc+2.08919e-9*Tc.^2-1.39872e-11*Tc.^3).*S.^1.5+......
% %            (-2.31065e-6-1.37674e-9*Tc-1.93316e-11*Tc.^2).*S.^2;
% 
% dlnawds = (-5.58651e-4+2.40452e-7*Tc-3.12165e-9*Tc.^2+2.40808e-11*Tc.^3)+...
%     1.5*(1.79613e-5-9.9422e-8*Tc+2.08919e-9*Tc.^2-1.39872e-11*Tc.^3).*S.^0.5+...
%     2*(-2.31065e-6-1.37674e-9*Tc-1.93316e-11*Tc.^2).*S;
% 
% % density derivative of refractive index from PMH model
% 
% end

% function n_density_derivative=PMH(n_wat)
% n_wat2 = n_wat.^2;
% n_density_derivative=(n_wat2-1).*(1+2/3*(n_wat2+2).*(n_wat/3-1/3./n_wat).^2);
% end





function [betasw,bsw]= betasw_HZP2019(lambda,theta,Tc,SA,p,delta)
% Scattering by pure seawater for all conditions: high salinity (e.g., Dead
% Sea) and high pressure (e.g., at depths > 1000 m)
%
% [betasw,bsw]= betasw_HZP2019(lambda,theta,Tc,SA,p,delta)
% lambda [nm]: wavelength; scalar or vector
% theta [deg]: scattering angle; scalar or vector
% Tc [deg C]: temperature in Celsius, scalar or vector
% SA [g kg^-1]: absolute salinity, SA = S*1.0047, where S is PSU. The
%   reason using absolute salinity is that PSU is not defined for S > 42.
%   Either scalar or vector
% p [dbar]: sea (gauge) pressure. Default = 0, which is under normal
%   condition at sea surface; scalar or vector
% delta: depolarization ratio. Default = 0.039. Must be a scalar
%
% Tc, SA, and p are either scalars or vectors of same size. If one or two
% of them are scalars and the other are vectors, all of them will be made
% into vectors. The idea is to represent a depth profile.
%
% betasw: spectral volume scattering at angle defined by theta. Its size is
%   [x y z], where x = length(Tc, SA or p), y = length(lambda)), and z =
%   length(theta)
% bw: total scattering coefficient. Its size is [x y]
%     for backscattering coefficients, divide total scattering by 2
%
% Revision history
% October 30, 2018
%   First version: only Tc, SA, and p are vectors
% December 20, 2018
%   Second version: allow lambda and theta to be vectors
% April 17, 2019
%   Third version: release to public
%
% Author: Xiaodong Zhang
%   Division of Marine Science
%   School of Ocean Science and Engineering
%   University of Southern Mississippi
%   xiaodong.zhang@ums.edu
%
% To cite:
%   L. Hu, X. Zhang, and M. J. Perry, "Light scattering by pure seawater:
%   Effect of pressure," Deep-Sea Res. I 146, 103-109 (2019).

if nargin == 4
    p = 0;
    delta = 0.039;
elseif nargin == 5
    delta = 0.039;
elseif nargin == 6
    % do nothing
else
    help betasw_all;
    error('The number of input argement is wrong');
end

% Make sure lambda, theta, Tc, SA and p are all vectors (or scalars)
if ~isvector(lambda) || ~isvector(theta) || ...
        ~isvector(Tc) || ~isvector(SA) || ~isvector(p)
    error('Wavelength, angel, temperature, salinity and pressure should be vector!');
end

% values of the Constants
Kbz = 1.3806503e-23;    %  Boltzmann constant
sal_con = 35.16504/35; % conversion, SA = sal_con x PSU

% Normal call is for theta to be a scalar
% If theta is a vector, we will first calculate beta90 and then calculate
% the VSF at other angles.
if numel(theta)>1
    [betasw90,bsw] = betasw_HZP2019(lambda,90,Tc,SA,p,delta);
    betasw = repmat(betasw90,1,1,numel(theta));
    for i=1:numel(theta)
        betasw(:,:,i) = ang_from_90(betasw90,theta(i),delta);
    end
    betasw = squeeze(betasw);
    return
else
    % do nothing
end

% Make sure the sizes of Tc, SA, and p are consistent
num = max([numel(Tc),numel(SA),numel(p)]);
ctd = zeros(num,3);
ctd(:,1) = Tc(:);
ctd(:,2) = SA(:);
ctd(:,3) = p(:);

Tc = ctd(:,1);
SA = ctd(:,2);
p = ctd(:,3);

lambda = lambda(:)'; % make it a row vector
% By now, Tc, SA and p are column vectors of the same length
Tk = Tc+273.15;         %  Absolute tempearture

% nsw: absolute refractive index of seawater
S = SA/sal_con;
[nsw, dnswds] = RIsw_MS(lambda/1000,Tc,S,p) ; % Millard and Seaver 1990
% Note that ds in dnswds from RIsw_MS is in unit of PSU whereas we need ds
% to be in absolute salinity in units of g/kg.
dnswds = dnswds/sal_con;

% density derivative of refractive index from PMH model
DFRI = PMH(nsw);  %% PMH model

% salinity derivative of refractive index. Note the salinity here is the
% absolute salinity with units of g/kg.
SFRI = 2*nsw.*dnswds;

% partial derivative of g w.r.t p
gp = gsw_gibbs(0,0,1,SA,Tc,p);

% partial derivative of g w.r.t p and p
gpp = gsw_gibbs(0,0,2,SA,Tc,p);

% partial derivative of g w.r.t. S and S
gss = gsw_gibbs(2,0,0,SA,Tc,p);

% total volume scattering at 90 degree
beta90 =  pi*pi*Kbz/2*((lambda*1e-9).^(-4)).*Tk.*...
    (DFRI.^2.*(-gpp./gp)+SFRI.^2.*(gp./gss))*(6+6*delta)/(6-7*delta);

% take care of 0 salinity, for which gss is not defined
ind = ~(SA>0);
if any(ind) % 
    beta90(ind,:) =  pi*pi*Kbz/2*((lambda*1e-9).^(-4)).*Tk(ind).*...
        (DFRI(ind,:).^2.*(-gpp(ind)./gp(ind)))*(6+6*delta)/(6-7*delta);
end
% volume scattering angle distribution
betasw = ang_from_90(beta90,theta,delta);

% total scattering coefficient
bsw=8*pi/3*beta90*(2+delta)/(1+delta);
end

function betasw = ang_from_90(beta90,theta,delta)
betasw = beta90*(1+(1-delta)/(1+delta)*cosd(theta)^2);
end

function [nsw, dnswds]= RIsw_MS(lambda,Tc,S,P)
% refractive index of seawater is from Millard and Seaver (1990, Deep Sea Research)
% lambda: wavelength in um (row vector)
% Tc: temperature in degrees Celsius (column vector)
% S: salinity in PSU (column vector)
% P: pressure in dbar (column vector)

A0 = 1.3280657;L2 = -0.0045536802; LM2 = 0.0025471707; LM4 = 7.501966E-6;
LM6 = 2.802632E-6; T1 = -5.2883907E-6; T2 = -3.0738272E-6;
T3 = 3.0124687E-8; T4 = -2.0883178E-10; TL = 1.0508621E-5;
T2L = 2.1282248E-7; T3L = -1.705881E-9;

S0 = 1.9029121E-4; S1LM2 = 2.4239607E-6; S1T = -7.3960297E-7;
S1T2 = 8.9818478E-9; S1T3 = 1.2078804E-10; STL = -3.589495E-7;

P1 = 1.5868383E-6; P2 = -1.574074E-11; PLM2 = 1.0712063E-8;
PT = -9.4834486E-9; PT2 = 1.0100326E-10; P2T2 = 5.8085198E-15;

P1S = -1.1177517E-9; PTS = 5.7311268E-11; PT2S = -1.5460458E-12;

n1 = A0+L2*lambda.^2+LM2./lambda.^2+LM4./lambda.^4+LM6./lambda.^6+...
    T1*Tc+T2*Tc.^2+T3*Tc.^3+T4*Tc.^4+...
    TL*Tc.*lambda+T2L*Tc.^2.*lambda+T3L*Tc.^3.*lambda;
n2 = S0*S+S1LM2*S./lambda.^2+S1T*S.*Tc+S1T2*S.*Tc.^2+...
    S1T3*S.*Tc.^3+STL*S.*Tc.*lambda;
n3 = P1*P+P2*P.^2+PLM2*P./lambda.^2+PT*P.*Tc+...
    PT2*P.*Tc.^2+P2T2*P.^2.*Tc.^2;
n4 = P1S*P.*S+PTS*P.*S.*Tc+PT2S*P.*Tc.^2.*S;

nsw = n1+n2+n3+n4;
dnswds = S0+S1LM2./lambda.^2+S1T*Tc+S1T2*Tc.^2+S1T3*Tc.^3+STL*Tc.*lambda+...
    P.*(P1S+PTS*Tc+PT2S*Tc.^2);
n_air = RIair_Ciddor(lambda);
nsw = nsw.*n_air; % absolute refractive index
dnswds = dnswds.*n_air;
end

function n_air = RIair_Ciddor(lambda)
% refractive index of standard air (T=15C,P=1013.25mb) is from Ciddor (1996,Applied Optics)
% lambda: wavelength in micron
n_air = 1.0+(5792105.0./(238.0185-1./lambda.^2)+...
    167917.0./(57.362-1./lambda.^2))/1e8;
end

function nsw = RInw(lambda,Tc,S)
% this function is used to calculate the refractive index of seawater 
% from Quan and Fry (1994, Applied Optics)
% lambda: wavelength (nm)
%     Tc: temperature in degree
%      S: salinity (pss)
%    nsw: refractive index of seawater

% refractive index of air is from Ciddor (1996,Applied Optics)
n_air = 1.0+(5792105.0./(238.0185-1./(lambda/1e3).^2)+...
    167917.0./(57.362-1./(lambda/1e3).^2))/1e8;

n0 = 1.31405; n1 = 1.779e-4 ; n2 = -1.05e-6 ; n3 = 1.6e-8 ; n4 = -2.02e-6 ;
n5 = 15.868; n6 = 0.01155;  n7 = -0.00423;  n8 = -4382 ; n9 = 1.1455e6;

nsw = n0+(n1+n2*Tc+n3*Tc.^2).*S+n4*Tc.^2+(n5+n6*S+n7*Tc)./lambda+...
    n8./lambda.^2+n9./lambda.^3; 
nsw = nsw.*n_air;
end

function n_density_derivative=PMH(nw)
% this function is used to calculate the density fluctuation of the
% refractive index (DFRI) based on Proutiere et al (1992,J. Phys. Chem.)
% nw: refractive index 

nw2 = nw.^2;
n_density_derivative=(nw2-1).*(1+2/3*(nw2+2).*(nw/3-1/3./nw).^2);
end

function gibbs = gsw_gibbs(ns,nt,np,SA,t,p)

% gsw_gibbs                                Gibbs energy and its derivatives
%==========================================================================
%
% USAGE:
%   gibbs = gsw_gibbs(ns,nt,np,SA,t,p)
%
% DESCRIPTION:
%  Calculates specific Gibbs energy and its derivatives up to order 3 for
%  seawater.  The Gibbs function for seawater is that of TEOS-10
%  (IOC et al., 2010), being the sum of IAPWS-08 for the saline part and
%  IAPWS-09 for the pure water part.  These IAPWS releases are the
%  officially blessed IAPWS descriptions of Feistel (2008) and the pure
%  water part of Feistel (2003).  Absolute Salinity, SA, in all of the GSW
%  routines is expressed on the Reference-Composition Salinity Scale of
%  2008 (RCSS-08) of Millero et al. (2008).
%
% INPUT:
%  ns  =  order of SA derivative                     [ integers 0, 1 or 2 ]
%  nt  =  order of t derivative                      [ integers 0, 1 or 2 ]
%  np  =  order of p derivative                      [ integers 0, 1 or 2 ]
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         (ie. absolute pressure - 10.1325 dbar)
%
%  SA, t and p need to have the same dimensions.
%
% OUTPUT:
%  gibbs  =  Specific Gibbs energy or its derivatives.
%
%            The Gibbs energy (when ns = nt = np = 0) has units of:
%                                                                  [ J/kg ]
%            The Absolute Salinity derivatives are output in units of:
%                                                   [ (J/kg) (g/kg)^(-ns) ]
%            The temperature derivatives are output in units of:
%                                                      [ (J/kg) (K)^(-nt) ]
%            The pressure derivatives are output in units of:
%                                                     [ (J/kg) (Pa)^(-np) ]
%            The mixed derivatives are output in units of:
%                              [ (J/kg) (g/kg)^(-ns) (K)^(-nt) (Pa)^(-np) ]
%  Note. The derivatives are taken with respect to pressure in Pa, not
%    withstanding that the pressure input into this routine is in dbar.
%
% AUTHOR:
%  David Jackett, Paul Barker and Trevor McDougall     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05 (27th January 2015)
%
% REFERENCES:
%  Feistel, R., 2003: A new extended Gibbs thermodynamic potential of
%   seawater,  Progr. Oceanogr., 58, 43-114.
%
%  Feistel, R., 2008: A Gibbs function for seawater thermodynamics
%   for -6 to 80°C and salinity up to 120 g kg1, Deep-Sea Res. I,
%   55, 1639-1671.
%
%  IAPWS, 2008: Release on the IAPWS Formulation 2008 for the
%   Thermodynamic Properties of Seawater. The International Association
%   for the Properties of Water and Steam. Berlin, Germany, September
%   2008, available from http://www.iapws.org.  This Release is referred
%   to as IAPWS-08.
%
%  IAPWS, 2009: Supplementary Release on a Computationally Efficient
%   Thermodynamic Formulation for Liquid Water for Oceanographic Use.
%   The International Association for the Properties of Water and Steam.
%   Doorwerth, The Netherlands, September 2009, available from
%   http://www.iapws.org.  This Release is referred to as IAPWS-09.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 2.6 and appendices A.6,  G and H of this TEOS-10 Manual.
%
%  Millero, F.J., R. Feistel, D.G. Wright, and T.J. McDougall, 2008:
%   The composition of Standard Seawater and the definition of the
%   Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72.
%
%  Reference page in Help browser
%       <a href="matlab:doc gsw_gibbs">doc gsw_gibbs</a>
%  Note that this reference page includes the code contained in gsw_gibbs.
%  We have opted to encode this programme as it is a global standard and 
%  such we cannot allow anyone to change it.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

% Set the upper and lower limits where the Gibbs function is defined.
SA(SA > 120 | t < -12 | t > 80 | p > 12000) = NaN;
t(SA > 120 | t < -12 | t > 80 | p > 12000) = NaN;
p(SA > 120 | t < -12 | t > 80 | p > 12000) = NaN;

sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).

x2 = sfac.*SA;
x = sqrt(x2);
y = t.*0.025;
z = p.*1e-4; %Note.The input pressure (p) is sea pressure in units of dbar.

if ns==0 & nt==0 & np==0
    
    g03 = 101.342743139674 + z.*(100015.695367145 + ...
        z.*(-2544.5765420363 + z.*(284.517778446287 + ...
        z.*(-33.3146754253611 + (4.20263108803084 - 0.546428511471039.*z).*z)))) + ...
        y.*(5.90578347909402 + z.*(-270.983805184062 + ...
        z.*(776.153611613101 + z.*(-196.51255088122 + (28.9796526294175 - 2.13290083518327.*z).*z))) + ...
        y.*(-12357.785933039 + z.*(1455.0364540468 + ...
        z.*(-756.558385769359 + z.*(273.479662323528 + z.*(-55.5604063817218 + 4.34420671917197.*z)))) + ...
        y.*(736.741204151612 + z.*(-672.50778314507 + ...
        z.*(499.360390819152 + z.*(-239.545330654412 + (48.8012518593872 - 1.66307106208905.*z).*z))) + ...
        y.*(-148.185936433658 + z.*(397.968445406972 + ...
        z.*(-301.815380621876 + (152.196371733841 - 26.3748377232802.*z).*z)) + ...
        y.*(58.0259125842571 + z.*(-194.618310617595 + ...
        z.*(120.520654902025 + z.*(-55.2723052340152 + 6.48190668077221.*z))) + ...
        y.*(-18.9843846514172 + y.*(3.05081646487967 - 9.63108119393062.*z) + ...
        z.*(63.5113936641785 + z.*(-22.2897317140459 + 8.17060541818112.*z))))))));
    
    g08 = x2.*(1416.27648484197 + z.*(-3310.49154044839 + ...
        z.*(384.794152978599 + z.*(-96.5324320107458 + (15.8408172766824 - 2.62480156590992.*z).*z))) + ...
        x.*(-2432.14662381794 + x.*(2025.80115603697 + ...
        y.*(543.835333000098 + y.*(-68.5572509204491 + ...
        y.*(49.3667694856254 + y.*(-17.1397577419788 + 2.49697009569508.*y))) - 22.6683558512829.*z) + ...
        x.*(-1091.66841042967 - 196.028306689776.*y + ...
        x.*(374.60123787784 - 48.5891069025409.*x + 36.7571622995805.*y) + 36.0284195611086.*z) + ...
        z.*(-54.7919133532887 + (-4.08193978912261 - 30.1755111971161.*z).*z)) + ...
        z.*(199.459603073901 + z.*(-52.2940909281335 + (68.0444942726459 - 3.41251932441282.*z).*z)) + ...
        y.*(-493.407510141682 + z.*(-175.292041186547 + (83.1923927801819 - 29.483064349429.*z).*z) + ...
        y.*(-43.0664675978042 + z.*(383.058066002476 + z.*(-54.1917262517112 + 25.6398487389914.*z)) + ...
        y.*(-10.0227370861875 - 460.319931801257.*z + y.*(0.875600661808945 + 234.565187611355.*z))))) + ...
        y.*(168.072408311545 + z.*(729.116529735046 + ...
        z.*(-343.956902961561 + z.*(124.687671116248 + z.*(-31.656964386073 + 7.04658803315449.*z)))) + ...
        y.*(880.031352997204 + y.*(-225.267649263401 + ...
        y.*(91.4260447751259 + y.*(-21.6603240875311 + 2.13016970847183.*y) + ...
        z.*(-297.728741987187 + (74.726141138756 - 36.4872919001588.*z).*z)) + ...
        z.*(694.244814133268 + z.*(-204.889641964903 + (113.561697840594 - 11.1282734326413.*z).*z))) + ...
        z.*(-860.764303783977 + z.*(337.409530269367 + ...
        z.*(-178.314556207638 + (44.2040358308 - 7.92001547211682.*z).*z))))));
    
    g08(x>0) = g08(x>0) + x2(x>0).*(5812.81456626732 + 851.226734946706.*y(x>0)).*log(x(x>0));
    
    gibbs = g03 + g08;
    
elseif ns==1 & nt==0 & np==0
    
    g08 = 8645.36753595126 + z.*(-6620.98308089678 + ...
        z.*(769.588305957198 + z.*(-193.0648640214916 + (31.6816345533648 - 5.24960313181984.*z).*z))) + ...
        x.*(-7296.43987145382 + x.*(8103.20462414788 + ...
        y.*(2175.341332000392 + y.*(-274.2290036817964 + ...
        y.*(197.4670779425016 + y.*(-68.5590309679152 + 9.98788038278032.*y))) - 90.6734234051316.*z) + ...
        x.*(-5458.34205214835 - 980.14153344888.*y + ...
        x.*(2247.60742726704 - 340.1237483177863.*x + 220.542973797483.*y) + 180.142097805543.*z) + ...
        z.*(-219.1676534131548 + (-16.32775915649044 - 120.7020447884644.*z).*z)) + ...
        z.*(598.378809221703 + z.*(-156.8822727844005 + (204.1334828179377 - 10.23755797323846.*z).*z)) + ...
        y.*(-1480.222530425046 + z.*(-525.876123559641 + (249.57717834054571 - 88.449193048287.*z).*z) + ...
        y.*(-129.1994027934126 + z.*(1149.174198007428 + z.*(-162.5751787551336 + 76.9195462169742.*z)) + ...
        y.*(-30.0682112585625 - 1380.9597954037708.*z + y.*(2.626801985426835 + 703.695562834065.*z))))) + ...
        y.*(1187.3715515697959 + z.*(1458.233059470092 + ...
        z.*(-687.913805923122 + z.*(249.375342232496 + z.*(-63.313928772146 + 14.09317606630898.*z)))) + ...
        y.*(1760.062705994408 + y.*(-450.535298526802 + ...
        y.*(182.8520895502518 + y.*(-43.3206481750622 + 4.26033941694366.*y) + ...
        z.*(-595.457483974374 + (149.452282277512 - 72.9745838003176.*z).*z)) + ...
        z.*(1388.489628266536 + z.*(-409.779283929806 + (227.123395681188 - 22.2565468652826.*z).*z))) + ...
        z.*(-1721.528607567954 + z.*(674.819060538734 + ...
        z.*(-356.629112415276 + (88.4080716616 - 15.84003094423364.*z).*z)))));
    
    g08(x>0) = g08(x>0) + (11625.62913253464 + 1702.453469893412.*y(x>0)).*log(x(x>0));
    g08(x==0) = nan;
    
    gibbs = 0.5.*sfac.*g08;
    
elseif ns==0 & nt==1 & np==0
    
    g03 = 5.90578347909402 + z.*(-270.983805184062 + ...
        z.*(776.153611613101 + z.*(-196.51255088122 + (28.9796526294175 - 2.13290083518327.*z).*z))) + ...
        y.*(-24715.571866078 + z.*(2910.0729080936 + ...
        z.*(-1513.116771538718 + z.*(546.959324647056 + z.*(-111.1208127634436 + 8.68841343834394.*z)))) + ...
        y.*(2210.2236124548363 + z.*(-2017.52334943521 + ...
        z.*(1498.081172457456 + z.*(-718.6359919632359 + (146.4037555781616 - 4.9892131862671505.*z).*z))) + ...
        y.*(-592.743745734632 + z.*(1591.873781627888 + ...
        z.*(-1207.261522487504 + (608.785486935364 - 105.4993508931208.*z).*z)) + ...
        y.*(290.12956292128547 + z.*(-973.091553087975 + ...
        z.*(602.603274510125 + z.*(-276.361526170076 + 32.40953340386105.*z))) + ...
        y.*(-113.90630790850321 + y.*(21.35571525415769 - 67.41756835751434.*z) + ...
        z.*(381.06836198507096 + z.*(-133.7383902842754 + 49.023632509086724.*z)))))));
    
    g08 = x2.*(168.072408311545 + z.*(729.116529735046 + ...
        z.*(-343.956902961561 + z.*(124.687671116248 + z.*(-31.656964386073 + 7.04658803315449.*z)))) + ...
        x.*(-493.407510141682 + x.*(543.835333000098 + x.*(-196.028306689776 + 36.7571622995805.*x) + ...
        y.*(-137.1145018408982 + y.*(148.10030845687618 + y.*(-68.5590309679152 + 12.4848504784754.*y))) - ...
        22.6683558512829.*z) + z.*(-175.292041186547 + (83.1923927801819 - 29.483064349429.*z).*z) + ...
        y.*(-86.1329351956084 + z.*(766.116132004952 + z.*(-108.3834525034224 + 51.2796974779828.*z)) + ...
        y.*(-30.0682112585625 - 1380.9597954037708.*z + y.*(3.50240264723578 + 938.26075044542.*z)))) + ...
        y.*(1760.062705994408 + y.*(-675.802947790203 + ...
        y.*(365.7041791005036 + y.*(-108.30162043765552 + 12.78101825083098.*y) + ...
        z.*(-1190.914967948748 + (298.904564555024 - 145.9491676006352.*z).*z)) + ...
        z.*(2082.7344423998043 + z.*(-614.668925894709 + (340.685093521782 - 33.3848202979239.*z).*z))) + ...
        z.*(-1721.528607567954 + z.*(674.819060538734 + ...
        z.*(-356.629112415276 + (88.4080716616 - 15.84003094423364.*z).*z)))));
    
    g08(x>0) = g08(x>0) + 851.226734946706.*x2(x>0).*log(x(x>0));
    
    gibbs = (g03 + g08).*0.025;
    
elseif ns==0 & nt==0 & np==1
    
    g03 = 100015.695367145 + z.*(-5089.1530840726 + ...
        z.*(853.5533353388611 + z.*(-133.2587017014444 + (21.0131554401542 - 3.278571068826234.*z).*z))) + ...
        y.*(-270.983805184062 + z.*(1552.307223226202 + ...
        z.*(-589.53765264366 + (115.91861051767 - 10.664504175916349.*z).*z)) + ...
        y.*(1455.0364540468 + z.*(-1513.116771538718 + ...
        z.*(820.438986970584 + z.*(-222.2416255268872 + 21.72103359585985.*z))) + ...
        y.*(-672.50778314507 + z.*(998.720781638304 + ...
        z.*(-718.6359919632359 + (195.2050074375488 - 8.31535531044525.*z).*z)) + ...
        y.*(397.968445406972 + z.*(-603.630761243752 + (456.589115201523 - 105.4993508931208.*z).*z) + ...
        y.*(-194.618310617595 + y.*(63.5113936641785 - 9.63108119393062.*y + ...
        z.*(-44.5794634280918 + 24.511816254543362.*z)) + ...
        z.*(241.04130980405 + z.*(-165.8169157020456 + 25.92762672308884.*z)))))));
    
    
    g08 = x2.*(-3310.49154044839 + z.*(769.588305957198 + ...
        z.*(-289.5972960322374 + (63.3632691067296 - 13.1240078295496.*z).*z)) + ...
        x.*(199.459603073901 + x.*(-54.7919133532887 + 36.0284195611086.*x - 22.6683558512829.*y + ...
        (-8.16387957824522 - 90.52653359134831.*z).*z) + ...
        z.*(-104.588181856267 + (204.1334828179377 - 13.65007729765128.*z).*z) + ...
        y.*(-175.292041186547 + (166.3847855603638 - 88.449193048287.*z).*z + ...
        y.*(383.058066002476 + y.*(-460.319931801257 + 234.565187611355.*y) + ...
        z.*(-108.3834525034224 + 76.9195462169742.*z)))) + ...
        y.*(729.116529735046 + z.*(-687.913805923122 + ...
        z.*(374.063013348744 + z.*(-126.627857544292 + 35.23294016577245.*z))) + ...
        y.*(-860.764303783977 + y.*(694.244814133268 + ...
        y.*(-297.728741987187 + (149.452282277512 - 109.46187570047641.*z).*z) + ...
        z.*(-409.779283929806 + (340.685093521782 - 44.5130937305652.*z).*z)) + ...
        z.*(674.819060538734 + z.*(-534.943668622914 + (176.8161433232 - 39.600077360584095.*z).*z)))));
    
    gibbs = (g03 + g08).*1e-8;
    % Note. This pressure derivative of the gibbs function is in units of (J/kg) (Pa^-1) = m^3/kg
    
elseif ns==1 & nt==1 & np==0
    
    g08 = 1187.3715515697959 + z.*(1458.233059470092 + ...
        z.*(-687.913805923122 + z.*(249.375342232496 + z.*(-63.313928772146 + 14.09317606630898.*z)))) + ...
        x.*(-1480.222530425046 + x.*(2175.341332000392 + x.*(-980.14153344888 + 220.542973797483.*x) + ...
        y.*(-548.4580073635929 + y.*(592.4012338275047 + y.*(-274.2361238716608 + 49.9394019139016.*y))) - ...
        90.6734234051316.*z) + z.*(-525.876123559641 + (249.57717834054571 - 88.449193048287.*z).*z) + ...
        y.*(-258.3988055868252 + z.*(2298.348396014856 + z.*(-325.1503575102672 + 153.8390924339484.*z)) + ...
        y.*(-90.2046337756875 - 4142.8793862113125.*z + y.*(10.50720794170734 + 2814.78225133626.*z)))) + ...
        y.*(3520.125411988816 + y.*(-1351.605895580406 + ...
        y.*(731.4083582010072 + y.*(-216.60324087531103 + 25.56203650166196.*y) + ...
        z.*(-2381.829935897496 + (597.809129110048 - 291.8983352012704.*z).*z)) + ...
        z.*(4165.4688847996085 + z.*(-1229.337851789418 + (681.370187043564 - 66.7696405958478.*z).*z))) + ...
        z.*(-3443.057215135908 + z.*(1349.638121077468 + ...
        z.*(-713.258224830552 + (176.8161433232 - 31.68006188846728.*z).*z))));
    
    g08(x>0) = g08(x>0) + 1702.453469893412.*log(x(x>0));
    g08(SA==0) = nan;
    
    gibbs = 0.5.*sfac.*0.025.*g08;
    
elseif ns==1 & nt==0 & np==1
    
    g08 =  -6620.98308089678 + z.*(1539.176611914396 + ...
        z.*(-579.1945920644748 + (126.7265382134592 - 26.2480156590992.*z).*z)) + ...
        x.*(598.378809221703 + x.*(-219.1676534131548 + 180.142097805543.*x - 90.6734234051316.*y + ...
        (-32.65551831298088 - 362.10613436539325.*z).*z) + ...
        z.*(-313.764545568801 + (612.4004484538132 - 40.95023189295384.*z).*z) + ...
        y.*(-525.876123559641 + (499.15435668109143 - 265.347579144861.*z).*z + ...
        y.*(1149.174198007428 + y.*(-1380.9597954037708 + 703.695562834065.*y) + ...
        z.*(-325.1503575102672 + 230.7586386509226.*z)))) + ...
        y.*(1458.233059470092 + z.*(-1375.827611846244 + ...
        z.*(748.126026697488 + z.*(-253.255715088584 + 70.4658803315449.*z))) + ...
        y.*(-1721.528607567954 + y.*(1388.489628266536 + ...
        y.*(-595.457483974374 + (298.904564555024 - 218.92375140095282.*z).*z) + ...
        z.*(-819.558567859612 + (681.370187043564 - 89.0261874611304.*z).*z)) + ...
        z.*(1349.638121077468 + z.*(-1069.887337245828 + (353.6322866464 - 79.20015472116819.*z).*z))));
    
    gibbs = g08.*sfac.*0.5e-8;
    % Note. This derivative of the Gibbs function is in units of (m^3/kg)/(g/kg) = m^3/g,
    % that is, it is the derivative of specific volume with respect to Absolute
    % Salinity measured in g/kg.
    
elseif ns==0 & nt==1 & np==1
    
    g03 = -270.983805184062 + z.*(1552.307223226202 + z.*(-589.53765264366 + ...
        (115.91861051767 - 10.664504175916349.*z).*z)) + ...
        y.*(2910.0729080936 + z.*(-3026.233543077436 + ...
        z.*(1640.877973941168 + z.*(-444.4832510537744 + 43.4420671917197.*z))) + ...
        y.*(-2017.52334943521 + z.*(2996.162344914912 + ...
        z.*(-2155.907975889708 + (585.6150223126464 - 24.946065931335752.*z).*z)) + ...
        y.*(1591.873781627888 + z.*(-2414.523044975008 + (1826.356460806092 - 421.9974035724832.*z).*z) + ...
        y.*(-973.091553087975 + z.*(1205.20654902025 + z.*(-829.084578510228 + 129.6381336154442.*z)) + ...
        y.*(381.06836198507096 - 67.41756835751434.*y + z.*(-267.4767805685508 + 147.07089752726017.*z))))));
    
    g08 = x2.*(729.116529735046 + z.*(-687.913805923122 + ...
        z.*(374.063013348744 + z.*(-126.627857544292 + 35.23294016577245.*z))) + ...
        x.*(-175.292041186547 - 22.6683558512829.*x + (166.3847855603638 - 88.449193048287.*z).*z + ...
        y.*(766.116132004952 + y.*(-1380.9597954037708 + 938.26075044542.*y) + ...
        z.*(-216.7669050068448 + 153.8390924339484.*z))) + ...
        y.*(-1721.528607567954 + y.*(2082.7344423998043 + ...
        y.*(-1190.914967948748 + (597.809129110048 - 437.84750280190565.*z).*z) + ...
        z.*(-1229.337851789418 + (1022.055280565346 - 133.5392811916956.*z).*z)) + ...
        z.*(1349.638121077468 + z.*(-1069.887337245828 + (353.6322866464 - 79.20015472116819.*z).*z))));
    
    gibbs = (g03 + g08).*2.5e-10;
    % Note. This derivative of the Gibbs function is in units of (m^3/(K kg)),
    % that is, the pressure of the derivative in Pa.
    
elseif ns==2 & nt==0 & np==0
    
    g08 = 5812.814566267320 + 851.2267349467060.*y + x.*(-3648.219935726910 ...
       + x.*(8103.204624147880 + x.*(-8187.513078222526 + x.*(4495.214854534080 ...
       - 850.3093707944657.*x))) + y.*(-740.1112652125230 + x.*(2175.341332000392 ...
       + x.*(-1470.212300173320 + 441.0859475949660.*x)) + y.*(-64.59970139670629 ...
       - 274.2290036817964.*x + y.*(-15.03410562928125 + 197.4670779425016.*x ...
       + y.*(1.313400992713418 - 68.55903096791521.*x + 9.987880382780312.*x.*y)))) ...
       + z.*(299.1894046108515 + x.*(-219.1676534131548 + 270.2131467083145.*x) ...
       + y.*(-262.9380617798205 - 90.67342340513160.*x + y.*(+ 574.5870990037140 ...
       + y.*(-690.4798977018854 + 351.8477814170325.*y))) + z.*(-78.44113639220025 ...
       - 16.32775915649044.*x + y.*(124.7885891702729 - 81.28758937756680.*y) ...
       + z.*(102.0667414089689 - 120.7020447884644.*x + y.*(-44.22459652414350 ...
       + 38.45977310848710.*y) - 5.118778986619230.*z))));
    
    g08(x>0) = g08(x>0)./(x2(x>0));
    g08(x==0) = NaN;
    
    gibbs = 0.5.*sfac.*sfac.*g08;
    
elseif ns==0 & nt==2 & np==0
    
    g03 = -24715.571866078 + z.*(2910.0729080936 + z.* ...
        (-1513.116771538718 + z.*(546.959324647056 + z.*(-111.1208127634436 + 8.68841343834394.*z)))) + ...
        y.*(4420.4472249096725 + z.*(-4035.04669887042 + ...
        z.*(2996.162344914912 + z.*(-1437.2719839264719 + (292.8075111563232 - 9.978426372534301.*z).*z))) + ...
        y.*(-1778.231237203896 + z.*(4775.621344883664 + ...
        z.*(-3621.784567462512 + (1826.356460806092 - 316.49805267936244.*z).*z)) + ...
        y.*(1160.5182516851419 + z.*(-3892.3662123519 + ...
        z.*(2410.4130980405 + z.*(-1105.446104680304 + 129.6381336154442.*z))) + ...
        y.*(-569.531539542516 + y.*(128.13429152494615 - 404.50541014508605.*z) + ...
        z.*(1905.341809925355 + z.*(-668.691951421377 + 245.11816254543362.*z))))));
    
    g08 = x2.*(1760.062705994408 + x.*(-86.1329351956084 + ...
        x.*(-137.1145018408982 + y.*(296.20061691375236 + y.*(-205.67709290374563 + 49.9394019139016.*y))) + ...
        z.*(766.116132004952 + z.*(-108.3834525034224 + 51.2796974779828.*z)) + ...
        y.*(-60.136422517125 - 2761.9195908075417.*z + y.*(10.50720794170734 + 2814.78225133626.*z))) + ...
        y.*(-1351.605895580406 + y.*(1097.1125373015109 + y.*(-433.20648175062206 + 63.905091254154904.*y) + ...
        z.*(-3572.7449038462437 + (896.713693665072 - 437.84750280190565.*z).*z)) + ...
        z.*(4165.4688847996085 + z.*(-1229.337851789418 + (681.370187043564 - 66.7696405958478.*z).*z))) + ...
        z.*(-1721.528607567954 + z.*(674.819060538734 + ...
        z.*(-356.629112415276 + (88.4080716616 - 15.84003094423364.*z).*z))));
    
    gibbs = (g03 + g08).*0.000625;
    
elseif ns==0 & nt==0 & np==2
    
    g03 = -5089.1530840726 + z.*(1707.1066706777221 + ...
        z.*(-399.7761051043332 + (84.0526217606168 - 16.39285534413117.*z).*z)) + ...
        y.*(1552.307223226202 + z.*(-1179.07530528732 + (347.75583155301 - 42.658016703665396.*z).*z) + ...
        y.*(-1513.116771538718 + z.*(1640.877973941168 + z.*(-666.7248765806615 + 86.8841343834394.*z)) + ...
        y.*(998.720781638304 + z.*(-1437.2719839264719 + (585.6150223126464 - 33.261421241781.*z).*z) + ...
        y.*(-603.630761243752 + (913.178230403046 - 316.49805267936244.*z).*z + ...
        y.*(241.04130980405 + y.*(-44.5794634280918 + 49.023632509086724.*z) + ...
        z.*(-331.6338314040912 + 77.78288016926652.*z))))));
    
    g08 = x2.*(769.588305957198 + z.*(-579.1945920644748 + (190.08980732018878 - 52.4960313181984.*z).*z) + ...
        x.*(-104.588181856267 + x.*(-8.16387957824522 - 181.05306718269662.*z) + ...
        (408.2669656358754 - 40.95023189295384.*z).*z + ...
        y.*(166.3847855603638 - 176.898386096574.*z + y.*(-108.3834525034224 + 153.8390924339484.*z))) + ...
        y.*(-687.913805923122 + z.*(748.126026697488 + z.*(-379.883572632876 + 140.9317606630898.*z)) + ...
        y.*(674.819060538734 + z.*(-1069.887337245828 + (530.4484299696 - 158.40030944233638.*z).*z) + ...
        y.*(-409.779283929806 + y.*(149.452282277512 - 218.92375140095282.*z) + ...
        (681.370187043564 - 133.5392811916956.*z).*z))));
    
    gibbs = (g03 + g08).*1e-16;
    % Note. This is the second derivative of the Gibbs function with respect to
    % pressure, measured in Pa.  This derivative has units of (J/kg) (Pa^-2).
    
elseif ns==0 & nt==3 & np==0
    
    g03 = 4420.44722490967251 + y.*(-3556.46247440779189 + y.*(3481.55475505542563 ...
        + y.*(-2278.12615817006417 + 640.671457624730749.*y))) + z.*(-4035.04669887042019 ...
        + y.*(9551.24268976732856 + y.*(-11677.0986370556993 + y.*(7621.36723970141975 ...
        - 2022.52705072543023.*y))) + z.*(2996.16234491491196 + y.*(-7243.56913492502372 ...
        + y.*(7231.23929412149937 - 2674.76780568550794.*y)) ...
        + z.*(-1437.27198392647188 + y.*(3652.71292161218389 + y.*(-3316.33831404091188 ...
        + 980.472650181734480.*y)) + z.*(292.807511156323187 + y.*(-632.996105358724890 ...
        + 388.914400846332569.*y) - 9.97842637253430098.*z))));
    
    g08 = x2.*(-1351.605895580406 + x.*(-60.1364225171250 + 296.200616913752.*x) ...
        + y.*(2194.22507460302 + y.*(-1299.61944525187 + 255.620365016620.*y) ...
        + x.*(21.0144158834147 + x.*(-411.354185807491 + 149.818205741705.*y))) ...
        + z.*(4165.46888479961 - 2761.91959080754.*x + y.*(-7145.48980769249 ...
        + 5629.56450267252.*x) + z.*(-1229.33785178942 + 1793.42738733014.*y ...
        + z.*(681.370187043564 - 875.695005603811.*y - 66.7696405958478.*z))));
    
    gibbs = (g03 + g08).*0.000625.*0.025;
    
elseif ns==0 & nt==2 & np==1
    
    g03 = 2910.07290809360 + y.*(-4035.04669887042 + y.*(4775.62134488366  ...
        + y.*(-3892.36621235190 + y.*(1905.34180992535 - 404.505410145086.*y))))  ...
        + z.*(-3026.23354307744 + y.*(5992.32468982982 + y.*(-7243.56913492502  ...
        + y.*(4820.82619608100 - 1337.38390284275.*y))) + z.*(1640.87797394117  ...
        + y.*(-4311.81595177942 + y.*(5479.06938241828 + y.*(-3316.33831404091  ...
        + 735.354487636301.*y)))+ z.*(-444.483251053774 + y.*(1171.23004462529  ...
        + y.*(-1265.99221071745 + 518.552534461777.*y))  ...
        + z.*(43.4420671917197 - 49.8921318626715.*y))));
    
    g08 = x2.*(-1721.52860756795 + 766.116132004952.*x - 2761.91959080754.*x.*y  ...
        + y.*(4165.46888479961 + y.*(-3572.74490384624 + 2814.78225133626.*x))  ...
        + z.*(1349.63812107747 - 216.766905006845.*x + y.*(-2458.67570357884 ...
        + 1793.42738733014.*y) + z.*(-1069.88733724583 + 153.839092433949.*x  ...
        + y.*(+ 2044.11056113069 - 1313.54250840572.*y) + z.*(353.632286646400  ...
        - 267.078562383391.*y - 79.2001547211682.*z))));
    
    gibbs = (g03 + g08).*2.5e-10.*0.025;
    
elseif ns==1 & nt==1 & np==1
    
    g08 =  1458.23305947009 + x.*(-525.876123559641 - 90.6734234051316.*x ...
        + y.*(2298.34839601486 + y.*(-4142.87938621131 + 2814.78225133626.*y))) ...
        + y.*(-3443.05721513591 + y.*(4165.46888479961 - 2381.82993589750.*y)) ...
        + z.*(-1375.82761184624 + x.*(499.154356681091 - 650.300715020534.*y) ...
        + y.*(2699.27624215494 + y.*(-2458.67570357884 + 1195.61825822010.*y)) ...
        + z.*(748.126026697488 + x.*(-265.347579144861 + 461.517277301845.*y) ...
        + y.*(-2139.77467449166 + y.*(2044.11056113069 - 875.695005603811.*y)) ...
        + z.*(-253.255715088584 + y.*(707.264573292800 - 267.078562383391.*y) ...
        + z.*(70.4658803315449 - 158.400309442336.*y))));
    
    gibbs = g08.*sfac.*0.5e-8.*0.025;
    
elseif ns==2 & nt==0 & np==1
    
    g08 = 299.1894046108515 + x.*(-219.1676534131548 + 270.2131467083145.*x) ...
        + y.*(-262.9380617798205 - 90.67342340513160.*x + y.*(+ 574.5870990037140 ...
        + y.*(-690.4798977018854 + 351.8477814170325.*y))) + z.*(- 156.8822727844005 ...
        - 32.65551831298088.*x + y.*( 249.5771783405457 - 162.5751787551336.*y) ...
        + z.*(306.2002242269066 - 362.1061343653932.*x + y.*(-132.6737895724305 ...
        + 115.3793193254613.*y) - 20.47511594647692.*z));

    g08(x>0) = g08(x>0)./x(x>0);
    g08(x==0) = nan;
    
    gibbs = 0.5.*sfac.*sfac.*g08.*1e-8;
    
elseif ns==1 & nt==2 & np==0
    
    g08 = 3520.125411988816 + x.*(-258.3988055868252 - 548.4580073635929.*x) ...
        + y.*(-2703.211791160812 + x.*(-180.409267551375 + 1184.8024676550094.*x) ...
        + y.*(2194.2250746030217 + x.*(31.52162382512202 + x.*(-822.70837161498247 ...
        + 199.7576076556064.*y)) + y.*(-866.4129635012441 + 127.8101825083098.*y))) ...
        + z.*(-3443.057215135908 + 2298.348396014856.*x + y.*(8330.937769599217 ...
        - 8285.758772422625.*x + y.*(-7145.4898076924878 + 8444.34675400878.*x)) ...
        + z.*(1349.638121077468 - 325.1503575102672.*x + y.*(-2458.675703578836 ...
        + 1793.427387330144.*y) + z.*(-713.258224830552 + 153.8390924339484.*x ...
        + y.*(1362.740374087128 - 875.6950056038112.*y) + z.*(176.8161433232 ...
        - 133.5392811916956.*y - 31.6800618884673.*z))));
    
    g08(x==0) = nan;
    
    gibbs = 0.5.*sfac.*0.025.*g08.*0.025;
    
elseif ns==2 & nt==1 & np==0
        
    g08 = 851.22673494670596 + x.*(-740.11126521252299 + x.*(2175.3413320003920 ...
        + x.*(-1470.2123001733199 + 441.08594759496600.*x)) + y.*(-129.19940279341259 ...
        - 548.45800736359286.*x + y.*(-45.102316887843749 + 592.40123382750477.*x ...
        + y.*(5.2536039708536704 + x.*(-274.23612387166082 + 49.939401913901600.*y)))) ...
        + z.*(-262.93806177982049 - 90.673423405131601.*x + y.*(1149.1741980074280 ...
        + y.*(- 2071.4396931056563 + 1407.3911256681299.*y)) + z.*( 124.78858917027286 ...
        - 162.57517875513361.*y + z.*(-44.224596524143500 + 76.919546216974197.*y))));

    g08(x>0) = g08(x>0)./x2(x>0);
    g08(x==0) = nan;
    
    gibbs = 0.5.*sfac.*sfac.*g08.*0.025;
    
elseif ns==1 & nt==0 & np==2
    
    g08 = 1539.17661191439606 + x.*(-313.764545568801 - 32.6555183129809.*x) ...
        + y.*(-1375.827611846244 + 499.15435668109143.*x + y.*(1349.638121077468 ...
        - 325.150357510267.*x + y.*(-819.558567859612 + 298.904564555024.*y))) ...
        + z.*(-1158.3891841289496 + x.*(1224.8008969076263 - 724.2122687307865.*x) ...
        + y.*(1496.252053394976 - 530.695158289722.*x + y.*(-2139.774674491656 ...
        + 461.5172773018452.*x + y.*(1362.740374087128 - 437.847502801906.*y))) ...
        + z.*(380.179614640378 - 122.850695678862.*x + y.*(-759.767145265752 ...
        + y.*(1060.8968599392 - 267.078562383391.*y)) + z.*(-104.992062636397 ...
        + y.*(281.86352132618 - 316.800618884673.*y))));
    
    gibbs = g08.*sfac.*0.5e-8.*1e-8;
    
elseif ns==0 & nt==1 & np==2
    
    g03 = 1552.307223226202 + y.*(-3026.233543077436 + y.*(2996.162344914912 ...
        + y.*(-2414.523044975008 + y.*(1205.20654902025 - 267.4767805685508.*y)))) ...
        + z.*(-1179.07530528732 + y.*(3281.755947882336 + y.*(-4311.815951779416 ...
        + y.*(3652.712921612184 + y.*(-1658.169157020456 + 294.141795054520.*y)))) ...
        + z.*(347.755831553010 + y.*(-1333.449753161323 + y.*(1756.84506693794 ...
        + y.*(-1265.99221071745 + 388.9144008463326.*y))) + z.*(-42.6580167036654 ...
        + y.*(173.7682687668788 - 99.7842637253430.*y))));
    
    g08 = x2.*(-687.9138059231220 + 166.3847855603638.*x + y.*(1349.638121077468 ...
        - 216.7669050068448.*x + y.*(-1229.337851789418 + 597.8091291100480.*y)) ...
        + z.*(748.1260266974880 - 176.8983860965740.*x + y.*(-2139.774674491656 ...
        + 307.6781848678968.*x + y.*(2044.110561130692 - 875.6950056038113.*y)) ...
        + z.*(-379.883572632876 + y.*(1060.89685993920 - 400.6178435750868.*y) ...
        + z.*(140.9317606630898 - 316.8006188846728.*y))));
    
    gibbs = (g03 + g08).*2.5e-10.*1e-8;
    % Note. This derivative of the Gibbs function is in units of (m^3/(K kg)),
    % that is, the pressure of the derivative in Pa.

end

end
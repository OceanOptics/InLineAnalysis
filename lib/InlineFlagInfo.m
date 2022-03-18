function InLineFlag = InlineFlagInfo(instrument)
% Read flags from binary flags 
%
% Author: Guillaume Bourdin
% Date: 2021-05-28
%%
% get rid of numbers in instrument name
instrument = instrument(isstrprop(instrument,'alpha'));
warning('off')

switch instrument
  case {'ACS', 'AC'}
    InLineFlag = table();
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'ap430_700_neg';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'ap 430-700 < -0.0015';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'cp_neg';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'cp < -0.0015';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'ap_shape';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'ap640 > ap676';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'cp_over10';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'cp > 10';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'noisy600_650';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'sum(abs(d(ap)/d(lambda(600-650)))) / ap450nm';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'ap460_640_04_450';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'd(ap)/d(lambda460-640) > 0.4 * ap_{450nm} | abs(d"(ap)/d(lambda460-640)) > 0.05)';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'positive_ap450_570';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = '4 consecutive d(ap)/d(lambda485-570) > 0'; 
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'poc_flag';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'POC cp < 0';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'chl_ap676lh_flag';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'ap676 line height chlorophyll < 0 | complex ap676';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'gamma_flag';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'gamma cp < 0';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'chl_Halh_flag';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'Housekeeper line height chlorophyll < 0';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'HH_mphi_flag';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'Housekeeper and Haentjens G50 > 0';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'HH_G50_flag';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'Housekeeper and Haentjens G50 < 0';
    
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'chlratio_flag';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'high variability between the two chlorophyll algorithms';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'gamma_suspicious';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'gamma < 0.2';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'poc_suspicious';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'POC suspicious';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'chl_ap676lh_suspicious';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'chlorophyll from ap676 line height suspicious';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'chl_Halh_suspicious';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'chlorophyll from Housekeeper line height suspicious';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'HH_G50_mphi_suspicious';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'Housekeeper and Haentjens G50 and mphi suspicious';
    
%   case {'BB', 'HBB'}
%   case 'TSG'
%   case 'PAR'
%   case 'WSCD'
%   case 'LISST'
  otherwise
    warning('%s not supported for flag reading', instrument)
end
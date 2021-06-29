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
    InLineFlag.name{size(InLineFlag,1)+1} = 'gamma_flag';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'gamma cp < 0';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'poc_flag';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'POC cp < 0';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'chl_ap676lh_flag';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'chlorophyll from ap676 line height < 0';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'chl_Halh_flag';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'chlorophyll from Housekeeper line height < 0';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'HH_G50_flag';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'Housekeeper and Haentjens G50 < 0';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'chlratio_flag';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'high variability between the two chlorophyll algorithms';
    
    InLineFlag.name{size(InLineFlag,1)+1} = 'gamma_suspicious';
    InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
    InLineFlag.description{size(InLineFlag,1)} = 'gamma suspicious';
    
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
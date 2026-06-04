function InLineFlag = InlineFlagInfo(instrument, product_type)
% Read flags from binary flags 
%
% Author: Guillaume Bourdin
% Date: 2021-05-28
%%
% get rid of numbers in instrument name
% instrument = instrument(isstrprop(instrument,'alpha'));
warning('off')

if nargin == 1
  product_type = 'particulate';
elseif nargin == 2
  if ~any(strcmp(product_type, {'particulate', 'part', 'dissolved', 'diw'}))
    error("Input 'product_type' not supported, should be either 'particulate' or 'dissolved'")
  end
end

% Built flag table
InLineFlag = table();
switch instrument
  case {'ACS', 'AC'}
    if any(strcmp(product_type, {'particulate', 'part'}))
      
      InLineFlag.name{size(InLineFlag,1)+1} = 'deltaTS_not_corrected_total';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'Total measurements not T/S corrected before substraction of filtered from total';
      
      InLineFlag.name{size(InLineFlag,1)+1} = 'deltaTS_not_corrected_filt0';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'Filtered preceeding current total event not T/S corrected before substraction of filtered from total';
      
      InLineFlag.name{size(InLineFlag,1)+1} = 'deltaTS_not_corrected_filt1';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'Filtered following current total event not T/S corrected before substraction of filtered from total';

      % InLineFlag.name{size(InLineFlag,1)+1} = 'fCDOM_mix_a_cluster0';
      % InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      % InLineFlag.description{size(InLineFlag,1)} = 'fCDOM interpolation: mixed clusters within start a filter event';
      % 
      % InLineFlag.name{size(InLineFlag,1)+1} = 'fCDOM_mix_a_cluster1';
      % InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      % InLineFlag.description{size(InLineFlag,1)} = 'fCDOM interpolation: mixed clusters within end a filter event';
      % 
      % InLineFlag.name{size(InLineFlag,1)+1} = 'fCDOM_mix_c_cluster0';
      % InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      % InLineFlag.description{size(InLineFlag,1)} = 'fCDOM interpolation: mixed clusters within start c filter event';
      % 
      % InLineFlag.name{size(InLineFlag,1)+1} = 'fCDOM_mix_c_cluster1';
      % InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      % InLineFlag.description{size(InLineFlag,1)} = 'fCDOM interpolation: mixed clusters within end c filter event';
      % 
      % InLineFlag.name{size(InLineFlag,1)+1} = 'fCDOM_a_cluster_chg';
      % InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      % InLineFlag.description{size(InLineFlag,1)} = 'fCDOM interpolation: change of a clusters between start and end filter events';
      % 
      % InLineFlag.name{size(InLineFlag,1)+1} = 'fCDOM_c_cluster_chg';
      % InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      % InLineFlag.description{size(InLineFlag,1)} = 'fCDOM interpolation: change of c clusters between start and end filter events';
      
      InLineFlag.name{size(InLineFlag,1)+1} = 'flag_linear_interp_a';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'a filter events linearly interpolated';

      InLineFlag.name{size(InLineFlag,1)+1} = 'flag_linear_interp_c';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'c filter events linearly interpolated';

      InLineFlag.name{size(InLineFlag,1)+1} = 'flag_Tresidual';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'Residual temperature correction failed';
  
      % InLineFlag.name{size(InLineFlag,1)+1} = 'a_filt0_not_clustered';
      % InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      % InLineFlag.description{size(InLineFlag,1)} = 'fCDOM interpolation: a_filt of start filter event not included in a_filt/fdom clustering';
      % 
      % InLineFlag.name{size(InLineFlag,1)+1} = 'a_filt1_not_clustered';
      % InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      % InLineFlag.description{size(InLineFlag,1)} = 'fCDOM interpolation: a_filt of end filter event not included in a_filt/fdom clustering';
      % 
      % InLineFlag.name{size(InLineFlag,1)+1} = 'c_filt0_not_clustered';
      % InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      % InLineFlag.description{size(InLineFlag,1)} = 'fCDOM interpolation: c_filt of start filter event not included in c_filt/fdom clustering';
      % 
      % InLineFlag.name{size(InLineFlag,1)+1} = 'c_filt1_not_clustered';
      % InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      % InLineFlag.description{size(InLineFlag,1)} = 'fCDOM interpolation: c_filt of end filter event not included in c_filt/fdom clustering';
      % 
      % InLineFlag.name{size(InLineFlag,1)+1} = 'flag_a_negative_slope';
      % InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      % InLineFlag.description{size(InLineFlag,1)} = 'fCDOM interpolation: fCDOM a filter interpolation with negative slope';
      % 
      % InLineFlag.name{size(InLineFlag,1)+1} = 'flag_c_negative_slope';
      % InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      % InLineFlag.description{size(InLineFlag,1)} = 'fCDOM interpolation: fCDOM c filter interpolation with negative slope';

      InLineFlag.name{size(InLineFlag,1)+1} = 'cp_neg';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'cp < -0.0015';
      
      InLineFlag.name{size(InLineFlag,1)+1} = 'ap_neg';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'ap 430-700 < -0.01';
      
      InLineFlag.name{size(InLineFlag,1)+1} = 'ap_shape';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'ap640 > ap676';
      
      InLineFlag.name{size(InLineFlag,1)+1} = 'ap_step';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'abs(d(ap)/d(lambda460-640)) > 3 * prctile(abs(d(ap)/d(lambda460-640)),95,2)';
  
      InLineFlag.name{size(InLineFlag,1)+1} = 'cp_step';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'abs(d(cp)/d(lambda460-640)) > 3 * prctile(abs(d(cp)/d(lambda460-640)),95,2)';
      
      InLineFlag.name{size(InLineFlag,1)+1} = 'ap430_700_neg';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'ap 430-700 < -0.0015';
      
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
      InLineFlag.description{size(InLineFlag,1)} = 'Housekeeper and Haentjens G50 < 0 or > 500';
      
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

    else

      InLineFlag.name{size(InLineFlag,1)+1} = 'deltaTSnotcorrected_filt';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'Filtered event not T/S corrected before before substraction of DIW from filtered';
      
    end

  case {'BB', 'HBB'}

    if any(strcmp(product_type, {'particulate', 'part'}))

      InLineFlag.name{size(InLineFlag,1)+1} = 'no_attenuation_correction';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'no attenuation correction applied';
      
      InLineFlag.name{size(InLineFlag,1)+1} = 'no_ag_attenuation_correction';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'no ag for attenuation correction';
      
      InLineFlag.name{size(InLineFlag,1)+1} = 'simplified_Doxaran_correction';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'Simplified Doxaran attenuation correction';

      InLineFlag.name{size(InLineFlag,1)+1} = 'ap_interpolated';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'ap missing and interpolated for attenuation correction';

      InLineFlag.name{size(InLineFlag,1)+1} = 'cp_interpolated';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'cp missing and interpolated for attenuation correction';

      InLineFlag.name{size(InLineFlag,1)+1} = 'ag_interpolated';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'ag missing and interpolated for attenuation correction';

      InLineFlag.name{size(InLineFlag,1)+1} = 'cg_interpolated';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'ap missing and interpolated for attenuation correction';
      
      InLineFlag.name{size(InLineFlag,1)+1} = 'ag_input_interpolated_linearly';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'input ap interpolated linearly for attenuation correction';
      
      InLineFlag.name{size(InLineFlag,1)+1} = 'cg_input_interpolated_linearly';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'input cp interpolated linearly for attenuation correction';
      
      InLineFlag.name{size(InLineFlag,1)+1} = 'fdom_input_interpolated_linearly';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'input fdom interpolated linearly for attenuation correction';
      
      InLineFlag.name{size(InLineFlag,1)+1} = 'deltaTSnotcorrected_filt0';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'Filtered preceeding current total event not T/S corrected before substraction of filtered from total';
      
      InLineFlag.name{size(InLineFlag,1)+1} = 'deltaTSnotcorrected_filt1';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'Filtered following current total event not T/S corrected before substraction of filtered from total';

      InLineFlag.name{size(InLineFlag,1)+1} = 'deltaTSnotcorrected_total';
      InLineFlag.bit(size(InLineFlag,1)) = size(InLineFlag,1)-1;
      InLineFlag.description{size(InLineFlag,1)} = 'Total measurements not T/S corrected before substraction of filtered from total';

    else
      
      warning('TODO: FINISH IMPLEMENTING betag flags')

    end

%   case 'TSG'
%   case 'PAR'
%   case 'WSCD'
%   case 'LISST'
  otherwise
    warning('%s not supported for flag reading', instrument)
end
function [betap_corr, bbp_corr, flags] = betap_Kcorrection(beta_total, beta_filt_interp, ...
  lambda_beta, k_exp, Chi, ac_p, ac_g, lambda_a, lambda_c, instrument_types, fdom_ag_correlation)
  % This function performs the attenuation correction on beta measured with backscattering sensors
  % This method is based on the HS4 correction devised by Doxaran et al., 2016 (10.1364/OE.24.003615)
  % but rather than iterating, we expand the exponential in that paper in a two terms Taylor series.
  % Author: Guillaume Bourdin & Emmanuel Boss
  % Date: 2024-08-29
  %
  % Inputs:
  % - beta_total: uncalibrated measured beta total
  % - beta_filt_interp: uncalibrated measured beta filtered
  % - lambda_beta: beta wavelength
  % - k_exp: optical pathlength
  % - Chi: factor needed to convert the particulate VSF to particulate backscattering from Sullivan et al. 2013
  % - ac_p: measured calibrated absorption and attenuation of particulate matter
  % - ac_g: measured calibrated absorption and attenuation of dissolved matter
  % - lambda_a: absorption channel wavelength
  % - lambda_c: absorption channel wavelength
  % - fdom_ag_correlation (option) if ag and cg are not computed from underway data
  %
  %%
  if nargin < 10
    instrument_types.c_type = 'ACS';
    instrument_types.bb_type = 'HyperBB';
  end
  if nargin < 11
    fdom_ag_correlation = [];
  end
  
  % round all timestamps to flag interpolated data
  beta_total = round_timestamp(beta_total, minutes(1));
  beta_filt_interp = round_timestamp(beta_filt_interp, minutes(1));
  
  %% create flag table
  betap_approx = beta_total.beta - beta_filt_interp.beta;
  bbp_approx = 2 * pi * Chi .* betap_approx;
  vnames = {'no_attenuation_correction', 'no_ag_attenuation_correction', ...
    'simplified_Doxaran_correction', 'ap_interpolated', 'cp_interpolated', ...
    'ag_interpolated', 'cg_interpolated', 'ag_input_interpolated_linearly', ...
    'cg_input_interpolated_linearly','fdom_input_interpolated_linearly'};
  flags = table('Size', [size(beta_filt_interp.dt, 1) size(vnames,2)], ...
    'VariableTypes', repmat({'logical'},1,size(vnames,2)), 'VariableNames', vnames);
  
  % Format ap, cp and compute bp if particulate data available
  agcg_modelled = false;
  if ~isempty(ac_p)
    ac_p = round_timestamp(ac_p, minutes(1));
    % if ag/cg modelled from FDOM exist, use them instead of DIW ag/cg
    if any(strcmp(ac_p.Properties.VariableNames, 'ag_modelled')) && any(strcmp(ac_p.Properties.VariableNames, 'cg_modelled'))
      ac_g = ac_p(:, strcmp(ac_p.Properties.VariableNames, 'dt') | strcmp(ac_p.Properties.VariableNames, 'ag_modelled')  | strcmp(ac_p.Properties.VariableNames, 'cg_modelled'));
      ac_g = renamevars(ac_g,{'ag_modelled','cg_modelled'},{'ag','cg'});
      agcg_modelled = true;
    else
      ac_g = round_timestamp(ac_g, minutes(1));
    end
    [ap_interp, flags.ap_interpolated] = interp_on_bb_lambda(ac_p.dt, ac_p.ap, lambda_a, beta_total.dt, lambda_beta);
    [cp_interp, flags.cp_interpolated] = interp_on_bb_lambda(ac_p.dt, ac_p.cp, lambda_c, beta_total.dt, lambda_beta);
    chl_interp = interp_on_bb_lambda(ac_p.dt, ac_p.chl_ap676lh, [], beta_total.dt, []);
    % compute bp
    bp_interp = cp_interp - ap_interp;
  else
    ap_interp = [];
    cp_interp = [];
    bp_interp = [];
  end

  % Format ag and cg if dissolved data available, otherwise use a function derived independently
  if ~isempty(ac_g) & ~agcg_modelled
    % rebuild ag from exponential fit parameters
    ag_rebuilt = NaN(size(ac_g.ag));
    cg_rebuilt = NaN(size(ac_g.cg));
    % define exponential function
    expfun = @(p, xd) p(1) * exp(p(2) * (xd - 440));
    for i = 1:size(ac_g, 1)
      ag_rebuilt(i, :) = expfun([ac_g.y_intercp_fit_ag(i) ac_g.base_fit_ag(i)], lambda_a);
      cg_rebuilt(i, :) = expfun([ac_g.y_intercp_fit_cg(i) ac_g.base_fit_cg(i)], lambda_c);
    end
    % merge ag and rebuilt ag
    ag_merged = ac_g.ag;
    ag_merged(isnan(ag_merged)) = ag_rebuilt(isnan(ag_merged));
    ag_merged(all(isnan(ac_g.ag), 2), :) = NaN;
    % merge cg and rebuilt cg
    cg_merged = ac_g.cg;
    cg_merged(isnan(cg_merged)) = cg_rebuilt(isnan(cg_merged));
    cg_merged(all(isnan(ac_g.ag), 2), :) = NaN;
    fprintf('betap attenuation correction: filter event ag and cg reconstructed from ag and cg exponential fit parameters\n')

    % Reconstruct ag and cg based on correlation with fdom if available otherwise interpolate linearly between filter events
    if any(strcmp(ac_p.Properties.VariableNames, 'fdom')) && any(strcmp(ac_g.Properties.VariableNames, 'fdom'))
      [beta_total.fdom, flags.fdom_input_interpolated_linearly] = interp_on_bb_lambda(ac_p.dt, ac_p.fdom, [], beta_total.dt, []);
      [a_reg, c_reg] = regress_ac_fdom(ag_merged, cg_merged, ac_g.fdom);
      % interpolate slope on lambda_beta
      a_slope_interp = interp1(lambda_a, a_reg.slope, lambda_beta, 'linear');
      c_slope_interp = interp1(lambda_c, c_reg.slope, lambda_beta, 'linear');
      % linear interpolation of a, c, and fdom
      lin_ag = interp_on_bb_lambda(ac_g.dt, ag_merged, lambda_a, beta_total.dt, lambda_beta);
      lin_cg = interp_on_bb_lambda(ac_g.dt, cg_merged, lambda_c, beta_total.dt, lambda_beta);
      lin_fdom = interp_on_bb_lambda(ac_g.dt, ac_g.fdom, [], beta_total.dt, []);
      % Re-build ag and cg interpolated between filter events
      ag_interp = lin_ag + (beta_total.fdom - lin_fdom) .* a_slope_interp; 
      cg_interp = lin_cg + (beta_total.fdom - lin_fdom) .* c_slope_interp;
      
      % visProd3D(lambda_beta, beta_total.dt, beta_total.ag_interp, false, 'Wavelength', false, 2);
      % title('ag fdom interpolated')
      % visProd3D(lambda_beta, beta_total.dt, beta_total.cg_interp, false, 'Wavelength', false, 3);
      % title('cg fdom interpolated')
    
      % Fill missing values (start/end of analyzed period) with nearest filt_interp data
      ag_interp = fillmissing(ag_interp, 'linear', 'SamplePoints', beta_total.dt, 'EndValues', 'nearest');
      cg_interp = fillmissing(cg_interp, 'linear', 'SamplePoints', beta_total.dt, 'EndValues', 'nearest');

      fprintf('betap attenuation correction: ag and cg interpolated between filter events based on FDOM data.\n')
    else
      % interpolate ag on BB lambda
      [ag_interp, flags.ag_interpolated] = interp_on_bb_lambda(ac_g.dt, ag_merged, lambda_a, beta_total.dt, lambda_beta);
      % interpolate cg on BB lambda
      [cg_interp, flags.cg_interpolated] = interp_on_bb_lambda(ac_g.dt, cg_merged, lambda_c, beta_total.dt, lambda_beta);
      fprintf('betap attenuation correction: ag and cg interpolated linearly between filter events.\n')
    end
  elseif any(strcmp(beta_filt_interp.Properties.VariableNames, 'fdom')) && ~isempty(fdom_ag_correlation) & ~agcg_modelled
    % Case when the correlation between fdom and ag is computed from data from a benchtop spectrophotometer
    a_slope_interp = interp1(fdom_ag_correlation.wl, fdom_ag_correlation.slope, lambda_beta, 'linear')';
    a_intercept_interp = interp1(fdom_ag_correlation.wl, fdom_ag_correlation.intercept, lambda_beta, 'linear')';
    % fit exponential to slope and intercept on wavelength < 680 nm
    [y_intercp_slope_fit, base_slope_fit] = FitExp(lambda_beta(lambda_beta <= 680), ...
        a_slope_interp(lambda_beta <= 680)', ones(1, sum(lambda_beta <= 680)));
    [y_intercp_intercept_fit, base_intercept_fit] = FitExp(lambda_beta(lambda_beta <= 680), ...
        a_intercept_interp(lambda_beta <= 680)', ones(1, sum(lambda_beta <= 680)));
    % reconstruct slope at all lambda_beta
    expfun = @(p, xd) p(1) * exp(p(2) * (xd - 440));
    slope_rebuilt = expfun([y_intercp_slope_fit base_slope_fit], lambda_beta);
    intercept_rebuilt = expfun([y_intercp_intercept_fit base_intercept_fit], lambda_beta);

    % figure; hold on
    % plot(lambda_beta, a_slope_interp);
    % plot(lambda_beta, slope_rebuilt);
    % ylabel('a_g slope')
    % legend('slope', 'slope fit')
    % figure; hold on
    % plot(lambda_beta, a_intercept_interp);
    % plot(lambda_beta, intercept_rebuilt);
    % ylabel('a_g intercept')
    % legend('intercept', 'intercept fit')
    
    ag_interp = beta_filt_interp.fdom .* slope_rebuilt + intercept_rebuilt;
    cg_interp = ag_interp;
    fprintf('betap attenuation correction: ag and cg reconstructed based on fdom_ag_parameters and FDOM data.\n')
  elseif agcg_modelled
    % interpolate ag on BB lambda
    [ag_interp, flags.ag_interpolated] = interp_on_bb_lambda(ac_g.dt, ac_g.ag, lambda_a, beta_total.dt, lambda_beta);
    % interpolate cg on BB lambda
    [cg_interp, flags.cg_interpolated] = interp_on_bb_lambda(ac_g.dt, ac_g.cg, lambda_c, beta_total.dt, lambda_beta);
    fprintf('betap attenuation correction: ag and cg modelled from FCDOM and R(a)_FCDOM(lambda) used.\n')
  else
    ag_interp = [];
    cg_interp = [];
  end
  
  % visProd3D(lambda_beta, beta_filt_interp.dt, ag_interp, false, 'Wavelength', false, 72);
  % title('ag_interp')
  
  if isempty(ap_interp) && isempty(cp_interp) && isempty(bp_interp)
    % Just use beta_total - beta_filt as betap and flag data
    warning('betap attenuation correction: No ap data available, attenuation correction not applied.')
    flags.no_attenuation_correction = true(size(beta_total, 1), 1);
    betap_corr = betap_approx;
    bbp_corr = bbp_approx;
  % elseif isempty(ag_interp) && isempty(cg_interp) % if condition need to be changed for method and compute should be checked
  %   % Approximate ag if dissolved not available and apply simplified Doxaran's correction and flag data
  %   warning("betap attenuation correction: No ag data available, simplified Doxaran's correction applied.")
  %   flags.no_ag_attenuation_correction = true(size(beta_total, 1), 1);
  %   flags.simplified_Doxaran_correction = true(size(beta_total, 1), 1);
  %   % estimate ag440 from chlorophyll based on Morel and Maritorena, 2001
  %   ag440 = 0.2 .* 0.06 .* chl_interp .^ 0.65;
  %   % % estimate entire ag spectra based on exponential equation from Bricaud et al., 1981
  %   ag_interp = ag440 .* exp(-0.014 .* (lambda_beta - 440));
  %   % Adapted Doxaran's method to work with only ap and bp as input
  %   k_total = ag_interp + ap_interp + 3.3 .* bp_interp;
  %   % apply correction
  %   betap_corr = (beta_total.beta - beta_filt_interp.beta) ./ (1 - k_total * k_exp);
  %   bbp_corr = 2 * pi * Chi .* betap_corr;
  % 
  %   error('Not implemented')
  % 
  %   % visProd3D(lambda_beta, beta_total.dt, ag_interp, false, 'Wavelength', false, 71);
  %   % title('ag interp')
  % 
  %   % visProd3D(lambda_beta, beta_total.dt, k_total, false, 'Wavelength', false, 72);
  %   % title('k total')
  % 
  %   % betap_notcorr = beta_total.beta - beta_filt_interp.beta;
  %   % visProd3D(lambda_beta, beta_total.dt, betap_notcorr, false, 'Wavelength', false, 73);
  %   % title('betap approx')
  % 
  %   % visProd3D(lambda_beta, beta_total.dt, bbp_corr, false, 'Wavelength', false, 74);
  %   % title('bbp corr')
  % 
  %   % visProd3D(lambda_beta, beta_total.dt, (bbp_corr-betap_notcorr)./betap_notcorr*100, false, 'Wavelength', false, 75);
  %   % title('% change Doxaran')
  else % apply Zhang's sigma correction
    if isempty(ag_interp)
      % estimate ag440 from chlorophyll based on Morel and Maritorena, 2001 (appendix B)
      ag440 = 0.2 .* 0.06 .* chl_interp .^ 0.65;
      % estimate entire ag spectra based on exponential equation from Bricaud et al., 1981
      ag_interp = ag440 .* exp(-0.014 .* (lambda_beta - 440));
      cg_interp = ag_interp;
    end
    % correct beta_filt
    beta_filt_corr = beta_filt_interp.beta .* exp(ag_interp * k_exp);
    % Estimate c
    c_estimated = cp_interp + cg_interp;
    % Estimate single scattering albedo (w)
    w_estimated = bp_interp ./ c_estimated;
    % Force w to be between 0 and 1 when 1 < w_estimated < 1.005 (uncertainty of ACS == 0.005 m^-1)
    w_estimated(w_estimated > 1 & w_estimated <= 1.005) = 1;
    
    % % load Zhang's LUT (1st method)
    % if isfile(fullfile('packages','Zhang_hyper_bb_sigma_correction_1st_method.mat'))
    %   load(fullfile('packages','Zhang_hyper_bb_sigma_correction_1st_method.mat'), 'B', 'bbr', 'c', 'w')
    % end
    % load Zhang's LUT (3rd method 2026)
    if isfile(fullfile('packages', 'BB_sigma_correction_2026_LUT.mat'))
      load(fullfile('packages', 'BB_sigma_correction_2026_LUT.mat'),'LUT_sigmacorr_2026')
    else
      LUT_sigmacorr_2026 = [];
    end

    % prepare auxilliary data for new Zhang correction
    env = table();
    env.T = beta_filt_interp.t; % replace missing T by 25
    env.T(isnan(env.T)) = 22;
    env.S = beta_filt_interp.s; % replace missing S by 33
    env.S(isnan(env.S)) = 33;
    
    fprintf("Computing sigma correction ...\n")
    % Find sigma correction by iteration starting with betap approximation for backscattering ratio computation
    keep_going = true;
    bbp_approx_change = NaN(5, size(beta_total.beta, 2));
    bbp_approx_new = bbp_approx;

    i = 1;
    while keep_going
      fprintf('Iteration #%i ...', i)
      % Estimate backscattering ratio
      bbr_estimated = bbp_approx_new ./ bp_interp ;
      % isolate bbp_approx to compare after iteration
      bbp_approx_old = bbp_approx_new;

      % % Find sigma in Zhang LUT (1st method)
      % sigma_corr_old = interp3(bbr, w, c, B, bbr_estimated, w_estimated, c_estimated, "linear");
      
      % Compute sigma with Zhang 3rd method (VECTORIZED)
      [sigma_corr, ~, LUT_sigmacorr_2026] = get_bb_sigmacorr_2026_vect(lambda_beta, ...
        c_estimated, bbr_estimated, w_estimated, env, instrument_types, LUT_sigmacorr_2026);

      sigma_corr = fillmissing(sigma_corr, 'nearest', 'SamplePoints', beta_total.dt, 'MaxGap',minutes(40));
      betap_approx_new = beta_total.beta .* sigma_corr - beta_filt_corr;
      bbp_approx_new = 2 * pi * Chi .* betap_approx_new;
      bbp_approx_change(i, :) = mean((bbp_approx_new - bbp_approx_old) ./ bbp_approx_old, 1, 'omitmissing');
      if max(bbp_approx_change(i, :)) < 0.000001 || i == size(bbp_approx_change, 1)
        keep_going = false;
        betap_corr = betap_approx_new;
        bbp_corr = bbp_approx_new;
      else
        i = i + 1;
      end
      fprintf(' done\n')
    end
    fprintf("betap attenuation correction: Zhang's sigma correction applied.\n")

    % visProd3D(lambda_beta, beta_total.dt, ag_interp, false, 'Wavelength', false, 68);
    % zlabel('a_g input')
    % 
    % visProd3D(lambda_beta, beta_total.dt, cg_interp, false, 'Wavelength', false, 69);
    % zlabel('c_g input')
    % 
    % visProd3D(lambda_beta, beta_total.dt, ap_interp, false, 'Wavelength', false, 70);
    % zlabel('a_p input')
    % 
    % visProd3D(lambda_beta, beta_total.dt, cp_interp, false, 'Wavelength', false, 71);
    % zlabel('c_p input')
    % 
    % visProd3D(lambda_beta, beta_total.dt, bbp_approx, false, 'Wavelength', false, 72);
    % zlabel('bbp approximation (2 * pi * Chi .* (beta_{total} - beta_{filt}))'); title('uncorrected')
    % 
    % visProd3D(lambda_beta, beta_total.dt, sigma_corr, false, 'Wavelength', false, 73);
    % zlabel('sigma correction'); title("Zhang's lookup table")
    % 
    % visProd3D(lambda_beta, beta_total.dt, bbp_corr, false, 'Wavelength', false, 74);
    % zlabel('bbp corrected'); title("Zhang's lookup table")
    % 
    % visProd3D(lambda_beta, beta_total.dt, (bbp_approx_new - bbp_approx) ./ bbp_approx*100, false, 'Wavelength', false, 75);
    % zlabel('% change bbp'); title("Zhang's lookup table")

  end
end

function [data_interp, flag] = interp_on_bb_lambda(data_dt, data, lambda, bb_dt, lambda_bb)
  flag = ~ismember(bb_dt, data_dt);
  % interpolate ap on BB lambda
  if ~isempty(lambda)
    data_bbwl = interp1(lambda, data', lambda_bb, 'linear')';
  else
    data_bbwl = data;
  end
  % interpolate ag on beta_total datetime
  data_interp = interp1(data_dt, data_bbwl, bb_dt, 'linear');
  data_interp = fillmissing(data_interp, 'nearest', 'SamplePoints', bb_dt);
end

%% regression between a/c and fdom
function [a_reg, c_reg] = regress_ac_fdom(a_filt, c_filt, fdom)
  % robust linear regression between filt a&c and fdom
  a_reg = array2table(NaN(size(a_filt, 2), 5), 'VariableNames', {'slope', 'intercept', 'RMSE', 'nRMSE', 'R2'});
  c_reg = array2table(NaN(size(c_filt, 2), 5), 'VariableNames', {'slope', 'intercept', 'RMSE', 'nRMSE', 'R2'});
  for i = 1:size(a_filt, 2)
    if sum(~isnan(a_filt(:, i))) > 3
      % regress variable with filtered water absorption
      [stats.b, stats.stats] = robustfit(fdom, a_filt(:, i));
      [~, MSGID] = lastwarn();
      warning('off', MSGID)
      a_reg.slope(i) = stats.b(2);
      a_reg.intercept(i) = stats.b(1);
      a_reg.RMSE(i) = stats.stats.robust_s;
      % calculated normalize RMSE by average and R2
      id_nonan = all(~isnan([a_filt(:, i) fdom]), 2);
      a_reg.nRMSE(i) = a_reg.RMSE(i)/abs(mean(a_filt(id_nonan, i), 'omitnan'))*100;
      a_reg.R2(i) = corr(a_filt(id_nonan, i), stats.b(1) + stats.b(2) * fdom(id_nonan))^2;
    end
    if sum(~isnan(c_filt(:, i))) > 3
      % regress variable with filtered water attenuation
      [stats.b, stats.stats] = robustfit(fdom, c_filt(:, i));
      [~, MSGID] = lastwarn();
      warning('off', MSGID)
      c_reg.slope(i) = stats.b(2);
      c_reg.intercept(i) = stats.b(1);
      c_reg.RMSE(i) = stats.stats.robust_s;
      % calculated normalize RMSE by average and R2
      id_nonan = all(~isnan([c_filt(:, i) fdom]), 2);
      c_reg.nRMSE(i) = c_reg.RMSE(i)/abs(mean(c_filt(id_nonan, i), 'omitnan'))*100;
      c_reg.R2(i) = corr(c_filt(id_nonan, i), stats.b(1) + stats.b(2) * fdom(id_nonan))^2;
    end
  end
end

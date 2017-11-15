function fh = visFlagParams(params, data, varname, varindex)
  % Define constants
  ColorSet = lines(5);
  fh = fig(51); hold('on');
  set(groot,'defaultAxesColorOrder',lines(7));
  % Wendy's Test
  plot(data.dt, abs(data.([varname '_dtc_md'])(:,varindex) - data.([varname '_dtc_mn'])(:,varindex)) * params.avg_sensitivity, '.', 'Color', ColorSet(1,:));
  % Uncertainty too large (based on variance)
  plot(data.dt, data.([varname '_dtc_unc'])(:,varindex) * params.unc1_sensitivity, '.', 'Color', ColorSet(2,:));
  % Uncertainty too large (based on standard deviation)
  plot(data.dt, data.([varname '_dtc_se'])(:,varindex) * params.unc2_sensitivity, '.', 'Color', ColorSet(3,:));
  % Max Absolute and Relative Uncertainty
  foo = max(params.abs_uncertainty, params.rel_uncertainty * data.(varname)(:,varindex));
  plot(data.dt, foo, 'Color', ColorSet(4,:));
  % Smoother Absolute and Relative Uncertainty
  plot(data.dt, filtfilt(ones(params.smooth_threshold,1), params.smooth_threshold, foo), 'Color', ColorSet(5,:));
  legend('|Median-Mean|', 'Uncertainty', 'Standard Error', 'Reference', 'Smoothed Ref');
  datetick2_doy();
  set(datacursormode(fh),'UpdateFcn',@data_cursor_display_date);
end
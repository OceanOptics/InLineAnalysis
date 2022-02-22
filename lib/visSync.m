function fh = visSync(fth, dt, obs, y_label)
  fh = fig(30);
  yyaxis('right');
  plot(fth.dt, fth.swt, 'k');
%   plot(fth.dt, fth.beta, 'k');
  ylabel('0 = Total | 1 = Filtered');
  ax = gca; ax.YAxis(2).Color = 'k';
  yyaxis('left'); hold('on');
  plot(dt, obs, '.');
  if nargin > 3; ylabel(y_label); end
  datetick2_doy();
  set(datacursormode(fh),'UpdateFcn',@data_cursor_display_date);
%   hp = pan(); hp.Motion = 'horizontal';
end
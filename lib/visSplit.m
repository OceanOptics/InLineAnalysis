function fh = visSplit(fth, dt_tot, tot, dt_filt, filt, dt_reject, reject, y_label)
  % Define color set
  % blue total | red fitlered | yellow is transition | green is flagged | purple target
  ColorSet=lines(3);
  col_tot = ColorSet(1,:);
  col_filt = ColorSet(2,:);
  col_switch = ColorSet(3, :);
%   col_flag = ColorSet(4, :);
%   col_target = ColorSet(5, :);
  
  fh = fig(40); hold('on'); ColorSet=lines(3);
  yyaxis('right');
  h(4) = plot(fth.dt, fth.swt, 'k', 'LineWidth', 1);
  ylabel('0 = Total | 1 = Filtered');
  ax = gca; ax.YAxis(2).Color = 'k';
  yyaxis('left');
  h(1) = plot(dt_tot, tot, '.', 'Color', col_tot);
  if ~isempty(dt_filt);  h(2) = plot(dt_filt, filt, '.', 'Color', col_filt); end
  h(3) = plot(dt_reject, reject, '.', 'Color', col_switch);
  if nargin > 7; ylabel(y_label); end
  if ~isempty(dt_filt);
    legend(h, 'Total', 'Filtered', 'Reject', 'Switch','AutoUpdate','off');
  else
    legend(h([1,3,4]), 'Good', 'Reject', 'Switch','AutoUpdate','off');
  end
  % Classic
  datetick2_doy();
  pan('xon');
  set(datacursormode(fh),'UpdateFcn',@data_cursor_display_date);
  % Interactive
  %   user_selection = guiSelectOnTimeSeries(fh);
end
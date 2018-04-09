function fh = visFlag(raw_tot, raw_filt,...
                      bin_tot_good, bin_tot_suspect,...
                      bin_filt_good, bin_filt_bad,...
                      varname, varindex)
  % Define constants
  ColorSet = lines(5);
  % . raw | o bin 
  sha_raw = '.';
  sha_bin = 'o';
  % blue total | red fitlered | yellow is transition | green is flagged | purple target
  col_tot = ColorSet(1,:);
  col_filt = ColorSet(2,:);
%   col_switch = ColorSet(3, :);
  col_flag = ColorSet(4, :);
%   col_target = ColorSet(5, :);
  
  fh = fig(52); hold('on');

  % Plot second bin data
  scatter(raw_tot.dt, raw_tot.(varname)(:,varindex), sha_raw, 'MarkerEdgeColor', col_tot, 'MarkerEdgeAlpha', 0.5);
  if ~isempty(raw_filt); scatter(raw_filt.dt, raw_filt.(varname)(:,varindex), sha_raw, 'MarkerEdgeColor', col_filt, 'MarkerEdgeAlpha', 0.5); end
  
  % Plot tot binned data
  if ~isempty(bin_tot_good); plot(bin_tot_good.dt, bin_tot_good.(varname)(:,varindex), sha_bin, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_tot); end
  if ~isempty(bin_tot_suspect); plot(bin_tot_suspect.dt, bin_tot_suspect.(varname)(:,varindex), sha_bin, 'MarkerEdgeColor', col_flag, 'MarkerFaceColor', col_flag); end
%   plot(bb3_tot_flag.dt(sel_tot_target), bb3_tot_flag.beta(sel_tot_target,1), 'o', 'MarkerEdgeColor', col_target);
  
  % Plot filt binned data
  if ~isempty(bin_filt_good); plot(bin_filt_good.dt, bin_filt_good.(varname)(:,varindex), sha_bin, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_filt); end
  if ~isempty(bin_filt_bad); plot(bin_filt_bad.dt, bin_filt_bad.(varname)(:,varindex), sha_bin, 'MarkerEdgeColor', col_flag, 'MarkerFaceColor', col_flag); end
%   plot(bb3_filt_flag.dt(sel_filt_target), bb3_filt_flag.beta(sel_filt_target,1), 'o', 'MarkerEdgeColor', col_target);  legend('pass 2', 'flag 16 & 32');
  
  ylabel(varname);
  datetick2_doy();
  set(datacursormode(fh),'UpdateFcn',@data_cursor_display_date);
end
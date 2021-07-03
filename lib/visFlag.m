function fh = visFlag(raw_tot, raw_filt, bin_tot_good, bin_tot_suspect, bin_filt_good, ...
                      bin_filt_bad, varname, varindex, raw_bad, fooflow)
if nargin < 9; raw_bad = []; end

% Define constants
ColorSet = lines(5);
% . raw | o bin 
sha_raw = '.';
sha_bin = 'o';
% blue total | red fitlered | yellow is transition | green is flagged | purple target
col_tot = ColorSet(1,:);
col_filt = ColorSet(2,:);
col_switch = ColorSet(3, :);
col_flag = ColorSet(4, :);
%   col_target = ColorSet(5, :);

fh = fig(52); hold('on');
yyaxis('left');

% Plot second bin data
%   plot(raw_tot.dt, raw_tot.(varname{1})(:,varindex), sha_raw, 'Color', col_tot);
%   plot(raw_tot.dt, raw_tot.(varname{2})(:,varindex), sha_raw, 'Color', col_tot);
if ~isempty(raw_tot);  scatter(raw_tot.dt, raw_tot.(varname)(:,varindex), sha_raw, 'MarkerEdgeColor', col_tot, 'MarkerEdgeAlpha', 0.5); end
if ~isempty(raw_filt); scatter(raw_filt.dt, raw_filt.(varname)(:,varindex), sha_raw, 'MarkerEdgeColor', col_filt, 'MarkerEdgeAlpha', 0.5); end
if ~isempty(raw_bad); scatter(raw_bad.dt, raw_bad.(varname)(:,varindex), sha_raw, 'MarkerEdgeColor', col_switch, 'MarkerEdgeAlpha', 0.5); end

% Plot tot binned data
if ~isempty(bin_tot_good); plot(bin_tot_good.dt, bin_tot_good.(varname)(:,varindex), sha_bin, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_tot); end
if ~isempty(bin_tot_suspect); plot(bin_tot_suspect.dt, bin_tot_suspect.(varname)(:,varindex), sha_bin, 'MarkerEdgeColor', col_flag, 'MarkerFaceColor', col_flag); end
%   plot(bb3_tot_flag.dt(sel_tot_target), bb3_tot_flag.beta(sel_tot_target,1), 'o', 'MarkerEdgeColor', col_target);

% Plot filt binned data
if ~isempty(bin_filt_good); plot(bin_filt_good.dt, bin_filt_good.(varname)(:,varindex), sha_bin, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_filt); end
if ~isempty(bin_filt_bad); plot(bin_filt_bad.dt, bin_filt_bad.(varname)(:,varindex), sha_bin, 'MarkerEdgeColor', col_flag, 'MarkerFaceColor', col_flag); end
%   plot(bb3_filt_flag.dt(sel_filt_target), bb3_filt_flag.beta(sel_filt_target,1), 'o', 'MarkerEdgeColor', col_target);  legend('pass 2', 'flag 16 & 32');

% set y limit to binned data in priority
if isempty(bin_tot_good) && ~isempty(bin_filt_good)
  ylim([min(bin_filt_good.(varname)(:,varindex)) max(bin_filt_good.(varname)(:,varindex))]); 
elseif ~isempty(bin_tot_good) && isempty(bin_filt_good)
  ylim([min(bin_tot_good.(varname)(:,varindex)) max(bin_tot_good.(varname)(:,varindex))]);
elseif ~isempty(bin_tot_good) && ~isempty(bin_filt_good)
  ylim([min([bin_tot_good.(varname)(:,varindex); bin_filt_good.(varname)(:,varindex)]) ...
    max([bin_tot_good.(varname)(:,varindex); bin_filt_good.(varname)(:,varindex)])]);
end

ylabel(varname);
datetick2_doy();

% plot flow right y axis
if ~strcmp(varname,'par')
    yyaxis('right'); 
    if ~isempty(fooflow) && any(~isnan(fooflow.spd))
      plot(fooflow.dt, fooflow.spd,'Color',[0 0.8 0]); popo=gca; 
      ylim([0 max(fooflow.spd)+4]); ylabel('Flow rate (lpm)'); popo.YColor = [0 0.8 0]; % green
      ylabel('flow rate (lpm)');
    end
end

set(datacursormode(fh),'UpdateFcn',@data_cursor_display_date);
end
function fh = visFlag(raw_tot, raw_filt, bin_tot_good, bin_tot_suspect, bin_filt_good, ...
                      bin_filt_bad, varname, varindex, raw_bad, fooflow, spd_var, autoscale)
if nargin < 9; raw_bad = []; end
if nargin < 10; fooflow = []; end
if nargin < 11; spd_var = 'spd'; end
if nargin < 12; autoscale = false; end
% Define constants
ColorSet = lines(5);
% . raw | o bin 
sha_raw = '.';
sha_bin = 'o';
% blue total | red fitlered | yellow bad (switch buffer + QCed)
col_tot = ColorSet(1,:);
col_filt = ColorSet(2,:);
col_switch = ColorSet(3, :);
col_flag = ColorSet(4, :);
%   col_target = ColorSet(5, :);

fh = fig(52); hold('on');
yyaxis('left');

% check if data is loaded
if isempty(raw_filt) && isempty(bin_tot_good) && isempty(bin_filt_good)
  error('No data loaded, check previous step')
end

% Plot second bin data
sc = [];
leg = {};
if ~isempty(raw_bad)
  sc(end+1) = scatter(raw_bad.dt, raw_bad.(varname)(:,varindex), sha_raw, 'MarkerEdgeColor', col_switch, 'MarkerEdgeAlpha', 0.5);
  leg{end+1} = 'raw bad';
end
if ~isempty(raw_tot)
  sc(end+1) = scatter(raw_tot.dt, raw_tot.(varname)(:,varindex), sha_raw, 'MarkerEdgeColor', col_tot, 'MarkerEdgeAlpha', 0.5);
  leg{end+1} = 'raw total';
end
if ~isempty(raw_filt)
  sc(end+1) = scatter(raw_filt.dt, raw_filt.(varname)(:,varindex), sha_raw, 'MarkerEdgeColor', col_filt, 'MarkerEdgeAlpha', 0.5);
  leg{end+1} = 'raw filtered';
end

% Plot tot binned data
if ~isempty(bin_tot_good)
  sc(end+1) = plot(bin_tot_good.dt, bin_tot_good.(varname)(:,varindex), sha_bin, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_tot);
  leg{end+1} = 'binned total';
end
if ~isempty(bin_tot_suspect)
  sc(end+1) = plot(bin_tot_suspect.dt, bin_tot_suspect.(varname)(:,varindex), sha_bin, 'MarkerEdgeColor', col_flag, 'MarkerFaceColor', col_flag);
  leg{end+1} = 'binned total suspect';
end

% Plot filt binned data
if ~isempty(bin_filt_good)
  sc(end+1) = plot(bin_filt_good.dt, bin_filt_good.(varname)(:,varindex), sha_bin, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_filt);
  leg{end+1} = 'binned filtered';
end
if ~isempty(bin_filt_bad)
  sc(end+1) = plot(bin_filt_bad.dt, bin_filt_bad.(varname)(:,varindex), sha_bin, 'MarkerEdgeColor', col_flag, 'MarkerFaceColor', col_flag);
  leg{end+1} = 'binned filtered bad';
end

% set y limit to binned data in priority
if ~autoscale
  if isempty(bin_tot_good) && ~isempty(bin_filt_good) && ~isempty(raw_filt)
    foo_sc = [bin_filt_good.(varname)(:,varindex); raw_filt.(varname)(:,varindex)];
  elseif ~isempty(bin_tot_good) && isempty(bin_filt_good)
    foo_sc = bin_tot_good.(varname)(:,varindex);
  elseif ~isempty(bin_tot_good) && ~isempty(bin_filt_good) && ~isempty(raw_filt)
    foo_sc = [bin_tot_good.(varname)(:,varindex); bin_filt_good.(varname)(:,varindex); raw_filt.(varname)(:,varindex)];
  elseif ~isempty(bin_tot_good) && ~isempty(bin_filt_good)
    foo_sc = [bin_tot_good.(varname)(:,varindex); bin_filt_good.(varname)(:,varindex)];
  end
  ylim([min(foo_sc) max(foo_sc)]);
end

ylabel(varname);
datetick2_doy();

% plot flow right y axis
if ~strcmp(varname,'par')
    yyaxis('right');
    if ~isempty(fooflow)
      if ~any(strcmp(fooflow.Properties.VariableNames, spd_var))
        spd_var = fooflow.Properties.VariableNames{contains(fooflow.Properties.VariableNames, 'spd') & ...
          ~contains(fooflow.Properties.VariableNames, 'avg')};
      end
      if any(~isnan(fooflow.(spd_var)))
        sc(end+1) = plot(fooflow.dt, fooflow.(spd_var),'Color',[0 0.8 0]); popo=gca; 
        ylim([0 max(fooflow.(spd_var))+4]); ylabel('Flow rate (lpm)'); popo.YColor = [0 0.8 0]; % green
        ylabel('flow rate (lpm)');
        leg{end+1} = 'flow';
      end
    end
end
if ~isempty(leg)
  legend(sc, leg, 'FontSize', 13, 'AutoUpdate','off')
end

set(datacursormode(fh),'UpdateFcn',@data_cursor_display_date);
end
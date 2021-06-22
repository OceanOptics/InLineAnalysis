function [ user_selection_t, user_selection_f, user_selection_s, user_selection_d] = ...
  guiSelectOnTimeSeries( figure_handler )
%GUISELECTONTIMESERIES Graphical User Interface Function o interact with a
%timeseries and select parts of it.
%
% INPUT:
%   figure_handler <figure> get it from figure(1);
%
%
% Requires:
%   datetick2_doy
%   data_cursor_display_date


% Set figure active
figure(figure_handler);

% Set auto datetick
% datetick2_doy();
% Set Zoom
hz = zoom();
% hz.Motion = 'horizontal';
% AutoZoom at the beginning
%   zoom((x_lim(2) - x_lim(1)) * 24/2);
% Set Pan
hp = pan();
hp.Motion = 'horizontal';
% check data in Z axis and set proper data cursor function
all_obj = findobj(figure_handler, '-property', 'XData','YData','ZData');
if isempty(all_obj(end).ZData)
  set(datacursormode(figure_handler),'UpdateFcn',@data_cursor_display_date);
else
  set(datacursormode(figure_handler),'UpdateFcn',@data_cursor_display_date_y);
end
hc = datacursormode();
% Set Select (user select a chunk of data)
user_selection_t = [];
user_selection_f = [];
user_selection_s = [];
user_selection_d = [];

hold('on'); % to overlay area
% User Interaction
fprintf('User Interface Modes:\n');
fprintf('\t [z] unconstrained zoom\n');
fprintf('\t [h] horizontal zoom\n');
fprintf('\t [v] vertical zoom\n');
fprintf('\t [p] pan\n');
fprintf('\t [c] cursor\n');
fprintf('\t [t] select red (total/trash)\n');
fprintf('\t [f] select green (filtered)\n');
fprintf('\t [s] select time point (clean)\n');
fprintf('\t [d] delete selected point(s)\n');
fprintf('\t [q] save & quit\n');
k='';
% Master Loop Listenning to commands
while ~strcmp(k, 'q')
%   k = input('>','s');
  w = waitforbuttonpress();
  if w
    k = get(gcf, 'CurrentCharacter'); 
    switch k
%       case 'z'
%         if strcmp(hz.Enable, 'on')
%           hz.Enable = 'off';
%         else
%           hz.Enable = 'on';
%         end
      case 'z'
          hz.Motion = 'both';
        if strcmp(hz.Enable, 'on')
          hz.Enable = 'off';
        else
          hz.Enable = 'on';
        end
      case 'h'
          hz.Motion = 'horizontal';
        if strcmp(hz.Enable, 'on')
          hz.Enable = 'off';
        else
          hz.Enable = 'on';
        end
      case 'v'
          hz.Motion = 'vertical';
        if strcmp(hz.Enable, 'on')
          hz.Enable = 'off';
        else
          hz.Enable = 'on';
        end
      case 'p'
        if strcmp(hp.Enable, 'on')
          hp.Enable = 'off';
        else
          hp.Enable = 'on';
        end
      case 'c'
        if strcmp(hc.Enable, 'on')
          hc.Enable = 'off';
        else
          hc.Enable = 'on';
        end
      case 't'
        foo = ginput(2);
        x_lim = xlim(); y_lim = ylim();
        if ~isdatetime(x_lim)
            dt_sel_t = foo(:,1)';
        else
            ax = gca; xdate = num2ruler(foo,ax.XAxis); dt_sel_t = xdate(:,1)';
        end
        user_selection_t = [user_selection_t; dt_sel_t];
        fprintf('%s - %s\n', datestr(dt_sel_t(1)), datestr(dt_sel_t(2)));
        % Plot selected area
        area(dt_sel_t,[y_lim(2), y_lim(2)], y_lim(1), 'FaceColor', [1 0.3 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

      case 'f'
        foo = ginput(2);
        x_lim = xlim(); y_lim = ylim();
        if ~isdatetime(x_lim)
            dt_sel_f = foo(:,1)';
        else
            ax = gca; xdate = num2ruler(foo,ax.XAxis); dt_sel_f = xdate(:,1)';
        end
        user_selection_f = [user_selection_f; dt_sel_f];
        fprintf('%s - %s\n', datestr(dt_sel_f(1)), datestr(dt_sel_f(2)));
        % Plot selected area
        area(dt_sel_f,[y_lim(2), y_lim(2)], y_lim(1), 'FaceColor', [0.3 1 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

      case 's'
        foo = ginput(1);
        x_lim = xlim(); y_lim = ylim();
        if ~isdatetime(x_lim)
            dt_sel_s = foo(:,1)';
        else
            ax = gca; xdate = num2ruler(foo,ax.XAxis); dt_sel_s = xdate(:,1)';
        end
        user_selection_s = [user_selection_s; dt_sel_s];
        fprintf('%s\n', datestr(dt_sel_s));
        % Plot selected area
        plot([dt_sel_s dt_sel_s], [y_lim(1), y_lim(2)], '-', 'Color', [1 0.2 0.2]);
      case 'd'
        % delete from plot directly
        datatip = struct2table(getCursorInfo(hc));
        if isempty(all_obj(end).ZData)
          user_selection_d = [user_selection_d; datatip.Position(:,1)];
          fprintf('%s\n', strjoin(cellstr(datestr(datatip.Position(:,1))), ' | '));
%           dt_dat = {all_obj.YData};
%           tomodify = {all_obj.ZData};
          for todel = 1:size(datatip.Position(:,1), 1)
            all_obj(end).YData(all_obj(end).XData == datatip.Position(todel,1)) = NaN;
%             tomodify{end}(dt_dat{end} == datatip.Position(todel,2), :) = NaN;
          end
%           set(all_obj, 'ZData', tomodify{end});
%           set(d, 'ZData', tomodify{end});
        else
          user_selection_d = [user_selection_d; datatip.Position(:,2)];
          fprintf('%s\n', strjoin(cellstr(datestr(datatip.Position(:,2))), ' | '));
%           dt_dat = {all_obj.YData};
%           tomodify = {all_obj.ZData};
          for todel = 1:size(datatip.Position(:,2), 1)
            all_obj(end).ZData(all_obj(end).YData == datatip.Position(todel,2), :) = NaN;
%             all_obj(end).ZData = fillmissing(all_obj(end).ZData, 'linear', 'SamplePoints', all_obj(end).YData);
%             tomodify{end}(dt_dat{end} == datatip.Position(todel,2), :) = NaN;
          end
%           set(all_obj, 'ZData', tomodify{end});
%           set(d, 'ZData', tomodify{end});
        end
        refresh(figure_handler)
      otherwise
        if ~strcmp(k,'q')
          fprintf('?\n');
        end
    end
  end
end
fprintf('Saved selection\n');

end


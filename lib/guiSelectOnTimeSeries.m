function [ user_selection_t, user_selection_f ] = guiSelectOnTimeSeries( figure_handler )
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
% Set Data Cursor
set(datacursormode(figure_handler),'UpdateFcn',@data_cursor_display_date);
hc = datacursormode();
% Set Select (user select a chunk of data)
user_selection_t = [];
user_selection_f = [];
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
      otherwise
        if ~strcmp(k,'q')
          fprintf('?\n');
        end
    end
  end
end
fprintf('Saved selection\n');

end


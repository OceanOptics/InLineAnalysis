function [ user_selection ] = guiSelectOnTimeSeries( figure_handler )
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
datetick2_doy();
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
user_selection = [];
hold('on'); % to overlay area
% User Interaction
fprintf('User Interface Modes:\n');
fprintf('\t [z] zoom\n');
fprintf('\t [p] pan\n');
fprintf('\t [c] cursor\n');
fprintf('\t [s] select\n');
fprintf('\t [q] save & quit\n');
k='';
% Master Loop Listenning to commands
while ~strcmp(k, 'q')
%   k = input('>','s');
  w = waitforbuttonpress();
  if w
    k = get(gcf, 'CurrentCharacter'); 
    switch k
      case 'z'
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
      case 's'
        foo = ginput(2); dt_sel = foo(:,1)';
        user_selection = [user_selection; dt_sel];
        fprintf('%s - %s\n', datestr(dt_sel(1)), datestr(dt_sel(2)));
        % Plot selected area
        x_lim = xlim(); y_lim = ylim();
        area(dt_sel,[y_lim(2), y_lim(2)], y_lim(1), 'FaceColor', [0.3 0.3 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
      otherwise
        if ~strcmp(k,'q')
          fprintf('?\n');
        end
    end
  end
end
fprintf('Saved selection\n');

end


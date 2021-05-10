function txt = data_cursor_display_date_y(empt,event_obj)
% DATA_CURSOR_DISPLAY_DATE 
%   Display date & time of cursor on timeseries figure
%   Add the following code after ploting your figure:
%--------------------------------------------------------------------
% set(datacursormode(figure(2)),'UpdateFcn',@data_cursor_display_date);
%--------------------------------------------------------------------

pos = get(event_obj,'Position');
txt = {datestr(pos(2), 'dd-mmm-yyyy HH:MM:SS.FFF'),['x: ',num2str(pos(1))],['z: ',num2str(pos(3))]};
end
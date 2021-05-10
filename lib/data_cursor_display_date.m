function txt = data_cursor_display_date(empt,event_obj)
% DATA_CURSOR_DISPLAY_DATE 
%   Display date & time of cursor on timeseries figure
%   Add the following code after ploting your figure:
%--------------------------------------------------------------------
% set(datacursormode(figure(2)),'UpdateFcn',@data_cursor_display_date);
%--------------------------------------------------------------------

pos = get(event_obj,'Position');
txt = {datestr(pos(1),'dd-mmm-yyyy HH:MM:SS.FFF'),['y: ',num2str(pos(2))]};
end
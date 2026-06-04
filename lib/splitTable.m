function [ data_total, data_filt, data_reject, data_lost ] = splitTable( ref, data, buffer, mode, constants, ref_varname, verbose )
%SPLITTABLE Split time series into total and filtered based on switch position
%   provided in the reference table.
% Note: the reference and data tables must be synchronised, the delay for the
% water between the reference and the data should be accounted for
% WARNING: only the intersect between the reference table and the data table it
% processed, any data recorded in data that does not has a timestamp in
% reference will be excluded from the output.
% WARNING: dt from the reference and data must be truncated to the second
% 
% INPUT:
%   - ref <Table>: must contain fields:
%         dt <Nx1 datenum> date and time precise to the second
%         ref_varname <Nx1 double> position of switch: 1 -> total; 0 -> filtered (see ref_var)
%   - data <Table>: table of data to split between two events,
%       must contain field:
%         dt <Nx1 datenum> date and time precise to the second
%   - buffer <1x2 double> in seconds
%         buffer(1): delay after switch event from total to filtered to get good data
%         buffer(2): delay after switch event from filtered to total to get good data
%   - mode <'split', 'rmBuffer'>
%         split: split data in total and filtered (remove buffered data)
%         rmBuffer: remove buffered data (all good data goes in data_total)
%   - constants <struct('SWITCH_FILTERED', 1, 'SWITCH_TOTAL', 0)>
%         change values assigned to position of switch
%   - ref_varname <char>: name of reference variable to use in ref table.
%
% OUTPUT:
%   - data_total <Table>: all total data from table data
%       contains all fields of table data
%   - data_filtered <Table>: all filtered data from table data
%       contains all fields of table data
%   - data_reject <Table>: all buffered data from table data
%   - data_lost <Table>: data from the data table that has no timestamp
%       matching in the reference table (does not include buffer data)
%
% author: Nils
% created: Sept 8, 2017

if nargin < 6; verbose = false; end
% Set switch position
if nargin < 5
  % Assume most recent FlowControl software
  SWITCH_FILTERED = 1;
  SWITCH_TOTAL = 0;
else
  SWITCH_FILTERED = constants.SWITCH_FILTERED;
  SWITCH_TOTAL = constants.SWITCH_TOTAL;
end
if isnumeric(buffer)
  buffer = seconds(buffer);
end

if verbose
  switch mode
    case 'split'
      fprintf('Splitting Filter|Total %s ... ', inputname(2));
    case 'rmBuffer'
      fprintf('Removing buffer %s ... ', inputname(2));
    otherwise
      error('Unknown mode');
  end
end

% Sort Rows of data by dt
data = sortrows(data, 'dt');

% Get intersect of reference and data
% Using intersect (not really good results, assume perfect logging of ref and data)
% [dt, i_ref, i_data] = intersect(ref.dt, data.dt);
% swt = ref.(ref_varname)(i_ref);
% d = data(i_data,:);

% Get intersect of reference and data
% Improved method by interpolation (assume switch has only zeros and ones)
[urefdt, i] = unique(ref.dt);
if islogical(ref.(ref_varname)); ref.(ref_varname) = double(ref.(ref_varname)); end
iref.swt = interp1(urefdt, ref.(ref_varname)(i), data.dt);
i_data = iref.swt == 0 | iref.swt == 1;
dt = data.dt(i_data);
swt = iref.swt(i_data);
d = data(i_data,:);

% Init few variables
n = size(swt,1);
flag_buffer = false(n,1);
sampling_frequency = 1 / seconds(median(diff(data.dt), 'omitnan')); % Hertz
step_buffer = round(seconds(buffer) * sampling_frequency);
% Find switch events from filtered to total
sel = find(swt(1:end-1) == SWITCH_FILTERED & swt(2:end) == SWITCH_TOTAL);
% Flag for buffer time after event
for i=1:size(sel, 1)
  flag_buffer(sel(i) - 1 + find(dt(sel(i):min(sel(i)+step_buffer(2)+10,n)) <= (dt(sel(i)) + seconds(buffer(2))/3600/24))) = true;
end
% % exact but much slower
% % Find switch events from filtered to total
% filt2tot = dt(swt(1:end-1) == SWITCH_FILTERED & swt(2:end) == SWITCH_TOTAL);
% % Flag for buffer time after event
% for i=progress(1:size(filt2tot, 1))
%   id = dt > filt2tot(i) & dt <= filt2tot(i) + buffer(2);
%   flag_buffer(id) = true;
% end

if strcmp(mode, 'split')
  % Query reference table
  sel_tot = swt == SWITCH_TOTAL & ~flag_buffer;
  % Return what we were asked
  data_total = d(sel_tot,:);
end

if nargout > 1 || strcmp(mode, 'rmBuffer')
  % Find switch events from total to filtered
  sel = find(swt(1:end-1) == SWITCH_TOTAL & swt(2:end) == SWITCH_FILTERED);
  % Flag for buffer time before event
  for i=1:size(sel, 1)
    flag_buffer(sel(i) - 1 + find(dt(sel(i):min(sel(i)+step_buffer(1)+10,n)) <= (dt(sel(i)) + seconds(buffer(1))/3600/24))) = true;
  end
  % % exact but much slower
  % % Find switch events from total to filtered
  % tot2filt = dt(swt(1:end-1) == SWITCH_TOTAL & swt(2:end) == SWITCH_FILTERED);
  % % Flag for buffer time before event
  % for i=progress(1:size(tot2filt, 1))
  %   id = dt > tot2filt(i) & dt <= tot2filt(i) + buffer(1);
  %   flag_buffer(id) = true;
  % end
  
  switch mode
    case 'split'
      % Query reference table
      sel_filt = swt == SWITCH_FILTERED & ~flag_buffer;
      % Return what we were asked
      data_filt = d(sel_filt,:);
    case 'rmBuffer'
      % Query reference table
      sel = ~flag_buffer;
      % Return what we were asked
      data_total = d(sel,:);
      data_filt = [];
  end
end

if nargout > 2
  data_reject = d(flag_buffer,:);
end

if nargout > 3
  % I don't want to work with that
  sel_lost = true(size(data,1),1);
  sel_lost(i_data) = false;
  data_lost = data(sel_lost,:);
end

if verbose; fprintf('Done\n'); end
% Enjoy the sea
end


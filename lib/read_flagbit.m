function flag = read_flagbit(flag_bit, instrument)
% Read flag bit of Inline and organize in a clean table
%
% Author: Guillaume Bourdin
% Date: 2021-05-28
%%
if ~isnumeric(flag_bit) || size(flag_bit,2) > 1
  error('flag_bit must be vector of double')
end

% get instrument flag info
FlagInfo = InlineFlagInfo(instrument);

% convert bit to table of logical
flag = logical(bitget(repmat(flag_bit, 1, size(FlagInfo,1)), ...
  repmat((FlagInfo.bit+1)', size(flag_bit,1), 1)));
flag = array2table(flag, 'VariableNames', FlagInfo.name);
flag.Properties.VariableDescriptions = FlagInfo.description;


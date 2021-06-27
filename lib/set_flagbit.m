function flag_bit = set_flagbit(flag)
% Set bit from matrix of double 0 and 1 or matrix of logical
%
% Author: Guillaume Bourdin
% Date: 2021-05-28
%%
if istable(flag)
  flag = table2array(flag);
end
if ~islogical(flag) && max(flag(:)) <= 1
  flag = logical(flag);
elseif max(flag(:)) > 1
  error('flag must be matrix of double 0 and 1 or matrix of logical')
end

p = repmat(1:size(flag, 2), size(flag, 1), 1);
flag_bit = zeros(size(flag));
flag_bit(flag) = 2.^(p(flag)-1);
flag_bit = sum(flag_bit, 2);

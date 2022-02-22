function [sel_good, sel_suspect] = selFlag(flag, varname, n_index)
  % Define constants 
%   FLAG_NO_DATA = 1;
%   FLAG_GOOD = 2;
%   FLAG_SUSPECT_MANUAL = 3;
%   FLAG_BAD_MANUAL = 4;
%   FLAG_VALUE_TOO_LARGE = 5;
%   FLAG_VARIANCE_TOO_LARGE = 6;
%   FLAG_AVG_TEST_FAIL = 8;
  FLAG_UNC_TEST_FAIL = 9;
  FLAG_SE_TEST_FAIL = 10;
  % Query
  sel_suspect = sum(bitget(flag.([varname '_flag']), FLAG_UNC_TEST_FAIL),2) >= n_index |...
                sum(bitget(flag.([varname '_flag']), FLAG_SE_TEST_FAIL),2) >= n_index;
              %sum(bitget(flag.([varname '_flag']), FLAG_VALUE_TOO_LARGE),2) >= 1 | ...
  sel_good = ~sel_suspect;
end
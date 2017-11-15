function S = datestr_doy(D,varargin)
%DATESTR String representation of date.
%   S = DATESTR_DOY(V) converts one or more date vectors V to date strings S.
%   Input V must be an M-by-6 matrix containing M full (six-element) date
%   vectors. Each element of V must be a positive double-precision number.
%   DATESTR_DOY returns a column vector of M date strings, where M is the total
%   number of date vectors in V. 
%
%   S = DATESTR_DOY(N) converts one or more serial date numbers N to date
%   strings S. Input argument N can be a scalar, vector, or
%   multidimensional array of positive double-precision numbers. DATESTR_DOY
%   returns a column vector of M date strings, where M is the total number
%   of date numbers in N. 
%
%   S = DATESTRN_DOY(D, F) converts one or more date vectors, serial date
%   numbers, or date strings D into the same number of date strings S.
%   Input argument F is a format number or string that determines the
%   format of the date string output. Valid values for F are given in Table
%   1, below. Input F may also contain a free-form date format string
%   consisting of format tokens as shown in Table 2, below. 
%
%   Date strings with 2-character years are interpreted to be within the
%   100 years centered around the current year. 
%
%   S = DATESTR_DOY(S1, F, P) converts date string S1 to date string S,
%   applying format F to the output string, and using pivot year P as the
%   starting year of the 100-year range in which a two-character year
%   resides. The default pivot year is the current year minus 50 years.
%   F = -1 uses the default format.
%
%	S = DATESTR_DOY(...,'local') returns the string in a localized format. The
%	default (which can be called with 'en_US') is US English. This argument 
%	must come last in the argument sequence.
%
%	Note:  The vectorized calling syntax can offer significant performance
%	improvement for large arrays.
%
%	Table 1: Standard MATLAB Date format definitions
%
%   Number           String                   Example
%   ===========================================================================
%      0             'dd-mmm-yyyy HH:MM:SS'   01-Mar-2000 15:45:17 
%      1             'dd-mmm-yyyy'            01-Mar-2000  
%      2             'mm/dd/yy'               03/01/00     
%      3             'mmm'                    Mar          
%      4             'm'                      M            
%      5             'mm'                     03            
%      6             'mm/dd'                  03/01        
%      7             'dd'                     01            
%      8             'ddd'                    Wed          
%      9             'd'                      W            
%     10             'yyyy'                   2000         
%     11             'yy'                     00           
%     12             'mmmyy'                  Mar00        
%     13             'HH:MM:SS'               15:45:17     
%     14             'HH:MM:SS PM'             3:45:17 PM  
%     15             'HH:MM'                  15:45        
%     16             'HH:MM PM'                3:45 PM     
%     17             'QQ-YY'                  Q1-96        
%     18             'QQ'                     Q1           
%     19             'dd/mm'                  01/03        
%     20             'dd/mm/yy'               01/03/00     
%     21             'mmm.dd,yyyy HH:MM:SS'   Mar.01,2000 15:45:17 
%     22             'mmm.dd,yyyy'            Mar.01,2000  
%     23             'mm/dd/yyyy'             03/01/2000 
%     24             'dd/mm/yyyy'             01/03/2000 
%     25             'yy/mm/dd'               00/03/01 
%     26             'yyyy/mm/dd'             2000/03/01 
%     27             'QQ-YYYY'                Q1-1996        
%     28             'mmmyyyy'                Mar2000        
%     29 (ISO 8601)  'yyyy-mm-dd'             2000-03-01
%     30 (ISO 8601)  'yyyymmddTHHMMSS'        20000301T154517 
%     31             'yyyy-mm-dd HH:MM:SS'    2000-03-01 15:45:17 
%     32             'yyyy/doy'               2000/365
%     33             'yy/doy'                 00/365
%     34             'yy/doy.HH:MM'           00/365.15:45
%     35             'doy'                    365
%     36             'doy.HH:MM'              365.15:45


%
%   Table 2: Free-form date format symbols
%   
%   Symbol  Interpretation of format symbol
%   ===========================================================================
%   yyyy    full year, e.g. 1990, 2000, 2002
%   yy      partial year, e.g. 90, 00, 02
%   mmmm    full name of the month, according to the calendar locale, e.g.
%           "March", "April" in the UK and USA English locales. 
%   mmm     first three letters of the month, according to the calendar 
%           locale, e.g. "Mar", "Apr" in the UK and USA English locales. 
%   mm      numeric month of year, padded with leading zeros, e.g. ../03/..
%           or ../12/.. 
%   m       capitalized first letter of the month, according to the
%           calendar locale; for backwards compatibility. 
%   dddd    full name of the weekday, according to the calendar locale, eg.
%           "Monday", "Tueday", for the UK and USA calendar locales. 
%   ddd     first three letters of the weekday, according to the calendar
%           locale, eg. "Mon", "Tue", for the UK and USA calendar locales. 
%   dd      numeric day of the month, padded with leading zeros, e.g. 
%           05/../.. or 20/../.. 
%   d       capitalised first letter of the weekday; for backwards 
%           compatibility
%   n       day of the year
%   HH      hour of the day, according to the time format. In case the time
%           format AM | PM is set, HH does not pad with leading zeros. In 
%           case AM | PM is not set, display the hour of the day, padded 
%           with leading zeros. e.g 10:20 PM, which is equivalent to 22:20; 
%           9:00 AM, which is equivalent to 09:00.
%   MM      minutes of the hour, padded with leading zeros, e.g. 10:15, 
%           10:05, 10:05 AM.
%   SS      second of the minute, padded with leading zeros, e.g. 10:15:30,
%           10:05:30, 10:05:30 AM. 
%   FFF     milliseconds field, padded with leading zeros, e.g.
%           10:15:30.015.
%   PM      set the time format as time of morning or time of afternoon. AM 
%           or PM is appended to the date string, as appropriate. 
%
%   Examples:
%	DATESTR_DOY(now) returns '24-Jan-2003 11:58:15' for that particular date,
%	on an US English locale DATESTRN(now,2) returns 01/24/03, the same as
%	for DATESTR_DOY(now,'mm/dd/yy') DATESTRN(now,'dd.mm.yyyy') returns
%	24.01.2003 To convert a non-standard date form into a standard MATLAB
%	dateform, first convert the non-standard date form to a date number,
%	using DATENUM, for example, 
%	DATESTR_DOY(DATENUM('24.01.2003','dd.mm.yyyy'),2) returns 01/24/03.
%
%	See also DATE, DATENUM, DATEVEC, DATETICK, DATESTR, FORMATDATE_DOY

%	Copyright 1984-2006 The MathWorks, Inc.
%	$Revision: 1.32.4.14 $  $Date: 2007/05/23 18:55:05 $

%   datestr.m modified by David Murphy, 11/20/09
%   Added day of year formats 
%   Called from datetick_doy.m

%==============================================================================
% handle input arguments
if (nargin<1) || (nargin>4)
    error('MATLAB:datestr:Nargin',nargchk(1,4,nargin));
end
last = nargin - 1;
islocal = 0;
if last > 0 && ischar(varargin{end})
    if strcmpi(varargin{end}, 'local')
        islocal = 1;
        last = last - 1;
    elseif strcmpi(varargin{end},'en_us')
        islocal = 0;
        last = last - 1;
    end
end

if last > 2
    error('MATLAB:datestr:Nargin','Too many arguments.');    
elseif last >= 1
    dateform = varargin{1};
    if last == 2
        pivotyear = varargin{2};
    end
end

isdatestr = ~isnumeric(D);
if last > 0
    if ~ischar(dateform);
        % lookup date form string on index
        dateformstr = getdateform(dateform);
    else
        dateformstr = dateform;
    end
else
    dateformstr = '';
end

if last == 2 && ischar(pivotyear)
    error('MATLAB:datestr:InputClass', 'Pivot year must be a number.');
end

% Convert strings and clock vectors to date numbers.
try
    if isdatestr || (size(D,2)==6 && all(all(D(:,1:5)==fix(D(:,1:5)))) &&...
        all(abs(sum(D,2)-2000)<500)) 
        if last <= 1 || ~isdatestr  %not a datestring or no pivot year.
            dtnumber = datenum(D);
        else %datestring and pivot year were passed in.
            dtnumber = datenum(D,pivotyear);
        end
    else %datenum was passed in
        dtnumber = D;
    end
catch
    error('MATLAB:datestr:ConvertToDateNumber',...
        'Cannot convert input into specified date string.\n%s.',lasterr);
end

% Determine format if none specified.  If all the times are zero,
% display using date only.  If all dates are all zero display time only.
% Otherwise display both time and date.
dtnumber = dtnumber(:);
if (last < 1) || (isnumeric(dateform) && (dateform == -1))
   if all(floor(dtnumber)==dtnumber)
      dateformstr = getdateform(1);
   elseif all(floor(dtnumber)==0)
      dateformstr = getdateform(16);
   else
      dateformstr = getdateform(0);
   end
end 

% Handle the empty case properly.  Return an empty which is the same
% length of the string that is normally returned for each dateform.
if isempty(dtnumber)
   S= reshape('', 0, length(dateformstr)); 
   return;
end

try
    if ~isfinite(dtnumber)
        %Don't bother to go through mex file, since datestr can not handle
        %non-finite dates.
        error('MATLAB:datestr', 'Date number out of range.');
    end
    % Obtain components using mex file
    [y,mo,d,h,minute,s] = datevecmx(dtnumber,true);  mo(mo==0) = 1;
catch
    err = lasterror;
    err.identifier = 'MATLAB:datestr:ConvertDateNumber';
    err.message = sprintf('%s\n%s\n%s',...
                          ['DATESTR failed converting date'...
                           ' number to date vector.'],...
                          err.message);
    rethrow(err);
end

% format date according to data format template
try
    S = formatdate_doy([y,mo,d,h,minute,s],dateformstr,islocal);
catch
    err=lasterror;
    err.message=sprintf('%s Format string %s.',...
                        err.message,dateformstr);
    rethrow(err);
end 
S = char(S);

%==============================================================================
function [formatstr] = getdateform(dateform)
% Determine date format string from date format index.
    switch dateform
        case -1, formatstr = 'dd-mmm-yyyy HH:MM:SS';
        case 0,  formatstr = 'dd-mmm-yyyy HH:MM:SS';
        case 1,  formatstr = 'dd-mmm-yyyy';
        case 2,  formatstr = 'mm/dd/yy';
        case 3,  formatstr = 'mmm';
        case 4,  formatstr = 'm';
        case 5,  formatstr = 'mm';
        case 6,  formatstr = 'mm/dd';
        case 7,  formatstr = 'dd';
        case 8,  formatstr = 'ddd';
        case 9,  formatstr = 'd';
        case 10, formatstr = 'yyyy';
        case 11, formatstr = 'yy';
        case 12, formatstr = 'mmmyy';
        case 13, formatstr = 'HH:MM:SS';
        case 14, formatstr = 'HH:MM:SS PM';
        case 15, formatstr = 'HH:MM';
        case 16, formatstr = 'HH:MM PM';
        case 17, formatstr = 'QQ-YY';
        case 18, formatstr = 'QQ';
        case 19, formatstr = 'dd/mm';
        case 20, formatstr = 'dd/mm/yy';
        case 21, formatstr = 'mmm.dd,yyyy HH:MM:SS';
        case 22, formatstr = 'mmm.dd,yyyy';
        case 23, formatstr = 'mm/dd/yyyy';
        case 24, formatstr = 'dd/mm/yyyy';
        case 25, formatstr = 'yy/mm/dd';
        case 26, formatstr = 'yyyy/mm/dd';
        case 27, formatstr = 'QQ-YYYY';
        case 28, formatstr = 'mmmyyyy'; 
        case 29, formatstr = 'yyyy-mm-dd';
        case 30, formatstr = 'yyyymmddHHMMSS';
        case 31, formatstr = 'yyyy-mm-dd HH:MM:SS';
        case 32, formatstr = 'yyyy/nnn';
        case 33, formatstr = 'yy/nnn';
        case 34, formatstr = 'yy/nnn.HH:MM';
        case 35, formatstr = 'nnn';
        case 36, formatstr = 'nnn.HH:MM'; 

        otherwise
            error('MATLAB:datestr:DateNumber',...
                'Unknown date format number: %s', dateform);
    end 

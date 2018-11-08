% Manually Quality Check a project
%	Script meant to be run after bin once the data was binned on a remote system
% author: Nils Haentjens
% created: Nov 4, 2018

% Intruction to manually QC the ACS data from a project
%   1. complete the section 0. Configuration
%   2. run the entire file (cmd + return)
%   3. wait for data to load and figure to display
%   4. look for instruction on the console and select bad section of data
%       When doing so make sure that if during extended period of time there is no filtered data
%       please remove those section especially when total data is collected and the instrument
%       cleanned before the next filtered seawater section
%   5. Once done the software will produice a folder named UI with a few files containing your work.
%       Send those to Nils



%% 0. Configuration
% Set the path to the configuration file of th project to run
path_to_configuration_file = 'cfg/default_cfg.m';

% Set path to data
path_to_data = '/Users/nils/Data/NAAMES/NAAMES3/';

% Days to QC
%   Set to NaN to QC the entire project
%   If you're computer is not able to load everything or you want to run day by day
%       set the date of beginning and end as folow
%       days_to_run = datenum(2015,11,10):datenum(2015,11,12)
days_to_qc = NaN;

% Instruments to QC
%   Set a cell array of instruments to QC, make sure you have the data for all of them
%   instruments_to_qc = {'TSG', 'ACS'};
instruments_to_qc = {'ACS'};

% Load raw data
%   Helps to see the trend as get full resolution in addition to minute binned data
%   Will require more RAM on your computer as 60 to 240 times more data to display
load_original_data = False;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  No modifications needed below here  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. Load configuration
global PATH_ROOT;
PATH_ROOT = path_to_data;
ila = InLineAnalysis(path_to_configuration_file);
% Overwrite configuration parameters
if ~isnan(days_to_qc); ila.cfg.days2run = days_to_qc; end
ila.cfg.instruments2run = instruments_to_qc;

%% 2. Load raw data [optional]
% Loads raw data that is already synchronized and split in total, filtered,
if load_raw_data; ila.Read('raw'); end

%% 3. Load binned data
ila.Read('bin');

%% 4. Skip automatic QC 
cfg.process.flag.skip = {'FTH', 'TSG', 'BB3', 'WSCD', 'ACS'};
ila.Flag();

%% 5. Manually QC the data
ila.cfg.qc.mode='ui';
ila.QC();

fprintf("Thank you !\n");
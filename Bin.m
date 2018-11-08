% Bin data
% 	Load raw data and run binning process
%	Script designed to run on super-computers
% author: Nils Haentjens
% created: Nov 5, 2018
% based on EXPORTS file


% Load InLineAnalysis
ila = InLineAnalysis('cfg/default_cfg.m');

% Quick Cfg update
ila.cfg.days2run = days2run;
ila.cfg.instruments2run = {'ACS'};

% Binning Method
% cfg.process.bin.method = '4flag'; % Method to use to flag automatically
cfg.process.bin.method = 'SB_IN_PRCTL'; % Faster method outputting only requirements for SeaBASS (not compatible with automatic flagging)
% cfg.process.bin.method = 'SB_ALL'; % Faster method outputting only requirements for SeaBASS (not compatible with automatic flagging)
% Binning mode
%     does not affect the outcome of the data but just the way the data is presented to the computer
% cfg.process.bin.mode = 'OneShot'; % Faster for small dataset fiting in the memory of the computer
cfg.process.bin.mode = 'ByDay'; % Slightly slower but can handle a lot more data at once as it will be binned one day at a time

%% 1. Read Sync & Split data
ila.Read('raw');

%% 2. Bin
ila.Bin()
ila.Write('bin')

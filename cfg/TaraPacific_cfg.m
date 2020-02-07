% TaraPacific Configuration file
% author: Nils
% created: Aug 16, 2018

cfg = struct('meta', struct(), 'instruments', struct(), 'process', struct());

%%%%%%%%%%%%%%
%% METADATA %%
%%%%%%%%%%%%%%
cfg.meta.investigators = 'Emmanuel_Boss,Guillaume_Bourdin';
cfg.meta.affiliations = 'University_of_Maine';
cfg.meta.emails = 'emmanuel.boss@maine.edu';
cfg.meta.experiment = 'Tara';
cfg.meta.cruise = 'TaraPacific';
cfg.meta.station = 'NA';
cfg.meta.documents = 'NA';
cfg.meta.calibration_files = 'NA';
cfg.meta.data_type = 'flow_thru';
cfg.meta.data_status = 'preliminary';
cfg.meta.measurement_depth = 1.5;


%%%%%%%%%%%%%%%%%
%% INSTRUMENTS %%
%%%%%%%%%%%%%%%%%

PATH_ROOT = 'C:\Users\Gui\Documents\MATLAB\InLineAnalysis\';

%%% TSG + GPS %%%
cfg.instruments.TSG = struct();
cfg.instruments.TSG.model = 'TSG';
cfg.instruments.TSG.boat = 'Tara';
cfg.instruments.TSG.path = struct('raw',  [PATH_ROOT 'raw' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'Underway' filesep]);
cfg.instruments.TSG.view = struct('varname', 't');

%%% FTH (FlowControl) %%%
cfg.instruments.FTH = struct();
cfg.instruments.FTH.model = 'FTH';
cfg.instruments.FTH.logger = 'FlowControl';
cfg.instruments.FTH.LoadPrevious = true;
cfg.instruments.FTH.path = struct('raw',  [PATH_ROOT 'raw' filesep 'flow' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'FlowControl' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'FlowControl' filesep]);
cfg.instruments.FTH.view = struct('varname', 'swt');

% cfg.instruments.FTH = struct();
% cfg.instruments.FTH.model = 'FTHv';
% cfg.instruments.FTH.logger = 'FlowControl';
% cfg.instruments.FTH.LoadPrevious = true;
% cfg.instruments.FTH.path = struct('raw',  [PATH_ROOT 'raw' filesep 'flow' filesep],...
%                                   'wk',   [PATH_ROOT 'wk' filesep 'FlowControl' filesep],...
%                                   'prod', [PATH_ROOT 'prod' filesep],...
%                                   'ui', [PATH_ROOT 'ui' filesep 'FlowControl' filesep]);
% cfg.instruments.FTH.view = struct('varname', 'swt');

%%% AC9 %%%
cfg.instruments.AC9 = struct();
cfg.instruments.AC9.model = 'AC9';
cfg.instruments.AC9.sn = '245';
cfg.instruments.AC9.logger = 'WetView';
% cfg.instruments.AC9.device_file = [PATH_ROOT 'acs091_20170410_20170904.dev']; % acs091_20180530_20180818
% cfg.instruments.AC9.lambda_reference = [400.9, 404.6, 408.4, 411.7, 415.5, 419.2, 423.7, 428.1, 432.2, 435.9, 440, 444, 449, 453.4, 457.6, 462, 466.2, 470.8, 475.7, 480.3, 485.2, 489.4, 493.7, 497.7, 502.1, 506.6, 511.4, 515.9, 520.6, 524.9, 529.2, 533.4, 537.5, 541.5, 545.7, 550, 554.3, 558.6, 562.8, 566.9, 570.7, 574.1, 578, 582.2, 586.1, 590.3, 594.5, 599.1, 603.8, 608.1, 612.7, 617.4, 621.7, 626, 630.5, 634.7, 639, 643.3, 647.8, 652.4, 656.8, 661.4, 665.9, 670.3, 674.4, 678.8, 682.8, 687.1, 691, 694.7, 698.4, 702.5, 706.2, 710, 713.7, 717.6, 721.4, 725.1, 728.8, 732.6, 736.2, 739.6];
% cfg.instruments.AC9.lambda_a = [400.9, 404.6, 408.4, 411.7, 415.5, 419.2, 423.7, 428.1, 432.2, 435.9, 440, 444, 449, 453.4, 457.6, 462, 466.2, 470.8, 475.7, 480.3, 485.2, 489.4, 493.7, 497.7, 502.1, 506.6, 511.4, 515.9, 520.6, 524.9, 529.2, 533.4, 537.5, 541.5, 545.7, 550, 554.3, 558.6, 562.8, 566.9, 570.7, 574.1, 578, 582.2, 586.1, 590.3, 594.5, 599.1, 603.8, 608.1, 612.7, 617.4, 621.7, 626, 630.5, 634.7, 639, 643.3, 647.8, 652.4, 656.8, 661.4, 665.9, 670.3, 674.4, 678.8, 682.8, 687.1, 691, 694.7, 698.4, 702.5, 706.2, 710, 713.7, 717.6, 721.4, 725.1, 728.8, 732.6, 736.2, 739.6];
% cfg.instruments.AC9.lambda_c = [400.9, 404.6, 408.4, 411.7, 415.5, 419.2, 423.7, 428.1, 432.2, 435.9, 440, 444, 449, 453.4, 457.6, 462, 466.2, 470.8, 475.7, 480.3, 485.2, 489.4, 493.7, 497.7, 502.1, 506.6, 511.4, 515.9, 520.6, 524.9, 529.2, 533.4, 537.5, 541.5, 545.7, 550, 554.3, 558.6, 562.8, 566.9, 570.7, 574.1, 578, 582.2, 586.1, 590.3, 594.5, 599.1, 603.8, 608.1, 612.7, 617.4, 621.7, 626, 630.5, 634.7, 639, 643.3, 647.8, 652.4, 656.8, 661.4, 665.9, 670.3, 674.4, 678.8, 682.8, 687.1, 691, 694.7, 698.4, 702.5, 706.2, 710, 713.7, 717.6, 721.4, 725.1, 728.8, 732.6, 736.2, 739.6];
cfg.instruments.AC9.path = struct('raw',  [PATH_ROOT 'raw' filesep 'AC9' filesep],...
                                  'di',  [PATH_ROOT 'raw' filesep 'AC9' filesep 'DI' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'AC9' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'AC9' filesep]);
cfg.instruments.AC9.view = struct('varname', 'a', 'varcol', 9);

%%% ACS 091 %%% (Aug 20 to ...)
cfg.instruments.ACS091 = struct();
cfg.instruments.ACS091.model = 'ACS';
cfg.instruments.ACS091.sn = '091';
cfg.instruments.ACS091.logger = 'Compass_2.1rc_scheduled_bin';
cfg.instruments.ACS091.device_file = [PATH_ROOT 'acs091_20180530_20180818.dev']; % acs091_20170410_20170904 acs091_20180530_20180818 
% cfg.instruments.ACS091.logger = 'Compass_2.1rc_scheduled';
% cfg.instruments.ACS091.lambda_reference = [400.9, 404.6, 408.4, 411.7, 415.5, 419.2, 423.7, 428.1, 432.2, 435.9, 440, 444, 449, 453.4, 457.6, 462, 466.2, 470.8, 475.7, 480.3, 485.2, 489.4, 493.7, 497.7, 502.1, 506.6, 511.4, 515.9, 520.6, 524.9, 529.2, 533.4, 537.5, 541.5, 545.7, 550, 554.3, 558.6, 562.8, 566.9, 570.7, 574.1, 578, 582.2, 586.1, 590.3, 594.5, 599.1, 603.8, 608.1, 612.7, 617.4, 621.7, 626, 630.5, 634.7, 639, 643.3, 647.8, 652.4, 656.8, 661.4, 665.9, 670.3, 674.4, 678.8, 682.8, 687.1, 691, 694.7, 698.4, 702.5, 706.2, 710, 713.7, 717.6, 721.4, 725.1, 728.8, 732.6, 736.2, 739.6];
% cfg.instruments.ACS091.lambda_a = [400.9, 404.6, 408.4, 411.7, 415.5, 419.2, 423.7, 428.1, 432.2, 435.9, 440, 444, 449, 453.4, 457.6, 462, 466.2, 470.8, 475.7, 480.3, 485.2, 489.4, 493.7, 497.7, 502.1, 506.6, 511.4, 515.9, 520.6, 524.9, 529.2, 533.4, 537.5, 541.5, 545.7, 550, 554.3, 558.6, 562.8, 566.9, 570.7, 574.1, 578, 582.2, 586.1, 590.3, 594.5, 599.1, 603.8, 608.1, 612.7, 617.4, 621.7, 626, 630.5, 634.7, 639, 643.3, 647.8, 652.4, 656.8, 661.4, 665.9, 670.3, 674.4, 678.8, 682.8, 687.1, 691, 694.7, 698.4, 702.5, 706.2, 710, 713.7, 717.6, 721.4, 725.1, 728.8, 732.6, 736.2, 739.6];
% cfg.instruments.ACS091.lambda_c = [400.9, 404.6, 408.4, 411.7, 415.5, 419.2, 423.7, 428.1, 432.2, 435.9, 440, 444, 449, 453.4, 457.6, 462, 466.2, 470.8, 475.7, 480.3, 485.2, 489.4, 493.7, 497.7, 502.1, 506.6, 511.4, 515.9, 520.6, 524.9, 529.2, 533.4, 537.5, 541.5, 545.7, 550, 554.3, 558.6, 562.8, 566.9, 570.7, 574.1, 578, 582.2, 586.1, 590.3, 594.5, 599.1, 603.8, 608.1, 612.7, 617.4, 621.7, 626, 630.5, 634.7, 639, 643.3, 647.8, 652.4, 656.8, 661.4, 665.9, 670.3, 674.4, 678.8, 682.8, 687.1, 691, 694.7, 698.4, 702.5, 706.2, 710, 713.7, 717.6, 721.4, 725.1, 728.8, 732.6, 736.2, 739.6];
cfg.instruments.ACS091.path = struct('raw',  [PATH_ROOT 'raw' filesep 'ACS' filesep],...
                                  'di',  [PATH_ROOT 'raw' filesep 'ACS' filesep 'DI' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'ACS' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'ACS' filesep]);
cfg.instruments.ACS091.view = struct('varname', 'a', 'varcol', 40);

%%% ACS 111 %%% (Aug 11 to Aug 20)
cfg.instruments.ACS111 = struct();
cfg.instruments.ACS111.model = 'ACS';
cfg.instruments.ACS111.sn = '111';
cfg.instruments.ACS111.ila_prefix = 'ACS';
cfg.instruments.ACS111.logger = 'Compass_2.1rc_scheduled_bin';
cfg.instruments.ACS111.device_file = [PATH_ROOT 'acs111_20171212_20180530.dev'];
% cfg.instruments.ACS111.logger = 'Compass_2.1rc_scheduled';
% cfg.instruments.ACS111.lambda_reference = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
% cfg.instruments.ACS111.lambda_a = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
% cfg.instruments.ACS111.lambda_c = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
cfg.instruments.ACS111.path = struct('raw',  [PATH_ROOT 'raw' filesep 'ACS' filesep],...
                                  'di',  [PATH_ROOT 'raw' filesep 'ACS' filesep 'DI' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'ACS' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'ACS' filesep]);
cfg.instruments.ACS111.view = struct('varname', 'a', 'varcol', 40);

%%% ACS 007 %%% (Aug 11 to Aug 20)
cfg.instruments.ACS007 = struct();
cfg.instruments.ACS007.model = 'ACS'; 
cfg.instruments.ACS007.sn = '007';
cfg.instruments.ACS007.ila_prefix = 'ACS';
cfg.instruments.ACS007.logger = 'WetView'; % 'WetView' 'Compass_2.1rc' 'Compass_2.1rc_scheduled' 'Compass_2.1rc_scheduled_bin'
cfg.instruments.ACS007.device_file = [PATH_ROOT 'acs007_20161101_20170220.dev']; % acs007_20160528_20160704  
% cfg.instruments.ACS007.logger = 'Compass_2.1rc_scheduled';
% cfg.instruments.ACS007.lambda_reference = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
% cfg.instruments.ACS007.lambda_a = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
% cfg.instruments.ACS007.lambda_c = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
cfg.instruments.ACS007.path = struct('raw',  [PATH_ROOT 'raw' filesep 'ACS' filesep],...
                                  'di',  [PATH_ROOT 'raw' filesep 'ACS' filesep 'DI' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'ACS' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'ACS' filesep]);
cfg.instruments.ACS007.view = struct('varname', 'a', 'varcol', 40);

%%% ACS 057 %%% (Aug 11 to Aug 20)
cfg.instruments.ACS057 = struct();
cfg.instruments.ACS057.model = 'ACS';
cfg.instruments.ACS057.sn = '057';
cfg.instruments.ACS057.ila_prefix = 'ACS';
cfg.instruments.ACS057.logger = 'Compass_2.1rc_scheduled_bin';
cfg.instruments.ACS057.device_file = [PATH_ROOT 'acs057_20160704_20160720.dev'];
% cfg.instruments.ACS057.logger = 'Compass_2.1rc_scheduled';
% cfg.instruments.ACS057.lambda_reference = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
% cfg.instruments.ACS057.lambda_a = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
% cfg.instruments.ACS057.lambda_c = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
cfg.instruments.ACS057.path = struct('raw',  [PATH_ROOT 'raw' filesep 'ACS' filesep],...
                                  'di',  [PATH_ROOT 'raw' filesep 'ACS' filesep 'DI' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'ACS' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'ACS' filesep]);
cfg.instruments.ACS057.view = struct('varname', 'a', 'varcol', 40);

%%% ACS 279 %%% (Aug 11 to Aug 20)
cfg.instruments.ACS279 = struct();
cfg.instruments.ACS279.model = 'ACS';
cfg.instruments.ACS279.sn = '279';
cfg.instruments.ACS279.ila_prefix = 'ACS';
cfg.instruments.ACS279.logger = 'Compass_2.1rc_scheduled_bin';
cfg.instruments.ACS279.device_file = [PATH_ROOT 'acs279_20180821_20180920.dev']; % acs279_20170902_20171213 acs279_20180821_20180920
% cfg.instruments.ACS279.logger = 'Compass_2.1rc_scheduled';
% cfg.instruments.ACS279.lambda_reference = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
% cfg.instruments.ACS279.lambda_a = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
% cfg.instruments.ACS279.lambda_c = [401 ,405 ,408.2 ,411.6 ,415.1 ,419.2 ,423.3 ,427.8 ,431.5 ,435.4 ,439.3 ,443.7 ,448.1 ,452.7 ,456.6 ,460.8 ,464.8 ,469.4 ,474 ,478.5 ,483.2 ,487.6 ,491.8 ,496.1 ,500.1 ,504.7 ,509.2 ,514 ,518.4 ,522.9 ,527.2 ,531.3 ,535.3 ,539.4 ,543.7 ,547.8 ,552.1 ,556.3 ,560.6 ,564.7 ,568.6 ,572.5 ,576.1 ,579.7 ,583.2 ,587.5 ,591.5 ,595.7 ,599.9 ,604.3 ,608.7 ,613 ,617.5 ,621.7 ,626 ,630.4 ,634.4 ,638.7 ,643 ,647.3 ,651.7 ,656.2 ,660.6 ,664.9 ,669.1 ,673.6 ,677.6 ,681.9 ,685.9 ,690 ,693.7 ,697.9 ,701.8 ,705.5 ,709.2 ,713.2 ,717.1 ,721.1 ,724.8 ,728.6 ,732.1 ,735.8 ,739.6 ,742.9];
cfg.instruments.ACS279.path = struct('raw',  [PATH_ROOT 'raw' filesep 'ACS' filesep],...
                                  'di',  [PATH_ROOT 'raw' filesep 'ACS' filesep 'DI' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'ACS' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'ACS' filesep]);
cfg.instruments.ACS279.view = struct('varname', 'a', 'varcol', 40);

%%% BB3 %%%
cfg.instruments.BB3 = struct();
cfg.instruments.BB3.model = 'BB';
cfg.instruments.BB3.sn = '1502';
cfg.instruments.BB3.ila_prefix = 'BB3';
cfg.instruments.BB3.logger = 'InlininoBB3';
cfg.instruments.BB3.lambda = [470,532,650];
cfg.instruments.BB3.theta = 124;
cfg.instruments.BB3.slope = [1.066E-05,7.076E-06,3.569E-06];
% cfg.instruments.BB3.slope = [8.407E-06,4.624E-06,4.090E-06];
cfg.instruments.BB3.dark = [50,44,45];
cfg.instruments.BB3.path = struct('raw',  [PATH_ROOT 'raw' filesep 'BB3' filesep],...
                                  'di',  [PATH_ROOT 'raw' filesep 'BB3' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'BB3' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'BB3' filesep]);
cfg.instruments.BB3.view = struct('varname', 'beta', 'varcol', 2);

%%% WSCD %%% (Mai 2016 to Oct 2018)
cfg.instruments.WSCD1082P = struct();
cfg.instruments.WSCD1082P.model = 'CD';
cfg.instruments.WSCD1082P.sn = '1082P';
cfg.instruments.WSCD1082P.ila_prefix = 'WSCD';
cfg.instruments.WSCD1082P.logger = 'InlininoWSCD';
cfg.instruments.WSCD1082P.slope = 62;
cfg.instruments.WSCD1082P.dark = 0.059;
cfg.instruments.WSCD1082P.path = struct('raw',  [PATH_ROOT 'raw' filesep 'WSCD' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'WSCD' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'WSCD' filesep]);
cfg.instruments.WSCD1082P.view = struct('varname', 'fdom');

% %%% WSCD %%% (Aug 11 to 17)
% cfg.instruments.WSCD859 = struct();
% cfg.instruments.WSCD859.model = 'CD';
% cfg.instruments.WSCD859.sn = '859';
% cfg.instruments.WSCD859.ila_prefix = 'WSCD';
% cfg.instruments.WSCD859.logger = 'InlininoWSCD';
% cfg.instruments.WSCD859.slope = 0.0909;
% cfg.instruments.WSCD859.dark = 44;
% cfg.instruments.WSCD859.path = struct('raw',  [PATH_ROOT 'raw/WSCD859/'],...
%                                   'wk',   [PATH_ROOT 'wk/WSCD859/'],...
%                                   'prod', [PATH_ROOT 'prod/'],...
%                                   'ui', [PATH_ROOT 'ui/WSCD859/']);
% cfg.instruments.WSCD859.view = struct('varname', 'fdom');
% 
% %%% LISST %%%
% cfg.instruments.LISST = struct();
% cfg.instruments.LISST.model = 'LISST';
% cfg.instruments.LISST.type = 'B'; 
% cfg.instruments.LISST.sn = '1183';
% cfg.instruments.LISST.logger = 'TeraTerm';
% cfg.instruments.LISST.zsc = [2.203500e+001, 2.568500e+001, 2.503000e+001, 2.986000e+001, 2.842500e+001, 3.283000e+001, 3.077000e+001, 3.659500e+001, 2.978000e+001, 3.552000e+001, 3.198000e+001, 4.216000e+001, 3.916500e+001, 4.662500e+001, 3.974000e+001, 4.454000e+001, 4.403500e+001, 4.604500e+001, 4.430000e+001, 4.510500e+001, 4.719500e+001, 3.850000e+001, 5.373000e+001, 2.664000e+001, 3.180500e+001, 1.655500e+001, 2.205500e+001, 1.554000e+001, 1.422000e+001, 1.123000e+001, 8.780000e+000, 8.555000e+000, 1.515000e+003, 1.167900e+003, 6.410000e+001, 1.055150e+003, 7.700000e+001, 2.116600e+003, 1.807000e+003, 5.476500e+003];
% % Original dcal file (ring area)
% % cfg.instruments.LISST.dcal = [1.0000000e+000, 1.0038000e+000, 9.9360000e-001, 1.0027000e+000, 9.9720000e-001, 9.9570000e-001, 9.9030000e-001, 9.9430000e-001, 9.9290000e-001, 9.9000000e-001, 9.9290000e-001, 9.9300000e-001, 9.9150000e-001, 9.9300000e-001, 9.9230000e-001, 9.9090000e-001, 1.1032000e+000, 1.1123000e+000, 1.2430000e+000, 1.1562000e+000, 1.3273000e+000, 1.1999000e+000, 1.0740000e+000, 1.7489000e+000, 1.5382000e+000, 2.5109000e+000, 2.5468000e+000, 3.5504000e+000, 3.9338000e+000, 5.1747342e+000, 7.5143548e+000, 1.2528083e+001];
% % dcal adhoc from email of Wayne Slade Dec 8, 2017
% cfg.instruments.LISST.dcal = [ 1.0179083e+00	   9.9213489e-01	   1.0108161e+00	   9.9492883e-01	   1.0043707e+00	   9.9891840e-01	   9.9859055e-01	   1.0042049e+00	   1.0000763e+00	   9.9889997e-01	   1.0009497e+00	   1.0004019e+00	   1.0011130e+00	   1.0004677e+00	   1.0213554e+00	   9.9990262e-01	   1.1115630e+00	   1.1206668e+00	   1.2493699e+00	   1.1643199e+00	   1.3355657e+00	   1.2090892e+00	   1.0781540e+00	   1.7620752e+00	   1.5508563e+00	   2.5304119e+00	   2.5638592e+00	   3.5757212e+00	   3.9631987e+00	   5.0166411e+00	   5.6381118e+00	   8.6881539e+00];
% % customize dcal of 1183 on Jan 16, 2018 due to bump in VSF at ring 30
% cfg.instruments.LISST.dcal(30) = 2.6; % Good for low values but bad for high values ??
% cfg.instruments.LISST.vcc = 48493;
% cfg.instruments.LISST.inversion = 'spherical';
% cfg.instruments.LISST.ds = [1.2500,1.4750,1.7405,2.0538,2.4235,2.8597,3.3744,3.9818,4.6986,5.5443,6.5423,7.7199,9.1095,10.7492,12.6841,14.9672,17.6613,20.8403,24.5916,29.0180,34.2413,40.4047,47.6776,56.2595,66.3863,78.3358,92.4362,109.0747,128.7082,151.8757,179.2133,211.4717,249.5366];
% cfg.instruments.LISST.theta = [0.082, 0.096, 0.114, 0.134, 0.158, 0.187, 0.221, 0.260, 0.307, 0.362, 0.428, 0.505, 0.596, 0.703, 0.829, 0.979, 1.155, 1.363, 1.609, 1.898, 2.240, 2.643, 3.119, 3.681, 4.344, 5.126, 6.049, 7.138, 8.424, 9.941, 11.73, 13.84];
% cfg.instruments.LISST.path = struct('raw',  [PATH_ROOT 'raw/LISST1183/'],...
%                                   'wk',   [PATH_ROOT 'wk/LISST1183/'],...
%                                   'prod', [PATH_ROOT 'prod/'],...
%                                   'ui', [PATH_ROOT 'ui/LISST1183/']);
% cfg.instruments.LISST.view = struct('varname', 'beta', 'varcol', 15);
% 
% %%% ALFA %%%
% cfg.instruments.ALFA = struct();
% cfg.instruments.ALFA.model = 'ALFA';
% cfg.instruments.ALFA.sn = '011';
% cfg.instruments.ALFA.logger = 'ALFA_LabView_m';
% cfg.instruments.ALFA.path = struct('raw',  [PATH_ROOT 'raw/ALFA011/'],...
%                                   'wk',   [PATH_ROOT 'wk/ALFA011/'],...
%                                   'prod', [PATH_ROOT 'prod/'],...
%                                   'ui', [PATH_ROOT 'ui/ALFA011/']);
% cfg.instruments.ALFA.view = struct('varname', 'FvFm');

%%% PAR %%% (Mai 2016 to Oct 2018)
cfg.instruments.PAR = struct();
cfg.instruments.PAR.model = 'PAR';
cfg.instruments.PAR.sn = '50168';
cfg.instruments.PAR.logger = 'Inlinino';
cfg.instruments.PAR.scale = 6.451E-04; % Volts/(uE/m²sec)
cfg.instruments.PAR.dark = 9.7E-03;
cfg.instruments.PAR.path = struct('raw',  [PATH_ROOT 'raw' filesep 'PAR' filesep],...
                                  'wk',   [PATH_ROOT 'wk' filesep 'PAR' filesep],...
                                  'prod', [PATH_ROOT 'prod' filesep],...
                                  'ui', [PATH_ROOT 'ui' filesep 'PAR' filesep]);
cfg.instruments.PAR.view = struct('varname', 'par');
%%%%%%%%%%%%%
%% PROCESS %%
%%%%%%%%%%%%%

%%% General parameters %%%
cfg.process.days2run = datenum(2018,5,30,0,0,0):datenum(2018,6,20,0,0,0);
cfg.process.instruments2run = {'FTH', 'TSG', 'ACS007','BB3','PAR'};
% cfg.process.instruments2run = {'FTH', 'PAR'};
cfg.process.write = true;
cfg.process.force_import = false;
cfg.process.parallel = Inf; % 0: disable parallel or Inf: as many thread available
cfg.process.di = struct();
cfg.process.di.skip = {'TSG'};
cfg.process.di.qc = struct('mode', 'ui');
cfg.process.di.bin = struct('bin_size', 30);

%%% Synchronization %%%
cfg.process.sync = struct();
cfg.process.sync.delay = struct();
cfg.process.sync.delay.FTH = 30;
cfg.process.sync.delay.AC9 = 60; % 
cfg.process.sync.delay.ACS007 = 60; % 
cfg.process.sync.delay.ACS057 = 60; % 
cfg.process.sync.delay.ACS091 = 60; % 280 95 65 82 60 57 52 -10 52 95
cfg.process.sync.delay.ACS111 = 60; % 280 95 65 82 60 57
cfg.process.sync.delay.ACS279 = 60; %70 128 75 62 58 10 60
cfg.process.sync.delay.BB3 = 0; % 55 45 22 20
% cfg.process.sync.delay.LISST = 1;
% cfg.process.sync.delay.WSCD859 = 5;
cfg.process.sync.delay.WSCD1082P = 40;
% cfg.process.sync.delay.ALFA = 15;
cfg.process.sync.skip = {'TSG','PAR'};

%%% QC Reference (Flow Control/FTH) %%%
cfg.process.qcref = struct();
cfg.process.qcref.reference = 'FTH';
cfg.process.qcref.view = 'WSCD1082P';
cfg.process.qcref.mode = 'ui'; % load or ui
cfg.process.qcref.MinFiltPeriod = 50; % filter even period in minute

%%% Split total/filtered %%%
cfg.process.split = struct();
cfg.process.split.reference = 'FTH';
cfg.process.split.buffer = struct();
cfg.process.split.buffer.AC9 = [180, 60];
cfg.process.split.buffer.ACS279 = [180, 60];
cfg.process.split.buffer.ACS091 = [180, 60];
cfg.process.split.buffer.ACS111 = [180, 60];
cfg.process.split.buffer.ACS007 = [180, 60];
cfg.process.split.buffer.ACS057 = [180, 60];
cfg.process.split.buffer.BB3 = [400, 220];
% cfg.process.split.buffer.LISST = [540, 360];
% cfg.process.split.buffer.WSCD859 = [310, 20];
cfg.process.split.buffer.WSCD1082P = [310, 20];
% cfg.process.split.buffer.ALFA = [180, 30];
% cfg.process.split.skip = {'FTH', 'TSG','WSCD1082P','PAR','ACS279'};
cfg.process.split.skip = {'FTH', 'TSG','WSCD1082P','PAR'};

%%% Binning %%%
cfg.process.bin = struct('bin_size', struct());
cfg.process.bin.prctile_detection = [2.5, 97.5];
% if StepQC with ACS: prctile_average = [2.5, 97.5]; otherwise prctile_average = [5, 75];
cfg.process.bin.prctile_average = [5, 75];
% Bin mode does not affect the outcome of the data but just the way the data is presented to the computer
% cfg.process.bin.mode = 'OneShot'; % Faster for small dataset fiting in the memory of the computer
cfg.process.bin.mode = 'ByDay'; % Slightly slower but can handle a lot more data at once as it will be binned one day at a time
cfg.process.bin.bin_size.FTH = 1;
cfg.process.bin.bin_size.AC9 = 1;
cfg.process.bin.bin_size.ACS007 = 1;
cfg.process.bin.bin_size.ACS057 = 1;
cfg.process.bin.bin_size.ACS091 = 1;
cfg.process.bin.bin_size.ACS279 = 1;
cfg.process.bin.bin_size.ACS111 = 1;
cfg.process.bin.bin_size.BB3 = 1;
cfg.process.bin.bin_size.PAR = 1;
% cfg.process.bin.bin_size.LISST = 10;
% cfg.process.bin.bin_size.WSCD859 = 1;
cfg.process.bin.bin_size.WSCD1082P = 1;
cfg.process.bin.bin_size.TSG = 1; % TSG does not need to be binned in most cases
% cfg.process.bin.bin_size.ALFA = 10;
cfg.process.bin.skip = {'TSG'};

%%% Automatically flagging %%%
cfg.process.flag = struct();
cfg.process.flag.skip = {'FTH', 'TSG', 'AC9', 'ACS007', 'ACS057','ACS091', 'ACS111', 'ACS279', 'BB3', 'LISST', 'WSCD859', 'WSCD1082P', 'ALFA','PAR'};
% Default: parameters set to all instruments if not specific parameters set
cfg.process.flag.default = struct();
cfg.process.flag.default.maximum_fudge_factor = 4;
cfg.process.flag.default.variance_fudge_factor = 3;
cfg.process.flag.default.avg_sensitivity = 1;
cfg.process.flag.default.unc1_sensitivity = 1;
cfg.process.flag.default.unc2_sensitivity = 2;
cfg.process.flag.default.smooth_threshold = 60;
cfg.process.flag.default.min_flag_n = 1;
cfg.process.flag.default.filt = struct('smooth_threshold', 2);

  
%%% Manually QC %%%
cfg.process.qc = struct();
cfg.process.qc.StepQCLim.a = 3;
cfg.process.qc.StepQCLim.c = 3;
cfg.process.qc.StepQCLim.bb = 3;
cfg.process.qc.Saturation_Threshold_bb = 4000; % (counts)
cfg.process.qc.mode = 'ui';
cfg.process.qc.global = struct();
cfg.process.qc.global.active = false;
cfg.process.qc.global.view = {'ACS279','BB3','PAR'};
cfg.process.qc.global.apply = {'ACS279','BB3','PAR'};
cfg.process.qc.specific = struct();
cfg.process.qc.specific.active = true;
cfg.process.qc.specific.run = {'ACS279','BB3','PAR'};
  
%%% Calibrate %%%
cfg.process.calibrate = struct();
cfg.process.calibrate.ACS007 = struct('compute_dissolved', false,...
                                  'interpolation_method', 'linear', ...
                                  'CDOM_source', 'WSCD1082P',... % Pay attention to serial number (could merge WSCD products first)
                                  'FTH_source', 'FTH');
cfg.process.calibrate.BB3 = struct('compute_dissolved', false,...
                                   'TSG_source', 'TSG',...
                                   'di_method', 'constant');
% cfg.process.calibrate.skip = {'FTH', 'TSG'};
cfg.process.calibrate.skip = {'FTH', 'TSG'};

%%% Write %%%
cfg.process.write = struct();
cfg.process.write.mode = 'One day one file'; % 'One file' or 'One day one file'
cfg.process.write.skip = {'TSG'};

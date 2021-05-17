function visProd_timeseries(data, instrument, lambda)
% Plot InLineAnalysis product time series to assess quality of processing
% and include in processing report
%
%% Author: Guillaume Bourdin
% Date: 03 Mar. 2021
%
% INPUT:
%   data: <NxM table> data containing:
%     - <1xM datenum> time vector
%     - <N-1xM double> data
%   instrument: <char> instrument name
%%
instrument = instrument(isstrprop(instrument,'alpha'));

% list_instru = {'ACS', 'AC', 'BB', 'TSG', 'PAR', 'WSCD', 'HBB'};
% idx = false(size(list_instru));
% for i = 1:max(size(list_instru)); idx(i) = contains(instrument, list_instru{i}); end

switch instrument
  case {'ACS', 'AC'}
    if all(contains({'Halh_chl', 'HH_G50'}, data.Properties.VariableNames))
      figure(77)
      clf
      subplot(3,1,1)
      hold on
      yyaxis('left')
      scatter(datetime(data.dt, 'ConvertFrom', 'datenum'), data.gamma, 6, 'filled')
      ylabel('gamma (unitless)')
      yyaxis('right')
      scatter(datetime(data.dt, 'ConvertFrom', 'datenum'), data.HH_G50, 6, 'filled')
      ylabel('H&H phytoplankton G50: cross-sectional area (\mum)')
      legend('gamma', 'H&H phytoplankton G50')
      hold off
      subplot(3,1,2)
      hold on
      scatter(datetime(data.dt, 'ConvertFrom', 'datenum'), data.Halh_chl, 6, 'filled')
      scatter(datetime(data.dt, 'ConvertFrom', 'datenum'), data.chl, 6, 'filled')
      ylabel('[chl] (mg.m^{-3})')
      legend('Houskeeper [chl]', 'a_{p676}[chl]')
      hold off
      subplot(3,1,3)
      scatter(datetime(data.dt, 'ConvertFrom', 'datenum'), data.poc, 6, 'filled')
      ylabel('[poc] (mg.m^{-3})')
      xlim([min(datetime(data.dt, 'ConvertFrom', 'datenum')) ...
        max(datetime(data.dt, 'ConvertFrom', 'datenum'))]);
      %
      figure(78);
      clf
      subplot(1,2,1)
      scatter(data.gamma, data.HH_G50, 6, 'filled')
      set(gca, 'XScale', 'log', 'YScale', 'log')
      xlabel('gamma (unitless)')
      ylabel('H&H phytoplankton G50: cross-sectional area (\mum)')
      subplot(1,2,2)
      scatter(data.chl, data.Halh_chl, 6, 'filled')
      set(gca, 'XScale', 'log', 'YScale', 'log')
      xlabel('a_{p676}[chl] (mg.m^{-3})')
      ylabel('Houskeeper [chl] (mg.m^{-3})')
    elseif contains({'gamma'}, data.Properties.VariableNames)
      figure(77)
      clf
      subplot(2,1,1)
      hold on
      yyaxis('left')
      scatter(datetime(data.dt, 'ConvertFrom', 'datenum'), data.gamma, 6, 'filled')
      ylabel('gamma (unitless)')
      yyaxis('right')
      scatter(datetime(data.dt, 'ConvertFrom', 'datenum'), data.chl, 6, 'filled')
      ylabel('[chl] (mg.m^{-3})')
      legend('gamma', 'a_{p676}[chl]')
      hold off
      subplot(2,1,2)
      scatter(datetime(data.dt, 'ConvertFrom', 'datenum'), data.poc, 6, 'filled')
      ylabel('[poc] (mg.m^{-3})')
      xlim([min(datetime(data.dt, 'ConvertFrom', 'datenum')) ...
        max(datetime(data.dt, 'ConvertFrom', 'datenum'))]);
    end
  case {'BB', 'HBB'}
    if size(lambda, 1) > size(lambda, 2); lambda = lambda'; end
    if any(contains(data.Properties.VariableNames, 'betap'))
      toplot = 'poc';
      unit = '(mg.m^{-3})';
      fignum = 79;
    elseif any(contains(data.Properties.VariableNames, 'betag'))
      toplot = 'betag';
      unit = 'm^{-1}';
      fignum = 80;
    end
    if contains(instrument, 'HBB')
      data.(toplot)(:, lambda ~= 430 & lambda ~= 550 & lambda ~= 660 & lambda ~= 680) = [];
      lambda(lambda ~= 430 & lambda ~= 550 & lambda ~= 660 & lambda ~= 680) = [];
    end
    C = reshape(spectrumRGB(lambda), max(size(lambda)),  3);
    figure(fignum);
    clf
    yyaxis('left')
    hold on
    scatter(datetime(data.dt, 'ConvertFrom', 'datenum'), data.(toplot), 10, C, 'filled');
    ylabel([toplot ' ' unit]);
    leg = cellfun(@(c) [toplot '_{' c 'nm}'], cellstr(num2str(lambda')), 'un', 0);
    if contains(instrument, 'HBB') && any(contains(data.Properties.VariableNames, 'gamma_bbp'))
      yyaxis('right')
      hold on
      scatter(datetime(data.dt, 'ConvertFrom', 'datenum'), data.gamma_bbp, 50, ...
        'k', 'filled', 'Marker', 'v');
      ylabel('Gamma bbp (unitless)');
      leg = [leg; {'gamma bbp'}];
    elseif contains(instrument, 'HBB') && any(contains(data.Properties.VariableNames, 'gamma_bbg'))
      yyaxis('right')
      scatter(datetime(data.dt, 'ConvertFrom', 'datenum'), data.gamma_bbg, 50, ...
        'k', 'filled', 'Marker', 'v');
      ylabel('Gamma bbg (unitless)');
      leg = [leg; {'gamma bbg'}];
    end
    hold off
    legend(leg)
    xlim([min(datetime(data.dt, 'ConvertFrom', 'datenum')) ...
      max(datetime(data.dt, 'ConvertFrom', 'datenum'))]);
  case 'TSG'
    figure(81);
    clf
    yyaxis('left')
    scatter(datetime(data.dt, 'ConvertFrom', 'datenum'), data.t, 6, 'filled'); ylabel('TSG T (°C)');
    yyaxis('right')
    scatter(datetime(data.dt, 'ConvertFrom', 'datenum'), data.s, 6, 'filled'); ylabel('TSG S (PSU)');
    xlim([min(datetime(data.dt, 'ConvertFrom', 'datenum')) ...
      max(datetime(data.dt, 'ConvertFrom', 'datenum'))]);
  case 'PAR'
    figure(82);
    clf
    scatter(datetime(data.dt, 'ConvertFrom', 'datenum'), data.par, 6, 'filled');
    ylabel('PAR (\muE.m^{-2}.s^{-1})');
    xlim([min(datetime(data.dt, 'ConvertFrom', 'datenum')) ...
      max(datetime(data.dt, 'ConvertFrom', 'datenum'))]);
  case 'WSCD'
    figure(83);
    clf
    scatter(datetime(data.dt, 'ConvertFrom', 'datenum'), data.fdom, 6, 'filled');
    ylabel('FDOM ppb');
    xlim([min(datetime(data.dt, 'ConvertFrom', 'datenum')) ...
      max(datetime(data.dt, 'ConvertFrom', 'datenum'))]);
%   case 'LISST'
    
  otherwise
    warning('%s not supported for product visualisation', instrument)
end
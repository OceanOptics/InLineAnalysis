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
if ~isdatetime(data.dt)
  data.dt = datetime(data.dt, 'ConvertFrom', 'datenum');
end

% list_instru = {'ACS', 'AC', 'BB', 'TSG', 'PAR', 'WSCD', 'HBB'};
% idx = false(size(list_instru));
% for i = 1:max(size(list_instru)); idx(i) = contains(instrument, list_instru{i}); end

switch instrument
  case {'ACS', 'AC'}
    if all(contains({'Halh_chl', 'HH_G50'}, data.Properties.VariableNames))
      vargam = data.Properties.VariableNames{contains(data.Properties.VariableNames, ...
        'gamma')};
      varHH_G50 = data.Properties.VariableNames{contains(data.Properties.VariableNames, ...
        'HH_G50')};
      varchl = data.Properties.VariableNames{contains(data.Properties.VariableNames, ...
        {'ap676lh_chl', 'Chl_lineheight'})};
      varHchl = data.Properties.VariableNames{contains(data.Properties.VariableNames, ...
        'Halh_chl')};
      varPOC = data.Properties.VariableNames{contains(data.Properties.VariableNames, ...
        {'poc','POC','POC_cp'})};
      fig(77);
      clf
      subplot(3,1,1)
      hold on
      yyaxis('left')
      scatter(data.dt, data.(vargam), 6, 'filled')
      ylabel('gamma (unitless)')
      yyaxis('right')
      scatter(data.dt, data.(varHH_G50), 6, 'filled')
      ylabel('H&H phytoplankton G50: cross-sectional area (\mum)')
      legend('gamma', 'H&H phytoplankton G50')
      hold off
      subplot(3,1,2)
      hold on
      scatter(data.dt, data.(varHchl), 6, 'filled')
      scatter(data.dt, data.(varchl), 6, 'filled')
      ylabel('[chl] (mg.m^{-3})')
      legend('Houskeeper [chl]', 'a_{p676}[chl]')
      hold off
      subplot(3,1,3)
      scatter(data.dt, data.(varPOC), 6, 'filled')
      ylabel('[poc] (mg.m^{-3})')
      xlim([min(data.dt) ...
        max(data.dt)]);
      %
      fig(78);
      clf
      subplot(1,2,1)
      scatter(data.(vargam), data.(varHH_G50), 6, 'filled')
      set(gca, 'XScale', 'log', 'YScale', 'log')
      xlabel('gamma (unitless)')
      ylabel('H&H phytoplankton G50: cross-sectional area (\mum)')
      subplot(1,2,2)
      scatter(data.(varchl), data.(varHchl), 6, 'filled')
      set(gca, 'XScale', 'log', 'YScale', 'log')
      xlabel('a_{p676}[chl] (mg.m^{-3})')
      ylabel('Houskeeper [chl] (mg.m^{-3})')
    elseif any(contains(data.Properties.VariableNames, 'gamma'))
      vargam = data.Properties.VariableNames{contains(data.Properties.VariableNames, ...
        'gamma')};
      varchl = data.Properties.VariableNames{contains(data.Properties.VariableNames, ...
        {'ap676lh_chl', 'Chl_lineheight'})};
      varPOC = data.Properties.VariableNames{contains(data.Properties.VariableNames, ...
        {'poc','POC','POC_cp'})};
      fig(77);
      clf
      subplot(2,1,1)
      hold on
      yyaxis('left')
      scatter(data.dt, data.(vargam), 6, 'filled')
      ylabel('gamma (unitless)')
      yyaxis('right')
      scatter(data.dt, data.(varchl), 6, 'filled')
      ylabel('[chl] (mg.m^{-3})')
      legend('gamma', 'a_{p676}[chl]')
      hold off
      subplot(2,1,2)
      scatter(data.dt, data.(varPOC), 6, 'filled')
      ylabel('[poc] (mg.m^{-3})')
      xlim([min(data.dt) ...
        max(data.dt)]);
    end
    if any(contains(data.Properties.VariableNames, 'base_fit_ag'))
      fig(79);
      clf
      hold on
      yyaxis('left')
      p = [];
      p1 = scatter(data.dt(~data.ag_fitflag), ...
        data.base_fit_ag(~data.ag_fitflag), 20, 'filled');
      p = [p, p1(1)];
      leg = {'a_g'};
      ylabel('a_g fit slope')
      if any(data.ag_fitflag)
        y_lim = ylim;
        p2 = plot([data.dt(data.ag_fitflag) ...
          data.dt(data.ag_fitflag)], ...
          [y_lim(1) y_lim(2)], '--b', 'LineWidth', 0.001);
        p = [p, p2(1)];
        leg = [leg, 'ag fit flagged'];
      end
      yyaxis('right')
      p3 = scatter(data.dt(~data.cg_fitflag), ...
        data.base_fit_cg(~data.cg_fitflag), 20, 'filled');
      p = [p, p3(1)];
      leg = [leg, 'c_g'];
      ylabel('c_g fit slope')
      if any(data.cg_fitflag)
        y_lim = ylim;
        p4 = plot([data.dt(data.cg_fitflag) ...
          data.dt(data.cg_fitflag)], ...
          [y_lim(1) y_lim(2)], '--r', 'LineWidth', 0.001);
        p = [p, p4(1)];
        leg = [leg, 'cg fit flagged'];
      end
      legend(p, leg)
      hold off
      xlim([min(data.dt) ...
        max(data.dt)]);
    end
      
  case {'BB', 'HBB'}
    if size(lambda, 1) > size(lambda, 2); lambda = lambda'; end
    if any(contains(data.Properties.VariableNames, 'betap'))
      toplot = 'poc';
      unit = '(mg.m^{-3})';
      fignum = 85;
    elseif any(contains(data.Properties.VariableNames, 'betag'))
      toplot = 'betag';
      unit = 'm^{-1}';
      fignum = 86;
    end
    if contains(instrument, 'HBB')
      data.(toplot)(:, lambda ~= 430 & lambda ~= 550 & lambda ~= 660 & lambda ~= 680) = [];
      lambda(lambda ~= 430 & lambda ~= 550 & lambda ~= 660 & lambda ~= 680) = [];
    end
    C = reshape(spectrumRGB(lambda), max(size(lambda)),  3);
    fig(fignum);
    clf
    yyaxis('left')
    hold on
    scatter(data.dt, data.(toplot), 10, C, 'filled');
    ylabel([toplot ' ' unit]);
    leg = cellfun(@(c) [toplot '_{' c 'nm}'], cellstr(num2str(lambda')), 'un', 0);
    if contains(instrument, 'HBB') && any(contains(data.Properties.VariableNames, 'gamma_bbp'))
      yyaxis('right')
      hold on
      scatter(data.dt, data.gamma_bbp, 50, ...
        'k', 'filled', 'Marker', 'v');
      ylabel('Gamma bbp (unitless)');
      leg = [leg; {'gamma bbp'}];
    elseif contains(instrument, 'HBB') && any(contains(data.Properties.VariableNames, 'gamma_bbg'))
      yyaxis('right')
      scatter(data.dt, data.gamma_bbg, 50, ...
        'k', 'filled', 'Marker', 'v');
      ylabel('Gamma bbg (unitless)');
      leg = [leg; {'gamma bbg'}];
    end
    hold off
    legend(leg)
    xlim([min(data.dt) ...
      max(data.dt)]);
  case 'TSG'
    fig(90);
    clf
    yyaxis('left')
    scatter(data.dt, data.t, 6, 'filled'); ylabel('TSG T (°C)');
    yyaxis('right')
    scatter(data.dt, data.s, 6, 'filled'); ylabel('TSG S (PSU)');
    xlim([min(data.dt) ...
      max(data.dt)]);
  case 'PAR'
    fig(92);
    clf
    scatter(data.dt, data.par, 6, 'filled');
    ylabel('PAR (\muE.m^{-2}.s^{-1})');
    xlim([min(data.dt) ...
      max(data.dt)]);
  case 'WSCD'
    fig(94);
    clf
    scatter(data.dt, data.fdom, 6, 'filled');
    ylabel('FDOM ppb');
    xlim([min(data.dt) ...
      max(data.dt)]);
%   case 'LISST'
  otherwise
    warning('%s not supported for product visualisation', instrument)
end
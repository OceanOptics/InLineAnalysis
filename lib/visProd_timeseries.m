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
% get rid of numbers in instrument name
instrument = instrument(isstrprop(instrument,'alpha'));
if ~isdatetime(data.dt)
  data.dt = datetime(data.dt, 'ConvertFrom', 'datenum');
end

% list_instru = {'ACS', 'AC', 'BB', 'TSG', 'PAR', 'WSCD', 'HBB','WS3S'};
% idx = false(size(list_instru));
% for i = 1:max(size(list_instru)); idx(i) = contains(instrument, list_instru{i}); end

switch instrument
  case {'ACS', 'AC'}
    if all(contains({'chl_Halh', 'HH_G50'}, data.Properties.VariableNames))
      vargam = data.Properties.VariableNames{contains(data.Properties.VariableNames, ...
        'gamma')};
      varHH_G50 = data.Properties.VariableNames{contains(data.Properties.VariableNames, ...
        'HH_G50')};
      varchl = data.Properties.VariableNames{contains(data.Properties.VariableNames, ...
        {'chl_ap676lh', 'Chl_lineheight'})};
      varHchl = data.Properties.VariableNames{contains(data.Properties.VariableNames, ...
        'chl_Halh')};
      varPOC = data.Properties.VariableNames{contains(data.Properties.VariableNames, ...
        {'poc','POC','POC_cp'})};
      wl550 = abs(lambda - 550) == min(abs(lambda - 550));
      fig(10);
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
      xlim([min(data.dt) max(data.dt)]);
      hold off
      subplot(3,1,2)
      hold on
      scatter(data.dt, data.(varHchl), 6, 'filled')
      scatter(data.dt, data.(varchl), 6, 'filled')
      ylabel('[chl] (mg.m^{-3})')
      legend('Houskeeper [chl]', 'a_{p676}[chl]')
      xlim([min(data.dt) max(data.dt)]);
      hold off
      subplot(3,1,3)
      scatter(data.dt, data.(varPOC), 6, 'filled')
      ylabel('[poc] (mg.m^{-3})')
      xlim([min(data.dt) max(data.dt)]);
      %
      if any(strcmp(data.Properties.VariableNames, 'cp'))
        nb_subplot = 3;
      else
        nb_subplot = 2;
      end
      fig(11);
      clf
      subplot(1,nb_subplot,1)
      binscatter(data.(varchl), data.(varHchl), 250);
      colormap('parula')
      set(gca, 'ColorScale', 'log')
      % scatter(data.(varchl), data.(varHchl), 6, 'filled')
      set(gca, 'XScale', 'log', 'YScale', 'log')
      xlabel('a_{p676}[chl] (mg.m^{-3})')
      ylabel('Houskeeper [chl] (mg.m^{-3})')
      subplot(1,nb_subplot,2)
      binscatter(data.(vargam), data.(varHH_G50), 250);
      colormap('parula')
      set(gca, 'ColorScale', 'log')
      % scatter(data.(vargam), data.(varHH_G50), 6, 'filled')
      set(gca, 'XScale', 'log', 'YScale', 'log')
      xlabel('gamma (unitless)')
      ylabel('H&H phytoplankton G50: cross-sectional area (\mum)')
      if nb_subplot == 3
        subplot(1,nb_subplot,3)
        binscatter(data.cp(:, wl550), data.(vargam), 250);
        colormap('parula')
        set(gca, 'ColorScale', 'log')
        % scatter(data.cp(:, wl550), data.(vargam), 6, 'filled')
        set(gca, 'XScale', 'log', 'YScale', 'log')
        xlabel('c_{p} 550 nm')
        ylabel('gamma c_p (unitless)')
      end
    elseif any(contains(data.Properties.VariableNames, 'gamma'))
      vargam = data.Properties.VariableNames{contains(data.Properties.VariableNames, ...
        'gamma')};
      varchl = data.Properties.VariableNames{contains(data.Properties.VariableNames, ...
        {'chl_ap676lh', 'Chl_lineheight'})};
      varPOC = data.Properties.VariableNames{contains(data.Properties.VariableNames, ...
        {'poc','POC','POC_cp'})};
      fig(12);
      clf
      subplot(2,1,1)
      hold on
      yyaxis('left')
      scatter(data.dt, data.(vargam), 7, 'filled')
      ylabel('gamma (unitless)')
      xlim([min(data.dt) max(data.dt)]);
      yyaxis('right')
      scatter(data.dt, data.(varchl), 7, 'filled')
      ylabel('[chl] (mg.m^{-3})')
      legend('gamma', 'a_{p676}[chl]')
      hold off
      subplot(2,1,2)
      scatter(data.dt, data.(varPOC), 7, 'filled')
      ylabel('[poc] (mg.m^{-3})')
      xlim([min(data.dt) max(data.dt)]);
    end
    if any(contains(data.Properties.VariableNames, 'base_fit_ag'))
      fig(15);
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
      xlim([min(data.dt) max(data.dt)]);
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
      xlim([min(data.dt) max(data.dt)]);
      legend(p, leg)
      hold off
    end
  case {'BB', 'HBB', 'HyperBB'}
    if size(lambda, 1) > size(lambda, 2); lambda = lambda'; end
    lambda_orig = lambda;
    if any(contains(data.Properties.VariableNames, 'betap'))
      toplot = {'poc', 'bbp'};
      unit = '(mg.m^{-3})';
      fignum = [20 21];
    elseif any(contains(data.Properties.VariableNames, 'betag'))
      toplot = {'betag'};
      unit = 'm^{-1}';
      fignum = 25;
    end
    for i = 1:size(toplot, 2)
      if contains(instrument, 'HBB')
        data.(toplot{i})(:, lambda_orig ~= 430 & lambda_orig ~= 550 & lambda_orig ~= 660 & lambda_orig ~= 680) = [];
        lambda(lambda ~= 430 & lambda ~= 550 & lambda ~= 660 & lambda ~= 680) = [];
      end
      C = reshape(spectrumRGB(lambda), max(size(lambda)),  3);
      fig(fignum(i));
      clf
      yyaxis('left')
      hold on
      scatter(data.dt, data.(toplot{i}), 10, C, 'filled');
      ylabel([toplot{i} ' ' unit]);
      leg = cellfun(@(c) [toplot{i} '_{' c 'nm}'], cellstr(num2str(lambda')), 'un', 0);
      if any(contains(data.Properties.VariableNames, 'gamma_bbp'))
        yyaxis('right')
        hold on
        scatter(data.dt, data.gamma_bbp, 50, 'k', 'filled', 'Marker', 'v');
        ylabel('Gamma bbp (unitless)');
        leg = [leg; {'gamma bbp'}];
        xlim([min(data.dt) max(data.dt)]);
      elseif any(contains(data.Properties.VariableNames, 'gamma_bbg'))
        yyaxis('right')
        scatter(data.dt, data.gamma_bbg, 50, 'k', 'filled', 'Marker', 'v');
        ylabel('Gamma bbg (unitless)');
        leg = [leg; {'gamma bbg'}];
        xlim([min(data.dt) max(data.dt)]);
      end
      hold off
      legend(leg)
      xlim([min(data.dt) max(data.dt)]);
    end
    if any(contains(data.Properties.VariableNames, 'gamma_bbp'))
      wl550 = abs(lambda - 550) == min(abs(lambda - 550));
      fig(11);
      clf
      binscatter(data.bbp(:, wl550), data.gamma_bbp, 250);
      colormap('parula')
      set(gca, 'ColorScale', 'log')
      % scatter(data.bbp(:, wl550), data.gamma_bbp, 7, 'filled')
      set(gca, 'XScale', 'log', 'YScale', 'log')
      xlabel('b_{bp} 550 nm')
      ylabel('gamma b_{bp} (unitless)')
    end
  case {'TSG','SBE','atlasTSG'}
    varT = data.Properties.VariableNames{strcmp(data.Properties.VariableNames, 't') | ...
      strcmp(data.Properties.VariableNames, 't2') | ...
      strcmp(data.Properties.VariableNames, 'sst') | ...
      strcmp(data.Properties.VariableNames, 'SST') | ...
      strcmp(data.Properties.VariableNames, 'SST') | ...
      strcmp(data.Properties.VariableNames, 'SST')};
    varS = data.Properties.VariableNames{strcmp(data.Properties.VariableNames, 's') | ...
      strcmp(data.Properties.VariableNames, 'sss') | ...
      strcmp(data.Properties.VariableNames, 'SSS') | ...
      strcmp(data.Properties.VariableNames, 'sss_adj')};
    fig(30);
    clf
    yyaxis('left')
    scatter(data.dt, data.(varT), 6, 'filled'); ylabel('TSG T (°C)');
    yyaxis('right')
    scatter(data.dt, data.(varS), 6, 'filled'); ylabel('TSG S (PSU)');
    xlim([min(data.dt) max(data.dt)]);
  case 'PAR'
    fig(40);
    clf
    scatter(data.dt, data.par, 6, 'filled');
    ylabel('PAR (\muE.m^{-2}.s^{-1})');
    xlim([min(data.dt) max(data.dt)]);
  case {'WSCD','WSCDP','SUVF'}
    fig(50);
    clf
    scatter(data.dt, data.fdom, 6, 'filled');
    ylabel('FDOM [v uncalibrated]');
    xlim([min(data.dt) max(data.dt)]);
  case 'WSS'
    fig(60);
    clf
    scatter(data.dt, data.chl, 6, 'filled');
    ylabel('Chl (mg.m^{-3})');
    xlim([min(data.dt) max(data.dt)]);
  case 'ALFA'
    fig(70);
    clf
    subplot(4, 1, 1);
    yyaxis('left'); hold on
    scatter(data.dt, data.Chlb, 6, 'filled');
    scatter(data.dt, data.Chlg, 6, 'filled');
    ylabel('Chl (mg.m^{-1}??)');
    yyaxis('right'); hold on
    scatter(data.dt, data.CDOMRb, 6, 'filled');
    ylabel('CDOMRb (ppb??)');
    xlim([min(data.dt) max(data.dt)]);
    legend('Chlb', 'Chlg', 'CDOMRb')
    subplot(4, 1, 2); hold on
    scatter(data.dt, data.R613Rb, 6, 'filled');
    scatter(data.dt, data.R625Rb, 6, 'filled');
    scatter(data.dt, data.R642Rb, 6, 'filled');
    scatter(data.dt, data.R662Rb, 6, 'filled');
    scatter(data.dt, data.R642Rg, 6, 'filled');
    scatter(data.dt, data.R662Rg, 6, 'filled');
    ylabel('R...Rb & R...Rg (mg.m^{-1}??)');
    xlim([min(data.dt) max(data.dt)]);
    legend('R613Rb', 'R625Rb', 'R642Rb', 'R662Rb', 'R642Rg', 'R662Rg')
    subplot(4, 1, 3); hold on
    scatter(data.dt, data.PE1Rg, 6, 'filled');
    scatter(data.dt, data.PE2Rg, 6, 'filled');
    scatter(data.dt, data.PE3Rg, 6, 'filled');
    scatter(data.dt, data.PE1CFg, 6, 'filled');
    scatter(data.dt, data.PE2CFg, 6, 'filled');
    scatter(data.dt, data.PE3CFg, 6, 'filled');
    scatter(data.dt, data.PE12Rg, 6, 'filled');
    scatter(data.dt, data.PE12CFg, 6, 'filled');
    ylabel('Phycoerythryn?? (mg.m^{-1}??)');
    xlim([min(data.dt) max(data.dt)]);
    legend('PE1Rg', 'PE2Rg', 'PE3Rg', 'PE1CFg', 'PE2CFg', 'PE3CFg', 'PE12Rg', 'PE12CFg')
    subplot(4, 1, 4); hold on
    scatter(data.dt, data.FvFm, 6, 'filled');
    scatter(data.dt, data.FvFmC, 6, 'filled');
    scatter(data.dt, data.FvFmG, 6, 'filled');
    scatter(data.dt, data.FvFmCG, 6, 'filled');
    ylabel('FvFm [??]');
    xlim([min(data.dt) max(data.dt)]);
    legend('FvFm', 'FvFmC', 'FvFmG', 'FvFmCG')
  case 'LISST'
    fig(80);
    clf
    scatter(data.dt, data.cp, 6, 'filled');
    ylabel('cp m^{-1}');
    xlim([min(data.dt) max(data.dt)]);
    visProd3D(data.Properties.UserData.diameters, data.dt, data.PSD, false, 'Intensity', false, 101);
    set(gca, 'ZScale', 'log', 'XScale', 'log');
    zlabel(['PSD (' data.Properties.VariableUnits{strcmp(data.Properties.VariableNames, 'PSD')} ')']);
    xlabel('Diameters \mum');
    ylabel('Time');
    visProd3D(data.Properties.UserData.diameters, data.dt, data.betap, false, 'Intensity', false, 102);
    set(gca, 'ZScale', 'log', 'XScale', 'log');
    zlabel(['VSF (' data.Properties.VariableUnits{strcmp(data.Properties.VariableNames, 'betap')} ')']);
    xlabel('Angle [degree]');
    ylabel('Time');
  case {'LISSTTau', 'LISSTTAU', 'TAU'}
    fig(90);
    clf
    scatter(data.dt, data.Beamcp, 6, 'filled'); ylabel('c_p [m^{-1}]');
    xlim([min(data.dt) max(data.dt)]);
  otherwise
    warning('%s not supported for product visualisation', instrument)
end
function [data_corrected, DIW_biofouling_correction] = BiofoulingCorrection(data, lambda, var2correct, varcol, flow, spd_variable, fdom, previous_DIW_biofouling_correction)
  % GUI to select cleaning datetime and reference sections for biofouling correction of dissolved ACS data
  % Author: Guillaume Bourdin
  % date: 2024-11-18
  %
  %%
  if ~isdatetime(data.dt)
    datenum_bool = true;
    data.dt = datetime(data.dt, 'ConvertFrom', 'datenum');
  else
    datenum_bool = true;
  end
  if nargin < 7
    fdom = [];
  end
  % correct until satisfied
  var2run_again = var2correct;
  if nargin < 8
    DIW_biofouling_correction = table();
  elseif isempty(previous_DIW_biofouling_correction)
    DIW_biofouling_correction = table();
  else
    DIW_biofouling_correction = previous_DIW_biofouling_correction;
  end

  while ~all(strcmp(var2run_again, 'none'))
    [data_corrected, correction_table] = DIW_biofouling_reference(data, lambda, var2run_again, varcol, flow, spd_variable, fdom);
    if isempty(DIW_biofouling_correction)
      DIW_biofouling_correction = correction_table;
    else
      % add missing variable in old table to merge with new
      missing_var_old = find(~ismember(correction_table.Properties.VariableNames, ...
        DIW_biofouling_correction.Properties.VariableNames));
      if ~isempty(missing_var_old)
        for i = 1:size(missing_var_old, 2)
          DIW_biofouling_correction = addvars(DIW_biofouling_correction, NaT(size(DIW_biofouling_correction,1), 2), ...
            'NewVariableNames', correction_table.Properties.VariableNames(missing_var_old(i)), ...
            'Before', DIW_biofouling_correction.Properties.VariableNames{missing_var_old(i)});
        end
      end
      % add missing variable in new table to merge with old
      missing_var_new = find(~ismember(DIW_biofouling_correction.Properties.VariableNames, ...
        correction_table.Properties.VariableNames));
      if ~isempty(missing_var_new)
        for i = 1:size(missing_var_new, 2)
          correction_table = addvars(correction_table, NaT(size(correction_table,1), 2), ...
            'NewVariableNames', DIW_biofouling_correction.Properties.VariableNames(missing_var_new(i)), ...
            'Before', correction_table.Properties.VariableNames{missing_var_new(i)});
        end
      end
      % merge tables
      DIW_biofouling_correction = [DIW_biofouling_correction; correction_table];
    end
    var2cor = [var2correct {'none'}];
    [idx, correction_ok] = listdlg('PromptString',{'Select variable(s) to correct again.', ...
      'Select none to continue calibrate.', '', "Select 'Cancel' to ignore correction and "}, ...
      'ListString', var2cor, 'ListSize',[250, 75], 'Name','Correct more?','InitialValue', size(var2cor, 2));
    if isempty(idx)
      idx = size(var2cor, 2);
    end
    % if correction ok, run next round on the same data, otherwise go back to previous round data and correct same variables
    if correction_ok
      data = data_corrected;
      var2run_again = var2cor(idx);
    end
  end

  if datenum_bool
    data_corrected.dt = datenum(data_corrected.dt);
  end
  % visProd3D(lambda.a, data.dt(id_start_period:id_end_period), biofouling_2rm, false, 'Wavelength', false, 43);
  % visProd3D(lambda.a, data.dt(id_start_period:id_end_period), biofouling_2rm_exp, false, 'Wavelength', false, 44);
  % % visProd3D(lambda.a, data.dt(id_start_period:id_end_period), data.ag(id_start_period:id_end_period, :), false, 'Wavelength', false, 45);
  % % visProd3D(lambda.a, data_corrected.dt(id_start_period:id_end_period), data_corrected.ag(id_start_period:id_end_period, :), false, 'Wavelength', false, 46);
  % visProd3D(lambda.a, data.dt, data.ag, false, 'Wavelength', false, 47);
  % visProd3D(lambda.c, data.dt, data.cg, false, 'Wavelength', false, 48);
  % visProd3D(lambda.a, data_corrected.dt, data_corrected.ag, false, 'Wavelength', false, 49);
  % visProd3D(lambda.c, data_corrected.dt, data_corrected.cg, false, 'Wavelength', false, 50);
  % fh = visFlag([], [], data_corrected, [], [], [], var2correct{1}, varcol, [], flow, spd_variable, false, [], [], fdom);
end

function [data_corrected, correction_table] = DIW_biofouling_reference(data, lambda, var2correct, varcol, flow, spd_variable, fdom)
  data.dt = dateshift(data.dt, 'Start', 'Seconds');
  correction_table = table();
  % UI to select instrument cleaning time
  fh = visFlag([], [], data, [], [], [], var2correct{1}, varcol, [], flow, spd_variable, false, [], [], fdom);
  legend('Flow', 'Binned fCDOM', 'Binned a_g','AutoUpdate','off','FontSize',14)
  title('Select cleaning event datetime to correct (press s)', 'FontSize', 14)
  [~, ~, cleaning_dt, ~, ~] = guiSelectOnTimeSeries(fh);
  data_corrected = data;
  if ~isempty(cleaning_dt)
    % add variables
    correction_table.cleaning_dt = [min(data.dt) - minutes(1); cleaning_dt];
    frq = median(diff(data.dt));
    for v = 1:length(var2correct)
      correction_table.([var2correct{v} '_before']) = NaT(size(correction_table.cleaning_dt, 1), 2);
      correction_table.([var2correct{v} '_after']) = NaT(size(correction_table.cleaning_dt, 1), 2);
      correction_table.([var2correct{v} '_shift_all_period']) = NaT(size(correction_table.cleaning_dt, 1), 2);
      % Select period before and after cleaning to get delta_cleaning
      fh = visFlag([], [], data, [], [], [], var2correct{v}, varcol, [], flow, spd_variable, false, [], [], fdom);
      legend('Flow', 'Binned fCDOM', 'Binned a_g','AutoUpdate','off','FontSize',14)
      for i = 2:size(correction_table, 1)
        % automatically select last event before cleaning and second after cleaning as reference for biofouling correction
        id550 = lambda.(strrep(var2correct{v}, 'g', '')) <= 550;
        foo = data.dt - correction_table.cleaning_dt(i);
        id_before = false(size(data.dt));
        if any(foo < minutes(0))
          id_before(data.dt == dateshift(max(foo(foo < minutes(0))) + correction_table.cleaning_dt(i), 'Start', 'Seconds')) = true;
        end
        id_after = false(size(data.dt));
        if any(foo > minutes(0))
          foo_after = foo(foo > minutes(0));
          id_after = data.dt == dateshift(min(foo_after(2:end)) + correction_table.cleaning_dt(i), 'Start', 'Seconds');
        end
        if any(id_before) && any(id_after)
          if any(~isnan(data.(var2correct{v})(id_before, id550))) && any(~isnan(data.(var2correct{v})(id_after, id550))) && ...
              sum(data.(var2correct{v})(id_before, id550) > data.(var2correct{v})(id_after, id550)) > sum(id550)/2
            correction_table.([var2correct{v} '_before'])(i, :) = [data.dt(id_before) - minutes(2) data.dt(id_before) + minutes(2)];
            correction_table.([var2correct{v} '_after'])(i, :) = [data.dt(id_after) - minutes(2) data.dt(id_after) + minutes(2)];
            foo_ylim = ylim;
            yyaxis('right')
            plot(data.dt(id_before | id_after), data.(var2correct{v})(id_before | id_after, varcol), 'o', 'Color', 'r', 'MarkerFaceColor', 'r')
            ylim(foo_ylim)
          end
        end
        title(['Cleaning period [' num2str(i) '/' num2str(size(correction_table, 1)) ']. Variable: ' var2correct{v} newline ...
          'Select period before (press F) and after (press T) instrument cleaning (red dots are automatically selected (press Q to save each entry), user input will overwrite automatic selection)' newline ...
          '\Delta_{cleaning} = x_{before cleaning} - x_{after cleaning}' newline ...
          'Select entire period (press X) to shift and match period selected after cleaning.'], 'FontSize', 12)
        % force axes limits to zoom on each cleaning event
        xlim([correction_table.cleaning_dt(i-1)-6*frq correction_table.cleaning_dt(i)+24*frq])
        ylim('auto')
        [aftr, bfr, ~, ~, shift_all_period] = guiSelectOnTimeSeries(fh);
        if ~isempty(bfr)
          correction_table.([var2correct{v} '_before'])(i, :) = bfr(end, :);
        end
        if ~isempty(aftr)
          correction_table.([var2correct{v} '_after'])(i, :) = aftr(end, :);
        end
        if ~isempty(shift_all_period)
          correction_table.([var2correct{v} '_shift_all_period'])(i, :) = shift_all_period(end, :);
        end
      end
      % fit biofouling with linear function
      for j = 2:size(correction_table, 1)
        if all(~isnat(correction_table.([var2correct{v} '_shift_all_period'])(j,:)))
          id_all_period = data.dt > correction_table.([var2correct{v} '_shift_all_period'])(j,1) & data.dt < correction_table.([var2correct{v} '_shift_all_period'])(j,2);
          if all(isnat(correction_table.([var2correct{v} '_before'])(j,:)) & isnat(correction_table.([var2correct{v} '_after'])(j,:)))
            foo = data.dt(id_all_period);
            correction_table.([var2correct{v} '_before'])(j,1) = foo(2)-minutes(2);
            correction_table.([var2correct{v} '_before'])(j,2) = foo(2)+minutes(2);
            foo = data.dt(data.dt < correction_table.([var2correct{v} '_shift_all_period'])(j,1));
            correction_table.([var2correct{v} '_after'])(j,1) = foo(end)-minutes(2);
            correction_table.([var2correct{v} '_after'])(j,2) = foo(end)+minutes(2);
          end
          % find sections before/after cleaning
          id_before = data.dt > correction_table.([var2correct{v} '_before'])(j,1) & data.dt < correction_table.([var2correct{v} '_before'])(j,2);
          id_after = data.dt > correction_table.([var2correct{v} '_after'])(j,1) & data.dt < correction_table.([var2correct{v} '_after'])(j,2);
          % compute delta biofouling as values before - after
          delta_cleaning = mean(data.(var2correct{v})(id_before, :), 1, 'omitnan') - mean(data.(var2correct{v})(id_after, :), 1, 'omitnan');
          % remove delta_cleaning
          data_corrected.(var2correct{v})(id_all_period, :) = data.(var2correct{v})(id_all_period, :) - delta_cleaning;
        elseif all(~isnat(correction_table.([var2correct{v} '_before'])(j,:)) & ~isnat(correction_table.([var2correct{v} '_after'])(j,:)))
          % find sections before/after cleaning
          id_before = data.dt > correction_table.([var2correct{v} '_before'])(j,1) & data.dt < correction_table.([var2correct{v} '_before'])(j,2);
          id_after = data.dt > correction_table.([var2correct{v} '_after'])(j,1) & data.dt < correction_table.([var2correct{v} '_after'])(j,2);
          if any(id_before) && any(id_after)
            % find start and end of section to correct for biofouling
            id_start_period = find(data.dt == min(data.dt(data.dt > correction_table.cleaning_dt(j-1))));
            id_end_period = find(data.dt == max(data.dt(data.dt < correction_table.cleaning_dt(j))));
            
            % EXPONENTIAL FIT VERSION, does not work for now
            % % fit exponential function to biofouling
            % biofouling_2rm_exp = NaN(size(data.(var2correct{v})(id_start_period:id_end_period, :)));
            % for kk = progress(1:size(data.(var2correct{v}), 2))
            %   foo_dt = data.dt(id_start_period:id_end_period) - min(data.dt(id_start_period:id_end_period));
            %   foo = data.(var2correct{v})(id_start_period:id_end_period, kk);
            %   foo_dt = minutes(foo_dt(~isnan(foo)));
            %   foo = foo(~isnan(foo));
            % 
            %   delta_neg = min(foo);
            %   foo = foo - delta_neg;
            %   f0 = fit(foo_dt, foo, 'exp2', 'Normalize', 'on', 'Robust', 'Bisquare');
            %   biof = f0(foo_dt);
            %   delta_start = interp1([foo_dt(1) foo_dt(end)], [biof(1) 0], foo_dt, 'linear'); 
            %   biof = biof - delta_start;
            %   biofouling_2rm_exp(:, kk) = biof + delta_neg;
            %   % figure(8); hold on;
            %   % plot(minutes(foo_dt(~isnan(foo))), biofouling_2rm_exp(:, kk));
            %   % scatter(minutes(foo_dt(~isnan(foo))), foo + delta_neg)
            %   % pause(1)
            % end
            
            % get section duration in minutes
            duration_minutes = minutes(data.dt(id_start_period:id_end_period) - data.dt(id_start_period));
            % compute delta_biofouling
            delta_biofouling = mean(data.(var2correct{v})(id_before, :), 1, 'omitnan') - mean(data.(var2correct{v})(id_after, :), 1, 'omitnan');
            % compute biofouling per minutes
            biofouling_per_min = delta_biofouling ./ max(duration_minutes);
            % reconstruct linear biofouling as biofouling_per_min .* duration_minutes
            biofouling_2rm = biofouling_per_min .* duration_minutes;
            % remove estimated biofouling
            data_corrected.(var2correct{v})(id_start_period:id_end_period, :) = data.(var2correct{v})(id_start_period:id_end_period, :) - biofouling_2rm;
            % data_corrected.(var2correct{v})(id_start_period:id_end_period, :) = data.(var2correct{v})(id_start_period:id_end_period, :) - biofouling_2rm_exp;
          end
        end
      end
    end
  end
end

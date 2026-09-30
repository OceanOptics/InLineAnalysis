function [ad, aph] = partition_ap_vec(ap_i, lambda, prcntl, verbose)
  % Description: This code partitions the spectral total particulate absorption coefficient, ap(lambda) into
  % phytoplankton, aph(lambda), and non-algal, ad(lambda), components.
  %
  % Vectorized version of partition_ap.m. Test against original code with
  % differences only within rounding noise range < 10e-15 1/m.
  % What makes this fast (in order of importance:
  %   1. FAST CROSSING LOCALIZATION (the dominant cost): finding where
  %      fres(s) = A*Bcurve(s) - D*Ccurve(s) crosses zero is equivalent to
  %      finding where the FIXED curve R(s) = Bcurve(s)/Ccurve(s) crosses
  %      the per-spectrum TARGET value D/A -- R only depends on the (x,y)
  %      grid point, not on the spectrum. When R is monotonic over s (true
  %      at ~84% of grid points, checked once per grid point up front),
  %      MATLAB's discretize() locates the bracketing interval for all K
  %      spectra in one vectorized call, instead of building a full K-by-131
  %      array and scanning it. The interpolation itself is then still done
  %      with the EXACT original fres-based formula (evaluated at just the
  %      2 bracketing points, gathered per spectrum), so results are
  %      bit-identical to the original method -- only the localization step
  %      changes, not the arithmetic that produces Ad/Sd. Grid points where
  %      R isn't monotonic (or Ccurve changes sign) fall back to the exact
  %      original per-spectrum scan.
  %   2. The 100-by-100 (x,y) grid search runs ONCE for ALL spectra together.
  %   3. STAGED / LAZY EVALUATION: 5 of the 6 stacked constraints only ever
  %      need 6 specific wavelengths, not the full hyperspectral curve.
  %      Only the (typically ~1-2%) survivors of those cheap checks get the
  %      expensive full-spectrum constraints (#11-13) and positivity check.
  %
  % Input:  ap_i, spectra of total particulate absorption coefficient,
  %               K-by-N matrix: K spectra (rows) by N wavelengths (columns).
  %               A single spectrum as a plain row vector also works (K=1).
  %         lambda, wavelengths of ap_i, must be a column vector, and must be hyperspectral from 400-700 nm
  %         prcntl, percentiles of feasible solutions, can be a scalar or vector
  %         verbose (optional): if true, prints a timing breakdown
  % Output: ad, non-algal (detrital) particulate absorption spectra, K-by-N-by-numel(prcntl)
  %         aph, phytoplankton specific absorption spectra, K-by-N-by-numel(prcntl)
  %         (same K-by-N row/column convention as ap_i, with percentile as the trailing dimension)
  %
  % If no feasible solutions are found for a given spectrum, its output values are NaNs.
  %
  % Performance note: the outer x-loop inside scm_ap (variable i, 100
  % iterations) is fully independent across i and can be changed to
  % "parfor i = 1:m" (requires Parallel Computing Toolbox)
  %
  % Author:  Guangming Zheng, gzheng@ucsd.edu (original single-spectrum algorithm)
  % Vectorized by: (completed from a stalled work-in-progress)
  %
  % An example to use this code:
  %
  % [ad_scm, aph_scm] = partition_ap_vec(ap_input, lambda_input, [10 50 90]);
  %             , where lambda_input is a column vector, [400:700]'
  %                     ap_input is K-by-N (K spectra, N wavelengths) at the same wavelengths as lambda_input
  %                     [10 50 90] specifies the output values as the 10th, 50th (median, i.e., optimal solution), and 90th percentiles of all feasible solutions.
  %                                If you only want the optimal solution, just ignore this input option.
  %
  % Reference:
  % Zheng, G., and D. Stramski (2013), A model based on stacked-constraints
  % approach for partitioning the light absorption coefficient of seawater
  % into phytoplankton and non-phytoplankton components, J. Geophys. Res. Oceans,
  % 118, 2155?2174, doi:10.1002/jgrc.20115.
  %
  % Zheng, G., and D. Stramski (2013), A model for partitioning the light
  % absorption coefficient of suspended marine particles into phytoplankton
  % and nonalgal components, J. Geophys. Res. Oceans, 118, 2977?2991,
  % doi:10.1002/jgrc.20206.
  %
  %%
  if nargin < 4; verbose = false; end
  [K, N] = size(ap_i);
  [Ad_c, Sd_c] = scm_ap(ap_i, lambda, verbose); % Ad_c, Sd_c: P-by-K (P = max feasible count across all spectra)
  P = size(Ad_c, 1);

  prcntl = prcntl(:)'; % row vector
  nP = numel(prcntl);
  ad  = NaN(K, N, nP);
  aph = NaN(K, N, nP);
  if P == 0
    return; % no spectrum has any feasible solution at all
  end

  if verbose; t0 = tic; end
  for w = 1:N
    adtemp_w  = Ad_c .* exp(-Sd_c * lambda(w));  % P-by-K, NaN where infeasible
    aphtemp_w = ap_i(:, w)' - adtemp_w;           % P-by-K
    % prctile ignores NaN values by default, so spectra with fewer
    % feasible grid points (or none at all) are handled automatically
    ad(:, w, :)  = permute(prctile(adtemp_w, prcntl, 1), [2 3 1]);  % nP-by-K -> K-by-1-by-nP
    aph(:, w, :) = permute(prctile(aphtemp_w, prcntl, 1), [2 3 1]);
  end
  if verbose; fprintf('percentile aggregation over %d wavelengths (%d feasible grid points per spectrum, max): %.2f s\n', N, P, toc(t0)); end
end


function [Ad_c, Sd_c] = scm_ap(ap_i, lambda, verbose)
  % Input:  ap_i, spectra of total particulate absorption coefficient, K-by-N
  %         lambda, wavelengths of ap_i, N-by-1
  % Output: Ad_c, Sd_c, feasible solutions, compacted to P-by-K (P = the
  %         largest number of feasible grid points found for any single
  %         spectrum), NaN-padded for spectra with fewer. This never
  %         materializes the full m-by-n-by-K speculative grid (100-by-100-by-K
  %         doubles, which would be tens of GB at K in the hundreds of
  %         thousands) since typically only ~1-3% of grid points are ever
  %         feasible for a given spectrum -- feasible (Ad,Sd,spectrum)
  %         triples are accumulated directly during the grid scan instead.
  %%
  % find the index of wavelengths involved in all constraints (numeric
  % indices throughout, matching partition_ap.m -- idx670 is the one
  % exception left as a logical mask, which is harmless for single-element
  % selection and matches the original)
  idx400 = find(lambda==400);
  idx412 = find(lambda==412);
  idx420 = find(lambda==420);
  idx430 = find(lambda==430);
  idx443 = find(lambda==443);
  idx450 = find(lambda==450);
  idx467 = find(lambda==467);
  idx490 = find(lambda==490);
  idx500 = find(lambda==500);
  idx510 = find(lambda==510);
  idx550 = find(lambda==550);
  idx555 = find(lambda==555);
  idx630 = find(lambda==630);
  idx650 = find(lambda==650);
  idx670 = find(lambda==670); % numeric here (needed as a landmark column index below)
  idx700 = find(lambda==700);
  
  % create vectors representing all possible values of aph ratios
  x = .01:.01:1;   m = numel(x);
  y = .01:.01:1;   n = numel(y);
  K = size(ap_i, 1);
  ss = 0.005:.0001:.018; % Constraint #3
  lam4 = [443 412 490 510];
  
  % upper bound of constraint #5
  ubrph12 = .38 + ap_i(:, idx467) ./ ap_i(:, idx412); % K-by-1
  
  % upper bound of constraint #7
  ubrph45 = 2.3 * ap_i(:, idx510) ./ ap_i(:, idx555) - 1.3; % K-by-1
  
  aps = ap_i(:, [idx443 idx412 idx490 idx510]); % K-by-4
  
  % --- precompute everything that depends on only x or only y, once ---
  A_all = aps(:, 2) - x .* aps(:, 1);         % K-by-m,  A_all(:,i) for x(i)
  D_all = aps(:, 4) - y .* aps(:, 3);         % K-by-n,  D_all(:,j) for y(j)
  Bcurve_all = exp(-ss * lam4(4)) - y' * exp(-ss * lam4(3)); % n-by-numel(ss)
  Ccurve_all = exp(-ss * lam4(2)) - x' * exp(-ss * lam4(1)); % m-by-numel(ss)
  
  % --- precompute, once, which grid points allow the fast crossing
  % localization: Ccurve must not change sign over s, and R=Bcurve/Ccurve
  % must be monotonic over s at that grid point (checked directly from the
  % already-precomputed curves, does not depend on K at all) ---
  if verbose; t0 = tic; end
  CsignOK = all(Ccurve_all > 0, 2) | all(Ccurve_all < 0, 2); % m-by-1
  fastOK = false(m, n);
  for i = 1:m
    if ~CsignOK(i); continue; end
    Ci = Ccurve_all(i, :);
    for j = 1:n
      Rij = Bcurve_all(j, :) ./ Ci;
      dR = diff(Rij);
      fastOK(i, j) = all(dR > 0) || all(dR < 0);
    end
  end
  if verbose; fprintf('fast-path eligibility precompute: %.2f s (%.1f%% of grid points eligible)\n', toc(t0), 100*sum(fastOK(:))/(m*n)); end
  
  % --- landmark wavelengths: constraints #4/5, #6/7, #8, #9, #10 only
  % ever need these 6 columns, never the full spectrum ---
  landIdx = [idx412 idx443 idx467 idx510 idx555 idx670];
  landLam = lambda(landIdx)'; % 1-by-6, actual wavelength values (for the exp() formula)
  apLand  = ap_i(:, landIdx); % K-by-6, precomputed once
  x0const = ap_i(:, idx412) ./ ap_i(:, idx443); % K-by-1, constraint #9's x0 -- independent of grid point
  
  % --- growable buffer for the (Ad, Sd, spectrum-index) triples of
  % feasible grid points, accumulated directly during the scan below,
  % instead of ever storing an m-by-n-by-K array ---
  bufCap = 1e6;
  feasAd = zeros(bufCap, 1);
  feasSd = zeros(bufCap, 1);
  feasK  = zeros(bufCap, 1);
  nFeas  = 0;
  
  % timing instrumentation (negligible overhead, only prints if verbose)
  t_cross = 0; t_stage1 = 0; t_stage2 = 0; nSurvTotal = 0; nFastUsed = 0;
  
  for i = 1:m
    xi = x(i);
    Ai = A_all(:, i);        % K-by-1
    Ci = Ccurve_all(i, :);   % 1-by-numel(ss)
    for j = 1:n
      yi = y(j);
      Dj = D_all(:, j);       % K-by-1
      Bj = Bcurve_all(j, :);  % 1-by-numel(ss)
  
      if verbose; t0 = tic; end
      Ad_xy = NaN(K, 1);
      Sd_xy = NaN(K, 1);
      usedFast = false;
  
      if fastOK(i, j)
        try
          % ---- FAST PATH: locate the (guaranteed <= 1) crossing via the
          % monotonic curve R = Bj./Ci, then interpolate with the EXACT
          % same fres-based formula as the original, evaluated only at the
          % 2 bracketing points (gathered per spectrum) ----
          Rij = Bj ./ Ci; % 1-by-numel(ss)
          target = Dj ./ Ai; % K-by-1
  
          if Rij(1) < Rij(end)
            c1 = discretize(target, Rij); % 1-based bracketing index, NaN if no crossing
          else
            Rrev = fliplr(Rij);
            c1rev = discretize(target, Rrev);
            c1 = numel(Rij) - c1rev;
          end
  
          validC = ~isnan(c1);
          if any(validC)
            c0 = c1(validC);
            Bc  = Bj(c0)';   Bc1 = Bj(c0 + 1)';
            Cc  = Ci(c0)';   Cc1 = Ci(c0 + 1)';
            Av  = Ai(validC); Dv = Dj(validC);
            fres0 = Av .* Bc  - Dv .* Cc;
            fres1 = Av .* Bc1 - Dv .* Cc1;
            y0v = ss(c0)'; y1v = ss(c0 + 1)';
            Sd_v = y0v + (0 - fres0) .* (y1v - y0v) ./ (fres1 - fres0);
            Ad_v = Av ./ (exp(-lam4(2) * Sd_v) - xi * exp(-lam4(1) * Sd_v));
            idxValid = find(validC);
            Sd_xy(idxValid) = Sd_v;
            Ad_xy(idxValid) = Ad_v;
          end
          usedFast = true;
        catch
          usedFast = false; % fall through to the exact method below
        end
      end
  
      if ~usedFast
        % ---- EXACT method (also the fallback if the fast path errors for
        % any reason) ----
        fres = Ai * Bj - Dj * Ci; % K-by-numel(ss)
  
        pos = fres > 0;
        foo = diff(pos, [], 2);
        nz = foo ~= 0;
        nCross = sum(nz, 2);
  
        single = nCross == 1;
        if any(single)
          rows = find(single);
          [~, c] = max(nz(rows, :), [], 2);
          lin0 = sub2ind(size(fres), rows, c);
          lin1 = sub2ind(size(fres), rows, c + 1);
          x0v = fres(lin0);  x1v = fres(lin1);
          y0v = ss(c)';      y1v = ss(c + 1)';
          Sd_s = y0v + (0 - x0v) .* (y1v - y0v) ./ (x1v - x0v);
          Ad_s = Ai(rows) ./ (exp(-lam4(2) * Sd_s) - xi * exp(-lam4(1) * Sd_s));
          Sd_xy(rows) = Sd_s;
          Ad_xy(rows) = Ad_s;
        end
  
        other = find(~single & nCross > 0);
        for k = other'
          idxc = find(nz(k, :));
          if numel(idxc) == 2 && abs(idxc(1) - idxc(2)) == 1
            continue
          end
          try
            Sd_k = interp1([fres(k, idxc) fres(k, idxc + 1)], [ss(idxc) ss(idxc + 1)], 0);
          catch
            fprintf("Warning: solver4s_ap's interpolation crashed for spectrum %i, replaced by NaN\n", k)
            continue
          end
          Sd_xy(k) = Sd_k;
          Ad_xy(k) = Ai(k) / (exp(-lam4(2) * Sd_k) - xi * exp(-lam4(1) * Sd_k));
        end
      end
  
      if verbose; t_cross = t_cross + toc(t0); nFastUsed = nFastUsed + usedFast; end
  
      valid = ~isnan(Ad_xy) & ~isnan(Sd_xy);
      if ~any(valid)
        continue
      end
  
      if verbose; t0 = tic; end
      % ==== STAGE 1 (cheap): constraints #4/5, #6/7, #8, #9, #10, using
      % ONLY the 6 landmark wavelengths, for every valid spectrum at once ====
      vrows = find(valid);
      ad_small  = Ad_xy(vrows) .* exp(-Sd_xy(vrows) .* landLam); % nValid-by-6
      aph_small = apLand(vrows, :) - ad_small;                    % nValid-by-6
      aph412 = aph_small(:,1); aph443c = aph_small(:,2); aph467 = aph_small(:,3);
      aph510 = aph_small(:,4); aph555  = aph_small(:,5); aph670  = aph_small(:,6);
  
      rph12 = aph467 ./ aph412;
      f1 = rph12 < 1.54 & rph12 < ubrph12(vrows) & rph12 > .74;
  
      rph45 = aph510 ./ aph555;
      f2 = rph45 < 10 & rph45 < ubrph45(vrows) & rph45 > 1.3;
  
      b2r = aph443c ./ aph670;
      f3 = b2r > 1.4 & b2r < 9.1;
  
      y0c = ad_small(:,1) ./ ap_i(vrows, idx412); % ad at 412 / ap at 412
      f4 = y0c < x0const(vrows) - .45 & y0c > x0const(vrows) - .91;
  
      naph510 = aph510 ./ aph443c; naph555 = aph555 ./ aph443c;
      sph45 = (naph510 - naph555) / (555 - 510);
      f5 = sph45 < 8.7e-3 & sph45 > 3e-3;
  
      % necessary (not sufficient) positivity check at the 6 landmarks --
      % the original requires positivity at ALL 301 wavelengths, so failing
      % here means it will definitely fail the full check later too
      posLand = all(aph_small > 0, 2) & all(ad_small > 0, 2);
  
      cheapPass = f1 & f2 & f3 & f4 & f5 & posLand;
      if verbose; t_stage1 = t_stage1 + toc(t0); end
      if ~any(cheapPass)
        continue
      end
  
      if verbose; t0 = tic; nSurvTotal = nSurvTotal + sum(cheapPass); end
      % ==== STAGE 2 (expensive): only for the cheap-stage survivors,
      % build the full hyperspectral ad/aph for constraint #11-13 (which
      % need wavelength ranges) and the full non-negativity requirement ====
      survivors = vrows(cheapPass);
      ad_temp  = Ad_xy(survivors) .* exp(-Sd_xy(survivors) .* lambda'); % nSurv-by-N
      aph_temp = ap_i(survivors, :) - ad_temp;
  
      g2rsum = sum(aph_temp(:, idx550:idx630), 2) ./ sum(aph_temp(:, idx650:idx700), 2); % Constraint parameter #11
      b2gsum = sum(aph_temp(:, idx450:idx490), 2) ./ sum(aph_temp(:, idx510:idx550), 2); % Constraint parameter #12
  
      naph = aph_temp ./ aph_temp(:, idx443);
      [rmax, irmax] = max(naph(:, idx450:idx500), [], 2);
      [rmin, irmin] = min(naph(:, idx500:idx550), [], 2);
      bgslp = (rmax - rmin) ./ (51 - irmax + irmin);
      [lmax, ilmax] = max(naph(:, idx430:idx450), [], 2);
      [lmin, ilmin] = min(naph(:, idx400:idx420), [], 2);
      bvslp = (lmax - lmin) ./ (30 + ilmax - ilmin);
  
      difslp = bgslp - bvslp; % Constraint parameter #13
      f6 = g2rsum < 1.28 & b2gsum > 2 & difslp > -.0047;
  
      posReq = all(aph_temp > 0, 2) & all(ad_temp > 0, 2);
  
      feasibleLocal = f6 & posReq;
      if any(feasibleLocal)
        flagRows = survivors(feasibleLocal);
        nNew = numel(flagRows);
        if nFeas + nNew > numel(feasAd)
          newCap = max(numel(feasAd) * 2, nFeas + nNew);
          feasAd(newCap) = 0; %#ok<AGROW> % auto-grows; new slots overwritten immediately below
          feasSd(newCap) = 0;
          feasK(newCap)  = 0;
        end
        feasAd(nFeas+1:nFeas+nNew) = Ad_xy(flagRows);
        feasSd(nFeas+1:nFeas+nNew) = Sd_xy(flagRows);
        feasK(nFeas+1:nFeas+nNew)  = flagRows;
        nFeas = nFeas + nNew;
      end
      if verbose; t_stage2 = t_stage2 + toc(t0); end
    end
  end
  
  % trim to actual size, then rebuild the compact P-by-K array from the
  % (Ad, Sd, spectrum) triples: sort by spectrum index (groups each
  % spectrum's own points contiguously), then compute each triple's
  % position within its own spectrum's group via a cumulative-offset trick
  % (no explicit per-spectrum loop needed)
  feasAd = feasAd(1:nFeas);
  feasSd = feasSd(1:nFeas);
  feasK  = feasK(1:nFeas);
  
  if nFeas == 0
    Ad_c = zeros(0, K);
    Sd_c = zeros(0, K);
  else
    countPerSpec = accumarray(feasK, 1, [K, 1]);
    P = max(countPerSpec);
    [sortedK, sortOrder] = sort(feasK);
    sortedAd = feasAd(sortOrder);
    sortedSd = feasSd(sortOrder);
    startOffset = [0; cumsum(countPerSpec(1:end-1))]; % cumulative count before spectrum k
    posInGroup = (1:nFeas)' - startOffset(sortedK);     % 1-based position within each spectrum's block
    linIdx = sub2ind([P, K], posInGroup, sortedK);
    Ad_c = NaN(P, K);
    Sd_c = NaN(P, K);
    Ad_c(linIdx) = sortedAd;
    Sd_c(linIdx) = sortedSd;
  end
  
  if verbose
    fprintf('--- scm_ap timing breakdown ---\n');
    fprintf('crossing localization + root solve: %.2f s (fast path used at %.1f%% of grid points)\n', t_cross, 100*nFastUsed/(m*n));
    fprintf('stage 1 (cheap, 6-wavelength constraints):          %.2f s\n', t_stage1);
    fprintf('stage 2 (expensive, full-spectrum constraints):     %.2f s (on %d total survivor-gridpoint pairs)\n', t_stage2, nSurvTotal);
    fprintf('sum of the above (loop overhead not included):      %.2f s\n', t_cross+t_stage1+t_stage2);
  end
end

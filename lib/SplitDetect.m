function [FTH] = SplitDetect (instrument, data, FTH, MinFiltPeriod, szFilt)
% author: Guillaume Bourdin
% created: Oct 10, 2019

%Detect filter events on ACS and BB3 inline
%
% INPUT:
%   - instrument cell string of instrument name e.g. 'BB3' or 'ACS007'
%   - data_in: instrument data <N Table> time series of data that must contains:
%         - dt <Nx1 datenum> date and time precise to the second
%         and either of the following variable must exist:
%         - beta <NxM double> in column M for BB3
%         - a <NxM double> in column M for ACS
%   - FTH log <Nx3 double> time series of data that must contains:
%         - dt <Nx1 datenum> date and time precise to the second
%         - swt <Nx1 double> 0 and 1 values for switch position
%                (0 = total, 1 = dissolved)
%         - spd <Nx1 double> flow rate (in lpm)
%
% OUTPUT:
%   - data_out: modified FTH log <Nx3 double> time series of data that must contains:
%         - dt <Nx1 datenum> date and time precise to the second
%         - swt <Nx1 double> 0 and 1 values for switch position
%                (0 = total, 1 = dissolved) automatically detected and
%                synchronised
%         - spd <Nx1 double> flow rate (in lpm)

if nargin < 2
  error('missing argument: instrument, data and FTH');
elseif nargin < 4
  warning('Optional arguments missing, set to default: MinFiltPeriod = 65 | szFilt = 10')
  MinFiltPeriod = 65;
  szFilt = 10;
elseif nargin < 5
  warning('Optional arguments missing, set to default: szFilt = 10')
  szFilt = 10;
end

% if contains(instrument, 'ACS')
%     sampling_freq = 4;
% elseif contains(instrument, 'AC9')
%     sampling_freq = 6;
% elseif contains(instrument, 'BB3')
%     sampling_freq = 1;
% end
% raw_swt = FTH.swt; % back-up FTH log

%% create new FTH vector
if ~isempty(FTH)
  % FTH.swt = zeros(size(FTH.swt,1),1);
  oldFTH = FTH;
  FTHstart = min([oldFTH.dt; data.dt]);
  FTHend = max([oldFTH.dt; data.dt]);
  FTHdt = (FTHstart:1/1/3600/24:FTHend)';
%   FTHdt = (oldFTH.dt(1):1/1/3600/24:oldFTH.dt(end))';

  % delete data with duplicats timestamp
  [~, I] = unique(oldFTH.dt, 'first');
  x = 1:length(oldFTH.dt);
  x(I) = [];
  oldFTH(x,:) = [];
  int_spd = interp1(oldFTH.dt, oldFTH.spd, FTHdt, 'linear');
  ism = ~ismember(FTHdt,oldFTH.dt);
  int_spd (ism) = NaN;
  FTH = table(FTHdt, zeros(size(FTHdt,1),1), int_spd, 'VariableNames', {'dt', 'swt', 'spd'});
else
  FTH = table(data.dt, zeros(size(data.dt,1),1), zeros(size(data.dt,1),1), ...
    'VariableNames', {'dt', 'swt', 'spd'});
end

%%
dt_ini = datetime(table2array(data(:,1)),'ConvertFrom','datenum');
if any(contains(instrument, 'ACS') | contains(instrument, 'AC9'))
  data = data.a; data(isinf(data)) = NaN;
  m_data_ini = mean(data, 2, 'omitnan');
elseif contains(instrument, 'BB')
  data = data.beta; data(isinf(data)) = NaN;
  m_data_ini = mean(data, 2, 'omitnan');
else
  error('Instrument not recognized, check "QC Reference" in cfg file');
end

d = designfilt('lowpassiir','FilterOrder',2, ...
  'PassbandFrequency',1/(60*5),'PassbandRipple',0.2, ...
  'SampleRate', 1);

sep = diff(dt_ini);
seg = find(sep > minutes(10)); seg = [1; seg; size(dt_ini,1)];

for i = progress(1:size(seg,1)-1)
  m_data = m_data_ini(seg(i):seg(i+1));
  dt = dt_ini(seg(i):seg(i+1));

  if any(isnan(m_data))
    filled = fillmissing(m_data,'linear','SamplePoints',dt);
  else
    filled = m_data;
  end
  m_M= filtfilt(d,filled);

    
%   try
%     if any(contains(instrument, 'ACS') | contains(instrument, 'AC9'))
% %       m_M = movmean(m_data,200);
%       m_M = sgolayfilt(m_data,1,201);
%     elseif contains(instrument, 'BB3')
%         m_M = movmean(m_data,200);
%     end
%   catch
%     if any(contains(instrument, 'ACS') | contains(instrument, 'AC9'))
%       m_M = movmean(m_data,200);
%     elseif contains(instrument, 'BB3')
%       m_M = movmean(m_data,200);
%     end
%   end

%   dbbp = NaN(size(m_M,1),1);
%   for j = 1:size(m_M,1)-1; dbbp (j) = m_M(j)-m_M(j+1); end
  dbbp = [0; -diff(m_M)];
  m_dbbp = movmean(dbbp,200);

  try
    if contains(instrument, 'ACS')
%       [~, x] = findpeaks(-movmean(m_M,1000),'MinPeakDistance',7000,'MinPeakWidth',1200); % 10000
      [~, x] = findpeaks(-sgolayfilt(m_M,1,401),'MinPeakDistance',MinFiltPeriod*0.6*240); % 10000  ,'MinPeakWidth',1200
    elseif contains(instrument, 'AC9')
%       [~, x] = findpeaks(-movmean(m_M,1000),'MinPeakDistance',7000,'MinPeakWidth',1200); % 10000
      [~, x] = findpeaks(-sgolayfilt(m_M,1,601),'MinPeakDistance',MinFiltPeriod*0.6*360); % 10000  ,'MinPeakWidth',1200
    elseif contains(instrument, 'BB3')
      [~, x] = findpeaks(-m_M,'MinPeakDistance',MinFiltPeriod*0.47*60,'MinPeakWidth',200); % 2000
    elseif contains(instrument, 'HBB')
      [~, x] = findpeaks(-m_M,'MinPeakDistance',MinFiltPeriod*2,'MinPeakWidth',20); % 2000
      smooth_m_data = movmedian(m_data,40);
    end
  catch
    if any(contains(instrument, 'ACS') | contains(instrument, 'AC9'))
%       [~, x] = findpeaks(-movmean(m_M,1000),'MinPeakDistance',size(m_data,1)-2,'MinPeakWidth',1200); % 10000
      [~, x] = findpeaks(-movmean(m_M,1000),'MinPeakDistance',size(m_data,1)-2); % 10000  ,'MinPeakWidth',1200
    elseif contains(instrument, 'AC9')
%       [~, x] = findpeaks(-movmean(m_M,1000),'MinPeakDistance',7000,'MinPeakWidth',1200); % 10000
      [~, x] = findpeaks(-movmean(m_M,1000),'MinPeakDistance',size(m_data,1)-2); % 10000  ,'MinPeakWidth',1200
    elseif contains(instrument, 'BB3')
      [~, x] = findpeaks(-m_M,'MinPeakDistance',size(m_data,1)-2,'MinPeakWidth',200); % 2000
    elseif contains(instrument, 'HBB')
      [~, x] = findpeaks(-m_M,'MinPeakDistance',size(m_data,1)-2,'MinPeakWidth',20); % 2000
      smooth_m_data = movmedian(m_data,40);
    end
  end

%     figure(62)
%     yyaxis('left')
%     hold on;
%     scatter(dt, m_data, 6, 'g', 'filled')
%     scatter(dt,-m_M,3,'r','filled');
% %     scatter(dt,-movmean(m_M,1000),2,'filled');
% %     scatter(dt,sgolayfilt(m_M,1,401),2,'g','filled');
%     scatter(dt(x), median(m_M, 'omitnan')*ones(size(x,1),1),30,'r','filled')
%     yyaxis('right')
%     scatter(dt,m_dbbp,3,'k','filled')
%     scatter(dt,dbbp,2,'b','filled')

% ACS conversion from szFilt to nb of points (initial method)
% szFilt*60*4/1.25 % 1920
% szFilt*60*4/1.92 % 1250

  filt_st = NaT(size(x,1),1);
  filt_end = NaT(size(x,1),1);
  if contains(instrument, 'ACS') % for ACS
    parfor k = 1:size(x,1)
      if all(x(k)<=round(szFilt*60*4/1.25) & x(k)>=size(m_data,1)-round(szFilt*60*4/1.25))% when seg is very small
        outlim = [1 size(m_M,1)];
      elseif x(k)<=round(szFilt*60*4/1.25) % when filter event in the beginning of the time series
        outlim = [1 x(k)+round(szFilt*60*4/1.25)];
      elseif x(k)>=size(m_data,1)-round(szFilt*60*4/1.25) % when filter event at the end of the time series
        outlim = [x(k)-round(szFilt*60*4/1.25) size(m_M,1)];
      else % when filter event in the middle of the time series
        outlim = [x(k)-round(szFilt*60*4/1.25) x(k)+round(szFilt*60*4/1.25)];
      end
      sel = m_M(outlim(1):outlim(2));
      seldt = dt(outlim(1):outlim(2));
      cent = median(seldt(sel<mean(sel, 'omitnan')));
      idx = find(abs(dt-cent)<seconds(1));
      if isempty(idx)
        idx = floor(median(find(abs(dt-cent)<seconds(3))));
      end

      if all(idx(1)<=round(szFilt*60*4/1.92) & idx(1)>=size(m_data,1)-round(szFilt*60*4/1.92))% when seg is very small
        inlim = [1 size(m_M,1)];
      elseif idx(1)<=round(szFilt*60*4/1.92) % when filter event in the beginning of the time series
        inlim = [1 idx(1)+round(szFilt*60*4/1.92)];
      elseif idx(1)>=size(m_data,1)-round(szFilt*60*4/1.92) % when filter event at the end of the time series
        inlim = [idx(1)-round(szFilt*60*4/1.92) size(m_M,1)];
      else % when filter event in the middle of the time series
        inlim = [idx(1)-round(szFilt*60*4/1.92) idx(1)+round(szFilt*60*4/1.92)];
      end    
      BEtemp = m_dbbp(inlim(1):inlim(2),:);
      AFtemp = m_dbbp(inlim(1):inlim(2),:);
      dtBEtmp = dt(inlim(1):inlim(2),:);
      dtAFtmp = dt(inlim(1):inlim(2),:);

      highpk = dtBEtmp(BEtemp > 0.7 * max(BEtemp)); % local high peak > 70% of max local peak of derivative
      highpk = highpk(highpk < cent); % keep only local high peak before filter event center
      popohpk = highpk(abs(dt(x(k))-highpk) == min(abs(dt(x(k))-highpk)));
      if all(isempty(popohpk) & x(k)<=round(szFilt*60*4/1.25))
        filt_st (k) = dtBEtmp(1);
      elseif isempty(popohpk) % find filter event longer than 10 in case flowcontrol is stuck
        foo = 2400;
        while all(max(BEtemp) < 0.03...% nanstd(BEtemp) < 12*nanstd(BEtemp(dtBEtmp > cent - minutes(0.5) & dtBEtmp < cent)) ...
                & foo < 30000 ...
                & size(dtBEtmp,1) < size(dt(1:idx(1)-foo),1))
          BEtemp = m_dbbp(inlim(1,1)-foo:inlim(1,2)-foo,:);
          dtBEtmp = dt(inlim(1,1)-foo:inlim(1,2)-foo,:);
          foo = foo + foo;
        end
%         highpk = dtBEtmp(BEtemp > 0.7 * max(BEtemp)); % local high peak > 70% of max local peak of derivative
        highpk = dtBEtmp(BEtemp > 0.7 * max(BEtemp) ...
          & max(BEtemp) > 0.03); % nanstd(BEtemp) > 12*nanstd(BEtemp(dtBEtmp > cent - minutes(0.5) & dtBEtmp < cent)));
        highpk = highpk(highpk < cent); % keep only local high peak before filter event center
        popohpk = highpk(abs(cent-highpk) == min(abs(cent-highpk)));
        if isempty(popohpk)
          filt_st (k) = NaT;
        else
          filt_st (k) = (popohpk - seconds(60));
        end
%         filt_st (k) = NaT;
      else
        filt_st (k) = (popohpk - seconds(60));
      end

      lowpk = dtAFtmp(AFtemp < 0.7 * min(AFtemp)); % local low peak < 70% of min local peak of derivative
      lowpk = lowpk(lowpk > cent); % keep only local low peak after filter event center
      popolpk = lowpk(abs(dt(x(k))-lowpk) == min(abs(dt(x(k))-lowpk)));
      if all(isempty(popolpk) & x(k)>=size(m_data,1)-round(szFilt*60*4/1.25))
        filt_st (k) = dtAFtmp(end);
      elseif isempty(popolpk) % find filter event longer than 10 in case flowcontrol is stuck
        foo = 2400;
        while all(min(AFtemp) > - 0.03...% nanstd(AFtemp) < 12*nanstd(AFtemp(dtAFtmp > cent & dtAFtmp < cent + minutes(0.5))) ...
                & foo < 30000 ...
                & size(dtAFtmp,1) < size(dt(idx(1):end-foo),1))
          AFtemp = m_dbbp(inlim(1,1)+foo:inlim(1,2)+foo,:);
          dtAFtmp = dt(inlim(1,1)+foo:inlim(1,2)+foo,:);
          foo = foo + foo;
        end
%         lowpk = dtAFtmp(AFtemp < 0.7 * min(AFtemp)); % local low peak < 70% of min local peak of derivative
        lowpk = dtAFtmp(AFtemp < 0.7 * min(AFtemp) ...
          & min(AFtemp) < - 0.03);% nanstd(AFtemp) > 12*nanstd(AFtemp(dtAFtmp > cent & dtAFtmp < cent + minutes(0.5))));
        lowpk = lowpk(lowpk > cent); % keep only local low peak after filter event "x"
        popolpk = lowpk(abs(cent-lowpk) == min(abs(cent-lowpk)));
        if isempty(popolpk)
          filt_end (k) = NaT;
        else
          filt_end (k) = (popolpk - seconds(0));
        end
%         filt_end (k) = NaT;
      else
        filt_end (k) = (popolpk - seconds(0));
      end
    end

% AC9 conversion from szFilt to nb of points (initial method)
% szFilt*60*4*1.2 % 2880
% szFilt*60*4/1.28 % 1875

  elseif contains(instrument, 'AC9') % for ACS
    parfor k = 1:size(x,1)
      if all(x(k)<=round(szFilt*60*4*1.2) & x(k)>=size(m_data,1)-round(szFilt*60*4*1.2))% when seg is very small
        outlim = [1 size(m_M,1)];
      elseif x(k)<=round(szFilt*60*4*1.2) % when filter event in the beginning of the time series
        outlim = [1 x(k)+round(szFilt*60*4*1.2)];
      elseif x(k)>=size(m_data,1)-round(szFilt*60*4*1.2) % when filter event at the end of the time series
        outlim = [x(k)-round(szFilt*60*4*1.2) size(m_M,1)];
      else % when filter event in the middle of the time series
        outlim = [x(k)-round(szFilt*60*4*1.2) x(k)+round(szFilt*60*4*1.2)];
      end
      sel = m_M(outlim(1):outlim(2));
      seldt = dt(outlim(1):outlim(2));
      cent = median(seldt(sel<mean(sel, 'omitnan')));
      idx = find(abs(dt-cent)<seconds(1));
      if isempty(idx)
        idx = floor(median(find(abs(dt-cent)<seconds(3))));
      end

      if all(idx(1)<=round(szFilt*60*4/1.28) & idx(1)>=size(m_data,1)-round(szFilt*60*4/1.28))% when seg is very small
        inlim = [1 size(m_M,1)];
      elseif idx(1)<=round(szFilt*60*4/1.28) % when filter event in the beginning of the time series
        inlim = [1 idx(1)+round(szFilt*60*4/1.28)];
      elseif idx(1)>=size(m_data,1)-round(szFilt*60*4/1.28) % when filter event at the end of the time series
        inlim = [idx(1)-round(szFilt*60*4/1.28) size(m_M,1)];
      else % when filter event in the middle of the time series
        inlim = [idx(1)-round(szFilt*60*4/1.28) idx(1)+round(szFilt*60*4/1.28)];
      end    
      BEtemp = m_dbbp(inlim(1):inlim(2),:);
      AFtemp = m_dbbp(inlim(1):inlim(2),:);
      dtBEtmp = dt(inlim(1):inlim(2),:);
      dtAFtmp = dt(inlim(1):inlim(2),:);

      highpk = dtBEtmp(BEtemp > 0.7 * max(BEtemp)); % local high peak > 70% of max local peak of derivative
      highpk = highpk(highpk < cent); % keep only local high peak before filter event center
      popohpk = highpk(abs(dt(x(k))-highpk) == min(abs(dt(x(k))-highpk)));
      if all(isempty(popohpk) & x(k)<=round(szFilt*60*4*1.2))
        filt_st (k) = dtBEtmp(1);
      elseif isempty(popohpk) % find filter event longer than 10 in case flowcontrol is stuck
        foo = 2400;
        while all(max(BEtemp) < 0.03...% nanstd(BEtemp) < 12*nanstd(BEtemp(dtBEtmp > cent - minutes(0.5) & dtBEtmp < cent)) ...
                & foo < 30000 ...
                & size(dtBEtmp,1) < size(dt(1:idx(1)-foo),1))
          BEtemp = m_dbbp(inlim(1,1)-foo:inlim(1,2)-foo,:);
          dtBEtmp = dt(inlim(1,1)-foo:inlim(1,2)-foo,:);
          foo = foo + foo;
        end
%         highpk = dtBEtmp(BEtemp > 0.7 * max(BEtemp)); % local high peak > 70% of max local peak of derivative
        highpk = dtBEtmp(BEtemp > 0.7 * max(BEtemp) ...
          & max(BEtemp) > 0.03); % nanstd(BEtemp) > 12*nanstd(BEtemp(dtBEtmp > cent - minutes(0.5) & dtBEtmp < cent)));
        highpk = highpk(highpk < cent); % keep only local high peak before filter event center
        popohpk = highpk(abs(cent-highpk) == min(abs(cent-highpk)));
        if isempty(popohpk)
          filt_st (k) = NaT;
        else
          filt_st (k) = (popohpk - seconds(55));
        end
%         filt_st (k) = NaT;
      else
        filt_st (k) = (popohpk - seconds(55));
      end

      lowpk = dtAFtmp(AFtemp < 0.7 * min(AFtemp)); % local low peak < 70% of min local peak of derivative
      lowpk = lowpk(lowpk > cent); % keep only local low peak after filter event cent
      popolpk = lowpk(abs(dt(x(k))-lowpk) == min(abs(dt(x(k))-lowpk)));
      if all(isempty(popolpk) & x(k)>=size(m_data,1)-round(szFilt*60*4*1.2))
        filt_st (k) = dtAFtmp(end);
      elseif isempty(popolpk) % find filter event longer than 10 in case flowcontrol is stuck
        foo = 2400;
        while all(min(AFtemp) > - 0.03...% nanstd(AFtemp) < 12*nanstd(AFtemp(dtAFtmp > cent & dtAFtmp < cent + minutes(0.5))) ...
                & foo < 30000 ...
                & size(dtAFtmp,1) < size(dt(idx(1):end-foo),1))
          AFtemp = m_dbbp(inlim(1,1)+foo:inlim(1,2)+foo,:);
          dtAFtmp = dt(inlim(1,1)+foo:inlim(1,2)+foo,:);
          foo = foo + foo;
        end
%         lowpk = dtAFtmp(AFtemp < 0.7 * min(AFtemp)); % local low peak < 70% of min local peak of derivative
        lowpk = dtAFtmp(AFtemp < 0.7 * min(AFtemp) ...
          & min(AFtemp) < - 0.03);% nanstd(AFtemp) > 12*nanstd(AFtemp(dtAFtmp > cent & dtAFtmp < cent + minutes(0.5))));
        lowpk = lowpk(lowpk > cent); % keep only local low peak after filter event "x"
        popolpk = lowpk(abs(cent-lowpk) == min(abs(cent-lowpk)));
        if isempty(popolpk)
          filt_end (k) = NaT;
        else
          filt_end (k) = (popolpk - seconds(0));
        end
%         filt_end (k) = NaT;
      else
        filt_end (k) = (popolpk - seconds(0));
      end
    end

% BB conversion from szFilt to nb of points (initial method)
% szFilt*60/1.5 % 400points ~6.7min
% szFilt*60*0.9167 % 550points ~9min
% szFilt*60*0.55 % 330points 5.5min

  elseif contains(instrument, 'BB3') % for BB3
    parfor k = 1:size(x,1)
%       for k = 1:size(x,1)

      if all(x(k)<=round(szFilt*60*0.9167) & x(k)>=size(m_data,1)-round(szFilt*60/1.5))% when seg is very small
        outlim = [1 size(m_M,1)];
      elseif x(k)<=round(szFilt*60*0.9167) % when filter event in the beginning of the time series
        outlim = [1 x(k)+round(szFilt*60/1.5)];
      elseif x(k)>=size(m_data,1)-round(szFilt*60/1.5) % when filter event at the end of the time series
        outlim = [x(k)-round(szFilt*60*0.9167) size(m_M,1)];
      else % when filter event in the middle of the time series
        outlim = [x(k)-round(szFilt*60*0.9167) x(k)+round(szFilt*60/1.5)];
      end
      sel = m_M(outlim(1,1):outlim(1,2));
      seldt = dt(outlim(1,1):outlim(1,2));
      cent = median(seldt(sel<mean(sel, 'omitnan')));
      idx = find(abs(dt-cent)<seconds(1));
      if isempty(idx)
        idx = floor(median(find(abs(dt-cent)<seconds(3))));
      end

      if all(idx(1)<=round(szFilt*60*0.55) & idx(1)>=size(m_data,1)-round(szFilt*60*0.55))% when seg is very small
        inlim = [1 size(m_M,1)];
      elseif idx(1)<=round(szFilt*60*0.55) % when filter event in the beginning of the time series
        inlim = [1 idx(1)+round(szFilt*60*0.55)];
      elseif idx(1)>=size(m_data,1)-round(szFilt*60*0.55) % when filter event at the end of the time series
        inlim = [idx(1)-round(szFilt*60*0.55) size(m_M,1)];
      else % when filter event in the middle of the time series
        inlim = [idx(1)-round(szFilt*60*0.55) idx(1)+round(szFilt*60*0.55)];
      end    
      BEtemp = m_dbbp(inlim(1,1):inlim(1,2),:);
      AFtemp = m_dbbp(inlim(1,1):inlim(1,2),:);
      dtBEtmp = dt(inlim(1,1):inlim(1,2),:);
      dtAFtmp = dt(inlim(1,1):inlim(1,2),:);

%       highpk = dtBEtmp(BEtemp > 0.7 * max(BEtemp)); % local high peak > 70% of max local peak of derivative
      highpk = dtBEtmp(BEtemp > 0.7 * max(BEtemp) ...
        & max(BEtemp) > 0.03); % nanstd(BEtemp) > 12*nanstd(BEtemp(dtBEtmp > cent - minutes(0.5) & dtBEtmp < cent)));
      highpk = highpk(highpk < cent); % keep only local high peak before filter event center
      popohpk = highpk(abs(cent-highpk) == min(abs(cent-highpk)));
      if all(isempty(popohpk) & idx<=round(szFilt*60*0.55))
        filt_st (k) = dtBEtmp(1);
      elseif isempty(popohpk) % find filter event longer than 10 in case flowcontrol is stuck
        foo = 600;
        while all(max(BEtemp) < 0.03...% nanstd(BEtemp) < 12*nanstd(BEtemp(dtBEtmp > cent - minutes(0.5) & dtBEtmp < cent)) ...
                & foo < 30000 ...
                & size(dtBEtmp,1) < size(dt(1:idx(1)-foo),1))
          BEtemp = m_dbbp(inlim(1,1)-foo:inlim(1,2)-foo,:);
          dtBEtmp = dt(inlim(1,1)-foo:inlim(1,2)-foo,:);
          foo = foo + foo;
        end
%         highpk = dtBEtmp(BEtemp > 0.7 * max(BEtemp)); % local high peak > 70% of max local peak of derivative
        highpk = dtBEtmp(BEtemp > 0.7 * max(BEtemp) ...
          & max(BEtemp) > 0.03); % nanstd(BEtemp) > 12*nanstd(BEtemp(dtBEtmp > cent - minutes(0.5) & dtBEtmp < cent)));
        highpk = highpk(highpk < cent); % keep only local high peak before filter event center
        popohpk = highpk(abs(cent-highpk) == min(abs(cent-highpk)));
        if isempty(popohpk)
          filt_st (k) = NaT;
        else
          filt_st (k) = (popohpk - seconds(190));
        end
%         filt_st (k) = NaT;
      else
        filt_st (k) = (popohpk - seconds(190));
      end

%       lowpk = dtAFtmp(AFtemp < 0.7 * min(AFtemp)); % local low peak < 70% of min local peak of derivative
      lowpk = dtAFtmp(AFtemp < 0.7 * min(AFtemp) ...
        & min(AFtemp) < - 0.03);% nanstd(AFtemp) > 12*nanstd(AFtemp(dtAFtmp > cent & dtAFtmp < cent + minutes(0.5))));
      lowpk = lowpk(lowpk > cent); % keep only local low peak after filter event center
      popolpk = lowpk(abs(cent-lowpk) == min(abs(cent-lowpk)));
      if all(isempty(popolpk) & idx>=size(m_data,1)-round(szFilt*60*0.55))
        filt_st (k) = dtAFtmp(end);
      elseif isempty(popolpk) % find filter event longer than 10 in case flowcontrol is stuck
        foo = 600;
        while all(min(AFtemp) > - 0.03...% nanstd(AFtemp) < 12*nanstd(AFtemp(dtAFtmp > cent & dtAFtmp < cent + minutes(0.5))) ...
                & foo < 30000 ...
                & size(dtAFtmp,1) < size(dt(idx(1):end-foo),1))
          AFtemp = m_dbbp(inlim(1,1)+foo:inlim(1,2)+foo,:);
          dtAFtmp = dt(inlim(1,1)+foo:inlim(1,2)+foo,:);
          foo = foo + foo;
        end
%         lowpk = dtAFtmp(AFtemp < 0.7 * min(AFtemp)); % local low peak < 70% of min local peak of derivative
        lowpk = dtAFtmp(AFtemp < 0.7 * min(AFtemp) ...
          & min(AFtemp) < - 0.03);% nanstd(AFtemp) > 12*nanstd(AFtemp(dtAFtmp > cent & dtAFtmp < cent + minutes(0.5))));
        lowpk = lowpk(lowpk > cent); % keep only local low peak after filter event "x"
        popolpk = lowpk(abs(cent-lowpk) == min(abs(cent-lowpk)));
        if isempty(popolpk)
          filt_end (k) = NaT;
        else
          filt_end (k) = (popolpk - seconds(-10));
        end
%         filt_end (k) = NaT;
      else
        filt_end (k) = (popolpk - seconds(-10));
      end
    end
    
% HBB conversion from szFilt to nb of points (initial method)
% szFilt*60*0.0283 % 21
% szFilt*60*0.038 % 21
% szFilt*60*0.046 % 28
% szFilt*60*0.0283 % 17
% szFilt*60*0.012 % 7.2

  elseif contains(instrument, 'HBB')
    for k = 1:size(x,1)
%       for k = 1:size(x,1)
      if all(x(k)<=round(szFilt*60*0.046) & x(k)>=size(m_data,1)-round(szFilt*60*0.038))% when seg is very small
        outlim = [1 size(m_M,1)];
      elseif x(k)<=round(szFilt*60*0.046) % when filter event in the beginning of the time series
        outlim = [1 round(x(k)+szFilt*60*0.038)];
      elseif x(k)>=size(m_data,1)-round(szFilt*60*0.038) % when filter event at the end of the time series
        outlim = [x(k)-round(szFilt*60*0.046) size(m_M,1)];
      else % when filter event in the middle of the time series
        outlim = [x(k)-round(szFilt*60*0.046) x(k)+round(szFilt*60*0.038)];
      end
      sel = m_M(outlim(1,1):outlim(1,2));
      seldt = dt(outlim(1,1):outlim(1,2));
      cent = median(seldt(sel<mean(sel, 'omitnan')));
      idx = find(abs(dt-cent)<seconds(10));
      if isempty(idx)
        idx = floor(median(find(abs(dt-cent)<seconds(30))));
      end

      if all(idx(1)<=round(szFilt*60*0.0283) & idx(1)>=size(m_data,1)-round(szFilt*60*0.012))% when seg is very small
        inlim = [1 size(m_M,1)];
      elseif idx(1)<=round(szFilt*60*0.0283) % when filter event in the beginning of the time series
        inlim = [1 idx(1)+round(szFilt*60*0.012)];
      elseif idx(1)>=size(m_data,1)-round(szFilt*60*0.012) % when filter event at the end of the time series
        inlim = [idx(1)-round(szFilt*60*0.012) size(m_M,1)];
      else % when filter event in the middle of the time series
        inlim = [idx(1)-round(szFilt*60*0.0283) idx(1)+round(szFilt*60*0.012)];
      end
      AFtemp = m_data(idx(1):outlim(2));
      dtAFtmp = dt(idx(1):outlim(2));
      
      lowpk = dtAFtmp(AFtemp == min(AFtemp) & movmedian(AFtemp,20) < median(AFtemp) / 1.2);
      if all(isempty(lowpk) & idx>=size(m_data,1)-round(szFilt*60*0.0283))
        filt_st (k) = dtAFtmp(end);
      elseif isempty(lowpk) % find filter event longer than 10 in case flowcontrol is stuck
        foo = 10;
        while ~any(AFtemp == min(AFtemp)...
                & movmedian(AFtemp,20) < median(AFtemp) / 2)  ...
                && foo < 30000 ...
                && size(dtAFtmp,1) < size(dt(idx(1):end-foo),1)
          AFtemp = m_data(idx(1):outlim(2)+foo);
          dtAFtmp = dt(idx(1):outlim(2)+foo,:);
          foo = foo + foo;
        end
        lowpk = dtAFtmp(AFtemp == min(AFtemp) & movmedian(AFtemp,20) < median(AFtemp) / 1.2);
        if isempty(lowpk)
          filt_end (k) = NaT;
        else
          filt_end (k) = (lowpk - seconds(5));
        end
%         filt_end (k) = NaT;
      else
        filt_end (k) = (lowpk - seconds(5));
      end
      
      if ~isempty(filt_end (k)) % if low pk is found then look for high pk
        idlowpk = find(abs(dt - filt_end (k)) == min(abs(dt - filt_end (k))));
        BEtemp = m_data(outlim(1):idlowpk);
        dtBEtmp = dt(outlim(1):idlowpk);

        highpk = dtBEtmp(BEtemp == max(BEtemp) & movmedian(BEtemp,20) > 1.2 * median(BEtemp));      
        if all(isempty(highpk) & idx<=round(szFilt*60*0.0283))
          filt_st (k) = dtBEtmp(1);
        elseif isempty(highpk) % find filter event longer than 10 in case flowcontrol is stuck
          foo = 10;
          while ~any(BEtemp == max(BEtemp)...
                  & movmedian(BEtemp,20) > 2 * median(BEtemp)) ...
                  && foo < 30000 ...
                  && size(dtBEtmp,1) < size(dt(1:idx(1)-foo),1)
            BEtemp = m_data(outlim(1)-foo:find(idlowpk));
            dtBEtmp = dt(outlim(1)-foo:find(idlowpk));
            foo = foo + foo;
          end
          highpk = dtBEtmp(BEtemp == max(BEtemp) & movmedian(BEtemp,20) > 2 * median(BEtemp)); 
          if isempty(highpk)
            filt_st (k) = NaT;
          else
            filt_st (k) = (highpk - seconds(50));
          end
  %         filt_st (k) = NaT;
        else
          filt_st (k) = (highpk - seconds(50));
        end
      else
        filt_st (k) = NaT;
      end
      
%       clf(63)
%       figure(63)
%       yyaxis('left'); hold on
%       plot(dtAFtmp, AFtemp,'k')
%       plot(dtBEtmp, BEtemp,'Color', [0.2 0.2 0.2])
%       scatter(cent, mean(sel),30,'r')
%       plot([dt(outlim(1)) dt(outlim(1))], [min(m_data(outlim(1):outlim(2))) max(m_data(outlim(1):outlim(2)))], ...
%         'g','Marker','none','LineStyle','-')
%       plot([dt(outlim(2)) dt(outlim(2))], [min(m_data(outlim(1):outlim(2))) max(m_data(outlim(1):outlim(2)))], ...
%         'g','Marker','none','LineStyle','-')
%       plot([dt(idx(1)) dt(idx(1))], [min(m_data(outlim(1):outlim(2))) max(m_data(outlim(1):outlim(2)))], ...
%         'b','Marker','none','LineStyle','-')
%       plot([filt_st(k) filt_st(k)], [min(m_data(outlim(1):outlim(2))) max(m_data(outlim(1):outlim(2)))], ...
%         'r','Marker','none','LineStyle','-')
%       plot([filt_end(k) filt_end(k)], [min(m_data(outlim(1):outlim(2))) max(m_data(outlim(1):outlim(2)))], ...
%         'r','Marker','none','LineStyle','-')
% %       yyaxis('right'); hold on
% %       plot(dt,dbbp,'r','Marker','none','LineStyle','-')
%       xlim([min(dtBEtmp) max(dtAFtmp)]);
%       
%       
%       scatter(seldt, sel,1,'b'); 
%       xlim([min(dtBEtmp) max(dtAFtmp)]);
%       yyaxis('right')
%       hold on
%       scatter(dt(dt > cent - minutes(0.8) & dt < cent + minutes(0.8)),m_dbbp(dt > cent - minutes(0.8) & dt < cent + minutes(0.8)),1,'o'); 
%       scatter(dt(inlim(1):inlim(2)),m_dbbp(inlim(1):inlim(2)),1,'r'); 
%       plot([filt_st(k) filt_st(k)],[min(m_dbbp(inlim(1):inlim(2))) max(m_dbbp(inlim(1):inlim(2)))],'g','Marker','none','LineStyle','-')
%       plot([filt_end(k) filt_end(k)],[min(m_dbbp(inlim(1):inlim(2))) max(m_dbbp(inlim(1):inlim(2)))],'g','Marker','none','LineStyle','-')

    end
      
  end

%     figure(64)
%     scatter(dt_ini,m_data_ini,1)
%     hold on
%     plot(datetime(FTH.dt,'ConvertFrom','datenum'),FTH.swt)
%     plot(dt,m_dbbp)
%     yyaxis('right')
%     plot([filt_st filt_st],[0 1],'g','Marker','none','LineStyle','-')
%     plot([filt_end filt_end],[0 1],'g','Marker','none','LineStyle','-')
% %     plot(datetime(FTH.dt,'ConvertFrom','datenum'),raw_swt,'r','LineStyle','-')
%     hold off

  filt = [filt_st,filt_end];

  %% keep ref only pair filt detected
  filt(any(isnat(filt),2),:) = [];

  % copy to FTH table
  for j=1:size(filt,1)
    FTH.swt(FTH.dt >= datenum(filt(j,1)) & FTH.dt <= datenum(filt(j,2))) = 1;
  end
end






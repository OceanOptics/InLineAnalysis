function [FTH] = SplitDetect (instrument, data, FTH, MinFiltPeriod)
% author: Guillaume Bourdin
% created: Oct 10, 2019

%Detect filter events on ACS and BB3 inline
%
% INPUT:
%   - instrument cell string of instrument name e.g. 'BB3' or 'ACS007'
%   - data_in: instrument data <N Table> time series of data that must contains:
%         - dt <Nx1 datenum> date and time precise to the second
%         for each other field the folowwing field must exist:
%         - bbp <NxM double> in column 2 for BB3
%         - a <NxM double> in column 2 for ACS
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
    error('missing argument: input instrument name, its data and FTH log');
end

% raw_swt = FTH.swt; % back-up FTH log
FTH.swt = zeros(size(FTH.swt,1),1);

dt_ini = datetime(table2array(data(:,1)),'ConvertFrom','datenum');
if any(contains(instrument, 'ACS') | contains(instrument, 'AC9'))
    data = table2array(data(:,3)); data(isinf(data)) = NaN;
    m_data_ini = nanmean(data,2);
elseif contains(instrument, 'BB3')
    data = table2array(data(:,2)); data(isinf(data)) = NaN;
    m_data_ini = nanmean(data,2);
else
    error('Instrument not recognized, check "QC Reference" in cfg file');
end

sep = diff(dt_ini);
seg = find(sep > minutes(10)); seg = [1; seg; size(dt_ini,1)];

for i = 1:size(seg,1)-1
m_data = m_data_ini(seg(i):seg(i+1));
dt = dt_ini(seg(i):seg(i+1));

try
if any(contains(instrument, 'ACS') | contains(instrument, 'AC9'))
%     m_M = movmean(m_data,200);
    m_M = sgolayfilt(m_data,1,201);
elseif contains(instrument, 'BB3')
    m_M = movmean(m_data,200);
end
catch
if any(contains(instrument, 'ACS') | contains(instrument, 'AC9'))
    m_M = movmean(m_data,200);
elseif contains(instrument, 'BB3')
    m_M = movmean(m_data,200);
end
end

dbbp = NaN(size(m_M,1),1);
for j = 1:size(m_M,1)-1; dbbp (j) = m_M(j)-m_M(j+1); end
m_dbbp = movmean(dbbp,100);

try
    if contains(instrument, 'ACS')
%         [~, x] = findpeaks(-movmean(m_M,1000),'MinPeakDistance',7000,'MinPeakWidth',1200); % 10000
        [~, x] = findpeaks(-sgolayfilt(m_M,1,401),'MinPeakDistance',MinFiltPeriod*0.6*240); % 10000  ,'MinPeakWidth',1200
    elseif contains(instrument, 'AC9')
%         [~, x] = findpeaks(-movmean(m_M,1000),'MinPeakDistance',7000,'MinPeakWidth',1200); % 10000
        [~, x] = findpeaks(-sgolayfilt(m_M,1,601),'MinPeakDistance',MinFiltPeriod*0.6*360); % 10000  ,'MinPeakWidth',1200
    elseif contains(instrument, 'BB3')
        [~, x] = findpeaks(-movmean(m_M,1000),'MinPeakDistance',MinFiltPeriod*0.47*60,'MinPeakWidth',1000); % 2000
    end
catch
    if any(contains(instrument, 'ACS') | contains(instrument, 'AC9'))
%         [~, x] = findpeaks(-movmean(m_M,1000),'MinPeakDistance',size(m_data,1)-2,'MinPeakWidth',1200); % 10000
        [~, x] = findpeaks(-movmean(m_M,1000),'MinPeakDistance',size(m_data,1)-2); % 10000  ,'MinPeakWidth',1200
    elseif contains(instrument, 'AC9')
%         [~, x] = findpeaks(-movmean(m_M,1000),'MinPeakDistance',7000,'MinPeakWidth',1200); % 10000
        [~, x] = findpeaks(-movmean(m_M,1000),'MinPeakDistance',size(m_data,1)-2); % 10000  ,'MinPeakWidth',1200
    elseif contains(instrument, 'BB3')
        [~, x] = findpeaks(-movmean(m_M,1000),'MinPeakDistance',size(m_data,1)-2,'MinPeakWidth',1000); % 2000
    end
end


% figure(62)
% % scatter(dt,m_Msg,15); hold on;
% scatter(dt,m_M,2,'filled'); hold
% scatter(dt,sgolayfilt(m_M,1,401),2,'g','filled');
% scatter(dt(x), nanmedian(m_M)*ones(size(x,1),1),10,'r','filled')
% yyaxis('right')
% scatter(dt,m_dbbp,2,'filled')

filt_st = NaT(size(x,1),1);
filt_end = NaT(size(x,1),1);
if contains(instrument, 'ACS') % for ACS
parfor k = 1:size(x,1)
    if all(x(k)<=1920 & x(k)>=size(m_data,1)-1920)% when seg is very small
        outlim = [1 size(m_M,1)];
    elseif x(k)<=1920 % when filter event in the beginning of the time series
        outlim = [1 x(k)+1920];
    elseif x(k)>=size(m_data,1)-1920 % when filter event at the end of the time series
        outlim = [x(k)-1920 size(m_M,1)];
    else % when filter event in the middle of the time series
        outlim = [x(k)-1920 x(k)+1920];
    end
    sel = m_M(outlim(1):outlim(2));
    seldt = dt(outlim(1):outlim(2));
    cent = median(seldt(sel<nanmean(sel)));
    idx = find(abs(dt-cent)<seconds(1));
    if isempty(idx)
        idx = floor(median(find(abs(dt-cent)<seconds(3))));
    end
    
    if all(idx(1)<=1250 & idx(1)>=size(m_data,1)-1250)% when seg is very small
        inlim = [1 size(m_M,1)];
    elseif idx(1)<=1250 % when filter event in the beginning of the time series
        inlim = [1 idx(1)+1250];
    elseif idx(1)>=size(m_data,1)-1250 % when filter event at the end of the time series
        inlim = [idx(1)-1250 size(m_M,1)];
    else % when filter event in the middle of the time series
        inlim = [idx(1)-1250 idx(1)+1250];
    end    
    BEtemp = m_dbbp(inlim(1):inlim(2),:);
    AFtemp = m_dbbp(inlim(1):inlim(2),:);
    dtBEtmp = dt(inlim(1):inlim(2),:);
    dtAFtmp = dt(inlim(1):inlim(2),:);
    
% % %     if x(k)<=1920 & x(k)>=size(m_data,1)-1920% when seg is very small
% % %     sel = m_M(1:end);
% % %     seldt = dt(1:end);
% % %     cent = median(seldt(sel<nanmean(sel)));
% % %     idx = find(abs(dt-cent)<seconds(1));
% % %     BEtemp = m_dbbp(1:idx(1),:);
% % %     AFtemp = m_dbbp(idx(1):end,:);
% % %     dtBEtmp = dt(1:idx(1),:);
% % %     dtAFtmp = dt(idx(1):end,:);
% % %     elseif x(k)<=1920 % when filter event in the beginning of the time series
% % %     sel = m_M(1:x(k)+1920);
% % %     seldt = dt(1:x(k)+1920);
% % %     cent = median(seldt(sel<nanmean(sel)));
% % %     idx = find(abs(dt-cent)<seconds(1));
% % %     BEtemp = m_dbbp(1:idx(1),:);
% % %     AFtemp = m_dbbp(idx(1):idx(1)+1250,:);
% % %     dtBEtmp = dt(1:idx(1),:);
% % %     dtAFtmp = dt(idx(1):idx(1)+1250,:);
% % %     elseif x(k)>=size(m_data,1)-1920 % when filter event at the end of the time series
% % %     sel = m_M(x(k)-1920:end);
% % %     seldt = dt(x(k)-1920:end);
% % %     cent = median(seldt(sel<nanmean(sel)));
% % %     idx = find(abs(dt-cent)<seconds(1));
% % %     BEtemp = m_dbbp(idx(1)-1250:idx(1),:);
% % %     AFtemp = m_dbbp(idx(1):end,:);
% % %     dtBEtmp = dt(idx(1)-1250:idx(1),:);
% % %     dtAFtmp = dt(idx(1):end,:);
% % %     else % when filter event in the middle of the time series
% % %     sel = m_M(x(k)-1920:x(k)+1920);
% % %     seldt = dt(x(k)-1920:x(k)+1920);
% % %     cent = median(seldt(sel<nanmean(sel)));
% % %     idx = find(abs(dt-cent)<seconds(1));
% % %     try
% % %     BEtemp = m_dbbp(idx(1)-1250:idx(1),:);
% % %     dtBEtmp = dt(idx(1)-1250:idx(1),:);
% % %     catch
% % %     end
% % %     try
% % %     AFtemp = m_dbbp(idx(1):idx(1)+1250,:);
% % %     dtAFtmp = dt(idx(1):idx(1)+1250,:);
% % %     catch
% % %     end
% % %     end
%     plot(seldt,sel); hold on; scatter(cent,nanmedian(sel),10,'r','filled')
%     plot(dtBEtmp,BEtemp); plot(dtAFtmp,AFtemp);

%     if std(BEtemp)<0.01*abs(sel(1000)) & std(AFtemp)<0.01*sel(1000)
%         continue
%     end
       
    highpk = dtBEtmp(BEtemp > 0.7 * max(BEtemp)); % local high peak > 10% of max local peak of derivative
    highpk = highpk(highpk < dt(x(k))); % keep only local high peak before filter event "x"
    popohpk = highpk(abs(dt(x(k))-highpk) == min(abs(dt(x(k))-highpk)));
    if all(isempty(popohpk) & x(k)<=1920)
        filt_st (k) = dtBEtmp(1);
    elseif isempty(popohpk)
        filt_st (k) = NaT;
    else
        filt_st (k) = (popohpk - seconds(55));
    end

    lowpk = dtAFtmp(AFtemp < 0.7 * min(AFtemp)); % local low peak < 10% of min local peak of derivative
    lowpk = lowpk(lowpk > dt(x(k))); % keep only local low peak after filter event "x"
    popolpk = lowpk(abs(dt(x(k))-lowpk) == min(abs(dt(x(k))-lowpk)));
    if all(isempty(popolpk) & x(k)>=size(m_data,1)-1920)
        filt_st (k) = dtAFtmp(end);
    elseif isempty(popolpk)
        filt_end (k) = NaT;
    else
        filt_end (k) = (popolpk - seconds(0));
    end
end

elseif contains(instrument, 'AC9') % for ACS
parfor k = 1:size(x,1)
    if all(x(k)<=2880 & x(k)>=size(m_data,1)-2880)% when seg is very small
        outlim = [1 size(m_M,1)];
    elseif x(k)<=2880 % when filter event in the beginning of the time series
        outlim = [1 x(k)+2880];
    elseif x(k)>=size(m_data,1)-2880 % when filter event at the end of the time series
        outlim = [x(k)-2880 size(m_M,1)];
    else % when filter event in the middle of the time series
        outlim = [x(k)-2880 x(k)+2880];
    end
    sel = m_M(outlim(1):outlim(2));
    seldt = dt(outlim(1):outlim(2));
    cent = median(seldt(sel<nanmean(sel)));
    idx = find(abs(dt-cent)<seconds(1));
    if isempty(idx)
        idx = floor(median(find(abs(dt-cent)<seconds(3))));
    end
    
    if all(idx(1)<=1875 & idx(1)>=size(m_data,1)-1875)% when seg is very small
        inlim = [1 size(m_M,1)];
    elseif idx(1)<=1875 % when filter event in the beginning of the time series
        inlim = [1 idx(1)+1875];
    elseif idx(1)>=size(m_data,1)-1875 % when filter event at the end of the time series
        inlim = [idx(1)-1875 size(m_M,1)];
    else % when filter event in the middle of the time series
        inlim = [idx(1)-1875 idx(1)+1875];
    end    
    BEtemp = m_dbbp(inlim(1):inlim(2),:);
    AFtemp = m_dbbp(inlim(1):inlim(2),:);
    dtBEtmp = dt(inlim(1):inlim(2),:);
    dtAFtmp = dt(inlim(1):inlim(2),:);
    
% % %     if x(k)<=1920 & x(k)>=size(m_data,1)-1920% when seg is very small
% % %     sel = m_M(1:end);
% % %     seldt = dt(1:end);
% % %     cent = median(seldt(sel<nanmean(sel)));
% % %     idx = find(abs(dt-cent)<seconds(1));
% % %     BEtemp = m_dbbp(1:idx(1),:);
% % %     AFtemp = m_dbbp(idx(1):end,:);
% % %     dtBEtmp = dt(1:idx(1),:);
% % %     dtAFtmp = dt(idx(1):end,:);
% % %     elseif x(k)<=1920 % when filter event in the beginning of the time series
% % %     sel = m_M(1:x(k)+1920);
% % %     seldt = dt(1:x(k)+1920);
% % %     cent = median(seldt(sel<nanmean(sel)));
% % %     idx = find(abs(dt-cent)<seconds(1));
% % %     BEtemp = m_dbbp(1:idx(1),:);
% % %     AFtemp = m_dbbp(idx(1):idx(1)+1250,:);
% % %     dtBEtmp = dt(1:idx(1),:);
% % %     dtAFtmp = dt(idx(1):idx(1)+1250,:);
% % %     elseif x(k)>=size(m_data,1)-1920 % when filter event at the end of the time series
% % %     sel = m_M(x(k)-1920:end);
% % %     seldt = dt(x(k)-1920:end);
% % %     cent = median(seldt(sel<nanmean(sel)));
% % %     idx = find(abs(dt-cent)<seconds(1));
% % %     BEtemp = m_dbbp(idx(1)-1250:idx(1),:);
% % %     AFtemp = m_dbbp(idx(1):end,:);
% % %     dtBEtmp = dt(idx(1)-1250:idx(1),:);
% % %     dtAFtmp = dt(idx(1):end,:);
% % %     else % when filter event in the middle of the time series
% % %     sel = m_M(x(k)-1920:x(k)+1920);
% % %     seldt = dt(x(k)-1920:x(k)+1920);
% % %     cent = median(seldt(sel<nanmean(sel)));
% % %     idx = find(abs(dt-cent)<seconds(1));
% % %     try
% % %     BEtemp = m_dbbp(idx(1)-1250:idx(1),:);
% % %     dtBEtmp = dt(idx(1)-1250:idx(1),:);
% % %     catch
% % %     end
% % %     try
% % %     AFtemp = m_dbbp(idx(1):idx(1)+1250,:);
% % %     dtAFtmp = dt(idx(1):idx(1)+1250,:);
% % %     catch
% % %     end
% % %     end
%     plot(seldt,sel); hold on; scatter(cent,nanmedian(sel),10,'r','filled')
%     plot(dtBEtmp,BEtemp); plot(dtAFtmp,AFtemp);

%     if std(BEtemp)<0.01*abs(sel(1000)) & std(AFtemp)<0.01*sel(1000)
%         continue
%     end
       
    highpk = dtBEtmp(BEtemp > 0.7 * max(BEtemp)); % local high peak > 10% of max local peak of derivative
    highpk = highpk(highpk < dt(x(k))); % keep only local high peak before filter event "x"
    popohpk = highpk(abs(dt(x(k))-highpk) == min(abs(dt(x(k))-highpk)));
    if all(isempty(popohpk) & x(k)<=2880)
        filt_st (k) = dtBEtmp(1);
    elseif isempty(popohpk)
        filt_st (k) = NaT;
    else
        filt_st (k) = (popohpk - seconds(55));
    end

    lowpk = dtAFtmp(AFtemp < 0.7 * min(AFtemp)); % local low peak < 10% of min local peak of derivative
    lowpk = lowpk(lowpk > dt(x(k))); % keep only local low peak after filter event "x"
    popolpk = lowpk(abs(dt(x(k))-lowpk) == min(abs(dt(x(k))-lowpk)));
    if all(isempty(popolpk) & x(k)>=size(m_data,1)-2880)
        filt_st (k) = dtAFtmp(end);
    elseif isempty(popolpk) 
        filt_end (k) = NaT;
    else
        filt_end (k) = (popolpk - seconds(0));
    end
end

elseif contains(instrument, 'BB3') % for BB3
parfor k = 1:size(x,1)
    
    if all(x(k)<=330 & x(k)>=size(m_data,1)-330)% when seg is very small
        outlim = [1 size(m_M,1)];
    elseif x(k)<=330 % when filter event in the beginning of the time series
        outlim = [1 x(k)+330];
    elseif x(k)>=size(m_data,1)-330 % when filter event at the end of the time series
        outlim = [x(k)-330 size(m_M,1)];
    else % when filter event in the middle of the time series
        outlim = [x(k)-330 x(k)+330];
    end
    sel = m_M(outlim(1):outlim(2));
    seldt = dt(outlim(1):outlim(2));
    cent = median(seldt(sel<nanmean(sel)));
    idx = find(abs(dt-cent)<seconds(1));
    if isempty(idx)
        idx = floor(median(find(abs(dt-cent)<seconds(3))));
    end
    
    if all(x(k)<=300 & x(k)>=size(m_data,1)-300)% when seg is very small
        inlim = [1 size(m_M,1)];
    elseif x(k)<=300 % when filter event in the beginning of the time series
        inlim = [1 x(k)+300];
    elseif x(k)>=size(m_data,1)-300 % when filter event at the end of the time series
        inlim = [x(k)-300 size(m_M,1)];
    else % when filter event in the middle of the time series
        inlim = [x(k)-300 x(k)+300];
    end    
    BEtemp = m_dbbp(inlim(1):inlim(2),:);
    AFtemp = m_dbbp(inlim(1):inlim(2),:);
    dtBEtmp = dt(inlim(1):inlim(2),:);
    dtAFtmp = dt(inlim(1):inlim(2),:);
    
    
% % %     if x(k)<=330 & x(k)>=size(m_data,1)-330% when seg is very small
% % %     sel = m_M(1:end);
% % %     seldt = dt(1:end);
% % %     cent = median(seldt(sel<nanmean(sel)));
% % %     idx = find(abs(dt-cent)<=seconds(1));
% % %     BEtemp = m_dbbp(1:idx(1),:);
% % %     AFtemp = m_dbbp(idx(1):end,:);
% % %     dtBEtmp = dt(1:idx(1),:);
% % %     dtAFtmp = dt(idx(1):end,:);
% % %     elseif x(k)<=330 % when filter event in the beginning of the time series
% % %     sel = m_M(1:x(k)+330);
% % %     seldt = dt(1:x(k)+330);
% % %     cent = median(seldt(sel<nanmean(sel)));
% % %     idx = find(abs(dt-cent)<=seconds(1));
% % %     BEtemp = m_dbbp(1:idx(1),:);
% % %     AFtemp = m_dbbp(idx(1):idx(1)+300,:);
% % %     dtBEtmp = dt(1:idx(1),:);
% % %     dtAFtmp = dt(idx(1):idx(1)+300,:);
% % %     elseif x(k)>=size(m_data,1)-330 % when filter event at the end of the time series
% % %     sel = m_M(x(k)-330:end);
% % %     seldt = dt(x(k)-330:end);
% % %     cent = median(seldt(sel<nanmean(sel)));
% % %     idx = find(abs(dt-cent)<seconds(1));
% % %     BEtemp = m_dbbp(idx(1)-300:idx(1),:);
% % %     AFtemp = m_dbbp(idx(1):end,:);
% % %     dtBEtmp = dt(idx(1)-300:idx(1),:);
% % %     dtAFtmp = dt(idx(1):end,:);
% % %     else % when filter event in the middle of the time series
% % %     sel = m_M(x(k)-330:x(k)+330);
% % %     seldt = dt(x(k)-330:x(k)+330);
% % %     cent = median(seldt(sel<nanmean(sel)));
% % %     idx = find(abs(dt-cent)<seconds(1));
% % %     BEtemp = m_dbbp(idx(1)-300:idx(1),:);
% % %     AFtemp = m_dbbp(idx(1):idx(1)+300,:);
% % %     dtBEtmp = dt(idx(1)-300:idx(1),:);
% % %     dtAFtmp = dt(idx(1):idx(1)+300,:);
% % %     end

%     if x(k)<=20 % when filter event in the beginning of the time series
%     BEtemp = m_dbbp(1:x(k),:);
%     AFtemp = m_dbbp(1:x(k)+300,:);
%     dtBEtmp = dt(1:x(k),:);
%     dtAFtmp = dt(1:x(k)+300,:);
%     elseif x(k)>20 & x(k)<=500 %#ok<*AND2> % when filter event in the beginning of the time series
%     BEtemp = m_dbbp(1:x(k),:);
%     AFtemp = m_dbbp(x(k)-20:x(k)+300,:);
%     dtBEtmp = dt(1:x(k),:);
%     dtAFtmp = dt(x(k)-20:x(k)+300,:);
%     elseif x(k)>=size(m_data,1)-300 % when filter event at the end of the time series
%     BEtemp = m_dbbp(x(k)-500:x(k),:);
%     AFtemp = m_dbbp(x(k)-20:size(m_data,1),:);
%     dtBEtmp = dt(x(k)-500:x(k),:);
%     dtAFtmp = dt(x(k)-20:size(m_data,1),:);
%     else % when filter event in the middle of the time series
%     BEtemp = m_dbbp(x(k)-500:x(k),:);
%     AFtemp = m_dbbp(x(k)-20:x(k)+300,:);
%     dtBEtmp = dt(x(k)-500:x(k),:);
%     dtAFtmp = dt(x(k)-20:x(k)+300,:);
%     end

    highpk = dtBEtmp(BEtemp > 0.7 * max(BEtemp)); % local high peak > 10% of max local peak of derivative
    highpk = highpk(highpk < dt(x(k))); % keep only local high peak before filter event "x"
    popohpk = highpk(abs(dt(x(k))-highpk) == min(abs(dt(x(k))-highpk)));
    if all(isempty(popohpk) & x(k)<=330)
        filt_st (k) = dtBEtmp(1);
    elseif isempty(popohpk)
        filt_st (k) = NaT;
    else
        filt_st (k) = (popohpk - seconds(3.5));
    end
        
%     if isempty(popohpk); filt_st (k) = NaT; else
%     filt_st (k) = (popohpk - minutes(3.5)); end

    lowpk = dtAFtmp(AFtemp < 0.7 * min(AFtemp)); % local low peak < 10% of min local peak of derivative
    lowpk = lowpk(lowpk > dt(x(k))); % keep only local low peak after filter event "x"
    popolpk = lowpk(abs(dt(x(k))-lowpk) == min(abs(dt(x(k))-lowpk)));
    if all(isempty(popolpk) & x(k)>=size(m_data,1)-330)
        filt_st (k) = dtAFtmp(end);
    elseif isempty(popolpk) 
        filt_end (k) = NaT;
    else
        filt_end (k) = (popolpk - seconds(-30));
    end
    
%     if isempty(popolpk); filt_end (k) = NaT; else
%     filt_end (k) = (popolpk - seconds(-30)); end
end
end

% figure(64)
% scatter(dt_ini,m_data_ini,1)
% hold on
% plot(datetime(FTH.dt,'ConvertFrom','datenum'),FTH.swt)
% plot(dt,m_dbbp)
% yyaxis('right')
% plot([filt_st filt_st],[0 1],'g','Marker','none','LineStyle','-')
% plot([filt_end filt_end],[0 1],'g','Marker','none','LineStyle','-')
% plot(datetime(FTH.dt,'ConvertFrom','datenum'),raw_swt,'r','LineStyle','-')
% hold off

filt = [filt_st,filt_end];

%% keep ref only pair filt detected
filt(any(isnat(filt),2),:) = [];

% copy to FTH table
for j=1:size(filt,1)
    FTH.swt(FTH.dt >= datenum(filt(j,1)) & FTH.dt <= datenum(filt(j,2))) = 1;
end
end
end






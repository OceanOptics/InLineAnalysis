function [FTH] = SplitDetect (instrument, data, FTH)
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

dt = datetime(table2array(data(:,1)),'ConvertFrom','datenum');
if contains(instrument, 'ACS')
    m_data = nanmean(table2array(data(:,3)),2);
elseif contains(instrument, 'BB3')
    m_data = nanmean(table2array(data(:,2)),2);
else
    error('Instrument not recognized, check "QC Reference" in cfg file');
end

if contains(instrument, 'ACS')
    m_M = movmean(m_data,200);
elseif contains(instrument, 'BB3')
    m_M = movmean(m_data,200);
end

dbbp = NaN(size(m_M,1),1);
for j = 1:size(m_M,1)-1; dbbp (j) = m_M(j)-m_M(j+1); end
m_dbbp = movmean(dbbp,100);

if contains(instrument, 'ACS')
    [~, x] = findpeaks(-m_M,'MinPeakDistance',10000);
elseif contains(instrument, 'BB3')
    [~, x] = findpeaks(-m_M,'MinPeakDistance',2000);
end

% figure(64)
% plot(dt(:,:),m_M(:,:)); hold on;
% scatter(dt(x), nanmean(m_M)*ones(size(x,1),1),10,'r','filled')
% yyaxis('right')
% plot(dt(:,:),m_dbbp(:,:))

filt_st = NaT(size(x,1),1);
filt_end = NaT(size(x,1),1);
for k = 1:size(x,1)-1
if contains(instrument, 'ACS') % for ACS
    
    if x(k)<=80 % when filter event in the beginning of the time series
    BEtemp = m_dbbp(1:x(k),:);
    AFtemp = m_dbbp(1:x(k)+2600,:);
    dtBEtmp = dt(1:x(k),:);
    dtAFtmp = dt(1:x(k)+2600,:);
    elseif x(k)>80 & x(k)<=2600 % when filter event in the beginning of the time series
    BEtemp = m_dbbp(1:x(k),:);
    AFtemp = m_dbbp(x(k)-80:x(k)+2600,:);
    dtBEtmp = dt(1:x(k),:);
    dtAFtmp = dt(x(k)-80:x(k)+2600,:);
    elseif x(k)>=size(m_data,1)-2600 % when filter event at the end of the time series
    BEtemp = m_dbbp(x(k)-2600:x(k),:);
    AFtemp = m_dbbp(x(k)-80:size(m_data,1),:);
    dtBEtmp = dt(x(k)-2600:x(k),:);
    dtAFtmp = dt(x(k)-80:size(m_data,1),:);
    else % when filter event in the middle of the time series
    BEtemp = m_dbbp(x(k)-2600:x(k),:);
    AFtemp = m_dbbp(x(k)-80:x(k)+2600,:);
    dtBEtmp = dt(x(k)-2600:x(k),:);
    dtAFtmp = dt(x(k)-80:x(k)+2600,:);
    end
    
    highpk = dtBEtmp(BEtemp > 0.7 * max(BEtemp)); % local high peak > 10% of max local peak of derivative
    highpk = highpk(highpk < dt(x(k))); % keep only local high peak before filter event "x"
    popohpk = highpk(abs(dt(x(k))-highpk) == min(abs(dt(x(k))-highpk)));
    filt_st (k) = (popohpk - seconds(45));

    lowpk = dtAFtmp(AFtemp < 0.7 * min(AFtemp)); % local low peak < 10% of min local peak of derivative
    lowpk = lowpk(lowpk > dt(x(k))); % keep only local low peak after filter event "x"
    popolpk = lowpk(abs(dt(x(k))-lowpk) == min(abs(dt(x(k))-lowpk)));
    filt_end (k) = (popolpk - seconds(0));

elseif contains(instrument, 'BB3') % for BB3
    
    if x(k)<=20 % when filter event in the beginning of the time series
    BEtemp = m_dbbp(1:x(k),:);
    AFtemp = m_dbbp(1:x(k)+300,:);
    dtBEtmp = dt(1:x(k),:);
    dtAFtmp = dt(1:x(k)+300,:);
    elseif x(k)>20 & x(k)<=500 %#ok<*AND2> % when filter event in the beginning of the time series
    BEtemp = m_dbbp(1:x(k),:);
    AFtemp = m_dbbp(x(k)-20:x(k)+300,:);
    dtBEtmp = dt(1:x(k),:);
    dtAFtmp = dt(x(k)-20:x(k)+300,:);
    elseif x(k)>=size(m_data,1)-300 % when filter event at the end of the time series
    BEtemp = m_dbbp(x(k)-500:x(k),:);
    AFtemp = m_dbbp(x(k)-20:size(m_data,1),:);
    dtBEtmp = dt(x(k)-500:x(k),:);
    dtAFtmp = dt(x(k)-20:size(m_data,1),:);
    else % when filter event in the middle of the time series
    BEtemp = m_dbbp(x(k)-500:x(k),:);
    AFtemp = m_dbbp(x(k)-20:x(k)+300,:);
    dtBEtmp = dt(x(k)-500:x(k),:);
    dtAFtmp = dt(x(k)-20:x(k)+300,:);
    end

    highpk = dtBEtmp(BEtemp > 0.7 * max(BEtemp)); % local high peak > 10% of max local peak of derivative
    highpk = highpk(highpk < dt(x(k))); % keep only local high peak before filter event "x"
    popohpk = highpk(abs(dt(x(k))-highpk) == min(abs(dt(x(k))-highpk)));
    filt_st (k) = (popohpk - minutes(3.5));

    lowpk = dtAFtmp(AFtemp < 0.7 * min(AFtemp)); % local low peak < 10% of min local peak of derivative
    lowpk = lowpk(lowpk > dt(x(k))); % keep only local low peak after filter event "x"
    popolpk = lowpk(abs(dt(x(k))-lowpk) == min(abs(dt(x(k))-lowpk)));
    filt_end (k) = (popolpk - seconds(-30));
end
end

filt_st(isnat(filt_st)) = [];
filt_end(isnat(filt_end)) = [];

% figure(64)
% scatter(dt,m_data,1)
% hold on
% plot(datetime(FTH.dt,'ConvertFrom','datenum'),FTH.swt)
% plot(dt,m_dbbp)
% yyaxis('right')
% plot([filt_st filt_st],[0 1],'g','Marker','none','LineStyle','-')
% plot([filt_end filt_end],[0 1],'g','Marker','none','LineStyle','-')
% plot(datetime(FTH.dt,'ConvertFrom','datenum'),raw_swt,'r','LineStyle','-')
% hold off

filt = table({filt_st},{filt_end});

%% keep ref only when the pair is in a 15 minutes interval
L = max([size(table2array(filt{:,1}),1) size(table2array(filt{:,2},1))]);
wh = [size(table2array(filt{:,1}),1) size(table2array(filt{:,2}),1)]==L;
if sum(wh)==2
    couple = NaT(size(table2array(filt{:,1}),1),2);
else
    couple = NaT(size(table2array(filt{:,~wh}),1),2);
end
for i=1:L
if sum(wh) < 2
    long = table2array(filt{:,wh});
    short = table2array(filt{:,~wh});
    cpl = find(abs(long(i)-short) < minutes(15));
    if ~isempty(cpl)
    couple (i,:) = [long(i) short(cpl(1))];
    end
else
    long = table2array(filt{:,1});
    short = table2array(filt{:,2});
    cpl = find(abs(long(i)-short) < minutes(15));
    if ~isempty(cpl)
    couple (i,:) = [long(i) short(cpl(1))];
    end
end
end

couple(any(isnat(couple),2),:) = [];

% copy to FTH table
FTH.swt = zeros(size(FTH.swt,1),1);
for j=1:size(couple,1)
if sum(wh) < 2
    if find(wh)==2
        FTH.swt(FTH.dt <= datenum(couple(j,1)) & FTH.dt >= datenum(couple(j,2))) = 1;
    elseif find(wh)==1
        FTH.swt(FTH.dt >= datenum(couple(j,1)) & FTH.dt <= datenum(couple(j,2))) = 1;
    end
else
    FTH.swt(FTH.dt >= datenum(couple(j,1)) & FTH.dt <= datenum(couple(j,2))) = 1;
end
end
end






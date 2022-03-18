function SimpleMap(data, station_info, colorbarlabel, varargin)
% author: Guillaume Bourdin
% created: August 13, 2020
%
% draws map with bubble map option if two variables are input
%
% INPUT:
%   - data: <Nx1 double> of values to map
%   - station_info: <Nx3 table> containing columns:
%       - station <Nx1 cellstr> of station ID
%       - lat <Nx1 double> of stations latitude
%       - lon <Nx1 double> of stations longitude
%   - colorbarlabel: <char> label for colorbar (with unit)
%
%   - optional argument:
%         - (1) <boolean> to activate bubble plot option (if true, another variable
%               is required as 5th input to the function)
%         - (2) <Nx1 double> of values to use as bubble size on map
%         - (3) <char> label for bubble legend (with unit)
%         - (4) <boolean> to print station name on map
%         - (5) <Nx3 double> text corresponding colors
%         - (6) <boolean> map onto biogeochemical provinces background
% OUTPUT: opens a figure
%
%%
switch nargin
    case 1
        error('not enough input arguments, provide station_info table with station# / lat / lon')
    case 2
        error('not enough input arguments, provide data label for colorbar (e.g. "Copepod abundance (#ind.m^{-3})"')
    case 3
        varargin{1} = false;
        varargin{2} = [];
        varargin{3} = [];
        varargin{4} = false;
        varargin{5} = [];
        varargin{6} = false;
    case 4
        varargin{3} = [];
        varargin{4} = false;
        varargin{5} = [];
        varargin{6} = false;
        if varargin{1} == true
            varargin{1} = false;
            varargin{2} = [];
            warning('Bubble plot boolean was true but biovolume data not provided, set to false (no bubble plot => standard scatter plot)')
        end
    case 5
        varargin{4} = false;
        varargin{5} = [];
        varargin{6} = false;
        if varargin{1} == true
            varargin{3} = 'You forgot to add the second variable label and unit';
            warning('variable 2 label and unit missing, check you variable 2 legend')
        end
    case 6
        varargin{4} = false;
        varargin{5} = [];
        varargin{6} = false;
    case 7
        varargin{6} = false;
        if varargin{4} == true
            varargin{5} = zeros(size(station_info,1), 3);
            warning('print text was true but no color provided, color set to black')
        end
    case 8
        varargin{6} = false;
    case 9
    otherwise
        error('Too many input arguments')
end

if varargin{4} == true && ~iscellstr(station_info.station)
    station_info.station = cellstr(station_info.station);
end

% if large data make edge transparent
if size(data, 1) > 300
  edgcol = 'none';
else
  edgcol = 'k';
end

figure('WindowState', 'maximize'); hold on
set(gca, 'Visible', 'off')

% wrap longitude over 360°
station_info.lon(station_info.lon < 0) = station_info.lon(station_info.lon < 0) + 360;
% get lat/lon limits
latlim = [min(station_info.lat) - (0.05 * (max(station_info.lat) - min(station_info.lat))) ...
    max(station_info.lat) + (0.05 * (max(station_info.lat) - min(station_info.lat)))];
lonlim = [min(station_info.lon) - (0.05 * (max(station_info.lon) - min(station_info.lon))) ...
    max(station_info.lon) + (0.05 * (max(station_info.lon) - min(station_info.lon)))];
% check if mapping toolbox is installed
if license('test','MAP_Toolbox')   
    axesm('robinson','MapLatLimit', latlim, 'MapLonLimit', lonlim,...
        'Frame','on','Grid','on', 'MLabelRound', 1,'PLabelRound', 1,...
        'MeridianLabel','on','ParallelLabel','on',...
        'MLineLocation',round((max(lonlim)-min(lonlim))/5, round(-log10((max(lonlim)-min(lonlim))/5))),...
        'PLineLocation',round((max(latlim)-min(latlim))/5, round(-log10((max(latlim)-min(latlim))/5))),...
        'MLabelLocation',round((max(lonlim)-min(lonlim))/5, round(-log10((max(lonlim)-min(lonlim))/5))),...
        'PLabelLocation',round((max(latlim)-min(latlim))/5, round(-log10((max(latlim)-min(latlim))/5))))
    
    % if biogeochemical provinces option == true
    if varargin{6}
        try 
            S=readtable('BIOGEOGRAPHY_05D.csv');
            lonbiogeo=reshape([S.Lon],360,720);
            lonbiogeo(:,1) = -180;
            lonbiogeo(:,720) = 180;
            latbiogeo=reshape([S.Lat],360,720);
            biogeo=reshape([S.BGCP],360,720);
            biogeo(:,1)=biogeo(:,2);
            biogeo(:,720)=biogeo(:,719);
            h=pcolorm(latbiogeo,lonbiogeo,biogeo);
            % shading flat
            set(h,'FaceAlpha',0.22)
            try 
                colormap(distinguishable_colors(56));
            catch
                warning('Distinguishable_colors function not found')
            end
        catch
            warning('Biogeochemical provinces file (BIOGEOGRAPHY_05D.csv) not found, biogeochemical provinces option set to false')
        end
    end
    
    % plot coastline
    land = shaperead('landareas.shp', 'UseGeoCoords', true);
    geoshow(land,'FaceColor', [0.9 0.9 0.9], 'EdgeColor', [0.6 0.6 0.6]);
    setm(gca,'ffacecolor',[1 1 1]);
    if ~varargin{1} && size(data,2) == 3 && all(~isnan(data),'all') && max(data(:)) <= 1
        scm(2) = scatterm(station_info.lat, ...
            station_info.lon, 150, ...
            data, ...
            'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,...
            'MarkerEdgeColor', edgcol, 'Marker', 'o');
            scm(2).Children.MarkerFaceAlpha = .5; %.75 .3
            scm(2).Children.MarkerEdgeAlpha = .8; %.15 .05
            scm(2).Children.LineWidth = 0.001;
    else
        % plot station with NaN in grey
        scm(1) = scatterm(station_info.lat(all(isnan(data),2)), ...
            station_info.lon(all(isnan(data),2)), 20, [0.65 0.65 0.65], ...
            'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1, ...
            'MarkerEdgeColor',[0.5 0.5 0.5], 'Marker', 'o');
        if varargin{1} && ~isempty(varargin{2}) % if bubble
            resized_data = exp(log10(varargin{2}) / max(log10(varargin{2}(:)))).^2.5 * 50 + 15;% normalize and adjust data range to bubble size range
            resized_data(varargin{2} == 0) = NaN;
            scm(2) = scatterm(station_info.lat(~all(isnan(data),2)), ...
                station_info.lon(~all(isnan(data),2)), ...
                resized_data(~all(isnan(data),2)), ...
                data(~all(isnan(data),2)), ...
                'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,...
                'MarkerEdgeColor',edgcol, 'Marker', 'o');
            scm(2).Children.MarkerFaceAlpha = .5; %.75 .3
            scm(2).Children.MarkerEdgeAlpha = .8; %.15 .05
            scm(2).Children.LineWidth = 0.001;
        else
%             data(data == 0) = NaN;
            scm(2) = scatterm(station_info.lat(~all(isnan(data),2)), ...
                station_info.lon(~all(isnan(data),2)), 150, ...
                data(~all(isnan(data),2)), ...
                'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1, ...
                'MarkerEdgeColor',edgcol, 'Marker', 'o');
            scm(2).Children.MarkerFaceAlpha = .75; %.75 .3 .5
            scm(2).Children.MarkerEdgeAlpha = .15; %.15 .05 .8
            scm(2).Children.LineWidth = 0.001;
        end
    end
    % print station name on map
    if varargin{4}
        distex_lat = (latlim(2)-latlim(1))*0.01;
        distex_lon = (lonlim(2)-lonlim(1))*0.005;
        for i = 1:size(station_info)
            textm(station_info.lat(i) + distex_lat, ...
                station_info.lon(i) + distex_lon, ...
                station_info.station(i), ...
                'Color', varargin{5}(i,:), 'FontWeight','bold','FontSize',14);
        end
    end
    
%     % scale
%     scaleruler('on');
%     setm(handlem('scaleruler1'), 'XLoc',0.9,'YLoc',-0.5, ... % 'MajorTick', [0,250], 'MinorTick', NaN,
%         'FontName', 'Times New Roman', 'FontSize', 10);

% if no mapping toolbox, check if m_map is installed
elseif exist('m_scatter', 'file')
    %% NEED TO be implemented
    
else
    error('Neither mapping toolbox, nor m_map were found')
end

if ~(~varargin{1} && size(data,2) == 3 && all(~isnan(data),'all') && max(data(:)) <= 1)
    colormap(brewermap(128, '*Spectral'));

    h = colorbar('eastoutside');
%     set(h, 'ylim', [min(data(:)) max(data(:))]);
    ylabel(h, colorbarlabel)
    if min(data(data>0)) / max(data(:)) > 0.01 % set logscale if data spans over 3 orders of magnitude
      extick = get(h, 'YTick');
      set(gca,'ColorScale','log')
      h.Ticks = extick;
    end
    set(gca,'FontSize',20)
end
% Legend bubble size map
if varargin{1} && ~isempty(varargin{1})% build logrithmic size scale
    n = floor(log(abs(varargin{2}))./log(10));
    n(isinf(n)) = 0;
    me = logspace(log10(10^min(n)), log10(10^max(n)), max(n)-min(n)+1);
    if size(me, 2) > 4 % reduce the legend when data spans over more than 5 order of magnitude
        foo = 1:ceil(size(me,2)/4):size(me,2);
        me = me(foo);
    end
    leg = legend('location', 'southeast'); % add legend to get position of the southeast corner
    legpos = get(leg,'position'); % get legend position
    delete(findobj('type','legend')) % delete legend
    southeast_corner = [legpos(1)+legpos(3) legpos(2)+legpos(4)];

    number_of_mark = 1;
    axes('Position',[southeast_corner(1)-0.03*number_of_mark*2 ...
        southeast_corner(2) .03*number_of_mark*2 0.045*size(me,2)])
    hold on
    y_space = 0.04;
    y = 0.12:y_space:0.12+y_space*(size(me,2)-1);

    x = repmat(0.07, 1, size(me,2));
    scatter(x,y, exp(log10(me) / max(log10(me))) .^ 2.5 * 50 + 15,'k','Marker','o');
%     x_spac = 0.095;
%     scatter(x,y, me/max(me)*400+100,'k','Marker','o'); % dynamic need to be implemented for multiple marker type (e.g. for 4 different marker)
%         scatter(x+x_spac,y,sqrt(me)*300+15,'k','Marker','*');
%         scatter(x+x_spac*2,y,sqrt(me)*300+15,'k','Marker','v');
%         scatter(x+x_spac*3,y,sqrt(me)*300+15,'k','Marker','o');
    xlim([-0.1 x(end)*10]); ylim([0.1 y(end) + y_space*3]); 
    for i = 1:size(me,2)
        text(x(i)*3.8,y(i),num2str(me(i)),'FontSize',12);
    end
    text(0,max(y)+0.04*1.7,varargin{3},'FontSize',14);
    set(gca,'xTick',[],'yTick',[]);
    box on
end

set(gcf,'PaperPositionMode','auto')


end


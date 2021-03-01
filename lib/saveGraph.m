function saveGraph(fig_name, savefmt, varargin)
% Author: Guillaume Bourdin
% Date: 2020-08-23
% 
% Save figure in different formats
%
%% Input:
%   - fig_name: <char> figure name
%   - savefmt: <char> save format ('jpeg', 'svg' or 'fig')
% Optional input:
%   - (1) <double> resolution in dpi (default = 200)
%
%%
if nargin < 2
    error('not enough input variable')
elseif nargin == 2
    varargin{1} = 200;
elseif nargin > 3
    error('Too many input variable')
end

vers = year(datetime(version('-date')));

switch savefmt
    case 'pdf'
        if vers >= 2020
            exportgraphics(gcf, [fig_name '.pdf'],'BackgroundColor','none','ContentType','vector')
        else
            error([savefmt ' format not supported'])
        end
    case {'jpg', 'jpeg'}
        if vers >= 2020
            exportgraphics(gcf, [fig_name '.jpg'],'Resolution',varargin{1})
        else
            set(gcf,'renderer','Painters');
            print(fig_name, '-djpeg', ['-r' num2str(varargin{1})]);
        end
    case 'svg'
        set(gcf,'renderer','Painters');
        print(fig_name, '-dsvg', ['-r' num2str(varargin{1})]);
    case 'fig'
        saveas(gcf, fig_name, savefmt);
    case 'png'
        if vers >= 2020
            exportgraphics(gcf, [fig_name '.png'])
        else
            error([savefmt ' format not supported'])
        end
    case {'tif','tiff'}
        if vers >= 2020
            exportgraphics(gcf, [fig_name '.tiff'])
        else
            error([savefmt ' format not supported'])
        end
    case 'emf'
        if vers >= 2020
            exportgraphics(gcf, [fig_name '.emf'])
        else
            error([savefmt ' format not supported'])
        end
    case 'eps'
        if vers >= 2020
            exportgraphics(gcf, [fig_name '.emf'])
        else
            error([savefmt ' format not supported'])
        end
    otherwise
        error([savefmt ' format not supported'])
end
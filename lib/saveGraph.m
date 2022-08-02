function saveGraph(fig_name, savefmt, res, hdl) % , legend_icons
% Author: Guillaume Bourdin
% Date: 2020-08-23
% 
% Harmonize command line to save graphics in different formats
%
%% Input:
%   - fig_name: <char> figure name
%   - savefmt: <char> save format ('pdf, 'jpeg', 'svg', 'fig', 'png', 'tiff', 'emf', or 'eps')
% Optional input:
%   - (1) <double> resolution in dpi (default = 200)
%   - (2) <figure handle> hdl (default = current handle)
%
%%
if nargin < 2
  error('not enough input variable')
elseif nargin == 2
  res = 200;
  hdl = gcf;
elseif nargin == 3
  hdl = gcf;
elseif nargin == 4 || nargin == 5
  if isempty(res)
    res = 200;
  end
else
  error('Too many input variable')
end

vers = year(datetime(version('-date')));
savefmt = strrep(savefmt, '-', '');
fig_name = strip(fig_name);

switch savefmt
  case 'pdf'
    if vers >= 2020
      set(hdl,'renderer','Painters');
      set(gcf,'PaperType','A4')
      exportgraphics(hdl, [fig_name '.pdf'],'BackgroundColor','none','ContentType','vector')
    else
      error([savefmt ' format not supported'])
    end
  case {'jpg', 'jpeg'}
    if vers >= 2020
      exportgraphics(hdl, [fig_name '.jpg'],'Resolution',res)
    else
      set(gcf,'renderer','Painters');
      print(fig_name, '-djpeg', ['-r' num2str(res)]);
    end
  case 'svg'
    set(hdl,'renderer','Painters');
%     if nargin == 5
%       fig2svg([fig_name '.' savefmt], '', '', legend_icons)
%     else
%       fig2svg([fig_name '.' savefmt])
%     end
    print(fig_name, '-dsvg', ['-r' num2str(res)]);
%     plot2svg([fig_name '.svg']);
  case 'fig'
    savefig(gcf, [fig_name '.' savefmt]);
  case 'png'
    if vers >= 2020
      exportgraphics(hdl, [fig_name '.png'])
    else
      error([savefmt ' format not supported'])
    end
  case {'tif','tiff'}
    if vers >= 2020
      exportgraphics(hdl, [fig_name '.tiff'])
    else
      error([savefmt ' format not supported'])
    end
  case 'emf'
    if vers >= 2020
      exportgraphics(hdl, [fig_name '.emf'])
    else
      error([savefmt ' format not supported'])
    end
  case 'eps'
    if vers >= 2020
      exportgraphics(hdl, [fig_name '.emf'])
    else
      error([savefmt ' format not supported'])
    end
  otherwise
    error([savefmt ' format not supported'])
end
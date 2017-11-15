function h = fig(id, title)
%FIG Select, reset and add a title to the figure corresponding to id
%
%Syntax:  fig(1, 'title');
%
%Inputs: 
%    Required:
%        id double containing the identification number of the figure
%    Optional:
%        title string containing the title of the figure    
%
%Outputs: go to the figure id
%
%Example: Plot Profiles in the Artic
% fig(1, 'title');
% plot(1:10, rand(10,1));
%
% Tested with: Matlab R2015a
%
% Author: Nils Haentjens, Ms, University of Maine
% Email: nils.haentjens@maine.edu
% Created: 14th July 2015
% Last update: 14th July 2015   

% Check Nargin
if nargin > 2
   error('Too many input arguments.')
elseif nargin < 1
   error('Not enough input arguments.')
end;

if nargin < 2
  fig_name = sprintf('%02d', id);
else
  fig_name = sprintf('%02d %s', id, title);
end;

% Select figure
h=figure(id);
% Reset/Empty figure
clf(id, 'reset');
% Set title
set(h, 'Name', fig_name, 'NumberTitle', 'off');


end
%  F = mkfig(num, width, height)
%
function F = mkfig(num, width, height, presentation)
if nargin <4
  presentation = false;
end
F = figure(num); clf
axes('FontName','Times New Roman') % Set axis font style
box('on'); % Define box around whole figure
hold on;

% First, set units to inches. Then get the (x0,y0) position, and only reset the
% width and height. Because it is annoying to always have the figure start at
% (0,0).
set(F, 'Units', 'inches');
pos = get(F, 'Position');

lft = pos(1);
bt = pos(2);

set(F, 'Units', 'Inches', 'Position', [lft, bt, width, height],...
    'PaperUnits', 'Inches', 'PaperSize', [width, height],...
    'InvertHardcopy', 'off');

if presentation
set(F, 'Color', fig_color);
else
  set(F, 'Color', 'w');
end



end
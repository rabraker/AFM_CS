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

set(F, 'Units', 'Inches', 'Position', [0, 0, width, height],...
    'PaperUnits', 'Inches', 'PaperSize', [width, height],...
    'InvertHardcopy', 'off');

if presentation
set(F, 'Color', fig_color);
else
  set(F, 'Color', 'w');
end



end
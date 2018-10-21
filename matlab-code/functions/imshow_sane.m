function imshow_sane(I, ax, width, height, lims)

I = flipud(I);
if ~exist('lims', 'var')
  lo = min(min(I));
  hi = max(max(I));
else
  lo = lims(1);
  hi = lims(2);
end
% hi = .6

xdata = [0, width];
ydata = [0, height];
imshow(flipud(I), [lo, hi], 'XData', xdata, 'YData', ydata, 'Parent', ax);
axis('on')
ax.YTickLabel = fliplr(yticks);

xlabel('x-dir')
ylabel('y-dir')




end


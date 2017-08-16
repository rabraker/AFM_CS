function imshow_sane(I, ax, width, height )

I = flipud(I);

lo = min(min(I));
hi = max(max(I));
xdata = [0, width];
ydata = [0, height];
imshow(flipud(I), [lo, hi], 'XData', xdata, 'YData', ydata, 'Parent', ax)
axis('on')
ax.YTickLabel = fliplr(yticks)

xlabel('x-dir')
ylabel('y-dir')




end

